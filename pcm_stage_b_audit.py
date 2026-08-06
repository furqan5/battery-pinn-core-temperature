"""Stage B audit -- the Stage B numbers looked too good, so they get audited.

Stage B said melted fraction is recoverable to ~0.1% of the full buffer and
latent heat to 0.015%.  Numbers that good on a 0.1 K sensor are a signal to look
for a leak, not to celebrate.  Three audits:

L1  HISTORY LEAK in S3/S4.  Stage B built the initial state from the reference
    trajectory at the NOMINAL q, then perturbed q only in the observation
    window.  A plant that had always run at a different q would have arrived at
    f = 0.5 with a different temperature field.  Holding the initial field fixed
    while varying q hands the estimator an anchor it does not have in the field.
    Fixed here by regenerating the trajectory at the trial q, so that q
    perturbs the history AND the window consistently.

L2  NOISE SCALING.  Conditioning is noise-independent; standard deviations scale
    linearly with sigma.  The brief fixes sigma = 0.1 K (instrument noise), but
    on a real cell the model-form error dominates and it is correlated, which
    costs effective sample size.  Reported, not buried.

L3  WHERE THE INFORMATION LIVES.  If dT_surf/df0 is essentially flat across the
    window, the trace carries one number (a level), not a rich time signature --
    and a monotone lookup from surface temperature to melted fraction would
    extract it without any inverse PDE solve.  This is a shape diagnostic on the
    sensitivity vector, not a fit; nothing is fitted in this session.
"""
import json
import numpy as np

import pcm_params as P
import pcm_identify as I

RESULTS = {}
T_WINDOW = 600.0
F0 = 0.50
T_REF_END = 9000.0     # enough to pass full melt at every q tested here


def banner(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)


# --------------------------------------------------------------------------- #
# L1 -- consistent-history parameterisation
# --------------------------------------------------------------------------- #

_CACHE = {}


def _ref_for(q, h):
    key = (round(q, 6), round(h, 9))
    if key not in _CACHE:
        _CACHE[key] = I.reference_trajectory(q=q, h=h, t_end=T_REF_END)
    return _CACHE[key]


def eval_consistent(theta, t_window=T_WINDOW, sensors=("surface",)):
    """Trial q and h drive BOTH the history that produced the state and the
    window that is observed.  This is the honest 'unknown constant load' case."""
    ref = _ref_for(theta["q"], theta["h"])
    T0, _ = I.state_at_f(ref, theta["f0"])
    y, _ = I.predict_window(ref["model"], T0, theta["q"], theta["h"],
                            t_window, P.SAMPLE_DT, sensors)
    return y


def l1_history_leak():
    banner("L1  History leak: does q perturb the state as well as the window?")
    base = {"f0": F0, "q": P.Q_NOMINAL, "h": P.H_CONV,
            "Lf": P.PCM_LF, "k_pcm": P.PCM_K_SOLID}
    ref_nom = _ref_for(P.Q_NOMINAL, P.H_CONV)

    ev_leaky = lambda th: I._eval_window(th, ref_nom, T_WINDOW, ("surface",))
    ev_fixed = lambda th: eval_consistent(th, T_WINDOW, ("surface",))

    out = {}
    for name, free in (("S3  {f, q}", ["f0", "q"]),
                       ("S4  {f, q, h}", ["f0", "q", "h"])):
        print(f"\n  --- {name} ---")
        row = {}
        for tag, ev in (("state fixed at nominal q  (Stage B, LEAKY)", ev_leaky),
                        ("state consistent with q   (corrected)", ev_fixed)):
            S, _, _ = I.sensitivity_matrix(base, free, ev)
            r = I.crlb_from_S(S, free, base)
            print(f"    {tag}")
            print(f"      cond(S) = {r['cond_scaled']:11.4e}   "
                  f"f sd = {100*r['sd_natural'][0]:8.3f} % of buffer")
            off = np.abs(r["corr"] - np.eye(len(free)))
            print(f"      max |off-diag corr| = {np.max(off):.4f}")
            row[tag.split("(")[0].strip()] = {
                "cond": r["cond_scaled"],
                "f_sd_buffer_pct": 100 * r["sd_natural"][0],
                "sd_relative": r["sd_relative"].tolist(),
                "corr": np.asarray(r["corr"]).tolist(),
                "singular": r["singular"],
            }
        out[name.split()[0]] = row
    RESULTS["L1"] = out
    return out


# --------------------------------------------------------------------------- #
# L2 -- noise scaling
# --------------------------------------------------------------------------- #

def l2_noise_scaling(ref):
    banner("L2  Noise level vs conditioning (trap 2)")
    base = {"f0": F0, "q": P.Q_NOMINAL, "h": P.H_CONV}
    ev = lambda th: eval_consistent(th, T_WINDOW, ("surface",))
    S, _, _ = I.sensitivity_matrix(base, ["f0", "q"], ev)

    print(f"  {'sigma (K)':>11s} {'cond(S)':>12s} {'f sd (% buffer)':>18s}"
          f"   note")
    rows = []
    notes = {0.1: "instrument noise, as specified in the brief",
             0.5: "plausible model-form error",
             1.0: "model-form error on a real cell"}
    for sig in (0.1, 0.5, 1.0):
        r = I.crlb_from_S(S, ["f0", "q"], base, sigma=sig)
        rows.append({"sigma": sig, "cond": r["cond_scaled"],
                     "f_sd_buffer_pct": 100 * r["sd_natural"][0]})
        print(f"  {sig:11.2f} {r['cond_scaled']:12.4e} "
              f"{100*r['sd_natural'][0]:18.4f}   {notes[sig]}")
    print("\n  cond(S) is IDENTICAL across rows -- conditioning is a property of")
    print("  the sensitivity matrix alone.  The standard deviations scale")
    print("  linearly with sigma.  These are different statements and the")
    print("  distinction decides how much of this survives on real hardware.")
    print()
    print("  Correlated residuals cost effective sample size on top of this:")
    n = int(T_WINDOW) + 1
    for tau in (1, 30, 120):
        neff = max(1.0, n / max(1.0, 2 * tau))
        infl = np.sqrt(n / neff)
        print(f"    correlation time {tau:4d} s -> N_eff ~ {neff:7.1f} "
              f"of {n}, bound inflated {infl:5.2f}x")
    RESULTS["L2"] = rows


# --------------------------------------------------------------------------- #
# L3 -- where the information lives
# --------------------------------------------------------------------------- #

def l3_information_shape(ref):
    banner("L3  Where does the information live? (shape of dT_surf/df0)")
    base = {"f0": F0, "q": P.Q_NOMINAL, "h": P.H_CONV}
    ev = lambda th: I._eval_window(th, ref, T_WINDOW, ("surface",))
    S, _, _ = I.sensitivity_matrix(base, ["f0"], ev)
    s = S[:, 0]

    mean = float(np.mean(s))
    var_frac = float(np.std(s) / abs(mean))
    # how much of the information survives if only the MEAN level is used
    n = len(s)
    info_full = float(s @ s)
    info_mean = float(n * mean ** 2)
    print(f"  window samples                 {n}")
    print(f"  mean dT_surf/df0               {mean:10.4f} K per unit fraction")
    print(f"  std / |mean| across window     {var_frac:10.4f}")
    print(f"  Fisher info, full trace        {info_full:12.4e}")
    print(f"  Fisher info, mean level only   {info_mean:12.4e}"
          f"   ({100*info_mean/info_full:.2f} % of full)")
    print()
    print("  A single 1 Hz sample, no window averaging:")
    sd1 = P.NOISE_SIGMA / abs(mean)
    print(f"    sd from one sample           {100*sd1:10.3f} % of buffer")
    print(f"    sd from {n} samples          "
          f"{100*sd1/np.sqrt(n):10.3f} % of buffer  (ideal 1/sqrt(N))")
    print()
    print("  If the mean level carries essentially all the information, the")
    print("  measurement is one number plus noise.  A monotone lookup from")
    print("  surface temperature to melted fraction would extract it with no")
    print("  inverse PDE solve at all -- which is the Part 7 finding again.")
    RESULTS["L3"] = {"n": n, "mean_sens": mean, "std_over_mean": var_frac,
                     "info_full": info_full, "info_mean_only": info_mean,
                     "info_mean_frac": info_mean / info_full,
                     "sd_one_sample_pct": 100 * sd1}


# --------------------------------------------------------------------------- #
# L4 -- is the signal the FRONT, or just the system warming up?
# --------------------------------------------------------------------------- #

def l4_front_vs_warming():
    banner("L4  Is the signal the melt front, or just monotone warming?")
    print("  Control: the same trajectory-position estimate on a cell whose PCM")
    print("  has NO latent heat (L_f -> ~0).  There is no front to track, so any")
    print("  remaining 'identifiability' is pure elapsed-time-on-a-trajectory.")
    print()
    ref_pcm = _ref_for(P.Q_NOMINAL, P.H_CONV)
    base = {"f0": F0, "q": P.Q_NOMINAL, "h": P.H_CONV}
    ev = lambda th: I._eval_window(th, ref_pcm, T_WINDOW, ("surface",))
    S, _, _ = I.sensitivity_matrix(base, ["f0"], ev)
    r = I.crlb_from_S(S, ["f0"], base)

    # sensible-only control: latent heat reduced 1000x, same mushy window, so
    # "melted fraction" degenerates to "position along a warming curve"
    mdl0 = P.build_model(Lf=P.PCM_LF / 1000.0)
    ref0 = I.reference_trajectory(model=mdl0, q=P.Q_NOMINAL, t_end=T_REF_END)
    ev0 = lambda th: I._eval_window(th, ref0, T_WINDOW, ("surface",))
    S0, _, _ = I.sensitivity_matrix(base, ["f0"], ev0)
    r0 = I.crlb_from_S(S0, ["f0"], base)

    f = ref_pcm["melted_fraction"]
    f2 = ref0["melted_fraction"]
    span = (float(ref_pcm["t"][int(np.searchsorted(f, 1e-9))]),
            float(ref_pcm["t"][int(np.searchsorted(f, 1 - 1e-9))]))
    span0 = (float(ref0["t"][int(np.searchsorted(f2, 1e-9))]),
             float(ref0["t"][int(np.searchsorted(f2, 1 - 1e-9))]))

    print(f"  {'case':28s} {'melt span (s)':>16s} {'dT/df0 (K)':>12s} "
          f"{'f sd (% buffer)':>17s}")
    print(f"  {'real PCM (L_f = 165 kJ/kg)':28s} "
          f"{span[1]-span[0]:16.0f} {float(np.mean(S[:,0])):12.4f} "
          f"{100*r['sd_natural'][0]:17.4f}")
    print(f"  {'no latent heat (L_f/1000)':28s} "
          f"{span0[1]-span0[0]:16.0f} {float(np.mean(S0[:,0])):12.4f} "
          f"{100*r0['sd_natural'][0]:17.4f}")
    print()
    print("  Read this as: how much of the melted-fraction information is")
    print("  specific to phase change, and how much is just 'the cell is")
    print("  warming monotonically so temperature encodes elapsed time'.")
    RESULTS["L4"] = {
        "pcm": {"melt_span_s": span[1] - span[0],
                "mean_sens": float(np.mean(S[:, 0])),
                "sd_buffer_pct": 100 * r["sd_natural"][0]},
        "no_latent": {"melt_span_s": span0[1] - span0[0],
                      "mean_sens": float(np.mean(S0[:, 0])),
                      "sd_buffer_pct": 100 * r0["sd_natural"][0]},
    }


if __name__ == "__main__":
    print("building nominal reference trajectory ...")
    ref = _ref_for(P.Q_NOMINAL, P.H_CONV)
    l1_history_leak()
    l2_noise_scaling(ref)
    l3_information_shape(ref)
    l4_front_vs_warming()
    with open("results/pcm_stage_b_audit.json", "w") as fh:
        json.dump(RESULTS, fh, indent=2, default=float)
    print("\nwrote results/pcm_stage_b_audit.json")
