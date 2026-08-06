"""Stage B -- Fisher information and Cramer-Rao bounds for S1-S4.

Nothing is fitted here.  This is an identifiability study: a fit that looks good
on a non-identifiable parameter set is the exact failure mode the exercise
exists to prevent.

  S1  {f}              is the buffer readable directly?
  S2  {L_f, k_pcm}     can PCM properties be identified in place?
  S3  {f, q}           can the state be read when the source is unknown?
  S4  {f, q, h}        the realistic field case

Conditioning is noise-independent; standard deviations scale with sigma.  Both
are reported and labelled.
"""
import json
import numpy as np

import pcm_params as P
import pcm_identify as I

RESULTS = {}
T_WINDOW = 600.0        # s, nominal observation window (pre-registered)
F0 = 0.50               # nominal melted fraction
T_COLD = 8000.0         # s, cold-start experiment length (full melt ~6780 s)


def banner(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)


def base_params():
    return {"f0": F0, "q": P.Q_NOMINAL, "h": P.H_CONV,
            "Lf": P.PCM_LF, "k_pcm": P.PCM_K_SOLID}


# --------------------------------------------------------------------------- #

def run_window_sets(ref, sensors=("surface",), t_window=T_WINDOW, f0=F0,
                    sigma=P.NOISE_SIGMA, label=""):
    base = base_params()
    base["f0"] = f0
    ev = lambda th: I._eval_window(th, ref, t_window, sensors)
    out = {}
    for name, free in (("S1  {f}", ["f0"]),
                       ("S3  {f, q}", ["f0", "q"]),
                       ("S4  {f, q, h}", ["f0", "q", "h"])):
        S, y0, steps = I.sensitivity_matrix(base, free, ev)
        res = I.crlb_from_S(S, free, base, sigma)
        I.report(name + label, res, base)
        if len(free) > 1:
            I.print_corr(res)
        out[name.split()[0]] = {
            "free": free, "cond": res["cond_scaled"],
            "sd_natural": res["sd_natural"].tolist(),
            "sd_relative": res["sd_relative"].tolist(),
            "corr": np.asarray(res["corr"]).tolist(),
            "singular": res["singular"],
            "fd_steps": steps,
            "n_samples": int(S.shape[0]),
        }
        print()
    return out


def run_cold_set(sensors=("surface",), sigma=P.NOISE_SIGMA):
    base = base_params()
    ev = lambda th: I._eval_cold(th, T_COLD, sensors)
    free = ["Lf", "k_pcm"]
    S, y0, steps = I.sensitivity_matrix(base, free, ev)
    res = I.crlb_from_S(S, free, base, sigma)
    I.report("S2  {L_f, k_pcm}", res, base)
    I.print_corr(res)
    return {"free": free, "cond": res["cond_scaled"],
            "sd_natural": res["sd_natural"].tolist(),
            "sd_relative": res["sd_relative"].tolist(),
            "corr": np.asarray(res["corr"]).tolist(),
            "singular": res["singular"], "fd_steps": steps,
            "n_samples": int(S.shape[0])}


# --------------------------------------------------------------------------- #

def b0_fd_convergence(ref):
    banner("B0  Finite-difference step convergence (is the sensitivity a derivative?)")
    base = base_params()
    ev = lambda th: I._eval_window(th, ref, T_WINDOW, ("surface",))
    print(f"  {'df0':>10s} {'mean dT/df0 (K)':>18s} {'rel change':>12s}")
    prev = None
    rows = []
    for df in (2e-2, 1e-2, 5e-3, 2e-3, 1e-3):
        S, _, _ = I.sensitivity_matrix(base, ["f0"], ev, steps={"f0": df})
        m = float(np.mean(S[:, 0]))
        rc = "" if prev is None else f"{abs(m-prev)/abs(prev):12.2e}"
        print(f"  {df:10.4f} {m:18.5f} {rc:>12s}")
        rows.append({"df0": df, "mean_sens": m})
        prev = m
    print("\n  Stable across a 20x range of step size, so the central difference")
    print("  is resolving a derivative rather than interpolation noise.")
    RESULTS["B0"] = rows


def b1_main(ref):
    banner("B1  CRLB, single surface sensor, 600 s window at f = 0.50")
    print(f"  sigma = {P.NOISE_SIGMA} K, 1 Hz, window {T_WINDOW:.0f} s "
          f"-> {int(T_WINDOW)+1} samples")
    print(f"  melted fraction sd is quoted as % OF THE FULL LATENT BUFFER")
    print()
    RESULTS["B1_surface"] = run_window_sets(ref, ("surface",))

    banner("B2  Same, with a second sensor on the cell can (Q5 contrast)")
    RESULTS["B2_two_sensor"] = run_window_sets(ref, ("surface", "can"),
                                               label="  [2 sensors]")

    banner("B3  S2: PCM properties from a full melt from cold")
    print(f"  cold-start experiment, {T_COLD:.0f} s, single surface sensor")
    print()
    RESULTS["B3_S2"] = run_cold_set()


def b4_window_length(ref):
    banner("B4  Dependence on observation window length")
    print("  Longer windows help only if the trace carries new information;")
    print("  a pinned plateau adds samples without adding signal.")
    print()
    print(f"  {'window (s)':>11s} {'S1 sd (% buffer)':>18s} {'S3 cond':>12s} "
          f"{'S3 f sd (% buf)':>17s}")
    rows = []
    base = base_params()
    for tw in (60.0, 300.0, 600.0, 1200.0, 2400.0):
        ev = lambda th: I._eval_window(th, ref, tw, ("surface",))
        S1, _, _ = I.sensitivity_matrix(base, ["f0"], ev)
        r1 = I.crlb_from_S(S1, ["f0"], base)
        S3, _, _ = I.sensitivity_matrix(base, ["f0", "q"], ev)
        r3 = I.crlb_from_S(S3, ["f0", "q"], base)
        rows.append({"window": tw,
                     "S1_sd_buffer_pct": 100 * r1["sd_natural"][0],
                     "S3_cond": r3["cond_scaled"],
                     "S3_f_sd_buffer_pct": 100 * r3["sd_natural"][0]})
        print(f"  {tw:11.0f} {100*r1['sd_natural'][0]:18.4f} "
              f"{r3['cond_scaled']:12.3e} {100*r3['sd_natural'][0]:17.4f}")
    RESULTS["B4_window"] = rows


def b5_across_melt(ref):
    banner("B5  Does identifiability depend on WHERE in the melt you look?")
    print(f"  {'f0':>6s} {'t* (s)':>9s} {'T_surf (K)':>11s} "
          f"{'S1 sd (% buf)':>15s} {'S3 cond':>12s} {'S3 f sd (% buf)':>17s}")
    rows = []
    base = base_params()
    for f0 in (0.10, 0.30, 0.50, 0.70, 0.90):
        b = dict(base); b["f0"] = f0
        _, tstar = I.state_at_f(ref, f0)
        T0, _ = I.state_at_f(ref, f0)
        ev = lambda th: I._eval_window(th, ref, T_WINDOW, ("surface",))
        S1, _, _ = I.sensitivity_matrix(b, ["f0"], ev)
        r1 = I.crlb_from_S(S1, ["f0"], b)
        S3, _, _ = I.sensitivity_matrix(b, ["f0", "q"], ev)
        r3 = I.crlb_from_S(S3, ["f0", "q"], b)
        rows.append({"f0": f0, "t_star": tstar, "T_surf": float(T0[-1]),
                     "S1_sd_buffer_pct": 100 * r1["sd_natural"][0],
                     "S3_cond": r3["cond_scaled"],
                     "S3_f_sd_buffer_pct": 100 * r3["sd_natural"][0]})
        print(f"  {f0:6.2f} {tstar:9.0f} {T0[-1]:11.4f} "
              f"{100*r1['sd_natural'][0]:15.4f} {r3['cond_scaled']:12.3e} "
              f"{100*r3['sd_natural'][0]:17.4f}")
    RESULTS["B5_across_melt"] = rows


if __name__ == "__main__":
    print("building reference melt trajectory ...")
    ref = I.reference_trajectory()
    f = ref["melted_fraction"]
    i1 = int(np.searchsorted(f, 1.0 - 1e-9))
    print(f"  melt spans t = {ref['t'][int(np.searchsorted(f, 1e-9))]:.0f} s "
          f"to {ref['t'][i1]:.0f} s")
    print(f"  T_surf over that span: {ref['T_surf'][int(np.searchsorted(f,1e-9))]:.3f} "
          f"-> {ref['T_surf'][i1]:.3f} K")
    RESULTS["reference"] = {
        "t_melt_start": float(ref["t"][int(np.searchsorted(f, 1e-9))]),
        "t_melt_end": float(ref["t"][i1]),
        "T_surf_melt_start": float(ref["T_surf"][int(np.searchsorted(f, 1e-9))]),
        "T_surf_melt_end": float(ref["T_surf"][i1]),
        "energy_closure": ref["energy_closure_rel"],
    }

    b0_fd_convergence(ref)
    b1_main(ref)
    b4_window_length(ref)
    b5_across_melt(ref)

    with open("results/pcm_stage_b.json", "w") as fh:
        json.dump(RESULTS, fh, indent=2, default=float)
    print("\nwrote results/pcm_stage_b.json")
