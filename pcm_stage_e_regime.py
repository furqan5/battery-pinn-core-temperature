"""Stage E -- the regime sweep.  The most transferable output of the session.

Vary the heating rate so the melting time spans more than an order of magnitude
around the PCM layer's diffusion time, and ask where -- if anywhere -- the melt
FRONT becomes identifiable during melting.  This is the direct analogue of the
timescale criterion that explained the battery result in Part 7.

The metric that matters is not raw separation.  Stage B/L4 and Stage C showed
that a cell with essentially no latent heat gives the SAME separation, so raw
separation measures monotone warming, not the front.  The front-specific metric
used here is the excess:

    front signal = (separation with real PCM) - (separation with L_f / 1000)

run at the same heating rate, same window fraction, same geometry.  If that
excess stays inside the noise at every heating rate, the front never becomes
observable at any rate, and the answer to the session's question is settled.

Two robustness sweeps follow, covering the two parameters flagged in
pcm_params.py as not fully sourced: the solid/liquid conductivity ratio, and the
mushy-zone width (trap 5).
"""
import json
import numpy as np

import pcm_params as P
import pcm_identify as I

RESULTS = {}
FRACTIONS = (0.30, 0.50, 0.70)
WINDOW_FRACTION = 0.10     # window as a fraction of the melt duration


def banner(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)


def melt_span(ref):
    f = ref["melted_fraction"]
    if f[-1] < 1.0 - 1e-9:
        return None
    i0 = int(np.searchsorted(f, 1e-9))
    i1 = int(np.searchsorted(f, 1.0 - 1e-9))
    return float(ref["t"][i0]), float(ref["t"][i1])


def build_ref(q, model=None, margin=2.0):
    """Reference trajectory with an end time adapted to the heating rate."""
    if model is None:
        model = P.build_model()
    ts = P.timescales(model, q=q)
    if not np.isfinite(ts["t_melt"]) or ts["Q_net"] <= 0:
        return None, ts
    t_end = float(margin * (ts["t_melt"] + 2000.0))
    t_end = min(max(t_end, 1500.0), 60000.0)
    ref = I.reference_trajectory(model=model, q=q, t_end=t_end)
    return ref, ts


def separation(ref, t_window, sensors=("surface",)):
    """Worst pairwise max|dT| across the 30/50/70 % cases."""
    model = ref["model"]
    tr = {}
    for f0 in FRACTIONS:
        T0, _ = I.state_at_f(ref, f0)
        y, _ = I.predict_window(model, T0, ref["q"], ref["h"], t_window,
                                P.SAMPLE_DT, sensors)
        tr[f0] = y
    pairs = [(0.30, 0.50), (0.50, 0.70), (0.30, 0.70)]
    seps = [float(np.max(np.abs(tr[a] - tr[b]))) for a, b in pairs]
    return min(seps), seps


def s1_bound(ref, t_window):
    base = {"f0": 0.50, "q": ref["q"], "h": ref["h"]}
    ev = lambda th: I._eval_window(th, ref, t_window, ("surface",))
    S, _, _ = I.sensitivity_matrix(base, ["f0"], ev)
    r = I.crlb_from_S(S, ["f0"], base)
    return 100.0 * r["sd_natural"][0], float(np.mean(S[:, 0]))


# --------------------------------------------------------------------------- #

def e1_regime_sweep():
    banner("E1  Heating-rate sweep: does the FRONT ever become observable?")
    print(f"  window = {WINDOW_FRACTION:.0%} of each melt duration, so every regime")
    print("  observes a comparable slice of its own melt rather than a fixed")
    print("  wall-clock time (which would confound window coverage with rate).")
    print()
    print(f"  {'q/q_3C':>7s} {'C-rate':>7s} {'t_melt':>9s} {'t_m/tau':>9s} "
          f"{'S1 sd':>8s} {'sep PCM':>9s} {'sep noLf':>9s} {'front':>9s} {'front/sig':>10s}")
    print(f"  {'':>7s} {'equiv':>7s} {'(s)':>9s} {'':>9s} "
          f"{'(%buf)':>8s} {'(K)':>9s} {'(K)':>9s} {'excess(K)':>9s} {'':>10s}")

    rows = []
    for mult in (1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 14.0, 24.0, 40.0):
        q = P.Q_NOMINAL * mult
        ref, ts = build_ref(q)
        if ref is None:
            print(f"  {mult:7.1f}  never melts (heat input below ambient loss)")
            continue
        span = melt_span(ref)
        if span is None:
            print(f"  {mult:7.1f}  did not reach full melt within the horizon")
            continue
        dur = span[1] - span[0]
        tw = max(30.0, round(WINDOW_FRACTION * dur))

        sep_pcm, _ = separation(ref, tw)
        sd1, _ = s1_bound(ref, tw)

        # Control at the same heating rate: no latent heat, hence no front.
        # Its window is set from ITS OWN melt duration, not the real PCM's.
        # Using the same wall-clock window for both would hand the control a
        # window covering ~40% of its (much faster) melt against 10% of the real
        # one, and the extra coverage alone would inflate its separation.
        mdl0 = P.build_model(Lf=P.PCM_LF / 1000.0)
        ref0, _ = build_ref(q, model=mdl0)
        sep0, tw0 = np.nan, np.nan
        if ref0 is not None:
            sp0 = melt_span(ref0)
            if sp0 is not None:
                tw0 = max(30.0, round(WINDOW_FRACTION * (sp0[1] - sp0[0])))
                sep0, _ = separation(ref0, tw0)

        excess = sep_pcm - sep0
        crate = np.sqrt(q * P.CELL_VOLUME / P.CELL_R_1KHZ) / P.CELL_CAPACITY_AH
        rows.append({
            "mult": mult, "q": q, "c_rate": float(crate),
            "t_melt": dur, "t_melt_over_tau": dur / ts["tau_diff"],
            "window": tw, "window_control": float(tw0), "S1_sd_buffer_pct": sd1,
            "sep_pcm": sep_pcm, "sep_no_latent": float(sep0),
            "front_excess": float(excess),
            "front_excess_over_noise": float(excess / P.NOISE_SIGMA),
        })
        print(f"  {mult:7.1f} {crate:7.1f} {dur:9.0f} {dur/ts['tau_diff']:9.2f} "
              f"{sd1:8.3f} {sep_pcm:9.4f} {sep0:9.4f} {excess:9.4f} "
              f"{excess/P.NOISE_SIGMA:10.2f}")

    print()
    print("  'front excess' is the separation attributable to phase change:")
    print("  what the real PCM gives MINUS what a latent-free cell gives at the")
    print("  same heating rate.  'front/sig' expresses it in units of the 0.1 K")
    print("  sensor noise.  Values below ~3 mean the front contributes nothing")
    print(f"  an estimator could use, whatever the raw separation looks like.")
    print()
    print(f"  NOTE: C-rate equivalents above {P.CELL_C_RATE_MAX:.0f}C exceed the")
    print("  manufacturer rating for this cell; they are included as regime")
    print("  probes, not as operating proposals.")
    RESULTS["E1"] = rows
    return rows


def e2_k_ratio():
    banner("E2  Robustness: solid/liquid conductivity ratio")
    print("  [S2] gives a single conductivity, so k_l/k_s is the least-sourced")
    print("  quantity in the study -- and it is the mechanism by which front")
    print("  POSITION changes the composite resistance.  Swept, not assumed.")
    print()
    print(f"  {'k_l/k_s':>9s} {'S1 sd (%buf)':>14s} {'sep PCM (K)':>13s} "
          f"{'dT/df0 (K)':>12s}")
    rows = []
    for kr in P.K_RATIO_SWEEP:
        mdl = P.build_model(k_ratio=kr)
        ref, _ = build_ref(P.Q_NOMINAL, model=mdl)
        if ref is None or melt_span(ref) is None:
            print(f"  {kr:9.2f}  did not melt")
            continue
        sd1, sens = s1_bound(ref, 600.0)
        sep, _ = separation(ref, 600.0)
        rows.append({"k_ratio": kr, "S1_sd_buffer_pct": sd1,
                     "sep": sep, "mean_sens": sens})
        print(f"  {kr:9.2f} {sd1:14.4f} {sep:13.4f} {sens:12.4f}")
    RESULTS["E2"] = rows
    return rows


def e3_mushy_width():
    banner("E3  Robustness: mushy-zone width (trap 5)")
    print("  For RT42 the 4 K width is largely PHYSICAL -- the datasheet melting")
    print("  range is 38-43 C -- but the apparent-heat-capacity dTm is also a")
    print("  numerical parameter, so the conclusions must be stable across it.")
    print()
    print(f"  {'dTm (K)':>9s} {'S1 sd (%buf)':>14s} {'sep PCM (K)':>13s} "
          f"{'dT/df0 (K)':>12s} {'melt span (s)':>15s}")
    rows = []
    for dTm in P.DTM_SWEEP:
        mdl = P.build_model(dTm=dTm)
        ref, _ = build_ref(P.Q_NOMINAL, model=mdl)
        sp = melt_span(ref) if ref is not None else None
        if sp is None:
            print(f"  {dTm:9.2f}  did not melt")
            continue
        sd1, sens = s1_bound(ref, 600.0)
        sep, _ = separation(ref, 600.0)
        rows.append({"dTm": dTm, "S1_sd_buffer_pct": sd1, "sep": sep,
                     "mean_sens": sens, "melt_span": sp[1] - sp[0]})
        print(f"  {dTm:9.2f} {sd1:14.4f} {sep:13.4f} {sens:12.4f} "
              f"{sp[1]-sp[0]:15.0f}")
    RESULTS["E3"] = rows
    return rows


if __name__ == "__main__":
    e1_regime_sweep()
    e2_k_ratio()
    e3_mushy_width()
    with open("results/pcm_stage_e.json", "w") as fh:
        json.dump(RESULTS, fh, indent=2, default=float)
    print("\nwrote results/pcm_stage_e.json")
