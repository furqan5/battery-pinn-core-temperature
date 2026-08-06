"""Stage C -- three-case discrimination, and Stage D -- time-resolved sensitivity.

Stage C asks the question the proposal actually turns on: are surface traces from
cells at 30 %, 50 % and 70 % melted distinguishable at 0.1 K noise?

Trap 6 says the total heat input must be held fixed or a different question is
being answered.  Here every case runs at the SAME q, for the SAME window
duration, through the SAME geometry, and differs only in how far along the melt
it starts.  The cumulative heat needed to REACH each state necessarily differs --
that is what melted fraction means -- but nothing in the compared window differs.

Stage D plots dT_surf/df0 through a full melt and past exhaustion, which is what
decides whether the honest claim is "we track the front" or "we detect
exhaustion".
"""
import json
import numpy as np

import pcm_params as P
import pcm_identify as I

RESULTS = {}
FRACTIONS = (0.30, 0.50, 0.70)
T_WINDOW = 600.0


def banner(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)


# --------------------------------------------------------------------------- #
# Stage C
# --------------------------------------------------------------------------- #

def stage_c(ref, t_window=T_WINDOW, sensors=("surface",), label=""):
    banner(f"STAGE C  Three-case discrimination at 30/50/70 % melted{label}")
    print(f"  identical q ({P.Q_NOMINAL:.0f} W/m3), identical {t_window:.0f} s window,")
    print(f"  identical geometry -- only the starting melted fraction differs")
    print(f"  noise sigma = {P.NOISE_SIGMA} K at 1 Hz")
    print()

    traces, t0s = {}, {}
    model = ref["model"]
    for f0 in FRACTIONS:
        T0, tstar = I.state_at_f(ref, f0)
        y, out = I.predict_window(model, T0, P.Q_NOMINAL, P.H_CONV,
                                  t_window, P.SAMPLE_DT, sensors)
        traces[f0] = y
        t0s[f0] = tstar
        idx = I.sensor_indices(model)
        print(f"  f = {f0:.2f}  starts at t = {tstar:7.1f} s, "
              f"T_surf {out['T'][0, idx['surface']]:8.4f} -> "
              f"{out['T'][-1, idx['surface']]:8.4f} K, "
              f"ends at f = {out['melted_fraction'][-1]:.4f}")

    print()
    print(f"  {'pair':>14s} {'max |dT| (K)':>14s} {'ratio to noise':>16s}  verdict")
    rows = []
    pairs = [(0.30, 0.50), (0.50, 0.70), (0.30, 0.70)]
    for a, b in pairs:
        d = float(np.max(np.abs(traces[a] - traces[b])))
        ratio = d / P.NOISE_SIGMA
        verdict = ("distinguishable" if ratio > 3 else
                   "marginal" if ratio > 1 else "BELOW NOISE")
        rows.append({"pair": [a, b], "max_dT": d, "ratio": ratio,
                     "verdict": verdict})
        print(f"  {f'{a:.2f} vs {b:.2f}':>14s} {d:14.4f} {ratio:16.2f}  {verdict}")

    worst = min(r["max_dT"] for r in rows)
    print()
    print(f"  smallest pairwise separation = {worst:.4f} K "
          f"({worst/P.NOISE_SIGMA:.2f} x noise)")
    return {"rows": rows, "t_starts": t0s, "worst_separation": worst,
            "window": t_window, "sensors": list(sensors)}


def stage_c_mechanism(ref):
    """Is the separation caused by the FRONT, or by the cell warming?

    Same three-case test on a cell whose PCM has almost no latent heat.  If the
    separation survives, it was never about the front.
    """
    banner("STAGE C control -- separation with the latent heat removed")
    mdl0 = P.build_model(Lf=P.PCM_LF / 1000.0)
    ref0 = I.reference_trajectory(model=mdl0, q=P.Q_NOMINAL, t_end=9000.0)
    r = stage_c(ref0, label="  [L_f / 1000, no real front]")
    return r


# --------------------------------------------------------------------------- #
# Stage D
# --------------------------------------------------------------------------- #

def stage_d(ref, f0=0.30, df=5e-3, t_end=9000.0):
    """dT_surf/df0 as a function of time, from f0 through exhaustion and beyond.

    Two trajectories differing only by df in initial melted fraction, propagated
    together.  During melting they are near time-shifted copies and differ by a
    roughly constant offset.  Around exhaustion one cell runs out of PCM before
    the other, so the traces separate sharply.  Where the sensitivity lives is
    what the method can honestly claim.
    """
    banner("STAGE D  Time-resolved dT_surf/d(melted fraction), through exhaustion")
    model = ref["model"]
    idx = I.sensor_indices(model)

    runs = {}
    for tag, ff in (("plus", f0 + df), ("minus", f0 - df)):
        T0, tstar = I.state_at_f(ref, ff)
        out = model.solve(T0, t_end, P.SAMPLE_DT, q_vol=P.Q_NOMINAL,
                          q_layer="core", bc_inner=("adiabatic",),
                          bc_outer=("robin", P.H_CONV, P.T_AMB),
                          scheme="enthalpy_newton", store_every=1)
        runs[tag] = out

    t = runs["plus"]["t"]
    sens = (runs["plus"]["T"][:, idx["surface"]]
            - runs["minus"]["T"][:, idx["surface"]]) / (2.0 * df)
    fbar = 0.5 * (runs["plus"]["melted_fraction"]
                  + runs["minus"]["melted_fraction"])

    # exhaustion time: when the mean melted fraction first reaches 1
    ex = int(np.searchsorted(fbar, 1.0 - 1e-9))
    t_ex = float(t[min(ex, len(t) - 1)])

    during = (fbar > 0.2) & (fbar < 0.8)
    after = (t > t_ex) & (t <= t_ex + 200.0)
    mean_during = float(np.mean(np.abs(sens[during]))) if during.any() else np.nan
    peak_after = float(np.max(np.abs(sens[after]))) if after.any() else np.nan
    ratio = peak_after / mean_during if mean_during else np.inf

    print(f"  started from f0 = {f0:.2f} +/- {df}")
    print(f"  exhaustion (f -> 1) at t = {t_ex:.0f} s after window start")
    print(f"  mean |dT/df0| during melting (0.2 < f < 0.8)  {mean_during:10.4f} K")
    print(f"  peak |dT/df0| in 200 s after exhaustion       {peak_after:10.4f} K")
    print(f"  ratio after / during                          {ratio:10.4f}")
    print()
    print(f"  {'t (s)':>8s} {'f':>8s} {'dT_surf/df0 (K)':>18s}")
    for tt in (0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 5000, 6000):
        j = int(np.searchsorted(t, tt))
        if j < len(t):
            print(f"  {t[j]:8.0f} {fbar[j]:8.4f} {sens[j]:18.4f}")

    RESULTS["D"] = {
        "f0": f0, "df": df, "t_exhaustion": t_ex,
        "mean_during": mean_during, "peak_after": peak_after, "ratio": ratio,
        "t": t.tolist(), "sens": sens.tolist(), "f": fbar.tolist(),
    }
    return t, sens, fbar, t_ex


if __name__ == "__main__":
    print("building reference trajectory ...")
    ref = I.reference_trajectory(t_end=14000.0)

    RESULTS["C_surface"] = stage_c(ref, sensors=("surface",))
    RESULTS["C_two_sensor"] = stage_c(ref, sensors=("surface", "can"),
                                      label="  [surface + can]")
    RESULTS["C_control_no_latent"] = stage_c_mechanism(ref)
    stage_d(ref)

    with open("results/pcm_stage_cd.json", "w") as fh:
        json.dump(RESULTS, fh, indent=2, default=float)
    print("\nwrote results/pcm_stage_cd.json")
