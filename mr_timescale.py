"""Multi-rate timescale test: where does the quasi-steady relation break down?

Chapter 6 established that on a drive cycle the radial profile shape is constant
to ~2 %, so one algebraic relation captures the internal field and a transient
solver adds nothing.  Chapter 7 wants to generalise that into a criterion.  This
script tests it across two and a half decades of forcing timescale on the SAME
CELL FAMILY (A123 ANR26650m1-b, per the dataset's own manufacturer sheet).

Design, and its honest limits:

  * h is fitted PER RECORD from its own 3600 s cooling tail, with the current
    off, so no heat-source model is involved and no fixture assumption is
    imported from Part 7.  (This is the "exploit the relaxation tail" route
    recommended in FINDINGS.)
  * R_eff is then fitted to the discharge surface rise, and the model is scored
    against the MEASURED surface at every rate.  That is a real validation.
  * The internal gradient is MODELLED, not measured — this dataset has no core
    sensor.  The model is the one Part 7 validated against a real core
    thermocouple on the same cell family, which is the strongest available
    warrant, but it is not an internal measurement and is not claimed as one.

Caveats carried into the write-up:
  * Cell variant differs: 2.50 Ah / 76 g here against the 2.3 Ah parameter set
    Richardson used.  rho*cp is taken from Part 7 and is (b).
  * Ambient is not logged; the nominal chamber temperature is used, cross-checked
    against the pre-discharge rest temperature.
"""
import os
import numpy as np
from scipy.optimize import curve_fit, minimize_scalar

from part7_lib import P, RadialFV
from mr_explore import load

fv = RadialFV(N=40)
RHO_CP = P.rho_cp()          # (b) from Part 7
K_T = P.k_t                  # (b) same cell family
TAU_DIFF = P.R_o ** 2 / (K_T / RHO_CP)

DISCHARGE_STEP = 5
REST_STEP = 6


def segments(f):
    """Discharge segment and its trailing cooling rest."""
    d = f[f["step"] == DISCHARGE_STEP]
    r = f[f["step"] == REST_STEP]
    return d.reset_index(drop=True), r.reset_index(drop=True)


def fit_h_from_tail(rest, T_amb):
    """h from the cooling decay alone. No source model, so no confounding.

    Lumped cooling: T - T_amb = dT0 exp(-t/tau), tau = rho_cp * R / (2h).
    """
    t = rest["t"].to_numpy() - rest["t"].to_numpy()[0]
    y = rest["Ts"].to_numpy() - T_amb
    if y[0] <= 0.3:
        return np.nan, np.nan, np.nan
    m = y > 0.05 * y[0]                       # stop before it buries in noise
    t, y = t[m], y[m]

    def mod(tt, dT0, tau, off):
        return dT0 * np.exp(-tt / tau) + off

    try:
        p, _ = curve_fit(mod, t, y, p0=(y[0], 600.0, 0.0),
                         bounds=([0, 30, -2], [50, 20000, 2]), maxfev=20000)
    except Exception:
        return np.nan, np.nan, np.nan
    tau = p[1]
    h = RHO_CP * P.R_o / (2.0 * tau)
    rmse = float(np.sqrt(np.mean((mod(t, *p) - y) ** 2)))
    return h, tau, rmse


def run_rate(rate, f, verbose=True):
    dis, rest = segments(f)
    if len(dis) < 20:
        return None
    t = dis["t"].to_numpy() - dis["t"].to_numpy()[0]
    I = dis["I"].to_numpy()
    Ts = dis["Ts"].to_numpy()
    dur = float(t[-1])

    # ambient: the plateau the cell relaxes to at the end of the cooling tail
    T_amb = float(np.median(rest["Ts"].to_numpy()[-200:])) if len(rest) > 220 \
        else float(Ts[0])

    h, tau_c, tail_rmse = fit_h_from_tail(rest, T_amb)
    if not np.isfinite(h):
        return None
    Bi = h * P.R_o / K_T

    # uniform 1 Hz grid for the solver
    n = int(np.floor(dur)) + 1
    tg = np.arange(n, dtype=float)
    Ig = np.interp(tg, t, I)
    Tsg = np.interp(tg, t, Ts)
    Tinf = np.full(n, T_amb)
    T0 = float(Tsg[0])

    def predict(R):
        return fv.solve(tg, (Ig ** 2) * R / P.V_b, Tinf, T0,
                        k=K_T, h=h, rho_cp=RHO_CP)

    r = minimize_scalar(lambda R: float(np.sum((predict(R)["T_surf"] - Tsg) ** 2)),
                        bounds=(1e-5, 1.0), method="bounded",
                        options={"xatol": 1e-12})
    o = predict(r.x)
    surf_rmse = float(np.sqrt(np.mean((o["T_surf"] - Tsg) ** 2)))

    # the quantity under test: transient gradient vs the quasi-steady relation
    grad_fv = o["T_core"] - o["T_surf"]
    rise = o["T_surf"] - T_amb
    grad_qs = (Bi / 2.0) * rise

    # score over the thermally active part only
    m = rise > 0.5 * rise.max()
    err = grad_qs[m] - grad_fv[m]
    rel = float(np.sqrt(np.mean(err ** 2)) / np.sqrt(np.mean(grad_fv[m] ** 2)))
    ratio_mean = float(np.mean((grad_fv / np.maximum(rise, 1e-9))[m]))

    return dict(rate=rate, dur=dur, tau_ratio=dur / TAU_DIFF, h=h, Bi=Bi,
                tau_cool=tau_c, tail_rmse=tail_rmse, R_eff=r.x,
                surf_rmse=surf_rmse, T_amb=T_amb,
                rise_max=float(rise.max()), grad_max=float(grad_fv.max()),
                qs_rel_err=rel, ratio_mean=ratio_mean, bi_over_2=Bi / 2.0,
                I_rms=float(np.sqrt((Ig ** 2).mean())))


if __name__ == "__main__":
    print("=" * 108)
    print("MULTI-RATE TIMESCALE TEST — A123 26650 LFP, 25 degC chamber")
    print("=" * 108)
    print(f"  tau_diff = R^2/alpha = {TAU_DIFF:.0f} s   "
          f"(k = {K_T} W/m/K, rho*cp = {RHO_CP:.4g} J/m^3/K, both (b) from Part 7)")
    print(f"  h fitted per record from its own cooling tail; no Part 7 fixture assumed.")
    print()

    d = load()
    rows = []
    for rate in sorted(d):
        out = run_rate(rate, d[rate])
        if out:
            rows.append(out)

    print(f"  {'C':>6s} {'t_dis':>8s} {'t/tau':>7s} {'h':>7s} {'Bi':>6s} "
          f"{'tau_cool':>8s} {'R_eff':>8s} {'surfRMSE':>9s} {'rise':>6s} "
          f"{'grad':>6s} {'QS err':>8s}")
    print(f"  {'':>6s} {'s':>8s} {'':>7s} {'W/m2K':>7s} {'':>6s} "
          f"{'s':>8s} {'mOhm':>8s} {'K':>9s} {'K':>6s} {'K':>6s} {'%':>8s}")
    print("  " + "-" * 100)
    for r in rows:
        print(f"  {r['rate']:>6.2f} {r['dur']:>8.0f} {r['tau_ratio']:>7.2f} "
              f"{r['h']:>7.2f} {r['Bi']:>6.3f} {r['tau_cool']:>8.0f} "
              f"{1000*r['R_eff']:>8.2f} {r['surf_rmse']:>9.4f} "
              f"{r['rise_max']:>6.2f} {r['grad_max']:>6.2f} "
              f"{100*r['qs_rel_err']:>8.2f}")

    np.savez("results/mr_timescale.npz",
             **{k: np.array([r[k] for r in rows]) for k in rows[0]})
    print()
    print("  Reading: QS err is the error of the quasi-steady relation against the")
    print("  transient solution, over the thermally active window. If the criterion")
    print("  is right it should be small when t_dis >> tau_diff and grow as the")
    print("  ratio falls through 1.")
    print()
    print("  saved results/mr_timescale.npz")
