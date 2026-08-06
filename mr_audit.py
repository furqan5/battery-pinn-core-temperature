"""Audit the multi-rate timescale result before it is written up.

Three things could make the monotonic trend an artefact rather than physics:

  1. Railed or ill-conditioned fits. At 0.05C the surface rise is 0.13 K and
     R_eff came back at the bound. A railed fit is not an estimate.
  2. h varies 22-39 W/m2K across records, but h is a fixture property and should
     be constant. If the trend is really a Bi trend in disguise, it says nothing
     about timescale.
  3. The surface fit degrades with rate (0.23 -> 1.63 K). If the FV model does
     not describe the cell at 20C, its internal gradient is not evidence.

Audit 2 re-runs the comparison at a FIXED Bi, removing h from the trend
entirely. Audit 3 perturbs R_eff to show the quasi-steady error does not depend
on how well the source was fitted.
"""
import numpy as np
from scipy.optimize import minimize_scalar

from part7_lib import P, RadialFV
from mr_explore import load
from mr_timescale import segments, fit_h_from_tail, RHO_CP, K_T, TAU_DIFF

fv = RadialFV(N=40)
BI_FIXED = 1.2142          # Part 7, record 2, from a surface-only fit on record 1
H_FIXED = BI_FIXED * K_T / P.R_o


def prep(rate, f):
    """Uniform 1 Hz discharge segment plus its ambient and cooling tail."""
    dis, rest = segments(f)
    t = dis["t"].to_numpy() - dis["t"].to_numpy()[0]
    I = dis["I"].to_numpy()
    Ts = dis["Ts"].to_numpy()
    T_amb = (float(np.median(rest["Ts"].to_numpy()[-200:])) if len(rest) > 220
             else float(Ts[0]))
    n = int(np.floor(t[-1])) + 1
    tg = np.arange(n, dtype=float)
    Ig = np.interp(tg, t, I)
    Tsg = np.interp(tg, t, Ts)
    return dict(tg=tg, Ig=Ig, Tsg=Tsg, Tinf=np.full(n, T_amb), T_amb=T_amb,
                T0=float(Tsg[0]), dur=float(t[-1]), rest=rest)


def solve(s, h, R):
    return fv.solve(s["tg"], (s["Ig"] ** 2) * R / P.V_b, s["Tinf"], s["T0"],
                    k=K_T, h=h, rho_cp=RHO_CP)


def fit_R(s, h):
    r = minimize_scalar(
        lambda R: float(np.sum((solve(s, h, R)["T_surf"] - s["Tsg"]) ** 2)),
        bounds=(1e-5, 1.0), method="bounded", options={"xatol": 1e-12})
    return float(r.x)


def qs_error(s, h, R):
    """Relative error of the quasi-steady relation vs the transient solution."""
    o = solve(s, h, R)
    Bi = h * P.R_o / K_T
    grad = o["T_core"] - o["T_surf"]
    rise = o["T_surf"] - s["T_amb"]
    m = rise > 0.5 * rise.max()
    if m.sum() < 5:
        return np.nan
    e = (Bi / 2.0) * rise[m] - grad[m]
    return float(np.sqrt((e ** 2).mean()) / np.sqrt((grad[m] ** 2).mean()))


d = load()
S = {r: prep(r, d[r]) for r in sorted(d)}
H = {}
R = {}

print("=" * 100)
print("AUDIT 1 — fit quality per record")
print("=" * 100)
print(f"  {'C':>6s} {'rise K':>8s} {'tailRMSE':>9s} {'tau_cool':>9s} {'h':>7s} "
      f"{'R_eff mOhm':>11s} {'surfRMSE':>9s}  flag")
keep = []
for r_ in sorted(S):
    s = S[r_]
    h, tau_c, tr = fit_h_from_tail(s["rest"], s["T_amb"])
    Re = fit_R(s, h)
    H[r_], R[r_] = h, Re
    sr = float(np.sqrt(np.mean((solve(s, h, Re)["T_surf"] - s["Tsg"]) ** 2)))
    rise = float(s["Tsg"].max() - s["T_amb"])
    flags = []
    if Re > 0.99 or Re < 2e-5:
        flags.append("R_eff RAILED")
    if rise < 1.0:
        flags.append("no thermal signal")
    if sr > 0.05 * rise:
        flags.append("surface fit >5% of rise")
    if not flags:
        keep.append(r_)
    print(f"  {r_:>6.2f} {rise:>8.2f} {tr:>9.4f} {tau_c:>9.0f} {h:>7.2f} "
          f"{1000*Re:>11.2f} {sr:>9.4f}  {'; '.join(flags) if flags else 'ok'}")
print(f"\n  fully clean rates: {keep}")
print(f"  h over the high-signal records (>=5C): "
      f"{[round(H[k],2) for k in sorted(H) if k >= 5]}")

print()
print("=" * 100)
print("AUDIT 2 — the same trend at a FIXED Bi (h removed from the comparison)")
print("=" * 100)
print(f"  Bi = {BI_FIXED:.4f} for every record (Part 7 value), so nothing below")
print(f"  can be an artefact of the per-record h fit.")
print()
print(f"  {'C':>6s} {'t_dis s':>8s} {'t/tau':>7s} {'QS err (fitted h)':>18s} "
      f"{'QS err (fixed Bi)':>18s}")
rows = []
for r_ in sorted(S):
    s = S[r_]
    e_fit = qs_error(s, H[r_], R[r_])
    e_fix = qs_error(s, H_FIXED, R[r_])
    rows.append((r_, s["dur"], s["dur"] / TAU_DIFF, e_fit, e_fix))
    print(f"  {r_:>6.2f} {s['dur']:>8.0f} {s['dur']/TAU_DIFF:>7.2f} "
          f"{100*e_fit:>17.2f}% {100*e_fix:>17.2f}%")

print()
print("=" * 100)
print("AUDIT 3 — sensitivity of the quasi-steady error to the source amplitude")
print("=" * 100)
print(f"  {'C':>6s} {'@0.7 R_eff':>12s} {'@1.0 R_eff':>12s} {'@1.3 R_eff':>12s} "
      f"{'spread':>8s}")
for r_ in sorted(S):
    s = S[r_]
    es = [qs_error(s, H[r_], R[r_] * f) for f in (0.7, 1.0, 1.3)]
    print(f"  {r_:>6.2f} {100*es[0]:>11.2f}% {100*es[1]:>11.2f}% "
          f"{100*es[2]:>11.2f}% {100*(max(es)-min(es)):>7.3f}%")
print()
print("  Near-identical columns mean the quasi-steady error is set by the")
print("  timescale and the geometry, not by how well the source was fitted.")

np.savez("results/mr_audit.npz",
         rate=np.array([r[0] for r in rows]),
         dur=np.array([r[1] for r in rows]),
         tau_ratio=np.array([r[2] for r in rows]),
         qs_fitted_h=np.array([r[3] for r in rows]),
         qs_fixed_bi=np.array([r[4] for r in rows]),
         h=np.array([H[k] for k in sorted(H)]),
         R_eff=np.array([R[k] for k in sorted(R)]))
print("\n  saved results/mr_audit.npz")
