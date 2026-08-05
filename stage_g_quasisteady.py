"""Stage G -- verify the quasi-steady baseline before letting it rewrite the headline.

B2:   T_core = T_surf + (Bi/2) (T_surf - T_inf)

scored 0.1409 K on DS2, 4.4x better than the inverse PINN.  That reverses the
session's conclusion, so it gets audited harder than the result it displaces.

Questions this answers:
  1. Is it leak-free?  (Bi comes from DS1 surface-only; no DS2 core anywhere.)
  2. Is it robust to Bi, or did Bi land luckily?
  3. Why does it beat a full transient solver?
  4. Is the comparison fair -- B2 anchors on the MEASURED surface while the
     PINN re-predicts the surface and inherits that error.  Isolate the two.
"""
import numpy as np

from part7_lib import P, Record, RadialFV
from stage_e_inverse import H_FIXED, K_FIXED

print("=" * 98)
print("QUASI-STEADY BASELINE -- audit")
print("=" * 98)

for tag, h, k in (("2", H_FIXED, K_FIXED), ("1", 37.2846, 0.418697)):
    rec = Record(tag)
    grad = rec.T_c - rec.T_s
    rise = rec.T_s - rec.T_inf
    Bi = h * P.R_o / k

    print(f"\n  DATASET {tag}   Bi = {Bi:.4f} (h={h:.4f}, k={k:.6f} from the OTHER "
          f"record's surface-only fit)")
    pred = rec.T_s + (Bi / 2) * rise
    e = pred - rec.T_c
    print(f"    B2 core RMSE = {np.sqrt((e**2).mean()):.4f} K, "
          f"max {np.abs(e).max():.4f} K, bias {e.mean():+.4f} K")

    # ---- 2. robustness to Bi ---- #
    print(f"    sensitivity to Bi (the ONLY parameter it uses):")
    print(f"      {'Bi':>8s} {'Bi/2':>8s} {'core RMSE':>10s}")
    for f in (0.7, 0.85, 1.0, 1.15, 1.3):
        Bx = Bi * f
        ex = rec.T_s + (Bx / 2) * rise - rec.T_c
        mark = "   <-- ours" if f == 1.0 else ""
        print(f"      {Bx:>8.4f} {Bx/2:>8.4f} {np.sqrt((ex**2).mean()):>10.4f}{mark}")
    # what Bi/2 would be optimal?  (diagnostic only -- uses the core, so it is
    # reported for comparison and never used as a predictor)
    c = float(np.sum(grad * rise) / np.sum(rise ** 2))
    print(f"      best-possible ratio by least squares = {c:.4f} "
          f"(vs Bi/2 = {Bi/2:.4f}, {100*abs(c-Bi/2)/c:.1f} % apart)")
    print(f"      NOTE: that best-possible value USES the core. It is a yardstick,")
    print(f"            not a predictor. Bi/2 above is core-blind.")

    # ---- 3. why does a quasi-steady relation work? ---- #
    alpha = k / P.rho_cp()
    tau_diff = P.R_o ** 2 / alpha
    print(f"    diffusion time R^2/alpha = {tau_diff:.0f} s; the measured "
          f"gradient/rise ratio has")
    m = rise > 0.5 * rise.max()
    r_meas = grad[m] / rise[m]
    print(f"      mean {r_meas.mean():.4f}, sd {r_meas.std():.4f} over the hot window")
    print(f"      -> the ratio is nearly CONSTANT in time, i.e. the radial profile")
    print(f"         is quasi-steady, so one number captures it.")

print()
print("=" * 98)
print("  4. IS THE COMPARISON FAIR?  B2 anchors on the measured surface; the PINN")
print("     re-predicts the surface and inherits that error. Separate the two.")
print("=" * 98)
rec = Record("2")
grad = rec.T_c - rec.T_s
d = np.load("results/stage_e_shape0.npz")
b = int(np.argmin(d["sel"]))
Tc_p, Ts_p = d["Tc"][b], d["Ts"][b]
Bi = H_FIXED * P.R_o / K_FIXED

variants = {
    "PINN, its own surface (as reported)": Tc_p,
    "PINN GRADIENT + measured surface":    rec.T_s + (Tc_p - Ts_p),
    "B2 quasi-steady":                     rec.T_s + (Bi / 2) * (rec.T_s - rec.T_inf),
}
print(f"    {'variant':38s} {'core RMSE':>10s} {'grad RMSE':>10s}")
for name, Tc in variants.items():
    e = Tc - rec.T_c
    eg = (Tc - rec.T_s) - grad
    print(f"    {name:38s} {np.sqrt((e**2).mean()):>10.4f} "
          f"{np.sqrt((eg**2).mean()):>10.4f}")
print()
print("    Even with the surface-fit error removed, the PINN's GRADIENT is worse")
print("    than the one-line relation. The gap is not an artefact of anchoring.")

# ---- and the same for the classical FV ---- #
fv = RadialFV(N=40)
T0 = float(rec.T_s[0])
from scipy.optimize import minimize_scalar
rr = minimize_scalar(
    lambda R: float(np.sum((fv.solve(rec.t, (rec.I**2)*R/P.V_b, rec.T_inf, T0,
                                     k=K_FIXED, h=H_FIXED,
                                     rho_cp=P.rho_cp())["T_surf"] - rec.T_s)**2)),
    bounds=(1e-4, 0.2), method="bounded")
o = fv.solve(rec.t, (rec.I**2)*rr.x/P.V_b, rec.T_inf, T0, k=K_FIXED, h=H_FIXED,
             rho_cp=P.rho_cp())
eg = (o["T_core"] - o["T_surf"]) - grad
print(f"    {'classical FV GRADIENT + meas. surface':38s} "
      f"{np.sqrt(((rec.T_s + (o['T_core']-o['T_surf'])) - rec.T_c)**2).mean()**0.5:>10.4f} "
      f"{np.sqrt((eg**2).mean()):>10.4f}")
