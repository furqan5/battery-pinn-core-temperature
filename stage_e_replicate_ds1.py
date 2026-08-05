"""Independent replication on DATASET 1, with the roles of the two records swapped.

This is the real test of whether the inverse PINN's advantage over the classical
control is CAPTURED PHYSICS or a fit to one record's idiosyncrasy.

Stage G's residual anatomy showed the classical model leaves a slow, strongly
autocorrelated (423 s, ~ the cell's own 405 s thermal time constant) drift whose
strongest correlate is TIME, not I^2.  A neural deviation field can absorb that
non-parametrically.  If that is all the PINN is doing, its advantage is a
property of DS2 and need not reappear on DS1.

Leak control, mirrored: h and k come from Stage B [B2] fitted to DS2's SURFACE
TRACE ONLY, and we fit DS1.  Nothing from DS1's core touches the fit.
"""
import sys
import time
import numpy as np
from scipy.optimize import minimize_scalar

from part7_lib import P, Record, RadialFV
import stage_e_inverse as E

# Stage B [B2] on DS2, SURFACE ONLY (clean initial condition)
H_DS2 = 37.2846
K_DS2 = 0.418697

n_seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 6
rec = Record("1")
fv = RadialFV(N=40)
T0 = float(rec.T_s[0])
grad = rec.T_c - rec.T_s

print("=" * 100)
print(f"REPLICATION -- inverse PINN on DATASET 1, surface only, {n_seeds} seeds")
print("=" * 100)
print(f"  h = {H_DS2:.4f} W/m2/K, k = {K_DS2:.6f} W/m/K  (Stage B [B2] on DS2 SURFACE ONLY)")
print(f"  Bi = {H_DS2*P.R_o/K_DS2:.4f}")
print(f"  measured core-surface max on DS1 = {grad.max():.4f} K")

# ---- classical control on the same footing, first ---- #
def pred(R, k=K_DS2, h=H_DS2):
    return fv.solve(rec.t, (rec.I ** 2) * R / P.V_b, rec.T_inf, T0,
                    k=k, h=h, rho_cp=P.rho_cp())


r = minimize_scalar(lambda R: float(np.sum((pred(R)["T_surf"] - rec.T_s) ** 2)),
                    bounds=(1e-4, 0.2), method="bounded")
o = pred(r.x)
cl_surf = float(np.sqrt(((o["T_surf"] - rec.T_s) ** 2).mean()))
cl_core = float(np.sqrt(((o["T_core"] - rec.T_c) ** 2).mean()))
print(f"\n  CLASSICAL control: R_eff {1000*r.x:.4f} mOhm, "
      f"surface {cl_surf:.4f} K, CORE {cl_core:.4f} K")
print(f"  electrical anchor on DS1: {1000*rec.R_ohmic_reg:.3f} mOhm")

runs = []
for sd in range(n_seeds):
    t0 = time.time()
    out = E.run_inverse(rec, h=H_DS2, k=K_DS2, n_shape=0, seed=sd)
    out["wall"] = time.time() - t0
    runs.append(out)
    print(f"    seed {sd}: {out['wall']:.0f} s, {out['closures']} closures, "
          f"R_eff {1000*out['R_eff']:.4f} mOhm, surface {out['surf_rmse']:.4f} K")

best, ok = E.summarise(rec, runs, "DS1 replication, order 0")

np.savez("results/stage_e_ds1_shape0.npz",
         Tc=np.array([x["Tc"] for x in runs]), Ts=np.array([x["Ts"] for x in runs]),
         R_eff=np.array([x["R_eff"] for x in runs]),
         sel=np.array([x["sel"] for x in runs]),
         core_rmse=np.array([x["core_rmse"] for x in runs]),
         surf_rmse=np.array([x["surf_rmse"] for x in runs]),
         T_c_meas=rec.T_c, T_s_meas=rec.T_s, t=rec.t)

print()
print("=" * 100)
print("  DOES THE PINN'S ADVANTAGE REPLICATE?")
print("=" * 100)
print(f"    {'':12s} {'classical':>12s} {'PINN':>12s} {'PINN/classical':>16s}")
print(f"    {'DS1 surface':12s} {cl_surf:>12.4f} {best['surf_rmse']:>12.4f} "
      f"{best['surf_rmse']/cl_surf:>16.2f}")
print(f"    {'DS1 core':12s} {cl_core:>12.4f} {best['core_rmse']:>12.4f} "
      f"{best['core_rmse']/cl_core:>16.2f}")
print(f"    (DS2 for comparison: classical 0.8943, PINN 0.6135, ratio 0.69)")
print()
print("    If the ratio here is also well below 1, the advantage replicates and is")
print("    a property of the METHOD.  If it is near or above 1, the DS2 advantage")
print("    was record-specific and should not be generalised.")
