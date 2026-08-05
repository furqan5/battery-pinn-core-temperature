"""Stage G diagnostic -- does the soft PDE constraint distort the recovered R_eff?

Seed 0 of Stage E fitted the surface to 0.196 K, better than the exact
finite-volume solver's best-possible 0.464 K.  The FV solution satisfies the PDE
by construction, so 0.464 K is the floor WITHIN the physics.  Beating it means
the PINN is spending PDE residual to buy data fit -- the constraint is soft.

If that is happening, the recovered R_eff should drift as the data weight rises.
Independent anchors for the true value:
    classical FV fit (physics exact)      14.30 mOhm
    electrical regression of V on I       14.49 mOhm

Cheaper settings than the headline run (this is about the TREND, not a final
number): fewer Adam epochs and L-BFGS steps, ~200 s per point.
"""
import time
import numpy as np
import torch

from part7_lib import P, Record, RadialFV
from stage_e_inverse import run_inverse, H_FIXED, K_FIXED

torch.set_num_threads(6)

rec = Record("2")
fv = RadialFV(N=40)
T0 = 0.5 * (rec.T_s[0] + rec.T_c[0])
grad = rec.T_c - rec.T_s

print("=" * 104)
print("STAGE G -- recovered R_eff vs data weight (soft-constraint distortion test)")
print("=" * 104)
print(f"  anchors:  classical FV fit 14.304 mOhm   |   "
      f"electrical regression {1000*rec.R_ohmic_reg:.3f} mOhm")
print(f"  physics floor on surface RMSE (exact FV, best R_eff) = 0.4637 K")
print()
print(f"  {'w_data':>8s} {'R_eff mOhm':>11s} {'surf RMSE':>10s} {'PDE+BC resid':>13s} "
      f"{'CORE RMSE':>10s} {'FV@sameR core':>14s} {'PDE violation K':>16s}")

for w_data in (5.0, 20.0, 200.0, 2000.0):
    t0 = time.time()
    r = run_inverse(rec, n_shape=0, seed=0, w_data=w_data,
                    adam_epochs=1000, lbfgs_steps=200)
    e = r["Tc"] - rec.T_c
    core = float(np.sqrt((e ** 2).mean()))
    o = fv.solve(rec.t, (rec.I ** 2) * r["R_eff"] / P.V_b, rec.T_inf, T0,
                 k=K_FIXED, h=H_FIXED, rho_cp=P.rho_cp())
    fv_core = float(np.sqrt(((o["T_core"] - rec.T_c) ** 2).mean()))
    viol = float(np.sqrt(((r["Ts"] - o["T_surf"]) ** 2).mean()))
    print(f"  {w_data:>8.0f} {1000*r['R_eff']:>11.4f} {r['surf_rmse']:>10.4f} "
          f"{r['sel']:>13.3e} {core:>10.4f} {fv_core:>14.4f} {viol:>16.4f}"
          f"   [{time.time()-t0:.0f} s]")

print()
print("  Reading: if R_eff drifts away from the ~14.3-14.5 mOhm anchors as w_data")
print("  rises, and surface RMSE drops below the 0.464 K physics floor, then the")
print("  network is buying data fit with PDE violation -- and the recovered")
print("  parameter is biased even though the fit looks excellent.")
