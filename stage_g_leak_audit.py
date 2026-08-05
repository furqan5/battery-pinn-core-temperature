"""Leak audit -- does ANY core-channel information reach the fit?

Triggered by the rule: if a result looks too good, check for a leak before
reporting it.  Order 0 came in at 0.6158 K core RMSE, better than the classical
control's 0.8895 K, so this ran before the number was reported.

FOUND ONE.  The initial condition was set as

    T0 = 0.5 * (T_s[0] + T_c[0])

which reads the measured CORE at t=0.  It is a single number, and the cell is
essentially isothermal at t=0 (core-surface = 0.038 K on DS2), so the leak is
small -- but it is real, and "small" is not "absent".

This script quantifies it and provides the clean replacement, T0 = T_s[0].
"""
import numpy as np
from scipy.optimize import minimize_scalar

from part7_lib import P, Record, RadialFV
from stage_e_inverse import H_FIXED, K_FIXED

fv = RadialFV(N=40)

print("=" * 96)
print("LEAK AUDIT")
print("=" * 96)

for tag in ("1", "2"):
    rec = Record(tag)
    leaky = 0.5 * (rec.T_s[0] + rec.T_c[0])
    clean = rec.T_s[0]
    print(f"\n  DATASET {tag}")
    print(f"    T_s(0) = {rec.T_s[0]:.4f} C, T_c(0) = {rec.T_c[0]:.4f} C, "
          f"core-surf at t=0 = {rec.T_c[0]-rec.T_s[0]:+.4f} K")
    print(f"    leaky T0 = {leaky:.4f} C   clean T0 = {clean:.4f} C   "
          f"difference {leaky-clean:+.4f} K")

    for lab, T0 in (("leaky (uses T_c[0])", leaky), ("clean (surface only)", clean)):
        def pred(R):
            return fv.solve(rec.t, (rec.I ** 2) * R / P.V_b, rec.T_inf, T0,
                            k=K_FIXED, h=H_FIXED, rho_cp=P.rho_cp())
        r = minimize_scalar(lambda R: float(np.sum((pred(R)["T_surf"] - rec.T_s) ** 2)),
                            bounds=(1e-4, 0.2), method="bounded")
        o = pred(r.x)
        es = np.sqrt(((o["T_surf"] - rec.T_s) ** 2).mean())
        ec = np.sqrt(((o["T_core"] - rec.T_c) ** 2).mean())
        print(f"      {lab:22s} R_eff {1000*r.x:7.4f} mOhm  "
              f"surf {es:.4f} K  CORE {ec:.4f} K")

print()
print("  Every other input is core-blind and was re-checked line by line:")
print("    h, k        Stage B [B2], fitted to DS1 SURFACE only, different record")
print("    R_eff       fitted to DS2 SURFACE only")
print("    dT_ref      T_s.max() - T_inf.mean()          surface + ambient")
print("    T_inf       measured ambient channel")
print("    I, V        measured electrical channels")
print("    selection   PDE + BC residual on a fixed collocation set, no data, no truth")
print("    seed choice lowest residual, NOT lowest core error")
print()
print("  VERDICT: the only core-channel contact was T0. Magnitude quantified above.")
print("           Stage E is re-run with the clean IC and the headline uses that.")
