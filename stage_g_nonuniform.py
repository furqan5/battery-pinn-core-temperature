"""Stage G diagnostic -- does NON-UNIFORM generation explain the PINN's behaviour?

Two facts need one explanation:
  * the inverse PINN recovers R_eff = 13.10 mOhm, ~9% BELOW both independent
    anchors (classical thermal fit 14.30, electrical regression 14.49);
  * yet its core prediction is 31% BETTER (0.614 K vs 0.894 K).

Under uniform generation those two should not go together: less heat means a
smaller gradient means a worse core.  So the uniform assumption is suspect.

Hypothesis: real generation is concentrated toward the axis (current collector,
tab and winding-core losses).  For a FIXED surface trace that requires LESS total
heat while producing a LARGER centre-to-surface gradient -- which is exactly the
pair of symptoms observed.

Test it classically, where nothing is opaque: redistribute generation with a
volume-preserving weight w(r) = 1 + beta(1 - 2(r/R)^2), refit R_eff to the
SURFACE at each beta, and watch R_eff and the core error move.

THIS IS A DIAGNOSTIC, NOT A FIT.  beta is scanned to explain a mechanism.  It is
NOT selected on core error and does NOT feed the headline, which keeps the
uniform-generation assumption.  Scanning beta against the core and then reporting
that core number would be textbook test-set selection.
"""
import numpy as np
from scipy.optimize import minimize_scalar

from part7_lib import P, Record, RadialFV
from stage_e_inverse import H_FIXED, K_FIXED

rec = Record("2")
fv = RadialFV(N=40)
T0 = float(rec.T_s[0])
grad = rec.T_c - rec.T_s

print("=" * 100)
print("STAGE G -- is non-uniform generation the explanation?")
print("=" * 100)
print(f"  targets to explain:  PINN R_eff 13.096 mOhm, core RMSE 0.614 K")
print(f"  uniform classical :       R_eff 14.301 mOhm, core RMSE 0.894 K")
print(f"  electrical anchor :       R_eff {1000*rec.R_ohmic_reg:.3f} mOhm")
print()
print(f"  {'beta':>6s} {'q(0)/q(R)':>10s} {'R_eff mOhm':>11s} {'surf RMSE':>10s} "
      f"{'CORE RMSE':>10s} {'grad RMSE':>10s}")

for beta in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6):
    w = fv.radial_weight(beta)

    def pred(R):
        return fv.solve(rec.t, (rec.I ** 2) * R / P.V_b, rec.T_inf, T0,
                        k=K_FIXED, h=H_FIXED, rho_cp=P.rho_cp(), q_weight=w)

    r = minimize_scalar(lambda R: float(np.sum((pred(R)["T_surf"] - rec.T_s) ** 2)),
                        bounds=(1e-4, 0.2), method="bounded")
    o = pred(r.x)
    es = np.sqrt(((o["T_surf"] - rec.T_s) ** 2).mean())
    ec = np.sqrt(((o["T_core"] - rec.T_c) ** 2).mean())
    eg = np.sqrt((((o["T_core"] - o["T_surf"]) - grad) ** 2).mean())
    print(f"  {beta:>6.2f} {w[0]/w[-1]:>10.3f} {1000*r.x:>11.4f} {es:>10.4f} "
          f"{ec:>10.4f} {eg:>10.4f}")

print()
print("  Reading:")
print("    - If R_eff FALLS toward ~13.1 mOhm as beta rises while the surface fit")
print("      stays flat, then centre-weighted generation reproduces the PINN's")
print("      lower R_eff without any neural network being involved.")
print("    - If the core error also falls toward ~0.61 K, the same mechanism")
print("      explains the better core prediction.")
print("    - The surface RMSE barely moving across beta is the point: the surface")
print("      CANNOT see the radial distribution. That is why beta is not")
print("      identifiable from this data, and why it stays out of the headline.")
