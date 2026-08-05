"""Like-for-like classical control for the PINN's core prediction (trap 5.8).

The PINN must be compared against a baseline solving the SAME problem with the
SAME information: surface trace only, h and k fixed from DS1's surface-only fit,
scalar R_eff, fitted on DS2.  Anything else flatters one side.

This also isolates how much of the headline number is the PINN and how much is
simply the physics plus a well-chosen k.
"""
import numpy as np
from scipy.optimize import minimize_scalar

from part7_lib import P, Record, RadialFV
from stage_e_inverse import H_FIXED, K_FIXED

rec = Record("2")
fv = RadialFV(N=40)
T0 = float(rec.T_s[0])   # surface only: no core leak
grad = rec.T_c - rec.T_s


def predict(R_eff, h=H_FIXED, k=K_FIXED, rho_cp=None):
    if rho_cp is None:
        rho_cp = P.rho_cp()
    q = (rec.I ** 2) * R_eff / P.V_b
    return fv.solve(rec.t, q, rec.T_inf, T0, k=k, h=h, rho_cp=rho_cp)


def sse(R):
    return float(np.sum((predict(R)["T_surf"] - rec.T_s) ** 2))


r = minimize_scalar(sse, bounds=(1e-4, 0.2), method="bounded",
                    options={"xatol": 1e-10})
out = predict(r.x)
e = out["T_core"] - rec.T_c
es = out["T_surf"] - rec.T_s
eg = (out["T_core"] - out["T_surf"]) - grad

print("=" * 92)
print("CLASSICAL CONTROL -- identical information, identical fixed parameters")
print("=" * 92)
print(f"  h = {H_FIXED:.4f}, k = {K_FIXED:.6f} (DS1 surface-only), fitted on DS2 surface")
print(f"  R_eff              = {1000*r.x:.4f} mOhm")
print(f"  surface RMSE       = {np.sqrt((es**2).mean()):.4f} K")
print(f"  CORE RMSE          = {np.sqrt((e**2).mean()):.4f} K   max {np.abs(e).max():.4f} K")
print(f"  gradient RMSE      = {np.sqrt((eg**2).mean()):.4f} K")
print(f"  core/surface ratio = {np.sqrt((e**2).mean())/np.sqrt((es**2).mean()):.2f}")

# ---- sensitivity of the CORE prediction to k, the parameter we cannot identify ----
print()
print("  Sensitivity of the core prediction to k (refitting R_eff at each k):")
print(f"    {'k':>8s} {'Bi':>7s} {'R_eff mOhm':>11s} {'surf RMSE':>10s} {'CORE RMSE':>10s}")
for k in (0.30, 0.35, K_FIXED, 0.404, 0.45, 0.55):
    rr = minimize_scalar(lambda R: float(np.sum((predict(R, k=k)["T_surf"] - rec.T_s) ** 2)),
                         bounds=(1e-4, 0.2), method="bounded")
    o = predict(rr.x, k=k)
    ee = o["T_core"] - rec.T_c
    ss = o["T_surf"] - rec.T_s
    print(f"    {k:>8.4f} {H_FIXED*P.R_o/k:>7.3f} {1000*rr.x:>11.4f} "
          f"{np.sqrt((ss**2).mean()):>10.4f} {np.sqrt((ee**2).mean()):>10.4f}")
print("  -> the surface fit barely moves while the core error swings widely.")
print("     That is the whole identifiability problem in one table: k is an INPUT.")

# ---- entropic-term sensitivity (headline keeps it off; this quantifies it) ----
print()
print("  Entropic-term sensitivity (Forgez et al. 2010, dU/dT = -0.5 mV/K at 50% SOC):")
Tc_K = rec.T_c + 273.15


def predict_with_rev(R_eff, dUdT):
    q_irr = (rec.I ** 2) * R_eff
    # I>0 = charge here, so Q_rev = +I T dU/dT (endothermic on charge)
    q_rev = rec.I * (rec.T_s + 273.15) * dUdT
    q = np.clip(q_irr + q_rev, 0.0, None) / P.V_b
    return fv.solve(rec.t, q, rec.T_inf, T0, k=K_FIXED, h=H_FIXED, rho_cp=P.rho_cp())


for dUdT in (0.0, P.dUdT_50):
    rr = minimize_scalar(
        lambda R: float(np.sum((predict_with_rev(R, dUdT)["T_surf"] - rec.T_s) ** 2)),
        bounds=(1e-4, 0.2), method="bounded")
    o = predict_with_rev(rr.x, dUdT)
    ee = o["T_core"] - rec.T_c
    ss = o["T_surf"] - rec.T_s
    lab = "off (headline)" if dUdT == 0.0 else f"on, {1000*dUdT:+.2f} mV/K"
    print(f"    {lab:<22s} R_eff {1000*rr.x:7.4f} mOhm  surf {np.sqrt((ss**2).mean()):.4f} K"
          f"  CORE {np.sqrt((ee**2).mean()):.4f} K")
print("  NOTE: T_s is used in the reversible term, not T_c -- using the core")
print("        temperature there would be a leak.  A fixed dU/dT at 50% SOC is an")
print("        approximation; the SOC window here is only 11% wide, which limits")
print("        the error that introduces.")
