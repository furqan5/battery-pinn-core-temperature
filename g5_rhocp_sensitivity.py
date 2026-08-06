"""G-5: how much does the assumed volumetric heat capacity actually matter?

Resolution of the flagged discrepancy. GAPS.md recorded a 5.1 % conflict between
the mass implied by Richardson's parameters (72.1 g) and a manufacturer sheet
giving 76 g. That comparison was wrong: 76 g belongs to the ANR26650M1-B
(2.6 Ah), a later variant, while Richardson's data is on the ANR26650M1
(2.3 Ah), whose datasheet gives a core cell weight of 70 g. Against the correct
figure the implied mass is +3.0 %, which is what an effective bulk density
should look like.

Two candidates therefore bracket the truth:
    rho*cp (Richardson)     = 2107 * 1171.6      = 2.4685e6 J/m^3/K
    rho*cp (70 g datasheet) = (0.070/V_b) * 1171.6 = 2.3966e6 J/m^3/K   (-2.9 %)

This script reports what changes, and what does not.
"""
import numpy as np
from scipy.optimize import minimize_scalar

from part7_lib import P, Record, RadialFV
from stage_e_inverse import H_FIXED, K_FIXED

M_DATASHEET = 0.070          # kg, A123 ANR26650M1 datasheet "core cell weight"
CP = P.cp                    # J/kg/K, Richardson

RHOCP = {
    "Richardson (rho=2107)": P.rho_cp(),
    "datasheet mass 70 g": (M_DATASHEET / P.V_b) * CP,
}

fv = RadialFV(N=40)
print("=" * 92)
print("G-5: sensitivity of every reported result to the assumed rho*cp")
print("=" * 92)
print(f"  cell volume V_b = {P.V_b:.5e} m^3")
for k, v in RHOCP.items():
    print(f"  {k:26s} rho*cp = {v:.4e} J/m^3/K   implied mass "
          f"{v/CP*P.V_b*1000:.2f} g   C_lump {v*P.V_b:.2f} J/K")
d = list(RHOCP.values())
print(f"  spread: {100*(d[0]-d[1])/d[1]:+.2f} %")
print()

for tag, h, k in (("2", H_FIXED, K_FIXED), ("1", 37.2846, 0.418697)):
    rec = Record(tag)
    T0 = float(rec.T_s[0])
    Bi = h * P.R_o / k
    print("-" * 92)
    print(f"  DATASET {tag}   Bi = {Bi:.4f}")
    print("-" * 92)

    # (a) the quasi-steady relation: does rho*cp appear at all?
    e = rec.T_s + (Bi / 2) * (rec.T_s - rec.T_inf) - rec.T_c
    qs = float(np.sqrt((e ** 2).mean()))
    print(f"    QUASI-STEADY relation core RMSE = {qs:.4f} K")
    print(f"      -> rho*cp does not appear in T_surf + (Bi/2)(T_surf - T_inf).")
    print(f"         The headline result is INDEPENDENT of this assumption.")

    # (b) the transient solver
    print(f"    {'rho*cp source':26s} {'R_eff mOhm':>11s} {'surf RMSE':>10s} "
          f"{'core RMSE':>10s} {'tau_diff s':>11s}")
    for name, rc in RHOCP.items():
        def pred(R):
            return fv.solve(rec.t, (rec.I ** 2) * R / P.V_b, rec.T_inf, T0,
                            k=k, h=h, rho_cp=rc)
        r = minimize_scalar(
            lambda R: float(np.sum((pred(R)["T_surf"] - rec.T_s) ** 2)),
            bounds=(1e-4, 0.2), method="bounded")
        o = pred(r.x)
        sr = float(np.sqrt(((o["T_surf"] - rec.T_s) ** 2).mean()))
        cr = float(np.sqrt(((o["T_core"] - rec.T_c) ** 2).mean()))
        tau = P.R_o ** 2 * rc / k
        print(f"    {name:26s} {1000*r.x:>11.4f} {sr:>10.4f} {cr:>10.4f} "
              f"{tau:>11.0f}")
    print()

# (c) effect on the timescale criterion
print("=" * 92)
print("  EFFECT ON THE TIMESCALE CRITERION")
print("=" * 92)
tau_a = P.R_o ** 2 * RHOCP["Richardson (rho=2107)"] / P.k_t
tau_b = P.R_o ** 2 * RHOCP["datasheet mass 70 g"] / P.k_t
print(f"    tau_diff: {tau_a:.0f} s (Richardson) vs {tau_b:.0f} s (datasheet), "
      f"{100*(tau_b-tau_a)/tau_a:+.2f} %")
print(f"    The criterion is eps = 7.48 (t/tau)^-1.18. A {100*(tau_b-tau_a)/tau_a:+.2f} %")
print(f"    shift in tau moves every plotted ratio by the same factor, so the")
print(f"    FITTED EXPONENT IS UNCHANGED and the crossings move by that amount:")
for tgt, ratio in ((2.0, 3.49), (5.0, 1.37), (10.0, 0.77), (25.0, 0.37)):
    print(f"      {tgt:5.1f} % error at t/tau = {ratio:.2f} "
          f"-> {ratio*tau_a/tau_b:.2f} under the datasheet value")

print()
print("=" * 92)
print("  INDEPENDENT ANCHOR FROM THE SAME DATASHEET")
print("=" * 92)
r2 = Record("2"); r1 = Record("1")
print(f"    datasheet internal impedance (1 kHz AC, 25 C): 8 mOhm typical")
print(f"    electrical regression of V on I, this work:    "
      f"{1000*r1.R_ohmic_reg:.2f} mOhm (DS1), {1000*r2.R_ohmic_reg:.2f} mOhm (DS2), DC at 8 C")
print(f"    recovered R_eff, classical:                    14.30 mOhm (DS2)")
print(f"    A 1 kHz AC value excludes charge-transfer and diffusion and is taken")
print(f"    at 25 C, so it must sit BELOW a DC value at 8 C. It does, by ~1.8x.")
print(f"    Consistent, and a third independent check on the absolute scale.")
