"""Stage G diagnostic -- how much physics is the PINN actually enforcing?

The finite-volume solver satisfies the PDE exactly by construction, so its
surface RMSE is the best achievable WITHIN the model class.  If the PINN fits
the surface better than that, it is not being cleverer about the physics -- it
is violating the PDE to chase the data, because the constraint is soft.

That matters here specifically: a network that bends the physics to fit the
surface has no reason to extrapolate correctly to the core, which is the one
thing we are asking it to do.

Test: take the PINN's recovered R_eff, push it through the EXACT solver, and
compare.  Divergence between the two is PINN physics violation, measured in K.
"""
import numpy as np

from part7_lib import P, Record, RadialFV
from stage_e_inverse import H_FIXED, K_FIXED

rec = Record("2")
fv = RadialFV(N=40)
T0 = float(rec.T_s[0])   # surface only: no core leak
grad = rec.T_c - rec.T_s

print("=" * 98)
print("STAGE G -- physics-consistency of the PINN solution")
print("=" * 98)

for n_shape in (0, 1):
    try:
        d = np.load(f"results/stage_e_shape{n_shape}.npz")
    except FileNotFoundError:
        continue
    b = int(np.argmin(d["sel"]))
    R = float(d["R_eff"][b])
    Tc_p, Ts_p = d["Tc"][b], d["Ts"][b]

    o = fv.solve(rec.t, (rec.I ** 2) * R / P.V_b, rec.T_inf, T0,
                 k=K_FIXED, h=H_FIXED, rho_cp=P.rho_cp())

    print(f"\n  --- shape order {n_shape}, selected seed {b}, "
          f"R_eff = {1000*R:.4f} mOhm ---")
    print(f"    {'':28s} {'surf RMSE':>10s} {'core RMSE':>10s} {'grad RMSE':>10s}")
    for lab, Tc, Ts in (("PINN solution", Tc_p, Ts_p),
                        ("same R_eff through exact FV", o["T_core"], o["T_surf"])):
        es = np.sqrt(((Ts - rec.T_s) ** 2).mean())
        ec = np.sqrt(((Tc - rec.T_c) ** 2).mean())
        eg = np.sqrt((((Tc - Ts) - grad) ** 2).mean())
        print(f"    {lab:28s} {es:10.4f} {ec:10.4f} {eg:10.4f}")

    dS = np.sqrt(((Ts_p - o["T_surf"]) ** 2).mean())
    dC = np.sqrt(((Tc_p - o["T_core"]) ** 2).mean())
    print(f"    PINN vs exact FV at the SAME R_eff: surface {dS:.4f} K, core {dC:.4f} K")
    print(f"    -> that difference IS the PDE violation, expressed in kelvin.")

    surf_rise = Ts_p - rec.T_inf
    m = surf_rise > 0.5 * surf_rise.max()
    ratio = float(((Tc_p - Ts_p)[m] / surf_rise[m]).mean())
    Bi = H_FIXED * P.R_o / K_FIXED
    print(f"    plausibility (trap 5.7): ratio {ratio:.4f} vs Bi/2 {Bi/2:.4f} "
          f"-> factor {ratio/(Bi/2):.3f} "
          f"({'PASS' if 1/3 < ratio/(Bi/2) < 3 else 'REJECT'})")

# measured reality check on the same statistic
m = (rec.T_s - rec.T_inf) > 0.5 * (rec.T_s - rec.T_inf).max()
print(f"\n  MEASURED ratio over the same hot window: "
      f"{float((grad[m]/(rec.T_s-rec.T_inf)[m]).mean()):.4f}   "
      f"(Bi/2 = {H_FIXED*P.R_o/K_FIXED/2:.4f})")
