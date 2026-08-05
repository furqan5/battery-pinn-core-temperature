"""Stage G -- skill against trivial baselines (trap 5.8, compare like with like).

The gradient figure shows the predicted core-minus-surface tracking the slow
envelope of the measured gradient but missing its fast oscillations almost
entirely.  A 0.61 K core RMSE sounds good, but "good" is meaningless without
asking: good compared to WHAT?

Three baselines that require no PINN, no PDE, and no fitting beyond a constant:

  B0  T_core = T_surf                      (assume no gradient at all)
  B1  T_core = T_surf + mean(measured gradient)
  B2  T_core = T_surf + Bi/2 * (T_surf - T_inf)   (the steady analytic relation,
                                                   using only surface + ambient
                                                   and the fixed Bi)

B1 uses the measured gradient's mean, so it is NOT core-blind -- it is an
unfairly strong baseline, deliberately.  If the PINN cannot beat a baseline that
was handed the answer's mean, the claim is weak.
"""
import numpy as np

from part7_lib import P, Record
from stage_e_inverse import H_FIXED, K_FIXED


def report(tag, Tc_pred, label_pred, h, k):
    rec = Record(tag)
    grad = rec.T_c - rec.T_s
    Bi = h * P.R_o / k

    cands = {
        label_pred:                       Tc_pred,
        "B0: T_core = T_surf":            rec.T_s.copy(),
        "B1: T_surf + mean(grad)*":       rec.T_s + grad.mean(),
        "B2: T_surf + Bi/2*(Ts-Tinf)":    rec.T_s + (Bi / 2) * (rec.T_s - rec.T_inf),
    }

    print(f"\n  DATASET {tag}   (measured gradient: mean {grad.mean():.3f} K, "
          f"sd {grad.std():.3f} K, max {grad.max():.3f} K)")
    print(f"    {'predictor':34s} {'core RMSE':>10s} {'grad RMSE':>10s} "
          f"{'skill vs B1':>12s}")
    b1 = float(np.sqrt(((cands['B1: T_surf + mean(grad)*'] - rec.T_c) ** 2).mean()))
    for name, Tc in cands.items():
        e = Tc - rec.T_c
        rm = float(np.sqrt((e ** 2).mean()))
        eg = float(np.sqrt((((Tc - rec.T_s) - grad) ** 2).mean()))
        skill = 1.0 - rm / b1
        print(f"    {name:34s} {rm:>10.4f} {eg:>10.4f} {100*skill:>11.1f} %")

    # How much of the gradient's VARIABILITY (not its level) is captured?
    gp = Tc_pred - rec.T_s
    ss_res = float(((gp - grad) ** 2).sum())
    ss_tot = float(((grad - grad.mean()) ** 2).sum())
    print(f"    variance of the gradient explained (R^2 vs its own mean): "
          f"{1 - ss_res/ss_tot:+.3f}")
    # correlation of the fluctuations
    print(f"    corr(predicted grad, measured grad) = "
          f"{np.corrcoef(gp, grad)[0,1]:+.4f}")
    print(f"    sd(predicted grad) = {gp.std():.3f} K vs measured {grad.std():.3f} K"
          f"  -> the prediction is {100*gp.std()/grad.std():.0f} % as variable")
    return b1


print("=" * 92)
print("STAGE G -- skill against trivial baselines")
print("=" * 92)
print("  * B1 is handed the measured gradient's MEAN, so it is not core-blind.")
print("    It is included as a deliberately unfair, strong baseline.")

d = np.load("results/stage_e_shape0.npz")
b = int(np.argmin(d["sel"]))
report("2", d["Tc"][b], "inverse PINN (order 0)", H_FIXED, K_FIXED)

try:
    d1 = np.load("results/stage_e_ds1_shape0.npz")
    b1i = int(np.argmin(d1["sel"]))
    report("1", d1["Tc"][b1i], "inverse PINN (order 0)", 37.2846, 0.418697)
except FileNotFoundError:
    print("\n  (DS1 replication not finished yet)")

print()
print("  Reading: a high core R^2 driven mostly by the gradient's MEAN LEVEL is a")
print("  much weaker claim than one driven by its DYNAMICS.  The honest statement")
print("  has to separate the two.")
