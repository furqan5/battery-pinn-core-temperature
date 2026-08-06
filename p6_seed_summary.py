"""Summarise the Part 6 seed sweep.

Separate from p6_seeds.py because the 12 trainings are expensive and must not be
repeated to fix a formatting bug. (The sweep itself completed; only the summary
line failed, on ndarray.ptp() which NumPy 2.0 removed.)
"""
import json
import numpy as np

TARGET_R0 = 144.352      # cycle-0 classical, the like-for-like comparator
TARGET_C1 = -1.305
QUOTED = 144.78          # the value Part 6 reports, scored as G2 HIT at +0.30 %

res = json.load(open("results/p6_seeds.json"))
print("=" * 96)
print("PART 6 SEED SWEEP — how stable is the fixed-residual fit?")
print("=" * 96)
print(f"  comparator: cycle-0 classical R0 = {TARGET_R0:.3f} mOhm, c1 = {TARGET_C1:.3f}")
print(f"  Part 6 quotes {QUOTED:.2f} mOhm (+{100*(QUOTED-TARGET_R0)/TARGET_R0:.2f} %) "
      f"and scores G2 a HIT on it.")
print()

for w in sorted({r["w_data"] for r in res}):
    g = sorted([r for r in res if r["w_data"] == w], key=lambda r: r["seed"])
    R = np.array([r["R0_mOhm"] for r in g])
    c1 = np.array([r["c1"] for r in g])
    rm = np.array([r["rmse"] for r in g])
    Lr = np.array([r["Lr"] for r in g])
    dev = 100 * (R - TARGET_R0) / TARGET_R0

    print("-" * 96)
    print(f"  w_data = {w:.0f}")
    print("-" * 96)
    print(f"    {'seed':>5s} {'R0 mOhm':>9s} {'dev %':>8s} {'c1':>8s} "
          f"{'RMSE K':>8s} {'Lr':>11s}  {'G2':>5s}")
    for r, d_ in zip(g, dev):
        hit = abs(d_) <= 3.0
        print(f"    {r['seed']:>5d} {r['R0_mOhm']:>9.3f} {d_:>+8.2f} "
              f"{r['c1']:>8.3f} {r['rmse']:>8.4f} {r['Lr']:>11.3e}  "
              f"{'HIT' if hit else 'MISS':>5s}")

    print(f"    all 6 : R0 {R.mean():.3f} +/- {R.std(ddof=1):.3f} mOhm "
          f"(range {R.min():.3f}-{R.max():.3f}, spread {np.ptp(R):.3f})")
    print(f"            deviation {dev.mean():+.2f} % "
          f"[{dev.min():+.2f}, {dev.max():+.2f}]")
    print(f"    G2 (|dev| <= 3 %) scores HIT on {int((np.abs(dev)<=3).sum())}/6 seeds")

    # truth-free selection: lowest PDE residual, the rule this project mandates
    b = int(np.argmin(Lr))
    print(f"    TRUTH-FREE PICK (lowest Lr): seed {g[b]['seed']} -> "
          f"R0 {R[b]:.3f} mOhm ({dev[b]:+.2f} %) -> "
          f"G2 {'HIT' if abs(dev[b]) <= 3 else 'MISS'}")

    # does Lr flag the outlier?
    worst = int(np.argmax(np.abs(dev)))
    print(f"    worst-parameter seed is {g[worst]['seed']} ({dev[worst]:+.2f} %); "
          f"its Lr is {Lr[worst]:.3e}, "
          f"{Lr[worst]/np.median(Lr):.1f}x the median -> "
          f"{'FLAGGED by the residual' if Lr[worst] == Lr.max() else 'NOT flagged'}")

    keep = np.arange(len(R)) != worst
    print(f"    excluding it: R0 {R[keep].mean():.3f} +/- {R[keep].std(ddof=1):.3f} "
          f"mOhm ({dev[keep].mean():+.2f} % +/- {dev[keep].std(ddof=1):.2f})")
    print(f"    c1 all negative: {bool((c1 < 0).all())}   "
          f"mean {c1.mean():+.3f} +/- {c1.std(ddof=1):.3f}")
    print(f"    RMSE {rm.mean():.4f} +/- {rm.std(ddof=1):.4f} K")
    print()

print("=" * 96)
print("VERDICT")
print("=" * 96)
w20 = [r for r in res if r["w_data"] == 20.0]
R20 = np.array([r["R0_mOhm"] for r in w20])
d20 = 100 * (R20 - TARGET_R0) / TARGET_R0
print(f"  The quoted {QUOTED:.2f} mOhm sits at the favourable edge of the distribution:")
print(f"  {int((R20 < QUOTED).sum())}/6 seeds fall below it, "
      f"{int((R20 > QUOTED).sum())}/6 above.")
print(f"  G2 as scored in the notebook (HIT at +0.30 %) is not representative;")
print(f"  it scores HIT on {int((np.abs(d20)<=3).sum())}/6 seeds.")
print()
print("  What DOES hold across every seed:")
print("    - c1 is negative in 12/12 runs (the broken residual gave +0.023)")
print("    - RMSE ~0.13 K at w=20 and ~0.07 K at w=200, against 0.62 K broken")
print("    - the shape is recovered; only the AMPLITUDE is seed-sensitive")
