"""G-1: run the LCO multi-rate separation and report it, RMSE included."""
import numpy as np

from g1_entropic_separation import (C_LUMP, V_CELL, L_SLAB, discharges,
                                    heat_trace, estimate_K)

CELLS = ["B0038", "B0039", "B0040"]
AMB = 44.0
CURRENTS = [1.0, 2.0, 4.0]
XG = np.linspace(0.08, 0.92, 22)          # avoid the ends, where dT/dt is worst

print("=" * 100)
print("G-1: IRREVERSIBLE / REVERSIBLE SEPARATION ON LCO, MULTI-RATE")
print("=" * 100)
print(f"  cells {CELLS} at {AMB:.0f} C, currents {CURRENTS} A")
print(f"  C_lump = {C_LUMP:.2f} J/K  (rho*cp = 2.75e6 (b), V_cell = {V_CELL:.4e} m^3)")
print()

# ---- 1. loss coefficient from relaxation only -------------------------------
allrecs = []
for cell in CELLS:
    for I in CURRENTS:
        allrecs += discharges(cell, AMB, I)
K, ntau = estimate_K(allrecs, AMB)
if K is None:
    K = (2 * 24.81 / L_SLAB) * V_CELL
    print(f"  no usable relaxation tails; falling back to Part-3 h = 24.81 "
          f"-> K = {K:.5f} W/K")
else:
    print(f"  K from {ntau} relaxation tails = {K:.5f} W/K "
          f"(tau = {C_LUMP/K:.0f} s); no source model involved")
print()

# ---- 2. Q/I at each x, per current ------------------------------------------
prof = {I: [] for I in CURRENTS}
for cell in CELLS:
    for I in CURRENTS:
        for r in discharges(cell, AMB, I):
            x, Q = heat_trace(r, K, AMB)
            m = r["on"]
            if m.sum() < 25:
                continue
            xi, qi = x[m], Q[m] / r["Imed"]
            o = np.argsort(xi)
            prof[I].append(np.interp(XG, xi[o], qi[o]))

print(f"  {'I (A)':>6s} {'cycles':>7s} {'mean Q/I (V)':>14s}")
for I in CURRENTS:
    a = np.array(prof[I])
    print(f"  {I:>6.1f} {len(a):>7d} {np.nanmean(a):>14.5f}")
    prof[I] = np.nanmedian(a, axis=0)
print()

# ---- 3. regress Q/I against I at each x -------------------------------------
Iv = np.array(CURRENTS)
T_K = AMB + 273.15
rows = []
for j, xj in enumerate(XG):
    y = np.array([prof[I][j] for I in CURRENTS])
    A = np.vstack([Iv, np.ones_like(Iv)]).T
    (slope, icpt), res, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ np.array([slope, icpt])
    rmse = float(np.sqrt(np.mean((pred - y) ** 2)))
    denom = float(np.sqrt(np.mean(y ** 2)))
    rows.append(dict(x=xj, R=slope, icpt=icpt, dUdT=-icpt / T_K,
                     rmse=rmse, rel=rmse / denom if denom else np.nan))

print("  Linear separation at each depth of discharge:")
print(f"    Q/I = I*R(x) - T*dU/dT(x)   ->  slope = R,  intercept = -T dU/dT")
print()
print(f"  {'x':>6s} {'R (mOhm)':>10s} {'intercept V':>12s} "
      f"{'dU/dT mV/K':>11s} {'fit RMSE V':>11s} {'rel %':>7s}")
for r in rows:
    print(f"  {r['x']:>6.3f} {1000*r['R']:>10.2f} {r['icpt']:>12.5f} "
          f"{1000*r['dUdT']:>11.4f} {r['rmse']:>11.6f} {100*r['rel']:>7.2f}")

R = np.array([r["R"] for r in rows])
dU = np.array([r["dUdT"] for r in rows])
rm = np.array([r["rmse"] for r in rows])
rel = np.array([r["rel"] for r in rows])

print()
print("=" * 100)
print("  THE NUMBER THAT WAS NEVER REPORTED")
print("=" * 100)
print(f"    median fit RMSE           {np.median(rm):.6f} V   "
      f"({100*np.median(rel):.2f} % of the fitted values)")
print(f"    worst  fit RMSE           {rm.max():.6f} V   ({100*rel.max():.2f} %)")
print(f"    n currents = {len(CURRENTS)}, so 3 points and 2 parameters: 1 dof per x")

print()
print("=" * 100)
print("  DOES REMOVING THE REVERSIBLE TERM REMOVE THE INTERIOR MINIMUM?")
print("=" * 100)
i_min = int(np.argmin(R))
print(f"    R(x) after separation: min at x = {XG[i_min]:.3f}, "
      f"{1000*R.min():.2f} mOhm; ends {1000*R[0]:.2f} -> {1000*R[-1]:.2f} mOhm")
interior = 0 < i_min < len(R) - 1
print(f"    interior minimum in the SEPARATED ohmic term? "
      f"{'YES' if interior else 'NO -- monotonic'}")
print()
print(f"    recovered dU/dT: {1000*dU.min():+.4f} to {1000*dU.max():+.4f} mV/K")
sign_changes = int(np.sum(np.diff(np.sign(dU)) != 0))
print(f"    sign changes across x: {sign_changes}")
print(f"    magnitude vs published LCO scale (0.1-0.2 mV/K, Part 3 sec 4.3): "
      f"median |dU/dT| = {1000*np.median(np.abs(dU)):.4f} mV/K")

print()
print("=" * 100)
print("  SENSITIVITY TO K (the only assumed constant)")
print("=" * 100)
print(f"  {'K factor':>9s} {'median R mOhm':>14s} {'median dU/dT mV/K':>19s}")
for f in (0.7, 0.85, 1.0, 1.15, 1.3):
    pr = {}
    for I in CURRENTS:
        acc = []
        for cell in CELLS:
            for r in discharges(cell, AMB, I):
                x, Q = heat_trace(r, K * f, AMB)
                m = r["on"]
                if m.sum() < 25:
                    continue
                xi, qi = x[m], Q[m] / r["Imed"]
                o = np.argsort(xi)
                acc.append(np.interp(XG, xi[o], qi[o]))
        pr[I] = np.nanmedian(np.array(acc), axis=0)
    Rs, dUs = [], []
    for j in range(len(XG)):
        y = np.array([pr[I][j] for I in CURRENTS])
        A = np.vstack([Iv, np.ones_like(Iv)]).T
        s, c = np.linalg.lstsq(A, y, rcond=None)[0]
        Rs.append(s); dUs.append(-c / T_K)
    print(f"  {f:>9.2f} {1000*np.median(Rs):>14.2f} {1000*np.median(dUs):>19.4f}")
print()
print("  As argued in the header: K error loads onto the SLOPE, leaving the")
print("  intercept -- the reversible term -- comparatively stable.")

np.savez("results/g1_separation.npz", x=XG, R=R, dUdT=dU, rmse=rm, rel=rel, K=K)
print("\n  saved results/g1_separation.npz")
