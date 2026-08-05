"""Verification of RadialFV against closed-form solutions.

If any of these fail, nothing downstream is trustworthy.  Run before Stage B.

Note on tolerances: the discretisation is constructed so that the analytic
uniform-generation steady state is reproduced EXACTLY (to round-off), because
  * interior conductance k*A_face/dr is exact for a parabolic-in-r^2 profile,
  * the area-exact half-cell conductance makes the wall flux exact,
  * the core extrapolation is exact for a profile linear in r^2.
So these are machine-precision tests, not loose ones.  A regression will show up
immediately rather than hiding under a 1% tolerance.
"""
import numpy as np
from part7_lib import (P, RadialFV, steady_core_minus_surf,
                       steady_surf_minus_inf, bi_over_2)

results = []


def check(name, got, want, tol, note=""):
    err = abs(got - want) / max(abs(want), 1e-30)
    ok = err <= tol
    results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:46s} got={got: .8g}  want={want: .8g}  "
          f"rel={err:.2e} (tol {tol:.0e}) {note}")


q = 5.0 / P.V_b          # 5 W over the cell volume, W/m^3
Tinf = 8.0

print("=" * 94)
print("TEST 1 -- steady state, uniform generation, vs exact solution (all N)")
print("=" * 94)
for N in (10, 20, 40, 80, 160):
    fv = RadialFV(N=N)
    t = np.arange(0, 20001, 1.0)
    out = fv.solve(t, np.full(len(t), q), np.full(len(t), Tinf), Tinf)
    d_cs = out["T_core"][-1] - out["T_surf"][-1]
    d_si = out["T_surf"][-1] - Tinf
    print(f"  N={N:4d}")
    check("core-surf   = q R^2 / 4k", d_cs, steady_core_minus_surf(q, P.k_t), 1e-9)
    check("surf-inf    = q R / 2h", d_si, steady_surf_minus_inf(q, P.h_lit), 1e-9)
    check("ratio       = Bi/2", d_cs / d_si, bi_over_2(P.h_lit, P.k_t), 1e-9)

print()
print("=" * 94)
print("TEST 2 -- global energy balance at steady state")
print("=" * 94)
for N in (20, 80):
    fv = RadialFV(N=N)
    t = np.arange(0, 20001, 1.0)
    out = fv.solve(t, np.full(len(t), q), np.full(len(t), Tinf), Tinf,
                   return_field=True)
    gen = q * fv.vol.sum()                              # W per unit height
    rej = fv.G_film * (out["T_surf"][-1] - Tinf)        # W per unit height
    check(f"N={N:3d}  generation == convection out", rej, gen, 1e-9)

print()
print("=" * 94)
print("TEST 3 -- profile shape is exactly parabolic in r^2")
print("=" * 94)
fv = RadialFV(N=40)
t = np.arange(0, 20001, 1.0)
o = fv.solve(t, np.full(len(t), q), np.full(len(t), Tinf), Tinf, return_field=True)
f = o["field"][-1]
# exact: T(r) = T_wall + q(R^2 - r^2)/(4k)
exact = o["T_surf"][-1] + q * (P.R_o ** 2 - fv.r ** 2) / (4 * P.k_t)
err = np.abs(f - exact).max()
results.append(err < 1e-9)
print(f"  [{'PASS' if err < 1e-9 else 'FAIL'}] max |T_cell - exact parabola| = {err:.3e} K")

print()
print("=" * 94)
print("TEST 4 -- adiabatic limit: all generated energy is stored")
print("=" * 94)
fv = RadialFV(N=80)
t = np.arange(0, 1001, 1.0)
o = fv.solve(t, np.full(len(t), q), np.full(len(t), Tinf), Tinf,
             h=1e-12, return_field=True)
Tbar = (o["field"] * fv.vol).sum(1) / fv.vol.sum()
check("mean rise = q t /(rho cp)", Tbar[-1] - Tinf, q * t[-1] / P.rho_cp(), 1e-9)

print()
print("=" * 94)
print("TEST 5 -- pure cooling at small Bi reproduces the lumped time constant")
print("=" * 94)
h_small = 0.5
fv = RadialFV(N=80)
t = np.arange(0, 4001, 1.0)
o = fv.solve(t, np.zeros(len(t)), np.zeros(len(t)), 10.0, h=h_small)
tau_pred = P.rho_cp() / (h_small * 2.0 / P.R_o)
Tm = o["T_surf"]
i1, i2 = 200, 3000
tau_fit = (t[i2] - t[i1]) / np.log(Tm[i1] / Tm[i2])
check("tau = rho cp V /(h A)", tau_fit, tau_pred, 3e-2,
      f"(Bi={P.biot(h_small):.4f}, lumped only approximate)")

print()
print("=" * 94)
print("TEST 6 -- symmetry at r=0 holds by construction")
print("=" * 94)
fv = RadialFV(N=40)
ok = fv.area[0] == 0.0
results.append(ok)
print(f"  [{'PASS' if ok else 'FAIL'}] inner face area identically zero -> no axis flux term")
_ts = np.arange(0, 20001, 1.0)
o = fv.solve(_ts, np.full(len(_ts), q), np.full(len(_ts), Tinf), Tinf,
             return_field=True)
f = o["field"][-1]
ok = np.all(np.diff(f) < 0)
results.append(ok)
print(f"  [{'PASS' if ok else 'FAIL'}] profile monotone decreasing outward")
# Discrete gradient at face 1 (r=dr) vs face N-1 (r=R-dr).  The identity
# ratio = dr/(R-dr) comes from the STEADY parabola dT/dr = -q r/(2k), so this
# must be evaluated at steady state -- at t=500 s (~1.2 tau) it is not yet valid.
g_in = (f[1] - f[0]) / fv.dr
g_out = (f[-1] - f[-2]) / fv.dr
check("grad(face 1)/grad(face N-1) = dr/(R-dr)", g_in / g_out,
      fv.dr / (P.R_o - fv.dr), 1e-9)

print()
print("=" * 94)
print("TEST 7 -- transient dt refinement (backward Euler, 1st order in time)")
print("=" * 94)
vals = []
for dt in (4.0, 2.0, 1.0, 0.5, 0.25, 0.125):
    fv = RadialFV(N=40)
    t = np.arange(0, 601, dt)
    qq = q * (1.0 + np.sin(2 * np.pi * t / 120.0))
    o = fv.solve(t, qq, np.full(len(t), Tinf), Tinf)
    v = o["T_core"][-1] - o["T_surf"][-1]
    vals.append(v)
    print(f"  dt={dt:6.3f}  core-surf at t=600 s = {v:.6f} K")
diffs = [abs(vals[i] - vals[-1]) for i in range(len(vals) - 1)]
print(f"  |value - finest| : {np.array(diffs).round(5)}")
rate = np.log2(diffs[0] / diffs[1]) if diffs[1] > 0 else np.nan
print(f"  refinement rate dt=4->2 : {rate:.2f}  (backward Euler => ~1)")
err_at_1s = diffs[2]
ok = err_at_1s < 0.05
results.append(ok)
print(f"  [{'PASS' if ok else 'FAIL'}] dt=1 s (the data's native rate) is within "
      f"{err_at_1s:.4f} K of the dt=0.125 s answer")

print()
print("=" * 94)
n_ok, n = sum(results), len(results)
print(f"SUMMARY: {n_ok}/{n} checks passed"
      f"{'  -- ALL GOOD' if n_ok == n else '  -- INVESTIGATE FAILURES'}")
print("=" * 94)
