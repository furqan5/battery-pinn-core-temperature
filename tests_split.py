"""Verification of the lumped base used by the split formulation.

The claim in pinn_split.py is that T = T_l + w is EXACT ALGEBRA, not an
approximation.  That claim rests on T_l being the true solution of the 0-D
model, so it gets tested rather than asserted.
"""
import numpy as np
import torch

from part7_lib import P, Record, RadialFV
from pinn_split import LumpedBase

results = []


def check(name, got, want, tol, note=""):
    err = abs(got - want) / max(abs(want), 1e-30)
    ok = err <= tol
    results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:52s} got={got: .8g} want={want: .8g} "
          f"rel={err:.2e} {note}")


rec = Record("2")
h, rho_cp, R, V_b = P.h_lit, P.rho_cp(), P.R_o, P.V_b
T0 = 0.5 * (rec.T_s[0] + rec.T_c[0])

print("=" * 96)
print("TEST 1 -- constant source: closed-form first-order response")
print("=" * 96)
# Build a record-like object with constant current so q is constant.
import copy
rc = copy.copy(rec)
rc.I = np.full(len(rec.t), 10.0)
rc.T_inf = np.full(len(rec.t), 8.0)
R_eff = 0.0143398
base = LumpedBase(rc, h, rho_cp, R, V_b, T0, n_shape=0)
Tl = base.T_l_np(R_eff)
tau = rho_cp * R / (2 * h)
q = 100.0 * R_eff / V_b                      # I^2 R / V
exact = 8.0 + (T0 - 8.0) * np.exp(-rec.t / tau) + \
        (q * tau / rho_cp) * (1 - np.exp(-rec.t / tau))
check("tau = rho cp R /(2h)", base.tau, tau, 1e-12)
check("T_l at t=500 s", Tl[500], exact[500], 1e-9)
check("T_l at t=2000 s", Tl[2000], exact[2000], 1e-9)
check("T_l steady limit = T_inf + q tau/(rho cp)",
      Tl[-1], exact[-1], 1e-9, f"(t/tau = {rec.t[-1]/tau:.1f})")

print()
print("=" * 96)
print("TEST 2 -- linearity in R_eff (the property superposition relies on)")
print("=" * 96)
b2 = LumpedBase(rec, h, rho_cp, R, V_b, T0, n_shape=0)
a, bb = 0.010, 0.020
Ta, Tb, Tm = b2.T_l_np(a), b2.T_l_np(bb), b2.T_l_np(0.5 * (a + bb))
dev = float(np.abs(Tm - 0.5 * (Ta + Tb)).max())
results.append(dev < 1e-12)
print(f"  [{'PASS' if dev < 1e-12 else 'FAIL'}] affine in R_eff (midpoint test): "
      f"max deviation {dev:.3e} K  (tolerance 1e-12)")

print()
print("=" * 96)
print("TEST 3 -- torch path agrees with the numpy path")
print("=" * 96)
Tl_np = b2.T_l_np(R_eff)
Tl_t = b2.T_l_torch(torch.tensor(R_eff)).numpy()
check("T_l_torch == T_l_np", float(np.abs(Tl_t - Tl_np).max()), 0.0, 1e-12,
      "(max abs difference)")
results.append(bool(np.abs(Tl_t - Tl_np).max() < 1e-12))
print(f"  [{'PASS' if results[-1] else 'FAIL'}] max difference "
      f"{np.abs(Tl_t-Tl_np).max():.3e} K")

print()
print("=" * 96)
print("TEST 4 -- lumped limit: as k -> inf (Bi -> 0), the FULL solution -> T_l")
print("=" * 96)
fv = RadialFV(N=40)
q_vol = (rec.I ** 2) * R_eff / V_b
devs = []
for k in (0.404, 10.0, 1000.0, 1e5):
    o = fv.solve(rec.t, q_vol, rec.T_inf, T0, k=k, h=h, rho_cp=rho_cp)
    d = float(np.sqrt(((o["T_surf"] - Tl_np) ** 2).mean()))
    devs.append(d)
    print(f"    k={k:9.1f}  Bi={h*R/k:9.5f}   rms(T_surf_FV - T_l) = {d:.6f} K")

# It converges to ~0.0025 K, NOT to zero, and that floor is explained:
# RadialFV steps with BACKWARD EULER while LumpedBase uses EXACT EXPONENTIAL
# integration.  Both solve the same ODE in this limit; they differ by the
# time-integration error at dt/tau = 1/405.  Confirm by reproducing the floor
# with a backward-Euler lumped integration.
tau = rho_cp * R / (2 * h)
Tbe = np.empty(len(rec.t)); Tbe[0] = T0
c = rec.dt / tau
for i in range(1, len(rec.t)):
    Tbe[i] = (Tbe[i-1] + c * rec.T_inf[i] + rec.dt * q_vol[i] / rho_cp) / (1 + c)
floor = float(np.sqrt(((Tbe - Tl_np) ** 2).mean()))
print(f"    backward-Euler lumped vs exponential lumped: {floor:.6f} K")
ok = abs(devs[-1] - floor) < 5e-4
results.append(ok)
print(f"  [{'PASS' if ok else 'FAIL'}] the {devs[-1]:.5f} K residual IS the "
      f"integrator difference, not a formulation error")
print(f"      (0.0025 K against a 6.5 K measured gradient = 0.04% -- negligible here,")
print(f"       but it is why FV and the split PINN can never agree to machine precision)")

print()
print("=" * 96)
print("TEST 5 -- gradient of T_l w.r.t. R_eff is exact (autograd vs analytic)")
print("=" * 96)
Rt = torch.tensor(R_eff, requires_grad=True)
out = b2.T_l_torch(Rt).sum()
out.backward()
analytic = b2.B[0].sum()
check("d(sum T_l)/d(R_eff)", float(Rt.grad), float(analytic), 1e-10)

print()
print("=" * 96)
n_ok, n = sum(results), len(results)
print(f"SUMMARY: {n_ok}/{n} checks passed"
      f"{'  -- ALL GOOD' if n_ok == n else '  -- INVESTIGATE'}")
print("=" * 96)
