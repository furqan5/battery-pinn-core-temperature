"""Stage D (second attempt) -- forward gate using the base-subtracted formulation.

The first attempt (stage_d_forward.py, full-field PINN with Fourier features)
FAILED at 39.2% core-surface error.  Two causes, both now measured rather than
guessed:

  1. The drive-cycle source is broadband -- 50% of its energy above harmonic 233
     (periods under 15 s).  A Fourier-feature MLP would need n_freq ~ 1500.
  2. Fourier features actively HURT.  On a constant-source control the
     core-surface error was 0.086% at n_freq=0 and 8.0% at n_freq=64,
     monotonically worse with more features.

Fix: split off the lumped solution exactly (the spiky source cancels), and drop
the Fourier features.  Scored, as before, on the CENTRE-TO-SURFACE DIFFERENCE
against the verified finite-volume solver.  No measured core data is used.
"""
import time
import numpy as np
import torch

from part7_lib import P, Record, RadialFV
from pinn_split import (LumpedBase, DeviationNet, split_pde_residual,
                        split_bc_residual, reconstruct)

torch.set_num_threads(6)


def run_forward_split(rec, R_eff, k=P.k_t, h=P.h_lit, rho_cp=None, shape=(),
                      seed=0, n_col=2000, adam_epochs=1500, lbfgs_steps=500,
                      width=48, depth=4, n_freq=0, w_bc=10.0, verbose=True):
    if rho_cp is None:
        rho_cp = P.rho_cp()
    torch.manual_seed(seed); np.random.seed(seed)

    t_ref = float(rec.t[-1])
    Fo = (k / rho_cp) * t_ref / P.R_o ** 2
    Bi = h * P.R_o / k
    T0 = 0.5 * (rec.T_s[0] + rec.T_c[0])
    dT_ref = float(rec.T_s.max() - rec.T_inf.mean())

    base = LumpedBase(rec, h, rho_cp, P.R_o, P.V_b, T0, n_shape=len(shape))
    net = DeviationNet(width, depth, n_freq)

    if verbose:
        print(f"    Fo={Fo:.4f}  Bi={Bi:.4f}  dT_ref={dT_ref:.4f} K  "
              f"tau_lumped={base.tau:.1f} s")
        print(f"    T_l(0) = {T0:.4f} C  (measured initial temperature, "
              f"{T0-rec.T_inf[0]:+.3f} K above ambient -- trap 5.2 carried here)")

    sh = torch.tensor(shape) if len(shape) else None
    opt = torch.optim.Adam(net.parameters(), lr=2e-3)
    sch = torch.optim.lr_scheduler.CosineAnnealingLR(opt, adam_epochs, eta_min=1e-5)
    t0 = time.time()
    for ep in range(adam_epochs):
        s = torch.rand(n_col, 1); t = torch.rand(n_col, 1)
        tb = torch.rand(n_col // 4, 1)
        rp = split_pde_residual(net, Fo, Bi, base, s, t, R_eff, sh, dT_ref)
        rb = split_bc_residual(net, Bi, base, tb, R_eff, sh, dT_ref)
        L = (rp ** 2).mean() + w_bc * (rb ** 2).mean()
        opt.zero_grad(); L.backward(); opt.step(); sch.step()
        if verbose and ep % 300 == 0:
            print(f"      adam {ep:5d} loss={L.item():.4e}")
    t_adam = time.time() - t0

    # L-BFGS on a FROZEN collocation set (trap 5.3)
    s_f = torch.rand(n_col * 2, 1); t_f = torch.rand(n_col * 2, 1)
    tb_f = torch.rand(n_col // 2, 1)
    calls = [0]
    lb = torch.optim.LBFGS(net.parameters(), max_iter=lbfgs_steps,
                           tolerance_grad=1e-14, tolerance_change=1e-16,
                           history_size=50, line_search_fn="strong_wolfe")

    def closure():
        calls[0] += 1
        lb.zero_grad()
        rp = split_pde_residual(net, Fo, Bi, base, s_f, t_f, R_eff, sh, dT_ref)
        rb = split_bc_residual(net, Bi, base, tb_f, R_eff, sh, dT_ref)
        L = (rp ** 2).mean() + w_bc * (rb ** 2).mean()
        L.backward()
        return L

    t0 = time.time(); lb.step(closure); t_lbfgs = time.time() - t0
    if verbose:
        print(f"      L-BFGS closures = {calls[0]} "
              f"({'OK' if calls[0] > 20 else 'SUSPICIOUS -- trap 5.3'}) "
              f"in {t_lbfgs:.0f} s")

    g = torch.Generator().manual_seed(999)
    s_e = torch.rand(20000, 1, generator=g); t_e = torch.rand(20000, 1, generator=g)
    tb_e = torch.rand(5000, 1, generator=g)
    rp = split_pde_residual(net, Fo, Bi, base, s_e, t_e, R_eff, sh, dT_ref)
    rb = split_bc_residual(net, Bi, base, tb_e, R_eff, sh, dT_ref)
    sel = float(((rp ** 2).mean() + w_bc * (rb ** 2).mean()).detach())

    Tc, Ts = reconstruct(net, base, rec.t / t_ref, R_eff, shape, dT_ref)
    return {"Tc": Tc, "Ts": Ts, "sel": sel, "closures": calls[0], "Bi": Bi,
            "t_adam": t_adam, "t_lbfgs": t_lbfgs, "base": base, "net": net}


if __name__ == "__main__":
    fv = RadialFV(N=40)
    R_EFF = 0.0143398

    for tag in ("2", "1"):
        rec = Record(tag)
        T0 = 0.5 * (rec.T_s[0] + rec.T_c[0])
        q = (rec.I ** 2) * R_EFF / P.V_b
        ref = fv.solve(rec.t, q, rec.T_inf, T0, k=P.k_t, h=P.h_lit,
                       rho_cp=P.rho_cp())
        d_ref = ref["T_core"] - ref["T_surf"]

        print("=" * 92)
        print(f"STAGE D (split) -- forward gate, dataset {tag}")
        print("=" * 92)
        print(f"  FD reference: core-surf max {d_ref.max():.4f} K, "
              f"rms {np.sqrt((d_ref**2).mean()):.4f} K")

        t0 = time.time()
        out = run_forward_split(rec, R_EFF, seed=0)
        print(f"  total {time.time()-t0:.0f} s")

        d = out["Tc"] - out["Ts"]
        e_s, e_c, e_d = out["Ts"] - ref["T_surf"], out["Tc"] - ref["T_core"], d - d_ref
        rel = np.sqrt((e_d ** 2).mean()) / np.sqrt((d_ref ** 2).mean())
        print(f"  PINN vs FD:")
        print(f"    surface   RMSE {np.sqrt((e_s**2).mean()):.5f} K  "
              f"max {np.abs(e_s).max():.5f} K")
        print(f"    core      RMSE {np.sqrt((e_c**2).mean()):.5f} K  "
              f"max {np.abs(e_c).max():.5f} K")
        print(f"    core-surf RMSE {np.sqrt((e_d**2).mean()):.5f} K  "
              f"max {np.abs(e_d).max():.5f} K   <-- SCORED")
        print(f"    relative core-surf error = {100*rel:.4f} %  "
              f"(P4 threshold 5 %)  -> {'PASS' if rel < 0.05 else 'FAIL'}")
        m = (out["Ts"] - rec.T_inf) > 0.5 * (out["Ts"] - rec.T_inf).max()
        ratio = float((d[m] / (out["Ts"] - rec.T_inf)[m]).mean())
        print(f"    plausibility ratio {ratio:.4f} vs Bi/2 {out['Bi']/2:.4f} "
              f"-> factor {ratio/(out['Bi']/2):.3f} "
              f"({'within' if 1/3 < ratio/(out['Bi']/2) < 3 else 'OUTSIDE'} 3x)")
        print()
