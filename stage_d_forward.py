"""Stage D -- forward gate.

Solve the FORWARD problem with KNOWN parameters and score the PINN against the
verified finite-volume solution.  Scoring is on the CENTRE-TO-SURFACE DIFFERENCE,
not just the surface trajectory: the trajectory is easy and carries almost no
spatial information, while the gradient is the entire quantity of interest.

Nothing here touches measured core data.  This is PINN-vs-FD only.
"""
import time
import numpy as np
import torch

from part7_lib import P, Record, RadialFV
from pinn_core import Scales, FieldNet, Source, pde_residual, bc_residual, predict_core_surf

torch.set_num_threads(6)          # i7-8850H: 6 physical cores, 12 logical


def run_forward(rec, R_eff, k=P.k_t, h=P.h_lit, rho_cp=None, seed=0,
                n_col=3000, adam_epochs=3000, lbfgs_steps=400,
                width=64, depth=4, n_freq=128, w_bc=10.0, verbose=True):
    if rho_cp is None:
        rho_cp = P.rho_cp()
    torch.manual_seed(seed)
    np.random.seed(seed)

    sc = Scales(rec, k, rho_cp, P.R_o)
    Bi = h * P.R_o / k
    src = Source(rec, P.V_b)

    theta0 = sc.nd_T(0.5 * (rec.T_s[0] + rec.T_c[0]))
    net = FieldNet(theta0, width, depth, n_freq)

    t_grid = torch.tensor(rec.t / rec.t[-1])
    th_inf = torch.tensor(sc.nd_T(rec.T_inf))
    theta_inf_t = (t_grid, th_inf)

    if verbose:
        print(f"    {sc}")
        print(f"    Bi = {Bi:.4f}   theta0 = {theta0:+.5f} "
              f"(= {sc.dim_T(theta0):.3f} C, measured, NOT ambient)")

    # ---- Adam with resampled collocation ---- #
    opt = torch.optim.Adam(net.parameters(), lr=1e-3)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, adam_epochs, eta_min=1e-5)
    t0 = time.time()
    for ep in range(adam_epochs):
        s = torch.rand(n_col, 1)
        t = torch.rand(n_col, 1)
        tb = torch.rand(n_col // 4, 1)
        rp = pde_residual(net, sc, src, s, t, R_eff)
        rb = bc_residual(net, sc, tb, Bi, theta_inf_t)
        loss = (rp ** 2).mean() + w_bc * (rb ** 2).mean()
        opt.zero_grad(); loss.backward(); opt.step(); sched.step()
        if verbose and (ep % 500 == 0 or ep == adam_epochs - 1):
            print(f"      adam {ep:5d}  loss={loss.item():.4e}  "
                  f"pde={(rp**2).mean().item():.3e}  bc={(rb**2).mean().item():.3e}")
    t_adam = time.time() - t0

    # ---- L-BFGS on a FROZEN collocation set ---- #
    # TRAP 5.3: if points are resampled inside the closure the objective changes
    # every call, the line search fails, and L-BFGS returns in under a second
    # REPORTING NO ERROR.  The closure counter below is the tell -- single
    # digits means it did nothing.
    s_f = torch.rand(n_col * 2, 1)
    t_f = torch.rand(n_col * 2, 1)
    tb_f = torch.rand(n_col // 2, 1)
    n_calls = [0]

    lb = torch.optim.LBFGS(net.parameters(), max_iter=lbfgs_steps,
                           tolerance_grad=1e-12, tolerance_change=1e-14,
                           history_size=50, line_search_fn="strong_wolfe")

    def closure():
        n_calls[0] += 1
        lb.zero_grad()
        rp = pde_residual(net, sc, src, s_f, t_f, R_eff)
        rb = bc_residual(net, sc, tb_f, Bi, theta_inf_t)
        loss = (rp ** 2).mean() + w_bc * (rb ** 2).mean()
        loss.backward()
        return loss

    t0 = time.time()
    lb.step(closure)
    t_lbfgs = time.time() - t0
    if verbose:
        print(f"      L-BFGS closure calls = {n_calls[0]}  "
              f"({'OK' if n_calls[0] > 20 else 'SUSPICIOUS -- see trap 5.3'})"
              f"  in {t_lbfgs:.1f} s")

    # ---- final residual on a fixed large set (truth-free selection metric) ---- #
    s_e = torch.rand(20000, 1)
    t_e = torch.rand(20000, 1)
    rp = pde_residual(net, sc, src, s_e, t_e, R_eff)
    rb = bc_residual(net, sc, torch.rand(5000, 1), Bi, theta_inf_t)
    res_metric = float((rp ** 2).mean() + w_bc * (rb ** 2).mean())

    Tc, Ts = predict_core_surf(net, sc, rec.t / rec.t[-1])
    return {"net": net, "sc": sc, "Tc": Tc, "Ts": Ts, "res": res_metric,
            "t_adam": t_adam, "t_lbfgs": t_lbfgs, "closures": n_calls[0], "Bi": Bi}


if __name__ == "__main__":
    fv = RadialFV(N=40)
    R_EFF = 0.0143398          # from Stage B [B5] on DS1 -- a KNOWN parameter here

    for tag in ("2",):
        rec = Record(tag)
        T0 = 0.5 * (rec.T_s[0] + rec.T_c[0])

        # reference solution from the verified FV solver
        q = (rec.I ** 2) * R_EFF / P.V_b
        ref = fv.solve(rec.t, q, rec.T_inf, T0, k=P.k_t, h=P.h_lit,
                       rho_cp=P.rho_cp())
        d_ref = ref["T_core"] - ref["T_surf"]

        print("=" * 92)
        print(f"STAGE D -- forward gate, dataset {tag}")
        print("=" * 92)
        print(f"  FD reference: surface rise {ref['T_surf'].max()-rec.T_inf.mean():.3f} K, "
              f"core-surf max {d_ref.max():.4f} K, mean {d_ref.mean():.4f} K")

        t0 = time.time()
        out = run_forward(rec, R_EFF, seed=0)
        print(f"  total {time.time()-t0:.1f} s  (adam {out['t_adam']:.1f} s, "
              f"lbfgs {out['t_lbfgs']:.1f} s)")

        d_pinn = out["Tc"] - out["Ts"]
        e_s = out["Ts"] - ref["T_surf"]
        e_c = out["Tc"] - ref["T_core"]
        e_d = d_pinn - d_ref

        print()
        print("  PINN vs FD (the gate):")
        print(f"    surface   RMSE {np.sqrt((e_s**2).mean()):.5f} K   "
              f"max {np.abs(e_s).max():.5f} K")
        print(f"    core      RMSE {np.sqrt((e_c**2).mean()):.5f} K   "
              f"max {np.abs(e_c).max():.5f} K")
        print(f"    core-surf RMSE {np.sqrt((e_d**2).mean()):.5f} K   "
              f"max {np.abs(e_d).max():.5f} K   <-- THE SCORED QUANTITY")
        rel = np.sqrt((e_d ** 2).mean()) / np.sqrt((d_ref ** 2).mean())
        print(f"    core-surf relative error = {100*rel:.3f} %   "
              f"(P4 threshold is 5 %)  -> {'PASS' if rel < 0.05 else 'FAIL'}")

        # physical plausibility gate (trap 5.7): ratio should be ~Bi/2
        surf_rise = out["Ts"] - rec.T_inf
        m = surf_rise > 0.5 * surf_rise.max()
        ratio = (d_pinn[m] / surf_rise[m]).mean()
        print(f"    plausibility: mean (core-surf)/(surf-rise) over the hot window "
              f"= {ratio:.4f}")
        print(f"                  Bi/2 = {out['Bi']/2:.4f}  -> ratio/(Bi/2) = "
              f"{ratio/(out['Bi']/2):.3f}  "
              f"({'within' if 1/3 < ratio/(out['Bi']/2) < 3 else 'OUTSIDE'} factor 3)")
        print(f"    (note: Bi/2 is the STEADY identity; this is a transient, so "
              f"agreement is indicative, not exact)")
