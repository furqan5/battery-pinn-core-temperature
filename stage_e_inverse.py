"""Stage E -- inverse PINN fitted to the SURFACE trace only, then Stage F scoring.

THE CORE CHANNEL IS HELD OUT ENTIRELY.  Not used for fitting, initialisation,
early stopping, or model selection.  It is read exactly once, at the end, to
score the prediction.

Leak control on the fixed parameters.  Richardson's published k = 0.404 and
h = 39.3 were identified on drive cycle 1 using BOTH thermocouples, so adopting
them would import core-derived information into a supposedly core-blind fit.
Instead h and k come from Stage B [B2], which fitted DS1's SURFACE TRACE ONLY:

    h = 37.0678 W/m^2/K     k = 0.390843 W/m/K      (core-free, other record)

and we fit DS2.  These land within 3.3% of Richardson's values, but his are not
used.

Formulation: base-subtracted split (see pinn_split.py).  The plain full-field
PINN failed the Stage D gate at 39.2%; the split passes it.

Model selection is truth-free: lowest PDE + BC residual on a large FIXED
collocation set, never the run closest to the core.
"""
import sys
import time
import numpy as np
import torch

from part7_lib import P, Record, RadialFV
from pinn_split import (LumpedBase, DeviationNet, split_pde_residual,
                        split_bc_residual, reconstruct)

torch.set_num_threads(6)

# Stage B [B2] fitted to DS1's SURFACE trace only, with the surface-only initial
# condition (an earlier version used 0.5*(T_s[0]+T_c[0]), which read the core at
# t=0; see stage_g_leak_audit.py).  Re-derived clean, these moved by 0.02% and
# 0.74% respectively -- small, but the headline uses the clean values.
H_FIXED = 37.0607         # W/m^2/K
K_FIXED = 0.393747        # W/m/K


def run_inverse(rec, h=H_FIXED, k=K_FIXED, rho_cp=None, n_shape=0, seed=0,
                n_col=2000, adam_epochs=1500, lbfgs_steps=500,
                width=48, depth=4, n_freq=0, w_bc=10.0, w_data=200.0,
                R_init=0.012, use_current=True, verbose=False):
    if rho_cp is None:
        rho_cp = P.rho_cp()
    torch.manual_seed(seed); np.random.seed(seed)

    t_ref = float(rec.t[-1])
    Fo = (k / rho_cp) * t_ref / P.R_o ** 2
    Bi = h * P.R_o / k
    # Initial condition from the SURFACE ONLY.  An earlier version used
    # 0.5*(T_s[0] + T_c[0]), which reads the measured core at t=0 -- a real, if
    # tiny, leak.  Quantified in stage_g_leak_audit.py: it moved the classical
    # core RMSE by 0.6% (0.8895 -> 0.8950 K) and R_eff by 0.02%.  Small is not
    # absent, so it is gone.  The cell is near-isothermal at t=0 anyway
    # (core-surface = 0.047 K on DS2), so a uniform T_s[0] is physically sound.
    T0 = float(rec.T_s[0])
    dT_ref = float(rec.T_s.max() - rec.T_inf.mean())

    base = LumpedBase(rec, h, rho_cp, P.R_o, P.V_b, T0, n_shape=n_shape)
    if not use_current:
        # TRAP 5.1 CONTROL: strip the I(t)^2 factor, keep its mean.  Deliberate.
        mean_I2 = float((rec.I ** 2).mean())
        rec2 = Record(rec.tag)
        rec2.I = np.full(len(rec.t), np.sqrt(mean_I2))
        base = LumpedBase(rec2, h, rho_cp, P.R_o, P.V_b, T0, n_shape=n_shape)

    net = DeviationNet(width, depth, n_freq)
    logR = torch.nn.Parameter(torch.tensor(float(np.log(R_init))))
    shp = torch.nn.Parameter(torch.zeros(n_shape)) if n_shape > 0 else None
    params = list(net.parameters()) + [logR] + ([shp] if shp is not None else [])

    t_dat = torch.tensor(rec.t / t_ref).reshape(-1, 1)
    s_dat = torch.ones_like(t_dat)
    Ts_dat = torch.tensor(rec.T_s)

    def terms(s, t, tb, idx):
        R = torch.exp(logR)
        rp = split_pde_residual(net, Fo, Bi, base, s, t, R, shp, dT_ref)
        rb = split_bc_residual(net, Bi, base, tb, R, shp, dT_ref)
        Tl = base.T_l_torch(R, shp)
        om = net(s_dat[idx], t_dat[idx]).squeeze(1)
        Ts_pred = Tl[idx] + dT_ref * om
        ld = ((Ts_pred - Ts_dat[idx]) ** 2).mean() / dT_ref ** 2
        return (rp ** 2).mean(), (rb ** 2).mean(), ld

    opt = torch.optim.Adam(params, lr=3e-3)
    sch = torch.optim.lr_scheduler.CosineAnnealingLR(opt, adam_epochs, eta_min=1e-5)
    t0 = time.time()
    n = len(t_dat)
    for ep in range(adam_epochs):
        s = torch.rand(n_col, 1); t = torch.rand(n_col, 1)
        tb = torch.rand(n_col // 4, 1)
        idx = torch.randint(0, n, (1024,))
        lp, lb_, ld = terms(s, t, tb, idx)
        L = lp + w_bc * lb_ + w_data * ld
        opt.zero_grad(); L.backward(); opt.step(); sch.step()
        if verbose and ep % 300 == 0:
            print(f"      adam {ep:5d} L={L.item():.3e} "
                  f"R={float(torch.exp(logR))*1000:.4f} mOhm")
    t_adam = time.time() - t0

    # ---- L-BFGS on a FROZEN collocation set (trap 5.3) ---- #
    s_f = torch.rand(n_col * 2, 1); t_f = torch.rand(n_col * 2, 1)
    tb_f = torch.rand(n_col // 2, 1)
    idx_f = torch.arange(n)
    calls = [0]
    lbf = torch.optim.LBFGS(params, max_iter=lbfgs_steps, tolerance_grad=1e-14,
                            tolerance_change=1e-16, history_size=50,
                            line_search_fn="strong_wolfe")

    def closure():
        calls[0] += 1
        lbf.zero_grad()
        lp, lb_, ld = terms(s_f, t_f, tb_f, idx_f)
        L = lp + w_bc * lb_ + w_data * ld
        L.backward()
        return L

    t0 = time.time(); lbf.step(closure); t_lbfgs = time.time() - t0

    # ---- truth-free selection metric, identical fixed set for every seed ---- #
    g = torch.Generator().manual_seed(12345)
    s_e = torch.rand(20000, 1, generator=g); t_e = torch.rand(20000, 1, generator=g)
    tb_e = torch.rand(5000, 1, generator=g)
    R = torch.exp(logR)
    rp = split_pde_residual(net, Fo, Bi, base, s_e, t_e, R, shp, dT_ref)
    rb = split_bc_residual(net, Bi, base, tb_e, R, shp, dT_ref)
    sel = float(((rp ** 2).mean() + w_bc * (rb ** 2).mean()).detach())

    R_eff = float(torch.exp(logR).detach())
    sh = shp.detach().numpy().copy() if shp is not None else np.array([])
    Tc, Ts = reconstruct(net, base, rec.t / t_ref, R_eff, sh, dT_ref)
    return {"R_eff": R_eff, "shape": sh, "Tc": Tc, "Ts": Ts, "sel": sel,
            "closures": calls[0], "t_adam": t_adam, "t_lbfgs": t_lbfgs,
            "surf_rmse": float(np.sqrt(np.mean((Ts - rec.T_s) ** 2))),
            "surf_max": float(np.abs(Ts - rec.T_s).max())}


def summarise(rec, runs, label):
    """Score against the measured core.  THE ONLY PLACE T_c IS READ."""
    print(f"\n  --- {label} ---")
    print(f"    {'seed':>4s} {'R_eff mOhm':>11s} {'shape':>18s} {'surfRMSE':>9s} "
          f"{'selection':>11s} {'clos':>5s} | {'CORE RMSE':>10s} {'coremax':>8s}")
    ok = []
    for i, r in enumerate(runs):
        e = r["Tc"] - rec.T_c
        r["core_rmse"] = float(np.sqrt((e ** 2).mean()))
        r["core_max"] = float(np.abs(e).max())
        conv = r["closures"] > 20 and np.isfinite(r["sel"])
        ok.append(conv)
        sh = ", ".join(f"{a:+.3f}" for a in r["shape"]) if len(r["shape"]) else "-"
        print(f"    {i:>4d} {1000*r['R_eff']:>11.4f} {sh:>18s} "
              f"{r['surf_rmse']:>9.4f} {r['sel']:>11.3e} {r['closures']:>5d} | "
              f"{r['core_rmse']:>10.4f} {r['core_max']:>8.4f}"
              f"{'' if conv else '  <-- NOT CONVERGED'}")
    R = np.array([r["R_eff"] for r in runs]) * 1000
    cr = np.array([r["core_rmse"] for r in runs])
    sr = np.array([r["surf_rmse"] for r in runs])
    print(f"    R_eff : mean {R.mean():.4f} mOhm  sd {R.std(ddof=1):.4f}  "
          f"spread {100*(R.max()-R.min())/R.mean():.3f} %")
    print(f"    surf  : mean {sr.mean():.4f} K  sd {sr.std(ddof=1):.4f}")
    print(f"    core  : mean {cr.mean():.4f} K  sd {cr.std(ddof=1):.4f}  "
          f"min {cr.min():.4f}  max {cr.max():.4f}")
    print(f"    non-convergence: {100*(1-np.mean(ok)):.1f} % "
          f"({len(ok)-sum(ok)}/{len(ok)})")
    b = int(np.argmin([r["sel"] for r in runs]))
    print(f"    SELECTED by lowest PDE+BC residual (truth-free): seed {b}")
    print(f"      -> core RMSE {runs[b]['core_rmse']:.4f} K, "
          f"max {runs[b]['core_max']:.4f} K, surface RMSE {runs[b]['surf_rmse']:.4f} K")
    bt = int(np.argmin([r["core_rmse"] for r in runs]))
    print(f"    (reference only, NOT used for selection: best-by-core is seed {bt} "
          f"at {runs[bt]['core_rmse']:.4f} K)")
    return runs[b], ok


if __name__ == "__main__":
    n_seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    rec = Record("2")
    print("=" * 100)
    print(f"STAGE E -- inverse PINN, DATASET 2, SURFACE ONLY, {n_seeds} seeds")
    print("=" * 100)
    print(f"  h = {H_FIXED:.4f} W/m2/K, k = {K_FIXED:.6f} W/m/K  "
          f"(Stage B [B2] on DS1 SURFACE ONLY)")
    print(f"  Bi = {H_FIXED*P.R_o/K_FIXED:.4f}   Bi/2 = {H_FIXED*P.R_o/K_FIXED/2:.4f}")
    print(f"  measured core-surf max on DS2 = {(rec.T_c-rec.T_s).max():.4f} K")

    allruns = {}
    for n_shape, lab in ((0, "order 0: scalar R_eff"),
                         (1, "order 1: R_eff (1 + a1 x)")):
        runs = []
        for sd in range(n_seeds):
            t0 = time.time()
            r = run_inverse(rec, n_shape=n_shape, seed=sd)
            r["wall"] = time.time() - t0
            runs.append(r)
            print(f"    seed {sd}: {r['wall']:.0f} s, {r['closures']} closures, "
                  f"R_eff {1000*r['R_eff']:.4f} mOhm, surf {r['surf_rmse']:.4f} K")
        best, ok = summarise(rec, runs, lab)
        allruns[n_shape] = runs
        np.savez(f"results/stage_e_shape{n_shape}.npz",
                 Tc=np.array([r["Tc"] for r in runs]),
                 Ts=np.array([r["Ts"] for r in runs]),
                 R_eff=np.array([r["R_eff"] for r in runs]),
                 sel=np.array([r["sel"] for r in runs]),
                 core_rmse=np.array([r["core_rmse"] for r in runs]),
                 surf_rmse=np.array([r["surf_rmse"] for r in runs]),
                 T_c_meas=rec.T_c, T_s_meas=rec.T_s, t=rec.t)

    # ---- trap 5.1 control: same machinery, current factor removed ---- #
    print("\n" + "=" * 100)
    print("  CONTROL (trap 5.1): identical fit with the I(t)^2 factor REMOVED")
    print("=" * 100)
    rc = run_inverse(rec, n_shape=0, seed=0, use_current=False)
    e = rc["Tc"] - rec.T_c
    print(f"    surface RMSE {rc['surf_rmse']:.4f} K   "
          f"CORE RMSE {np.sqrt((e**2).mean()):.4f} K")
    good = allruns[0][0]
    print(f"    vs with I^2: surface {good['surf_rmse']:.4f} K, "
          f"core {good['core_rmse']:.4f} K")
    print(f"    -> removing the current factor costs "
          f"{rc['surf_rmse']/good['surf_rmse']:.1f}x on the surface and "
          f"{np.sqrt((e**2).mean())/good['core_rmse']:.1f}x on the core")
