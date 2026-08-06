"""Stage A -- verify the phase-change solver BEFORE believing anything it says.

A phase-change solver that has not been checked against an analytical benchmark
is not evidence of anything.  Five checks:

  A1  one-phase Stefan (Neumann) analytical solution, all three schemes
  A2  convergence in dx, dt and mushy-zone width
  A3  global energy closure on the ACTUAL battery configuration
  A4  latent closure: latent absorbed == melted mass * L_f
  A5  the trap: pointwise apparent heat capacity silently losing latent heat

Run:  python pcm_stage_a_verify.py
"""
import json
import numpy as np

from pcm_solver import Material, Model, stefan_front, stefan_profile, stefan_lambda
import pcm_params as P

RESULTS = {}


def banner(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)


# --------------------------------------------------------------------------- #
# Stefan benchmark harness
# --------------------------------------------------------------------------- #

SUPERHEAT = 20.0          # K above T_m at the wall
L_DOMAIN = 0.060          # m -- semi-infinite for the durations used here
T_END = 3000.0            # s


def run_stefan(N=600, dt=1.0, dTm=0.5, scheme="enthalpy_newton", t_end=T_END,
               max_iter=200):
    """One-phase Stefan: solid initially at the solidus, wall held at T_m + 20 K.

    Only the liquid conducts in the analytical problem, so the benchmark is run
    with k_liquid == k_solid and alpha_l = k / (rho cp).
    """
    pcm = Material("pcm", rho=P.PCM_RHO, cp=P.PCM_CP, k=P.PCM_K_SOLID,
                   Lf=P.PCM_LF, Tm=P.PCM_TM, dTm=dTm, k_liquid=P.PCM_K_SOLID)
    mdl = Model([(pcm, L_DOMAIN, N)], geom="planar", height=1.0, r_inner=0.0)

    Tw = P.PCM_TM + SUPERHEAT
    T0 = np.full(N, pcm.Ts)          # fully solid, at the bottom of the mush

    out = mdl.solve(T0, t_end, dt,
                    q_vol=0.0,
                    bc_inner=("dirichlet", Tw),
                    bc_outer=("adiabatic",),
                    scheme=scheme, max_iter=max_iter,
                    store_every=max(1, int(round(t_end / dt / 60))))

    alpha_l = P.PCM_K_SOLID / (P.PCM_RHO * P.PCM_CP)
    Ste = P.PCM_CP * SUPERHEAT / P.PCM_LF

    s_num = out["front"]
    s_ana = stefan_front(out["t"], Ste, alpha_l)
    ok = out["t"] > 0
    front_rel_err = np.abs(s_num[ok] - s_ana[ok]) / s_ana[ok]

    x = mdl.rc
    T_ana = stefan_profile(x, out["t"][-1], Tw, P.PCM_TM, Ste, alpha_l)
    prof_rmse = float(np.sqrt(np.mean((out["T"][-1] - T_ana) ** 2)))

    return {
        "mdl": mdl, "out": out, "Ste": Ste, "alpha_l": alpha_l,
        "s_num": s_num, "s_ana": s_ana, "t": out["t"],
        "front_rel_err_final": float(front_rel_err[-1]),
        "front_rel_err_max": float(np.max(front_rel_err)),
        "profile_rmse": prof_rmse,
        "energy_closure_rel": out["energy_closure_rel"],
        "converged": out["converged"],
    }


# --------------------------------------------------------------------------- #
# A1 -- benchmark, three schemes
# --------------------------------------------------------------------------- #

def a1_stefan_schemes():
    banner("A1  One-phase Stefan benchmark (Neumann analytical solution)")
    Ste = P.PCM_CP * SUPERHEAT / P.PCM_LF
    lam = stefan_lambda(Ste)
    print(f"  wall superheat        {SUPERHEAT:.1f} K above T_m")
    print(f"  Stefan number         {Ste:.4f}")
    print(f"  lambda (root)         {lam:.6f}")
    print(f"  analytical s(3000 s)  {stefan_front(T_END, Ste, P.PCM_K_SOLID/(P.PCM_RHO*P.PCM_CP))*1e3:.4f} mm")
    print()
    print(f"  {'scheme':18s} {'front err':>11s} {'profile RMSE':>13s} "
          f"{'energy closure':>15s}  conv")
    rows = {}
    # The chord scheme needs a high iteration cap near the liquidus (see the
    # note in pcm_solver); 600 is enough for machine tolerance here.
    caps = {"ahc_pointwise": 200, "ahc_chord": 600, "enthalpy_newton": 200}
    for scheme in ("ahc_pointwise", "ahc_chord", "enthalpy_newton"):
        r = run_stefan(N=600, dt=1.0, dTm=0.5, scheme=scheme, max_iter=caps[scheme])
        rows[scheme] = {k: v for k, v in r.items() if k not in ("mdl", "out",
                                                               "s_num", "s_ana", "t")}
        print(f"  {scheme:18s} {100*r['front_rel_err_final']:10.4f}% "
              f"{r['profile_rmse']:12.3e} K {r['energy_closure_rel']:14.3e}  "
              f"{'yes' if r['converged'] else 'NO'}")
    RESULTS["A1"] = rows
    return rows


# --------------------------------------------------------------------------- #
# A2 -- convergence
# --------------------------------------------------------------------------- #

def a2_convergence():
    banner("A2  Convergence in dx, dt and mushy-zone width")

    print("  refining mesh (dt=1 s, dTm=0.5 K, chord scheme)")
    print(f"    {'N':>6s} {'dx (mm)':>9s} {'front err':>11s} {'profile RMSE':>13s}")
    mesh = []
    for N in (150, 300, 600, 1200):
        r = run_stefan(N=N, dt=1.0, dTm=0.5)
        mesh.append({"N": N, "dx_mm": L_DOMAIN / N * 1e3,
                     "front_err": r["front_rel_err_final"],
                     "profile_rmse": r["profile_rmse"]})
        print(f"    {N:6d} {L_DOMAIN/N*1e3:9.4f} {100*r['front_rel_err_final']:10.4f}% "
              f"{r['profile_rmse']:12.3e} K")

    print("\n  refining timestep (N=600, dTm=0.5 K, chord scheme)")
    print(f"    {'dt (s)':>8s} {'front err':>11s} {'profile RMSE':>13s}")
    time = []
    for dt in (4.0, 2.0, 1.0, 0.5):
        r = run_stefan(N=600, dt=dt, dTm=0.5)
        time.append({"dt": dt, "front_err": r["front_rel_err_final"],
                     "profile_rmse": r["profile_rmse"]})
        print(f"    {dt:8.2f} {100*r['front_rel_err_final']:10.4f}% "
              f"{r['profile_rmse']:12.3e} K")

    print("\n  narrowing the mushy zone toward the sharp-interface limit")
    print("  (the analytical solution IS the dTm -> 0 limit, so this must converge)")
    print(f"    {'dTm (K)':>9s} {'front err':>11s} {'profile RMSE':>13s} {'energy':>12s}")
    mush = []
    for dTm in (4.0, 2.0, 1.0, 0.5, 0.25, 0.1):
        r = run_stefan(N=600, dt=1.0, dTm=dTm)
        mush.append({"dTm": dTm, "front_err": r["front_rel_err_final"],
                     "profile_rmse": r["profile_rmse"],
                     "energy": r["energy_closure_rel"]})
        print(f"    {dTm:9.2f} {100*r['front_rel_err_final']:10.4f}% "
              f"{r['profile_rmse']:12.3e} K {r['energy_closure_rel']:11.2e}")

    RESULTS["A2"] = {"mesh": mesh, "time": time, "mushy": mush}
    return mesh, time, mush


# --------------------------------------------------------------------------- #
# A3 / A4 -- energy and latent closure on the real configuration
# --------------------------------------------------------------------------- #

def a3_energy_closure():
    banner("A3/A4  Energy and latent closure on the battery configuration")
    mdl = P.build_model()
    T0 = np.full(mdl.N, P.T_INIT)
    t_end = 9000.0
    out = mdl.solve(T0, t_end, dt=1.0, q_vol=P.Q_NOMINAL, q_layer="core",
                    bc_inner=("adiabatic",),
                    bc_outer=("robin", P.H_CONV, P.T_AMB),
                    scheme="ahc_chord", store_every=10)

    f_end = mdl.melted_fraction(out["T_final"])
    latent_state = mdl.latent_stored(out["T_final"])
    latent_expect = f_end * mdl.pcm_mass() * P.PCM_LF

    print(f"  heat in                    {out['E_in']:12.2f} J")
    print(f"  heat out (convection)      {out['E_out']:12.2f} J")
    print(f"  enthalpy stored            {out['E_stored']:12.2f} J")
    print(f"  residual (in-out-stored)   {out['energy_residual']:12.4e} J")
    print(f"  relative closure error     {out['energy_closure_rel']:12.4e}"
          f"   <-- tolerance 1e-10")
    print()
    print(f"  melted fraction at {t_end:.0f} s   {f_end:12.4f}")
    print(f"  latent held in field       {latent_state:12.2f} J")
    print(f"  melted mass * L_f          {latent_expect:12.2f} J")
    rel = abs(latent_state - latent_expect) / max(latent_expect, 1e-12)
    print(f"  relative difference        {rel:12.4e}")
    print(f"  max Picard iterations      {out['max_picard_iters']:6d}"
          f"   converged: {out['converged']}")

    RESULTS["A3"] = {
        "E_in": out["E_in"], "E_out": out["E_out"], "E_stored": out["E_stored"],
        "residual": out["energy_residual"],
        "closure_rel": out["energy_closure_rel"],
        "melted_fraction": f_end,
        "latent_state": latent_state, "latent_expect": latent_expect,
        "latent_rel_diff": rel, "converged": out["converged"],
    }
    return out


def a4_cross_scheme():
    banner("A4b  Cross-check: chord AHC vs enthalpy-Newton on the real problem")
    mdl = P.build_model()
    T0 = np.full(mdl.N, P.T_INIT)
    res = {}
    for scheme in ("ahc_chord", "enthalpy_newton"):
        out = mdl.solve(T0, 6000.0, dt=1.0, q_vol=P.Q_NOMINAL, q_layer="core",
                        bc_outer=("robin", P.H_CONV, P.T_AMB),
                        scheme=scheme, store_every=10,
                        max_iter=600 if scheme == "ahc_chord" else 200)
        res[scheme] = out
        print(f"  {scheme:18s} T_surf(end)={out['T_surf'][-1]:8.4f} K   "
              f"melt frac={out['melted_fraction'][-1]:.6f}   "
              f"closure={out['energy_closure_rel']:.2e}")
    d = np.max(np.abs(res["ahc_chord"]["T_surf"] - res["enthalpy_newton"]["T_surf"]))
    df = np.max(np.abs(res["ahc_chord"]["melted_fraction"]
                       - res["enthalpy_newton"]["melted_fraction"]))
    print(f"\n  max |T_surf| difference between schemes  {d:.3e} K")
    print(f"  max |melted fraction| difference         {df:.3e}")
    print("  (agreement verifies the implementation, not the model form --")
    print("   the model form is what A1 checks against the analytical solution)")
    RESULTS["A4b"] = {"max_dT_surf": float(d), "max_dfrac": float(df)}
    return d, df


# --------------------------------------------------------------------------- #
# A5 -- the trap, demonstrated rather than asserted
# --------------------------------------------------------------------------- #

def a5_pointwise_trap():
    banner("A5  Trap 1: pointwise apparent heat capacity stepping over the latent plateau")
    print("  Narrow mushy zone + coarse timestep.  A cell can cross the whole")
    print("  mushy range within one step, so no iterate ever sees c_eff spike")
    print("  and the latent heat is silently skipped.  Energy closure is the")
    print("  only check that catches it -- the temperature field looks plausible.")
    print()
    print(f"  {'dTm (K)':>9s} {'dt (s)':>8s} {'scheme':>16s} {'front err':>11s} "
          f"{'energy closure':>15s}")
    rows = []
    for dTm in (2.0, 0.5, 0.1):
        for dt in (1.0, 10.0):
            for scheme in ("ahc_pointwise", "ahc_chord"):
                r = run_stefan(N=600, dt=dt, dTm=dTm, scheme=scheme)
                rows.append({"dTm": dTm, "dt": dt, "scheme": scheme,
                             "front_err": r["front_rel_err_final"],
                             "closure": r["energy_closure_rel"]})
                print(f"  {dTm:9.2f} {dt:8.1f} {scheme:>16s} "
                      f"{100*r['front_rel_err_final']:10.4f}% "
                      f"{r['energy_closure_rel']:14.3e}")
    RESULTS["A5"] = rows
    return rows


def verdict():
    banner("STAGE A VERDICT")
    a1 = RESULTS["A1"]
    a3 = RESULTS["A3"]
    checks = []

    e = a1["ahc_chord"]["front_rel_err_final"]
    checks.append(("Stefan front position, chord AHC", e < 0.01, f"{100*e:.3f}% (tol 1%)"))
    e = a1["enthalpy_newton"]["front_rel_err_final"]
    checks.append(("Stefan front position, enthalpy Newton", e < 0.01, f"{100*e:.3f}% (tol 1%)"))
    e = a3["closure_rel"]
    checks.append(("Global energy closure, battery config", e < 1e-10, f"{e:.2e} (tol 1e-10)"))
    e = a3["latent_rel_diff"]
    checks.append(("Latent == melted mass * L_f", e < 1e-10, f"{e:.2e} (tol 1e-10)"))
    e = RESULTS["A4b"]["max_dT_surf"]
    checks.append(("AHC vs enthalpy agreement", e < 1e-6, f"{e:.2e} K (tol 1e-6)"))

    allok = True
    for name, ok, detail in checks:
        allok = allok and ok
        print(f"  [{'PASS' if ok else 'FAIL'}]  {name:42s} {detail}")
    print()
    print(f"  Stage A {'PASSES -- the solver may be trusted for Stages B-E.' if allok else 'FAILS -- do not trust Stages B-E.'}")
    RESULTS["verdict"] = {"pass": bool(allok),
                          "checks": [(n, bool(o), d) for n, o, d in checks]}
    return allok


if __name__ == "__main__":
    a1_stefan_schemes()
    a2_convergence()
    a3_energy_closure()
    a4_cross_scheme()
    a5_pointwise_trap()
    ok = verdict()
    with open("results/pcm_stage_a.json", "w") as fh:
        json.dump(RESULTS, fh, indent=2, default=float)
    print("\n  wrote results/pcm_stage_a.json")
    raise SystemExit(0 if ok else 1)
