"""Stage B -- classical baseline fitted to the SURFACE trace only.

Purpose (per the project's hard-won lesson): the classical fitter is the
instrument that detects PINN bugs.  Comparing the PINN's recovered source
against this one is what catches a missing I(t)^2 factor, which no loss curve
will reveal.

The core channel is NOT used anywhere here except to report error at the very
end.  Nothing is fitted, initialised, or selected using it.
"""
import time
import numpy as np
from scipy.optimize import least_squares, minimize_scalar

from part7_lib import P, Record, RadialFV

np.set_printoptions(suppress=True)

# Parameter bounds.  Unbounded fits on this problem ran h to 2.5e8 W/m^2/K when
# fed a bad source -- bounds turn a silent divergence into a visible pin.
BOUNDS = {
    "k":      (0.05, 5.0),        # W/m/K   radial, jellyroll across layers
    "h":      (1.0, 500.0),       # W/m^2/K natural conv ~5-25, forced/liquid higher
    "rho_cp": (0.5e6, 6.0e6),     # J/m^3/K
    "R_eff":  (1e-4, 0.2),        # Ohm
    "R_slope": (-5.0, 5.0),       # - per unit DoD
}


def fit(rec, fv, free, fixed, q_kind, T0, verbose=True):
    """Bounded least-squares fit of `free` parameters to the SURFACE trace."""
    names = list(free)
    x0 = np.array([free[n] for n in names], float)
    lo = np.array([BOUNDS[n][0] for n in names], float)
    hi = np.array([BOUNDS[n][1] for n in names], float)

    def params_of(x):
        p = dict(fixed)
        p.update(dict(zip(names, x)))
        return p

    def q_of(p):
        if q_kind == "measured":
            return rec.q_measured() / P.V_b
        R = p["R_eff"]
        if "R_slope" in p:
            R = R * (1.0 + p["R_slope"] * (rec.dod - rec.dod.mean()))
        return (rec.I ** 2) * R / P.V_b

    def predict(x):
        p = params_of(x)
        out = fv.solve(rec.t, q_of(p), rec.T_inf, T0,
                       k=p["k"], h=p["h"], rho_cp=p["rho_cp"])
        return out["T_surf"], out["T_core"]

    def resid(x):
        return predict(x)[0] - rec.T_s

    t0 = time.time()
    sol = least_squares(resid, x0, bounds=(lo, hi), method="trf",
                        x_scale="jac", ftol=1e-12, xtol=1e-12)
    el = time.time() - t0

    p = params_of(sol.x)
    Ts_p, Tc_p = predict(sol.x)
    rmse = float(np.sqrt(np.mean((Ts_p - rec.T_s) ** 2)))
    mx = float(np.abs(Ts_p - rec.T_s).max())

    J = sol.jac
    dof = max(len(rec.t) - len(names), 1)
    s2 = float(2 * sol.cost / dof)
    try:
        sd = np.sqrt(np.diag(np.linalg.inv(J.T @ J) * s2))
    except np.linalg.LinAlgError:
        sd = np.full(len(names), np.nan)
    # scale-normalised condition number of the sensitivity matrix
    Js = J * np.abs(sol.x)
    cond = np.linalg.cond(Js) if Js.shape[1] > 1 else 1.0

    if verbose:
        pinned = [n for n, v in zip(names, sol.x)
                  if abs(v - BOUNDS[n][0]) < 1e-8 or abs(v - BOUNDS[n][1]) < 1e-8]
        print(f"    {el:5.1f} s, {sol.nfev:3d} evals | surface RMSE {rmse:.4f} K, "
              f"max {mx:.4f} K | cond(J_scaled) = {cond:.3e}")
        for n, v, s in zip(names, sol.x, sd):
            print(f"      {n:8s} = {v:12.6g}  +/- {s:10.4g}  ({100*s/abs(v):8.2f} %)")
        if pinned:
            print(f"      !! PINNED AT BOUND: {pinned}  -- not an estimate")

    return {"p": p, "x": sol.x, "names": names, "Ts": Ts_p, "Tc": Tc_p,
            "rmse": rmse, "max": mx, "sd": dict(zip(names, sd)), "cond": cond}


def core_score(rec, r):
    e = r["Tc"] - rec.T_c
    return float(np.sqrt((e ** 2).mean())), float(np.abs(e).max())


if __name__ == "__main__":
    fv = RadialFV(N=40)
    store = {}

    for tag in ("1", "2"):
        rec = Record(tag)
        # Initial condition: the cell is essentially isothermal at t=0 but it is
        # NOT at ambient.  Using ambient here is the classic silent error (5.2).
        T0 = float(rec.T_s[0])   # surface only: no core leak

        print("=" * 94)
        print(f"DATASET {tag}   n={len(rec.t)}  duration={rec.t[-1]:.0f} s")
        print("=" * 94)
        print(f"  T_s(0)={rec.T_s[0]:.3f} C  T_c(0)={rec.T_c[0]:.3f} C  "
              f"T_inf(0)={rec.T_inf[0]:.3f} C")
        print(f"  uniform T0 = {T0:.4f} C, sitting {T0-rec.T_inf[0]:+.3f} K above "
              f"ambient  <-- printed per trap 5.2")
        print(f"  fitted plateau OCV U = {rec.U_ocv_reg:.4f} V, "
              f"ohmic regression slope = {1000*rec.R_ohmic_reg:.3f} mOhm")

        Q = rec.q_measured()
        print(f"  measured Q = I(V-U): mean {Q.mean():.4f} W, max {Q.max():.3f} W, "
              f"negative fraction {np.mean(Q < 0):.3f}")

        t0 = time.time(); fv.solve(rec.t, Q / P.V_b, rec.T_inf, T0)
        print(f"  one FV solve (N=40, {len(rec.t)} steps): "
              f"{1000*(time.time()-t0):.0f} ms")
        print()

        print("  [B1] source = MEASURED I(V-U) | fit {h}, k & rho_cp at literature")
        r1 = fit(rec, fv, {"h": 39.3}, {"k": P.k_t, "rho_cp": P.rho_cp()},
                 "measured", T0)
        print("  [B2] source = MEASURED I(V-U) | fit {h, k}")
        r2 = fit(rec, fv, {"h": 39.3, "k": 0.404}, {"rho_cp": P.rho_cp()},
                 "measured", T0)
        print("  [B3] source = MEASURED I(V-U) | fit {h, k, rho_cp}  <- expect trouble")
        r3 = fit(rec, fv, {"h": 39.3, "k": 0.404, "rho_cp": P.rho_cp()}, {},
                 "measured", T0)
        print("  [B4] source = I^2 R_eff       | fit {R_eff, h}, k & rho_cp fixed")
        r4 = fit(rec, fv, {"R_eff": 0.012, "h": 39.3},
                 {"k": P.k_t, "rho_cp": P.rho_cp()}, "ohmic", T0)
        print("  [B5] source = I^2 R_eff       | fit {R_eff} only, h at literature")
        r5 = fit(rec, fv, {"R_eff": 0.012},
                 {"k": P.k_t, "h": P.h_lit, "rho_cp": P.rho_cp()}, "ohmic", T0)
        print("  [B6] source = I^2 R_eff(DoD)  | fit {R_eff, R_slope}, h at literature")
        r6 = fit(rec, fv, {"R_eff": 0.012, "R_slope": 0.0},
                 {"k": P.k_t, "h": P.h_lit, "rho_cp": P.rho_cp()}, "ohmic", T0)

        # ---- trap 5.1 control: same model with the current factor removed ---- #
        print("  [B7] CONTROL: source with NO I(t) factor at all (constant power)")

        def pred_const(Pw):
            q = np.full(len(rec.t), Pw) / P.V_b
            return fv.solve(rec.t, q, rec.T_inf, T0, k=P.k_t, h=P.h_lit,
                            rho_cp=P.rho_cp())
        rr = minimize_scalar(
            lambda Pw: float(np.sum((pred_const(Pw)["T_surf"] - rec.T_s) ** 2)),
            bounds=(1e-4, 50.0), method="bounded")
        oc = pred_const(rr.x)
        rmse_c = float(np.sqrt(np.mean((oc["T_surf"] - rec.T_s) ** 2)))
        ec = oc["T_core"] - rec.T_c
        print(f"      best constant power = {rr.x:.4f} W | surface RMSE {rmse_c:.4f} K "
              f"({rmse_c/r5['rmse']:.1f}x worse than B5) | "
              f"CORE RMSE {np.sqrt((ec**2).mean()):.4f} K")

        print()
        print("  ---- CORE channel: scoring only, never used in any fit above ----")
        rows = [("B1", r1), ("B2", r2), ("B3", r3), ("B4", r4), ("B5", r5), ("B6", r6)]
        print(f"    {'fit':4s} {'surf RMSE':>10s} {'core RMSE':>10s} {'core max':>9s}"
              f"  {'ratio core/surf':>15s}")
        for nm, r in rows:
            cr, cm = core_score(rec, r)
            print(f"    {nm:4s} {r['rmse']:10.4f} {cr:10.4f} {cm:9.4f}"
                  f"  {cr/r['rmse']:15.2f}")
        store[tag] = {nm: r for nm, r in rows}
        print()

    # ---- cross-record: does h identified on DS1 transfer to DS2? ---- #
    print("=" * 94)
    print("CROSS-RECORD CHECK (this is what makes h independent for Stage E)")
    print("=" * 94)
    for nm in ("B1", "B2", "B4", "B5"):
        h1 = store["1"][nm]["p"]["h"]
        h2 = store["2"][nm]["p"]["h"]
        print(f"  {nm}: h(DS1) = {h1:8.3f}   h(DS2) = {h2:8.3f}   "
              f"difference {100*abs(h1-h2)/h1:6.2f} %")
    for nm in ("B4", "B5", "B6"):
        R1 = store["1"][nm]["p"]["R_eff"]
        R2 = store["2"][nm]["p"]["R_eff"]
        print(f"  {nm}: R_eff(DS1) = {1000*R1:7.3f} mOhm  R_eff(DS2) = {1000*R2:7.3f} mOhm"
              f"   difference {100*abs(R1-R2)/R1:6.2f} %")
