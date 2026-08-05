"""Stage C -- Cramer-Rao identifiability analysis, BEFORE any inverse fit.

Rule for this project: if a parameter set is not identifiable, do not fit it.
A fit that looks fine and means nothing is the failure mode this whole exercise
exists to avoid -- and Stage B [B3] already showed it can produce the BEST
surface fit while being 20x worse on the core.

Sensitivities are finite-differenced through the ACTUAL discretised FV model,
not an analytic surrogate, so the numbers describe the estimator we will
really run.
"""
import numpy as np
from part7_lib import P, Record, RadialFV

np.set_printoptions(suppress=True, linewidth=140)


# --------------------------------------------------------------------------- #
# Noise level -- needed to turn sensitivities into a bound
# --------------------------------------------------------------------------- #

def estimate_sigma(y):
    """Robust high-frequency noise estimate from the second difference.

    For white noise, Var(y[i+1] - 2y[i] + y[i-1]) = 6 sigma^2.  Uses the MAD so
    that genuine fast transients in the drive cycle do not inflate it.

    IMPORTANT: run this on the RAW 1.1 s trace.  Applying it after
    interpolation onto a 1 s grid correlates neighbouring samples and reports a
    noise level several times too small, which in turn makes every parameter
    set look identifiable.
    """
    d2 = y[2:] - 2 * y[1:-1] + y[:-2]
    mad = np.median(np.abs(d2 - np.median(d2)))
    return 1.4826 * mad / np.sqrt(6.0)


def raw_surface(tag):
    """Un-interpolated surface trace, for an honest instrument-noise estimate."""
    import scipy.io as sio
    from part7_lib import _RH
    return sio.loadmat(_RH + r"\Temperature_data.mat")["T_data_" + tag][:, 1]


# --------------------------------------------------------------------------- #
# Forward model in terms of a named parameter vector
# --------------------------------------------------------------------------- #

def build_predictor(rec, fv, T0, source="ohmic", n_shape=0):
    """Return (predict(theta_dict) -> T_surf, T_core)."""

    def predict(th):
        if source == "measured":
            q = rec.q_measured() / P.V_b
        else:
            R = th["R_eff"] * np.ones(len(rec.t))
            xc = rec.dod - rec.dod.mean()
            span = max(rec.dod.max() - rec.dod.min(), 1e-12)
            for j in range(1, n_shape + 1):
                R = R * (1.0 + th[f"a{j}"] * (2 * xc / span) ** j)
            q = (rec.I ** 2) * R / P.V_b
        out = fv.solve(rec.t, q, rec.T_inf, T0,
                       k=th["k"], h=th["h"], rho_cp=th["rho_cp"])
        return out["T_surf"], out["T_core"]

    return predict


def crlb(rec, fv, T0, base, free, sigma, source="ohmic", n_shape=0, rel_step=1e-4):
    """Scale-normalised CRLB for the parameters named in `free`.

    Returns dict with relative standard deviations (as fractions), the
    condition number of the scale-normalised sensitivity matrix, and the
    correlation matrix of the estimates.
    """
    predict = build_predictor(rec, fv, T0, source, n_shape)
    n = len(rec.t)
    S = np.zeros((n, len(free)))

    for j, name in enumerate(free):
        th_p, th_m = dict(base), dict(base)
        step = abs(base[name]) * rel_step
        if step == 0:
            step = rel_step
        th_p[name] = base[name] + step
        th_m[name] = base[name] - step
        yp, _ = predict(th_p)
        ym, _ = predict(th_m)
        # scale-normalised: d T / d ln(theta)  (use the step, not theta, if theta=0)
        S[:, j] = (yp - ym) / (2 * step) * (abs(base[name]) if base[name] != 0 else 1.0)

    F = S.T @ S / sigma ** 2
    cond = np.linalg.cond(S)
    # Eigen-decomposition rather than a bare inverse: a singular Fisher matrix
    # is the RESULT here, not an error, so report it instead of raising.
    w = np.linalg.eigvalsh(F)
    w = np.sort(w)[::-1]
    if w[-1] <= 0 or (w[0] / max(w[-1], 1e-300)) > 1e14:
        rel_sd = np.full(len(free), np.inf)
        corr = np.full((len(free), len(free)), np.nan)
        singular = True
    else:
        C = np.linalg.inv(F)
        rel_sd = np.sqrt(np.diag(C))
        corr = C / np.outer(rel_sd, rel_sd)
        singular = False
    return {"rel_sd": rel_sd, "cond": cond, "corr": corr, "S": S, "free": free,
            "eig": w, "singular": singular}


def report(name, res, threshold=0.25):
    worst = np.max(res["rel_sd"])
    if res["singular"]:
        verdict = "RANK DEFICIENT"
    elif worst > threshold:
        verdict = "NOT IDENTIFIABLE"
    else:
        verdict = "identifiable"
    print(f"  {name:32s} cond(S)={res['cond']:9.3e}  worst rel sd="
          f"{'    singular' if res['singular'] else f'{100*worst:9.3f} %'}"
          f"  -> {verdict}")
    for nm, s in zip(res["free"], res["rel_sd"]):
        flag = "  <-- unconstrained" if s > threshold else ""
        val = "     singular" if not np.isfinite(s) else f"{100*s:12.3f} %"
        print(f"      {nm:9s} rel sd = {val}{flag}")
    if len(res["free"]) > 1 and not res["singular"]:
        mx = np.max(np.abs(res["corr"] - np.eye(len(res["free"]))))
        print(f"      max |off-diagonal correlation| = {mx:.4f}"
              f"{'   <-- near-collinear' if mx > 0.95 else ''}")


if __name__ == "__main__":
    fv = RadialFV(N=40)

    # Stage B best-model residual RMSE, the honest model-form error scale.
    SIGMA_MODEL = {"1": 0.2200, "2": 0.2133}

    for tag in ("1", "2"):
        rec = Record(tag)
        T0 = 0.5 * (rec.T_s[0] + rec.T_c[0])
        sig_raw = estimate_sigma(raw_surface(tag))
        sig_interp = estimate_sigma(rec.T_s)
        sigma = SIGMA_MODEL[tag]

        print("=" * 96)
        print(f"DATASET {tag}")
        print("=" * 96)
        print(f"  instrument noise on RAW 1.1 s trace   sigma = {sig_raw:.5f} K")
        print(f"  same estimator after 1 Hz interpolation      = {sig_interp:.5f} K"
              f"   <-- {sig_raw/max(sig_interp,1e-12):.1f}x smaller, an artefact")
        print(f"  Stage B best-model residual RMSE             = {sigma:.4f} K")
        print()
        print(f"  CRLB is computed at sigma = {sigma:.4f} K, the MODEL-FORM error,")
        print(f"  not the instrument noise.  Rationale: Stage B residuals are ~40x")
        print(f"  the thermocouple noise and visibly structured, so model error --")
        print(f"  not sensor error -- sets how well parameters can be recovered.")
        print(f"  Even this is OPTIMISTIC: the residuals are correlated, and")
        print(f"  correlated errors carry less information than white ones.")
        print(f"  DoD window visited: {rec.dod.min():+.4f} to {rec.dod.max():+.4f}"
              f"  (span {rec.dod.max()-rec.dod.min():.4f} -- only "
              f"{100*(rec.dod.max()-rec.dod.min()):.0f}% of the x axis)")
        print()

        base = {"R_eff": 0.0145, "h": P.h_lit, "k": P.k_t, "rho_cp": P.rho_cp()}

        print("  --- thermal parameters, source = I^2 R_eff ---")
        report("{R_eff}", crlb(rec, fv, T0, base, ["R_eff"], sigma))
        report("{R_eff, h}", crlb(rec, fv, T0, base, ["R_eff", "h"], sigma))
        report("{R_eff, h, k}", crlb(rec, fv, T0, base, ["R_eff", "h", "k"], sigma))
        report("{R_eff, h, k, rho_cp}",
               crlb(rec, fv, T0, base, ["R_eff", "h", "k", "rho_cp"], sigma))
        report("{h, k, rho_cp}  (Stage B B3)",
               crlb(rec, fv, T0, base, ["h", "k", "rho_cp"], sigma))
        print()

        print("  --- P3's set: R_eff plus a SHAPE in DoD, h FREE ---")
        for order in (1, 2, 3, 4):
            b = dict(base)
            for j in range(1, order + 1):
                b[f"a{j}"] = 0.10
            free = ["R_eff"] + [f"a{j}" for j in range(1, order + 1)] + ["h"]
            report(f"{{R_eff, a1..a{order}, h}}",
                   crlb(rec, fv, T0, b, free, sigma, n_shape=order))
        print()

        print("  --- same shapes with h FIXED (the Stage E configuration) ---")
        for order in (1, 2, 3, 4):
            b = dict(base)
            for j in range(1, order + 1):
                b[f"a{j}"] = 0.10
            free = ["R_eff"] + [f"a{j}" for j in range(1, order + 1)]
            report(f"{{R_eff, a1..a{order}}}  h fixed",
                   crlb(rec, fv, T0, b, free, sigma, n_shape=order))
        print()
