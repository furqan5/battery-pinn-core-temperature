"""Shared identifiability machinery for Stages B-E.

HOW "MELTED FRACTION" IS TREATED AS A PARAMETER
-----------------------------------------------
Melted fraction is a STATE, not a model coefficient, so it needs a
parameterisation before anything can be differentiated with respect to it.

The one used here: run a reference melt from cold, then define state(f0) as the
field on that trajectory at the instant its melted fraction equals f0.  The
estimation problem is then "given the next T_window seconds of surface
temperature, how well is f0 determined?" -- which is the actual field question.

Two consequences worth being explicit about, because they cut in opposite
directions:

  1. Along a single trajectory, f0 and the temperature field co-vary.  So
     dT_surf/df0 picks up both the direct effect of front position AND the
     indirect effect of the correlated core temperature.  That is MORE
     information than melted fraction alone would carry, so this framing is
     GENEROUS to identifiability.  A negative result under it is conservative.

  2. It follows that the traces for different f0 are time-shifted copies of one
     reference trace: T_surf(t | f0) = T_ref(t + t*(f0)).  This is not an
     artefact of the framing, it is the physics -- with everything else held
     fixed, "how much has melted" IS "how far along" -- and it means the
     sensitivity is the plateau slope multiplied by dt*/df0.  When the plateau
     is flat because the surface is pinned at the melt point, the sensitivity is
     small no matter how good the estimator.

An alternative framing would hold the core temperature fixed and move only the
front.  That state is not reachable by any heating history, so it would answer a
question no estimator ever faces.  Not used.
"""
import numpy as np

import pcm_params as P


# --------------------------------------------------------------------------- #
# Reference trajectory and state extraction
# --------------------------------------------------------------------------- #

def reference_trajectory(model=None, q=None, h=None, t_end=14000.0, dt=1.0,
                         T_init=None):
    """Full melt from cold.  Stores every step so states can be interpolated."""
    if model is None:
        model = P.build_model()
    q = P.Q_NOMINAL if q is None else q
    h = P.H_CONV if h is None else h
    T0 = np.full(model.N, P.T_INIT) if T_init is None else np.asarray(T_init)
    out = model.solve(T0, t_end, dt, q_vol=q, q_layer="core",
                      bc_inner=("adiabatic",),
                      bc_outer=("robin", h, P.T_AMB),
                      scheme="enthalpy_newton", store_every=1)
    out["model"] = model
    out["q"] = q
    out["h"] = h
    return out


def state_at_f(ref, f0):
    """Temperature field on the reference trajectory where melted fraction = f0.

    Linear interpolation in f between stored snapshots.  With dt = 1 s and a
    melt lasting thousands of seconds, a finite-difference step of df = 0.005
    spans ~20 stored points, so the piecewise-linear kinks average out and the
    central difference is a clean derivative (checked in fd_convergence).
    """
    f = ref["melted_fraction"]
    if not (0.0 < f0 < 1.0):
        raise ValueError(f"f0 must lie strictly inside (0,1), got {f0}")
    if f[-1] < f0:
        raise ValueError(f"reference trajectory only reaches f={f[-1]:.4f} < {f0}")

    # restrict to the strictly-melting portion; f is flat at 0 during warm-up
    i0 = int(np.searchsorted(f, 1e-9, side="right")) - 1
    i0 = max(i0, 0)
    i1 = int(np.searchsorted(f, 1.0 - 1e-9, side="left")) + 1
    i1 = min(i1, len(f) - 1)
    fs, ts = f[i0:i1 + 1], ref["t"][i0:i1 + 1]

    j = int(np.searchsorted(fs, f0))
    j = min(max(j, 1), len(fs) - 1)
    denom = fs[j] - fs[j - 1]
    w = 0.0 if denom <= 0 else (f0 - fs[j - 1]) / denom
    k = i0 + j
    T = (1.0 - w) * ref["T"][k - 1] + w * ref["T"][k]
    t_star = (1.0 - w) * ts[j - 1] + w * ts[j]
    return T, float(t_star)


# --------------------------------------------------------------------------- #
# Observation model
# --------------------------------------------------------------------------- #

def sensor_indices(model):
    """Outer PCM surface, and the cell can (the two-sensor contrast)."""
    return {
        "surface": model.N - 1,
        "can": model.layer_slices["can"].stop - 1,
    }


def predict_window(model, T0, q, h, t_window, dt=P.SAMPLE_DT, sensors=("surface",)):
    """Forward-simulate from state T0 and return the sampled sensor traces.

    Returns a flat vector: sensors concatenated, which is the correct stacking
    for a joint Fisher information over multiple independent sensors.
    """
    out = model.solve(T0, t_window, dt, q_vol=q, q_layer="core",
                      bc_inner=("adiabatic",),
                      bc_outer=("robin", h, P.T_AMB),
                      scheme="enthalpy_newton", store_every=1)
    idx = sensor_indices(model)
    return np.concatenate([out["T"][:, idx[s]] for s in sensors]), out


def predict_cold(model, q, h, t_end, dt=P.SAMPLE_DT, sensors=("surface",)):
    """Full melt from ambient -- the experiment for identifying PCM properties."""
    T0 = np.full(model.N, P.T_INIT)
    return predict_window(model, T0, q, h, t_end, dt, sensors)


# --------------------------------------------------------------------------- #
# Parameter handling
# --------------------------------------------------------------------------- #

# scale used to normalise each sensitivity column.  Positive-definite physical
# parameters get log scaling (scale = value); melted fraction gets scale 1
# because its natural range is [0,1].  Fixed in PCM_PREDICTIONS.md in advance.
def param_scale(name, base):
    return 1.0 if name == "f0" else abs(base[name])


DEFAULT_STEPS = {"f0": 5e-3, "q": 1e-3, "h": 1e-3, "Lf": 1e-3, "k_pcm": 1e-3}


def _eval_window(theta, ref, t_window, sensors, dt=P.SAMPLE_DT):
    """Predicted trace for the window experiment (params f0, q, h)."""
    model = ref["model"]
    T0, _ = state_at_f(ref, theta["f0"])
    y, _ = predict_window(model, T0, theta["q"], theta["h"], t_window, dt, sensors)
    return y


def _eval_cold(theta, t_end, sensors, dt=P.SAMPLE_DT):
    """Predicted trace for the cold-start experiment (params Lf, k_pcm, q, h).

    The model is rebuilt because L_f and k_pcm are material properties.
    """
    model = P.build_model(Lf=theta["Lf"], k_pcm=theta["k_pcm"],
                          k_ratio=theta.get("k_ratio", P.K_RATIO),
                          dTm=theta.get("dTm", P.PCM_DTM))
    y, _ = predict_cold(model, theta["q"], theta["h"], t_end, dt, sensors)
    return y


def sensitivity_matrix(base, free, evaluate, steps=None):
    """Central-difference sensitivities, scale-normalised.

    Differencing runs through the ACTUAL discretised solver, not an analytic
    surrogate, so the numbers describe the estimator that would really be run.
    """
    steps = dict(DEFAULT_STEPS if steps is None else steps)
    y0 = evaluate(base)
    S = np.zeros((len(y0), len(free)))
    used = {}
    for j, name in enumerate(free):
        st = steps.get(name, 1e-3)
        step = st if name == "f0" else abs(base[name]) * st
        tp, tm = dict(base), dict(base)
        tp[name] = base[name] + step
        tm[name] = base[name] - step
        S[:, j] = (evaluate(tp) - evaluate(tm)) / (2.0 * step) * param_scale(name, base)
        used[name] = step
    return S, y0, used


def crlb_from_S(S, free, base, sigma=P.NOISE_SIGMA):
    """Cramer-Rao bound from a scale-normalised sensitivity matrix.

    Conditioning is a property of S alone and is NOISE-INDEPENDENT.  The
    standard deviations scale linearly with sigma.  Both are returned, labelled.
    """
    F = S.T @ S / sigma ** 2
    cond = float(np.linalg.cond(S))
    w = np.sort(np.linalg.eigvalsh(F))[::-1]
    singular = bool(w[-1] <= 0 or (w[0] / max(w[-1], 1e-300)) > 1e14)

    if singular:
        sd_scaled = np.full(len(free), np.inf)
        corr = np.full((len(free), len(free)), np.nan)
    else:
        C = np.linalg.inv(F)
        sd_scaled = np.sqrt(np.diag(C))
        corr = C / np.outer(sd_scaled, sd_scaled)

    # back out standard deviations in natural units
    scales = np.array([param_scale(n, base) for n in free])
    sd_nat = sd_scaled * scales
    rel = np.array([sd_nat[j] / abs(base[n]) if base[n] != 0 else np.inf
                    for j, n in enumerate(free)])

    return {
        "free": list(free), "cond_scaled": cond, "singular": singular,
        "sd_scaled": sd_scaled,       # sd in scale-normalised units
        "sd_natural": sd_nat,         # sd in the parameter's own units
        "sd_relative": rel,           # sd / |base value|
        "corr": corr, "eig": w, "sigma": sigma, "S": S,
    }


def report(name, res, base, extra=""):
    tag = "RANK DEFICIENT" if res["singular"] else ""
    print(f"  {name:26s} cond(S) = {res['cond_scaled']:10.3e}  {tag}{extra}")
    for j, n in enumerate(res["free"]):
        if n == "f0":
            pct_buffer = 100.0 * res["sd_natural"][j]
            pct_value = 100.0 * res["sd_relative"][j]
            val = ("        singular" if not np.isfinite(pct_buffer)
                   else f"{pct_buffer:9.3f} % of buffer")
            extra2 = ("" if not np.isfinite(pct_value)
                      else f"   ({pct_value:8.2f} % of f={base['f0']:.2f})")
            print(f"      {n:7s} sd = {val}{extra2}")
        else:
            v = res["sd_relative"][j]
            val = "     singular" if not np.isfinite(v) else f"{100*v:9.3f} %"
            print(f"      {n:7s} sd = {val} of nominal")
    if len(res["free"]) > 1 and not res["singular"]:
        off = np.abs(res["corr"] - np.eye(len(res["free"])))
        i, j = np.unravel_index(np.argmax(off), off.shape)
        mx = off[i, j]
        flag = "   <-- near-collinear" if mx > 0.95 else ""
        print(f"      max |off-diag correlation| = {mx:.4f}"
              f"  ({res['free'][i]},{res['free'][j]}){flag}")


def print_corr(res):
    if res["singular"]:
        print("      correlation matrix: undefined (Fisher matrix singular)")
        return
    names = res["free"]
    print("      correlation matrix:")
    print("          " + "".join(f"{n:>10s}" for n in names))
    for i, n in enumerate(names):
        print(f"      {n:>8s}" + "".join(f"{res['corr'][i, j]:10.4f}"
                                         for j in range(len(names))))
