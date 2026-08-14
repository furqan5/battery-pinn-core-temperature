"""G-1: separate irreversible from reversible heat on LCO, using multi-rate data.

THE POINT. Part 3 recovered a non-monotonic heat-generation profile on B0005 with
a minimum near x = 0.32. Whether that minimum is a resistance feature or an
artefact of absorbing the reversible term has been open, because B0005 is a
single-rate record and the two contributions are perfectly confounded in it.
MANIFEST.md states the resolution requires "an LCO entropy-profile subtraction or
a same-chemistry replication", and notes that the LFP run previously relied on
cannot test an LCO hypothesis.

The NASA archive contains LCO cells cycled at 1, 2 and 4 A within one cell at a
single ambient (B0038/39/40 at 44 C). That makes the separation possible on the
chemistry actually in question.

METHOD, and why it needs no source model at all. Bernardi with discharge-positive
current gives

    Q(t) = I (U_oc - V) - I T dU/dT  =  I^2 R(x) - I T (dU/dT)(x)

so, dividing by current,

    Q/I  =  I R(x)  -  T (dU/dT)(x)                                        (*)

At fixed depth of discharge x, (*) is a straight line in I whose SLOPE is the
ohmic coefficient and whose INTERCEPT is the reversible term. Two or more
currents therefore separate them.

Q itself is obtained from the lumped energy balance by direct differentiation,

    Q(t) = C dT/dt + K (T - T_amb),

which involves no assumed source shape whatsoever -- the quantity being
decomposed is measured, not fitted. C and K are the only inputs.

ROBUSTNESS -- AND A PREDICTION OF MINE THAT THE CONTROL REFUTED. I argued that an
error in K would enter Q as K*dT, that dT grows with current, and that it would
therefore contaminate the SLOPE more than the INTERCEPT, leaving the reversible
term well conditioned. The sensitivity sweep in g1_run.py says otherwise:
scaling K from 0.70 to 1.30 leaves the ohmic slope at 58-62 mOhm (+/-3 %) while
the reversible term moves from -0.163 to -0.416 mV/K, a factor of 2.6, roughly
proportional to K.

So the conditioning is the reverse of what I claimed. R(x) is the trustworthy
output here and dU/dT is not: it inherits the uncertainty in K almost one for
one. That matters for how the result may be quoted, and it is the reason the
linear-fit RMSE reported below UNDERSTATES the real uncertainty -- the dominant
error is systematic in K and leaves no trace in the fit residual.
"""
import glob
import os

import numpy as np
from scipy.io import loadmat
from scipy.optimize import curve_fit

from paths import NASA_ROOT as ROOT  # resolved at import, not hardcoded

# Part-2/3 locked thermal constants for this 18650 (b)
RHO_CP = 2500.0 * 1100.0
V_CELL = np.pi * 0.009 ** 2 * 0.065
C_LUMP = RHO_CP * V_CELL                    # J/K
L_SLAB = 0.018


def load_cycles(name):
    f = glob.glob(os.path.join(ROOT, "**", name + ".mat"), recursive=True)[0]
    m = loadmat(f, simplify_cells=True)
    k = [x for x in m if not x.startswith("__")][0]
    return np.atleast_1d(m[k]["cycle"])


def discharges(name, amb, cur, tol=0.35):
    out = []
    for c in load_cycles(name):
        if str(c.get("type", "")).lower() != "discharge":
            continue
        if abs(float(c.get("ambient_temperature", -999)) - amb) > 0.6:
            continue
        d = c["data"]
        I = np.abs(np.asarray(d["Current_measured"], float))
        on = I > 0.5
        if on.sum() < 25:
            continue
        Imed = float(np.median(I[on]))
        if abs(Imed - cur) > tol:
            continue
        t = np.asarray(d["Time"], float)
        T = np.asarray(d["Temperature_measured"], float)
        V = np.asarray(d["Voltage_measured"], float)
        out.append(dict(t=t - t[0], I=I, T=T, V=V, Imed=Imed, on=on))
    return out


def smooth(y, w):
    """Centred moving average; the NASA series carries no white noise, so a
    light window is enough to make differentiation stable."""
    if w < 3:
        return y.copy()
    k = np.ones(w) / w
    return np.convolve(np.pad(y, (w // 2, w // 2), mode="edge"), k, "valid")[:len(y)]


def heat_trace(rec, K, amb, w=5):
    """Q(t) from the lumped balance, and coulomb-counted x(t)."""
    t, T, I = rec["t"], rec["T"], rec["I"]
    Ts = smooth(T, w)
    dTdt = np.gradient(Ts, t)
    Q = C_LUMP * dTdt + K * (Ts - amb)
    ah = np.concatenate([[0], np.cumsum(0.5 * (I[1:] + I[:-1]) * np.diff(t))])
    x = ah / ah[-1] if ah[-1] > 0 else ah
    return x, Q


def estimate_K(recs, amb):
    """K from the post-discharge relaxation, where the source is off.

    C dT/dt = -K (T - T_amb) => exponential decay with tau = C/K. No source model
    enters, so this cannot be confounded with the heat-generation terms.
    """
    taus = []
    for r in recs:
        off = ~r["on"]
        # trailing rest only
        i = len(off)
        while i > 0 and off[i - 1]:
            i -= 1
        if len(off) - i < 12:
            continue
        t = r["t"][i:] - r["t"][i]
        y = r["T"][i:] - amb
        if y[0] < 1.0:
            continue
        m = y > 0.15 * y[0]
        if m.sum() < 8:
            continue
        try:
            p, _ = curve_fit(lambda tt, a, tau, c: a * np.exp(-tt / tau) + c,
                             t[m], y[m], p0=(y[0], 1500.0, 0.0),
                             bounds=([0, 60, -3], [80, 40000, 3]), maxfev=20000)
            taus.append(p[1])
        except Exception:
            pass
    if not taus:
        return None, 0
    return C_LUMP / float(np.median(taus)), len(taus)
