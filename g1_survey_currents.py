"""G-1 step 1: does the NASA archive contain LCO cells discharged at DIFFERENT currents?

The entropic separation needs q/I plotted against I at fixed depth of discharge:
the slope is the ohmic coefficient, the intercept is -T dU/dT. That requires two
or more discharge currents on the SAME chemistry. B0005 is all at 2 A, which is
why Part 3 had to absorb the reversible term into effective coefficients.

If other NASA cells ran at other currents, the separation becomes possible on
LCO -- the chemistry whose interior minimum is actually in question. If not, G-1
cannot be closed with this archive and the honest answer is to say so.
"""
import glob
import os

import numpy as np
from scipy.io import loadmat

from paths import NASA_ROOT as ROOT  # resolved at import, not hardcoded

files = sorted(glob.glob(os.path.join(ROOT, "**", "B*.mat"), recursive=True))
print(f"{len(files)} .mat files found\n")
print(f"{'cell':8s} {'n_dis':>6s} {'I mean':>8s} {'I set':>26s} "
      f"{'T_amb set':>14s} {'rise K':>8s}")

summary = {}
for f in files:
    name = os.path.basename(f).replace(".mat", "")
    try:
        m = loadmat(f, simplify_cells=True)
    except Exception as e:
        print(f"{name:8s} load failed: {type(e).__name__}")
        continue
    key = [k for k in m if not k.startswith("__")]
    if not key:
        continue
    cyc = m[key[0]]["cycle"]
    Is, Ts, rises = [], [], []
    for c in np.atleast_1d(cyc):
        if str(c.get("type", "")).lower() != "discharge":
            continue
        d = c["data"]
        I = np.abs(np.asarray(d["Current_measured"], float))
        T = np.asarray(d["Temperature_measured"], float)
        on = I > 0.5
        if on.sum() < 20:
            continue
        Is.append(float(np.median(I[on])))
        Ts.append(float(c.get("ambient_temperature", np.nan)))
        rises.append(float(T.max() - T.min()))
    if not Is:
        continue
    Iset = sorted({round(v, 1) for v in Is})
    Tset = sorted({t for t in Ts if np.isfinite(t)})
    summary[name] = dict(n=len(Is), I=Iset, T=Tset,
                         rise=float(np.median(rises)))
    print(f"{name:8s} {len(Is):>6d} {np.mean(Is):>8.2f} "
          f"{str(Iset)[:26]:>26s} {str(Tset)[:14]:>14s} {np.median(rises):>8.2f}")

allI = sorted({i for s in summary.values() for i in s["I"]})
print(f"\ndistinct discharge currents across the archive: {allI}")
if len(allI) > 1:
    print("MULTIPLE CURRENTS PRESENT -> an LCO separation is possible in principle.")
    print("Caveat to check next: are they on the same cell, or only across cells?")
    for name, s in summary.items():
        if len(s["I"]) > 1:
            print(f"   {name}: MULTIPLE currents within one cell -> {s['I']}")
else:
    print("SINGLE current only -> the separation cannot be done on this archive.")
