"""G-1 step 2: are the multi-current LCO cells usable for an entropic separation?

Wanted per cell: several discharges at each of 1, 2 and 4 A, at ONE ambient, each
carrying enough temperature rise to invert. Mixed ambients confound the entropic
term because dU/dT is multiplied by absolute temperature.
"""
import glob
import os

import numpy as np
from scipy.io import loadmat

from paths import NASA_ROOT as ROOT  # resolved at import, not hardcoded

CELLS = ["B0038", "B0039", "B0040", "B0042", "B0043", "B0044"]


def load(name):
    f = glob.glob(os.path.join(ROOT, "**", name + ".mat"), recursive=True)[0]
    m = loadmat(f, simplify_cells=True)
    k = [x for x in m if not x.startswith("__")][0]
    return np.atleast_1d(m[k]["cycle"])


for name in CELLS:
    cyc = load(name)
    rows = []
    for i, c in enumerate(cyc):
        if str(c.get("type", "")).lower() != "discharge":
            continue
        d = c["data"]
        I = np.abs(np.asarray(d["Current_measured"], float))
        T = np.asarray(d["Temperature_measured"], float)
        t = np.asarray(d["Time"], float)
        on = I > 0.5
        if on.sum() < 20:
            continue
        rows.append(dict(idx=i, I=round(float(np.median(I[on])), 1),
                         amb=float(c.get("ambient_temperature", np.nan)),
                         rise=float(T.max() - T.min()),
                         dur=float(t[-1] - t[0]),
                         ah=float(np.trapezoid(I[on], t[on]) / 3600.0),
                         n=len(t)))
    print("=" * 92)
    print(f"{name}: {len(rows)} discharges")
    print("=" * 92)
    combos = {}
    for r in rows:
        combos.setdefault((r["amb"], r["I"]), []).append(r)
    print(f"  {'T_amb':>7s} {'I (A)':>6s} {'n cyc':>6s} {'rise K':>8s} "
          f"{'dur s':>8s} {'Ah':>6s} {'samples':>8s}")
    for (amb, I), g in sorted(combos.items()):
        print(f"  {amb:>7.1f} {I:>6.1f} {len(g):>6d} "
              f"{np.median([x['rise'] for x in g]):>8.2f} "
              f"{np.median([x['dur'] for x in g]):>8.0f} "
              f"{np.median([x['ah'] for x in g]):>6.2f} "
              f"{int(np.median([x['n'] for x in g])):>8d}")
    ambs = {a for a, _ in combos}
    for a in sorted(ambs):
        Is = sorted({I for aa, I in combos if aa == a})
        if len(Is) >= 3:
            n = min(len(combos[(a, I)]) for I in Is)
            print(f"  --> USABLE at {a:.0f} C: currents {Is}, "
                  f"min {n} cycles per current")
    print()
