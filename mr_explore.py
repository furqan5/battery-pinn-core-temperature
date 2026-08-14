"""Multi-rate LFP: understand the records before modelling anything.

Same cell family as Part 7 (A123 ANR26650m1-b per the dataset's own
manufacturer sheet), eight C-rates, surface temperature only.

Questions this answers, in order:
  1. What is the current sign convention here?  (Part 7 taught us not to assume.)
  2. Where is the discharge segment in each record?
  3. How long does discharge last, against the ~1000 s diffusion time?
  4. Is there a rest tail (needed to see cooling)?
"""
import glob
import os
import re

import numpy as np
import pandas as pd

#: Root of the Catenaro & Onori galvanostatic discharge set
#: (Mendeley Data, doi:10.17632/kxsbr4x3j2.1). Not redistributed with this
#: repository -- see DATA.md for how to obtain it.
#:
#: Override with the BATTERY_PINN_DATA environment variable, which should point
#: at the directory containing `kxsbr4x3j2-2`. Defaults to `Data_Sets/` beside
#: this file, which is where DATA.md tells you to unpack it.
_DATA_ROOT = os.environ.get(
    "BATTERY_PINN_DATA",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "Data_Sets"),
)
BASE = os.path.join(_DATA_ROOT, "kxsbr4x3j2-2", "galvanostatic_discharge_test")

CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results",
                     "mr_cache")
os.makedirs(CACHE, exist_ok=True)


def rate_of(name):
    m = re.search(r"_(\d+(?:_\d+)?)C_", name)
    return float(m.group(1).replace("_", ".")) if m else None


def load(chem="LFP", amb="25degC"):
    """Return {C-rate: DataFrame}, cached as .npz so the xlsx read happens once."""
    out = {}
    for p in sorted(glob.glob(os.path.join(BASE, chem, amb, "*.xlsx"))):
        r = rate_of(os.path.basename(p))
        if r is None or r in out:
            continue
        cf = os.path.join(CACHE, f"{chem}_{amb}_{r}.npz")
        if os.path.exists(cf):
            d = np.load(cf)
            out[r] = pd.DataFrame({k: d[k] for k in d.files})
        else:
            df = pd.read_excel(p)[["Test_Time(s)", "Step_Index", "Voltage(V)",
                                   "Current(A)", "Surface_Temp(degC)"]]
            df.columns = ["t", "step", "V", "I", "Ts"]
            np.savez_compressed(cf, **{c: df[c].to_numpy() for c in df.columns})
            out[r] = df
    return out


if __name__ == "__main__":
    d = load()
    print(f"LFP @ 25 degC: {len(d)} C-rates -> {sorted(d)}")
    print()
    print(f"{'C':>6s} {'n':>6s} {'span s':>8s} {'step ids':>22s} "
          f"{'I range':>16s} {'corr(I,V)':>10s}")
    for r in sorted(d):
        f = d[r]
        m = f["I"].abs() > 0.05
        cc = np.corrcoef(f["I"][m], f["V"][m])[0, 1] if m.sum() > 10 else np.nan
        steps = sorted(f["step"].unique())
        print(f"{r:>6.2f} {len(f):>6d} {f['t'].max():>8.0f} "
              f"{str(steps)[:22]:>22s} "
              f"{f['I'].min():>7.2f}..{f['I'].max():<7.2f} {cc:>10.3f}")

    print()
    print("per-step breakdown for the 5C record (identify the discharge):")
    f = d[5.0]
    for s in sorted(f["step"].unique()):
        g = f[f["step"] == s]
        print(f"  step {s}: n={len(g):6d}  dur={g['t'].max()-g['t'].min():8.1f} s  "
              f"I mean={g['I'].mean():+8.3f}  V {g['V'].min():.3f}->{g['V'].max():.3f}  "
              f"Ts {g['Ts'].iloc[0]:.2f}->{g['Ts'].iloc[-1]:.2f}")
