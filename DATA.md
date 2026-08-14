# Data required to reproduce, and where to get it

**No measurement data is redistributed here.** Every dataset below belongs to its
original authors and must be obtained from the authoritative source. This is both
a licensing decision and a scientific one: you should get the numbers from the
people who measured them, not from a cache in someone else's repository.

Place the downloads under `Data_Sets/` at the repository root. That directory is
gitignored.

---

## 1. A123 26650 LFP cell with an internal core thermocouple

**Used by:** Part 7 core-temperature validation (`part7_lib.py`, `stage_*.py`,
`Part7_CoreValidation.ipynb`), and paper 1.

This is the dataset the whole project exists to test against — surface *and*
measured core temperature on the same cell, over two drive cycles.

> R. R. Richardson and D. A. Howey, "Sensorless battery internal temperature
> estimation using a Kalman filter with impedance measurement," *IEEE Trans.
> Sustainable Energy*, vol. 6, no. 4, pp. 1190–1199, Oct. 2015.
> doi: [10.1109/TSTE.2015.2420375](https://doi.org/10.1109/TSTE.2015.2420375)

The measurement records are distributed with the associated
`EKF-Battery-Impedance-Temperature` code release by the same authors. Obtain it
from their repository and unpack it to:

```
Data_Sets/EKF-Battery-Impedance-Temperature-master/
```

## 2. Multi-rate galvanostatic discharge, three chemistries

**Used by:** the timescale criterion (`mr_explore.py`, `mr_timescale.py`,
`mr_audit.py`, `mr_figure.py`) and the entropic-separation test
(`g1_*.py`). This supplies the eight discharge rates behind
`epsilon ~ 7.5 (t_f/tau_d)^-1.18`.

> E. Catenaro and S. Onori, "Experimental data of three lithium-ion batteries
> under galvanostatic discharge tests at different C-rates and operating
> temperatures," Mendeley Data, 2021.
> doi: [10.17632/kxsbr4x3j2.1](https://doi.org/10.17632/kxsbr4x3j2.1)

Download and unpack so the layout is:

```
Data_Sets/kxsbr4x3j2-2/galvanostatic_discharge_test/<chem>/<amb>/*.xlsx
```

`mr_explore.py` reads the `.xlsx` files once and writes an `.npz` cache to
`results/mr_cache/`. **That cache is deliberately not committed** — it contains
the raw time series (test time, step index, voltage, current, surface
temperature) rather than derived results, so committing it would redistribute the
source dataset. It regenerates automatically on first run.

The path is currently hardcoded near the top of `mr_explore.py`. Override it with
the `BATTERY_PINN_DATA` environment variable, or edit `BASE` directly.

## 3. NASA Ames battery ageing set (cell B0005)

**Used by:** the CRLB-gated heat-generation recovery over 168 discharge cycles
(`verify/Part3_*.ipynb`, `nasa_*_fits.csv`), and paper 2.

> B. Saha and K. Goebel, "Battery data set," NASA Ames Prognostics Data
> Repository, NASA Ames Research Center, Moffett Field, CA, USA, 2007.

Unpack to `Data_Sets/Li-ion Battery Dataset from NASA PCoE/`.

One provenance note that matters for reading paper 2: the `.mat` mirror carries
no `metadata.csv`, so the cycles **cannot** be stratified by ambient temperature.
The set analysed is *all 168 discharge cycles in `B0005.mat`*, not an
ambient-selected subset. The manuscript describes it that way. See G-6 in
`GAPS.md`.

---

## What *is* committed, and why that is not redistribution

| Path | Contents | Why it is safe |
|---|---|---|
| `nasa_*_fits.csv` | fitted coefficients per cycle | analysis output, not measurements |
| `results/*.json`, `results/*.log` | summary statistics, CRLB bounds, audit logs | analysis output |
| `results/*.npz` (except `mr_cache/`) | saved model fields and sensitivities | analysis output |
| `figures/*.png` | plots | analysis output |

The line drawn is: **anything that would let you skip the original measurement is
not committed; anything that is a conclusion about it is.**

## Citation obligation

If you use this code with these datasets, cite the datasets. The analysis is
worth nothing without them, and their authors did the expensive part.
