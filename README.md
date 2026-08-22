# Core temperature from the outside of a lithium-ion cell — and when the machinery is worth it

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22043513.svg)](https://doi.org/10.5281/zenodo.22043513)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg)](https://www.python.org/)
[![Tests 24/24](https://img.shields.io/badge/verification-24%2F24%20passing-brightgreen.svg)](tests_fd.py)

Analysis code and results for two manuscripts on inverse estimation of internal
cell state from a single external temperature channel.

**The short version of what is in here: the internal temperature of a
cylindrical LFP cell is recoverable from surface data to 0.14 K, verified against
a real core thermocouple — and the inverse physics-informed neural network built
to do it is the worst of the three methods tried. One line of algebra beats it by
4.4×.** That negative result is the most useful thing in the repository and it is
reported as the headline rather than a footnote.

![Internal temperature reconstructed from surface data alone. Black is the
measured core thermocouple, held out from every fit. Blue is the measured
surface — the only input. Red is the neural estimator's prediction.
](figures/core_validation.png)

*The blue trace is all the estimator sees. The black trace is what it has to
predict, and it never sees it. The gap between them peaks at 6.54 K. Shown here
is the neural estimator at 0.614 K RMSE; the one-line algebraic relation reaches
0.141 K on the same data, which is the finding.*

---

## The two manuscripts

| | Subject |
|---|---|
| **Paper 1** | Core temperature estimation validated against an internal thermocouple, and a timescale criterion for when a transient inverse solver is warranted |
| **Paper 2** | Identifiability-gated inverse estimation of heat generation from one surface channel — Cramér–Rao bounds computed *before* fitting |

**The manuscripts themselves are not in this repository — this is a code
release.** Everything needed to reproduce every number they report is here.

---

## Headline results

Fitting the **surface trace only**, scored against a measured core thermocouple
that no estimator sees, on two independent drive cycles:

| method | DS2 core RMSE | DS1 core RMSE |
|---|---|---|
| **quasi-steady one-liner** `T_core = T_surf + (Bi/2)(T_surf − T_inf)` | **0.1409 K** | **0.3366 K** |
| classical transient finite volume | 0.8943 K | 1.0544 K |
| inverse PINN (the method under test) | 0.6135 K | 0.5024 K |

against measured core-to-surface gradients peaking at 6.54 K and 7.13 K.

The reason the one-liner wins is quasi-steadiness, and it is demonstrated rather
than asserted: the ratio (core − surface)/(surface − ambient) is constant to
**2.0 %** across a 3541 s drive cycle with 30 A peaks. One number therefore
carries the entire spatial structure, and there is nothing dynamic left for a PDE
solver to add.

Two qualifications that belong beside the headline:

- **The PINN does not solve the PDE it claims to.** Its field departs from a true
  solution by 0.80 K on the surface — larger than its own core error. Its
  recovered parameter, put through an exact solver, gives 1.52 K, *worse* than the
  classical fit.
- **Its target hit depends on a hyperparameter that truth-free selection gets
  wrong.** Choosing the data weight by lowest PDE residual — the criterion this
  project mandates — gives 1.05 K and misses the pre-registered < 1 K target.

Full accounting, including the six scored predictions (four hits, two misses) and
three mid-session predictions that were wrong, is in **`FINDINGS.md`**.

The companion PCM study reaches the same shape of conclusion by a different
route: the remaining latent buffer is identifiable to 0.1 % of capacity from one
surface sensor, but the **melt-front position is not identifiable at any heating
rate tested** across a 40× sweep, because 99 % of the Fisher information sits in
the mean level of the trace. See **`PCM_FINDINGS.md`**.

---

## Reproducing a specific claim

Data first — see **`DATA.md`**. No measurement data is redistributed here.

```bash
pip install -r requirements.txt
```

| Claim | Run this |
|---|---|
| Cylindrical FV solver correct to ~1e-13 against the analytic steady state | `python tests_fd.py` |
| Split/lumped-base machinery, 9 checks | `python tests_split.py` |
| Classical surface-only baseline | `python stage_b_classical.py` |
| CRLB identifiability, computed before any inverse fit | `python stage_c_crlb.py` |
| Forward PINN gate | `python stage_d_forward_split.py` |
| Inverse PINN, headline run | `python stage_e_inverse.py` |
| DS1 replication with record roles swapped | `python stage_e_replicate_ds1.py` |
| Core-channel leak audit | `python stage_g_leak_audit.py` |
| Quasi-steady one-liner (the result that wins) | `python stage_g_quasisteady.py` |
| Data-weight sweep showing the P5 hit is fragile | `python stage_g_weight_sweep.py` |
| Timescale criterion across eight discharge rates | `python mr_timescale.py` |
| Entropic vs irreversible separation | `python g1_run.py` |
| PCM solver verification against the Neumann Stefan solution | `python pcm_stage_a_verify.py` |
| PCM CRLB and the leak audit that reinterprets it | `python pcm_stage_b_crlb.py`, `python pcm_stage_b_audit.py` |
| PCM regime sweep across heating rate | `python pcm_stage_e_regime.py` |
| Everything end to end, as notebooks | `Part7_CoreValidation.ipynb`, `PCM_Identifiability.ipynb` |

CPU only. No GPU is needed anywhere; `torch` runs on CPU throughout. The full
Part 7 notebook takes about 40 minutes.

`results_inventory.md` traces **every quantitative claim to the file that
produced it**, with provenance tiers. Start there if you want to check a specific
number rather than run everything.

---

## How to read the repository

The documents are the deliverable as much as the code is.

- **`FINDINGS.md`** — Part 7. The core-temperature result, scored predictions,
  and the sections on what is *not* established. Read §5.5 and §5.6 if you only
  read two things; they are the parts unfavourable to the method.
- **`PCM_FINDINGS.md`** — the phase-change identifiability study. Five
  independent lines of evidence that the melt front is not what is being
  observed.
- **`GAPS.md`** — every `[VERIFY]` marker raised during the project and what
  closed it, including one claim that was **dropped** because its supporting run
  could not be verified (G-1).
- **`results_inventory.md`** — claim-to-file provenance map.
- **`PCM_PREDICTIONS.md`** — registered before the corresponding results existed.
  The git history is the evidence: commit `101bf7a` registers Q1–Q6 and precedes
  `e1d2b02`, which produces the first Stage B number, both dated 2026-08-06.
- **`thesis/ch6_validation.md`** — longer-form validation chapter.

## Methodological commitments

These are the rules the project was run under, and they are the reason several
results here are negative:

1. **Compare against the trivial baseline before claiming a method works.** The
   entire PINN apparatus is beaten by one line of algebra. Without that baseline
   in the table, this repository would have reported a 0.61 K success and buried a
   0.14 K result that needed no network.
2. **Report the parameter correlation matrix, not just marginal standard
   deviations.** The CRLB rated `{h, k, rho·c_p}` identifiable at 3.9 %; the fit
   then drove `k` 1140 % away, roughly 300× outside the CRLB ellipse. The 0.99
   `k`–`rho·c_p` correlation was the warning that a marginal bound hid.
3. **Pre-register predictions and score the misses.** Four hits and two misses on
   Part 7; three and three on the PCM study. The misses are recorded with why they
   were wrong.
4. **Audit for leaks before reporting a good number.** An earlier version leaked
   the core channel through the initial condition. It was found by audit,
   quantified at 0.6 % on core RMSE, removed, and every downstream constant
   re-derived.
5. **Test discretisations against closed-form solutions at machine precision.**
   Two silent bugs — a boundary source broadcast to every cell, and an
   O(dr/4R) wall-conductance bias — were caught only because the analytic steady
   state reproduces to ~1e-13.

## What is not established

Stated here because it belongs in the README, not only in the papers.

- One cell, two drive cycles, one ambient (~8 °C). Nothing here speaks to
  temperature dependence, ageing, or cell-to-cell variation.
- **The Biot number is an input, not an output.** Everything rests on it and it is
  not identifiable from a single surface channel; its correlation with `k` is
  0.998.
- The core prediction is essentially a function of `k`, and surface data cannot
  identify `k`. The headline therefore rests on an input.
- `R_eff` is a blend, not an internal resistance. The entropic term is off in the
  headline; switching it on moves core RMSE by 0.4 %.
- No shape `R(x)` was recovered over `[0,1]`. The depth-of-discharge window spans
  11–14 % of the axis, so order 1 is a narrow-band interpolant.
- A plain full-field PINN does not work here at all — it failed the forward gate
  at 39.2 %. Everything reported depends on the base-subtracted formulation.

## Citation

See `CITATION.cff`. If you use the datasets, cite the datasets — `DATA.md` lists
them. Their authors did the expensive part.

## Licence

MIT for code, CC BY 4.0 for documents, figures and results. No third-party
measurement data is redistributed. See `LICENSE` and `DATA.md`.

## What this repository is an export of

This is a filtered export of a private working repository. All 29 commits, their
messages, authors and dates are preserved, so the pre-registration ordering is
verifiable here. Commit **hashes differ** from the private original because the
export was filtered; `PCM_FINDINGS.md` carries a note on that.

Four things were removed or held back, each for a stated reason:

- **Unpublished proposal material** for a third-party collaboration. Not mine to
  publish.
- **`results/mr_cache/*.npz`** — the raw Catenaro & Onori time series. Source
  measurement data rather than analysis output, and it regenerates on first run.
- **Absolute developer paths**, replaced by `paths.py`. Seven files carried a
  hardcoded path that made them unrunnable anywhere else.
- **The manuscripts.** This repository is the code and the research record; the
  papers are published through the journal, not here.

Author identity was also consolidated. Some commits were made on a machine whose
git config carried a family member's identity; all 29 are my own work and are now
attributed accordingly.
