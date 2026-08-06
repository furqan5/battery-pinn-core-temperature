# Results inventory — every quantitative claim, with its source

Compiled 6 Aug 2026, before any thesis text was drafted. Nothing here is a number I have
supplied; every entry points at a file. Where a number is quoted in conversation or in the
session brief but does not appear in a file, it is listed in §7 rather than used.

---

## 1. Provenance tiers

Not all numbers in this project stand on the same footing, and the thesis must not flatten them.

| Tier | Meaning | Applies to |
|---|---|---|
| **V1 — re-derived** | Recomputed from committed code + committed data during this inventory. Highest confidence. | All Part 7 results; **Part 3 (re-executed 6 Aug 2026, see §3.0)** |
| **V2 — re-executable, not re-executed** | Code and input data are both present in this repository; the numbers live in notebook prose because the cells carry no stored outputs. | Part 6 (re-execution in progress at time of writing) |
| **V3 — self-contained, not re-executed** | Code present, inputs synthetic (generated in-notebook), no stored outputs. | Parts 2, 4, 5, 5b |
| **V4 — external, unverifiable** | Referenced in project documents; code and raw output never available to this repository. | The Iontech LFP multi-C-rate run |

**A caveat that belongs in the methodology chapter.** The Part 2–6 notebooks have **no stored cell
outputs**. Their reported values are the author's transcription of runs on other machines
(`MANIFEST.md` states all six "have been executed at least once with outputs reproduced across two
independent environments"). Part 7, by contrast, was executed end to end in this repository
(commit `addd906`: 54 cells, 24 code, 0 errors, 39.9 min) with outputs embedded. That asymmetry is
worth one honest sentence in the thesis rather than concealment.

---

## 2. Part 7 — validation against measured core (tier V1)

All re-derived from `results/stage_e_*.npz` + `part7_lib.py` during this inventory. Agreed with
`FINDINGS.md` to all quoted digits.

### 2.1 Headline

| Method | DS2 core RMSE | DS1 core RMSE | Source |
|---|---|---|---|
| Quasi-steady one-liner `T_core = T_surf + (Bi/2)(T_surf − T_inf)` | **0.1409 K** (max 0.3742) | **0.3366 K** (max 0.6900) | `stage_g_quasisteady.py`, re-derived |
| Classical transient finite-volume | 0.8943 K | 1.0544 K | `stage_f_classical_control.py`; `stage_e_replicate_ds1.py` |
| Inverse PINN, order 0 | 0.6135 K (max 1.9093, bias −0.3577) | 0.5024 K (max 1.3208, bias −0.3245) | `results/stage_e_shape0.npz`, `stage_e_ds1_shape0.npz` |
| Inverse PINN, order 1 | 0.5710 K (max 1.4981, bias −0.2976) | — | `results/stage_e_shape1.npz` |

Measured core-minus-surface gradient peaks: **6.5426 K** (DS2), **7.1314 K** (DS1).

### 2.2 The quasi-steadiness mechanism

| Quantity | DS2 | DS1 | Source |
|---|---|---|---|
| Measured ratio (core−surf)/(surf−ambient), hot window | **0.6161 ± 0.0126** (2.0 %) | **0.6218 ± 0.0103** (1.7 %) | re-derived |
| Bi (from other record's surface-only fit) | 1.2142 | 1.1487 | `stage_e_inverse.py:H_FIXED,K_FIXED` |
| Bi/2 | 0.6071 | 0.5744 | — |
| Diffusion time R²/α | 1043 s | 981 s | `stage_g_quasisteady.py` |

### 2.3 Estimator internals

| Quantity | Value | Source |
|---|---|---|
| R_eff, DS2 order 0 | 13.1081 mΩ (seed spread 0.489 %) | npz, re-derived |
| R_eff, DS1 order 0 | 13.2651 mΩ (spread 1.089 %) | npz, re-derived |
| Core RMSE across 6 seeds, DS2 order 0 | mean 0.6274, sd 0.0086, range 0.6135–0.6369 | re-derived |
| **Non-convergence rate, Part 7** | **0 / 6 on every configuration** | `results/stage_e.log`, `stage_e_ds1.log` |
| Electrical anchor (regression of V on I) | 14.492 mΩ (DS2), 13.872 mΩ (DS1) | `part7_lib.py:Record.R_ohmic_reg` |
| Error concentration | first 10 min = 17 % of record, **56 % of squared error**; RMSE 1.1146 K early vs 0.4470 K after | `results/stage_f.log` |

### 2.4 The two damaging results about the PINN

**PDE violation** (`results/consistency.log`):

| | surface RMSE | core RMSE |
|---|---|---|
| PINN's own field | 0.1974 K | 0.6135 K |
| Same R_eff, physics enforced exactly | 0.8395 K | 1.5238 K |
| **Difference = PDE violation** | **0.8029 K** | **1.0570 K** |

**Data-weight fragility** (`results/weight_sweep.log`):

| w_data | R_eff | surface RMSE | PDE+BC residual | core RMSE |
|---|---|---|---|---|
| 5 | 13.32 mΩ | 0.4827 K | **0.1015 (lowest)** | **1.0500 K** |
| 20 | 13.62 mΩ | 0.3025 K | 0.1136 | 0.7403 K |
| **200 (used)** | 13.05 mΩ | 0.1963 K | 0.1436 | 0.6332 K |
| 2000 | 11.26 mΩ | 0.1374 K | 0.2665 | 0.7721 K |

Truth-free selection (lowest PDE+BC residual) picks w_data = 5 → 1.05 K → **P5 would have been a MISS**.

### 2.5 Supporting Part 7 numbers

| Quantity | Value | Source |
|---|---|---|
| k-sensitivity: core RMSE at k = 0.30 / 0.3937 / 0.55 | 1.8209 / 0.8943 / 1.6983 K | `results/classical_control.log` |
| Surface RMSE at same k | 0.5433 / 0.4652 / 0.4007 K (monotone) | same |
| Forward gate (split formulation) | 4.66 % (DS2), 1.39 % (DS1) core-surface rel. error | `results/stage_d_split.log` |
| Forward gate, naive full-field PINN | **39.2 % — FAILED** | commit `943357a` |
| Fourier-feature sweep (constant-source control) | 0.086 % (n_freq 0) → 4.883 % (n_freq 64), monotone | commit `943357a` |
| Trap 5.1 control (no I(t)²) | 2.7× worse surface, 1.8× worse core | `results/stage_e.log` |
| FV solver verification | 24/24 checks to ~1e-13 | `tests_fd.py` |
| Lumped-base verification | 9/9 | `tests_split.py` |
| Leak: core RMSE with / without | 0.8895 / 0.8950 K (0.6 %) | `results/leak_audit.log` |
| Leak propagation: h, k re-derived | h 37.0678→37.0607; k 0.390843→**0.393747** | commit `531ed0c` |
| Skill vs baselines, DS2 | B0 (T_core=T_surf) 5.1855 K; B1 (+mean grad) 1.2101 K | `results/skill.log` |
| Gradient correlation, PINN vs measured | 0.9330 (DS2), 0.9882 (DS1) | `results/skill.log` |
| Entropic sensitivity | core RMSE 0.8895 → 0.8927 K (0.4 %) | `results/classical_control.log` |
| Non-uniform generation hypothesis | **refuted** — R_eff rose 14.30→14.36 mΩ with β | `results/nonuniform.log` |
| Classical residual anatomy | 100 % model-form; autocorr 1/e at 423 s; corr with time r = +0.51, with I² r = +0.10 | `results/residual_anatomy.log` |

### 2.6 Cell and data provenance (Part 7)

| Item | Value | Source |
|---|---|---|
| Cell | A123 ANR26650M1, LFP, 2.3 Ah | Richardson & Howey 2015 |
| Data | 2 HEV drive cycles, core + surface thermocouples | `Data_Sets\EKF-Battery-Impedance-Temperature-master\...\Temperature_data.mat` |
| Channel identities | `T_data(:,1)=t, (:,2)=T_surf, (:,3)=T_core, (:,4)=T_inf` | `MainScript.m:136-139` |
| R_o, V_b, ρ, cp, k, h | 0.0129 m, 3.4219e-5 m³, 2107 kg/m³, 1171.6 J/kgK, 0.404 W/mK, 39.3 W/m²K | `MainScript.m:87-97` |
| Current sign | **I > 0 = charge**; corr(I,V) = +0.90; regression slope +13.9 mΩ | `part7_lib.py:Record`, session diagnostics |
| Records | DS1 5981.8 s / 5439 pts; DS2 3543.1 s / 3222 pts, 1.1 s raw | `part7_lib.py` |
| DoD window spanned | 0.138 (DS1), 0.112 (DS2) — 11–14 % of axis | Stage A |

---

## 3. Part 3 — identifiability and classical inverse (tier V1, re-executed)

The intellectual spine. `Part3_CRLB_Classical_Inverse.ipynb`.

### 3.0 Re-execution outcome — 6 Aug 2026

Executed in this repository against the local NASA archive: **24/24 code cells, 0 errors, 109 s**
(`verify/Part3_CRLB_Classical_Inverse_executed.ipynb`, `verify/run.log`). **Every transcribed
value reproduced**, to the last quoted digit in all cases checked, including the CRLB table, the
degeneracy eigenvector, the lumped-inverse bias, the profile likelihood, the final 168-cycle fit,
and the EIS comparison. The chapter's numbers can now be cited as verified here, not merely
reported.

Three things the re-execution added that the prose did not state:

1. **The 15.87 K / 14.3 K inconsistency is resolved.** Cell 36 prints, verbatim: *"relative figure
   of merit: RMSE 0.228 K over a rise of 15.87 K = 1.44 %"*, and cell 17 gives "temperature rise
   across cycles: median 15.87 K, min 13.74, max 17.31". **15.87 K is correct**; the "14.3 K, 1.6 %"
   in the §4.10 prose is a transcription slip. Use 1.44 %.
2. **The `.mat` mirror carries no `metadata.csv`**, so cell 18 prints "No metadata — cannot
   stratify". The "B0005 @ 24 °C stratum" is therefore *all 168 discharge cycles in `B0005.mat`*,
   not a subset selected on logged ambient. The label is nominal. This does not change any number
   but it should be described accurately in the thesis. `[VERIFY: whether B0005 cycles are all at
   24 °C nominal — the CSV mirror's metadata would settle it]`
3. **The source-order limit is directly printed** (cell 29), which is stronger evidence than the
   prose summary:

| Model | scaled cond(F) | worst relative CRLB sd |
|---|---|---|
| order 1, h fixed | 2.412×10² | 1.63 % |
| **order 2, h fixed** | **1.070×10⁴** | **7.48 %** |
| order 3, h fixed | 6.136×10⁵ | 52.69 % |
| order 4, h fixed | 3.269×10⁷ | 362.97 % |

"Highest order supported with h fixed (cond < 1e5 and all sd < 25 %): **2**" — printed by the
notebook, not asserted in prose.

### 3.1 CRLB, synthetic slab, σ = 0.10 K

| Configuration | scaled cond(F) | relative CRLB sd | Cell |
|---|---|---|---|
| {R_int, h, ρc_p}, 2 sensors | **8.800×10⁹** | 421.658 / 421.649 / 416.666 % | 9 |
| {R_int, h}, 2 sensors | 2.542×10¹ | 0.026 / 0.031 % | 9 |
| {R_int, h}, 1 sensor | 2.542×10¹ | 0.037 / 0.044 % | 9 |
| {R_int}, 1 sensor | 1.000 | 0.014 % | 9 |

Weakest eigenvector **(+0.580, +0.580, +0.573)** ≈ (1,1,1)/√3. Observables: ΔT_ss = 10.2022 K,
τ = 2475.0 s — two observables, three unknowns.

### 3.2 CRLB for the enriched source (real-data configuration)

| Parameter set | scaled cond(F) | relative CRLB sd | Cell |
|---|---|---|---|
| {R₀, h} | 2.43×10² | R₀ 0.16 %, h 0.44 % | 24 |
| **{R₀, β, h}** | **3.58×10⁶** | **β 36.8 %, h 38.9 %** | 24 |
| {R₀, β}, h known | 2.26×10² | R₀ 0.17 %, β 0.42 % | 24 |

Mock confirming the bound: true β = 1.500, h = 11.0 → unconstrained fit returned **β = 2.305
(+54 %), h = 17.15 (+56 %)** at **RMSE 0.051 K**. Excellent trajectory, both parameters wrong.

### 3.3 Final real-data fit, B0005 @ 24 °C, 168 cycles

| Quantity | h = 29.57 (truncated) | **h = 24.81 (profile likelihood)** |
|---|---|---|
| usable | 159/168 (94.6 %) | **168/168 (100 %)**, 0 railed |
| RMSE median | 0.260 K | **0.228 K** |
| R₀ median | 157.5 mΩ | **153.0 mΩ** |
| c₁, c₂ | −1.255, +2.276 | **−1.352, +2.192** |
| R_end median | 306.8 mΩ | **280.4 mΩ** |
| corr(capacity, R₀) | −0.906 | **−0.908** |
| corr(capacity, R_end) | −0.986 | **−0.981** |
| source minimum | x = 0.276, 0.827× | **x = 0.308, 0.792×** |
| rise to cutoff | 2.02× | **1.84×** |

Bootstrap x_min: median **0.318**, CI [0.304, 0.324], 6.3 % — the most stable quantity in the
analysis. (The 0.308 vs 0.318 difference is explained in cell 38: median of a ratio ≠ ratio of
medians. Not an inconsistency.)

Shape stability across ageing: c₁ IQR 0.167 (12.4 %), c₂ IQR 0.144 (6.6 %), R₀ IQR 14.8 %.

### 3.4 External validation — the only one in the project

| Source | Value |
|---|---|
| EIS R_e median | 55.9 mΩ |
| EIS R_ct median | 77.5 mΩ |
| **EIS R_e + R_ct** | **134.3 mΩ** |
| **Thermal R₀** | **153.0 mΩ** |
| **Ratio** | **1.14×** |

278 impedance records, same battery and ambient. `R_end/(R_e+R_ct) = 2.09×`.

### 3.5 The noise-floor result

`ac1(d2) = +0.3190` against a white-noise reference of **−4/6 = −0.6667**; smallest temperature
step 0.00135 K; σ̂ = 0.0124 K. Reference behaviour: white noise −0.670/−0.610/−0.689 at
σ = 0.10/0.02/0.005 K; 3-point boxcar ≈ −0.50; linear interpolation ≈ 0.00; noiseless signal +0.99.
**Conclusion: no white-noise floor exists in this dataset**, so any "N× the noise floor" statistic
is void. Part 6 re-ran it on the raw `.mat` mirror and got identical values to five significant
figures — the smoothing is in the original NASA files, not the curation.

Honest figure of merit: **0.228 K over a 15.87 K rise = 1.44 %**.

> ✅ **Resolved by re-execution.** The prose carried two versions ("15.87 K rise = 1.44 %" and
> "14.3 K rise, 1.6 %"). The notebook prints 15.87 K, and cell 17 confirms it as the median rise
> across the 168 cycles. **Quote 1.44 %.**

---

## 4. Part 2, 4, 5, 5b (tier V3)

### 4.1 Part 2 — FD truth model

Slab, L = 0.018 m; ρ = 2500, cp = 1100, k = 2.5, h = 10 (all **(b)**); V_cell = 1.654×10⁻⁵ m³;
I = 2.5 A, R_int = 30 mΩ → q''' = 1.134×10⁴ W/m³. **Bi = 0.036.**

τ_diff ≈ 89 s; τ_lump = 2475 s; **τ₁ (first Robin eigenmode) = 2506 s**.

10/10 pre-registered checks pass: B1 7.72e-3 K, B1b 7.86e-3 K, B2 0.1835 K, **B3 2504.7 s**,
B4 307.626 K, B5 300.287 K, B6 3.05e-5 K, B7 6.22e-7, B8 3.20e-5 K, B9 0.00 K.

**The volume trap**: Bernardi gives watts, the PDE needs W/m³; dividing by a unit-area slab volume
rather than the real cell volume is wrong by **~1000×**.

### 4.2 Part 4 — forward PINN gate

| Metric | Value |
|---|---|
| max \|T_PINN − T_FD\|, whole field | **0.00100 K** |
| field RMSE | 0.00046 K |
| **core-surface bow, PINN vs FD** | **0.17358 vs 0.17324 K → 0.196 %** |
| Adam 3000 + L-BFGS | 115 s + 17 s (435 closure calls) |

**Ablation 1 — scaling** (this is the "37× smaller loss, 120× worse answer" claim, verified):

| Formulation | final residual loss | max abs error |
|---|---|---|
| Physical units | **1.699×10⁻⁵** | **9.7997 K** |
| Scaled | 6.272×10⁻⁴ | 0.0812 K |

Loss ratio 6.272e-4 / 1.699e-5 = **36.9×**; error ratio 9.7997 / 0.0812 = **120.7×**. Both confirmed.

**Ablation 2 — L-BFGS frozen collocation set:**

| Metric | broken (1 call) | fixed (435 calls) | factor |
|---|---|---|---|
| max abs error | 0.02587 K | 0.00100 K | **26×** |
| field RMSE | 0.01282 K | 0.00046 K | 28× |

**Ablation 3 — loss weighting:** w_bc = 1 → 1.2783 K; w_bc = 100 → 0.0264 K (**48×**).

### 4.3 Part 5 / 5b — inverse PINN on synthetic data

CRLB for this configuration: scaled cond 25.38, R_int 0.116 %, h 0.138 %.

| Quantity | Lumped rival | Inverse PINN | CRLB sd |
|---|---|---|---|
| R_int error | −1.813 % | **−0.307 % ± 0.115 %** | 0.116 % |
| ratio R/h error | −0.070 % | +0.062 % | — |
| core field RMSE | **not available** | 0.009–0.028 K | — |

**Data-weight sweep (Part 5b) — directly relevant to Part 7:**

| w_data | mean R_err | sd | mean \|bow\| | mean Lr |
|---|---|---|---|---|
| 2 | −0.309 % | 0.253 | 6.93 % | 7.10×10⁻³ |
| **20** | **−0.307 %** | **0.115** | **3.34 %** | **4.34×10⁻³** |
| 200 | −1.293 % | 0.648 | 3.93 % | 4.49×10⁻³ |
| 2000 | −2.706 % | 1.951 | 2.71 % | 1.05×10⁻² |

**Non-convergence: 1 in 6.** Seed 3 returned R_err = **−4.271 %** (20× the typical error) with
**Lr = 5.048×10⁻³ sitting inside the range of the good runs (2.97–5.05×10⁻³)** — the truth-free
selection criterion cannot see the failure. Converged runs: **−0.113 % ± 0.104 %** (n = 5) against
CRLB sd 0.116 %, ratio 0.90.

**The metric trap (Part 5, §5.1):** core max error 0.0492 K looked like a pass, while the *bow* at
7200 s was **0.23251 vs 0.17325 K = +34.2 % wrong**. A lenient metric hid a large spatial error.

CPU vs GPU drift: worst relative difference **5.67×10⁻²** on the bow error → GPU is a separate
baseline.

---

## 5. Part 6 — inverse PINN on real data (tier V2 → **partially V1, with a caveat**)

### 5.0 Re-execution outcome — 6 Aug 2026

Executed in this repository: **9 code cells, 0 errors, 1781 s**
(`verify/Part6_Real_Data_PINN_executed.ipynb`). Both residual variants are present in the archived
notebook — the broken one at cell 6 (`mult = 1 + c1·x + c2·x²`) and the corrected one at cell 14
(`mult = (1 + c1·x + c2·x²)·Ir²`) — so the fix *is* in the artifact, not only in the prose.

**The qualitative conclusion reproduces. The quoted parameter value does not.**

| Variant | R₀ (mΩ) | c₁ | c₂ | RMSE (K) | Lr | bow ratio |
|---|---|---|---|---|---|---|
| **re-executed, fixed, w = 20** | **152.83** | −1.532 | 2.141 | 0.1269 | 2.012×10⁻³ | 0.74 |
| *prose, fixed, w = 20* | *144.78* | *−1.232* | *+1.845* | *0.1402* | *8.56×10⁻³* | *0.73* |
| **re-executed, fixed, w = 200** | **149.06** | −1.416 | 2.040 | 0.0704 | 5.486×10⁻³ | 0.86 |
| *prose, fixed, w = 200* | *149.48* | *−1.441* | *+2.069* | *0.0713* | *5.48×10⁻³* | *0.86* |
| re-executed, broken, w = 20 | 143.86 | +0.023 | 0.430 | 0.6157 | 2.749×10⁻³ | 1.11 |
| cycle-0 classical target | 144.35 | −1.305 | 1.934 | 0.1297 | — | — |

**Reproduces robustly:** the defect diagnosis. Without the current factor c₁ sits at +0.023 — the
wrong sign — and RMSE is 0.62 K; with it, c₁ is firmly negative and RMSE falls to ~0.13 K. The
w = 200 row matches the prose to 0.3 %.

**Does not reproduce:** the headline parameter at w = 20. A fresh run returns **152.83 mΩ**, which
against the cycle-0 classical of 144.35 mΩ is **+5.87 %**, not the +0.30 % quoted. Note the fresh
run has *better* physics (Lr 2.0×10⁻³ against 8.6×10⁻³) while sitting further from the classical
value — the same failure mode Part 5b documented, where the truth-free selection criterion cannot
distinguish a good run from a badly-parameterised one.

> ⚠ **Consequence for the prediction ledger.** Part 6 scores **G2** ("R₀ within ±3 % of classical")
> as a **HIT** at +0.30 %. On the re-executed run it is **+5.87 %, a MISS**. The verdict depends on
> which training run is taken, and no seed spread was ever reported for this fit.
>
> **Recommendation for the thesis:** quote Part 6's R₀ as a range across runs, not to five
> significant figures, and either re-score G2 as indeterminate or run the 5+ seeds that the
> project's own methodology requires elsewhere. `[VERIFY: seed spread for the Part 6 fixed-residual
> fit at w = 20 — currently n = 1 in the notebook and n = 1 here, and they disagree by 5.6 %]`

This is not a defect in the Part 6 argument, which is about a residual defect and survives intact.
It is a defect in the *precision* with which a single stochastic fit was reported, and it is
exactly the kind of thing re-execution exists to catch.

### 5.1 The missing-I(t)² story

| Variant | R₀ | c₁ | c₂ | R_end | x_min | RMSE | bow ratio |
|---|---|---|---|---|---|---|---|
| broken (no I²), w = 20 | 143.86 | **+0.023** | +0.430 | 209.0 | — | 0.6157 | 1.11 |
| **fixed, w = 20** | **144.78** | **−1.232** | **+1.845** | 233.5 | **0.334** | 0.1402 | 0.73 |
| fixed, w = 200 | 149.48 | −1.441 | +2.069 | 243.4 | 0.348 | 0.0713 | 0.86 |
| *cycle-0 classical target* | *144.35* | *−1.305* | *+1.934* | *235.1* | *0.337* | *0.1297* | — |

**The defect's signature**: without the I(t)² factor, c₁ flips sign (+0.023 vs −1.232). The
amplitude R₀ was barely affected; the *shape* was destroyed. This is the single clearest argument
in the project for keeping a classical rival beside the neural method.

Quoted result: R₀ **144.78 mΩ** (+0.30 % vs cycle-0 classical), x_min agreeing to 1 %.

Lumped τ = 997.6 s vs 1-D mode-1 τ = 1027.5 s → **+2.99 %** (small-Bi theory Bi_a/3 = +2.98 %).

---

## 6. The full prediction ledger

**Note: the brief says "roughly a dozen scored predictions". The actual count is 35 scored
predictions plus 9 verification checks** — considerably stronger than claimed.

| Series | Where | Hits | Misses | Partial | Void |
|---|---|---|---|---|---|
| B1–B9 | Part 2 | 9 (verification checks, not predictions) | 0 | — | — |
| C1–C5 | Part 3 | 3 | 1 (C3) | 1 (C5) | — |
| D1–D5 | Part 3 | 3 | 2 (D2, D5) | — | (D1, D3 void on first pass) |
| E1–E7 (E6 split) | Part 4 | 7 | 1 (E6b) | — | — |
| F1–F6 | Part 5 | 4 | 2 (F1, F5-as-run) | — | — |
| G1–G5 | Part 6 | 3 | 1 (G5) | 1 (G4 band retracted) | — |
| P1–P6 | Part 7 | 4 | 2 (P1, P3) | — | — |
| **Total** | | **24** | **9** | **2** | — |

Plus **six self-logged mid-session misses** (three in Part 7, one in Part 3 §4.10 on R₀ ∝ h
scaling, one in Part 5b on drift, one in Part 6 on the comparator).

Most informative misses, for the honesty chapter:
- **C3** — confused *shape* redundancy with *statistical* redundancy; second sensor gave the
  standard √2, not <5 %.
- **D5** — predicted thermal/EIS 1.5–5×, got 1.14×. Mis-specified prediction, favourable result.
- **E6b** — Gate B error was **0.0189 K near the switch vs 0.0540 K away**, the opposite of predicted.
- **P1** — gradient 43 % above the predicted interval.
- **P3** — predicted unidentifiable at >25 %; observed 1.3–10.6 % at orders 1–3.

---

## 7. Claims I could not source, and where I disagree with the brief

These are the items §0 asked me to flag rather than write up.

### 7.1 "The reversible-heat contamination diagnosed in the Part 6 work"

**Not supported as stated.** Part 6 does not diagnose reversible-heat contamination. It
*absorbs* the entropic term into the effective coefficients c₁, c₂ and explicitly declines to
model it (Part 3 §4.3: "the effect is absorbed, not modelled"; Part 6 cell: "Adding it explicitly
would be degenerate with c₁, c₂").

The contamination hypothesis comes from the **Iontech LFP multi-C-rate run**, which `MANIFEST.md`
states plainly: *"Its code was never shared with me — `train_pinn.py` was referenced but not
uploaded — so it is not in this package and I cannot verify it."* Tier V4.

### 7.2 "The multi-rate evidence points to a reversible-heat artifact"

**The project's own manifest contradicts this.** `MANIFEST.md` §"Two corrections to carry forward":

> *"LFP cannot test an LCO hypothesis. The false minimum in LCO comes from a sign change in
> dU_oc/dT near 66 % SOC. LFP does not share that structure. The Iontech test shows the method can
> **detect** entropic contamination; it says nothing about whether B0005's specific minimum is
> entropic. That still needs either an LCO entropy-profile subtraction or a same-chemistry
> replication."*

And the multi-rate fit itself is flagged as unusable: reported c₁ = −3.01, c₂ = 0.07 gives
R(x)/R₀ crossing zero at x ≈ 0.34 and reaching **−1.94** at full discharge — negative irreversible
heat, thermodynamically forbidden — **and its RMSE was never reported**.

**Recommendation.** The thesis should say the interior minimum's origin is **unresolved**, list
the entropic hypothesis as the leading candidate, and state the two experiments that would settle
it. Writing it as diagnosed would be the one indefensible claim in the document, and a viva
examiner who reads MANIFEST.md would find it immediately.

`[VERIFY: Iontech multi-C-rate run — RMSE, code, and raw output, or drop the claim entirely]`

### 7.3 The stated mechanism for why the one-liner wins

The brief says the mechanism is that the relation *"re-anchors on the measured surface at every
step and never needs the heat-generation term, so heat-source error cannot accumulate."*

That is **half right and worth stating carefully**. Both parts are supported:

- *Never needs the source* — correct, and it is why the one-liner beats the transient FV, which
  must get R_eff right.
- *Re-anchoring* — correct but **not sufficient**, and this was tested. Giving the PINN the same
  anchor (its gradient + the measured surface) still yields 0.5990 K against the one-liner's
  0.1409 K (`results/quasisteady.log`). So anchoring alone does not close the gap.

The dominant mechanism documented in `FINDINGS.md` is **quasi-steadiness** — the ratio holds to
~2 %, so the algebraic conversion of surface rise to gradient is simply accurate. The thesis should
lead with quasi-steadiness and cite source-independence as the second, reinforcing reason.

### 7.4 Scope drift between chapters — needs an explicit statement

The chapters do not share a cell, a chemistry, or a geometry:

| Part | Geometry | Cell | Chemistry | Bi |
|---|---|---|---|---|
| 2, 4, 5, 5b | 1-D slab | generic 18650 | agnostic | **0.036** |
| 3, 6 | lumped / slab | NASA B0005 18650 | **LCO** | ≈ 0.07 |
| 7 | **1-D cylinder** | A123 26650 | **LFP** | **1.21** |

**Bi differs by a factor of ~34 between the early chapters and Part 7.** The thesis argument (a
timescale/Biot criterion for when a transient solver earns its place) actually *needs* this spread
— it is a strength, not an embarrassment — but it must be stated in Chapter 1, not discovered by
the examiner in Chapter 6.

### 7.5 Smaller items

- **"1.44 % relative accuracy"** appears with two different denominators (§3.5). Resolve before use.
- **Part 5b's w_data sweep found 20 optimal and 200 markedly worse** (−0.307 % vs −1.293 %). Part 7
  used w_data = 200. The two are on different problems, but a viva examiner may well connect them;
  the thesis should address it in the Part 7 chapter rather than leave it.
- **`FINDINGS.md` §6 has two items numbered 6.** Cosmetic, but fix before submission.
- No number in the brief was found to be fabricated; the issues above are framing and sourcing.

---

## 8. Figure list

| # | Figure | Produced by | Exists now? |
|---|---|---|---|
| 6.1 | Measured core, predicted core, measured surface + error trace | `stage_f_figure.py` | ✅ `figures/core_validation.png` |
| 6.2 | Measured vs predicted core−surface gradient | `stage_f_figure.py` | ✅ `figures/gradient_validation.png` |
| 3.1 | CRLB condition number vs parameter set | — | ❌ table only; would need new code |
| 3.2 | Recovered R(x)/R₀ with cycle-to-cycle band | Part 3 cell 7m | ❌ code present, not executed here |
| 3.3 | corr(capacity, R₀) ageing scatter | Part 3 | ❌ same |
| 4.1 | PINN vs FD field error map | Part 4 | ❌ same |
| 4.2 | Scaled vs physical-units ablation | Part 4 | ❌ same |
| 5.1 | Seed scatter vs CRLB ellipse | Part 5b | ❌ same |
| 6.3 | k-sensitivity: surface RMSE and core RMSE vs k | `stage_f_classical_control.py` | ❌ data in log, plot not written |
| 6.4 | Data-weight sweep: accuracy vs PDE violation | `stage_g_weight_sweep.py` | ❌ data in log, plot not written |
| 1.1 | Bi across the three cell/geometry cases | — | ❌ new |

Two of eleven exist. Six more are cheap (data already computed, plotting code only). Three need
notebook re-execution.

---

## 9. Recommended reading order for the examiner-facing argument

1. §2.2 quasi-steadiness — the physical claim
2. §2.1 three-method comparison — the measurement
3. §3.1–3.2 CRLB — why the parameter sets were chosen before fitting
4. §2.4 PDE violation and weight fragility — the honesty core
5. §7 — what the thesis must not claim
