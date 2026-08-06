# GAPS — every `[VERIFY: …]` marker, and what closes it

Collected 6 Aug 2026. Ordered by how much damage each would do if an examiner found it first.

---

## Tier 1 — must be resolved before submission

### G-1. The reversible-heat claim — **CLOSED, and the claim is NOT supported**

**Done:** the separation was implemented, run on LCO, and its RMSE reported. Code:
`g1_survey_currents.py`, `g1_inspect_multirate.py`, `g1_entropic_separation.py`, `g1_run.py`,
`g1_figure.py`. Output: `results/g1_separation.log`, `results/g1_separation.npz`,
`figures/g1_separation.png`.

**The data that made it possible.** MANIFEST.md said the question needed "an LCO entropy-profile
subtraction or a same-chemistry replication", and that LFP could not test an LCO hypothesis. It
can be done: the NASA archive contains LCO cells cycled at **1, 2 and 4 A within one cell at one
ambient** — B0038/B0039/B0040 at 44 °C. Same chemistry as B0005.

**Method, which assumes no source shape at all.** With discharge-positive current,
Q/I = I·R(x) − T·(dU/dT)(x), so at fixed x a plot of Q/I against I has slope R and intercept
−T·dU/dT. Q itself comes from the lumped balance Q = C·dT/dt + K(T − T_amb) by direct
differentiation, and K comes from 60 post-discharge relaxation tails where the source is off
(τ = 800 s, K = 0.05684 W/K). Nothing is fitted to a heat-generation model.

**The number that was never reported.** Median linear-fit RMSE **0.01249 V (4.89 % of the fitted
values)**, worst 0.02660 V (11.49 %). Three currents, two parameters — one degree of freedom per x.

**Result 1 — the interior minimum SURVIVES the separation.** After removing the reversible term,
R(x) still falls to a minimum (x ≈ 0.72 in the well-behaved region) and rises again. On this
evidence the minimum is **not** a reversible-heat artefact.

**Result 2 — and the reason the first result is the only usable one.** Scaling K by ±30 %:

| | ±30 % in K | verdict |
|---|---|---|
| ohmic term R(x) | 58–62 mΩ, band **2 %** of value | well conditioned |
| reversible term dU/dT | −0.163 to −0.416 mV/K, band **112 %** of value | **not determined** |

**A prediction of mine, refuted by the control I wrote to test it.** I argued in the module
docstring that K-error would load onto the slope and leave the intercept well conditioned. The
sweep shows the reverse, almost exactly proportional in K. The docstring now records the
correction. This also means **the 4.89 % fit RMSE understates the true uncertainty**: the dominant
error is systematic in K and leaves no trace in the residual — the same lesson as everywhere else
in this project.

**What this does and does not settle.**

- Settles: the claim no longer rests on an unverifiable external run. There is local, runnable
  code and a reported RMSE.
- Settles: on the one multi-rate LCO test available, the entropic-artefact explanation is **not
  supported**.
- Does **not** settle B0005 specifically. Different cell, different ambient (44 °C against 24 °C,
  where resistance is much higher), only 6 cycles at the 4 A leverage point, one degree of freedom
  per fit, and a visible artefact band at x < 0.28 where coulomb normalisation across cells is
  least reliable.

**How to write it.** The interior minimum is a feature of the recovered *effective* coefficient,
and single-rate data cannot decompose it. The one multi-rate test on the right chemistry does not
support the entropic explanation, and the reversible term that test recovers is itself undetermined
to better than a factor of two. Say that, and claim no more.

---

### G-1 (original text, retained for the record). The reversible-heat claim has no verifiable source
**Where:** the thesis brief asks for "the reversible-heat contamination diagnosed in the Part 6
work"; the intended supporting evidence is the Iontech LFP multi-C-rate run.

**Problem:** Part 6 does not diagnose it — it *absorbs* the entropic term into the effective
coefficients c₁, c₂ and explicitly declines to model it. The contamination hypothesis comes from a
run whose code `MANIFEST.md` states was never shared ("`train_pinn.py` was referenced but not
uploaded — so it is not in this package and I cannot verify it"). The reported fit
(c₁ = −3.01, c₂ = 0.07) implies R(x)/R₀ = −1.94 at full discharge — **negative irreversible heat,
thermodynamically forbidden** — and **its RMSE was never reported**.

**To close:** either (a) obtain the code and raw output of that run, or (b) drop the claim and
write the interior minimum's origin as unresolved, naming the entropic hypothesis as leading
candidate. **Recommendation: (b).** MANIFEST.md itself says "LFP cannot test an LCO hypothesis…
it says nothing about whether B0005's specific minimum is entropic."

`[VERIFY: Iontech multi-C-rate run — code, RMSE, raw output, or drop]`

### G-2. Part 6's headline R₀ needs a range, and G2 must be re-scored — **RESOLVED**
**Where:** Part 6 quotes R₀ = 144.78 mΩ, "+0.30 % vs cycle-0 classical", scoring **G2 a HIT**.

**Closed by a 6-seed × 2-weight sweep** (`p6_seeds.py`, `results/p6_seed_summary.log`, 12 runs,
~23 min):

| w_data | R₀ across 6 seeds | deviation | G2 HIT rate |
|---|---|---|---|
| 20 (quoted) | 147.05 ± 6.32 mΩ | +1.87 % [−6.71, +5.87] | **2/6** |
| 200 | 146.50 ± 6.11 mΩ | +1.48 % [−7.11, +3.58] | 1/6 |

**The quoted 144.78 mΩ sits at the favourable edge — 5 of 6 seeds land above it.** As scored in
the notebook G2 is not representative.

**But the project's own selection rule rescues it.** Picking the run with the lowest PDE residual
— the truth-free criterion mandated throughout — selects seed 5 at w = 20, giving 147.93 mΩ
(+2.48 %), **a legitimate HIT**. So the correct statement is: *G2 passes under truth-free
selection, though only 2 of 6 individual seeds satisfy it, and the value originally quoted was at
the favourable extreme.*

**A genuinely favourable finding fell out of this.** Seed 2 is the outlier at both weights
(−6.7 %, −7.1 %), and its PDE residual is **8.8× the median at w = 20 and 4.6× at w = 200** — the
truth-free criterion **flags it correctly**. That is the opposite of Part 5b, where the bad run's
residual sat inside the good range and the criterion was blind. Worth stating as a contrast: the
criterion's reliability is problem-dependent, and on real data here it worked.

Excluding the flagged run, w = 200 gives **148.98 ± 0.65 mΩ** — a 0.4 % spread. The estimator is
stable when it converges; the issue is the 1-in-6 failure, which **replicates Part 5b's rate on
real data**.

**Robust across all 12 runs:** c₁ negative in 12/12 (the broken residual gave +0.023), RMSE
0.132 ± 0.009 K at w = 20 and 0.074 ± 0.003 K at w = 200 against 0.62 K broken. **The shape is
recovered every time; only the amplitude is seed-sensitive.** The chapter's argument is intact.

**Remaining action:** quote R₀ as a range, re-score G2 with the distribution stated, and add the
1-in-6 failure rate wherever the Part 6 accuracy figure appears.

### G-3. University of the Punjab thesis template not supplied
Front matter, chapter numbering, heading levels, equation numbering and citation style are all
provisional. Chapter 6 carries a banner saying so.

`[VERIFY: UoP template]`

---

## Tier 2 — affects specific numbers

### G-4. Reference metadata verified against the publisher record — **RESOLVED**
Every reference in both manuscripts was checked against the publisher's own record on 6 Aug 2026.
One carried a wrong title and several were incomplete:

| Reference | Correction |
|---|---|
| **perez2012** | **The published title contains "cylindrical"** — the source archive's readme omits it. Added pp. 41–50, paper DSCC2012-MOVIC2012-8782, venue detail, doi 10.1115/DSCC2012-MOVIC2012-8782 |
| **lin2013** | Vol. 21, no. 5, pp. 1745–1755, Sep. 2013 and the full eight-author list were all missing |
| **catenaro2021** | Companion article was not cited: *Data in Brief* **35**, art. 106894, doi 10.1016/j.dib.2021.106894 — verified to reference Mendeley `kxsbr4x3j2` specifically, and it lists exactly the eight LFP C-rates used here |
| richardson2015 | doi 10.1109/TSTE.2015.2420375 added |
| forgez2010 | doi 10.1016/j.jpowsour.2009.10.105 added |
| bernardi1985 | already correct |
| saha2007 | repository entry, no DOI issued |

Applied to both manuscripts by `fix_references.py`; citation order re-verified as IEEE-compliant
and both recompiled clean.

### G-5. Cell mass — **RESOLVED, and the flagged discrepancy was an error of mine**
I compared Richardson's cell against the wrong datasheet. The 76 g figure belongs to the
**ANR26650M1-B** (2.6 Ah), a later variant, which is the cell in the multi-rate dataset.
Richardson's data is on the **ANR26650M1**, whose manufacturer datasheet gives:

```
Nominal capacity and voltage      2.3 Ah, 3.3 V
Internal impedance (1kHz AC)      8 mOhm typical
Core cell weight                  70 grams
```

Against the correct figure, Richardson's ρ×V = **72.10 g is +3.0 %**, not −5.1 %, which is what an
*effective* bulk density for a thermal model should look like. There is no discrepancy.

**Sensitivity quantified anyway** (`g5_rhocp_sensitivity.py`, `results/g5_rhocp.log`), because
ρc_p is assumed either way:

| | Richardson (72.10 g) | datasheet (70.00 g) | change |
|---|---|---|---|
| **Quasi-steady core RMSE** | **0.1409 / 0.3366 K** | **identical** | **none** — ρc_p is absent from the relation |
| Classical R_eff (DS2) | 14.3009 mΩ | 14.2414 mΩ | 0.4 % |
| Classical core RMSE (DS2) | 0.8943 K | 0.8300 K | 7 % |
| τ_diff | 1017 s | 987 s | −2.9 % |
| Criterion exponent | −1.18 | **unchanged** | — |
| 5 % crossing | t/τ = 1.37 | t/τ = 1.41 | 2.9 % |

**No reported conclusion changes.** The headline is structurally immune because ρc_p does not
appear in the quasi-steady relation, and the criterion's exponent is unaffected because a uniform
shift in τ moves every plotted point together.

**A third independent anchor fell out of it.** The same datasheet gives 8 mΩ at 1 kHz, 25 °C. That
excludes charge-transfer and diffusion and is taken at 25 °C, so it must sit *below* a DC value at
8 °C — and it does, by about 1.8× against the recovered 14.30 mΩ. Now cited in the paper.

### G-7. A123 datasheet — **RESOLVED**
Consulted directly and now cited as `\bibitem{a123ds}` in the manuscript. Supplies nominal
capacity, core cell weight and 1 kHz impedance, all used above.

### G-6. "B0005 @ 24 °C stratum" is a nominal label
The `.mat` mirror carries no `metadata.csv`; the notebook prints "No metadata — cannot stratify".
The stratum is therefore *all 168 discharge cycles in `B0005.mat`*, not an ambient-selected subset.
No number changes, but the description must be accurate.

`[VERIFY: whether all B0005 cycles are at 24 °C nominal — the CSV mirror's metadata would settle it]`

*(G-7 resolved — moved up beside G-5, which it settled.)*

---

## Tier 3 — noted, low impact

### G-8. Multi-rate internal gradient is modelled, not measured
The timescale criterion in Chapter 7 rests on gradients computed with the verified finite-volume
solver, because that dataset has no core sensor. The solver is the one validated against a real
core in Chapter 6 on the same cell family, which is the strongest available warrant — but it is not
an internal measurement and is not claimed as one. Stated in the chapter.

### G-9. Ambient not logged in the multi-rate dataset
The nominal chamber temperature is used, cross-checked against the plateau at the end of each
cooling tail.

### G-10. `FINDINGS.md` has two items numbered 6
Cosmetic.

### G-11. Part 5b found w_data = 20 optimal; Part 7 used 200
Different problems, but an examiner may connect them. Addressed in Chapter 6 §6.6.1.

---

## Resolved during this work

| Was | Resolution |
|---|---|
| "1.44 % or 1.6 %?" — two rise values for the same RMSE | **1.44 %.** Re-execution prints "0.228 K over a rise of 15.87 K = 1.44 %"; cell 17 confirms 15.87 K as the median rise across 168 cycles. |
| Part 3 numbers were transcriptions from another machine | **Re-executed here**: 24/24 cells, 0 errors, every value reproduced. Tier V1. |
| Was the timescale criterion just an argument? | **Measured** across eight C-rates: ε = 7.48 (t/τ)^−1.18, r² = 0.996. |
| Did any core data leak into Part 7? | One leak found (initial condition), quantified at 0.6 %, removed, all downstream constants re-derived, contaminated results deleted. |
