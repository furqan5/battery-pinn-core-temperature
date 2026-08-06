# GAPS — every `[VERIFY: …]` marker, and what closes it

Collected 6 Aug 2026. Ordered by how much damage each would do if an examiner found it first.

---

## Tier 1 — must be resolved before submission

### G-1. The reversible-heat claim has no verifiable source
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

### G-4. Reference metadata taken from code archives, not articles
Richardson & Howey (2015) and Forgez et al. (2010) are cited from the README and readme-text of the
code releases that use them, not from the publisher record.

`[VERIFY: page numbers, volume, DOI for both]`

### G-5. Cell mass discrepancy, ~5 %
The manufacturer sheet inside the previously unused multi-rate dataset gives **0.076 kg** for the
A123 ANR26650m1-b. Part 7 uses ρ = 2107 kg m⁻³ with V = 3.4219×10⁻⁵ m³, implying **72.1 g** — 5.1 %
lower. ρc_p propagates directly into every recovered coefficient.

**To close:** decide which mass applies to the Richardson cell (2.3 Ah variant) versus the
multi-rate cell (2.5 Ah "-b" variant) and state the assumption. This is cheap and converts part of
ρc_p from **(b)** toward **(a)**.

`[VERIFY: mass and capacity of the exact cell in the Richardson archive]`

### G-6. "B0005 @ 24 °C stratum" is a nominal label
The `.mat` mirror carries no `metadata.csv`; the notebook prints "No metadata — cannot stratify".
The stratum is therefore *all 168 discharge cycles in `B0005.mat`*, not an ambient-selected subset.
No number changes, but the description must be accurate.

`[VERIFY: whether all B0005 cycles are at 24 °C nominal — the CSV mirror's metadata would settle it]`

### G-7. A123 datasheet not consulted directly
Nominal capacity, DC resistance and dimensions come from the Richardson archive and the multi-rate
dataset's summary sheet, not the manufacturer document.

`[VERIFY: A123 ANR26650M1 / -B datasheet]`

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
