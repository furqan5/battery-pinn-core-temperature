# PCM melt-front identifiability — findings

Two-day de-risk run before the WUT collaboration meeting. Question:

> Is the PCM melt-front position — equivalently the remaining latent buffer —
> identifiable from external temperature measurement alone?

**Answer: the buffer is identifiable; the front is not.** Those turn out to be
different statements, and the difference is what the proposal has to be
rewritten around.

---

## 1. Verdict, in the language of the brief's decision thresholds

*(References to "the brief's §5" below mean the three decision outcomes in the
session brief, not section 5 of this document.)*

The brief's thresholds for the **strongest** outcome are met on the numbers:
S1 standard deviation 0.10 % of the full buffer (threshold < 10 %), and
three-case discrimination at 7.5× the sensor noise (threshold > 3×).

**But the mechanism is not the one the proposal claims, and by mechanism this is
the brief's second outcome: identifiable by integration, not by front
tracking — so the deck's melt-front-tracking language must change before the
meeting.**

Melted fraction is recoverable because it is a monotone function of the net heat
the cell has accumulated, and accumulated heat shows up in the surface
temperature. It is *not* recoverable from the front's effect on the thermal
resistance of the PCM layer. Five independent lines say so:

| evidence | result |
|---|---|
| Latent heat cut 1000× (no front exists) | buffer still recovered as well: **0.114 %** vs 0.101 % |
| Stage C control, latent-free cell | **more** separation than real PCM: 0.832 K vs 0.751 K |
| Stage E, front-attributable excess separation | **negative at every heating rate**, −0.86σ at 3C to −32.6σ at 19C |
| k_liquid/k_solid swept 1.0 → 0.5 (the front-changes-resistance mechanism) | bound moves only 0.102 % → 0.118 % |
| Fisher information in the mean level of the trace | **99.06 %** — the measurement is one number, not a time signature |

The last one is the Part 7 finding again in a new setting: if 99 % of the
information is in the mean level, a monotone lookup from surface temperature to
melted fraction extracts it, and an inverse PDE solve adds nothing.

Phase change does not merely fail to help. It actively **suppresses**
discriminability — latent heat compresses the temperature excursion, which is
why the latent-free control beats the real PCM at every heating rate and by a
growing margin as melting accelerates.

### What this licenses the proposal to claim

- **Supported:** remaining latent buffer estimated to ~0.1 % of capacity from a
  single surface sensor, *given* the heat load and boundary coefficient.
- **Supported:** PCM exhaustion detected as a distinct event — sensitivity jumps
  **13.3×** within 200 s of exhaustion (4.13 K → 55.0 K per unit melted
  fraction).
- **Supported:** PCM properties identifiable in place from one full melt —
  L_f to 0.015 %, k_pcm to 0.087 %, correlation only 0.18.
- **Not supported:** continuous tracking of melt-front *position* as a
  distinguishable physical quantity. The front is not what is being observed.
- **Not supported:** the claim that phase change creates observability that a
  quasi-steady field lacks. The history-carrying argument is correct in
  principle — the buffer *is* an integral of past heat — but the integral is
  read off the temperature level, not off the front.

### The load-dependence caveat, which decides how much survives

Because the information is accumulated heat, the estimate is only as good as
your knowledge of the heat load and the boundary condition. In the realistic
field case **S4 {f, q‴, h}**, once the assumed history is made consistent with
the trial heat load, **corr(f, h) = 0.9993** and the bound degrades from 0.56 %
to 2.73 %. That is the same near-collinear structure that warned in Part 7
(k vs h at 0.998) before the CRLB-blessed parameter set failed.

And σ = 0.1 K is *instrument* noise. On real hardware the model-form error
dominates and is correlated. Scaling is linear in σ, and correlated residuals
cost effective sample size on top: at a 30 s correlation time the bound inflates
7.8×, at 120 s it inflates 15.5×. Compounding those on S4: 2.73 % × 5 (for
σ = 0.5 K) × 7.75 (30 s correlation) ≈ **106 %** — i.e. no information at all
about the buffer in the realistic field case, on plausible real-hardware
assumptions. **The headline 0.1 % is a best-case number and should not be quoted
at the meeting without the caveat.**

---

## 2. Pre-registered predictions, scored

Registered in `PCM_PREDICTIONS.md`, committed before any Stage B number existed
(git `101bf7a`, ahead of results commit `e1d2b02`, both dated 2026-08-06).
**Three hits, three misses.**

*Hash note:* this public repository is a filtered export of the private working
repository, so hashes differ from the originals (`868b4ef` and `62cbedb` there).
Commit order, messages and dates are unchanged, and those are what carry the
claim — the registration commit precedes the first result commit.

| # | Prediction | Threshold | Actual | |
|---|---|---|---|---|
| Q1 | Melted fraction alone (S1) identifiable | sd < 10 % | **0.101 %** of buffer | ✅ **HIT** |
| Q2 | {f, q‴} poorly conditioned | cond > 10⁴ | **60.5** as computed, **14.6** leak-corrected | ❌ **MISS** — by ~3 orders of magnitude |
| Q3 | Discrimination below 0.5 K during melting | max ΔT < 0.5 K | **0.751 K** worst pair (7.5σ) | ❌ **MISS** |
| Q4 | Sensitivity ≥ 5× larger after full melt | ratio ≥ 5 | **13.3×** | ✅ **HIT** |
| Q5 | Second sensor improves buffer bound < 2× | factor < 2 | **1.997×** on S1 | ✅ **HIT** — but marginal, and **4.81× on S4** |
| Q6 | A high heating rate exists where the front is identifiable | exists | **No rate, in the opposite direction** | ❌ **MISS** |

### On the misses

**Q2** was wrong about the mechanism. I predicted "more heat into more PCM looks
like less heat into less PCM". Over a 600 s window that degeneracy does not bite,
because f sets the *level* of the surface trace and q‴ sets its *slope*, and
those are separable. The degeneracy I should have predicted is the one that
actually appeared — f against **h**, not against q‴ — and it only appeared after
the history leak was fixed. Fixing the leak *improved* S3 (cond 60.5 → 14.6)
while making S4 far worse (corr 0.969 → 0.9993). I had the right worry attached
to the wrong pair of parameters.

**Q3** was wrong because I assumed the surface would be pinned near the melt
point. It is not: RT42 melts over a 38–43 °C range and there is a real gradient
across a 5 mm layer of k = 0.2 W/m·K material, so the surface climbs 7.7 K over
the melt rather than sitting flat. The pinning argument in §1 of the brief —
which I accepted when writing the predictions — is qualitatively wrong for a
technical-grade paraffin in this geometry. It would be closer to right for a
sharp-melting PCM, and the dTm sweep shows that direction (separation 0.751 K at
the physical 4 K width, 0.692 K as the width → 0.25 K) — but even a perfectly
sharp PCM still gives 0.69 K, so pinning was never the binding constraint.

**Q6** was wrong in *direction*, which is the most useful miss of the six. I
predicted the front would become identifiable at high heating rates, where the
mushy zone spans a real temperature range. The opposite happens: the
front-attributable signal is negative at every rate and becomes *more* negative
as melting accelerates, because the latent plateau suppresses exactly the
temperature excursion that carries the information. There is no regime — none
was found across a 40× span of heating rate, t_melt/τ_diff from 29.6 down to
1.3 — in which the front contributes anything an estimator could use.

Q5 is scored a hit on S1 at 1.997×, which is close enough to the threshold to be
luck rather than insight; on S4 the same second sensor buys 4.81×, so the
underlying claim ("the limit is physical pinning, not sensor count") is **not**
supported even where the number passed.

---

## 3. The sentence for the meeting

> We can estimate the remaining latent buffer from a single surface sensor to
> about 0.1 % of capacity and flag PCM exhaustion as a thirteen-fold jump in
> sensitivity — but that information comes from integrated heat, not from the
> melt front, so what we should propose is exhaustion detection and buffer
> estimation by integration, and the twelve months are better spent on making
> that robust to an unknown heat load than on tracking a front the surface
> cannot see.

---

## 4. Which deck claims this supports

**The proposal deck was not on disk, so this is written against the claim-types
named in the brief rather than against your actual slides — send me the deck and
I will mark it up claim by claim.**

Against what the brief describes:

Claims of the form *"we reconstruct the melt front position from external
sensors"* and *"we track the melt front continuously"* are **not supported** and
should be cut or rewritten. The front's position is not observable in this
configuration at any heating rate tested; what varies with the surface trace is
accumulated heat, of which melted fraction is a bookkeeping consequence.

The motivating argument in the deck — that a melting PCM carries history a
quasi-steady field does not, so the Part 7 negative result does not transfer — is
**correct in principle and survives**, but it does not do the work it was meant
to. The buffer really is an integral of past heat and really cannot be expressed
by a quasi-steady relation. Yet the integral is legible in the temperature
*level*, so a lookup table reaches it and the PDE solver is again unnecessary.
The Part 7 pattern repeats: the physics is richer, and the estimator still is not
needed.

Claims of the form *"we detect thermal-runaway precursors via loss of PCM
buffer"* are **supported and are the strongest thing here** — exhaustion is a
sharp, unambiguous, 13× event.

The apparent-heat-capacity method named in the deck is **verified and sound**
(§5 below), with one caveat worth a footnote: implemented literally, pointwise,
it loses 4.7–33 % of the energy while producing a plausible-looking temperature
field. Anyone reproducing this needs the enthalpy-chord or Newton form.

---

## 5. Stage A — solver verification

The identifiability results are only worth as much as the solver. All checks pass.

| check | result | tolerance |
|---|---|---|
| One-phase Stefan (Neumann) front position | 0.537 % error | < 1 % |
| Convergence to sharp-interface limit | 2.81 % → 0.155 % as ΔT_m 4 K → 0.1 K, first order | — |
| Mesh convergence | converged by N = 300 (0.53 %, stable to N = 1200) | — |
| Timestep convergence | 0.537 % at every dt from 4 s to 0.5 s | — |
| Global energy closure, battery configuration | 5.30 × 10⁻¹³ | < 10⁻¹⁰ |
| Latent absorbed = melted mass × L_f | 2.13 × 10⁻¹⁶ | < 10⁻¹⁰ |
| Apparent heat capacity vs enthalpy method | agree to 2.95 × 10⁻¹⁰ K | < 10⁻⁶ |

Residual benchmark error is entirely mushy-zone width, converging linearly — the
correct behaviour, since the analytical solution *is* the ΔT_m → 0 limit.

**Trap 1 confirmed, not assumed.** The literal pointwise apparent heat capacity
method loses 4.7 % to 33 % of the energy depending on mushy width and timestep,
worst at narrow ΔT_m and coarse dt, exactly as warned. The temperature field
looks entirely plausible throughout; only the energy closure catches it.

---

## 6. Regime and parameters

| quantity | value | |
|---|---|---|
| t_melt / τ_diff | 29.6 (nominal 3C) → 1.3 (19C equivalent) | swept |
| Stefan number, c_p ΔT_m / L_f | 0.048 | latent-dominated |
| Biot number, PCM layer | 0.25 | |
| Latent buffer | 4274 J (25.9 g RT42) | |
| Melt duration at 3C-equivalent | 5618 s | |

Cell: HAKADI/Selian 26700 sodium-ion, 3.5 Ah, 84.2 g, ≤20 mΩ at 1 kHz — Bischof
et al., *Journal of Power Sources Advances* **27** (2024) 100148. PCM: Rubitherm
RT42, L_f 165 kJ/kg, 311.5–315.5 K melting range, k 0.2 W/m·K, c_p 2000 J/kg·K —
Hammoodi et al., *Scientific Reports* **15** (2025) 31308, Table 1, cross-checked
against the Rubitherm RT product table (42 °C, 165 kJ/kg). Cell radial
conductivity from Bhundiya, Hunt & Drolen, TFAWS 2018 (0.43 ± 0.07 W/m·K for
18650 Li-ion, transferred to Na-ion — tagged as an estimate). Full provenance
tags in `pcm_params.py`.

**Two provenance points that cut against the result being too pessimistic:**

1. Rubitherm's 165 kJ/kg is latent *plus* sensible by their convention; used
   here as pure latent, following the cited paper. This **overstates** the
   buffer, which is generous to the proposal.
2. The cited source gives one conductivity, not separate solid and liquid
   values. Rather than silently set k_l = k_s — which would have manufactured
   the null result — the ratio was swept 1.0 → 0.5 and moves the bound by 16 %.

Neither assumption is doing the work. Nor is the mushy-zone width: from 0.25 K
to 2 K the separation is flat at 0.69 K.

---

## 7. Limitations — what could still overturn this

**Natural convection in the molten PCM is the real one, and it was tested.**
A 1-D model cannot carry buoyancy-driven flow in the liquid region, which in a
real melting paraffin raises effective heat transfer substantially. The standard
proxy is an enhanced conductivity in the liquid only, and it is the one direction
the E2 sweep did not cover (E2 went *downward*, k_l/k_s = 0.5–1.0). Running it
upward:

| k_liquid/k_solid | 1 (molecular) | 2 | 3 | 5 | 10 |
|---|---|---|---|---|---|
| buffer sd (% of buffer) | 0.102 | 0.083 | 0.078 | 0.081 | 0.089 |
| worst separation (K) | 0.751 | 0.931 | 1.052 | 1.012 | 0.922 |

Convective enhancement helps modestly and **non-monotonically**, peaking near
k_l/k_s ≈ 3, and at best moves the separation from 0.75 K to 1.05 K and the bound
from 0.102 % to 0.078 %. Same order of magnitude, same conclusion. It does not
rescue front tracking. (`results/pcm_convection_check.json`.)

Remaining limitations, none of which were tested:

- **1-D radial.** No axial gradient, no tab heating, no end effects. The in-plane
  tab-heat localisation idea in the proposal is a genuinely different problem and
  nothing here speaks to it either way.
- **Constant density.** No volume change on melting. The real 880 → 760 kg/m³
  contraction is what drives the convection above and also causes void formation
  at the can wall, which would change the contact resistance — plausibly a larger
  effect than anything modelled here, and not capturable in 1-D.
- **Sustained constant heat load.** Represents repeated cycling at an average
  dissipation, not a single 3C discharge (which lasts only 20 minutes and would
  not melt the PCM on its own).
- **Model assumed exact, noise assumed white.** The L2 audit quantifies how much
  this flatters the result; it is the largest single caveat on the headline number.

---

## 8. What I would do with the twelve months

Not a deliverable, but it follows directly and you will be asked.

The PCM inverse problem as posed is not where the value is. Three options that
survive this result:

1. **Exhaustion detection and buffer estimation by integration**, made robust to
   an unknown and time-varying heat load. The 0.9993 f–h collinearity is the
   real technical problem, and it is a genuine one — worth twelve months, and
   honestly framed.
2. **In-situ PCM property identification.** S2 came out clean (cond 6.1, L_f to
   0.015 %, correlation 0.18) and is the least-oversold result in this study.
   Measuring how a PCM's properties drift with thermal cycling is a real
   question that this configuration genuinely can answer.
3. **In-plane tab-heat localisation**, as the brief already anticipated — not
   tested here, and the one remaining place where a spatial inverse problem
   might have something a lookup table does not.

Reproduce with: `pcm_stage_a_verify.py`, `pcm_stage_b_crlb.py`,
`pcm_stage_b_audit.py`, `pcm_stage_cd.py`, `pcm_stage_e_regime.py`, or
`PCM_Identifiability.ipynb` top to bottom.
