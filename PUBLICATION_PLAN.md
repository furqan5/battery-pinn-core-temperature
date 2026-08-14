# Publication plan — two papers, zero publication fee

Verified 2026-08-15. Supersedes the cost section of `paper/SUBMISSION_NOTES.md`,
whose figures came from IEEE's 2020 schedule.

> **On the file paths below.** The manuscript sources (`paper/`, `paper2/` and
> `paper/SUBMISSION_NOTES.md`) are not in this public repository yet — they are
> held until the author list is final. Paths are kept so the plan stays readable,
> and the files will appear at them when the manuscripts are added.

---

## 1. The filter rule that answers "no fee at all"

There is one structural rule and it removes most of the guesswork:

> **A hybrid journal is free to publish in on the subscription route. A fully
> open-access journal always charges.**

Hybrid means the journal is subscription-funded and *offers* open access as a
paid option. Decline the option and the cost is zero. Fully open access means
the article-processing charge is the business model and there is no free route.

So the filter is: **choose hybrid, decline OA, and then check only one thing —
whether the publisher levies page charges on top.** IEEE societies do. Elsevier
and IOP generally do not.

### Excluded outright — these always charge

IEEE Access, IEEE Open Journal of Instrumentation and Measurement, every MDPI
title (Batteries, Energies, Sensors), Frontiers, Scientific Reports,
eTransportation. No amount of venue-fit argues around a mandatory APC.

---

## 2. IEEE TIM — verified fee position for 2026

| Item | Value | Status |
|---|---|---|
| Submission fee | none | verified |
| Mandatory publication charge, subscription route | **none** | verified — TIM is hybrid, OA is optional |
| **Free page allowance, regular paper** | **8 pages** | verified |
| Minimum length, regular paper | 5 pages | from IMS author information |
| Overlength charge, non-member | $265 per page or part page | verified |
| **Overlength charge, IMS member** | **$220 per page** | verified — this is new information, the old notes had only the $265 rate |
| Voluntary page charge | optional — decline it | verified |
| Optional OA APC | **$2,645** | verified; the notes' $2,045 is out of date. Irrelevant if OA is declined |
| Lower-middle-income discount | 25–50 % on overlength, by country GDP | verified; Pakistan is lower-middle-income for FY2026 |

**Conclusion: at 8 pages or fewer, publishing in TIM costs exactly $0.** The
judgment in `SUBMISSION_NOTES.md` was right and is now verified against 2026
figures rather than 2020 ones.

Two notes worth acting on:

- **IMS membership pays for itself only if you overrun.** At $220 vs $265 per
  page the saving is $45/page, so membership is not worth buying as insurance
  unless the paper is already over. Keep both papers under 8 pages instead.
- The lower-middle-income discount applies to *overlength* charges. It is
  insurance, not a plan.

---

## 3. Page-count estimates — both papers fit, one may be too short

Estimated with `scratchpad/pagecount.py` (IEEEtran two-column heuristic: ~975
body words per page, a single-column float ~0.33 page, an equation ~28 words).
**This is an estimate, not a compile.** Neither paper has been compiled; there is
no LaTeX toolchain on this machine.

| | body words | figs | tables | eqs | refs | estimated pages | fee |
|---|---|---|---|---|---|---|---|
| Paper 1 — core temperature validation | 3329 | 3 | 3 | 7 | 8 | **5.7** | $0 |
| Paper 2 — identifiability-gated inverse estimation | 2393 | 3 | 3 | 3 | 5 | **4.5** | $0 |

Both are comfortably inside the 8-page allowance. Neither needs trimming, so
`SUBMISSION_NOTES.md` item 2 ("trim Section III or VII first") is not required.

**The live issue is the opposite one.** Paper 2 estimates at 4.5 pages against a
stated 5-page minimum for regular papers. That is inside the estimator's error
bar, but it is the wrong side of a threshold, and a paper submitted below the
minimum risks being redirected to a short-paper category with a different review
track. There is ready material to expand it with, all of it already established
and currently held back only for length:

- The full CRLB condition-number table across parameter sets, not just the
  reported worst-relative-bound figures.
- The profile-likelihood construction for the pooled convection coefficient,
  currently compressed.
- The white-noise diagnostic that detects the series has been filtered — this is
  a genuinely reusable contribution and is worth its own subsection with the
  method stated fully, rather than a sentence in the abstract.
- The per-cycle quality-gate definitions, so the "all 168 cycles inside
  pre-registered gates" claim is checkable by a reader.

Expanding to 6 pages is still free and removes the threshold risk.

---

## 4. Zero-fee fallback ladder

TIM's acceptance rate makes a rejection likely enough to plan for. Both papers
target TIM first; the fallbacks below are deliberately **different for each
paper** so that one rejection does not queue both behind the same second choice.

### Paper 1 — core temperature, estimator comparison, timescale criterion

| Rank | Venue | Fee on subscription route | Fit |
|---|---|---|---|
| 1 | **IEEE TIM** | $0 at ≤8 pp | measurement validated against a reference sensor; a negative result is publishable here |
| 2 | **Measurement Science and Technology** (IOP) | $0 — hybrid `[verify per-journal page charges]` | instrumentation and measurement science; comfortable with method-limitation papers |
| 3 | **Measurement** (Elsevier) | $0 — hybrid | broad measurement scope, high throughput |
| 4 | **Journal of Energy Storage** (Elsevier) | $0 — hybrid, 24-month embargo | good scope fit, but the audience expects methods that win |

### Paper 2 — CRLB-gated model selection, heat-generation recovery

| Rank | Venue | Fee on subscription route | Fit |
|---|---|---|---|
| 1 | **IEEE TIM** | $0 at ≤8 pp | parameter estimation from a measurement channel, plus the filtered-noise diagnostic |
| 2 | **Journal of Power Sources** (Elsevier) | $0 — hybrid | heat generation vs SOC with an EIS cross-check is squarely in scope |
| 3 | **Applied Thermal Engineering** (Elsevier) | $0 — hybrid | thermal inverse problem |
| 4 | **IEEE Trans. Industrial Electronics** | $0, but check page allowance | estimation and identifiability for a power-electronics audience |

`[VERIFY at submission]` Elsevier and IOP page-charge policies per title. The
subscription route being free is a general publisher policy; a specific journal
levying page or colour charges on top is the exception to check. Colour charges
have largely been abolished but confirm rather than assume.

---

## 5. Preprints — do this first, it is free and it helps

Post both to **arXiv** (`eess.SY` or `physics.ins-det`) before or at submission.
IEEE permits preprints on arXiv and permits updating them to the accepted
version. Three reasons this matters here:

1. It establishes priority on the negative result, which is the part most likely
   to be scooped by someone reporting the positive half.
2. It gives a citable anchor so **paper 2 can cite paper 1** rather than
   describing its result inline. Paper 2's identifiability framing reads better
   standing on paper 1.
3. It gives the code release something to point at, and gives the papers
   something to point back to — see `CODE_RELEASE.md`.

**Submission order:** paper 1 first. Post paper 2's preprint at the same time but
submit it two to four weeks later, after paper 1 has an arXiv identifier to cite.

---

## 6. Blocking items before either submission

Everything substantive is closed. What remains is mechanical.

1. **Compile both papers.** No LaTeX toolchain locally. Use Overleaf's *IEEE
   Journal Paper* template, which supplies `IEEEtran.cls`. This is the only item
   that can change the fee conclusion, because the page estimate is a heuristic.
2. **Confirm the real page count** against the 8-page ceiling and the 5-page
   floor. Expand paper 2 if it lands under 5.
3. **ORCID.** Required by IEEE at submission. Free to register.
4. **Decline the voluntary page charge and decline open access** in the IEEE
   submission workflow. These are opt-in prompts; accepting either is how a
   free submission stops being free.
5. **Sign the overlength agreement only if actually over 8 pages.** Do not sign
   it pre-emptively.

### Already closed — the old notes are stale on these

- `SUBMISSION_NOTES.md` says to switch `\documentclass` from `onecolumn` to
  `[journal]`. Both papers already declare `\documentclass[journal]{IEEEtran}`,
  which is two-column submission format. Nothing to do.
- G-4 (reference metadata and DOIs) and G-5 (cell mass) are both marked
  **RESOLVED** in `GAPS.md`. The notes list them as outstanding.
- G-1 (the unverifiable Iontech reversible-heat claim) is handled the way
  `GAPS.md` recommends: paper 2 §"Is the minimum irreversible or entropic?"
  reports the multi-rate separation, states that the minimum survives it, and
  then explicitly bounds the claim — the reversible term is undetermined to
  better than a factor of two, and the test "weakens the entropic explanation
  without settling it." No unsupported claim is made.
- G-6 (the "24 °C stratum" is a nominal label) is handled: paper 2 describes the
  data as "168 discharge cycles of cell B0005" with no stratification claim.

---

## 7. What this costs, end to end

| | cost |
|---|---|
| arXiv preprints, both papers | $0 |
| IEEE TIM submission, both papers | $0 |
| IEEE TIM publication at ≤8 pages, both papers | $0 |
| ORCID | $0 |
| Overleaf (free tier is sufficient) | $0 |
| **Total** | **$0** |

The only route to a non-zero cost is accepting an optional prompt in the IEEE
workflow, or overrunning 8 pages. Both are under your control.
