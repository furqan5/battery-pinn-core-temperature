# FINDINGS — Part 7: internal temperature validated against a measured core

## The headline, stated the way the evidence supports it

Fitting to the **surface trace only**, the **measured** core temperature of an A123 26650 LFP
cell is reproduced to:

| method | DS2 core RMSE | DS1 core RMSE |
|---|---|---|
| **quasi-steady one-liner** `T_core = T_surf + (Bi/2)(T_surf − T_inf)` | **0.1409 K** | **0.3366 K** |
| classical transient finite-volume | 0.8943 K | 1.0544 K |
| **inverse PINN** (the method under test) | **0.6135 K** | 0.5024 K |

against measured core-to-surface gradients peaking at **6.54 K** (DS2) and
**7.13 K** (DS1).

**So the internal field IS reconstructable from surface data to a fraction of a kelvin — but
the PINN is not what does it.** A single algebraic line beats it by 4.4× on DS2.
This is a negative result for the method and the most useful thing in the notebook.

The PINN still passes its pre-registered target (P5, < 1 K), and every reported number is
leak-free. It is simply not the best tool for this problem, and saying so is the point.

---

## 1. Why the one-liner wins

Because **the problem is quasi-steady**. The measured ratio (core − surface)/(surface − ambient)
over the hot part of each record is

* DS2: **0.6161 ± 0.0126** (2.0 % variability)
* DS1: **0.6218 ± 0.0103** (1.7 % variability)

Constant to about 2 % across a 3541 s drive cycle with 30 A peaks. The radial profile shape
never really changes; only its amplitude does. One number therefore captures the entire spatial
structure, and there is nothing dynamic left for a PDE solver to add.

That constant is the steady analytic result Bi/2 = 0.6071 (DS2), which matches the measured
0.6161 to 1.5 %. The relation quoted in the
project's own plausibility gate turns out to be a better *predictor* than the machinery it was
meant to check.

**The comparison is fair.** The one-liner anchors on the measured surface while the PINN
re-predicts it, so the PINN's surface error could have been inherited. Removing that — taking
only the PINN's *gradient* and adding it to the measured surface — still gives 0.5990 K
against the one-liner's 0.1409 K. The gap is in the gradient itself, not in the anchoring.

**What the transient models are still for.** The one-liner needs Bi, and Bi is *not* free: a
±15 % error in it costs 5–10× in core RMSE. Bi = hR/k came from a **transient** surface-only fit
on the other record. So the honest division of labour is: **use a transient model to identify Bi;
use algebra to predict.** The PINN adds nothing to either half.

---

## 2. Scored predictions, misses included

| # | Prediction | Outcome | Observed |
|---|---|---|---|
| P1 | core−surface max between 0.5 and 5 K | **MISS** | 7.13 K (DS1), 6.54 K (DS2) — 43 % above the interval |
| P2 | classical surface RMSE ≤ 0.3 K | **HIT** | 0.220 K (DS1), 0.213 K (DS2) with the measured source |
| P3 | {R₀, shape, h} not identifiable, worst rel sd > 25 % | **MISS** | 1.3–10.6 % at shape orders 1–3; only order 4 breaches 25 % |
| P4 | forward PINN within 5 % on centre-surface difference | **HIT** | 4.66 % (DS2), 1.39 % (DS1) — DS2 close to the line |
| P5 | **predicted core < 1 K RMSE** | **HIT** | **0.6135 K** (order 0), 0.5710 K (order 1) |
| P6 | core error ≥ 2× surface-fit error | **HIT** | ratio 3.11 (order 0), 3.15 (order 1) |

Four hits, two misses. Plus **three predictions I made mid-session and got wrong**, logged
because they cost real debugging:

1. **"Fourier time features are needed for the broadband source."** Backwards. On a controlled
   problem they made things monotonically worse — 0.086 % error at n_freq = 0 rising to 4.88 %
   at n_freq = 64.
2. **"The PINN's core error will be worse than the classical control's, because its R_eff is
   9 % low."** It was 31 % *better*. The mechanical chain — weak source ⇒ small gradient ⇒ bad
   core — simply did not hold.
3. **"Centre-weighted generation explains the low R_eff and the good core together."** Refuted:
   R_eff moved the wrong way (14.30 → 14.36 mΩ) as the weighting increased.

---

## 3. Full results

| | R_eff | surface RMSE | **core RMSE** |
|---|---|---|---|
| **DS2** | | | |
| quasi-steady one-liner | — | — | **0.1409 K** (max 0.3742) |
| classical FV | 14.301 mΩ | 0.4652 K | 0.8943 K |
| inverse PINN, order 0 | 13.108 mΩ | 0.1974 K | 0.6135 K (max 1.9093) |
| inverse PINN, order 1 | 13.058 mΩ | 0.1815 K | 0.5710 K (max 1.4981) |
| **DS1** (roles swapped) | | | |
| quasi-steady one-liner | — | — | **0.3366 K** |
| classical FV | 13.649 mΩ | 0.5613 K | 1.0544 K |
| inverse PINN (order 0) | 13.265 mΩ | 0.1238 K | 0.5024 K |

Across 6 seeds (order 0, DS2): core RMSE sd 0.0086 K, R_eff spread
0.489 %, **non-convergence 0 %**. Selection was on the PDE + BC residual over a
fixed 20 000-point collocation set, never on proximity to the core.

**Where the PINN's error lives.** It is not uniform in time: the **first 10 minutes is 17 % of
the record but 56 % of the squared error** (RMSE 1.1146 K there,
0.4470 K afterwards). The start-up transient has the sharpest radial gradients and the
model assumes a perfectly uniform initial field — the one regime where "quasi-steady" is false.

**Trap 5.1 control.** Removing the I(t)² factor costs 2.7× on the surface fit and 1.8× on the
core. The trap would have been caught here.

---

## 4. The dominant lever: k

Refitting R_eff at each k against the same surface data:

| k (W/m/K) | Bi | surface RMSE | core RMSE |
|---|---|---|---|
| 0.3000 | 1.594 | 0.5433 | 1.8209 |
| 0.3500 | 1.366 | 0.4961 | 1.0998 |
| 0.3937 **(ours)** | 1.214 | 0.4652 | 0.8943 |
| 0.4040 *(Richardson)* | 1.183 | 0.4590 | 0.9027 |
| 0.4500 | 1.062 | 0.4355 | 1.0941 |
| 0.5500 | 0.869 | 0.4007 | 1.6983 |

The surface fit improves **monotonically** with k while the core error has a minimum near
k ≈ 0.39 and roughly doubles either side. At k = 0.55 the surface fit is the best in the table
and the core prediction is nearly twice as bad.

The core prediction is essentially a function of k; surface data cannot identify k (CRLB
correlation with h is 0.998). **So the headline rests on an input, not an output.** Our
leak-free k = 0.39375, obtained from a surface-only fit on the *other* record, lands within
about 1 % of the core-optimal value — either good transfer or luck, and one record cannot tell
the difference.

---

## 5. What is established

1. **The internal temperature of this cell is recoverable from surface data to
   0.14–0.34 K**, verified against a real core thermocouple rather than a model. That
   is the gap this project existed to close, and it is closed.
2. **The mechanism is quasi-steadiness**, demonstrated directly: the gradient/rise ratio is
   constant to ~2 % across both records.
3. **No core data entered any fit.** h and k came from a *surface-only* fit on the *other*
   record. An earlier version leaked the core through the initial condition
   `0.5(T_s[0]+T_c[0])`; it was found by audit, quantified (0.6 % on core RMSE), removed, and all
   downstream constants re-derived.
4. **The inverse PINN works** — it passes its forward gate, converges on 6/6 seeds with 0.8 %
   parameter spread, and meets its < 1 K target — **and is still the worst of the three methods
   tried.**

## 6. What is NOT established

1. **This is one cell, two drive cycles, one ambient (~8 °C).** Nothing here speaks to
   temperature dependence, ageing, or cell-to-cell variation.
2. **Bi is an input.** Everything rests on it and it is not identifiable from a single surface
   channel.
3. **R_eff is a blend, not an internal resistance.** The entropic term is off in the headline.
   mean(|Q_rev|) is 57–63 % of mean ohmic, but the **signed** mean is only −1.4 % because the
   term alternates with charge and discharge. Switching it on (Forgez et al. 2010, −0.5 mV/K at
   50 % SOC) moves the core RMSE by 0.4 %. Real, but immaterial here.
4. **No shape R(x) was recovered.** The DoD window spans 11–14 % of the axis. Order 1 is a
   narrow-band interpolant, not a curve over [0,1].
5. **The PINN is not solving the PDE it claims to solve, and its accuracy comes from that.**
   This is the sharpest result in the notebook and it is not favourable. Take the PINN's own
   recovered R_eff and put it through the exact solver:

   | | surface RMSE | core RMSE |
   |---|---|---|
   | PINN's own field | 0.1974 K | 0.6135 K |
   | same R_eff, physics enforced exactly | 0.8395 K | 1.5238 K |
   | **difference = PDE violation** | **0.8029 K** | **1.0570 K** |

   The PINN's field departs from a true solution by ~0.8 K on the surface —
   larger than its own core error. Solved correctly, its recovered parameter gives
   1.52 K, *worse* than the classical fit's 0.89 K. So the PINN is
   functioning as a flexible curve-fitter with a physics-flavoured regulariser, and its good core
   number is substantially a consequence of bending the physics, not obeying it.

   What it is bending toward is real, though: the classical residual is 100 % model-form error,
   autocorrelated over 423 s (≈ the 405 s thermal time constant), correlating most with *time*
   (r = 0.51) rather than I² (r = 0.10). A smooth neural deviation field absorbs that drift
   non-parametrically. It buys accuracy without insight, and costs a 9 % bias in R_eff against
   two independent anchors.

   The advantage over the classical solver does **replicate** (core ratio 0.69 on DS2, 0.48 on
   DS1), so it is a property of the method rather than one record's quirk. But "reliably better
   by breaking the constraint" is a weaker claim than it first appears, and it does not survive
   comparison with the one-liner.
6. **A plain full-field PINN does not work here at all** — it failed the forward gate at 39.2 %.
   Everything reported depends on the base-subtracted formulation.

---

## 7. Method lessons

1. **Compare against the trivial baseline before claiming a method works.** The whole PINN
   apparatus is beaten by one line of algebra. Without B0/B1/B2 in the table, this notebook would
   have reported a 0.61 K success and buried a 0.14 K result that needed no network.
2. **Report the parameter correlation matrix, not just marginal standard deviations.** The CRLB
   rated {h, k, ρc_p} identifiable at 3.9 %; the fit then drove k 1140 % away, ~300× outside
   the CRLB ellipse, for 4 K of core error. The 0.99 k–ρc_p correlation was the warning.
3. **A better surface fit can mean a worse internal answer.** Stage B [B3] had the best surface
   RMSE of any model and the worst core error.
4. **Fourier time features can hurt badly.** Monotonically worse with n_freq on a controlled
   problem.
5. **Write heat sources sign-explicitly.** `abs(I*(V-U))` hides the current convention; the
   explicit form exposed that I > 0 means charge here, via a 55 % negative-heat fraction.
6. **Test discretisations against closed-form solutions at machine precision.** Two silent bugs
   — a boundary source broadcast to every cell, and an O(dr/4R) wall-conductance bias — were
   caught only because the analytic steady state is reproduced to ~1e-13.

---

## 8. What the next experiment should be

**Stop improving the estimator; go after Bi.** Everything now rests on it and nothing measures
it.

1. **Joint fit across both records with a shared k and per-record h.** DS1 and DS2 have different
   current statistics and DS1 has a full relaxation tail. That may break the k–h collinearity
   (0.998) that defeats a single record. Costs nothing — the data is already here.
2. **Exploit the relaxation tail deliberately.** With the current off, the decay *shape* — not
   its amplitude — constrains k independently of any source model. DS1 already contains one.
3. **Test the quasi-steady relation's limits.** It should fail when the forcing timescale
   approaches the diffusion time (~1000 s). A step or a fast pulse would find that boundary, and
   *that* is the regime where a transient solver would finally earn its place.

If the quasi-steady relation holds up under (3), the practical conclusion for battery thermal
management is that core-temperature estimation on cells like this needs a good Bi and a
thermocouple — not a neural network.
