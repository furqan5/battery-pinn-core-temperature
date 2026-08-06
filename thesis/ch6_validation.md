# Chapter 6 — Validation against an internal measurement

> **Formatting provisional.** Chapter numbering, heading levels, equation numbering and
> reference style follow a conventional engineering-thesis convention pending the University of
> the Punjab template. `[VERIFY: UoP thesis template — front matter, chapter numbering, citation
> style]`

Provenance labels used throughout: **(a)** verified against a named source or an analytic result;
**(b)** engineering estimate with stated assumptions; **(c)** inference drawn in this work.

---

## 6.1 The claim this chapter tests

Every result in Chapters 3 to 5 shares one weakness. The reconstructed internal temperature field
was scored against a finite-difference model, against a classical estimator, and against its own
convergence diagnostics, but never against a thermometer inside a cell. A method validated only
against another model has not been validated as a measurement. On a sealed 18650 there is no way
to do better, because the cell has no core sensor.

This chapter closes that gap. The experiment is narrow and its statement is short: fit the
estimator to a **surface** temperature trace alone, predict the **core** temperature, and score
the prediction against a core thermocouple that the fit never saw. The deliverable is an error in
kelvin.

The answer is that the internal temperature of an A123 26650 cell is recoverable from surface data
to a few tenths of a kelvin, and that the inverse physics-informed neural network is not what
recovers it. A one-line algebraic relation does the same job about four times more accurately. The
negative result is the more useful half, and Section 6.7 gives its mechanism.

---

## 6.2 Cell, instrumentation and data provenance

The experiment requires a cell with simultaneous core and surface thermometry. That rules out the
18650 geometry used in Chapters 3 and 6 of the earlier work and requires a larger format with an
instrumented centre.

The data are two hybrid-electric drive cycles measured on an **A123 ANR26650M1** lithium iron
phosphate cell, 2.3 A h nominal, published by Richardson and Howey with the code accompanying
their impedance-based estimator [1]. Each record carries time, surface temperature, core
temperature, ambient temperature, current and terminal voltage.

Channel identities were not inferred from filenames. They were read from the reference
implementation's own variable assignment (`MainScript.m`, lines 136–139), which maps columns one
to four of the temperature array to time, surface, core and ambient respectively **(a)**. This
matters more than it sounds: mistaking the core channel for the surface channel would invert the
entire experiment and would still produce plausible-looking numbers.

**Table 6.1** — Record summary. Values re-derived in this work from the source archive.

| Quantity | Record 1 | Record 2 |
|---|---|---|
| Duration | 5981.8 s | 3543.1 s |
| Samples (raw) | 5439 at 1.1 s | 3222 at 1.1 s |
| Surface temperature range | 8.13 – 19.53 °C | 8.20 – 18.55 °C |
| Core temperature range | 8.24 – 26.59 °C | 8.24 – 25.03 °C |
| Ambient (mean) | 7.92 °C | 7.88 °C |
| Current range | −23 to +30 A | −23 to +30 A |
| **Peak core − surface** | **7.131 K** | **6.543 K** |
| Coulomb-counted depth-of-discharge window | 0.138 | 0.112 |

Two features of Table 6.1 govern everything that follows. The peak core-to-surface difference is
**7.13 K**, roughly forty times the 0.18 K bow available in the slab study of Chapter 3 **(a)**,
because the cell runs at an 8 °C ambient with 30 A peaks (about 13 C on a 2.3 A h cell). There is
a real spatial signal to reconstruct. Against that, the depth-of-discharge window spans only 11 to
14 per cent of the axis, because these are charge-sustaining hybrid cycles. Recovering a
heat-generation *shape* across the full state-of-charge range would therefore extrapolate over
86 per cent of a domain with no observations behind it, and Section 6.5 restricts the source model
accordingly.

### 6.2.1 The current sign convention, established rather than assumed

The reference implementation computes irreversible heat as `abs(I*(V-3.3))`. The absolute value
conceals the sign convention, so it was re-derived here. Writing the Bernardi irreversible term
sign-explicitly as

$$Q_\mathrm{irr} = I_\mathrm{dis}\,(U_\mathrm{oc} - V)$$

and evaluating it with the discharge-positive convention gave **negative heat in 55 per cent of
samples**, which is thermodynamically impossible. Three independent checks resolve it **(a)**:
the correlation between current and voltage over samples with |I| > 1 A is **+0.90**, meaning
voltage rises when current is positive, which happens only on charge; regressing V on I gives a
slope of **+13.9 mΩ**, consistent with the direct-current resistance of this cell at 8 °C; and the
record's energy balance closes only in that orientation, requiring a mean source of 1.226 W
against the 1.197 W the corrected expression supplies, a 2.4 per cent agreement.

In this dataset **I > 0 denotes charge**. The irreversible heat is therefore

$$Q_\mathrm{irr}(t) = I(t)\,\bigl[V(t) - U_\mathrm{oc}\bigr] \ \ge 0 .$$

This is a small point with a large moral, and it belongs in the record: an `abs()` that makes a
quantity look well-behaved also makes a sign error unfindable.

---

## 6.3 Hold-out protocol

The core channel is the test set. It is not used for fitting, for initialisation, for early
stopping, for hyperparameter choice, or for selecting among random seeds. It is read once, after
every estimator is frozen, to compute the reported errors.

Two thermal parameters cannot be identified from a single surface trace (Chapter 4 establishes this
quantitatively), so they must be fixed from outside the record being fitted. The published values
from [1] were **not** used, because they were identified on drive cycle 1 using *both*
thermocouples; adopting them would import core-derived information into a supposedly core-blind
fit. Instead, each record is fitted using the convection coefficient and radial conductivity
obtained from a **surface-only** fit on the *other* record:

**Table 6.2** — Cross-record parameter transfer. Neither column uses any core data.

| Fitted record | h (W m⁻²K⁻¹) | k (W m⁻¹K⁻¹) | Bi = hR/k | Source of h, k |
|---|---|---|---|---|
| Record 2 | 37.0607 | 0.393747 | 1.2142 | surface-only fit on record 1 |
| Record 1 | 37.2846 | 0.418697 | 1.1487 | surface-only fit on record 2 |

These land within 3.3 per cent of the published values of 39.3 and 0.404 **(a)** [1], which is
reassuring, but the published pair is not used anywhere in this chapter.

### 6.3.1 A leak found by the protocol, and what it cost

The first complete run produced a core error better than the classical control, which triggered
the project rule that an unexpectedly good result is audited before it is reported. The audit
found a leak. The initial condition had been set to the mean of the first surface and first core
sample,

$$T(r,0) = \tfrac{1}{2}\bigl[T_\mathrm{surf}(0) + T_\mathrm{core}(0)\bigr],$$

which reads the held-out channel. The cell is nearly isothermal at t = 0, where the measured
core-to-surface difference is 0.047 K on record 2, so the leak is small. It is not zero.
Quantified with the classical estimator, where the arithmetic is transparent, removing it moves
the core RMSE from 0.8895 K to **0.8950 K**, a 0.6 per cent effect, and R_eff by 0.02 per cent.

The consequence that mattered was downstream. The convection coefficient and conductivity in
Table 6.2 come from a surface-only fit that had used the *same* contaminated initial condition, so
the "independent" parameters carried the leak too. Re-deriving them cleanly moved k from 0.390843
to **0.393747**, a 0.74 per cent change. The initial condition is now the first surface sample
alone, every downstream constant was recomputed, and the contaminated results were deleted rather
than adjusted.

The episode is reported in full because the alternative, quietly patching a 0.6 per cent effect,
would have left an unquantified dependency in the headline. It also shows the hold-out rule doing
the work it exists for.

---

## 6.4 Numerical foundation and the forward gate

### 6.4.1 Governing equation

Radial conduction in a cylinder with volumetric generation:

$$\rho c_p \frac{\partial T}{\partial t}
= \frac{k}{r}\frac{\partial}{\partial r}\!\left(r \frac{\partial T}{\partial r}\right) + q'''(t)$$

with a convective outer boundary and symmetry at the axis:

$$-k \left.\frac{\partial T}{\partial r}\right|_{r=R} = h\bigl[T(R,t) - T_\infty\bigr],
\qquad \left.\frac{\partial T}{\partial r}\right|_{r=0} = 0 .$$

In words: stored heat per unit volume on the left is supplied by conduction from neighbouring
shells plus internal generation. The 1/r factor distinguishes the cylinder from the slab of
Chapter 3. An outward shell has more area than the one inside it, so the same temperature gradient
carries more power outward as r grows. The second boundary condition is less a physical constraint
than a statement that the axis has no neighbour on its inner side.

### 6.4.2 Discretisation and its verification

The truth model is a conservative finite-volume discretisation on a cell-centred radial grid,
marched with backward Euler. Finite volume rather than finite difference for one specific reason:
the inner face of the first cell has zero area, so the symmetry condition holds **identically by
construction** rather than being imposed as an extra equation that can be got wrong.

The scheme reproduces the analytic uniform-generation steady state, centre-minus-surface
$q'''R^2/4k$ and surface-minus-ambient $q'''R/2h$, to within **1×10⁻¹³ relative** at every grid
resolution tested from N = 10 to N = 160 **(a)**. Twenty-four such checks pass, covering global
energy conservation, the parabolic profile shape, the adiabatic limit, the lumped cooling time
constant, axis symmetry and time-step refinement.

Two defects were caught by those checks and would have been invisible otherwise. The convective
boundary source had been added to every control volume rather than only the outermost one,
injecting heat throughout the domain; the signature was an error that grew *linearly with mesh
count*, and the global energy balance was in error by a factor of 2.5. Separately, the half-cell
wall conductance neglected the area variation between the last node and the wall, an O(Δr/4R) bias
sitting in precisely the centre-to-surface quantity under test. Both are now fixed and both are
covered by regression tests.

### 6.4.3 Forward gate: the network must solve the easy direction first

Before attempting the inverse, the network was required to reproduce the finite-volume solution
with every parameter known, scored on the **centre-to-surface difference** rather than the surface
trajectory. The trajectory is easy and carries almost no spatial information; the gradient is the
entire quantity of interest.

The first attempt failed, and instructively. A full-field network with Fourier time features
reached only **39.2 per cent** relative error on the gradient against a 5 per cent target. The
partial differential equation residual stalled at 2.68 while the dimensionless source term has a
mean of 2.3. The residual was the same size as the forcing, so the equation was not being
satisfied at all.

Two causes were measured rather than assumed. First, the drive-cycle source is genuinely
broadband: half of its spectral energy lies above harmonic 233 of the record, corresponding to
periods under 15 s, and 99 per cent below 2.3 s. Representing that with Fourier features would
require roughly 1500 of them. Second, and less expected, the Fourier features were actively
harmful. On a constant-source control problem with a known answer, the gradient error rose
monotonically with the number of features: 0.086 per cent with none, 0.150 per cent with four,
0.493 per cent with sixteen, and 4.88 per cent with sixty-four **(a)**. They could not reach the
bandwidth the problem needed and meanwhile degraded the optimisation landscape for the smooth
solution that remained.

The remedy is exact algebra rather than a tuning change. Writing the field as the sum of the
lumped zero-dimensional solution and a deviation,

$$T(r,t) = T_l(t) + w(r,t),
\qquad \rho c_p \frac{\mathrm{d}T_l}{\mathrm{d}t} = q'''(t) - \frac{2h}{R}\bigl[T_l - T_\infty\bigr],$$

and substituting into the governing equation, the spiky source **cancels identically**, leaving

$$\rho c_p \frac{\partial w}{\partial t} = k\,\nabla^2 w + \frac{2h}{R}\bigl[T_l(t) - T_\infty\bigr],
\qquad w(r,0) = 0 .$$

The physical reading is that the fast, spatially uniform response to a spiky load is something we
can integrate exactly, so there is no reason to ask a network to learn it. What remains is driven
by $T_l - T_\infty$, the output of a first-order filter with a 405 s time constant, which
attenuates 15 s content by a factor of about 169. The deviation $w$ is smooth.

This does not reduce the network's role to bookkeeping. The lumped solution $T_l$ carries **no
radial information whatsoever**, so the entire core-to-surface gradient, the only quantity this
chapter is about, lives in $w$. The split removes the part that can be integrated exactly and
leaves the network exactly the part that cannot.

With that formulation the gate passes: **4.66 per cent** on record 2 and **1.39 per cent** on
record 1, against the 5 per cent criterion, and roughly six times faster than the failed attempt.
Record 2 sits close to the line and that is reported rather than smoothed.

---

## 6.5 The three estimators

All three see the same information: the surface temperature trace, the measured current and
voltage, the ambient trace, and the convection coefficient and conductivity of Table 6.2. None
sees the core.

**Estimator A — classical transient.** The finite-volume model of Section 6.4.2 with a single
free parameter, an effective heat-generation coefficient $R_\mathrm{eff}$ entering through

$$q'''(t) = \frac{I(t)^2\,R_\mathrm{eff}}{V_\mathrm{cell}},
\qquad V_\mathrm{cell} = 3.4219\times10^{-5}\ \mathrm{m^3}\ \textbf{(a)}\ [1].$$

The current-squared factor is not decorative. These cycles contain genuine zero-current rests, and
a source written without it keeps injecting heat through them. Removing that factor as a
deliberate control costs a factor of 2.7 on the surface fit and 1.8 on the core — the same class
of defect that Chapter 5 diagnosed in the earlier inverse network.

The recovered quantity is called an **effective heat-generation coefficient**, never an internal
resistance, because it blends irreversible and reversible contributions that a single temperature
channel cannot separate.

**Estimator B — inverse PINN.** The base-subtracted formulation of Section 6.4.3, with
$R_\mathrm{eff}$ trained jointly with the network weights. Symmetry at the axis is enforced by
construction by taking $s = (r/R)^2$ as the spatial input, under which the radial Laplacian becomes
$4\,\partial_s w + 4s\,\partial_{ss} w$ with no 1/r term at all — the singularity is removed
analytically rather than masked. The initial condition is exact because $w(r,0) = 0$ by
construction, with the measured initial offset carried entirely by $T_l(0)$. Six random seeds were
run for every configuration; L-BFGS used a frozen collocation set, with the closure-call count
reported as a convergence check.

**Estimator C — quasi-steady algebraic relation.** For steady uniform generation in a cylinder the
centre-to-surface difference and the surface rise above ambient are, exactly,

$$T(0) - T(R) = \frac{q''' R^2}{4k},
\qquad T(R) - T_\infty = \frac{q''' R}{2h},$$

whose ratio is $hR/2k = \mathrm{Bi}/2$, independent of the source. Eliminating $q'''$ between them
gives a predictor that never references the heat generation at all:

$$\boxed{\;T_\mathrm{core}(t) = T_\mathrm{surf}(t) + \frac{\mathrm{Bi}}{2}\bigl[T_\mathrm{surf}(t) - T_\infty(t)\bigr]\;}$$

Read plainly: if the radial profile keeps the same shape and only changes amplitude, then measuring
how far the surface sits above ambient tells you how far the core sits above the surface, through
one dimensionless number. This relation appears in the earlier chapters only as a plausibility
check on other estimators. Here it is used as a predictor, and Section 6.6 shows it outperforms
both of them.

---

## 6.6 Results

**Table 6.3** — Predicted core against measured core. Surface-only fits; the core channel entered
only at scoring. All values re-derived from the archived arrays for this chapter.

| Estimator | Record 2 core RMSE | Record 2 max | Record 1 core RMSE | Record 1 max |
|---|---|---|---|---|
| **C — quasi-steady relation** | **0.1409 K** | 0.3742 K | **0.3366 K** | 0.6900 K |
| B — inverse PINN, scalar source | 0.6135 K | 1.9093 K | 0.5024 K | 1.3208 K |
| B — inverse PINN, linear shape | 0.5710 K | 1.4981 K | — | — |
| A — classical transient | 0.8943 K | 2.7087 K | 1.0544 K | — |
| *measured peak gradient, for scale* | *6.5426 K* | | *7.1314 K* | |

The internal temperature is reconstructed from surface data to **0.14 K** on record 2 and 0.34 K
on record 1, against gradients peaking at 6.5 and 7.1 K. Expressed against the quantity being
reconstructed, the best estimator's error is 2.2 per cent of the peak gradient.

The ordering is identical on two independent records, with the roles of the records swapped
between them: the algebraic relation is best, the network second, the transient solver third.

### 6.6.1 What the network's number does and does not mean

The network meets its pre-registered target of better than 1 K, converged on six of six seeds for
every configuration, and shows a seed-to-seed standard deviation of 0.0086 K on the core RMSE with
a 0.49 per cent spread in the recovered coefficient. Those are the favourable facts, and they are
real. Two others belong beside them rather than in a later section, because quoting 0.6135 K
without them overstates what was demonstrated.

**The network does not satisfy the equation it is built on.** Taking its own recovered coefficient
and solving the same physics exactly:

**Table 6.4** — Physics consistency of the network solution, record 2.

| | Surface RMSE | Core RMSE |
|---|---|---|
| Network's own field | 0.1974 K | 0.6135 K |
| Same $R_\mathrm{eff}$, equation enforced exactly | 0.8395 K | 1.5238 K |
| **Difference — the residual violation** | **0.8029 K** | **1.0570 K** |

The network's field departs from a true solution of its own governing equation by 0.80 K on the
surface, which is larger than its entire core error. Its recovered coefficient, used correctly,
gives 1.52 K — worse than the classical estimator's 0.89 K. The network is therefore behaving as a
flexible interpolator with a physics-flavoured penalty, and its good core number is substantially
a consequence of relaxing the constraint rather than obeying it.

A supporting observation makes this concrete. The classical estimator reaches a surface RMSE of
0.4652 K, and because it satisfies the equation exactly, that is the best achievable *within the
physics* for this model class. The network reaches 0.1974 K. Fitting the surface better than the
physics permits is only possible by not solving the physics.

**Its pre-registered success depends on a hyperparameter that truth-free selection gets wrong.**
Sweeping the weight on the data term, with everything else held:

**Table 6.5** — Data-weight sweep, record 2. The residual column is the selection criterion this
project mandates; the core column was not visible during selection.

| $w_\mathrm{data}$ | $R_\mathrm{eff}$ (mΩ) | Surface RMSE | PDE + BC residual | Core RMSE | That $R_\mathrm{eff}$ solved exactly |
|---|---|---|---|---|---|
| 5 | 13.32 | 0.4827 K | **0.1015 (lowest)** | **1.0500 K** | 1.3683 K |
| 20 | 13.62 | 0.3025 K | 0.1136 | 0.7403 K | 1.1695 K |
| **200 (used)** | 13.05 | 0.1963 K | 0.1436 | **0.6332 K** | 1.5650 K |
| 2000 | 11.26 | 0.1374 K | 0.2665 | 0.7721 K | 3.1037 K |

Every trend is monotonic. More weight on the data buys surface accuracy and pays for it in
equation residual, in physics violation, and in bias on the recovered coefficient, which reaches
11.26 mΩ — 22 per cent below both independent anchors. The independent anchors are the classical
fit at 14.30 mΩ and the electrical regression of terminal voltage on current at 14.49 mΩ, which
agree with each other to 1.3 per cent and disagree with the network's 13.05 mΩ by about 9 per cent.

The difficulty is in the residual column. Selecting the weight by the lowest equation residual, the
truth-free criterion this project requires precisely so that selection cannot see the answer, picks
$w_\mathrm{data} = 5$, which yields **1.05 K and would have failed** the pre-registered target. The
reported value used 200, fixed before any core scoring and never revised, which happens to sit near
the core-error minimum. That is a fortunate choice, not a defensible method, and it is recorded as
such. Chapter 5's synthetic sweep had independently found 20 optimal and 200 markedly worse on a
different problem, which sharpens the point rather than softening it.

### 6.6.2 Scoring against trivial baselines

A method is only as impressive as the baseline it beats, so three predictors requiring no
estimation at all were scored on the same data.

**Table 6.6** — Skill against baselines, record 2.

| Predictor | Core RMSE |
|---|---|
| $T_\mathrm{core} = T_\mathrm{surf}$ (no gradient at all) | 5.1855 K |
| $T_\mathrm{core} = T_\mathrm{surf} + \overline{\Delta T}$ (given the measured mean gradient) | 1.2101 K |
| Inverse PINN | 0.6135 K |
| **Quasi-steady relation** | **0.1409 K** |

The second baseline is deliberately unfair — it is handed the mean of the quantity being predicted,
which no core-blind method could know. The network beats it, which is a genuine result. The
algebraic relation beats the network by a factor of 4.4 while using strictly less information than
the second baseline, since Bi comes from a fit on a different record and no core statistic enters
at any point.

---

## 6.7 Why the algebraic relation wins

The mechanism is that **the problem is quasi-steady**, and this is measurable rather than argued.
Over the thermally active part of each record, the ratio of the measured core-to-surface difference
to the measured surface rise above ambient is

$$\frac{T_\mathrm{core} - T_\mathrm{surf}}{T_\mathrm{surf} - T_\infty}
= 0.6161 \pm 0.0126 \ \ (\text{record 2}), \qquad
0.6218 \pm 0.0103 \ \ (\text{record 1}),$$

a variability of 2.0 and 1.7 per cent across a 3541 s cycle carrying 30 A peaks **(a)**. The shape
of the radial profile does not change; only its amplitude does. One number therefore captures the
entire spatial structure, and a transient solver has nothing left to contribute.

Those measured ratios are to be compared with $\mathrm{Bi}/2 = 0.6071$ and $0.5744$ from Table 6.2,
agreeing to 1.5 and 7.6 per cent. Both came from fits on the other record, with no core data
anywhere in the chain.

There is a second, reinforcing reason, and the two should not be conflated. The relation **never
references the volumetric source**. Estimators A and B must both recover $R_\mathrm{eff}$ and then
propagate it through a thermal history, so any error in the heat-generation model accumulates into
the predicted temperature. Estimator C sidesteps the source entirely by reading the surface
temperature at every instant, which is the quantity the source would have been used to predict.

That second reason invites an obvious objection: the algebraic relation is anchored on the measured
surface at every step, while the network re-predicts the surface and therefore inherits its own
surface error. If that were the whole story the comparison would be unfair. It was tested.
Removing the network's surface error, by taking only its predicted *gradient* and adding it to the
measured surface, gives **0.5990 K** against the relation's 0.1409 K. Anchoring is not what
separates them. The same test applied to the classical estimator gives 0.5622 K. The gap lives in
the gradient itself.

**What the transient models are still for.** The relation requires Bi, and Bi is not free. Varying
it by ±15 per cent degrades the core RMSE by a factor of five to ten (Section 6.8). It was obtained
from a **transient** surface-only fit on an independent record — which is to say, the transient
modelling of Chapters 3 to 5 is what makes the algebraic predictor usable. The honest division of
labour is to identify Bi with a transient model and to predict with algebra. What this chapter does
not find is any role for the neural formulation in either half.

---

## 6.8 Sensitivity to the Biot number

The result rests on a parameter that the data cannot determine. Refitting the source coefficient at
each assumed conductivity, against the same surface trace:

**Table 6.7** — Conductivity sensitivity, record 2, classical estimator.

| k (W m⁻¹K⁻¹) | Bi | Surface RMSE | Core RMSE |
|---|---|---|---|
| 0.30 | 1.594 | 0.5433 K | 1.8209 K |
| 0.35 | 1.366 | 0.4961 K | 1.0998 K |
| **0.3937 (used)** | 1.214 | 0.4652 K | **0.8943 K** |
| 0.404 [1] | 1.183 | 0.4590 K | 0.9027 K |
| 0.45 | 1.062 | 0.4355 K | 1.0941 K |
| 0.55 | 0.869 | **0.4007 K (best)** | 1.6983 K |

This table is the most important diagnostic in the chapter. The **surface fit improves
monotonically** as conductivity rises, while the **core error passes through a minimum near
k ≈ 0.39 and roughly doubles either side**. At k = 0.55 the surface fit is the best in the table
and the core prediction is nearly twice as bad. Judged on what is observable, the worst entry looks
like the best.

The Cramér–Rao analysis of Chapter 4 shows why: conductivity and convection coefficient carry a
correlation of 0.998 in this configuration, so the surface trace constrains their combination and
not their separation. The core prediction is essentially a function of Bi, and **Bi is an input to
this work, not an output of it**. The value used here, obtained from a surface-only fit on the
other record, lands within about 1 per cent of the core-optimal value. That is either good transfer
between records or luck, and two records cannot distinguish the two possibilities **(c)**.

For completeness, switching on the reversible term with the published lithium iron phosphate
entropy coefficient of −0.5 mV K⁻¹ at 50 per cent state of charge **(a)** [2] moves the core RMSE
from 0.8895 K to 0.8927 K, a 0.4 per cent change. The instantaneous reversible heat is not small:
its mean absolute value is 57 to 63 per cent of the mean ohmic heat. But it alternates with charge
and discharge, so its signed mean is only −1.4 per cent, and its effect on the thermal prediction
is immaterial on this duty cycle **(c)**.

---

## 6.9 Where the error lives

The reported RMSE is an average over conditions that are not equivalent, and the average conceals
the structure.

**Table 6.8** — Temporal decomposition of the network's core error, record 2.

| Window | Core RMSE | Max | Bias | Share of squared error |
|---|---|---|---|---|
| First 10 min | 1.1146 K | 1.9093 K | −0.9973 K | **56 %** |
| After first 10 min | 0.4470 K | 1.2363 K | −0.2273 K | 44 % |
| Whole record | 0.6135 K | 1.9093 K | −0.3577 K | 100 % |

The first ten minutes occupy 17 per cent of the record and carry 56 per cent of the squared error.
The interpretation is physically clean and it also bounds the chapter's central claim: start-up is
the one regime in which quasi-steadiness is false. The cell begins nearly isothermal, so the radial
profile has not yet developed, the ratio of Section 6.7 has not yet reached its plateau, and every
estimator that assumes a developed profile is wrong in the same direction. Excluding start-up, the
network reaches 0.45 K.

This is also where the two records differ. The correlation between predicted and measured gradient
is 0.933 on record 2 and 0.988 on record 1, and record 1's longer duration gives the profile more
time to develop relative to the record length.

---

## 6.10 Summary

Fitting only a surface temperature trace, the internal temperature of an A123 26650 cell was
predicted and scored against a core thermocouple that no estimator saw. The best error is
**0.1409 K on record 2 and 0.3366 K on record 1**, against measured gradients peaking at 6.54 and
7.13 K. This is the first internal-temperature validation in this project, and it closes the gap
identified at the end of Chapter 5.

The estimator that achieves it is a single algebraic relation, not the inverse network. The network
meets its pre-registered target at 0.6135 K, converges reliably, and beats the classical transient
solver on both records — but it does not satisfy its own governing equation at the operating point
where it performs best, and its success depends on a weighting that the project's own truth-free
selection rule would have chosen badly.

The mechanism is quasi-steadiness: the radial profile shape is constant to about 2 per cent across
both records, so one dimensionless number captures the spatial structure.

That raises the question this chapter cannot answer from two drive cycles alone — whether
quasi-steadiness is a property of this cell or of this *duty cycle*. Chapter 7 settles it on an
independent multi-rate dataset covering the same cell family across eight discharge rates, and
finds the answer is a timescale criterion: the quasi-steady error follows

$$\varepsilon \approx 7.5\,\left(\frac{t_\mathrm{forcing}}{\tau_\mathrm{diff}}\right)^{-1.18}
\ \ \text{per cent}, \qquad r^2 = 0.996,$$

crossing 5 per cent at a timescale ratio of 1.4 and 25 per cent at 0.37. The drive cycles analysed
here sit at a ratio of 3.4, comfortably inside the quasi-steady regime — which is why the algebraic
relation wins on this data, and where it would stop winning.

---

## References cited in this chapter

[1] R. R. Richardson and D. A. Howey, "Sensorless battery internal temperature estimation using a
Kalman filter with impedance measurement," *IEEE Transactions on Sustainable Energy*, vol. 6,
no. 4, pp. 1190–1199, 2015. *(Verified: repository `README.md` of the accompanying code archive,
which also supplies the parameter values and channel assignments used here.)*

[2] C. Forgez, D. V. Do, G. Friedrich, M. Morcrette and C. Delacourt, "Thermal modeling of a
cylindrical LiFePO₄/graphite lithium-ion battery," *Journal of Power Sources*, vol. 195,
pp. 2961–2968, 2010. *(Verified: cited as the source of the entropy coefficient in the University
of Michigan electro-thermal model release, `00 read me first.txt`, item 10.)*

`[VERIFY: page numbers and DOI for both references against the publisher record before submission —
metadata above is taken from the code archives, not from the articles themselves.]`

`[VERIFY: A123 ANR26650M1 datasheet for nominal capacity, mass and DC resistance — currently taken
from the Richardson archive and from typical values, not from the manufacturer document.]`
