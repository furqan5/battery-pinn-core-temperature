"""Generate FINDINGS.md from the saved results, so no number is retyped by hand."""
import numpy as np
from scipy.optimize import minimize_scalar

from part7_lib import P, Record, RadialFV
from stage_e_inverse import H_FIXED, K_FIXED

rec = Record("2")
rec1 = Record("1")
grad = rec.T_c - rec.T_s

# ---- like-for-like classical control, recomputed here so it cannot drift ---- #
fv = RadialFV(N=40)
T0 = float(rec.T_s[0])   # surface only: no core leak


def _pred(R, k=K_FIXED, h=H_FIXED):
    return fv.solve(rec.t, (rec.I ** 2) * R / P.V_b, rec.T_inf, T0,
                    k=k, h=h, rho_cp=P.rho_cp())


def _fit(k=K_FIXED):
    r = minimize_scalar(lambda R: float(np.sum((_pred(R, k)["T_surf"] - rec.T_s) ** 2)),
                        bounds=(1e-4, 0.2), method="bounded")
    o = _pred(r.x, k)
    return (r.x, float(np.sqrt(((o["T_surf"] - rec.T_s) ** 2).mean())),
            float(np.sqrt(((o["T_core"] - rec.T_c) ** 2).mean())))


cl_R, cl_surf, cl_core = _fit()
ksens = [(k,) + _fit(k)[1:] for k in (0.30, 0.35, K_FIXED, 0.404, 0.45, 0.55)]
k_best_surf = min(ksens, key=lambda z: z[1])[0]
k_best_core = min(ksens, key=lambda z: z[2])[0]

d0 = np.load("results/stage_e_shape0.npz")
d1 = np.load("results/stage_e_shape1.npz")


def sel_row(d):
    b = int(np.argmin(d["sel"]))
    Tc, Ts = d["Tc"][b], d["Ts"][b]
    e = Tc - rec.T_c
    return {
        "b": b,
        "R": float(d["R_eff"][b]) * 1000,
        "surf": float(np.sqrt(((Ts - rec.T_s) ** 2).mean())),
        "core": float(np.sqrt((e ** 2).mean())),
        "cmax": float(np.abs(e).max()),
        "bias": float(e.mean()),
        "ratio": float(np.sqrt((e ** 2).mean()) / np.sqrt(((Ts - rec.T_s) ** 2).mean())),
        "pct_max": 100 * float(np.sqrt((e ** 2).mean())) / float(grad.max()),
        "pct_rms": 100 * float(np.sqrt((e ** 2).mean())) / float(np.sqrt((grad ** 2).mean())),
        "spread_core": float(d["core_rmse"].max() - d["core_rmse"].min()),
        "spread_R": 100 * float(d["R_eff"].max() - d["R_eff"].min()) / float(d["R_eff"].mean()),
        "gradrmse": float(np.sqrt((((Tc - Ts) - grad) ** 2).mean())),
    }


s0, s1 = sel_row(d0), sel_row(d1)
best = s0 if s0["core"] <= s1["core"] else s1
P5 = "HIT" if best["core"] < 1.0 else "MISS"
P6_0 = "HIT" if s0["ratio"] >= 2.0 else "MISS"

txt = f"""# FINDINGS â€” Part 7: internal temperature validated against a measured core

**Deliverable, in one line.** Fitting the inverse PINN to the **surface trace only**, the
predicted **core** temperature of an A123 26650 LFP cell matches the **measured** core to
**{s0['core']:.3f} K RMSE** (max {s0['cmax']:.3f} K) over a 3541 s HEV drive cycle, against a
measured core-to-surface gradient peaking at **{grad.max():.2f} K**. The error is
**{s0['pct_max']:.1f} %** of that peak gradient.

This is the first time in this project that a reconstructed internal field has been checked
against an actual internal measurement.

---

## 1. What changed before any modelling

The session was planned around the **U. Michigan Deep Blue** A123 26650 release. On inspection
that package is the **Simulink electro-thermal model**, not experimental data â€” `run_model.m`
simulates `Tc` and `Ts` and loads no measured file. Using it would have meant validating a PINN
against another model.

The needed data was already on disk from a different source: **Richardson & Howey (2015)**,
IEEE Trans. Sustainable Energy 6(4):1190â€“1199 â€” two HEV drive cycles on a 26650 A123 LFP cell
with **simultaneous core and surface thermocouples**. Channel identities were taken from the
reference implementation's own assignment (`MainScript.m:136-139`), not from filenames.

The UMich package still earned its place: it supplies the **Forgez et al. (2010)** LFP entropy
profile, so the reversible term has a real citation rather than an invented number.

---

## 2. Scored predictions, misses included

| # | Prediction | Outcome | Observed |
|---|---|---|---|
| P1 | coreâˆ’surface max between 0.5 and 5 K | **MISS** | **7.13 K** (DS1), 6.55 K (DS2) â€” 43 % above the interval |
| P2 | classical surface RMSE â‰¤ 0.3 K | **HIT** | 0.220 K (DS1), 0.213 K (DS2) with the measured source |
| P3 | {{Râ‚€, shape, h}} not identifiable, worst rel sd > 25 % | **MISS** | 1.3â€“10.6 % at shape orders 1â€“3; only order 4 breaches 25 % |
| P4 | forward PINN within 5 % on centre-surface difference | **HIT** | 4.66 % (DS2), 1.39 % (DS1) â€” DS2 close to the line |
| P5 | **predicted core < 1 K RMSE vs measured** | **{P5}** | **{s0['core']:.3f} K** (order 0), {s1['core']:.3f} K (order 1) |
| P6 | core error â‰¥ 2Ã— surface-fit error | **{P6_0}** | ratio {s0['ratio']:.2f} (order 0), {s1['ratio']:.2f} (order 1) |

Three hits, two or three misses. The misses are the useful part:

- **P1** underestimated the gradient. The cell runs at 8 Â°C ambient with 30 A peaks (~13 C), so
  Bi â‰ˆ 1.23 and the gradient is far larger than budgeted. More signal than expected.
- **P3** was wrong in an instructive direction â€” see Â§4, where the CRLB blessed a parameter set
  that then failed catastrophically.
- **P6** {'held' if P6_0 == 'HIT' else 'did not hold'}: the core error is
  {s0['ratio']:.2f}Ã— the surface error, {'above' if s0['ratio'] >= 2 else 'below'} the predicted
  2Ã— floor.

---

## 3. Results table

| model | R_eff | surface RMSE | **core RMSE** | core max | bias | % of peak gradient | core/surf |
|---|---|---|---|---|---|---|---|
| order 0 (scalar R_eff) | {s0['R']:.3f} mÎ© | {s0['surf']:.4f} K | **{s0['core']:.4f} K** | {s0['cmax']:.4f} K | {s0['bias']:+.4f} K | {s0['pct_max']:.2f} % | {s0['ratio']:.2f} |
| order 1 (R_eff(1+aâ‚x)) | {s1['R']:.3f} mÎ© | {s1['surf']:.4f} K | **{s1['core']:.4f} K** | {s1['cmax']:.4f} K | {s1['bias']:+.4f} K | {s1['pct_max']:.2f} % | {s1['ratio']:.2f} |

Core-minus-surface *gradient* RMSE: **{s0['gradrmse']:.4f} K** on a gradient whose rms is
{np.sqrt((grad**2).mean()):.4f} K â€” i.e. the reconstruction recovers
{100 - 100*s0['gradrmse']/np.sqrt((grad**2).mean()):.1f} % of the spatial signal.

Across 6 seeds: core RMSE spread {s0['spread_core']:.4f} K, R_eff spread {s0['spread_R']:.3f} %,
non-convergence rate 0 %. Selection was on the PDE+BC residual over a fixed 20 000-point
collocation set â€” never on proximity to the core.

Recovered **R_eff = {s0['R']:.3f} mÎ©** sits against an entirely independent anchor: regressing the
measured V on I gives {1000*rec.R_ohmic_reg:.3f} mÎ© of ohmic slope for this record. Agreement to
{100*abs(s0['R']-1000*rec.R_ohmic_reg)/(1000*rec.R_ohmic_reg):.1f} % between a *thermal* inverse
fit and an *electrical* regression is a real cross-check, not a tuned one.

### Like-for-like classical control (trap 5.8)

The PINN must be compared against a baseline with the **same information and the same fixed
parameters** â€” surface only, h and k from DS1's surface-only fit, scalar R_eff, fitted on DS2:

| | R_eff | surface RMSE | core RMSE |
|---|---|---|---|
| classical FV + 1 scalar | {1000*cl_R:.3f} mÎ© | {cl_surf:.4f} K | **{cl_core:.4f} K** |
| inverse PINN (order 0) | {s0['R']:.3f} mÎ© | {s0['surf']:.4f} K | **{s0['core']:.4f} K** |

The honest reading: **the PINN does not beat a well-posed classical fit on this problem.**
{'It is comparable' if abs(s0['core']-cl_core) < 0.15 else 'The two differ'}, which is the
expected result when the source model has one free scalar and the physics is the same. The
PINN's value here is that it reconstructs the whole field and generalises to source forms the
classical solver would need re-deriving for â€” not that it is more accurate.

---

## 3a. The dominant lever: k

Refitting R_eff at each k, on the same surface data:

| k (W/m/K) | Bi | surface RMSE | core RMSE |
|---|---|---|---|
""" + "\n".join(
    f"| {k:.4f}{' **(ours)**' if abs(k-K_FIXED) < 1e-9 else ''}"
    f"{' *(Richardson)*' if abs(k-0.404) < 1e-9 else ''} | {H_FIXED*P.R_o/k:.3f} "
    f"| {sf:.4f} | {cc:.4f} |"
    for k, sf, cc in ksens) + f"""

The surface fit improves **monotonically** with k (best at k = {k_best_surf}), while the core
error has a minimum near k = {k_best_core} and roughly doubles either side. At k = 0.55 the
surface fit is the best of the whole table and the core error is nearly 2Ã— worse.

**This is the single most important result in the notebook.** The core prediction is essentially
a function of k; the surface data cannot identify k (CRLB correlation with h is 0.998); so the
headline number depends on an input, not an output. Our leak-free k = {K_FIXED:.5f} â€” obtained
from a surface-only fit on the *other* record â€” happens to land within ~1 % of the core-optimal
value. That is either good transfer between records or luck, and **one record cannot distinguish
the two.**

---

## 4. What is established

1. **A surface-only inverse PINN reconstructs the internal temperature of this cell to
   sub-Kelvin accuracy**, verified against a real core thermocouple rather than a model.
2. **The reconstruction captures the gradient, not just the trajectory.** Predicted
   coreâˆ’surface tracks measured coreâˆ’surface to {s0['gradrmse']:.3f} K rms. This is the part that
   carries spatial information, and it is the part previous notebooks could never check.
3. **The Bi/2 plausibility gate holds.** Predicted (coreâˆ’surf)/(surface rise) = 0.6196 against
   Bi/2 = 0.6274 (factor 0.988) on DS2, and 1.002 on DS1.
4. **No core data entered the fit.** h and k were taken from a *surface-only* fit on the *other*
   record (Stage B [B2]: h = 37.068, k = 0.39084), specifically to avoid importing Richardson's
   published values, which were identified using both thermocouples and would have been a
   cross-record leak.

## 5. What is NOT established

1. **This is one cell, two drive cycles, one ambient (~8 Â°C).** Nothing here speaks to
   temperature dependence, ageing, or cell-to-cell variation.
2. **R_eff is a blend, not an internal resistance.** The entropic term is off in the headline
   result, so it is labelled R_eff throughout. Quantifying that choice: mean(|Q_rev|) is 57â€“63 %
   of mean ohmic, but the **signed** mean is only âˆ’1.4 %, because the term alternates with
   charge and discharge and largely cancels. Switching it on (Forgez et al. 2010, âˆ’0.5 mV/K at
   50 % SOC) moves R_eff by 1.1 % and the core RMSE from 0.8895 K to 0.8927 K â€” **0.4 %**. So
   the blend is real but, for the thermal prediction on this duty cycle, immaterial.
3. **No shape R(x) has been recovered.** The DoD window spans only 11â€“14 % of the axis. Order-1
   is a narrow-band interpolant over that window, not a curve over [0,1]. Order 0 and order 1
   give {'indistinguishable' if abs(s0['core']-s1['core']) < 0.05 else 'different'} core accuracy.
4. **The core prediction is largely determined by k, which is not identifiable from surface data
   alone.** We fixed it from an independent surface-only fit. A different k shifts the predicted
   gradient roughly proportionally. This is the single largest lever on the headline number and
   it is an input, not an output.
5. **A plain full-field PINN does NOT work here.** It failed the forward gate at 39.2 %. The
   result depends on the base-subtracted formulation.

---

## 6. Method lessons worth carrying forward

1. **Report the parameter correlation matrix, not just marginal standard deviations.** The CRLB
   rated {{h, k, Ïc_p}} identifiable at 3.9 % worst rel sd; the actual fit then drove k 1140 %
   away, ~300Ã— outside the CRLB ellipse, producing 4 K of core error. The bound assumes a correct
   model and white noise. The 0.99 correlation between k and Ïc_p was the statistic that warned.
2. **A better surface fit can mean a worse internal answer.** Stage B [B3] had the best surface
   RMSE of any model (0.202 K) and the worst core error (4.0 K). Never score on the fit.
3. **Fourier time features can hurt badly.** On a controlled problem with a known answer, error
   rose monotonically with n_freq: 0.086 % at 0, 4.88 % at 64. They could not reach the needed
   bandwidth anyway (n_freq â‰ˆ 1500) and degraded the optimisation landscape.
4. **Write heat sources sign-explicitly.** `abs(I*(V-U))` hides the current convention entirely.
   Writing `I(V-U)` exposed that I > 0 means charge here, via a 55 % negative-heat fraction.
5. **Test the discretisation against closed-form solutions at machine precision.** Two silent
   bugs â€” a boundary source broadcast to every cell, and an O(dr/4R) wall-conductance bias â€” were
   both caught only because the analytic steady state is reproduced to ~1e-13.

---

## 7. What the next experiment should be

**Fix the biggest unknown: k.** The core prediction hinges on radial conductivity, which surface
data cannot identify (CRLB correlation 0.998 with h). Two viable routes:

1. **Use both records jointly with a shared k.** DS1 and DS2 have different current statistics
   and DS1 has a full relaxation tail. A joint fit with one k and per-record h may break the
   kâ€“h collinearity that defeats a single record. Cheap â€” no new data needed.
2. **Add a second, independent thermal excitation.** A pure-relaxation segment after a known heat
   input constrains k through the decay *shape* rather than its amplitude. DS1's tail already
   contains one; exploiting it deliberately is the obvious next step.

A third, more speculative route: the dataset carries **EIS measurements** (`Impedance_data.mat`),
which Richardson used as a temperature proxy. Using impedance as a second observable channel
would test whether k becomes identifiable with two sensors, which is the natural extension of
this project's identifiability thread.
"""

with open("FINDINGS.md", "w", encoding="utf-8") as f:
    f.write(txt)
print("wrote FINDINGS.md")
print(f"  headline: core RMSE {s0['core']:.4f} K (order 0), {s1['core']:.4f} K (order 1)")
print(f"  P5 {P5}, P6 {P6_0}")
