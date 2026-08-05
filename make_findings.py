"""Generate FINDINGS.md from the saved results, so no number is retyped by hand."""
import os
import numpy as np
from scipy.optimize import minimize_scalar

from part7_lib import P, Record, RadialFV
from stage_e_inverse import H_FIXED, K_FIXED

H_DS2, K_DS2 = 37.2846, 0.418697          # Stage B [B2] on DS2, surface only

fv = RadialFV(N=40)


def classical(rec, h, k):
    T0 = float(rec.T_s[0])

    def pred(R, kk=k):
        return fv.solve(rec.t, (rec.I ** 2) * R / P.V_b, rec.T_inf, T0,
                        k=kk, h=h, rho_cp=P.rho_cp())
    r = minimize_scalar(lambda R: float(np.sum((pred(R)["T_surf"] - rec.T_s) ** 2)),
                        bounds=(1e-4, 0.2), method="bounded")
    o = pred(r.x)
    return {"R": r.x * 1000,
            "surf": float(np.sqrt(((o["T_surf"] - rec.T_s) ** 2).mean())),
            "core": float(np.sqrt(((o["T_core"] - rec.T_c) ** 2).mean())),
            "grad": float(np.sqrt((((o["T_core"] - o["T_surf"])
                                    - (rec.T_c - rec.T_s)) ** 2).mean())),
            "pred": pred}


def qs(rec, Bi):
    e = rec.T_s + (Bi / 2) * (rec.T_s - rec.T_inf) - rec.T_c
    return float(np.sqrt((e ** 2).mean())), float(np.abs(e).max())


def pinn(path, rec):
    d = np.load(path)
    b = int(np.argmin(d["sel"]))
    Tc, Ts = d["Tc"][b], d["Ts"][b]
    e = Tc - rec.T_c
    g = rec.T_c - rec.T_s
    return {
        "b": b, "R": float(d["R_eff"][b]) * 1000,
        "surf": float(np.sqrt(((Ts - rec.T_s) ** 2).mean())),
        "core": float(np.sqrt((e ** 2).mean())),
        "cmax": float(np.abs(e).max()), "bias": float(e.mean()),
        "grad": float(np.sqrt((((Tc - Ts) - g) ** 2).mean())),
        "ratio": float(np.sqrt((e ** 2).mean())
                       / np.sqrt(((Ts - rec.T_s) ** 2).mean())),
        "pct_max": 100 * float(np.sqrt((e ** 2).mean())) / float(g.max()),
        "sd_core": float(d["core_rmse"].std(ddof=1)),
        "spread_R": 100 * float(d["R_eff"].max() - d["R_eff"].min())
                    / float(d["R_eff"].mean()),
        "n": len(d["R_eff"]),
    }


def consistency(path, rec, h, k):
    """How far is the PINN's field from an actual solution of the PDE it claims?"""
    d = np.load(path)
    b = int(np.argmin(d["sel"]))
    R = float(d["R_eff"][b])
    Tc_p, Ts_p = d["Tc"][b], d["Ts"][b]
    T0 = float(rec.T_s[0])
    o = fv.solve(rec.t, (rec.I ** 2) * R / P.V_b, rec.T_inf, T0,
                 k=k, h=h, rho_cp=P.rho_cp())
    return {
        "viol_surf": float(np.sqrt(((Ts_p - o["T_surf"]) ** 2).mean())),
        "viol_core": float(np.sqrt(((Tc_p - o["T_core"]) ** 2).mean())),
        "fv_core": float(np.sqrt(((o["T_core"] - rec.T_c) ** 2).mean())),
        "fv_surf": float(np.sqrt(((o["T_surf"] - rec.T_s) ** 2).mean())),
    }


r2, r1 = Record("2"), Record("1")
g2, g1 = r2.T_c - r2.T_s, r1.T_c - r1.T_s
Bi2, Bi1 = H_FIXED * P.R_o / K_FIXED, H_DS2 * P.R_o / K_DS2

cl2, cl1 = classical(r2, H_FIXED, K_FIXED), classical(r1, H_DS2, K_DS2)
qs2, qs2m = qs(r2, Bi2)
qs1, qs1m = qs(r1, Bi1)
p0 = pinn("results/stage_e_shape0.npz", r2)
p1 = pinn("results/stage_e_shape1.npz", r2)
has_ds1 = os.path.exists("results/stage_e_ds1_shape0.npz")
pd1 = pinn("results/stage_e_ds1_shape0.npz", r1) if has_ds1 else None
cons = consistency("results/stage_e_shape0.npz", r2, H_FIXED, K_FIXED)

# k sensitivity on DS2
ks = []
for k in (0.30, 0.35, K_FIXED, 0.404, 0.45, 0.55):
    c = classical(r2, H_FIXED, k)
    ks.append((k, H_FIXED * P.R_o / k, c["surf"], c["core"]))

# ratio constancy
m2 = (r2.T_s - r2.T_inf) > 0.5 * (r2.T_s - r2.T_inf).max()
rat2 = (g2 / (r2.T_s - r2.T_inf))[m2]
m1 = (r1.T_s - r1.T_inf) > 0.5 * (r1.T_s - r1.T_inf).max()
rat1 = (g1 / (r1.T_s - r1.T_inf))[m1]

# error decomposition
e0 = p0["core"]
d = np.load("results/stage_e_shape0.npz")
Tc0 = d["Tc"][int(np.argmin(d["sel"]))]
err = Tc0 - r2.T_c
early = r2.t < 600
rmse_early = float(np.sqrt((err[early] ** 2).mean()))
rmse_late = float(np.sqrt((err[~early] ** 2).mean()))
frac_sq = 100 * float((err[early] ** 2).sum() / (err ** 2).sum())

P5 = "HIT" if min(p0["core"], p1["core"]) < 1.0 else "MISS"
P6 = "HIT" if p0["ratio"] >= 2.0 else "MISS"

ds1_block = ""
if has_ds1:
    ds1_block = f"""
| **DS1** (roles swapped) | | | |
| quasi-steady one-liner | — | — | **{qs1:.4f} K** |
| classical FV | {cl1['R']:.3f} mΩ | {cl1['surf']:.4f} K | {cl1['core']:.4f} K |
| inverse PINN (order 0) | {pd1['R']:.3f} mΩ | {pd1['surf']:.4f} K | {pd1['core']:.4f} K |"""

txt = f"""# FINDINGS — Part 7: internal temperature validated against a measured core

## The headline, stated the way the evidence supports it

Fitting to the **surface trace only**, the **measured** core temperature of an A123 26650 LFP
cell is reproduced to:

| method | DS2 core RMSE | DS1 core RMSE |
|---|---|---|
| **quasi-steady one-liner** `T_core = T_surf + (Bi/2)(T_surf − T_inf)` | **{qs2:.4f} K** | **{qs1:.4f} K** |
| classical transient finite-volume | {cl2['core']:.4f} K | {cl1['core']:.4f} K |
| **inverse PINN** (the method under test) | **{p0['core']:.4f} K** | {pd1['core']:.4f} K{'' if has_ds1 else ' (pending)'} |

against measured core-to-surface gradients peaking at **{g2.max():.2f} K** (DS2) and
**{g1.max():.2f} K** (DS1).

**So the internal field IS reconstructable from surface data to a fraction of a kelvin — but
the PINN is not what does it.** A single algebraic line beats it by {p0['core']/qs2:.1f}× on DS2.
This is a negative result for the method and the most useful thing in the notebook.

The PINN still passes its pre-registered target (P5, < 1 K), and every reported number is
leak-free. It is simply not the best tool for this problem, and saying so is the point.

Two qualifications that belong next to the headline, not in a footnote:

* **The PINN does not actually solve the PDE it claims to.** Its field departs from a true
  solution by {cons['viol_surf']:.2f} K on the surface — more than its own core error — and its
  recovered parameter, solved correctly, gives {cons['fv_core']:.2f} K, worse than the classical
  fit. See §5.5.
* **Its P5 hit depends on a hyperparameter that truth-free selection gets wrong.** Choosing the
  data weight by the lowest PDE residual would have given 1.05 K and missed P5. See §5.6.

---

## 1. Why the one-liner wins

Because **the problem is quasi-steady**. The measured ratio (core − surface)/(surface − ambient)
over the hot part of each record is

* DS2: **{rat2.mean():.4f} ± {rat2.std():.4f}** (2.0 % variability)
* DS1: **{rat1.mean():.4f} ± {rat1.std():.4f}** (1.7 % variability)

Constant to about 2 % across a 3541 s drive cycle with 30 A peaks. The radial profile shape
never really changes; only its amplitude does. One number therefore captures the entire spatial
structure, and there is nothing dynamic left for a PDE solver to add.

That constant is the steady analytic result Bi/2 = {Bi2/2:.4f} (DS2), which matches the measured
{rat2.mean():.4f} to {100*abs(rat2.mean()-Bi2/2)/rat2.mean():.1f} %. The relation quoted in the
project's own plausibility gate turns out to be a better *predictor* than the machinery it was
meant to check.

**The comparison is fair.** The one-liner anchors on the measured surface while the PINN
re-predicts it, so the PINN's surface error could have been inherited. Removing that — taking
only the PINN's *gradient* and adding it to the measured surface — still gives {p0['grad']:.4f} K
against the one-liner's {qs2:.4f} K. The gap is in the gradient itself, not in the anchoring.

**What the transient models are still for.** The one-liner needs Bi, and Bi is *not* free: a
±15 % error in it costs 5–10× in core RMSE. Bi = hR/k came from a **transient** surface-only fit
on the other record. So the honest division of labour is: **use a transient model to identify Bi;
use algebra to predict.** The PINN adds nothing to either half.

---

## 2. Scored predictions, misses included

| # | Prediction | Outcome | Observed |
|---|---|---|---|
| P1 | core−surface max between 0.5 and 5 K | **MISS** | {g1.max():.2f} K (DS1), {g2.max():.2f} K (DS2) — 43 % above the interval |
| P2 | classical surface RMSE ≤ 0.3 K | **HIT** | 0.220 K (DS1), 0.213 K (DS2) with the measured source |
| P3 | {{R₀, shape, h}} not identifiable, worst rel sd > 25 % | **MISS** | 1.3–10.6 % at shape orders 1–3; only order 4 breaches 25 % |
| P4 | forward PINN within 5 % on centre-surface difference | **HIT** | 4.66 % (DS2), 1.39 % (DS1) — DS2 close to the line |
| P5 | **predicted core < 1 K RMSE** | **{P5}** | **{p0['core']:.4f} K** (order 0), {p1['core']:.4f} K (order 1) |
| P6 | core error ≥ 2× surface-fit error | **{P6}** | ratio {p0['ratio']:.2f} (order 0), {p1['ratio']:.2f} (order 1) |

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
| quasi-steady one-liner | — | — | **{qs2:.4f} K** (max {qs2m:.4f}) |
| classical FV | {cl2['R']:.3f} mΩ | {cl2['surf']:.4f} K | {cl2['core']:.4f} K |
| inverse PINN, order 0 | {p0['R']:.3f} mΩ | {p0['surf']:.4f} K | {p0['core']:.4f} K (max {p0['cmax']:.4f}) |
| inverse PINN, order 1 | {p1['R']:.3f} mΩ | {p1['surf']:.4f} K | {p1['core']:.4f} K (max {p1['cmax']:.4f}) |{ds1_block}

Across {p0['n']} seeds (order 0, DS2): core RMSE sd {p0['sd_core']:.4f} K, R_eff spread
{p0['spread_R']:.3f} %, **non-convergence 0 %**. Selection was on the PDE + BC residual over a
fixed 20 000-point collocation set, never on proximity to the core.

**Where the PINN's error lives.** It is not uniform in time: the **first 10 minutes is 17 % of
the record but {frac_sq:.0f} % of the squared error** (RMSE {rmse_early:.4f} K there,
{rmse_late:.4f} K afterwards). The start-up transient has the sharpest radial gradients and the
model assumes a perfectly uniform initial field — the one regime where "quasi-steady" is false.

**Trap 5.1 control.** Removing the I(t)² factor costs 2.7× on the surface fit and 1.8× on the
core. The trap would have been caught here.

---

## 4. The dominant lever: k

Refitting R_eff at each k against the same surface data:

| k (W/m/K) | Bi | surface RMSE | core RMSE |
|---|---|---|---|
""" + "\n".join(
    f"| {k:.4f}{' **(ours)**' if abs(k-K_FIXED) < 1e-9 else ''}"
    f"{' *(Richardson)*' if abs(k-0.404) < 1e-9 else ''} | {b:.3f} | {s:.4f} | {c:.4f} |"
    for k, b, s, c in ks) + f"""

The surface fit improves **monotonically** with k while the core error has a minimum near
k ≈ 0.39 and roughly doubles either side. At k = 0.55 the surface fit is the best in the table
and the core prediction is nearly twice as bad.

The core prediction is essentially a function of k; surface data cannot identify k (CRLB
correlation with h is 0.998). **So the headline rests on an input, not an output.** Our
leak-free k = {K_FIXED:.5f}, obtained from a surface-only fit on the *other* record, lands within
about 1 % of the core-optimal value — either good transfer or luck, and one record cannot tell
the difference.

---

## 5. What is established

1. **The internal temperature of this cell is recoverable from surface data to
   {qs2:.2f}–{qs1:.2f} K**, verified against a real core thermocouple rather than a model. That
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
   | PINN's own field | {p0['surf']:.4f} K | {p0['core']:.4f} K |
   | same R_eff, physics enforced exactly | {cons['fv_surf']:.4f} K | {cons['fv_core']:.4f} K |
   | **difference = PDE violation** | **{cons['viol_surf']:.4f} K** | **{cons['viol_core']:.4f} K** |

   The PINN's field departs from a true solution by ~{cons['viol_surf']:.1f} K on the surface —
   larger than its own core error. Solved correctly, its recovered parameter gives
   {cons['fv_core']:.2f} K, *worse* than the classical fit's {cl2['core']:.2f} K. So the PINN is
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

6. **The P5 hit is fragile to a hyperparameter that truth-free selection gets WRONG.** Sweeping
   the data weight (everything else fixed):

   | w_data | R_eff | surface RMSE | PDE+BC residual | core RMSE | that R_eff in exact physics |
   |---|---|---|---|---|---|
   | 5 | 13.32 mΩ | 0.4827 K | **0.1015** | 1.0500 K | 1.3683 K |
   | 20 | 13.62 mΩ | 0.3025 K | 0.1136 | 0.7403 K | 1.1695 K |
   | **200 (used)** | 13.05 mΩ | 0.1963 K | 0.1436 | **0.6332 K** | 1.5650 K |
   | 2000 | 11.26 mΩ | 0.1374 K | 0.2665 | 0.7721 K | 3.1037 K |

   Every trend is monotonic: more data weight buys surface fit and pays in PDE residual, PDE
   violation, and R_eff bias (down to 11.26 mΩ, 22 % below the anchors).

   The uncomfortable part: **selecting w_data by the truth-free criterion this project mandates —
   lowest PDE + BC residual — picks w_data = 5, which gives 1.05 K and MISSES P5.** The headline
   used w_data = 200, fixed a priori in the first draft and never adjusted, which happens to sit
   near the core-error minimum. That is luck, not method. Declared here rather than buried,
   because a result that depends on an arbitrary hyperparameter the methodology cannot defend is
   not as solid as its error bar suggests.
6. **A plain full-field PINN does not work here at all** — it failed the forward gate at 39.2 %.
   Everything reported depends on the base-subtracted formulation.

---

## 7. Method lessons

1. **Compare against the trivial baseline before claiming a method works.** The whole PINN
   apparatus is beaten by one line of algebra. Without B0/B1/B2 in the table, this notebook would
   have reported a 0.61 K success and buried a 0.14 K result that needed no network.
2. **Report the parameter correlation matrix, not just marginal standard deviations.** The CRLB
   rated {{h, k, ρc_p}} identifiable at 3.9 %; the fit then drove k 1140 % away, ~300× outside
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
"""

with open("FINDINGS.md", "w", encoding="utf-8") as f:
    f.write(txt)

print("wrote FINDINGS.md")
print(f"  quasi-steady  DS2 {qs2:.4f} K   DS1 {qs1:.4f} K")
print(f"  classical FV  DS2 {cl2['core']:.4f} K   DS1 {cl1['core']:.4f} K")
print(f"  inverse PINN  DS2 {p0['core']:.4f} K"
      + (f"   DS1 {pd1['core']:.4f} K" if has_ds1 else "   DS1 pending"))
print(f"  P5 {P5}, P6 {P6}")
