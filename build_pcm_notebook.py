"""Assemble PCM_Identifiability.ipynb from the stage modules.

The notebook imports the same modules the command-line scripts use, so there is
one implementation, not two.  Runnable top to bottom with no hidden state.
"""
import nbformat as nbf

nb = nbf.v4.new_notebook()
C = []
md = lambda s: C.append(nbf.v4.new_markdown_cell(s.strip()))
code = lambda s: C.append(nbf.v4.new_code_cell(s.strip()))


md(r"""
# PCM melt-front identifiability — two-day de-risk

**Question.** Is the PCM melt-front position — equivalently the remaining latent
buffer — identifiable from external temperature measurement alone?

**Why this exists.** A twelve-month collaboration is about to be proposed on
reconstructing the internal thermal state of a sodium-ion cell under
phase-change-material cooling. The preceding study on a lithium-ion cell ended
in a negative result: the problem was quasi-steady, a PDE solver added nothing
over a one-line algebraic relation, and a month went into an estimator the
physics did not require. The lesson was that identifiability must be checked in
the actual regime before anything is built. This notebook applies that lesson.

**Runtime** about 6 minutes on a 6-core CPU laptop. No GPU. Nothing is fitted.

---

## Why this regime might differ from the battery case

This is the reason the question is open rather than already answered.

In the battery case the internal state was a temperature field, and with a small
Biot number that field was slaved to the surface — quasi-steady, no dynamics,
nothing for a solver to recover. A melting PCM is different in kind. The
**remaining latent buffer is an integral of past heat absorbed**, so it carries
history that the instantaneous temperature field does not. A quasi-steady
relation cannot express it in principle.

That is the argument for why the proposal might survive the battery result.
This notebook tests whether the argument holds numerically.

There is a countervailing effect and it is the crux: during melting the PCM
absorbs heat at nearly constant temperature, so the external surface is *pinned*
near the melt point. That pinning is exactly why external sensing would be
uninformative during melting — the proposal's own motivating argument, turned
against it. The question is whether the front leaves any usable signature before
the PCM is exhausted, or whether the only observable event is the temperature
knee when melting completes.
""")

md(r"""
## Pre-registered predictions

Written and committed **before** Stage B was run (git `868b4ef`, ahead of the
results commit `62cbedb`), with the scoring conventions fixed in advance so that
a threshold could not be reinterpreted after seeing the numbers. Full text in
`PCM_PREDICTIONS.md`. Scored at the end of this notebook, misses included.

| # | Prediction | Threshold |
|---|---|---|
| Q1 | Melted fraction alone (S1) is identifiable | sd < 10 % of full buffer |
| Q2 | {f, q‴} is poorly conditioned | scaled cond > 10⁴ |
| Q3 | Three-case discrimination below 0.5 K during melting | max ΔT < 0.5 K |
| Q4 | Sensitivity ≥ 5× larger after full melt than during | ratio ≥ 5 |
| Q5 | A second sensor improves the buffer bound < 2× | factor < 2 |
| Q6 | Some high heating rate makes the front identifiable | exists |
""")

code("""
import json, numpy as np, matplotlib.pyplot as plt
import pcm_params as P, pcm_identify as I
np.set_printoptions(suppress=True, linewidth=140)
%matplotlib inline
""")

md(r"""
## Parameters and provenance

Every value is tagged **(a) sourced**, **(b) estimated with the assumption
stated**, or **(c) derived**. Nothing is invented. Two provenance points are
flagged in `pcm_params.py` because they could bias the answer, and both are
handled so they cut *against* a convenient result:

1. Rubitherm's 165 kJ/kg is latent **plus sensible** by their convention, used
   here as pure latent following the cited paper. This overstates the buffer,
   which is generous to the proposal.
2. The cited source gives **one** conductivity, not separate solid and liquid
   values. Setting `k_l = k_s` silently would manufacture a null result, since
   that contrast is the mechanism by which front position changes the composite
   resistance. It is swept in Stage E instead.
""")

code("""
mdl = P.build_model()
ts = P.timescales(mdl)
print(f"cell rho            {P.CELL_RHO:9.1f} kg/m3   (c) 84.2 g / 26700 volume")
print(f"q''' at 3C          {P.Q_3C:9.0f} W/m3    ({P.Q_3C*P.CELL_VOLUME:.3f} W)")
print(f"PCM mass            {ts['pcm_mass']*1e3:9.2f} g")
print(f"latent buffer       {ts['E_latent']:9.0f} J")
print(f"net melt power      {ts['Q_net']:9.3f} W   (input minus loss at T_m)")
print(f"tau_diff (PCM)      {ts['tau_diff']:9.0f} s")
print(f"Stefan number       {ts['Stefan']:9.4f}   latent-dominated")
print(f"Biot (PCM layer)    {ts['Biot_pcm']:9.4f}")
""")

md(r"""
## Stage A — build and verify the forward solver

A phase-change solver that has not been checked against an analytical benchmark
is not evidence of anything. If Stage A does not pass, nothing below it counts.

Three schemes solve the *same* discrete system: the literal pointwise apparent
heat capacity named in the proposal, an enthalpy-chord form, and Newton on the
enthalpy residual. Checks: the one-phase Stefan (Neumann) analytical solution,
convergence in dx / dt / mushy width, global energy closure, and latent closure.
""")

code("""
import pcm_stage_a_verify as A
A.a1_stefan_schemes()
_ = A.a2_convergence()
""")

md(r"""
Front error converges linearly to the sharp-interface limit as the mushy zone
narrows, and is flat in both mesh and timestep — so the residual benchmark error
is *entirely* mushy-zone width, which is the correct behaviour: the analytical
solution **is** the ΔT_m → 0 limit.
""")

code("""
A.a3_energy_closure()
_ = A.a4_cross_scheme()
""")

md(r"""
### Trap 1, demonstrated rather than assumed

Apparent heat capacity with too narrow a mushy zone steps over the latent
plateau and silently conserves no energy. The temperature field looks entirely
plausible while this happens; only the energy closure catches it.
""")

code("""
A.a5_pointwise_trap()
ok = A.verdict()
assert ok, "Stage A failed -- do not trust anything below this cell"
""")

md(r"""
## Stage B — the identifiability computation

Fisher information and Cramér–Rao bounds from central-difference sensitivities
through the **actual discretised solver**, not an analytic surrogate.

Melted fraction is a *state*, not a coefficient, so it needs a parameterisation
before it can be differentiated. The one used: run a reference melt from cold,
define `state(f0)` as the field on that trajectory where melted fraction equals
`f0`, and ask how well the next window of surface data determines `f0`.

Along a single trajectory `f0` co-varies with the whole temperature field, so
this framing carries **more** information than melted fraction alone — it is
generous to identifiability, and a negative result under it is conservative.
""")

code("""
import pcm_stage_b_crlb as B
ref = I.reference_trajectory()
B.b0_fd_convergence(ref)
""")

code("""
B.b1_main(ref)
""")

code("""
B.b4_window_length(ref)
B.b5_across_melt(ref)
""")

md(r"""
## Stage B audit — the numbers looked too good

Melted fraction to ~0.1 % of the buffer and latent heat to 0.015 % on a 0.1 K
sensor is a reason to look for a leak, not to celebrate. Four audits:

- **L1** — a history leak. Stage B built the initial state from the *nominal*-q
  trajectory then perturbed q only inside the observation window, so a plant
  that had always run at a different heat load was never represented.
- **L2** — conditioning is noise-independent; standard deviations are not.
- **L3** — where the information actually lives in the trace.
- **L4** — the control that matters: is the signal the front, or just warming?
""")

code("""
import pcm_stage_b_audit as AUD
_ = AUD.l1_history_leak()
""")

code("""
AUD.l2_noise_scaling(ref)
AUD.l3_information_shape(ref)
AUD.l4_front_vs_warming()
""")

md(r"""
**This is the pivot of the study.** Fixing the history leak takes S4's worst
correlation from 0.969 to **0.9993** — the same near-collinear structure that
warned in Part 7 (k vs h at 0.998) before a CRLB-blessed parameter set failed.

And L4 settles the mechanism: with latent heat cut 1000×, so that **no melt
front exists at all**, the melted fraction is recovered just as well. The
information is monotone warming encoding elapsed time on a known trajectory —
not the front. L3 says the same thing from the other side: 99 % of the Fisher
information is in the *mean level* of the trace, which a lookup table extracts
with no inverse PDE solve.
""")

md(r"""
## Stage C — the discrimination test

Three cases identical in every respect except melted fraction — 30 %, 50 %,
70 % — at 0.1 K noise.

Trap 6 requires the total heat input be held fixed. Every case here runs at the
same q‴, for the same window duration, through the same geometry, and differs
only in how far along the melt it starts. The cumulative heat needed to *reach*
each state necessarily differs — that is what melted fraction means — but
nothing inside the compared window differs.
""")

code("""
import pcm_stage_cd as CD
r_surface = CD.stage_c(ref, sensors=("surface",))
r_two     = CD.stage_c(ref, sensors=("surface", "can"), label="  [surface + can]")
""")

code("""
r_control = CD.stage_c_mechanism(ref)
""")

md(r"""
The control is the finding: a cell whose PCM has **no latent heat**, and
therefore no front, gives **more** separation than the real one. Whatever the
0.75 K separation is measuring, it is not the melt front.
""")

md(r"""
## Stage D — time-resolved sensitivity

Where the sensitivity lives determines what the method can honestly claim.
""")

code("""
t, sens, fbar, t_ex = CD.stage_d(ref)
""")

code("""
fig, ax = plt.subplots(figsize=(9, 4.2))
ax.plot(t, np.abs(sens), color="#C44E52", lw=1.6)
ax.axvline(t_ex, color="k", ls="--", lw=1.2)
ax.text(t_ex - 150, 1.02*np.max(np.abs(sens)), "PCM exhausted", ha="right", fontsize=9)
ax.set_yscale("log"); ax.grid(alpha=0.3)
ax.set_xlabel("time from window start (s)")
ax.set_ylabel(r"$|\\partial T_{surf}/\\partial f|$  (K per unit melted fraction)")
ax2 = ax.twinx(); ax2.plot(t, fbar, color="#8172B2", lw=1.1, alpha=0.65)
ax2.set_ylabel("melted fraction", color="#8172B2"); ax2.grid(False)
ax.set_title("Stage D -- sensitivity is flat through the melt, then jumps at exhaustion")
plt.show()
""")

md(r"""
## Stage E — the regime sweep

The most transferable output. Vary the heating rate so melting time spans more
than an order of magnitude around the PCM layer's diffusion time, and ask where
— if anywhere — the **front** becomes identifiable.

Raw separation will not answer this, because Stage C showed a latent-free cell
gives the same separation. The front-specific metric is the *excess*:
separation with real PCM minus separation with the latent heat removed, at the
same heating rate and the same fraction of each system's own melt duration.
""")

code("""
import pcm_stage_e_regime as E
rows = E.e1_regime_sweep()
""")

code("""
E.e2_k_ratio()
_ = E.e3_mushy_width()
""")

md(r"""
The front-attributable excess is **negative at every heating rate** and grows
more negative as melting accelerates. Latent heat compresses the temperature
excursion, so phase change makes the state *less* discriminable, never more.
Q6 predicted the opposite.

The two robustness sweeps close the flagged provenance gaps: the solid/liquid
conductivity ratio moves the bound only 0.102 % → 0.118 %, and the mushy-zone
width leaves separation flat at 0.69 K from 0.25 K to 2 K. Neither assumption is
carrying the conclusion.
""")

md(r"""
## Scored predictions and verdict

| # | Threshold | Actual | |
|---|---|---|---|
| Q1 | sd < 10 % | **0.101 %** of buffer | ✅ HIT |
| Q2 | cond > 10⁴ | **60.5** as computed, **14.6** leak-corrected | ❌ MISS, by ~3 orders |
| Q3 | max ΔT < 0.5 K | **0.751 K** worst pair (7.5σ) | ❌ MISS |
| Q4 | ratio ≥ 5 | **13.3×** | ✅ HIT |
| Q5 | factor < 2 | **1.997×** on S1, but 4.81× on S4 | ✅ HIT, marginal |
| Q6 | some rate works | **no rate, opposite direction** | ❌ MISS |

Three hits, three misses.

**Verdict.** The §5 thresholds for the strongest outcome are met on the numbers
(S1 sd 0.10 % < 10 %, discrimination 7.5σ > 3σ). But by mechanism this is the
second outcome: **identifiable by integration, not by front tracking.** The
buffer is recoverable because melted fraction is a monotone function of
accumulated net heat; it is not recoverable from the front's effect on thermal
resistance. Phase change does not merely fail to help — it actively suppresses
discriminability.

The honest claim is *"we detect PCM exhaustion and estimate the buffer by
integration"*, not *"we track the melt front continuously"*. The full write-up,
the load-dependence caveat that decides how much of this survives on real
hardware, and the note on which deck claims are supported, are in
`PCM_FINDINGS.md`.
""")

nb["cells"] = C
nb["metadata"] = {
    "kernelspec": {"display_name": "Python 3", "language": "python",
                   "name": "python3"},
    "language_info": {"name": "python", "version": "3.14"},
}
with open("PCM_Identifiability.ipynb", "w", encoding="utf-8") as fh:
    nbf.write(nb, fh)
print(f"wrote PCM_Identifiability.ipynb  ({len(C)} cells)")
