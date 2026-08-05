"""Generate Part7_CoreValidation.ipynb.

The notebook is the deliverable; this builder is the single source of truth for
its cells, so the two cannot drift.  Run, then execute the notebook to embed
outputs.
"""
import nbformat as nbf

nb = nbf.v4.new_notebook()
C = []


def md(s):
    C.append(nbf.v4.new_markdown_cell(s.strip("\n")))


def code(s):
    C.append(nbf.v4.new_code_cell(s.strip("\n")))


# --------------------------------------------------------------------------- #
md(r"""
# Part 7 — Validating a reconstructed internal temperature against a **measured** core

Six earlier notebooks validated this inverse PINN against a finite-difference truth model and
against a classical estimator on real data. None of them ever checked the reconstructed
**internal** field against an **actual internal measurement**, because sealed 18650 cells have
no core sensor.

This notebook closes that gap:

> **Fit the inverse PINN using the SURFACE trace only. Predict the CORE temperature.
> Compare against the MEASURED core temperature. Report the error in Kelvin.**

**Runtime warning.** The PINN stages take roughly 35–50 minutes total on a 6-core CPU laptop.
Stages A–C run in under a minute. Every cell is runnable top to bottom with no hidden state.

---

### A correction to the plan, made before any fitting

The session was designed around the University of Michigan Deep Blue release
*"An Electro-thermal Model for the A123 26650 LiFePO4 Battery"*. On inspection that package is
the **Simulink model**, not experimental data — `run_model.m` *simulates* `Tc` and `Ts` from a
block diagram and loads no measured file. Building on it would have meant validating a PINN
against another model, which is Part 2 again.

The data actually needed was already on disk from a different source: **Richardson & Howey
(2015)**, *"Sensorless battery internal temperature estimation using a Kalman filter with
impedance measurement"*, IEEE Trans. Sustainable Energy **6**(4):1190–1199 — two HEV drive
cycles on a **26650 A123 LFP cell with simultaneous core and surface thermocouples**.

The UMich package is still used, for one thing it genuinely provides: the **Forgez et al. (2010)
LFP entropy profile** (`dUdT.mat`), which lets the reversible term be handled with a real
citation instead of an invented number.
""")

code(r"""
import time, numpy as np, torch, matplotlib.pyplot as plt
import scipy.io as sio

from part7_lib import P, Record, RadialFV, load_forgez_dudt, bi_over_2
torch.set_num_threads(6)   # i7-8850H: 6 physical cores (12 logical)

plt.rcParams.update({"figure.dpi": 120, "font.size": 9, "axes.grid": True,
                     "grid.alpha": 0.25, "axes.axisbelow": True, "legend.frameon": False})
np.set_printoptions(suppress=True, linewidth=130)
print("torch", torch.__version__, "| threads", torch.get_num_threads())
""")

# ---------------------------------------------------------------- Stage A --- #
md(r"""
## Stage A — inspect before assuming

Load both records and print every channel's range, the sampling rate, and the record length.
Column meanings are **not** inferred from filenames: they are taken from the reference
implementation's own assignment (`MainScript.m:136-139`), which reads
`T_data(:,1)=t, (:,2)=T_surf, (:,3)=T_core, (:,4)=T_inf`.

**The gate:** if the measured core-minus-surface maximum is below about 0.3 K, the dataset
cannot support the claim and we stop.
""")

code(r"""
recs = {tag: Record(tag) for tag in ("1", "2")}
for tag, r in recs.items():
    s = r.summary()
    d = r.T_c - r.T_s
    print(f"=== DATASET {tag} " + "="*76)
    print(f"  duration {s['duration_s']:.0f} s, n={s['n']} on a {r.dt:g} s grid "
          f"(raw 1.1 s: {r.n_raw_T} temperature, {r.n_raw_VC} voltage/current samples)")
    print(f"  T_surf  {s['T_s_min']:7.3f} -> {s['T_s_max']:7.3f} degC")
    print(f"  T_core  {s['T_c_min']:7.3f} -> {s['T_c_max']:7.3f} degC")
    print(f"  T_inf   mean {s['T_inf_mean']:.3f} degC")
    print(f"  I       {s['I_min']:+7.2f} -> {s['I_max']:+7.2f} A  (rms {s['I_rms']:.2f})")
    print(f"  gross throughput {s['gross_Ah']:.2f} Ah = {s['gross_Ah']/P.Cap_Ah:.2f} "
          f"equivalent full cycles, but net DoD span only {s['dod_span']:.4f}")
    print(f"  >> core - surface: max {s['core_surf_max']:.4f} K, "
          f"mean {s['core_surf_mean']:.4f} K, final {s['core_surf_final']:.4f} K")
    print(f"  >> GATE (>= 0.3 K): {'PASS' if s['core_surf_max'] >= 0.3 else 'FAIL - STOP'}")
    rise = s['T_s_max'] - s['T_inf_mean']
    print(f"  >> max(core-surf)/surface-rise = {s['core_surf_max']/rise:.4f} "
          f"-> implies Bi ~ {2*s['core_surf_max']/rise:.3f}")
""")

md(r"""
Both records pass by a wide margin — 7.13 K and 6.55 K, not 0.3 K. Ambient is ~8 °C and peak
current is 30 A (≈13 C on a 2.3 Ah cell), which is why the gradient is so much larger than on
the earlier 18650 work.

Note the ratio `max(core−surf)/surface-rise` is **0.614 and 0.613** — two independent records
agreeing to three decimals. Read as Bi/2 that implies Bi ≈ 1.23, against Bi = 1.255 from
Richardson's independently identified `h`, `R`, `k`. **(c) inference**, and encouraging — but
Bi/2 is a *steady-state* identity and these are transients, so it is indicative, not exact.
""")

code(r"""
fig, axes = plt.subplots(2, 2, figsize=(11, 5.2), sharex="col")
for j, tag in enumerate(("1", "2")):
    r = recs[tag]; t = r.t/60
    ax = axes[0, j]
    ax.plot(t, r.T_c, "k-", lw=1.3, label="core (measured)")
    ax.plot(t, r.T_s, color="#1565c0", lw=1.2, label="surface (measured)")
    ax.plot(t, r.T_inf, ":", color="#777", lw=0.9, label="ambient")
    ax.set_title(f"dataset {tag}"); ax.set_ylabel("T / $^\\circ$C")
    if j == 0: ax.legend(fontsize=8)
    ax = axes[1, j]
    ax.plot(t, r.T_c - r.T_s, color="#c62828", lw=1.2)
    ax.set_ylabel("core $-$ surface / K"); ax.set_xlabel("time / min")
    ax.axhline(0.3, color="#888", ls="--", lw=0.8)
    ax.text(0.02, 0.9, f"max {np.max(r.T_c-r.T_s):.2f} K", transform=ax.transAxes, fontsize=8)
plt.tight_layout(); plt.show()
""")

md(r"""
### Assumed material properties, with provenance

Every number below is tagged **(a)** verified with source, **(b)** engineering estimate with
stated assumptions, or **(c)** inference. Nothing here is invented.
""")

code(r"""
print(f"(a) R_o    = {P.R_o*1000:.2f} mm          Richardson & Howey 2015, MainScript.m:87")
print(f"(a) V_cell = {P.V_b:.5e} m^3   MainScript.m:88")
print(f"(a) rho    = {P.rho:.1f} kg/m^3      MainScript.m:90")
print(f"(a) cp     = {P.cp:.1f} J/kg/K     MainScript.m:94 (identified on drive cycle 1)")
print(f"(a) k      = {P.k_t:.3f} W/m/K       MainScript.m:95 (radial, effective)")
print(f"(a) h      = {P.h_lit:.1f} W/m^2/K     MainScript.m:97")
print(f"(c) rho*cp = {P.rho_cp():.4e} J/m^3/K")
print(f"(c) C_lump = {P.C_lump():.2f} J/K")
print(f"(c) implied cell mass = rho*V = {P.rho*P.V_b*1000:.1f} g  "
      f"(datasheet ~76 g -> {100*(P.rho*P.V_b*1000-76)/76:+.1f} %)")
print(f"(c) Bi at literature h = {P.biot(P.h_lit):.4f}   Bi/2 = {bi_over_2(P.h_lit, P.k_t):.4f}")

soc, dudt = load_forgez_dudt()
print(f"\n(a) Forgez et al. 2010 LFP entropy profile, {len(soc)} points, "
      f"{1000*dudt.min():+.3f} to {1000*dudt.max():+.3f} mV/K")
print(f"    value at 50% SOC used by Richardson: {1000*P.dUdT_50:+.2f} mV/K")
""")

md(r"""
### Two traps caught here, both silent

**Unit trap.** `MainScript.m:121` documents the current channel as `I(mA)`. **That comment is
wrong.** Line 194 computes `Q = abs(I*(V-3.3))` with no factor of 1000, and the energy balance
settles it: at amperes the mean source is ~1.2 W against a lumped capacity of 84.5 J/K, which
matches the observed 18 K rise. At milliamps it would be 1.2 mW and the cell would not warm.

**Sign trap.** The current sign convention is established from the data, not assumed.
""")

code(r"""
for tag, r in recs.items():
    m = np.abs(r.I) > 1.0
    cc = np.corrcoef(r.I[m], r.V[m])[0, 1]
    Q = r.q_measured()
    print(f"dataset {tag}: corr(I,V)={cc:+.4f} over |I|>1A -> "
          f"I>0 means {'CHARGE' if cc > 0 else 'discharge'}")
    print(f"   regression V = {r.R_ohmic_reg*1000:+.3f} mOhm * I + {r.U_ocv_reg:.4f} V "
          f"(A123 26650 DC resistance ~10 mOhm at 25 C; higher at 8 C is expected)")
    print(f"   Q = I(V-U): mean {Q.mean():.4f} W, max {Q.max():.3f} W, "
          f"negative fraction {np.mean(Q<0):.3f}")
    # energy balance closes the argument
    A_out = 2*np.pi*P.R_o*P.height()
    need = (P.C_lump()*(0.5*(r.T_c+r.T_s)[-1]-0.5*(r.T_c+r.T_s)[0])
            + (P.h_lit*A_out*(r.T_s-r.T_inf)).sum()*r.dt) / r.t[-1]
    print(f"   energy balance needs mean {need:.3f} W; measured source gives "
          f"{Q.mean():.3f} W  ({100*abs(Q.mean()-need)/need:.1f} % apart)\n")
""")

md(r"""
Richardson writes `abs(I*(V-3.3))`, which hides the sign entirely. Writing it sign-explicitly as
`I(V − U)` exposed it: with the discharge-positive convention, 55 % of samples came out with
*negative* irreversible heat, which is physically impossible. `I > 0` is **charge**.

Because `V` is measured, `Q = I(V − U_ocv)` is a **directly measured** heat source with no free
parameter. It is not the core channel, so using it is not a leak — it is an independent
instrument, and it is what makes the energy-balance check above possible.
""")

# ---------------------------------------------------------------- Stage A2 -- #
md(r"""
### One finding that reshapes Stage E

These are **charge-sustaining** HEV cycles. Gross throughput is 3.4–3.6 equivalent full cycles,
but the net DoD excursion spans only **0.138 and 0.112** — the cell never leaves an ~11–14 %
window around its starting SOC.

So recovering a *shape* `R_eff(x)` across x ∈ [0,1] is **not supportable**: 86 % of the axis has
no observations behind it. This is confirmed quantitatively in Stage C rather than assumed.
""")

code(r"""
fig, ax = plt.subplots(1, 2, figsize=(11, 2.7))
for j, tag in enumerate(("1", "2")):
    r = recs[tag]
    ax[j].plot(r.t/60, r.dod, color="#6a1b9a", lw=1.2)
    ax[j].axhspan(r.dod.min(), r.dod.max(), color="#6a1b9a", alpha=0.12)
    ax[j].set_ylim(-0.55, 0.55)
    ax[j].set_title(f"dataset {tag}: DoD window = {r.dod.max()-r.dod.min():.3f} "
                    f"of the full [0,1] axis")
    ax[j].set_xlabel("time / min"); ax[j].set_ylabel("DoD (relative to t=0)")
plt.tight_layout(); plt.show()
""")

# --------------------------------------------------------------- Stage B ---- #
md(r"""
## Stage B — classical baseline first

The classical fitter is not busywork: on the previous cell, comparing the PINN's source term
against the classical fitter's is what exposed a missing `I(t)²` factor that no loss curve
revealed. **It is the instrument that detects PINN bugs.**

The forward model is a conservative finite-volume discretisation of radial conduction. Symmetry
at r = 0 holds *by construction* — the inner face of the first cell has zero area, so there is no
axis flux term to get wrong. It is verified against closed-form solutions to ~1e-13 before use.
""")

code(r"""
import subprocess, sys
print(subprocess.run([sys.executable, "tests_fd.py"], capture_output=True,
                     text=True).stdout[-1500:])
""")

md(r"""
Two bugs were caught by those checks, both of which would have been silent:

1. **The convective boundary source was broadcast to every cell** instead of only the outermost
   one, injecting heat throughout the domain. Global energy balance was off by 2.5×, and the
   error grew *linearly with N* — the signature that located it.
2. **The half-cell wall conductance ignored area variation**, an O(dr/4R) bias sitting in exactly
   the core-minus-surface quantity under test.

Now the fits. All of them use the **surface trace only**; the core column is read at the end
purely to score.
""")

code(r"""
print(subprocess.run([sys.executable, "stage_b_classical.py"], capture_output=True,
                     text=True).stdout)
""")

md(r"""
**The headline warning is [B3].** Fitting {h, k, ρc_p} gives the **best surface RMSE of any
model** (0.202 K, better than B1's 0.220 K) and the **worst core prediction** — 4.0 K against
B1's 0.40 K, a factor of 20. It buys 8 % on the surface and costs 10× on the core, with `k`
pinned at 5 W/m/K, physically absurd for a jellyroll.

That is the "smaller loss, worse answer" failure mode in its natural habitat, and it is
**completely invisible from the surface alone**. It is also the reason Stage C exists.

Two more results that matter:

- **[B7] control**: removing the `I(t)` factor entirely costs 7.1× on the surface fit. Trap 5.1
  would be caught on this dataset.
- **Cross-record**: `h` = 37.08 (DS1) vs 37.25 (DS2), 0.47 % apart. That reproducibility is what
  makes `h` legitimately independent when we fix it in Stage E.
""")

# --------------------------------------------------------------- Stage C ---- #
md(r"""
## Stage C — identifiability *before* fitting

If a parameter set is not identifiable, we do not fit it. Sensitivities are finite-differenced
through the **actual discretised model**, so the bound describes the estimator we really run.

The noise level matters and the first attempt got it wrong: the high-frequency estimator applied
to the 1 Hz-interpolated trace reports 0.005 K because interpolation correlates neighbours. On
the raw 1.1 s trace it is 0.0087 K. **Neither is the right σ here** — Stage B's best model leaves
0.22 K of visibly structured residual, ~25× the instrument noise. Model-form error, not sensor
error, sets recoverability.
""")

code(r"""
print(subprocess.run([sys.executable, "stage_c_crlb.py"], capture_output=True,
                     text=True).stdout)
""")

md(r"""
**{R_eff, h, k, ρc_p} is rank deficient** at cond ≈ 1.4×10⁹ / 1.7×10⁹, echoing the 8.8×10⁹ found
in earlier work. Fitting all four is ruled out.

**The CRLB is necessary but not sufficient**, and Stage B proves it. The bound rates
{h, k, ρc_p} as identifiable at 3.9 % worst relative sd — yet the actual B3 fit drove `k` from
0.404 to 5.0, a 1140 % move, roughly 300× outside the CRLB ellipse, and produced 4 K of core
error. A variance bound assumes the model is correct and the noise white; these residuals are
structured, so the optimiser chases *bias*, which no variance bound can see. The statistic that
**did** warn was the parameter correlation: 0.99 between `k` and `ρc_p`.

This is worth carrying forward as a method lesson: **report the correlation matrix, not just the
marginal standard deviations.**
""")

# ------------------------------------------------------- pre-registration --- #
md(r"""
## Pre-registered predictions

Recorded here **before** Stages D–F are run, and scored in Stage G including the misses.

| # | Prediction |
|---|---|
| P1 | Measured core-minus-surface maximum is between 0.5 K and 5 K |
| P2 | The classical fitter reaches surface RMSE ≤ 0.3 K |
| P3 | {R₀, shape, h} jointly is **not** identifiable from one channel — worst relative CRLB sd > 25 % |
| P4 | With h fixed, the forward PINN recovers the centre-surface difference to within 5 % of the finite-difference value |
| P5 | **Predicted core matches measured core to better than 1 K RMSE** — the headline |
| P6 | Predicted-core error exceeds surface-fit error by at least 2×, because the core is an extrapolation away from the data |
""")

# --------------------------------------------------------------- Stage D ---- #
md(r"""
## Stage D — forward gate

Solve the forward problem with **known** parameters and score the PINN against the verified FV
solution, specifically on the **centre-to-surface difference**. The trajectory is easy and
carries almost no spatial information; the gradient is the whole point.

### First attempt: plain full-field PINN — **FAILED at 39.2 %**

The PDE loss stalled at 2.68 while the dimensionless source has mean ≈ 2.3, so the residual was
the same size as the forcing: the PDE was essentially not being satisfied. Two causes, both
measured rather than guessed.
""")

code(r"""
# Why it failed, cause 1: the drive-cycle source is genuinely broadband.
r = recs["2"]
Qh = (r.I**2) * 0.0143398 / P.V_b * P.R_o**2 / (P.k_t * (r.T_s.max()-r.T_inf.mean()))
A = np.abs(np.fft.rfft(Qh - Qh.mean())); f = np.fft.rfftfreq(len(Qh), 1.0)
c = np.cumsum(A**2)/np.sum(A**2)
print("SOURCE bandwidth (this is what a Fourier-feature MLP would have to represent):")
for frac in (0.5, 0.9, 0.99):
    i = int(np.searchsorted(c, frac))
    print(f"  {100*frac:5.1f}% of source energy below harmonic #{i:5d} "
          f"(period {1/max(f[min(i,len(f)-1)],1e-9):6.1f} s)")
print(f"  n_freq=128 reaches only period {r.t[-1]/128:.1f} s -> needs n_freq ~ 1500. Infeasible.")

print("\nCause 2, measured on a constant-source control (sweep already run):")
print("  n_freq :     0      2      4      8     16     64")
print("  err %  : 0.086  0.247  0.150  0.231  0.493  4.883")
print("  -> Fourier features make it MONOTONICALLY WORSE. Zero features wins by 57x.")
print("     They cannot reach the needed bandwidth anyway, and they wreck the")
print("     optimisation landscape for the smooth solution we actually have.")
""")

md(r"""
### The fix: split off the analytically-known stiff part

This is **exact algebra, not an approximation**. Write `T(r,t) = T_l(t) + w(r,t)` where `T_l` is
the lumped 0-D solution. Substituting into the PDE, the spiky source **cancels identically**:

$$\rho c_p \frac{\partial w}{\partial t} = k\,\nabla^2 w + \frac{2h}{R}\,[T_l(t) - T_\infty]$$

with `w(r,0) = 0` exactly. The remaining forcing is proportional to `T_l − T_∞`, the output of a
first-order low-pass with τ = 405 s, which attenuates the 15 s content by ~169×. So `w` is smooth
and a tiny network represents it well.

**Is the PINN then doing trivial work?** No. `T_l` is a 0-D model containing *no radial
information whatsoever*; the entire core-to-surface gradient — the only quantity this project
cares about — lives in `w`. The split removes the stiff part we can integrate exactly and leaves
the PINN precisely the part we cannot.

The measured initial offset (trap 5.2) is carried by `T_l(0)`, and printed.
""")

code(r"""
print(subprocess.run([sys.executable, "stage_d_forward_split.py"], capture_output=True,
                     text=True).stdout)
""")

# --------------------------------------------------------------- Stage E ---- #
md(r"""
## Stage E — inverse on surface data only

**The core channel is held out entirely** — not for fitting, not for initialisation, not for
early stopping, not for model selection. It is read exactly once, in Stage F, to score.

### Leak control on the fixed parameters

Richardson's published `k = 0.404` and `h = 39.3` were identified on drive cycle 1 using **both
thermocouples**. Adopting them would import core-derived information into a supposedly
core-blind fit — a cross-record leak. Instead we take both from **Stage B [B2], which fitted
DS1's surface trace only**:

```
h = 37.0678 W/m²/K      k = 0.390843 W/m/K       (core-free, independent record)
```

and we fit **DS2**. These land within 3.3 % of Richardson's values, which is reassuring, but his
are not used.

Model selection is **truth-free**: the retained run is the one with the lowest PDE + BC residual
on a large *fixed* collocation set, never the one closest to the core. Six seeds; spread and
non-convergence rate both reported.
""")

md(r"""
**Cost control, stated explicitly.** A full 6-seed fit of both source orders takes roughly
**80 minutes** on this laptop (~6 min per seed). The two switches below are the only knobs; both
print what they did, so nothing is hidden.
""")

code(r"""
N_SEEDS      = 6      # >= 5 required. Set to 2 for a ~15 min smoke run.
REUSE_CACHED = True   # Load results/stage_e_shape*.npz if present instead of refitting.
                      # Set False to force a fresh fit. Either way it says which it did.

import os
from stage_e_inverse import run_inverse, summarise, H_FIXED, K_FIXED

rec2 = recs["2"]
print(f"h = {H_FIXED:.4f} W/m2/K, k = {K_FIXED:.6f} W/m/K  (Stage B [B2], DS1 surface only)")
print(f"Bi = {H_FIXED*P.R_o/K_FIXED:.4f}   Bi/2 = {H_FIXED*P.R_o/K_FIXED/2:.4f}")
print(f"measured core-surface max on DS2 = {(rec2.T_c-rec2.T_s).max():.4f} K\n")

stage_e = {}
for n_shape, lab in ((0, "order 0: scalar R_eff"), (1, "order 1: R_eff (1 + a1 x)")):
    path = f"results/stage_e_shape{n_shape}.npz"
    if REUSE_CACHED and os.path.exists(path):
        d = np.load(path)
        print(f"[LOADED cached fit from {path} -- set REUSE_CACHED=False to refit]")
        runs = [{"R_eff": float(d["R_eff"][i]), "shape": np.array([]),
                 "Tc": d["Tc"][i], "Ts": d["Ts"][i], "sel": float(d["sel"][i]),
                 "closures": 999, "surf_rmse": float(d["surf_rmse"][i])}
                for i in range(len(d["R_eff"]))]
    else:
        print(f"[FITTING {N_SEEDS} seeds for shape order {n_shape} -- this is the slow part]")
        runs = []
        for sd in range(N_SEEDS):
            t0 = time.time()
            r = run_inverse(rec2, n_shape=n_shape, seed=sd)
            print(f"   seed {sd}: {time.time()-t0:.0f} s, R_eff "
                  f"{1000*r['R_eff']:.4f} mOhm, surface {r['surf_rmse']:.4f} K")
            runs.append(r)
    summarise(rec2, runs, lab)
    stage_e[n_shape] = runs
""")

# --------------------------------------------------------------- Stage F ---- #
md(r"""
## Stage F — the validation

Predicted core against measured core: RMSE, maximum error, and error as a fraction of the
measured core-to-surface gradient. Plus the required figure.
""")

code(r"""
import stage_f_figure
stage_f_figure.main()
from IPython.display import Image, display
display(Image("figures/core_validation.png"))
display(Image("figures/gradient_validation.png"))
""")

md(r"""
### The like-for-like control, and the parameter that actually decides the answer

Trap 5.8 says compare like with like. The PINN must be measured against a baseline given the
**same information and the same fixed parameters** — not against a weaker one. The cell below
does that, then sweeps `k`, the parameter the surface data cannot identify.
""")

code(r"""
print(subprocess.run([sys.executable, "stage_f_classical_control.py"],
                     capture_output=True, text=True).stdout)
""")

md(r"""
Read the `k` sweep carefully — it is the most important table in this notebook. The **surface fit
improves monotonically with k**, while the **core error has a minimum and roughly doubles either
side**. At k = 0.55 the surface fit is the best in the table and the core prediction is nearly
twice as bad.

So the headline number rests on `k`, which is an **input**, not something this experiment
measured. Our leak-free k = 0.3908 (from a surface-only fit on the *other* record) lands within
about 1 % of the core-optimal value. That is either good transfer between records or luck, and a
single record cannot tell the difference. It is the first thing the next experiment should attack.
""")

# --------------------------------------------------------------- Stage G ---- #
md(r"""
## Stage G — honesty pass

### Is the PINN actually enforcing its own physics?

The finite-volume solver satisfies the PDE exactly by construction, so **its** surface RMSE is
the best achievable *within the physics*. If the PINN fits the surface better than that, it is
not being cleverer — it is spending PDE residual to buy data fit, because the constraint is soft.

That matters here specifically: a network that bends the physics to fit the surface has no reason
to extrapolate correctly to the core, which is the one thing we are asking of it. The test below
takes the PINN's own recovered `R_eff`, pushes it through the exact solver, and reports the
divergence in kelvin.
""")

code(r"""
print(subprocess.run([sys.executable, "stage_g_consistency.py"],
                     capture_output=True, text=True).stdout)
""")

md(r"""
### Does the data weight distort the recovered parameter?

Two independent anchors on the true `R_eff`: the classical fit with physics enforced exactly
(14.30 mΩ) and the electrical regression of V on I (14.49 mΩ). If the PINN's estimate drifts away
from those as the data weight rises, the soft constraint is biasing it.

**Note on selection.** If a weighting is preferred on the basis of this sweep, it is chosen on
*physics consistency* — PDE residual and agreement with the electrical anchor — both of which are
**core-blind**. Choosing a weighting because it improves the core error would be test-set
selection, and it would invalidate the headline.
""")

code(r"""
print(subprocess.run([sys.executable, "stage_g_weight_sweep.py"],
                     capture_output=True, text=True).stdout)
""")

md(r"""
### Verification suites

Both discretisations are checked against closed-form solutions before anything depends on them.
""")

code(r"""
for suite in ("tests_fd.py", "tests_split.py"):
    out = subprocess.run([sys.executable, suite], capture_output=True, text=True).stdout
    print(f"--- {suite} ---")
    print(out.strip().split("\n")[-2])
""")

md(r"""
### The written record

The scored prediction table, seed spread, non-convergence rate, and a plain statement of what is
and is not established are in `FINDINGS.md`, generated from the saved results so no number is
retyped by hand.
""")

code(r"""
import subprocess, sys
subprocess.run([sys.executable, "make_findings.py"], check=True)
print(open("FINDINGS.md", encoding="utf-8").read())
""")

nb["cells"] = C
nb.metadata.kernelspec = {"display_name": "Python 3", "language": "python", "name": "python3"}
nbf.write(nb, "Part7_CoreValidation.ipynb")
print(f"wrote Part7_CoreValidation.ipynb with {len(C)} cells")
