# Pre-registered predictions — PCM melt-front identifiability

Written and committed **before** Stage B was run. Stage A (solver verification)
was complete at the time of writing; no identifiability number existed yet.
Scored in `PCM_FINDINGS.md`, misses included.

| # | Prediction | Threshold as scored |
|---|---|---|
| Q1 | Melted fraction alone (S1) is identifiable | relative CRLB sd < 10 % |
| Q2 | {melted fraction, q‴} (S3) is **poorly** conditioned — more heat into more PCM looks like less heat into less PCM | scaled condition number > 10⁴ |
| Q3 | Three-case discrimination (30/50/70 %) gives max separation **below 0.5 K** during melting, rising sharply after exhaustion | max |ΔT_surf| < 0.5 K during melt |
| Q4 | Sensitivity to melted fraction is **at least 5× larger** after full melt than during melting | ratio ≥ 5 |
| Q5 | A second sensor at the cell can improves the melted-fraction bound by **less than a factor of 2** — the limit is physical pinning, not sensor count | improvement factor < 2 |
| Q6 | There exists a heating-rate regime where the front is identifiable during melting, and it is at **high** rates where the mushy zone spans a real temperature range | exists q‴ with S1 sd < 10 % |

## Scoring conventions fixed in advance

So that a threshold cannot be reinterpreted after seeing the numbers:

- **"Relative CRLB sd" for melted fraction** is the standard deviation expressed
  as a **percentage of the full latent buffer** (i.e. sd in melted-fraction
  units × 100 %). This is the operationally meaningful quantity — it is how well
  the remaining buffer is known. The sd relative to the current value of *f* is
  also reported, but Q1 and Q6 are scored on percentage-of-full-buffer.
- **Scaled condition number** is cond(S) where column *j* of S is
  ∂T/∂θⱼ × scale(θⱼ), with scale = θⱼ for strictly positive parameters
  (q‴, h, L_f, k_pcm) and scale = 1 for the melted fraction, whose natural
  range is [0, 1]. Conditioning is noise-independent; the standard deviations
  scale as σ. Both are reported separately.
- **Nominal observation window** is 600 s at 1 Hz with σ = 0.1 K, starting from
  a state at f = 0.50 taken from the reference melt trajectory.
- **Q4's "during melting"** is the sensitivity averaged over the plateau
  (0.2 < f < 0.8); **"after full melt"** is the peak sensitivity in the
  200 s following exhaustion.

## Stated expectation of the mechanism

The argument the proposal rests on is that the remaining latent buffer is an
*integral of past heat absorbed*, so it carries history that an instantaneous
temperature field does not, and a quasi-steady relation cannot express it in
principle. That argument is correct as stated.

The countervailing effect, which is the crux: during melting the PCM absorbs
heat at nearly constant temperature, so the outer surface is **pinned** near the
melt point. Pinning is what makes external sensing uninformative *while* melting
proceeds. The regime numbers computed in Stage A say pinning should be strong
here — melting is ~23× slower than PCM diffusion (t_melt/τ_diff = 22.7) and the
Stefan number is 0.048, both deep in the latent-dominated, near-isothermal
limit. So Q3 and Q4 are expected to hold, and Q6 is the genuinely open one.
