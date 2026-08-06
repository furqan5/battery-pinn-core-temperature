"""Figure and criterion for the multi-rate timescale test."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"figure.dpi": 130, "savefig.dpi": 170, "font.size": 9,
                     "axes.grid": True, "grid.alpha": 0.25,
                     "axes.axisbelow": True, "legend.frameon": False})

d = np.load("results/mr_audit.npz")
rate, tau_ratio = d["rate"], d["tau_ratio"]
e_fit, e_fix = d["qs_fitted_h"], d["qs_fixed_bi"]

# 0.05C is excluded: R_eff rails at the bound and the modelled rise (0.13 K)
# cannot reproduce the measured 2.39 K, which at 0.125 A is not ohmic heating.
m = rate >= 1.0
x, yf, yx, rr = tau_ratio[m], e_fit[m] * 100, e_fix[m] * 100, rate[m]

print("=" * 88)
print("CRITERION — quasi-steady error against forcing/diffusion timescale ratio")
print("=" * 88)
print(f"  {'C-rate':>7s} {'t_dis/tau_diff':>15s} {'QS error %':>11s}")
for a, b, c in zip(rr, x, yf):
    print(f"  {a:>7.1f} {b:>15.2f} {c:>11.2f}")

# log-log fit: is the error a power law in the ratio?
p = np.polyfit(np.log(x), np.log(yf), 1)
print(f"\n  power-law fit: QS error % = {np.exp(p[1]):.2f} * (t/tau)^({p[0]:.3f})")
print(f"  r^2 = {np.corrcoef(np.log(x), np.log(yf))[0,1]**2:.4f}")


def ratio_at(target):
    """Interpolate the timescale ratio at which the error crosses `target` %."""
    o = np.argsort(x)
    return float(np.interp(target, yf[o][::-1], x[o][::-1]))


for tgt in (2.0, 5.0, 10.0, 25.0):
    print(f"  QS error = {tgt:5.1f} %  at  t_dis/tau_diff = {ratio_at(tgt):.2f}")

# Part 7 comparison, computed on the same normalisation
P7_RATIO = 3541.0 / 1043.0            # record length / tau_diff, record 2
P7_QS = 0.1409 / 5.1855 * 100         # QS core RMSE / rms measured gradient
print(f"\n  Part 7 drive cycle: t/tau = {P7_RATIO:.2f}, "
      f"measured QS error = {P7_QS:.2f} %")
print(f"  criterion predicts {np.exp(p[1])*P7_RATIO**p[0]:.2f} % pure approximation "
      f"error at that ratio.")
print(f"  Part 7's figure is larger because it is scored against a REAL core and")
print(f"  therefore also carries model-form and sensor error, not only the")
print(f"  quasi-steady approximation. The ordering is the right way round.")

fig, ax = plt.subplots(figsize=(7.2, 4.4))
ax.loglog(x, yf, "o-", color="#c62828", lw=1.6, ms=7,
          label="quasi-steady error, h fitted per record")
ax.loglog(x, yx, "s--", color="#1565c0", lw=1.1, ms=5, alpha=0.85,
          label="same, at a fixed Bi = 1.214 (h removed)")
for a, b, c in zip(rr, x, yf):
    ax.annotate(f"{a:g}C", (b, c), textcoords="offset points",
                xytext=(6, -11), fontsize=7.5, color="#555")

ax.axvline(1.0, color="#444", lw=0.9, ls=":")
# NOTE: ylim is set further down, so this y must sit inside the final limits.
# An earlier version placed it at y=0.55, below the axis, where it collided with
# the x-axis label in the compiled proof.
ax.text(1.06, 12.0, "forcing time = diffusion time", rotation=90,
        fontsize=7.5, color="#444", va="bottom")
ax.axhline(5.0, color="#2e7d32", lw=0.9, ls="--")
ax.text(x.min() * 1.05, 5.4, "5 % error", fontsize=7.5, color="#2e7d32")

ax.plot([P7_RATIO], [P7_QS], "*", ms=16, color="#000000", zorder=5,
        label="Part 7 drive cycle, scored against a MEASURED core")

ax.set_xlabel(r"forcing timescale / diffusion timescale,  $t_\mathrm{dis}/\tau_\mathrm{diff}$")
ax.set_ylabel("quasi-steady relation error  /  %")
ax.set_title("When does a transient solver earn its place?\n"
             "A123 26650 LFP, eight C-rates, same cell family as the core validation",
             fontsize=9.5)
ax.legend(loc="upper right", fontsize=8)
ax.set_ylim(1.4, 90)
fig.savefig("figures/timescale_criterion.png", bbox_inches="tight")
print("\n  wrote figures/timescale_criterion.png")
