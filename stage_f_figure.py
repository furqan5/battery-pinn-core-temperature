"""Stage F -- score the core prediction and draw the one required figure.

Reads the Stage E results saved to results/, scores the selected run against the
MEASURED core, and produces:
    figures/core_validation.png  -- measured core, predicted core, measured
                                    surface on shared axes, plus the error trace.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from part7_lib import P, Record

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 160, "font.size": 9,
    "axes.grid": True, "grid.alpha": 0.25, "axes.axisbelow": True,
    "legend.frameon": False,
})


def load(n_shape):
    d = np.load(f"results/stage_e_shape{n_shape}.npz")
    return {k: d[k] for k in d.files}


def score(rec, Tc_pred):
    e = Tc_pred - rec.T_c
    grad = rec.T_c - rec.T_s
    return {
        "rmse": float(np.sqrt((e ** 2).mean())),
        "max": float(np.abs(e).max()),
        "bias": float(e.mean()),
        "rmse_pct_of_grad": 100 * float(np.sqrt((e ** 2).mean())) / float(grad.max()),
        "rmse_pct_of_grad_rms": 100 * float(np.sqrt((e ** 2).mean()))
                                / float(np.sqrt((grad ** 2).mean())),
    }


def main():
    rec = Record("2")
    grad = rec.T_c - rec.T_s

    rows = []
    picked = None
    for n_shape in (0, 1):
        try:
            d = load(n_shape)
        except FileNotFoundError:
            continue
        b = int(np.argmin(d["sel"]))              # truth-free selection
        Tc, Ts = d["Tc"][b], d["Ts"][b]
        s = score(rec, Tc)
        s["surf_rmse"] = float(np.sqrt(((Ts - rec.T_s) ** 2).mean()))
        s["ratio"] = s["rmse"] / s["surf_rmse"]
        s["R_eff"] = float(d["R_eff"][b]) * 1000
        s["n_shape"] = n_shape
        s["spread"] = float(d["core_rmse"].max() - d["core_rmse"].min())
        s["seed"] = b
        rows.append(s)
        if n_shape == 0:
            picked = (Tc, Ts, d, b)

    print("=" * 104)
    print("STAGE F -- PREDICTED CORE vs MEASURED CORE  (dataset 2, surface-only fit)")
    print("=" * 104)
    print(f"  measured core-minus-surface: max {grad.max():.4f} K, "
          f"rms {np.sqrt((grad**2).mean()):.4f} K, mean {grad.mean():.4f} K")
    print()
    hdr = (f"  {'model':<26s} {'R_eff':>9s} {'surfRMSE':>9s} {'CORE RMSE':>10s} "
           f"{'core max':>9s} {'bias':>8s} {'% of max grad':>14s} {'core/surf':>10s}")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for s in rows:
        name = f"order {s['n_shape']} (seed {s['seed']})"
        print(f"  {name:<26s} {s['R_eff']:>8.3f}m {s['surf_rmse']:>9.4f} "
              f"{s['rmse']:>10.4f} {s['max']:>9.4f} {s['bias']:>+8.4f} "
              f"{s['rmse_pct_of_grad']:>13.2f}% {s['ratio']:>10.2f}")
    print()
    for s in rows:
        print(f"  order {s['n_shape']}: P5 (core RMSE < 1 K) -> "
              f"{'HIT' if s['rmse'] < 1.0 else 'MISS'}   "
              f"P6 (core err >= 2x surface err) -> "
              f"{'HIT' if s['ratio'] >= 2.0 else 'MISS'} (ratio {s['ratio']:.2f})")

    # ---- where does the error actually live? ---- #
    print()
    print("  Error is NOT uniform in time -- decomposition of the order-0 result:")
    Tc0 = picked[0]
    e = Tc0 - rec.T_c
    for lo, hi, lab in ((0, 600, "first 10 min (start-up transient)"),
                        (600, 1200, "10-20 min"),
                        (600, int(rec.t[-1]), "after the first 10 min"),
                        (0, int(rec.t[-1]), "whole record")):
        m = (rec.t >= lo) & (rec.t < hi)
        print(f"    {lab:34s} RMSE {np.sqrt((e[m]**2).mean()):.4f} K   "
              f"max {np.abs(e[m]).max():.4f} K   bias {e[m].mean():+.4f} K")
    m_early = rec.t < 600
    print(f"    -> the first 10 min is {100*m_early.mean():.0f} % of the record but "
          f"{100*(e[m_early]**2).sum()/(e**2).sum():.0f} % of the squared error.")
    print("    Physically reasonable: the start-up transient has the sharpest radial")
    print("    gradients and the model assumes a perfectly uniform initial field.")

    # ---------------- figure ---------------- #
    Tc, Ts, d, b = picked
    t = rec.t / 60.0
    fig, (ax, ax2) = plt.subplots(
        2, 1, figsize=(7.6, 6.2), sharex=True,
        gridspec_kw={"height_ratios": [3, 1.15], "hspace": 0.08})

    # spread across seeds, shown as a band
    lo, hi = d["Tc"].min(axis=0), d["Tc"].max(axis=0)
    ax.fill_between(t, lo, hi, color="#c62828", alpha=0.18, lw=0,
                    label="predicted core, spread over 6 seeds")

    ax.plot(t, rec.T_c, color="#111111", lw=1.6, label="measured CORE (held out)")
    ax.plot(t, Tc, color="#c62828", lw=1.4, ls="--", label="predicted core (PINN)")
    ax.plot(t, rec.T_s, color="#1565c0", lw=1.3, label="measured SURFACE (the only input)")
    ax.plot(t, rec.T_inf, color="#777777", lw=0.9, ls=":", label="ambient")

    ax.set_ylabel("temperature  /  $^\\circ$C")
    # headroom so the legend never sits on the traces
    ylo, yhi = ax.get_ylim()
    ax.set_ylim(ylo, yhi + 0.30 * (yhi - ylo))
    ax.legend(loc="upper left", fontsize=8, ncol=2)
    s0 = rows[0]
    ax.set_title(
        "A123 26650 LFP: internal temperature reconstructed from surface data alone\n"
        f"core RMSE {s0['rmse']:.3f} K, max {s0['max']:.3f} K  "
        f"(measured core-surface gradient peaks at {grad.max():.2f} K)",
        fontsize=9.5)

    ax2.axhline(0, color="#444444", lw=0.8)
    ax2.plot(t, Tc - rec.T_c, color="#c62828", lw=1.1,
             label="predicted $-$ measured core")
    ax2.plot(t, Ts - rec.T_s, color="#1565c0", lw=1.0, alpha=0.85,
             label="surface fit residual")
    ax2.set_xlabel("time  /  minutes")
    ax2.set_ylabel("error  /  K")
    ax2.legend(loc="upper left", fontsize=8, ncol=2)
    ax2.set_xlim(t[0], t[-1])

    fig.savefig("figures/core_validation.png", bbox_inches="tight")
    print("  wrote figures/core_validation.png")

    # second figure: the gradient itself, which is the physics under test
    fig2, ax3 = plt.subplots(figsize=(7.6, 3.4))
    ax3.plot(t, grad, color="#111111", lw=1.5, label="measured core $-$ surface")
    ax3.plot(t, Tc - Ts, color="#c62828", lw=1.3, ls="--",
             label="predicted core $-$ surface")
    ax3.set_xlabel("time  /  minutes"); ax3.set_ylabel("$T_{core}-T_{surf}$  /  K")
    ax3.legend(loc="upper left", fontsize=8); ax3.set_xlim(t[0], t[-1])
    ax3.set_title("The quantity that carries the spatial information", fontsize=9.5)
    fig2.savefig("figures/gradient_validation.png", bbox_inches="tight")
    print("  wrote figures/gradient_validation.png")

    eg = (Tc - Ts) - grad
    print(f"  gradient RMSE {np.sqrt((eg**2).mean()):.4f} K "
          f"({100*np.sqrt((eg**2).mean())/np.sqrt((grad**2).mean()):.2f} % of rms gradient)")


if __name__ == "__main__":
    main()
