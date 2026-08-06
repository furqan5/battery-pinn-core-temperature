"""Figures for the PCM identifiability study.

  fig_pcm_crlb.png       CRLB across S1-S4, with the leak correction shown
  fig_pcm_sensitivity.png Stage D time-resolved dT_surf/d(melted fraction)
  fig_pcm_regime.png     Stage E heating-rate sweep (the transferable one)
"""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.size": 9, "axes.grid": True, "grid.alpha": 0.3,
    "figure.dpi": 150, "savefig.dpi": 150, "savefig.bbox": "tight",
})

B = json.load(open("results/pcm_stage_b.json"))
A = json.load(open("results/pcm_stage_b_audit.json"))
CD = json.load(open("results/pcm_stage_cd.json"))
E = json.load(open("results/pcm_stage_e.json"))

NOISE = 0.1


# --------------------------------------------------------------------------- #

def fig_crlb():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.2))

    labels = ["S1\n{f}", "S3\n{f, q}", "S4\n{f, q, h}"]
    one = [100 * B["B1_surface"][k]["sd_natural"][0] for k in ("S1", "S3", "S4")]
    two = [100 * B["B2_two_sensor"][k]["sd_natural"][0] for k in ("S1", "S3", "S4")]
    corr = [100 * B["B1_surface"]["S1"]["sd_natural"][0],
            A["L1"]["S3"]["state consistent with q"]["f_sd_buffer_pct"],
            A["L1"]["S4"]["state consistent with q"]["f_sd_buffer_pct"]]

    x = np.arange(3)
    w = 0.27
    ax1.bar(x - w, one, w, label="1 sensor (surface)", color="#4C72B0")
    ax1.bar(x, two, w, label="2 sensors (+ can)", color="#55A868")
    ax1.bar(x + w, corr, w, label="1 sensor, history leak fixed",
            color="#C44E52")
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels)
    ax1.set_ylabel("melted-fraction sd  (% of full latent buffer)")
    ax1.set_yscale("log")
    ax1.axhline(10.0, color="k", ls="--", lw=1)
    ax1.text(2.42, 11.0, "Q1 threshold 10 %", ha="right", fontsize=8)
    ax1.set_title("Cramér–Rao bound on the remaining buffer\n"
                  "σ = 0.1 K, 1 Hz, 600 s window at f = 0.50", fontsize=9.5)
    ax1.legend(fontsize=7.5, loc="upper left")
    for xi, v in zip(x - w, one):
        ax1.text(xi, v * 1.15, f"{v:.2f}", ha="center", fontsize=7)
    for xi, v in zip(x + w, corr):
        ax1.text(xi, v * 1.15, f"{v:.2f}", ha="center", fontsize=7)

    conds = [B["B1_surface"][k]["cond"] for k in ("S1", "S3", "S4")]
    conds_c = [B["B1_surface"]["S1"]["cond"],
               A["L1"]["S3"]["state consistent with q"]["cond"],
               A["L1"]["S4"]["state consistent with q"]["cond"]]
    ax2.bar(x - 0.17, conds, 0.34, label="as computed in Stage B",
            color="#4C72B0")
    ax2.bar(x + 0.17, conds_c, 0.34, label="history leak fixed", color="#C44E52")
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels)
    ax2.set_yscale("log")
    ax2.set_ylabel("scale-normalised cond(S)")
    ax2.axhline(1e4, color="k", ls="--", lw=1)
    ax2.text(2.45, 1.3e4, "Q2 predicted > 10⁴", ha="right", fontsize=8)
    ax2.set_ylim(0.5, 1e5)
    ax2.set_title("Conditioning — noise-independent\n"
                  "(the sd panel scales with σ; this one does not)", fontsize=9.5)
    ax2.legend(fontsize=7.5, loc="upper left")
    for xi, v in zip(x - 0.17, conds):
        ax2.text(xi, v * 1.3, f"{v:.0f}", ha="center", fontsize=7)
    for xi, v in zip(x + 0.17, conds_c):
        ax2.text(xi, v * 1.3, f"{v:.0f}", ha="center", fontsize=7)

    s4 = A["L1"]["S4"]["state consistent with q"]
    mx = max(abs(v) for row in s4["corr"] for v in row if abs(v) < 0.999999)
    fig.text(0.5, -0.06,
             f"S4 is the realistic field case. Once the history is made "
             f"consistent with the trial heat load, corr(f, h) = {mx:.4f} — the "
             f"same near-collinear structure that warned in Part 7.",
             ha="center", fontsize=8.5, style="italic")
    fig.tight_layout()
    fig.savefig("figures/fig_pcm_crlb.png")
    print("  wrote figures/fig_pcm_crlb.png")
    plt.close(fig)


def fig_sensitivity():
    d = CD["D"]
    t = np.array(d["t"])
    s = np.array(d["sens"])
    f = np.array(d["f"])
    t_ex = d["t_exhaustion"]

    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    ax.plot(t, np.abs(s), color="#C44E52", lw=1.6,
            label=r"$|\partial T_{surf}/\partial f|$")
    ax.axvline(t_ex, color="k", ls="--", lw=1.2)
    ax.text(t_ex - 150, 1.02 * np.max(np.abs(s)), "PCM exhausted",
            fontsize=8.5, ha="right")

    ax.axhline(d["mean_during"], color="#4C72B0", ls=":", lw=1.4)
    ax.text(120, d["mean_during"] * 1.25,
            f"mean during melting = {d['mean_during']:.2f} K", fontsize=8,
            color="#4C72B0")
    ax.annotate(f"peak after exhaustion = {d['peak_after']:.1f} K\n"
                f"ratio {d['ratio']:.1f}×",
                xy=(t_ex, d["peak_after"]),
                xytext=(t_ex + 1500, d["peak_after"] * 0.55),
                fontsize=8.5, ha="left",
                arrowprops=dict(arrowstyle="->", lw=1))

    ax.set_xlabel("time from window start (s)")
    ax.set_ylabel(r"$|\partial T_{surf}/\partial f|$   (K per unit melted fraction)")
    ax.set_yscale("log")
    ax.set_title("Stage D — where the sensitivity to melted fraction lives\n"
                 "started from f = 0.30, propagated through exhaustion",
                 fontsize=10)

    ax2 = ax.twinx()
    ax2.plot(t, f, color="#8172B2", lw=1.1, alpha=0.65)
    ax2.set_ylabel("melted fraction", color="#8172B2")
    ax2.tick_params(axis="y", labelcolor="#8172B2")
    ax2.grid(False)
    ax2.set_ylim(0, 1.05)

    ax.legend(loc="upper left", fontsize=8.5)
    fig.text(0.5, -0.04,
             "The sensitivity is flat and modest for the whole melt, then jumps "
             "an order of magnitude at exhaustion. Exhaustion is an event the "
             "surface reports clearly; the front's position, during melting, is not.",
             ha="center", fontsize=8.5, style="italic")
    fig.tight_layout()
    fig.savefig("figures/fig_pcm_sensitivity.png")
    print("  wrote figures/fig_pcm_sensitivity.png")
    plt.close(fig)


def fig_regime():
    rows = E["E1"]
    tt = np.array([r["t_melt_over_tau"] for r in rows])
    sep = np.array([r["sep_pcm"] for r in rows])
    sep0 = np.array([r["sep_no_latent"] for r in rows])
    exc = np.array([r["front_excess_over_noise"] for r in rows])
    cr = np.array([r["c_rate"] for r in rows])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.2))

    # Axis is inverted so that FASTER melting is to the right, which is the
    # direction Q6 predicted the front would become observable in.
    def _regime_axis(ax):
        ax.set_xscale("log")
        ax.invert_xaxis()
        ax.set_xticks([30, 10, 3, 1])
        ax.set_xticklabels(["30", "10", "3", "1"])
        ax.minorticks_off()
        ax.set_xlabel(r"$t_{melt}/\tau_{diff}$"
                      "\n" r"faster melting $\rightarrow$")

    ax1.plot(tt, sep, "o-", color="#4C72B0", label="real PCM (RT42)")
    ax1.plot(tt, sep0, "s--", color="#937860",
             label="control: latent heat removed")
    ax1.axhline(3 * NOISE, color="k", ls=":", lw=1)
    ax1.text(tt.min(), 3 * NOISE * 1.12, "3σ", ha="left", fontsize=8)
    ax1.set_ylabel("worst pairwise separation, 30/50/70 % (K)")
    ax1.set_title("Discrimination vs heating regime", fontsize=9.5)
    ax1.legend(fontsize=8, loc="upper left")
    _regime_axis(ax1)

    ax2.plot(tt, exc, "o-", color="#C44E52")
    ax2.axhline(0, color="k", lw=1)
    ax2.axhspan(-3, 3, color="grey", alpha=0.18)
    ax2.set_ylim(min(exc.min() * 1.15, -5), 8)
    ax2.text(tt.max(), 4.4, "±3σ band: front contributes nothing usable",
             ha="left", fontsize=8)
    ax2.set_ylabel("front-attributable excess separation / σ")
    ax2.set_title("Signal attributable to phase change\n"
                  "(real PCM minus latent-free control)", fontsize=9.5)
    _regime_axis(ax2)
    for x, y, c in zip(tt, exc, cr):
        if c in (cr[0], cr[-1]):
            ax2.annotate(f"{c:.0f}C equiv.", xy=(x, y), xytext=(0, 10),
                         textcoords="offset points", fontsize=8, ha="center")

    fig.text(0.5, -0.05,
             "The excess is negative at every heating rate and grows more "
             "negative as melting speeds up: latent heat compresses the "
             "temperature excursion, so phase change makes the state less "
             "discriminable, never more. Q6 predicted the opposite.",
             ha="center", fontsize=8.5, style="italic")
    fig.tight_layout()
    fig.savefig("figures/fig_pcm_regime.png")
    print("  wrote figures/fig_pcm_regime.png")
    plt.close(fig)


if __name__ == "__main__":
    import os
    os.makedirs("figures", exist_ok=True)
    fig_crlb()
    fig_sensitivity()
    fig_regime()
