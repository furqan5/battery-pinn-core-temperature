"""Figure for G-1: the separated ohmic and reversible terms, with the K band."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from g1_entropic_separation import discharges, heat_trace, estimate_K

plt.rcParams.update({"figure.dpi": 130, "savefig.dpi": 250, "font.size": 8.5,
                     "axes.grid": True, "grid.alpha": 0.25,
                     "axes.axisbelow": True, "legend.frameon": False})

CELLS = ["B0038", "B0039", "B0040"]
AMB, CURRENTS = 44.0, [1.0, 2.0, 4.0]
XG = np.linspace(0.08, 0.92, 22)
T_K = AMB + 273.15
Iv = np.array(CURRENTS)
A = np.vstack([Iv, np.ones_like(Iv)]).T

recs = [r for c in CELLS for I in CURRENTS for r in discharges(c, AMB, I)]
K0, _ = estimate_K(recs, AMB)


def separate(Kf):
    pr = {}
    for I in CURRENTS:
        acc = []
        for c in CELLS:
            for r in discharges(c, AMB, I):
                x, Q = heat_trace(r, K0 * Kf, AMB)
                m = r["on"]
                if m.sum() < 25:
                    continue
                o = np.argsort(x[m])
                acc.append(np.interp(XG, x[m][o], (Q[m] / r["Imed"])[o]))
        pr[I] = np.nanmedian(np.array(acc), axis=0)
    R, dU = [], []
    for j in range(len(XG)):
        y = np.array([pr[I][j] for I in CURRENTS])
        s, c = np.linalg.lstsq(A, y, rcond=None)[0]
        R.append(s); dU.append(-c / T_K)
    return np.array(R), np.array(dU)


base_R, base_dU = separate(1.0)
band = [separate(f) for f in (0.7, 1.3)]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 2.7))

ax1.fill_between(XG, 1000 * np.minimum(band[0][0], band[1][0]),
                 1000 * np.maximum(band[0][0], band[1][0]),
                 color="#1565c0", alpha=0.18, lw=0, label="K varied $\\pm$30 %")
ax1.plot(XG, 1000 * base_R, "o-", color="#1565c0", ms=3, lw=1.4, label="separated")
i = int(np.argmin(base_R[4:])) + 4
ax1.plot(XG[i], 1000 * base_R[i], "*", ms=11, color="#c62828", zorder=5)
ax1.set_xlabel("depth of discharge  $x$")
ax1.set_ylabel("ohmic coefficient  $R(x)$  /  m$\\Omega$")
ax1.set_title("Irreversible term: stable", fontsize=9)
ax1.legend(fontsize=6.5, loc="upper right")

ax2.fill_between(XG, 1000 * np.minimum(band[0][1], band[1][1]),
                 1000 * np.maximum(band[0][1], band[1][1]),
                 color="#c62828", alpha=0.18, lw=0, label="K varied $\\pm$30 %")
ax2.plot(XG, 1000 * base_dU, "s-", color="#c62828", ms=3, lw=1.4, label="separated")
ax2.axhline(0, color="#444", lw=0.8)
ax2.set_xlabel("depth of discharge  $x$")
ax2.set_ylabel("$\\mathrm{d}U/\\mathrm{d}T$  /  mV K$^{-1}$")
ax2.set_title("Reversible term: not determined", fontsize=9)
ax2.legend(fontsize=6.5, loc="lower left")

fig.tight_layout()
fig.savefig("figures/g1_separation.png", bbox_inches="tight")

print(f"K0 = {K0:.5f} W/K")
print(f"R  : {1000*base_R.min():.1f} to {1000*base_R.max():.1f} mOhm, "
      f"median {1000*np.median(base_R):.1f}")
print(f"     band width at median x: "
      f"{1000*abs(band[0][0][11]-band[1][0][11]):.1f} mOhm "
      f"({100*abs(band[0][0][11]-band[1][0][11])/base_R[11]:.0f} % of value)")
print(f"dU/dT: {1000*base_dU.min():+.3f} to {1000*base_dU.max():+.3f} mV/K")
print(f"     band width at median x: "
      f"{1000*abs(band[0][1][11]-band[1][1][11]):.3f} mV/K "
      f"({100*abs(band[0][1][11]-band[1][1][11])/abs(base_dU[11]):.0f} % of value)")
print(f"interior minimum in R at x = {XG[i]:.2f} ({1000*base_R[i]:.1f} mOhm)")
print("wrote figures/g1_separation.png")
