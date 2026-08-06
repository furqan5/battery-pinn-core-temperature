"""Figure: the surface fit and the core error move in OPPOSITE directions with k.

This is Table 6.7 as a figure. The crossing behaviour is the point and it reads
far better visually: the observable quantity improves monotonically while the
quantity of interest passes through a minimum.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

from part7_lib import P, Record, RadialFV
from stage_e_inverse import H_FIXED, K_FIXED

plt.rcParams.update({"figure.dpi": 130, "savefig.dpi": 300, "font.size": 8.5,
                     "axes.grid": True, "grid.alpha": 0.25,
                     "axes.axisbelow": True, "legend.frameon": False})

rec = Record("2")
fv = RadialFV(N=40)
T0 = float(rec.T_s[0])

ks = np.linspace(0.28, 0.60, 17)
surf, core = [], []
for k in ks:
    def pred(R):
        return fv.solve(rec.t, (rec.I ** 2) * R / P.V_b, rec.T_inf, T0,
                        k=k, h=H_FIXED, rho_cp=P.rho_cp())
    r = minimize_scalar(lambda R: float(np.sum((pred(R)["T_surf"] - rec.T_s) ** 2)),
                        bounds=(1e-4, 0.2), method="bounded")
    o = pred(r.x)
    surf.append(float(np.sqrt(((o["T_surf"] - rec.T_s) ** 2).mean())))
    core.append(float(np.sqrt(((o["T_core"] - rec.T_c) ** 2).mean())))
surf, core = np.array(surf), np.array(core)

fig, ax = plt.subplots(figsize=(3.5, 2.7))
ax2 = ax.twinx()
l1, = ax.plot(ks, surf, "o-", color="#1565c0", ms=3.5, lw=1.4,
              label="surface RMSE (observable)")
l2, = ax2.plot(ks, core, "s-", color="#c62828", ms=3.5, lw=1.4,
               label="core RMSE (what we want)")
ax.axvline(K_FIXED, color="#2e7d32", lw=1.0, ls="--")
ax.text(K_FIXED + 0.006, ax.get_ylim()[1] * 0.97, "k used\n(leak-free)",
        fontsize=6.5, color="#2e7d32", va="top")

ax.set_xlabel(r"radial conductivity $k$  /  W m$^{-1}$K$^{-1}$")
ax.set_ylabel("surface RMSE  /  K", color="#1565c0")
ax2.set_ylabel("core RMSE  /  K", color="#c62828")
ax.tick_params(axis="y", labelcolor="#1565c0")
ax2.tick_params(axis="y", labelcolor="#c62828")
ax2.grid(False)
ax.legend(handles=[l1, l2], loc="upper center", fontsize=6.5)
fig.tight_layout()
fig.savefig("figures/biot_sensitivity.png", bbox_inches="tight")

i = int(np.argmin(core))
print(f"surface RMSE is monotone decreasing: {surf[0]:.4f} -> {surf[-1]:.4f} K")
print(f"core RMSE minimum at k = {ks[i]:.4f} ({core[i]:.4f} K)")
print(f"  at k = 0.28: core {core[0]:.4f} K ; at k = 0.60: core {core[-1]:.4f} K")
print(f"  leak-free k = {K_FIXED:.5f} gives core "
      f"{np.interp(K_FIXED, ks, core):.4f} K, "
      f"{100*abs(K_FIXED-ks[i])/ks[i]:.1f} % from the optimum")
print("wrote figures/biot_sensitivity.png")
