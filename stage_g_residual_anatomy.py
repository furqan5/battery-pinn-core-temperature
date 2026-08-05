"""Stage G -- what structure does the classical model leave behind, and does the
PINN capture it?

The classical FV fit leaves 0.465 K of surface residual against 0.0097 K of
thermocouple noise, so ~98% of it is model-form error, not measurement error.
The PINN gets that down to 0.197 K.  Whatever the PINN is representing, it lives
in the difference.

Correlating the classical residual against candidate physical drivers says WHICH
idealisation is being violated.  This is diagnosis, not fitting.
"""
import numpy as np
from scipy.optimize import minimize_scalar

from part7_lib import P, Record, RadialFV
from stage_e_inverse import H_FIXED, K_FIXED

rec = Record("2")
fv = RadialFV(N=40)
T0 = float(rec.T_s[0])


def pred(R):
    return fv.solve(rec.t, (rec.I ** 2) * R / P.V_b, rec.T_inf, T0,
                    k=K_FIXED, h=H_FIXED, rho_cp=P.rho_cp())


r = minimize_scalar(lambda R: float(np.sum((pred(R)["T_surf"] - rec.T_s) ** 2)),
                    bounds=(1e-4, 0.2), method="bounded")
o = pred(r.x)
res = o["T_surf"] - rec.T_s

print("=" * 96)
print("CLASSICAL SURFACE RESIDUAL -- anatomy")
print("=" * 96)
print(f"  rms {res.std():.4f} K, mean {res.mean():+.4f} K, "
      f"range {res.min():+.3f} to {res.max():+.3f} K")
print(f"  thermocouple noise ~0.0097 K -> {100*(1-0.0097**2/res.var()):.1f} % of the "
      f"residual VARIANCE is model-form error, not measurement")

# autocorrelation time: white noise decorrelates in one sample
a = res - res.mean()
ac = np.correlate(a, a, "full")[len(a)-1:] / np.dot(a, a)
tau_ac = int(np.argmax(ac < 1/np.e))
print(f"  residual autocorrelation falls to 1/e at lag {tau_ac} s "
      f"-> strongly correlated, NOT white")

print()
print("  Correlation of the residual against candidate drivers:")
drivers = {
    "surface rise (T_s - T_inf)": rec.T_s - rec.T_inf,
    "instantaneous heat I^2":     rec.I ** 2,
    "|current|":                  np.abs(rec.I),
    "ambient T_inf":              rec.T_inf,
    "time":                       rec.t,
    "dT_s/dt":                    np.gradient(rec.T_s, rec.dt),
    "measured core-surface":      rec.T_c - rec.T_s,
}
for name, d in drivers.items():
    c = float(np.corrcoef(res, d)[0, 1])
    flag = "   <-- strong" if abs(c) > 0.5 else ""
    print(f"    {name:28s} r = {c:+.4f}{flag}")

print()
print("  Interpretation guide:")
print("    correlation with SURFACE RISE  -> h is not constant (temperature- or")
print("                                      flow-dependent convection)")
print("    correlation with I^2           -> the source model is wrong in shape")
print("    correlation with dT_s/dt       -> a missing thermal mass, e.g. the can")
print("                                      or the thermocouple's own dynamics")
print("    correlation with TIME          -> a slow drift (ageing, flow change)")

# Does allowing h to depend on surface rise absorb it?  Cheap classical test.
print()
print("  Test: let h vary linearly with surface rise, h = h0 (1 + g*(T_s - T_inf)).")
print("  This stays a CLASSICAL fit -- physics enforced exactly -- and uses only")
print("  surface and ambient data.")
print(f"    {'g':>8s} {'h0':>8s} {'R_eff mOhm':>11s} {'surf RMSE':>10s} {'CORE RMSE':>10s}")


class VarH(RadialFV):
    """FV with h depending on the measured surface rise (still surface-only data)."""

    def solve_varh(self, t, q_vol, T_inf, T0, h_series, k, rho_cp):
        import scipy.linalg as sla
        n, N, dr = len(t), self.N, self.dr
        dt = float(t[1] - t[0])
        Tc_o = np.empty(n); Ts_o = np.empty(n)
        T = np.full(N, T0, float)
        cap = rho_cp * self.vol / dt
        cond = k * self.area[1:N] / dr
        G_half = 4.0 * np.pi * k * self.R ** 2 / (self.R ** 2 - self.r[-1] ** 2)
        Tc_o[0] = self._core(T)
        self.G_half = G_half
        self.G_film = h_series[0] * self.area[N]
        Ts_o[0] = self._surf(T, T_inf[0])
        for i in range(1, n):
            G_film = h_series[i] * self.area[N]
            G_out = 1.0 / (1.0 / G_half + 1.0 / G_film)
            main = cap.copy()
            main[:-1] += cond; main[1:] += cond; main[-1] += G_out
            ab = np.zeros((3, N))
            ab[0, 1:] = -cond; ab[1, :] = main; ab[2, :-1] = -cond
            rhs = cap * T + q_vol[i] * self.vol
            rhs[-1] += G_out * T_inf[i]
            T = sla.solve_banded((1, 1), ab, rhs)
            self.G_film = G_film
            Tc_o[i] = self._core(T); Ts_o[i] = self._surf(T, T_inf[i])
        return {"T_core": Tc_o, "T_surf": Ts_o}


vh = VarH(N=40)
rise = rec.T_s - rec.T_inf
for g in (0.0, 0.01, 0.02, 0.03, 0.05):
    hs = H_FIXED * (1.0 + g * rise)

    def p2(R):
        return vh.solve_varh(rec.t, (rec.I ** 2) * R / P.V_b, rec.T_inf, T0,
                             hs, K_FIXED, P.rho_cp())
    rr = minimize_scalar(lambda R: float(np.sum((p2(R)["T_surf"] - rec.T_s) ** 2)),
                         bounds=(1e-4, 0.2), method="bounded")
    oo = p2(rr.x)
    es = np.sqrt(((oo["T_surf"] - rec.T_s) ** 2).mean())
    ec = np.sqrt(((oo["T_core"] - rec.T_c) ** 2).mean())
    print(f"    {g:>8.3f} {H_FIXED:>8.3f} {1000*rr.x:>11.4f} {es:>10.4f} {ec:>10.4f}")
print("  (g is scanned to DIAGNOSE, not selected on the core -- the headline keeps")
print("   constant h.)")
