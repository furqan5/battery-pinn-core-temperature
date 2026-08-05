"""
part7_lib.py -- shared physics + data plumbing for Part 7 (core-temperature validation).

Everything here is imported by Part7_CoreValidation.ipynb.  Kept in a module so the
notebook stays readable and so the numerics can be unit-tested independently.

Geometry/physics: 1-D radial conduction in a cylinder,

    rho*cp dT/dt = (k/r) d/dr( r dT/dr ) + q'''(t)
    dT/dr|_{r=0} = 0                       (symmetry)
    -k dT/dr|_{r=R} = h [ T(R,t) - T_inf ] (Robin)

Data: Richardson & Howey (2015), 26650 A123 LFP, two HEV drive cycles with
simultaneous core and surface thermocouples.
"""

from __future__ import annotations

import numpy as np
import scipy.io as sio
import scipy.linalg as sla

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #

REPO = r"C:\Users\Nouman\Desktop\Furqan's Docs\battery_pinn"

_RH = (REPO + r"\Data_Sets\EKF-Battery-Impedance-Temperature-master"
              r"\EKF-Battery-Impedance-Temperature-master"
              r"\Temperature Estimation using Impedance with Matlab - Github\Data"
              r"\MainScriptData")

_UM = (REPO + r"\Data_Sets\Battery Electrothermal Model-1"
              r"\Battery Electrothermal Model")


# --------------------------------------------------------------------------- #
# Cell / material parameters
# --------------------------------------------------------------------------- #

class P:
    """A123 ANR26650M1 parameters.

    Provenance tags follow the convention used throughout this project:
      (a) verified with source   (b) engineering estimate   (c) inference
    """
    # (a) Richardson & Howey 2015, MainScript.m lines 87-97
    R_o    = 0.0129        # m      outer radius
    V_b    = 3.4219e-5     # m^3    cell volume
    rho    = 2107.0        # kg/m^3 effective bulk density
    cp     = 1171.6        # J/kg/K effective bulk specific heat   (identified on DS1)
    k_t    = 0.404         # W/m/K  effective radial conductivity  (identified on DS1)
    h_lit  = 39.3          # W/m^2/K convection coefficient        (identified on DS1)

    # (a) A123 datasheet / Richardson
    Cap_Ah = 2.3           # Ah
    U_plat = 3.3           # V  LFP plateau OCV used by Richardson (MainScript.m:194)

    # (a) Forgez et al. 2010, J. Power Sources 195:2961-2968, via Richardson
    dUdT_50 = -0.5e-3      # V/K  entropic coefficient at 50% SOC

    @classmethod
    def rho_cp(cls) -> float:
        """(c) Volumetric heat capacity, J/m^3/K."""
        return cls.rho * cls.cp

    @classmethod
    def C_lump(cls) -> float:
        """(c) Lumped cell heat capacity, J/K."""
        return cls.rho_cp() * cls.V_b

    @classmethod
    def height(cls) -> float:
        """(c) Effective height implied by V_b and R_o, m."""
        return cls.V_b / (np.pi * cls.R_o ** 2)

    @classmethod
    def biot(cls, h: float) -> float:
        """(c) Biot number hR/k."""
        return h * cls.R_o / cls.k_t


# --------------------------------------------------------------------------- #
# Data loading
# --------------------------------------------------------------------------- #

class Record:
    """One HEV drive cycle resampled onto a uniform 1 Hz grid.

    Attributes are all 1-D arrays of equal length unless noted.
      t      s        time from start
      T_s    degC     measured SURFACE temperature      (the only T used for fitting)
      T_c    degC     measured CORE temperature         (HELD OUT -- validation only)
      T_inf  degC     measured ambient/coolant
      I      A        current, + = discharge
      V      V        terminal voltage
      dod    -        coulomb-counted depth of discharge relative to t=0
    """

    def __init__(self, tag: str, dt: float = 1.0):
        T = sio.loadmat(_RH + r"\Temperature_data.mat")["T_data_" + tag]
        VC = sio.loadmat(_RH + r"\Voltage_current_data.mat")["VC_data_" + tag]

        t_raw, Ts_raw, Tc_raw, Tinf_raw = T[:, 0], T[:, 1], T[:, 2], T[:, 3]
        ti, I_raw, V_raw = VC[:, 0], VC[:, 1], VC[:, 2]

        t_end = np.floor(min(t_raw.max(), ti.max()))
        self.tag = tag
        self.dt = dt
        self.t = np.arange(0.0, t_end + dt, dt)

        self.T_s = np.interp(self.t, t_raw, Ts_raw)
        self.T_c = np.interp(self.t, t_raw, Tc_raw)
        self.T_inf = np.interp(self.t, t_raw, Tinf_raw)
        self.I = np.interp(self.t, ti, I_raw)
        self.V = np.interp(self.t, ti, V_raw)

        # SIGN CONVENTION -- established empirically, not assumed:
        #   corr(I, V) = +0.90 over |I| > 1 A, i.e. voltage RISES with positive
        #   current, so I > 0 is CHARGE.  Confirmed by regressing V on I, which
        #   gives +13.9 mOhm (right for an A123 26650 at 8 C), and by the energy
        #   balance below.  Richardson's abs(I*(V-3.3)) hides this.
        # Therefore charge throughput is +int(I dt) and DEPTH OF DISCHARGE falls
        # as the cell charges.
        self.soc_shift = np.cumsum(self.I) * dt / 3600.0 / P.Cap_Ah
        self.dod = -self.soc_shift

        # measured plateau OCV for this record: intercept of V = R*I + U.
        # The SOC window here is only ~11% wide and LFP is flat, so a single
        # fitted plateau value is defensible; the slope doubles as an
        # independent anchor on R_eff.
        m = np.abs(self.I) > 1.0
        A = np.vstack([self.I[m], np.ones(int(m.sum()))]).T
        self.R_ohmic_reg, self.U_ocv_reg = np.linalg.lstsq(A, self.V[m], rcond=None)[0]

        self.n_raw_T = len(t_raw)
        self.n_raw_VC = len(ti)

    # -- derived heat-source estimates ------------------------------------- #

    def q_measured(self, U_ocv: float | None = None) -> np.ndarray:
        """Irreversible heat from measured terminal quantities, W.  No free parameter.

        With I > 0 = charge (see __init__), the Bernardi irreversible term is
            Q_irr = I_discharge (U_ocv - V) = -I (U_ocv - V) = I (V - U_ocv),
        which is >= 0 whenever the overpotential shares the sign of the current.
        We do NOT wrap this in abs(): a negative value is a real signal that the
        sign convention or the OCV is wrong, and abs() would swallow it.

        Validated by energy balance: mean here is 1.197 W (DS1) / 1.934 W (DS2)
        against the 1.226 W / 1.942 W required to explain the measured heating.
        """
        U = self.U_ocv_reg if U_ocv is None else U_ocv
        return self.I * (self.V - U)

    def q_reversible(self, T_K: np.ndarray, dUdT: float = P.dUdT_50) -> np.ndarray:
        """(a) Bernardi reversible term, W.  Forgez et al. 2010 via Richardson.

        Standard form is Q_rev = -I_discharge * T * dU/dT.  With I > 0 = charge,
        I_discharge = -I, so Q_rev = +I * T * dU/dT.  For LFP near 50% SOC
        dU/dT < 0, so this is endothermic on charge and exothermic on
        discharge, which is the expected behaviour.
        """
        return self.I * T_K * dUdT

    def summary(self) -> dict:
        d = self.T_c - self.T_s
        return {
            "tag": self.tag,
            "duration_s": float(self.t[-1]),
            "n": int(len(self.t)),
            "T_s_min": float(self.T_s.min()), "T_s_max": float(self.T_s.max()),
            "T_c_min": float(self.T_c.min()), "T_c_max": float(self.T_c.max()),
            "T_inf_mean": float(self.T_inf.mean()),
            "core_surf_max": float(d.max()),
            "core_surf_mean": float(d.mean()),
            "core_surf_final": float(d[-1]),
            "I_min": float(self.I.min()), "I_max": float(self.I.max()),
            "I_rms": float(np.sqrt((self.I ** 2).mean())),
            "dod_span": float(self.dod.max() - self.dod.min()),
            "gross_Ah": float(np.abs(self.I).sum() * self.dt / 3600.0),
        }


def load_forgez_dudt() -> tuple[np.ndarray, np.ndarray]:
    """(a) Full SOC-dependent LFP entropy profile, Forgez et al. 2010.

    Shipped inside the U. Michigan electro-thermal model release
    (dUdT.mat), which credits it to Forgez.  Returns (SOC, dUdT [V/K]).
    """
    m = sio.loadmat(_UM + r"\dUdT.mat")
    soc = m["SOC_ent"].ravel().astype(float)
    dudt = m["dUdt"].ravel().astype(float) * 1e-3     # mV/K -> V/K
    order = np.argsort(soc)
    return soc[order], dudt[order]


# --------------------------------------------------------------------------- #
# Finite-volume radial solver  (the truth model / bug detector)
# --------------------------------------------------------------------------- #

class RadialFV:
    """Conservative finite-volume discretisation of radial conduction in a cylinder.

    Cell-centred grid, N cells, node i centred at r_i = (i+1/2)*dr.

    Why finite volume rather than plain finite difference: the inner face of
    cell 0 has zero area, so the symmetry condition dT/dr|_0 = 0 is satisfied
    *identically* by construction rather than being imposed as an extra equation.
    That removes one whole class of silent error at r = 0.

    Time stepping is backward Euler (unconditionally stable).  At dt = 1 s with
    N = 40 an explicit scheme would be unstable -- the explicit limit here is
    about 0.32 s -- so this is not optional.
    """

    def __init__(self, N: int = 40, R: float = P.R_o):
        self.N, self.R = N, R
        self.dr = R / N
        self.r = (np.arange(N) + 0.5) * self.dr          # cell centres
        self.r_face = np.arange(N + 1) * self.dr         # faces, r_face[0] = 0

        # per-unit-height volumes and areas (H cancels in the balance)
        self.vol = np.pi * (self.r_face[1:] ** 2 - self.r_face[:-1] ** 2)
        self.area = 2.0 * np.pi * self.r_face            # area[0] = 0 -> symmetry

    # -- operator assembly -------------------------------------------------- #

    def _matrices(self, k: float, h: float, rho_cp: float, dt: float):
        """Build backward-Euler system  M T^{n+1} = rho_cp*vol/dt * T^n + b."""
        N, dr = self.N, self.dr
        main = np.zeros(N)
        lower = np.zeros(N - 1)
        upper = np.zeros(N - 1)

        cap = rho_cp * self.vol / dt
        main += cap

        # interior faces 1..N-1
        cond = k * self.area[1:N] / dr          # W/K per unit height, faces 1..N-1
        main[:-1] += cond
        main[1:] += cond
        upper[:] = -cond
        lower[:] = -cond

        # outer boundary: half-cell conduction in series with the convective film.
        #
        # The naive half-cell conductance k*A_R/(dr/2) ignores the area variation
        # between the last cell centre and the wall and carries an O(dr/4R) error.
        # Using the area-exact form below makes the scheme reproduce the analytic
        # uniform-generation steady state to machine precision, which turns the
        # steady-state test into a sharp instrument instead of a loose one.
        G_half = 4.0 * np.pi * k * self.R ** 2 / (self.R ** 2 - self.r[-1] ** 2)
        G_film = h * self.area[N]
        self.G_half = G_half
        self.G_film = G_film
        self.G_out = 1.0 / (1.0 / G_half + 1.0 / G_film)   # W/K per unit height
        main[-1] += self.G_out

        return lower, main, upper, cap

    # -- solve -------------------------------------------------------------- #

    def solve(self, t, q_vol, T_inf, T0, k=P.k_t, h=P.h_lit,
              rho_cp=None, return_field=False):
        """March the field forward.

        t        (n,)   uniform time grid, s
        q_vol    (n,)   volumetric heat generation, W/m^3
        T_inf    (n,)   ambient, degC
        T0       float or (N,)  initial temperature, degC
        returns  dict with T_core, T_surf (both degC), optionally full field
        """
        if rho_cp is None:
            rho_cp = P.rho_cp()
        n, N, dr = len(t), self.N, self.dr
        dt = float(t[1] - t[0])

        lower, main, upper, cap = self._matrices(k, h, rho_cp, dt)
        ab = np.zeros((3, N))
        ab[0, 1:] = upper
        ab[1, :] = main
        ab[2, :-1] = lower

        T = np.full(N, T0, dtype=float) if np.isscalar(T0) else np.array(T0, float)

        Tc = np.empty(n)
        Ts = np.empty(n)
        field = np.empty((n, N)) if return_field else None

        q_vol = np.asarray(q_vol, float)
        T_inf = np.asarray(T_inf, float)
        if T_inf.ndim == 0:
            T_inf = np.full(n, float(T_inf))

        Tc[0], Ts[0] = self._core(T), self._surf(T, T_inf[0])
        if return_field:
            field[0] = T

        for i in range(1, n):
            # NOTE: the convective source enters ONLY the outermost cell's equation.
            # Adding it to every cell (scalar broadcast) silently injects heat
            # everywhere and breaks the global energy balance -- error then grows
            # linearly with N, which is how this was caught.
            rhs = cap * T + q_vol[i] * self.vol
            rhs[-1] += self.G_out * T_inf[i]
            T = sla.solve_banded((1, 1), ab, rhs)
            Tc[i] = self._core(T)
            Ts[i] = self._surf(T, T_inf[i])
            if return_field:
                field[i] = T

        out = {"T_core": Tc, "T_surf": Ts, "G_out": self.G_out}
        if return_field:
            out["field"] = field
            out["r"] = self.r
        return out

    # -- observables -------------------------------------------------------- #

    def _core(self, T):
        """Centreline temperature, r = 0.

        Cell 0 is centred at dr/2, not at 0.  Near the axis the profile is
        locally parabolic, T(r) ~ T(0) + c r^2, so fit c from the first two
        cells and extrapolate inward.  With N = 40 the correction is small but
        it is not zero, and it is exactly the quantity under test here.
        """
        r0, r1 = self.r[0], self.r[1]
        c = (T[1] - T[0]) / (r1 ** 2 - r0 ** 2)
        return T[0] - c * r0 ** 2

    def _surf(self, T, T_inf):
        """Outer wall temperature T(R) -- where the thermocouple actually sits.

        Found by matching the area-exact half-cell conduction to the convective
        film.  This is NOT the last cell centre; reporting the cell centre as
        'surface' would bias the core-minus-surface difference low.
        """
        return ((self.G_half * T[-1] + self.G_film * T_inf)
                / (self.G_half + self.G_film))


# --------------------------------------------------------------------------- #
# Analytic checks
# --------------------------------------------------------------------------- #

def steady_core_minus_surf(q_vol: float, k: float, R: float = P.R_o) -> float:
    """(a) Exact steady solution for uniform generation in a cylinder.

    T(0) - T(R) = q''' R^2 / (4k).  Independent of h.
    """
    return q_vol * R ** 2 / (4.0 * k)


def steady_surf_minus_inf(q_vol: float, h: float, R: float = P.R_o) -> float:
    """(a) Steady surface rise: total generation over the convective area.

    q''' * pi R^2 H = h * 2 pi R H * (T(R)-T_inf)  ->  q''' R / (2h).
    """
    return q_vol * R / (2.0 * h)


def bi_over_2(h: float, k: float, R: float = P.R_o) -> float:
    """(a) Steady ratio [T(0)-T(R)] / [T(R)-T_inf] = hR/(2k) = Bi/2."""
    return h * R / (2.0 * k)
