"""1D composite phase-change conduction solver.

Geometry: core -> can wall -> PCM -> ambient, planar (m=0) or cylindrical (m=1).

    rho c_p(T) dT/dt = (1/r^m) d/dr ( r^m k(T) dT/dr ) + q'''

Phase change is carried in the volumetric enthalpy

    H(T) = rho * [ cp * (T - T_ref) + L_f * g(T) ]

with g the liquid fraction ramping linearly over the mushy zone
[T_m - dTm/2, T_m + dTm/2].  Differentiating gives the apparent heat capacity
named in the proposal,

    c_p,eff = cp + L_f / dTm     inside the mushy zone,

so the three schemes below are three ways of solving the SAME continuous model,
not three different models.  They differ only in how the nonlinearity is handled:

  'ahc_pointwise'   c_eff evaluated at the current iterate.  This is the literal
                    textbook apparent-heat-capacity method.  It LOSES ENERGY when
                    a cell steps across the whole mushy zone in one timestep,
                    because no iterate ever sees the latent spike.  Kept
                    deliberately so the failure can be shown rather than assumed.
  'ahc_chord'       c_eff as the enthalpy chord (H(T^k)-H(T^n))/(T^k-T^n).  At
                    Picard convergence this reproduces H(T^{n+1})-H(T^n) exactly,
                    so it conserves energy to solver tolerance.
  'enthalpy_newton' Newton on the enthalpy residual directly.  An independent
                    algorithm (different iteration matrix and convergence path)
                    reaching the same discrete nonlinear system -- this
                    cross-checks the IMPLEMENTATION, not the model form.  The
                    load-bearing checks are the Stefan benchmark and energy
                    closure in pcm_stage_a_verify.py.

'enthalpy_newton' is the DEFAULT and is what Stages B-E run on.  Not because the
chord scheme is wrong -- the two discretise the identical system and agree to
<1e-6 K when both are run to convergence -- but because Picard convergence
degrades badly right at the liquidus: a cell that has just finished melting has
an enthalpy chord ~100x its sensible capacity, giving a contraction ratio near
0.94 and needing ~300 iterations per timestep.  Newton gets there in 3-4.

Everything is finite volume with harmonic-mean face conductivities, which is the
correct averaging across a material jump (core/can/PCM) and matters here because
the conductivity contrast is ~100x at the can/PCM interface.
"""
import numpy as np
from scipy.linalg import solve_banded

T_REF = 298.15  # enthalpy datum, K.  Arbitrary; cancels in all energy balances.


# --------------------------------------------------------------------------- #
# Materials
# --------------------------------------------------------------------------- #

class Material:
    """A homogeneous material, optionally undergoing phase change.

    k_liquid defaults to k (solid).  Giving them different values is what makes
    the melt-front position affect the composite thermal resistance at all --
    with k_s == k_l a moving front is very nearly invisible from outside, so
    this is not a cosmetic option.
    """

    def __init__(self, name, rho, cp, k, Lf=0.0, Tm=None, dTm=None, k_liquid=None):
        self.name = name
        self.rho = rho
        self.cp = cp
        self.k = k
        self.Lf = Lf
        self.Tm = Tm
        self.dTm = dTm
        self.k_liquid = k if k_liquid is None else k_liquid
        self.pcm = Lf > 0.0
        if self.pcm:
            if Tm is None or dTm is None or dTm <= 0:
                raise ValueError(f"{name}: phase-change material needs Tm and dTm>0")
            self.Ts = Tm - 0.5 * dTm
            self.Tl = Tm + 0.5 * dTm

    def g(self, T):
        """Liquid fraction in [0, 1]."""
        if not self.pcm:
            return np.zeros_like(T)
        return np.clip((T - self.Ts) / (self.Tl - self.Ts), 0.0, 1.0)

    def H(self, T):
        """Volumetric enthalpy, J/m^3."""
        h = self.rho * self.cp * (T - T_REF)
        if self.pcm:
            h = h + self.rho * self.Lf * self.g(T)
        return h

    def c_app(self, T):
        """Apparent volumetric heat capacity dH/dT, J/(m^3 K)."""
        c = np.full_like(np.asarray(T, dtype=float), self.rho * self.cp)
        if self.pcm:
            inside = (T > self.Ts) & (T < self.Tl)
            c = c + inside * self.rho * self.Lf / (self.Tl - self.Ts)
        return c

    def k_of_T(self, T):
        """Conductivity, interpolated across the mushy zone."""
        if not self.pcm or self.k_liquid == self.k:
            return np.full_like(np.asarray(T, dtype=float), self.k)
        return self.k + (self.k_liquid - self.k) * self.g(T)


# --------------------------------------------------------------------------- #
# Mesh + model
# --------------------------------------------------------------------------- #

class Model:
    """Composite 1D domain.

    layers : list of (Material, outer_coordinate, n_cells)
             Inner coordinate of the first layer is `r_inner` (0 by default).
    geom   : 'planar' or 'cylindrical'
    height : cylinder height, m.  For 'planar' this is the cross-sectional area.
             Only scales absolute energies, never temperatures.
    """

    def __init__(self, layers, geom="cylindrical", height=1.0, r_inner=0.0):
        if geom not in ("planar", "cylindrical"):
            raise ValueError("geom must be 'planar' or 'cylindrical'")
        self.geom = geom
        self.m = 1 if geom == "cylindrical" else 0
        self.height = height

        edges = [r_inner]
        mat_idx = []
        self.materials = []
        self.layer_slices = {}
        for mat, r_out, n in layers:
            i0 = len(mat_idx)
            self.materials.append(mat)
            mi = len(self.materials) - 1
            sub = np.linspace(edges[-1], r_out, n + 1)[1:]
            edges.extend(sub.tolist())
            mat_idx.extend([mi] * n)
            self.layer_slices[mat.name] = slice(i0, i0 + n)

        self.re = np.asarray(edges, dtype=float)          # N+1 cell edges
        self.mat_idx = np.asarray(mat_idx, dtype=int)     # N
        self.N = len(self.mat_idx)
        self.rc = 0.5 * (self.re[:-1] + self.re[1:])      # N cell centres
        self.dr = np.diff(self.re)

        if self.m == 1:
            self.V = np.pi * (self.re[1:] ** 2 - self.re[:-1] ** 2) * height
            self.A = 2.0 * np.pi * self.re * height       # N+1 face areas
        else:
            self.V = self.dr * height
            self.A = np.full(self.N + 1, height)

        # per-cell constant properties, gathered once
        self.rho = np.array([self.materials[i].rho for i in self.mat_idx])
        self.cp = np.array([self.materials[i].cp for i in self.mat_idx])
        self.is_pcm = np.array([self.materials[i].pcm for i in self.mat_idx])

        # distance from each cell centre to its left / right face
        self.d_left = self.rc - self.re[:-1]
        self.d_right = self.re[1:] - self.rc

    # ---- per-cell material functions, vectorised over the whole field ---- #

    def _per_cell(self, T, fn):
        out = np.empty(self.N)
        for mi, mat in enumerate(self.materials):
            sel = self.mat_idx == mi
            if sel.any():
                out[sel] = getattr(mat, fn)(T[sel])
        return out

    def H(self, T):
        return self._per_cell(T, "H")

    def c_app(self, T):
        return self._per_cell(T, "c_app")

    def k_of_T(self, T):
        return self._per_cell(T, "k_of_T")

    def g(self, T):
        return self._per_cell(T, "g")

    def energy(self, T):
        """Total enthalpy of the domain, J."""
        return float(np.sum(self.V * self.H(T)))

    def pcm_mass(self):
        sel = self.is_pcm
        return float(np.sum(self.V[sel] * self.rho[sel]))

    def melted_fraction(self, T):
        """Volume-weighted mean liquid fraction over all PCM cells."""
        sel = self.is_pcm
        if not sel.any():
            return 0.0
        return float(np.sum(self.V[sel] * self.g(T)[sel]) / np.sum(self.V[sel]))

    def latent_stored(self, T):
        """Latent heat currently held in the PCM, J."""
        sel = self.is_pcm
        Lf = np.array([self.materials[i].Lf for i in self.mat_idx])
        return float(np.sum(self.V[sel] * self.rho[sel] * Lf[sel] * self.g(T)[sel]))

    def front_position(self, T):
        """Melt-front coordinate, defined by the integral of liquid fraction.

        Energy-consistent (equals the sharp-interface position when the mushy
        zone is thin) and smooth, which matters because Stage B differentiates
        through it.
        """
        sel = self.is_pcm
        if not sel.any():
            return np.nan
        gg = self.g(T)[sel]
        if self.m == 1:
            # find r such that the melted annulus volume matches
            r0 = self.re[:-1][sel][0]
            Vmelt = np.sum(self.V[sel] * gg) / self.height
            return float(np.sqrt(r0 ** 2 + Vmelt / np.pi))
        x0 = self.re[:-1][sel][0]
        return float(x0 + np.sum(self.dr[sel] * gg))

    # ------------------------------------------------------------------ #
    # Face conductances
    # ------------------------------------------------------------------ #

    def _conductances(self, T):
        """Interior face conductances G (W/K), length N-1.

        Harmonic mean weighted by the two half-cell distances: exact for a
        series resistance across a material jump, unlike an arithmetic mean.
        """
        k = self.k_of_T(T)
        dR = self.d_right[:-1] / k[:-1]     # resistance, centre -> face
        dL = self.d_left[1:] / k[1:]        # resistance, face -> next centre
        return self.A[1:-1] / (dR + dL)

    def _outer_conductance(self, T, bc):
        k = self.k_of_T(T)
        if bc[0] == "robin":
            h = bc[1]
            return self.A[-1] / (self.d_right[-1] / k[-1] + 1.0 / h)
        if bc[0] == "dirichlet":
            return self.A[-1] * k[-1] / self.d_right[-1]
        return 0.0

    def _inner_conductance(self, T, bc):
        k = self.k_of_T(T)
        if bc[0] == "robin":
            h = bc[1]
            return self.A[0] / (self.d_left[0] / k[0] + 1.0 / h)
        if bc[0] == "dirichlet":
            return self.A[0] * k[0] / self.d_left[0]
        return 0.0

    # ------------------------------------------------------------------ #
    # Time integration
    # ------------------------------------------------------------------ #

    def solve(self, T_init, t_end, dt, q_vol=0.0, q_layer=None,
              bc_inner=("adiabatic",), bc_outer=("robin", 10.0, 298.15),
              scheme="enthalpy_newton", store_every=1, tol=1e-9, max_iter=200,
              t_start=0.0):
        """Backward-Euler march.  Returns a dict of histories and diagnostics.

        q_vol   scalar W/m^3, or callable q(t) -> W/m^3.
        q_layer name of the layer the source sits in (default: first layer).
        scheme  'ahc_chord' | 'ahc_pointwise' | 'enthalpy_newton'

        Backward Euler is used throughout: it is unconditionally stable, which
        matters because the apparent heat capacity varies by ~50x between mushy
        and non-mushy cells and an explicit step would be limited by the
        smallest.  It is first-order in time, so dt convergence is checked in
        Stage A rather than assumed.
        """
        T = np.array(T_init, dtype=float).copy()
        if T.shape != (self.N,):
            raise ValueError(f"T_init must have shape ({self.N},)")

        sl = self.layer_slices[q_layer] if q_layer else slice(0, 0)
        if q_layer is None:
            sl = list(self.layer_slices.values())[0]
        q_mask = np.zeros(self.N)
        q_mask[sl] = 1.0

        nsteps = int(round(t_end / dt))
        ts, Ts, fs, fronts = [t_start], [T.copy()], [self.melted_fraction(T)], \
                             [self.front_position(T)]

        E0 = self.energy(T)
        E_in = 0.0
        E_out = 0.0
        iters_max = 0
        converged_all = True

        for n in range(1, nsteps + 1):
            t_new = t_start + n * dt
            qv = q_vol(t_new) if callable(q_vol) else q_vol
            src = qv * q_mask * self.V                      # W per cell

            Told = T.copy()
            Hold = self.H(Told)

            T, nit, ok = self._step(Told, Hold, dt, src, bc_inner, bc_outer,
                                    scheme, tol, max_iter)
            iters_max = max(iters_max, nit)
            converged_all = converged_all and ok

            E_in += float(np.sum(src)) * dt
            E_out += self._boundary_loss(T, bc_inner, bc_outer) * dt

            if n % store_every == 0 or n == nsteps:
                ts.append(t_new)
                Ts.append(T.copy())
                fs.append(self.melted_fraction(T))
                fronts.append(self.front_position(T))

        E1 = self.energy(T)
        denom = max(abs(E_in), abs(E1 - E0), 1e-30)
        return {
            "t": np.array(ts),
            "T": np.array(Ts),                   # (nstore, N)
            "T_surf": np.array([x[-1] for x in Ts]),
            "T_core": np.array([x[0] for x in Ts]),
            "melted_fraction": np.array(fs),
            "front": np.array(fronts),
            "T_final": T,
            "E_in": E_in, "E_out": E_out,
            "E_stored": E1 - E0,
            "energy_residual": (E_in - E_out - (E1 - E0)),
            "energy_closure_rel": abs(E_in - E_out - (E1 - E0)) / denom,
            "max_picard_iters": iters_max,
            "converged": converged_all,
            "scheme": scheme, "dt": dt,
        }

    def _boundary_loss(self, T, bc_inner, bc_outer):
        """Net heat leaving the domain through both boundaries, W."""
        loss = 0.0
        Go = self._outer_conductance(T, bc_outer)
        if bc_outer[0] == "robin":
            loss += Go * (T[-1] - bc_outer[2])
        elif bc_outer[0] == "dirichlet":
            loss += Go * (T[-1] - bc_outer[1])
        Gi = self._inner_conductance(T, bc_inner)
        if bc_inner[0] == "robin":
            loss += Gi * (T[0] - bc_inner[2])
        elif bc_inner[0] == "dirichlet":
            loss += Gi * (T[0] - bc_inner[1])
        return loss

    def _assemble(self, C, T, dt, src, bc_inner, bc_outer, Told, extra_rhs=None):
        """Tridiagonal system for one Picard/Newton iterate.

        C is the capacitance (J/m^3/K) used on the diagonal; which C is passed
        is exactly what distinguishes the three schemes.
        """
        N = self.N
        G = self._conductances(T)
        ab = np.zeros((3, N))
        rhs = np.zeros(N)

        cap = self.V * C / dt
        diag = cap.copy()
        rhs += cap * Told + src

        diag[:-1] += G
        diag[1:] += G
        ab[0, 1:] = -G          # upper
        ab[2, :-1] = -G         # lower

        Go = self._outer_conductance(T, bc_outer)
        if bc_outer[0] == "robin":
            diag[-1] += Go
            rhs[-1] += Go * bc_outer[2]
        elif bc_outer[0] == "dirichlet":
            diag[-1] += Go
            rhs[-1] += Go * bc_outer[1]

        Gi = self._inner_conductance(T, bc_inner)
        if bc_inner[0] == "robin":
            diag[0] += Gi
            rhs[0] += Gi * bc_inner[2]
        elif bc_inner[0] == "dirichlet":
            diag[0] += Gi
            rhs[0] += Gi * bc_inner[1]

        if extra_rhs is not None:
            rhs += extra_rhs
        ab[1, :] = diag
        return ab, rhs

    def _step(self, Told, Hold, dt, src, bc_inner, bc_outer, scheme, tol, max_iter):
        T = Told.copy()

        if scheme == "enthalpy_newton":
            for it in range(1, max_iter + 1):
                G = self._conductances(T)
                # residual R(T) = V(H(T)-Hold)/dt - conduction - source
                R = self.V * (self.H(T) - Hold) / dt - src
                flux = np.zeros(self.N)
                flux[:-1] += G * (T[1:] - T[:-1])
                flux[1:] += G * (T[:-1] - T[1:])
                Go = self._outer_conductance(T, bc_outer)
                if bc_outer[0] == "robin":
                    flux[-1] += Go * (bc_outer[2] - T[-1])
                elif bc_outer[0] == "dirichlet":
                    flux[-1] += Go * (bc_outer[1] - T[-1])
                Gi = self._inner_conductance(T, bc_inner)
                if bc_inner[0] == "robin":
                    flux[0] += Gi * (bc_inner[2] - T[0])
                elif bc_inner[0] == "dirichlet":
                    flux[0] += Gi * (bc_inner[1] - T[0])
                R -= flux

                # Jacobian: tridiagonal, conductivity dependence lagged
                ab = np.zeros((3, self.N))
                diag = self.V * self.c_app(T) / dt
                diag[:-1] += G
                diag[1:] += G
                ab[0, 1:] = -G
                ab[2, :-1] = -G
                if bc_outer[0] in ("robin", "dirichlet"):
                    diag[-1] += Go
                if bc_inner[0] in ("robin", "dirichlet"):
                    diag[0] += Gi
                ab[1, :] = diag

                dT = solve_banded((1, 1), ab, -R)
                # c_app is discontinuous at the mushy-zone edges, so a cell
                # sitting exactly on one can two-cycle forever at an amplitude
                # of a few nK.  Damping past 30 iterations breaks the cycle.
                T = T + (dT if it <= 30 else 0.5 * dT)
                if np.max(np.abs(dT)) < tol:
                    return T, it, True
            return T, max_iter, False

        # --- Picard schemes ------------------------------------------------ #
        for it in range(1, max_iter + 1):
            if scheme == "ahc_pointwise":
                C = self.c_app(T)
            elif scheme == "ahc_chord":
                dT = T - Told
                C = self.c_app(Told).copy()
                big = np.abs(dT) > 1e-10
                if big.any():
                    C[big] = (self.H(T)[big] - Hold[big]) / dT[big]
            else:
                raise ValueError(f"unknown scheme {scheme!r}")

            ab, rhs = self._assemble(C, T, dt, src, bc_inner, bc_outer, Told)
            Tnew = solve_banded((1, 1), ab, rhs)
            err = np.max(np.abs(Tnew - T))
            T = Tnew
            if err < tol:
                return T, it, True
        # Not converged.  For 'ahc_chord' this is usually not oscillation but
        # slow LINEAR convergence: a cell that has just crossed the liquidus has
        # an enthalpy chord ~100x the sensible capacity, which makes the Picard
        # map contract at a ratio near 0.94, so ~300 iterations are needed for
        # machine tolerance.  Damping makes it worse, not better.  Raise
        # max_iter, or use 'enthalpy_newton', which solves the same discrete
        # system quadratically.
        return T, max_iter, False


# --------------------------------------------------------------------------- #
# Analytical benchmark
# --------------------------------------------------------------------------- #

def stefan_lambda(Ste):
    """Root of  lam * exp(lam^2) * erf(lam) = Ste / sqrt(pi).

    One-phase Stefan problem (Neumann solution): semi-infinite solid initially
    exactly at T_m, wall at x=0 raised to T_w > T_m at t=0.
    """
    from scipy.optimize import brentq
    from scipy.special import erf
    f = lambda lam: lam * np.exp(lam ** 2) * erf(lam) - Ste / np.sqrt(np.pi)
    return brentq(f, 1e-10, 10.0, xtol=1e-14, rtol=1e-15)


def stefan_front(t, Ste, alpha_l):
    """Analytical melt-front position s(t) = 2 lam sqrt(alpha_l t)."""
    return 2.0 * stefan_lambda(Ste) * np.sqrt(alpha_l * np.asarray(t, dtype=float))


def stefan_profile(x, t, Tw, Tm, Ste, alpha_l):
    """Analytical temperature profile of the one-phase Stefan problem."""
    from scipy.special import erf
    lam = stefan_lambda(Ste)
    x = np.asarray(x, dtype=float)
    s = 2.0 * lam * np.sqrt(alpha_l * t)
    T = np.where(x < s,
                 Tw - (Tw - Tm) * erf(x / (2.0 * np.sqrt(alpha_l * t))) / erf(lam),
                 Tm)
    return T
