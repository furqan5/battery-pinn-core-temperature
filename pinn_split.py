"""Base-subtracted (lumped-split) formulation of the same inverse problem.

WHY THIS EXISTS.  Stage D's first attempt -- a plain full-field PINN -- failed
its gate at 39% core-surface error.  Diagnosis: the drive-cycle source is
genuinely broadband (50% of its energy above harmonic 233, i.e. periods under
15 s; 99% below 2.3 s), while a Fourier-feature MLP with n_freq=128 resolves
only down to 27.7 s.  No amount of training fixes a basis that cannot represent
the forcing.

THE FIX IS EXACT ALGEBRA, NOT AN APPROXIMATION.  Split

    T(r,t) = T_l(t) + w(r,t)

where T_l is the LUMPED (0-D) solution

    rho cp dT_l/dt = q(t) - (2h/R)(T_l - T_inf).

Substituting into the full PDE, the spiky source q CANCELS IDENTICALLY:

    rho cp dw/dt = k lap(w) + (2h/R)(T_l - T_inf)

with  dw/dr|_R = -(h/k)(w(R) + T_l - T_inf)  and  w(r,0) = 0.

The remaining forcing is proportional to T_l - T_inf, which is the output of a
first-order low-pass with tau = rho cp R/(2h) = 405 s.  At the 15 s timescale
that filter attenuates by ~169x, so w is SMOOTH and a small Fourier basis
suffices.

Is the PINN then doing trivial work?  No.  T_l is a 0-D model and contains NO
radial information whatsoever; the entire core-to-surface gradient -- the only
quantity this project cares about -- lives in w.  The split removes the stiff
part we can integrate exactly and leaves the PINN exactly the part we cannot.

The initial offset (trap 5.2) is carried by T_l(0), which is set to the MEASURED
initial temperature.  Because w(r,0) = 0 exactly, the hard IC w = t_hat * N is
now exactly right rather than fighting a boundary layer.
"""

from __future__ import annotations
import numpy as np
import torch
import torch.nn as nn

from pinn_core import FourierTime, torch_interp

torch.set_default_dtype(torch.float64)


# --------------------------------------------------------------------------- #
# Lumped base, precomputed by superposition
# --------------------------------------------------------------------------- #

class LumpedBase:
    """Exact exponential integration of the 0-D model, linear in the source.

    Because the lumped ODE is LINEAR in q, and q = I^2 R_eff (1 + sum_j a_j x^j),
    we can precompute one response per basis term and then form

        T_l(t) = A(t) + R_eff * [ B_0(t) + sum_j a_j B_j(t) ]

    as a cheap differentiable combination.  No ODE solve inside the training
    loop.

      A    : q = 0, actual T_inf(t), initial condition = measured T0
      B_j  : q = I^2 x^j / V_b, T_inf = 0, zero initial condition
    """

    def __init__(self, rec, h, rho_cp, R, V_b, T0, n_shape=0):
        self.tau = rho_cp * R / (2.0 * h)
        dt = rec.dt
        self.dt = dt
        n = len(rec.t)
        e = np.exp(-dt / self.tau)
        g = self.tau / rho_cp * (1.0 - e)      # gain on q over one step

        # homogeneous + ambient response
        A = np.empty(n)
        A[0] = T0
        for i in range(1, n):
            A[i] = rec.T_inf[i] + (A[i - 1] - rec.T_inf[i]) * e
        self.A = A

        # source basis responses
        xn = self._xn(rec)
        self.B = []
        for j in range(n_shape + 1):
            qj = (rec.I ** 2) * (xn ** j) / V_b
            B = np.zeros(n)
            for i in range(1, n):
                B[i] = B[i - 1] * e + g * qj[i]
            self.B.append(B)

        self.n_shape = n_shape
        self.t_hat = torch.tensor(rec.t / rec.t[-1])
        self.A_t = torch.tensor(A)
        self.B_t = [torch.tensor(b) for b in self.B]
        self.T_inf_t = torch.tensor(rec.T_inf)
        self.xn = xn

    @staticmethod
    def _xn(rec):
        c = float(rec.dod.mean())
        span = float(max(rec.dod.max() - rec.dod.min(), 1e-12))
        return 2.0 * (rec.dod - c) / span            # in [-1, 1]

    def T_l_np(self, R_eff, shape=()):
        out = self.A + R_eff * self.B[0]
        for j, a in enumerate(shape, start=1):
            out = out + R_eff * a * self.B[j]
        return out

    def T_l_torch(self, R_eff, shape=None):
        """Lumped base on the native data grid, differentiable in R_eff/shape.

        No interpolation: the data grid and the base grid are the same, so this
        is exact and cheap -- which is what makes the inverse fit affordable.
        """
        out = self.A_t + R_eff * self.B_t[0]
        if shape is not None and len(shape) > 0:
            for j, a in enumerate(shape, start=1):
                out = out + R_eff * a * self.B_t[j]
        return out

    def excess(self, t_hat, R_eff, shape=None):
        """(T_l - T_inf) at arbitrary t_hat, differentiable in R_eff/shape."""
        A = torch_interp(t_hat, self.t_hat, self.A_t)
        Ti = torch_interp(t_hat, self.t_hat, self.T_inf_t)
        tot = A + R_eff * torch_interp(t_hat, self.t_hat, self.B_t[0])
        if shape is not None and len(shape) > 0:
            for j, a in enumerate(shape, start=1):
                tot = tot + R_eff * a * torch_interp(t_hat, self.t_hat, self.B_t[j])
        return tot - Ti


# --------------------------------------------------------------------------- #
# Network for the deviation w
# --------------------------------------------------------------------------- #

class DeviationNet(nn.Module):
    """omega(s, t_hat) = t_hat * N(s, t_hat).

    w(r,0) = 0 EXACTLY by construction, because the measured initial offset is
    carried entirely by the lumped base T_l(0).  Symmetry at the axis is still
    by construction via s = r_hat^2.

    MEASURED DEFAULT: n_freq = 0.  A sweep on the constant-source control gave
    core-surface errors of 0.086% (n_freq=0), 0.15% (4), 0.49% (16) and 8.0%
    (64) -- monotonically worse with more Fourier features.  They cannot reach
    the bandwidth the raw problem needs (that would take n_freq ~ 1500) and they
    badly degrade the optimisation landscape for the smooth solution we actually
    have after splitting.  So: none.
    """

    def __init__(self, width=48, depth=4, n_freq=0):
        super().__init__()
        self.n_freq = n_freq
        self.ff = FourierTime(n_freq) if n_freq > 0 else None
        d_in = 2 if n_freq == 0 else 1 + (1 + 2 * n_freq)
        layers, d = [], d_in
        for _ in range(depth):
            layers += [nn.Linear(d, width), nn.Tanh()]
            d = width
        layers += [nn.Linear(d, 1)]
        self.net = nn.Sequential(*layers)
        for m in self.net:
            if isinstance(m, nn.Linear):
                nn.init.xavier_normal_(m.weight)
                nn.init.zeros_(m.bias)

    def forward(self, s, t):
        z = torch.cat([s, t], dim=1) if self.ff is None \
            else torch.cat([s, self.ff(t)], dim=1)
        return t * self.net(z)


# --------------------------------------------------------------------------- #
# Residuals for the split problem
# --------------------------------------------------------------------------- #

def split_pde_residual(net, Fo, Bi, base, s, t, R_eff, shape=None, dT_ref=1.0):
    """(1/Fo) om_t - (4 om_s + 4 s om_ss) - 2 Bi * g,  g = (T_l - T_inf)/dT_ref."""
    s = s.requires_grad_(True)
    t = t.requires_grad_(True)
    om = net(s, t)
    om_s = torch.autograd.grad(om, s, torch.ones_like(om), create_graph=True)[0]
    om_ss = torch.autograd.grad(om_s, s, torch.ones_like(om_s), create_graph=True)[0]
    om_t = torch.autograd.grad(om, t, torch.ones_like(om), create_graph=True)[0]
    g = base.excess(t, R_eff, shape) / dT_ref
    return (1.0 / Fo) * om_t - (4.0 * om_s + 4.0 * s * om_ss) - 2.0 * Bi * g


def split_bc_residual(net, Bi, base, t, R_eff, shape=None, dT_ref=1.0):
    """2 om_s + Bi (om + g) = 0 at s = 1."""
    s = torch.ones_like(t).requires_grad_(True)
    t = t.requires_grad_(True)
    om = net(s, t)
    om_s = torch.autograd.grad(om, s, torch.ones_like(om), create_graph=True)[0]
    g = base.excess(t, R_eff, shape) / dT_ref
    return 2.0 * om_s + Bi * (om + g)


def reconstruct(net, base, t_hat_np, R_eff, shape=None, dT_ref=1.0):
    """Dimensional core and surface temperatures T = T_l + dT_ref * omega."""
    t = torch.tensor(t_hat_np).reshape(-1, 1)
    with torch.no_grad():
        oc = net(torch.zeros_like(t), t).squeeze(1).numpy()
        os_ = net(torch.ones_like(t), t).squeeze(1).numpy()
    sh = () if shape is None else tuple(np.atleast_1d(shape))
    Tl = base.T_l_np(R_eff, sh)
    return Tl + dT_ref * oc, Tl + dT_ref * os_
