"""Cylindrical inverse PINN -- network, physics residual, and training loops.

Non-dimensionalisation (mandatory here, per the project's own history: training
in physical units once gave a residual loss 37x smaller while the answer was
120x worse).

    r_hat = r/R          in [0,1]
    t_hat = t/t_ref      in [0,1]
    theta = (T - T_ref)/dT_ref

PDE becomes
    (1/Fo) d(theta)/d(t_hat) = laplacian(theta) + Q_hat
    Fo    = alpha * t_ref / R^2
    Q_hat = q''' R^2 / (k dT_ref)

SYMMETRY AT THE AXIS IS HANDLED BY CONSTRUCTION, not by a penalty.  The network
takes s = r_hat^2 as its spatial input, so

    d(theta)/d(r_hat) = 2 r_hat * theta_s   ->  vanishes at r_hat = 0 identically,

and the radial Laplacian written in s carries NO 1/r term at all:

    (1/r) d/dr ( r dT/dr )  ==  4 theta_s + 4 s theta_ss

so the r=0 singularity is removed analytically rather than masked.
"""

from __future__ import annotations
import math
import numpy as np
import torch
import torch.nn as nn

torch.set_default_dtype(torch.float64)


# --------------------------------------------------------------------------- #
# Scaling bundle
# --------------------------------------------------------------------------- #

class Scales:
    def __init__(self, rec, k, rho_cp, R):
        self.R = R
        self.t_ref = float(rec.t[-1])
        self.T_ref = float(rec.T_inf.mean())
        self.dT_ref = float(rec.T_s.max() - self.T_ref)
        self.alpha = k / rho_cp
        self.Fo = self.alpha * self.t_ref / R ** 2
        self.k = k
        self.rho_cp = rho_cp

    def nd_T(self, T):
        return (T - self.T_ref) / self.dT_ref

    def dim_T(self, th):
        return self.T_ref + th * self.dT_ref

    def q_hat(self, q_vol):
        """W/m^3 -> dimensionless source."""
        return q_vol * self.R ** 2 / (self.k * self.dT_ref)

    def __repr__(self):
        return (f"Scales(t_ref={self.t_ref:.1f}s, T_ref={self.T_ref:.3f}C, "
                f"dT_ref={self.dT_ref:.3f}K, Fo={self.Fo:.4f}, "
                f"alpha={self.alpha:.4e} m2/s)")


# --------------------------------------------------------------------------- #
# Network
# --------------------------------------------------------------------------- #

class FourierTime(nn.Module):
    """Fixed Fourier features in t_hat.

    Justified by measurement, not habit: 99.9% of the surface-signal energy sits
    below 0.024 Hz, i.e. harmonic ~86 of this record, because the cell itself is
    a low-pass filter with a 405 s convective time constant.  n_freq=128 covers
    that with margin while keeping the input width modest.
    """

    def __init__(self, n_freq=128):
        super().__init__()
        self.register_buffer("kf", torch.arange(1, n_freq + 1, dtype=torch.float64))

    def forward(self, t):
        a = 2.0 * math.pi * t * self.kf
        return torch.cat([t, torch.sin(a), torch.cos(a)], dim=1)


class FieldNet(nn.Module):
    """theta(s, t_hat) with the initial condition imposed HARD.

        theta = theta0 + t_hat * N(s, t_hat)

    theta0 is the MEASURED initial offset, not zero.  Forcing theta(.,0)=0 would
    assert the cell starts at ambient; here it starts 0.43 K above it (DS1).
    That offset is printed at construction so it can never be silently lost.
    """

    def __init__(self, theta0, width=64, depth=4, n_freq=128):
        super().__init__()
        self.theta0 = float(theta0)
        self.ff = FourierTime(n_freq)
        d_in = 1 + (1 + 2 * n_freq)
        layers, d = [], d_in
        for _ in range(depth):
            layers += [nn.Linear(d, width), nn.Tanh()]
            d = width
        layers += [nn.Linear(d, 1)]
        self.net = nn.Sequential(*layers)
        for m in self.net:
            if isinstance(m, nn.Linear):
                nn.init.xavier_normal_(m.weight, gain=1.0)
                nn.init.zeros_(m.bias)

    def forward(self, s, t):
        z = torch.cat([s, self.ff(t)], dim=1)
        return self.theta0 + t * self.net(z)


# --------------------------------------------------------------------------- #
# Source term
# --------------------------------------------------------------------------- #

def torch_interp(x, xp, fp):
    """Linear interpolation, differentiable in fp but not needed in x."""
    i = torch.clamp(torch.searchsorted(xp, x.contiguous()), 1, len(xp) - 1)
    x0, x1 = xp[i - 1], xp[i]
    y0, y1 = fp[i - 1], fp[i]
    w = (x - x0) / (x1 - x0)
    return y0 + w * (y1 - y0)


class Source:
    """q'''(t) = I(t)^2 * R_eff(dod) / V_cell, in W/m^3.

    TRAP 5.1 LIVES HERE.  The I(t)^2 factor is the whole point: this drive cycle
    contains genuine zero-current rests, and a source written without the current
    factor keeps injecting heat through them.  It still fits the mean acceptably
    and corrupts the recovered coefficient completely.  `use_current=False`
    exists only to demonstrate that failure deliberately.
    """

    def __init__(self, rec, V_b, use_current=True):
        self.t_grid = torch.tensor(rec.t / rec.t[-1])
        self.I2 = torch.tensor(rec.I ** 2)
        self.dod = torch.tensor(rec.dod)
        self.dod_c = float(rec.dod.mean())
        self.dod_span = float(max(rec.dod.max() - rec.dod.min(), 1e-12))
        self.V_b = V_b
        self.use_current = use_current
        self.I2_mean = float((rec.I ** 2).mean())

    def q_vol(self, t_hat, R_eff, shape=None):
        I2 = (torch_interp(t_hat, self.t_grid, self.I2) if self.use_current
              else torch.full_like(t_hat, self.I2_mean))
        R = R_eff
        if shape is not None and len(shape) > 0:
            x = torch_interp(t_hat, self.t_grid, self.dod)
            xn = 2.0 * (x - self.dod_c) / self.dod_span      # in [-1,1]
            f = torch.ones_like(xn)
            for j, a in enumerate(shape, start=1):
                f = f + a * xn ** j
            R = R_eff * f
        return I2 * R / self.V_b


# --------------------------------------------------------------------------- #
# Physics residual
# --------------------------------------------------------------------------- #

def pde_residual(net, sc, src, s, t, R_eff, shape=None):
    """(1/Fo) theta_t - (4 theta_s + 4 s theta_ss) - Q_hat."""
    s = s.requires_grad_(True)
    t = t.requires_grad_(True)
    th = net(s, t)
    th_s = torch.autograd.grad(th, s, torch.ones_like(th), create_graph=True)[0]
    th_ss = torch.autograd.grad(th_s, s, torch.ones_like(th_s), create_graph=True)[0]
    th_t = torch.autograd.grad(th, t, torch.ones_like(th), create_graph=True)[0]
    lap = 4.0 * th_s + 4.0 * s * th_ss
    Qh = sc.q_hat(src.q_vol(t, R_eff, shape))
    return (1.0 / sc.Fo) * th_t - lap - Qh


def bc_residual(net, sc, t, Bi, theta_inf_t):
    """Robin at r_hat = 1:  d(theta)/d(r_hat) + Bi (theta - theta_inf) = 0.

    At s = 1, d(theta)/d(r_hat) = 2 theta_s.
    """
    s = torch.ones_like(t).requires_grad_(True)
    t = t.requires_grad_(True)
    th = net(s, t)
    th_s = torch.autograd.grad(th, s, torch.ones_like(th), create_graph=True)[0]
    th_inf = torch_interp(t, theta_inf_t[0], theta_inf_t[1])
    return 2.0 * th_s + Bi * (th - th_inf)


# --------------------------------------------------------------------------- #
# Observables
# --------------------------------------------------------------------------- #

@torch.no_grad()
def predict_profile(net, sc, t_hat_np, s_np):
    t = torch.tensor(t_hat_np).reshape(-1, 1)
    out = np.empty((len(t_hat_np), len(s_np)))
    for j, sv in enumerate(s_np):
        s = torch.full_like(t, float(sv))
        out[:, j] = net(s, t).squeeze(1).numpy()
    return out


def predict_core_surf(net, sc, t_hat_np):
    """Dimensional core (s=0) and surface (s=1) temperatures."""
    t = torch.tensor(t_hat_np).reshape(-1, 1)
    with torch.no_grad():
        thc = net(torch.zeros_like(t), t).squeeze(1).numpy()
        ths = net(torch.ones_like(t), t).squeeze(1).numpy()
    return sc.dim_T(thc), sc.dim_T(ths)
