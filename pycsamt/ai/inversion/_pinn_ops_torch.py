# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
PyTorch implementations of the PINN optimisation
kernels used by pinn1d / pinn2d / pinn3d and the
corresponding Hybrid inverters.

All ``import torch`` calls are deferred to function
bodies so that importing this module never fails on
systems without PyTorch.
"""

from __future__ import annotations

import math

import numpy as np

__all__ = [
    "_mt1d_torch",
    "_mt1d_torch_batch",
    "_fit_station_torch",
    "_fit_2d_joint_torch",
    "_fit_3d_joint_torch",
]

_MU0 = 4e-7 * math.pi
_PI = math.pi


# ── single-station Wait recursion ──────────────────


def _mt1d_torch(freqs_t, log_rho, log_thick):
    r"""
    Wait (1954) 1-D MT recursion in PyTorch.

    Parameters
    ----------
    freqs_t : Tensor (F,) float64
        Frequencies in Hz.
    log_rho : Tensor (L,) float64
        Base-10 log-resistivity per layer.
    log_thick : Tensor (L-1,) float64
        Base-10 log-thickness per layer in m.

    Returns
    -------
    rho_a : Tensor (F,) float64
        Apparent resistivity in Ohm.m.
    phase : Tensor (F,) float64
        Impedance phase in degrees.
    """
    import torch

    rho = 10.0**log_rho
    thick = 10.0**log_thick
    omega = 2.0 * _PI * freqs_t

    j_mu_w = torch.complex(
        torch.zeros_like(omega),
        _MU0 * omega,
    )
    rho_c = rho.to(torch.cdouble)
    k_n = torch.sqrt(j_mu_w / rho_c[-1])
    Z = j_mu_w / k_n

    L = rho.shape[0]
    j_unit = torch.tensor(
        0 + 1j,
        dtype=torch.cdouble,
        device=freqs_t.device,
    )
    for idx in range(L - 2, -1, -1):
        k_j = torch.sqrt(j_mu_w / rho_c[idx])
        Z0_j = j_mu_w / k_j
        h = thick[idx].to(torch.cdouble)
        kh = k_j * h
        tanh_v = torch.tanh(j_unit * kh)
        denom = Z0_j + Z * tanh_v
        Z = Z0_j * (Z + Z0_j * tanh_v) / denom

    rho_a = torch.abs(Z) ** 2 / (_MU0 * omega)
    phase = torch.angle(Z) * (180.0 / _PI)
    return rho_a, phase


# ── batched Wait recursion (S stations) ────────────


def _mt1d_torch_batch(freqs_t, log_rho, log_thick):
    r"""
    Batched Wait (1954) MT-1D forward over S stations.

    Parameters
    ----------
    freqs_t : Tensor (F,) float64
        Common frequency grid in Hz.
    log_rho : Tensor (S, L) float64
        Base-10 log-resistivity per station/layer.
    log_thick : Tensor (S, L-1) float64
        Base-10 log-thickness (m) per station.

    Returns
    -------
    rho_a : Tensor (S, F) float64
    phase : Tensor (S, F) float64
    """
    import torch

    S, L = log_rho.shape
    rho = 10.0**log_rho
    thick = 10.0**log_thick
    omega = 2.0 * _PI * freqs_t

    j_mu_w = torch.complex(torch.zeros_like(omega), _MU0 * omega,).unsqueeze(
        0
    )  # (1, F)

    rho_c = rho.to(torch.cdouble)

    k_n = torch.sqrt(j_mu_w / rho_c[:, -1:])
    Z = j_mu_w / k_n  # (S, F)

    j_unit = torch.tensor(
        0 + 1j,
        dtype=torch.cdouble,
        device=freqs_t.device,
    )
    for idx in range(L - 2, -1, -1):
        k_j = torch.sqrt(j_mu_w / rho_c[:, idx : idx + 1])
        Z0_j = j_mu_w / k_j
        h = thick[:, idx].to(torch.cdouble).unsqueeze(1)
        kh = k_j * h
        tanh_v = torch.tanh(j_unit * kh)
        denom = Z0_j + Z * tanh_v
        Z = Z0_j * (Z + Z0_j * tanh_v) / denom

    rho_a = torch.abs(Z) ** 2 / (_MU0 * omega.unsqueeze(0))
    phase = torch.angle(Z) * (180.0 / _PI)
    return rho_a, phase


# ── 1-D per-station optimisation ───────────────────


def _fit_station_torch(
    obs,
    n_layers: int,
    depth_max: float,
    lam: float,
    lr: float,
    epochs: int,
    device: str,
    log_every: int,
    init_log_rho: np.ndarray | None = None,
    init_log_thick: np.ndarray | None = None,
) -> dict:
    """
    Optimise log-rho and log-thick for one station.

    Parameters
    ----------
    obs : SiteObs1D
        Observed data container.
    n_layers : int
        Number of earth layers.
    depth_max : float
        Approximate investigation depth in m.
    lam : float
        Smoothness regularisation weight.
    lr : float
        Adam learning rate.
    epochs : int
        Gradient-descent iterations.
    device : str
        PyTorch device string.
    log_every : int
        Print-loss interval (0 = silent).
    init_log_rho : ndarray (n_layers,) or None
    init_log_thick : ndarray (n_layers-1,) or None

    Returns
    -------
    dict : log_rho, log_thick, history
    """
    import torch

    freq = obs.freq
    rho_obs = obs.rho_obs
    ph_obs = obs.phase_obs

    freqs_t = torch.tensor(
        freq,
        dtype=torch.float64,
        device=device,
    )
    log_rho_obs_t = torch.tensor(
        np.log10(np.maximum(rho_obs, 1e-12)),
        dtype=torch.float64,
        device=device,
    )
    ph_obs_t = torch.tensor(
        ph_obs,
        dtype=torch.float64,
        device=device,
    )

    if init_log_rho is None:
        init_lr = float(np.mean(np.log10(np.maximum(rho_obs, 1e-3))))
        init_lr_v = np.full(n_layers, init_lr, dtype=np.float64)
    else:
        init_lr_v = np.asarray(init_log_rho, dtype=np.float64)

    if init_log_thick is None:
        dz = depth_max / n_layers
        init_lt_v = np.full(
            n_layers - 1,
            np.log10(max(dz, 1.0)),
            dtype=np.float64,
        )
    else:
        init_lt_v = np.asarray(init_log_thick, dtype=np.float64)

    log_rho = torch.tensor(
        init_lr_v,
        dtype=torch.float64,
        device=device,
        requires_grad=True,
    )
    log_thick = torch.tensor(
        init_lt_v,
        dtype=torch.float64,
        device=device,
        requires_grad=True,
    )

    opt = torch.optim.Adam([log_rho, log_thick], lr=lr)
    history: list[float] = []

    valid = torch.isfinite(log_rho_obs_t) & torch.isfinite(ph_obs_t)

    for ep in range(1, epochs + 1):
        opt.zero_grad()
        rho_pred, ph_pred = _mt1d_torch(freqs_t, log_rho, log_thick)
        lr_pred = torch.log10(rho_pred.clamp(min=1e-12))

        # NaN-safe loss: detach pred at missing obs
        lr_safe = torch.where(valid, log_rho_obs_t, lr_pred.detach())
        ph_safe = torch.where(valid, ph_obs_t, ph_pred.detach())
        n_v = valid.sum().clamp(min=1)
        data_loss = (
            (lr_pred - lr_safe) ** 2 + ((ph_pred - ph_safe) / 90.0) ** 2
        ).sum() / n_v
        smooth = ((log_rho[1:] - log_rho[:-1]) ** 2).mean()
        loss = data_loss + lam * smooth
        loss.backward()
        opt.step()
        history.append(loss.detach().item())

        if log_every > 0 and (ep == 1 or ep % log_every == 0):
            print(f"  {obs.name}  ep {ep:>5d}/{epochs}  loss={loss.item():.5f}")

    return {
        "log_rho": (log_rho.detach().cpu().numpy()),
        "log_thick": (log_thick.detach().cpu().numpy()),
        "history": history,
    }


# ── 2-D joint optimisation ─────────────────────────


def _fit_2d_joint_torch(
    log_rho_obs,
    ph_obs,
    freqs_t,
    n_layers,
    depth_max,
    lam_z,
    lam_x,
    lr,
    epochs,
    device,
    log_every,
    init_log_rho=None,
    init_log_thick=None,
):
    r"""
    Jointly optimise a pseudo-2D section via Adam.

    Parameters
    ----------
    log_rho_obs : ndarray (S, F)
        Observed log10(rho_a), NaN for missing.
    ph_obs : ndarray (S, F)
        Observed phase in degrees, NaN for missing.
    freqs_t : Tensor (F,) float64
        Common frequency grid.
    n_layers : int
    depth_max : float
    lam_z : float
        Vertical smoothness weight.
    lam_x : float
        Lateral smoothness weight.
    lr : float
    epochs : int
    device : str
    log_every : int
    init_log_rho : ndarray (S, L) or None
    init_log_thick : ndarray (S, L-1) or None

    Returns
    -------
    dict : log_rho (S,L), log_thick (S,L-1), history
    """
    import torch
    import torch.nn as nn

    # Accept numpy or tensor for freqs_t
    freqs_t = torch.as_tensor(
        np.asarray(freqs_t),
        dtype=torch.float64,
        device=device,
    )

    S, F = log_rho_obs.shape
    L = n_layers

    lr_obs_t = torch.tensor(
        log_rho_obs,
        dtype=torch.float64,
        device=device,
    )
    ph_obs_t = torch.tensor(
        ph_obs,
        dtype=torch.float64,
        device=device,
    )
    valid = torch.isfinite(lr_obs_t) & torch.isfinite(ph_obs_t)

    if init_log_rho is None:
        mean_lr = float(np.nanmean(log_rho_obs))
        if not np.isfinite(mean_lr):
            mean_lr = 2.0
        init_lr_v = np.full((S, L), mean_lr)
    else:
        init_lr_v = np.asarray(init_log_rho, dtype=float)

    if init_log_thick is None:
        dz = max(depth_max / L, 1.0)
        init_lt_v = np.full((S, L - 1), np.log10(dz))
    else:
        init_lt_v = np.asarray(init_log_thick, dtype=float)

    log_rho = torch.tensor(
        init_lr_v,
        dtype=torch.float64,
        device=device,
        requires_grad=True,
    )
    log_thick = torch.tensor(
        init_lt_v,
        dtype=torch.float64,
        device=device,
        requires_grad=True,
    )

    opt = torch.optim.Adam([log_rho, log_thick], lr=lr)
    history: list[float] = []

    for ep in range(1, epochs + 1):
        opt.zero_grad()
        rho_pred, ph_pred = _mt1d_torch_batch(freqs_t, log_rho, log_thick)
        lr_pred = torch.log10(rho_pred.clamp(min=1e-12))

        lr_safe = torch.where(valid, lr_obs_t, lr_pred.detach())
        ph_safe = torch.where(valid, ph_obs_t, ph_pred.detach())
        res_rho = (lr_pred - lr_safe) ** 2
        res_ph = ((ph_pred - ph_safe) / 90.0) ** 2
        n_v = valid.sum().clamp(min=1)
        data_loss = (res_rho + res_ph).sum() / n_v

        vert = ((log_rho[:, 1:] - log_rho[:, :-1]) ** 2).mean()

        if S > 1:
            lat = ((log_rho[1:, :] - log_rho[:-1, :]) ** 2).mean()
        else:
            lat = torch.tensor(
                0.0,
                dtype=torch.float64,
                device=device,
            )

        loss = data_loss + lam_z * vert + lam_x * lat
        loss.backward()
        nn.utils.clip_grad_norm_([log_rho, log_thick], max_norm=5.0)
        opt.step()
        with torch.no_grad():
            log_thick.clamp_(0.0, 5.0)
        history.append(loss.detach().item())

        if log_every > 0 and ep % log_every == 0:
            print(f"  epoch {ep:>5d}/{epochs}  loss={history[-1]:.5f}")

    return {
        "log_rho": (log_rho.detach().cpu().numpy()),
        "log_thick": (log_thick.detach().cpu().numpy()),
        "history": history,
    }


# ── 3-D joint optimisation ─────────────────────────


def _fit_3d_joint_torch(
    log_rho_obs,
    ph_obs,
    freqs_t,
    adjacency,
    n_layers,
    depth_max,
    lam_z,
    lam_g,
    lr,
    epochs,
    device,
    log_every,
    init_log_rho=None,
    init_log_thick=None,
):
    r"""
    Jointly optimise a quasi-3D section via Adam.

    Uses a graph-Laplacian spatial smoothness term
    to couple adjacent stations.

    Parameters
    ----------
    log_rho_obs : ndarray (S, F)
    ph_obs : ndarray (S, F)
    freqs_t : Tensor (F,) float64
    adjacency : ndarray (S, S)
        Station adjacency matrix.
    n_layers : int
    depth_max : float
    lam_z : float
        Vertical smoothness weight.
    lam_g : float
        Graph spatial smoothness weight.
    lr : float
    epochs : int
    device : str
    log_every : int
    init_log_rho : ndarray (S, L) or None
    init_log_thick : ndarray (S, L-1) or None

    Returns
    -------
    dict : log_rho (S,L), log_thick (S,L-1), history
    """
    import torch
    import torch.nn as nn

    # Accept numpy or tensor for freqs_t
    freqs_t = torch.as_tensor(
        np.asarray(freqs_t),
        dtype=torch.float64,
        device=device,
    )

    S, F = log_rho_obs.shape
    L = n_layers

    lr_obs_t = torch.tensor(
        log_rho_obs,
        dtype=torch.float64,
        device=device,
    )
    ph_obs_t = torch.tensor(
        ph_obs,
        dtype=torch.float64,
        device=device,
    )
    valid = torch.isfinite(lr_obs_t) & torch.isfinite(ph_obs_t)

    A_t = torch.tensor(
        adjacency,
        dtype=torch.float64,
        device=device,
    )
    d_vec = A_t.sum(dim=1)

    if init_log_rho is None:
        mean_lr = float(np.nanmean(log_rho_obs))
        if not np.isfinite(mean_lr):
            mean_lr = 2.0
        init_lr_v = np.full((S, L), mean_lr)
    else:
        init_lr_v = np.asarray(init_log_rho, dtype=float)

    if init_log_thick is None:
        dz = max(depth_max / L, 1.0)
        init_lt_v = np.full((S, L - 1), np.log10(dz))
    else:
        init_lt_v = np.asarray(init_log_thick, dtype=float)

    log_rho = torch.tensor(
        init_lr_v,
        dtype=torch.float64,
        device=device,
        requires_grad=True,
    )
    log_thick = torch.tensor(
        init_lt_v,
        dtype=torch.float64,
        device=device,
        requires_grad=True,
    )

    opt = torch.optim.Adam([log_rho, log_thick], lr=lr)
    history: list[float] = []

    for ep in range(1, epochs + 1):
        opt.zero_grad()
        rho_pred, ph_pred = _mt1d_torch_batch(freqs_t, log_rho, log_thick)
        lr_pred = torch.log10(rho_pred.clamp(min=1e-12))

        lr_safe = torch.where(valid, lr_obs_t, lr_pred.detach())
        ph_safe = torch.where(valid, ph_obs_t, ph_pred.detach())
        res_rho = (lr_pred - lr_safe) ** 2
        res_ph = ((ph_pred - ph_safe) / 90.0) ** 2
        n_v = valid.sum().clamp(min=1)
        data_loss = (res_rho + res_ph).sum() / n_v

        vert = ((log_rho[:, 1:] - log_rho[:, :-1]) ** 2).mean()

        # Graph Laplacian: L = D - A
        Lm = d_vec.unsqueeze(1) * log_rho - A_t @ log_rho
        spatial = (Lm * log_rho).sum() / log_rho.numel()

        loss = data_loss + lam_z * vert + lam_g * spatial
        loss.backward()
        nn.utils.clip_grad_norm_([log_rho, log_thick], max_norm=5.0)
        opt.step()
        with torch.no_grad():
            log_thick.clamp_(0.0, 5.0)
        history.append(loss.detach().item())

        if log_every > 0 and ep % log_every == 0:
            print(f"  epoch {ep:>5d}/{epochs}  loss={history[-1]:.5f}")

    return {
        "log_rho": (log_rho.detach().cpu().numpy()),
        "log_thick": (log_thick.detach().cpu().numpy()),
        "history": history,
    }
