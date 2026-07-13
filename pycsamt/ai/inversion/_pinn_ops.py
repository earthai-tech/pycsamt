# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Backend-dispatching PINN optimisation kernels.

Each public function probes the active backend via
:func:`pycsamt.backends.get_backend` and forwards the
call to either :mod:`_pinn_ops_torch` or
:mod:`_pinn_ops_tf`.

User code (pinn1d / pinn2d / pinn3d and the Hybrid
classes) should import from here rather than from the
backend-specific modules.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "fit_station",
    "fit_2d_joint",
    "fit_3d_joint",
]


def _backend() -> str:
    """Return active backend name or raise."""
    from pycsamt.backends import get_backend

    name = get_backend()
    if name == "none":
        raise ImportError(
            "PINNInverter requires PyTorch or "
            "TensorFlow.\n"
            "Install one with:\n"
            "    pip install torch\n"
            "    pip install tensorflow"
        )
    return name


# ── 1-D per-station ────────────────────────────────


def fit_station(
    obs,
    *,
    n_layers: int,
    depth_max: float,
    lam: float,
    lr: float,
    epochs: int,
    device: str,
    log_every: int,
    init_log_rho: np.ndarray | None = None,
    init_log_thick: np.ndarray | None = None,
    backend: str | None = None,
) -> dict:
    """
    Optimise log-rho/log-thick for one station.

    Dispatches to the PyTorch or TensorFlow kernel
    based on the active backend.

    Parameters
    ----------
    obs : SiteObs1D
        Observed data container.
    n_layers : int
    depth_max : float
    lam : float
        Smoothness weight.
    lr : float
        Adam learning rate.
    epochs : int
    device : str
        Device string (``'cpu'``, ``'cuda'``,
        ``'/GPU:0'``, …).
    log_every : int
        Print-loss interval (0 = silent).
    init_log_rho : ndarray (n_layers,) or None
    init_log_thick : ndarray (n_layers-1,) or None
    backend : str or None
        Force ``'torch'`` or ``'tensorflow'``.
        Auto-detected when ``None``.

    Returns
    -------
    dict
        Keys: ``log_rho``, ``log_thick``, ``history``.
    """
    bk = backend or _backend()
    if bk == "torch":
        from ._pinn_ops_torch import (
            _fit_station_torch as _impl,
        )
    else:
        from ._pinn_ops_tf import (
            _fit_station_tf as _impl,
        )
    return _impl(
        obs,
        n_layers=n_layers,
        depth_max=depth_max,
        lam=lam,
        lr=lr,
        epochs=epochs,
        device=device,
        log_every=log_every,
        init_log_rho=init_log_rho,
        init_log_thick=init_log_thick,
    )


# ── 2-D joint ──────────────────────────────────────


def fit_2d_joint(
    log_rho_obs: np.ndarray,
    ph_obs: np.ndarray,
    freqs_t,
    *,
    n_layers: int,
    depth_max: float,
    lam_z: float,
    lam_x: float,
    lr: float,
    epochs: int,
    device: str,
    log_every: int,
    init_log_rho: np.ndarray | None = None,
    init_log_thick: np.ndarray | None = None,
    backend: str | None = None,
) -> dict:
    """
    Jointly optimise a pseudo-2D section via Adam.

    Parameters
    ----------
    log_rho_obs : ndarray (S, F)
        Observed log10(rho_a), NaN for missing.
    ph_obs : ndarray (S, F)
        Observed phase in degrees, NaN for missing.
    freqs_t : Tensor or ndarray (F,)
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
    backend : str or None

    Returns
    -------
    dict
        Keys: ``log_rho``, ``log_thick``, ``history``.
    """
    bk = backend or _backend()
    if bk == "torch":
        from ._pinn_ops_torch import (
            _fit_2d_joint_torch as _impl,
        )
    else:
        from ._pinn_ops_tf import (
            _fit_2d_joint_tf as _impl,
        )
    return _impl(
        log_rho_obs,
        ph_obs,
        freqs_t,
        n_layers=n_layers,
        depth_max=depth_max,
        lam_z=lam_z,
        lam_x=lam_x,
        lr=lr,
        epochs=epochs,
        device=device,
        log_every=log_every,
        init_log_rho=init_log_rho,
        init_log_thick=init_log_thick,
    )


# ── 3-D joint ──────────────────────────────────────


def fit_3d_joint(
    log_rho_obs: np.ndarray,
    ph_obs: np.ndarray,
    freqs_t,
    adjacency: np.ndarray,
    *,
    n_layers: int,
    depth_max: float,
    lam_z: float,
    lam_g: float,
    lr: float,
    epochs: int,
    device: str,
    log_every: int,
    init_log_rho: np.ndarray | None = None,
    init_log_thick: np.ndarray | None = None,
    backend: str | None = None,
) -> dict:
    """
    Jointly optimise a quasi-3D section via Adam.

    Parameters
    ----------
    log_rho_obs : ndarray (S, F)
    ph_obs : ndarray (S, F)
    freqs_t : Tensor or ndarray (F,)
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
    backend : str or None

    Returns
    -------
    dict
        Keys: ``log_rho``, ``log_thick``, ``history``.
    """
    bk = backend or _backend()
    if bk == "torch":
        from ._pinn_ops_torch import (
            _fit_3d_joint_torch as _impl,
        )
    else:
        from ._pinn_ops_tf import (
            _fit_3d_joint_tf as _impl,
        )
    return _impl(
        log_rho_obs,
        ph_obs,
        freqs_t,
        adjacency,
        n_layers=n_layers,
        depth_max=depth_max,
        lam_z=lam_z,
        lam_g=lam_g,
        lr=lr,
        epochs=epochs,
        device=device,
        log_every=log_every,
        init_log_rho=init_log_rho,
        init_log_thick=init_log_thick,
    )
