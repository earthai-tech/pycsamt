# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

from __future__ import annotations

from typing import Iterable, Tuple
import numpy as np

__all__ = [
    "MU0",
    "rotate_impedance",
    "rotate_tipper",
    "z_to_rho_phi",
    "rho_phi_to_z",
    "time_vector",
    "coherence_ms",
    "pack_hermitian",
    "unpack_hermitian",
    "rotate_spectra",
    "amp_or_psd",
]

# magnetic permeability of free space [H/m]
MU0: float = 4.0e-7 * np.pi


# ----------------------- rotations ----------------------------

def _rot_mat(theta_deg: float | np.ndarray) -> np.ndarray:
    """2-D rotation matrix R(θ) with our SEG axis convention.

    R = [[ cosθ,  sinθ],
         [-sinθ,  cosθ]]
    """
    th = np.asarray(theta_deg, float) * np.pi / 180.0
    c = np.cos(th)
    s = np.sin(th)
    R = np.empty((2, 2) + th.shape, float)
    R[0, 0] = c
    R[0, 1] = s
    R[1, 0] = -s
    R[1, 1] = c
    return np.moveaxis(R, (2, 3), (0, 1)) if R.ndim > 2 else R


def rotate_impedance(
    z: np.ndarray,
    theta_deg: float | np.ndarray,
) -> np.ndarray:
    """Rotate 2x2 impedance per freq: Z' = R Z R^T.

    Shapes supported:
      - z: (nfreq, 2, 2) or (2, 2)
      - theta: scalar or (nfreq,)
    """
    Z = np.asarray(z)
    if Z.ndim == 2:
        Z = Z[np.newaxis, ...]
    n = Z.shape[0]
    th = np.asarray(theta_deg, float)
    if th.ndim == 0:
        th = np.full(n, float(th))
    R = np.stack([_rot_mat(t) for t in th], axis=0)
    Rt = np.swapaxes(R, -1, -2)
    return (R @ Z @ Rt).squeeze()


def rotate_tipper(
    t: np.ndarray,
    theta_deg: float | np.ndarray,
) -> np.ndarray:
    """Rotate tipper horizontal vector(s): t' = R t.

    Accepts:
      - t: (nfreq, 2) or (nfreq, 1, 2) or (2,)
      - theta: scalar or (nfreq,)
    """
    T = np.asarray(t)
    if T.ndim == 3 and T.shape[1] == 1:
        T = T[:, 0, :]
    if T.ndim == 1:
        T = T[np.newaxis, ...]
    n = T.shape[0]
    th = np.asarray(theta_deg, float)
    if th.ndim == 0:
        th = np.full(n, float(th))
    R = np.stack([_rot_mat(t0) for t0 in th], axis=0)
    out = (R @ T[..., np.newaxis]).squeeze(-1)
    return out.squeeze()


# ------------------- Z <-> (rho, phi) -------------------------

def z_to_rho_phi(
    z: np.ndarray,
    freq: np.ndarray | Iterable[float],
    *,
    mu0: float = MU0,
) -> Tuple[np.ndarray, np.ndarray]:
    """Component-wise apparent resistivity and phase from Z.

    ρ = |Z|^2 / (μ0 ω) ;  φ = atan2(Im(Z), Re(Z)) [deg]
    """
    Z = np.asarray(z)
    if Z.ndim == 2:
        Z = Z[np.newaxis, ...]
    f = np.asarray(freq, float).reshape(-1)
    if f.size != Z.shape[0]:
        raise ValueError("freq and Z length mismatch")
    w = 2.0 * np.pi * f
    mag2 = (Z.real ** 2 + Z.imag ** 2)
    rho = mag2 / (mu0 * w)[:, None, None]
    phi = np.degrees(np.arctan2(Z.imag, Z.real))
    return rho.squeeze(), phi.squeeze()


def rho_phi_to_z(
    rho: np.ndarray,
    phi_deg: np.ndarray,
    freq: np.ndarray | Iterable[float],
    *,
    mu0: float = MU0,
) -> np.ndarray:
    """Build Z from ρ, φ, and freq (per component)."""
    R = np.asarray(rho, float)
    P = np.asarray(phi_deg, float)
    if R.ndim == 2:
        R = R[np.newaxis, ...]
        P = P[np.newaxis, ...]
    f = np.asarray(freq, float).reshape(-1)
    if f.size != R.shape[0]:
        raise ValueError("freq and rho length mismatch")
    w = 2.0 * np.pi * f
    mag = np.sqrt(R * (mu0 * w)[:, None, None])
    ph = np.radians(P)
    Z = mag * (np.cos(ph) + 1j * np.sin(ph))
    return Z.squeeze()


# ---------------------- time-series ---------------------------

def time_vector(npts: int, dt: float) -> np.ndarray:
    """Evenly spaced time vector [s] from N and Δt."""
    n = int(max(0, npts))
    return np.arange(n, dtype=float) * float(dt)


# ------------------------ coherency ---------------------------

def coherence_ms(
    Sxy: np.ndarray,
    Sxx: np.ndarray,
    Syy: np.ndarray,
) -> np.ndarray:
    """Magnitude-squared coherency γ^2."""
    Sxy = np.asarray(Sxy)
    Sxx = np.asarray(Sxx)
    Syy = np.asarray(Syy)
    num = np.abs(Sxy) ** 2
    den = (np.abs(Sxx) * np.abs(Syy))
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.where(den > 0.0, num / den, 0.0)
    return np.clip(out, 0.0, 1.0)


# -------------------- spectra packing -------------------------

def pack_hermitian(C: np.ndarray) -> np.ndarray:
    """Pack Hermitian NxN into upper-triangular vector.

    Output length = N*(N+1)/2, complex dtype preserved.
    """
    A = np.asarray(C)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("C must be square (N x N)")
    n = A.shape[0]
    idx = np.triu_indices(n)
    return A[idx]


def unpack_hermitian(v: np.ndarray, n: int) -> np.ndarray:
    """Unpack upper-triangular vector into Hermitian matrix."""
    v = np.asarray(v)
    tri_len = n * (n + 1) // 2
    if v.size != tri_len:
        raise ValueError("vector size incompatible with n")
    A = np.zeros((n, n), complex)
    iu = np.triu_indices(n)
    A[iu] = v
    # mirror to lower triangle (conjugate symmetry)
    il = (iu[1], iu[0])
    A[il] = np.conjugate(A[iu])
    # ensure diagonal is real-valued
    d = np.diag_indices(n)
    A[d] = np.real(A[d])
    return A


def rotate_spectra(
    C: np.ndarray,
    theta_deg: float,
) -> np.ndarray:
    """Rotate cross-power matrix: C' = R C Rᵀ (per freq)."""
    A = np.asarray(C)
    if A.ndim == 2:
        A = A[np.newaxis, ...]
    n = A.shape[0]
    R = _rot_mat(theta_deg)
    if R.ndim == 2:
        R = np.broadcast_to(R, (n, 2, 2))
    Rt = np.swapaxes(R, -1, -2)
    return (R @ A @ Rt).squeeze()


# ------------------------ signals -----------------------------

def amp_or_psd(
    x: np.ndarray,
    fs: float,
    *,
    mode: str = "amp",
) -> Tuple[np.ndarray, np.ndarray]:
    """Return (f, y) where y is amplitude or PSD.

    mode="amp": |X(f)| / N
    mode="psd": (|X(f)|^2) / (fs * N)
    """
    sig = np.asarray(x, float).ravel()
    N = sig.size
    if N == 0:
        return np.asarray([]), np.asarray([])
    X = np.fft.rfft(sig)
    f = np.fft.rfftfreq(N, d=1.0 / float(fs))
    if mode.lower() == "psd":
        y = (np.abs(X) ** 2) / (float(fs) * N)
    else:
        y = np.abs(X) / max(1, N)
    return f, y
