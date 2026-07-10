# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
General utilities for working with impedance tensors (Z) and
tippers.

This module provides small, composable helpers that complement
the object APIs in :mod:`pycsamt.z`.
"""
from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from ..exceptions import ZError
from ..log.logger import get_logger
from ..utils.zmath import (
    invertmatrix_incl_errors,
    rotatematrix_incl_errors,
    z_error2r_phi_error,
)

logger = get_logger(__name__)

# Orientation correction
def correct_for_sensor_orientation(
    z_prime: np.ndarray,
    bx: float = 0.0,
    by: float = 90.0,
    ex: float = 0.0,
    ey: float = 90.0,
    z_prime_err: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray | None]:
    """
    Correct an impedance tensor for sensor misorientation.

    Given measured fields in non-standard sensor frames,

    .. math::

       \\mathbf{E}' = \\mathbf{Z}' \\, \\mathbf{B}',

    and change-of-basis matrices :math:`\\mathbf{T}` (for E) and
    :math:`\\mathbf{U}` (for B), the tensor in standard coordinates
    is

    .. math::

       \\mathbf{Z} = \\mathbf{T} \\, \\mathbf{Z}' \\, \\mathbf{U}^{-1}.

    This function builds :math:`\\mathbf{T}` and :math:`\\mathbf{U}`
    from the sensor orientations and applies the transform.

    Parameters
    ----------
    z_prime : ndarray
        Input impedance tensor(s). Shape ``(2, 2)`` or
        ``(n_freq, 2, 2)``. Complex valued.
    bx, by : float, default 0, 90
        Orientations (deg) of Bx, By relative to geographic
        North (0°). Positive clockwise.
    ex, ey : float, default 0, 90
        Orientations (deg) of Ex, Ey relative to geographic
        North (0°). Positive clockwise.
    z_prime_err : ndarray, optional
        Standard deviations for :math:`\\mathbf{Z}'`, same shape as
        ``z_prime`` and real-valued. If provided, a conservative
        error propagation is applied (see Notes).

    Returns
    -------
    z : ndarray
        Orientation-corrected tensor(s), same shape as ``z_prime``.
    z_err : ndarray or None
        Propagated standard deviations in the default orientation.
        ``None`` if ``z_prime_err`` was not given.

    Notes
    -----
    - Matrices :math:`\\mathbf{T}` and :math:`\\mathbf{U}` are real and
      built from unit vectors at angles ``ex, ey`` and ``bx, by``
      respectively:

      .. code-block:: text

         T = [[cos(ex), cos(ey)],
              [sin(ex), sin(ey)]],
         U = [[cos(bx), cos(by)],
              [sin(bx), sin(by)]]

    - Error propagation uses a simple 1-norm–like bound for
      :math:`Z = T Z' U^{-1}`:

      .. math::

         \\Delta Z \\approx |T| \\, \\Delta Z' \\, |U^{-1}|,

      applied element-wise by summation. This is conservative and
      ignores covariances.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.z.utils import correct_for_sensor_orientation
    >>> Zp = np.array([[0+0j, 1+1j], [-1-1j, 0+0j]])
    >>> Z, Zerr = correct_for_sensor_orientation(Zp, bx=5, by=95)
    """

    Zp = np.asarray(z_prime)
    if Zp.ndim == 2:
        if Zp.shape != (2, 2):
            raise ZError(
                "For 2-D input, 'z_prime' shape must be (2, 2); "
                f"got {Zp.shape!r}."
            )
        Zp = Zp[None, ...]  # promote to (n,2,2)
        squeeze_output = True
    elif Zp.ndim == 3 and Zp.shape[1:] == (2, 2):
        squeeze_output = False
    else:
        raise ZError(
            "'z_prime' must have shape (2, 2) or (n_freq, 2, 2); "
            f"got {Zp.shape!r}."
        )

    if not np.issubdtype(Zp.dtype, np.complexfloating):
        # accept real, but cast to complex for safety
        Zp = Zp.astype(complex, copy=False)

    if z_prime_err is not None:
        Ze = np.asarray(z_prime_err)
        if Ze.ndim == 2:
            if Ze.shape != (2, 2):
                raise ZError(
                    "For 2-D input, 'z_prime_err' must be (2, 2); "
                    f"got {Ze.shape!r}."
                )
            Ze = Ze[None, ...]
        elif Ze.ndim == 3 and Ze.shape[1:] == (2, 2):
            pass
        else:
            raise ZError(
                "'z_prime_err' must match (2, 2) or (n, 2, 2); "
                f"got {Ze.shape!r}."
            )
        if Ze.shape[0] != Zp.shape[0]:
            raise ZError(
                "'z_prime_err' and 'z_prime' must have matching "
                f"n_freq: {Ze.shape[0]} vs {Zp.shape[0]}."
            )
        if not np.issubdtype(Ze.dtype, np.floating):
            Ze = Ze.astype(float, copy=False)
    else:
        Ze = None

    # build T and U transformation
    def _rotmat(ax: float, ay: float) -> np.ndarray:
        # angles in degrees (clockwise positive)
        a = np.deg2rad(ax)
        b = np.deg2rad(ay)
        return np.array([[np.cos(a), np.cos(b)], [np.sin(a), np.sin(b)]],
                        dtype=float)

    T = _rotmat(ex, ey)
    U = _rotmat(bx, by)

    # Inverse of U
    try:
        Uinv = np.linalg.inv(U)
    except np.linalg.LinAlgError as exc:
        raise ZError(
            "Magnetic sensor basis 'U' is singular; cannot invert. "
            "Check angles 'bx' and 'by'."
        ) from exc

    # apply transformation
    # Z = T * Z' * U^{-1}
    # vectorized over frequencies
    Z_left = np.einsum("ik,nkj->nij", T, Zp, optimize=True)
    Z = np.einsum("nij,jl->nil", Z_left, Uinv, optimize=True)

    # error propagation
    if Ze is not None:
        # ΔZ ≈ |T| ΔZ' |U^{-1}|  (elementwise bound by summation)
        Tabs = np.abs(T)
        Uabs = np.abs(Uinv)
        # First contraction over k: (i,k) x (n,k,l) -> (i,n,l)
        tmp = np.tensordot(Tabs, Ze, axes=([1], [1]))  # (2, n, 2)
        # Then over l: (i,n,l) x (l,j) -> (i,n,j)
        dZ = np.tensordot(tmp, Uabs, axes=([2], [0]))  # (2, n, 2)
        Zerr = np.transpose(dZ, (1, 0, 2))             # (n, 2, 2)
    else:
        Zerr = None

    # finish
    if squeeze_output:
        Z = Z[0]
        if Zerr is not None:
            Zerr = Zerr[0]

    return Z, Zerr


# Frequency / period helpers
def periods_from_freq(
    freq: Sequence[float], *, log10: bool = False
) -> np.ndarray:
    """
    Convert frequency (Hz) to period (s).

    Parameters
    ----------
    freq : array-like of float
        Frequencies in Hz.
    log10 : bool, default False
        If ``True``, return ``log10(period)``.

    Returns
    -------
    ndarray
        Periods with the same length as ``freq``.
    """
    f = np.asarray(freq, dtype=float).ravel()
    if np.any(f <= 0.0):
        raise ZError("Frequencies must be strictly positive.")
    T = 1.0 / f
    return np.log10(T) if log10 else T


def freq_from_periods(
    period: Sequence[float], *, log10: bool = False
) -> np.ndarray:
    """
    Convert period (s) to frequency (Hz).

    Parameters
    ----------
    period : array-like of float
        Periods in seconds. If ``log10=True``, the input is
        interpreted as ``log10(period)``.
    log10 : bool, default False
        Whether input is in base-10 log scale.

    Returns
    -------
    ndarray
        Frequencies in Hz.
    """
    p = np.asarray(period, dtype=float).ravel()
    P = np.power(10.0, p) if log10 else p
    if np.any(P <= 0.0):
        raise ZError("Periods must be strictly positive.")
    return 1.0 / P



# Shape / validation helpers
def ensure_z3(z: np.ndarray) -> np.ndarray:
    """
    Normalize an impedance array to shape ``(n_freq, 2, 2)``.

    Parameters
    ----------
    z : ndarray
        Array with shape ``(2, 2)`` or ``(n_freq, 2, 2)`` (real or
        complex).

    Returns
    -------
    ndarray
        Complex array of shape ``(n_freq, 2, 2)``.

    Raises
    ------
    ZError
        If the shape is not supported.
    """
    arr = np.asarray(z)
    if arr.ndim == 2 and arr.shape == (2, 2):
        arr = arr[None, ...]
    elif not (arr.ndim == 3 and arr.shape[1:] == (2, 2)):
        raise ZError(
            "Z must be shape (2, 2) or (n_freq, 2, 2); got "
            f"{arr.shape!r}."
        )
    return arr.astype(complex, copy=False)


# Stateless ρ–φ from Z (optionally with errors)
def rho_phi_from_z(
    z: np.ndarray,
    freq: Sequence[float],
    z_err: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None, np.ndarray | None]:
    """
    Compute apparent resistivity and phase from Z.

    Parameters
    ----------
    z : ndarray
        Impedance tensor(s), shape ``(2, 2)`` or ``(n, 2, 2)``.
    freq : array-like of float
        Frequency vector (Hz). Length must be ``n``.
    z_err : ndarray, optional
        Absolute errors for Z, same shape as ``z``. If given,
        uncertainties are propagated component-wise.

    Returns
    -------
    rho : ndarray
        Apparent resistivity (Ohm·m), shape ``(n, 2, 2)``.
    phi : ndarray
        Phase (deg), shape ``(n, 2, 2)``.
    rho_err : ndarray or None
        Error on resistivity (Ohm·m) or ``None`` if ``z_err`` is
        ``None``.
    phi_err : ndarray or None
        Error on phase (deg) or ``None`` if ``z_err`` is ``None``.
    """
    Z = ensure_z3(z)
    f = np.asarray(freq, dtype=float).ravel()
    if f.size != Z.shape[0]:
        raise ZError(
            "Length of 'freq' must match number of Z stacks: "
            f"{f.size} vs {Z.shape[0]}."
        )
    # ρ = |Z|^2 / (μ0 ω) ; use 0.2 / f (μ0≈4πe-7, ω=2πf)
    rho = (np.abs(Z) ** 2) * (0.2 / f)[:, None, None]
    phi = np.degrees(np.angle(Z))

    if z_err is None:
        return rho, phi, None, None

    E = ensure_z3(z_err).astype(float, copy=False)
    rho_e = np.zeros_like(rho, dtype=float)
    phi_e = np.zeros_like(phi, dtype=float)
    for k in range(Z.shape[0]):
        for i in range(2):
            for j in range(2):
                r_rel, p_err = z_error2r_phi_error(
                    Z[k, i, j].real, Z[k, i, j].imag, E[k, i, j]
                )
                rho_e[k, i, j] = rho[k, i, j] * r_rel
                phi_e[k, i, j] = p_err
    return rho, phi, rho_e, phi_e



# Antisymmetry enforcement (Zxy ~ -Zyx)
def enforce_offdiag_antisymmetry(
    z: np.ndarray, z_err: np.ndarray | None = None
) -> tuple[np.ndarray, np.ndarray | None]:
    """
    Enforce the magnetotelluric property ``Zxy = -Zyx``.

    Parameters
    ----------
    z : ndarray
        Z tensor(s), shape ``(2, 2)`` or ``(n, 2, 2)``.
    z_err : ndarray, optional
        Errors with same shape as ``z``. If provided, off-diagonal
        errors are combined via RMS.

    Returns
    -------
    z_new : ndarray
        Copy of Z with enforced antisymmetry on off-diagonals.
    z_err_new : ndarray or None
        Propagated errors (off-diagonals via RMS), or ``None`` if
        ``z_err`` was not given.

    Notes
    -----
    We set

    .. math::

       m = \\tfrac{1}{2}(Z_{xy} - Z_{yx}),\\; Z_{xy}=m,\\; Z_{yx}=-m,

    which minimally enforces antisymmetry.
    """
    Z = ensure_z3(z)
    Zc = Z.copy()

    xy = Zc[:, 0, 1]
    yx = Zc[:, 1, 0]
    m = 0.5 * (xy - yx)
    Zc[:, 0, 1] = m
    Zc[:, 1, 0] = -m

    if z_err is None:
        return (Zc[0] if z.ndim == 2 else Zc), None

    E = ensure_z3(z_err).astype(float, copy=False)
    Ec = E.copy()
    # combine off-diagonal errors conservatively (RMS)
    exy = Ec[:, 0, 1]
    eyx = Ec[:, 1, 0]
    e_m = np.sqrt(0.5 * (exy**2 + eyx**2))
    Ec[:, 0, 1] = e_m
    Ec[:, 1, 0] = e_m

    if z.ndim == 2:
        return Zc[0], Ec[0]
    return Zc, Ec



# Vectorized rotate / invert with optional errors
def rotate_z(
    z: np.ndarray,
    angle_deg: float | Sequence[float],
    z_err: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray | None]:
    """
    Rotate Z by ``angle_deg`` (CW positive), vectorized over freq.

    Parameters
    ----------
    z : ndarray
        Z tensor(s), shape ``(2, 2)`` or ``(n, 2, 2)``.
    angle_deg : float or sequence of float
        Single angle for all stacks or one per frequency.
    z_err : ndarray, optional
        Error tensor matching ``z`` shape.

    Returns
    -------
    z_rot : ndarray
        Rotated Z with same shape as ``z``.
    z_err_rot : ndarray or None
        Rotated error, if ``z_err`` was provided.
    """
    Z = ensure_z3(z)
    n = Z.shape[0]
    if np.isscalar(angle_deg) or (
        isinstance(angle_deg, (list, tuple)) and len(angle_deg) == 1
    ):
        alphas = np.full(n, float(np.asarray(angle_deg).ravel()[0]), dtype=float)
    else:
        a = np.asarray(angle_deg, dtype=float).ravel()
        if a.size != n:
            raise ZError(
                f"Expected {n} angles, got {a.size}."
            )
        alphas = a

    E = None if z_err is None else ensure_z3(z_err).astype(float, copy=False)

    Zr = np.empty_like(Z)
    Er = None if E is None else np.empty_like(E)
    for k in range(n):
        if E is None:
            Zr[k], _ = rotatematrix_incl_errors(Z[k], alphas[k])
        else:
            Zr[k], Er[k] = rotatematrix_incl_errors(Z[k], alphas[k], E[k])

    if z.ndim == 2:
        return (Zr[0], None if Er is None else Er[0])
    return Zr, Er


def invert_z(
    z: np.ndarray, z_err: np.ndarray | None = None
) -> tuple[np.ndarray, np.ndarray | None]:
    """
    Invert Z (or each Z slice) with optional error propagation.

    Parameters
    ----------
    z : ndarray
        Z tensor(s), shape ``(2, 2)`` or ``(n, 2, 2)``.
    z_err : ndarray, optional
        Error tensor matching ``z`` shape.

    Returns
    -------
    z_inv : ndarray
        Inverse of Z with same shape as ``z``.
    z_inv_err : ndarray or None
        Propagated error of Z inverse, or ``None`` if no input
        error given.

    Raises
    ------
    ZError
        If any slice is singular.
    """
    Z = ensure_z3(z)
    E = None if z_err is None else ensure_z3(z_err).astype(float, copy=False)

    Zi = np.empty_like(Z)
    Ei = None if E is None else np.empty_like(E)

    for k in range(Z.shape[0]):
        try:
            if E is None:
                Zi[k], _ = invertmatrix_incl_errors(Z[k], np.zeros((2, 2)))
            else:
                Zi[k], Ei[k] = invertmatrix_incl_errors(Z[k], E[k])
        except np.linalg.LinAlgError as exc:
            raise ZError("Singular impedance tensor; cannot invert.") from exc

    if z.ndim == 2:
        return (Zi[0], None if Ei is None else Ei[0])
    return Zi, Ei


# Frequency alignment (generic)
def align_frequency_stack(
    ref_freq: Sequence[float],
    freq: Sequence[float],
    z: np.ndarray,
    *,
    fill_value: complex | float = np.nan,
) -> np.ndarray:
    """
    Align a Z component/vector to a reference frequency grid.

    This is a generalization of a "fit tensor" operation: values at
    ``freq`` are placed at the matching indices of ``ref_freq`` and
    missing locations are filled with ``fill_value``.

    Parameters
    ----------
    ref_freq : array-like of float
        Target frequency grid (Hz), length ``N``.
    freq : array-like of float
        Source frequency grid (Hz), length ``M`` with values
        contained in ``ref_freq`` (within exact match).
    z : ndarray
        Data to align. Shape must be compatible with the first
        dimension being ``M`` (e.g., ``(M,)`` or ``(M, 2, 2)``).
    fill_value : scalar, default NaN
        Value used to fill missing entries.

    Returns
    -------
    ndarray
        Aligned array with shape ``(N, ...)``.

    Raises
    ------
    ZError
        If shapes are incompatible.
    """
    f_ref = np.asarray(ref_freq, dtype=float).ravel()
    f = np.asarray(freq, dtype=float).ravel()
    Z = np.asarray(z)

    if f.size != Z.shape[0]:
        raise ZError(
            "Data first dimension must match 'freq' length: "
            f"{Z.shape[0]} vs {f.size}."
        )

    out_shape = (f_ref.size,) + Z.shape[1:]
    out = np.full(out_shape, fill_value, dtype=Z.dtype)

    # exact match indices (assumes identical frequency values)
    idx = {val: i for i, val in enumerate(f_ref)}
    try:
        dst = np.array([idx[val] for val in f], dtype=int)
    except KeyError as exc:
        raise ZError(
            "All 'freq' values must be present in 'ref_freq'."
        ) from exc

    out[dst, ...] = Z
    return out


# Simple finite/outlier masks
def finite_mask(a: np.ndarray) -> np.ndarray:
    """
    Boolean mask where all real and imaginary parts are finite.

    Parameters
    ----------
    a : ndarray
        Input array.

    Returns
    -------
    ndarray of bool
        Mask with ``True`` for finite entries.
    """
    A = np.asarray(a)
    if np.iscomplexobj(A):
        return np.isfinite(A.real) & np.isfinite(A.imag)
    return np.isfinite(A)


def sigma_clip_mask(
    a: np.ndarray, *, nsigma: float = 3.0, axis: int = 0
) -> np.ndarray:
    """
    Sigma-clip mask for outlier rejection.

    Parameters
    ----------
    a : ndarray
        Input array.
    nsigma : float, default 3.0
        Threshold in units of standard deviation.
    axis : int, default 0
        Axis along which statistics are computed.

    Returns
    -------
    ndarray of bool
        Mask of inliers (``True``) and outliers (``False``).
    """
    A = np.asarray(a, dtype=float)
    mu = np.nanmean(A, axis=axis, keepdims=True)
    sd = np.nanstd(A, axis=axis, keepdims=True)
    return np.abs(A - mu) <= nsigma * (sd + 1e-12)
