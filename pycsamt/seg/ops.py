# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

from __future__ import annotations

from collections.abc import Iterable

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
    "synthesize_spectra_from_z",
]

# magnetic permeability of free space [H/m]
MU0: float = 4.0e-7 * np.pi


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


def z_to_rho_phi(
    z: np.ndarray,
    freq: np.ndarray | Iterable[float],
    *,
    mu0: float = MU0,
) -> tuple[np.ndarray, np.ndarray]:
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
    mag2 = Z.real**2 + Z.imag**2
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


def time_vector(npts: int, dt: float) -> np.ndarray:
    """Evenly spaced time vector [s] from N and Δt."""
    n = int(max(0, npts))
    return np.arange(n, dtype=float) * float(dt)


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
    den = np.abs(Sxx) * np.abs(Syy)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.where(den > 0.0, num / den, 0.0)
    return np.clip(out, 0.0, 1.0)


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


def amp_or_psd(
    x: np.ndarray,
    fs: float,
    *,
    mode: str = "amp",
) -> tuple[np.ndarray, np.ndarray]:
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


def synthesize_spectra_from_z(
    Z: np.ndarray,
    *,
    S_HH: np.ndarray | None = None,
    H_psd: tuple[np.ndarray, np.ndarray, np.ndarray | None] | None = None,
    tipper: np.ndarray | None = None,  # shape (nf,1,2) or (nf,2) accepted
    include_hz: bool = False,
    chan_order: tuple[str, ...] = ("HX", "HY", "EX", "EY"),
    e_noise: float | np.ndarray | None = None,
    h_noise: float | np.ndarray | None = None,
) -> tuple[np.ndarray, list[str]]:
    r"""
    Synthesize a full Hermitian cross–spectral tensor from an
    impedance tensor ``Z(f)`` and optional tipper.

    Given electric–magnetic relation ``E = Z H``, blocks are:

    * ``S_EH = Z S_HH``
    * ``S_EE = Z S_HH Z^H``

    If a horizontal tipper ``T`` (TX, TY) is supplied and HZ is
    requested, additional blocks are:

    * ``S_ZH = T S_HH``, ``S_ZZ = T S_HH T^H``
    * ``S_EZ = Z S_HH T^H``

    If magnetic spectra are not provided, a unit–power model is
    used (``S_HH = I``).  The result is explicitly symmetrized
    to be Hermitian.

    Parameters
    ----------
    Z : ndarray
        Impedance per frequency.  Shape ``(n, 2, 2)`` or
        ``(2, 2)``.  If 2–D, it is broadcast to one
        frequency.
    S_HH : ndarray, optional
        Magnetic spectra per frequency, shape ``(n, 2, 2)``.
        Must be Hermitian.  Overrides ``H_psd`` if both are
        given.
    H_psd : tuple of ndarray, optional
        ``(Pxx, Pyy, Pxy)`` where ``Pxx, Pyy`` are real
        ``(n,)`` and ``Pxy`` is optional complex ``(n,)``.
        Used to assemble ``S_HH`` when it is not supplied.
    tipper : ndarray, optional
        Horizontal tipper as ``(n, 1, 2)`` or ``(n, 2)``.
        Required when ``include_hz`` is ``True`` or when
        ``"HZ"`` is in ``chan_order``.
    include_hz : bool, default False
        If ``True``, add HZ and its cross–terms.  Requires
        ``tipper``.
    chan_order : tuple of str, optional
        Channel order of the output spectra.  Defaults to
        ``("HX","HY","EX","EY")``.  Include ``"HZ"`` to
        place the vertical magnetic channel.
    e_noise : float or ndarray, optional
        Diagonal noise power added to the electric block
        ``S_EE``.  Scalar or length ``n``.
    h_noise : float or ndarray, optional
        Diagonal noise power added to the magnetic block
        ``S_HH``.  Scalar or length ``n``.

    Returns
    -------
    Sfull : ndarray
        Cross–spectral tensor, shape ``(n, m, m)``, complex,
        Hermitian.  ``m = len(chan_order)``.
    chan_ids : list of str
        Channel names in the order used for ``Sfull``.

    Raises
    ------
    ValueError
        If shapes are incompatible, the tipper is missing
        when HZ is requested, or dimensions do not match.

    Notes
    -----
    Absolute scaling of spectra depends on ``S_HH``.  If
    neither ``S_HH`` nor ``H_psd`` is provided, the unit–
    power model is used.  This is convenient for tests and
    structure checks, but does not preserve physical power.

    The function enforces Hermitian symmetry by conjugating
    off–diagonal blocks and real–ifying the diagonal.

    Examples
    --------
    Minimal, unit–power synthesis:

    >>> S, order = synthesize_spectra_from_z(Z_arr)

    With magnetic PSDs and tipper, including HZ:

    >>> Pxx = np.full(n, 1e-2)
    >>> Pyy = np.full(n, 1e-2)
    >>> S, order = synthesize_spectra_from_z(
    ...     Z_arr,
    ...     H_psd=(Pxx, Pyy, None),
    ...     tipper=T_arr,  # (n,1,2) or (n,2)
    ...     include_hz=True,
    ...     chan_order=("HX", "HY", "HZ", "EX", "EY"),
    ... )

    See Also
    --------
    spectra_from_Z
        Builds a :class:`~pycsamt.seg.spectra.Spectra` using
        this tensor.
    Spectra.from_Z
        Class method wrapper returning a spectra object.
    Spectra.to_Z
        Inverse operation (spectra → transfer functions).

    References
    ----------
    .. [1] Chave, A. D., & Jones, A. G. (2012). *The
           Magnetotelluric Method: Theory and Practice*.
           Cambridge Univ. Press.
    .. [2] Bendat, J. S., & Piersol, A. G. (2011).
           *Random Data: Analysis and Measurement
           Procedures*. Wiley.
    .. [3] SEG EDI MT/EMAP standard (1987). MTNet.
    """

    Zm = np.asarray(Z, complex)
    if Zm.ndim == 2:
        Zm = Zm[np.newaxis, ...]
    nf = Zm.shape[0]

    # ---- S_HH assembly
    if S_HH is None:
        if H_psd is not None:
            Pxx, Pyy, Pxy = H_psd
            Pxx = np.asarray(Pxx, float).reshape(nf)
            Pyy = np.asarray(Pyy, float).reshape(nf)
            if Pxy is None:
                Pxy = np.zeros(nf, complex)
            else:
                Pxy = np.asarray(Pxy, complex).reshape(nf)
            S_HH = np.zeros((nf, 2, 2), complex)
            S_HH[:, 0, 0] = Pxx
            S_HH[:, 1, 1] = Pyy
            S_HH[:, 0, 1] = Pxy
            S_HH[:, 1, 0] = np.conjugate(Pxy)
        else:
            # unit-power default (good for tests, not absolute scaling)
            S_HH = np.zeros((nf, 2, 2), complex)
            S_HH[:, 0, 0] = 1.0
            S_HH[:, 1, 1] = 1.0
    else:
        S_HH = np.asarray(S_HH, complex)
        if S_HH.shape != (nf, 2, 2):
            raise ValueError("S_HH must have shape (nf,2,2)")

    # ---- optional diag noise
    def _add_diag_noise(Sblk: np.ndarray, noise):
        if noise is None:
            return
        if np.isscalar(noise):
            Sblk[:, 0, 0] += float(noise)
            Sblk[:, 1, 1] += float(noise)
        else:
            v = np.asarray(noise, float).reshape(nf)
            Sblk[:, 0, 0] += v
            Sblk[:, 1, 1] += v

    # E/H and E/E
    S_EH = np.einsum("fij,fjk->fik", Zm, S_HH)  # Z S_HH
    S_EE = np.einsum("fij,fjk,flk->fil", Zm, S_HH, np.conjugate(Zm))  # Z S_HH Z^H
    _add_diag_noise(S_EE, e_noise)
    _add_diag_noise(S_HH, h_noise)

    # Tipper / HZ (optional)
    T = None
    if include_hz:
        if tipper is None:
            raise ValueError("include_hz=True requires tipper.")
        T = np.asarray(tipper, complex)
        if T.ndim == 2 and T.shape[1] == 2:
            T = T[:, np.newaxis, :]  # (nf,1,2)
        if T.shape != (nf, 1, 2):
            raise ValueError("tipper must have shape (nf,1,2) or (nf,2)")
        S_ZH = np.einsum("fik,fkj->fij", T, S_HH)  # (nf,1,2)
        S_ZZ = np.einsum("fik,fkj,flk->fil", T, S_HH, np.conjugate(T))  # (nf,1,1)
        S_EZ = np.einsum("fij,fkj->fik", Zm, np.conjugate(S_ZH))  # (nf,2,1)

    order = [s.upper() for s in chan_order]
    nchan = len(order)
    Sfull = np.zeros((nf, nchan, nchan), complex)
    idx = {lbl: i for i, lbl in enumerate(order)}

    def _put(lbl_i, lbl_j, arr):
        i = idx[lbl_i]
        j = idx[lbl_j]
        Sfull[:, i, j] = arr

    # H/H
    _put("HX", "HX", S_HH[:, 0, 0])
    _put("HX", "HY", S_HH[:, 0, 1])
    _put("HY", "HX", np.conjugate(S_HH[:, 0, 1]))
    _put("HY", "HY", S_HH[:, 1, 1])
    # E/H
    _put("EX", "HX", S_EH[:, 0, 0])
    _put("EX", "HY", S_EH[:, 0, 1])
    _put("EY", "HX", S_EH[:, 1, 0])
    _put("EY", "HY", S_EH[:, 1, 1])
    # H/E (Hermitian)
    _put("HX", "EX", np.conjugate(S_EH[:, 0, 0]))
    _put("HY", "EX", np.conjugate(S_EH[:, 0, 1]))
    _put("HX", "EY", np.conjugate(S_EH[:, 1, 0]))
    _put("HY", "EY", np.conjugate(S_EH[:, 1, 1]))
    # E/E
    _put("EX", "EX", S_EE[:, 0, 0])
    _put("EX", "EY", S_EE[:, 0, 1])
    _put("EY", "EX", np.conjugate(S_EE[:, 0, 1]))
    _put("EY", "EY", S_EE[:, 1, 1])

    # Z blocks if requested
    if "HZ" in order:
        if T is None:
            raise ValueError("HZ requested in chan_order but no tipper provided.")
        _put("HZ", "HX", S_ZH[:, 0, 0])
        _put("HZ", "HY", S_ZH[:, 0, 1])
        _put("HX", "HZ", np.conjugate(S_ZH[:, 0, 0]))
        _put("HY", "HZ", np.conjugate(S_ZH[:, 0, 1]))
        _put("HZ", "HZ", S_ZZ[:, 0, 0])
        _put("EX", "HZ", S_EZ[:, 0, 0])
        _put("EY", "HZ", S_EZ[:, 1, 0])
        _put("HZ", "EX", np.conjugate(S_EZ[:, 0, 0]))
        _put("HZ", "EY", np.conjugate(S_EZ[:, 1, 0]))

    # ensure exact Hermitian symmetry on the diagonal
    d = np.diag_indices(nchan)
    Sfull[:, d[0], d[1]] = np.real(Sfull[:, d[0], d[1]])
    return Sfull, order


def effective_dof_from_meta(
    *,
    segnum: int | np.ndarray | None = None,
    avgt: float | np.ndarray | None = None,
    bw: float | np.ndarray | None = None,
    min_dof: int = 1,
) -> int | np.ndarray | None:
    r"""
    Estimate effective DoF (independent averages) from
    per-frequency metadata. If ``segnum`` is given, use it.
    Else, if both ``avgt`` and ``bw`` are finite, use
    ``round(avgt * bw)``. Returns ``None`` if nothing can be
    inferred.
    """
    if segnum is not None:
        return np.asarray(segnum).astype(int)
    if avgt is None or bw is None:
        return None
    a = np.asarray(avgt, float)
    b = np.asarray(bw, float)
    with np.errstate(invalid="ignore"):
        m = np.rint(a * b)
    m = np.where(np.isfinite(m), m, 0.0).astype(int)
    if np.ndim(m) == 0:
        return int(max(min_dof, m))
    return np.maximum(m, int(min_dof))


def _safe_inv2(
    a: np.ndarray,
    ridge: float | None = None,
) -> np.ndarray:
    r"""
    Invert 2x2 complex matrix with optional Tikhonov
    regularization. Adds ``ridge * I`` if ``ridge`` > 0.
    """
    a = np.asarray(a, complex)
    if ridge is not None and ridge > 0.0:
        a = a + float(ridge) * np.eye(2, dtype=float)
    return np.linalg.inv(a)


def z_error_from_blocks(
    S_EE: np.ndarray,
    S_EH: np.ndarray,
    S_HH: np.ndarray,
    *,
    M: float | None,
    ridge: float | None = None,
) -> np.ndarray:
    r"""
    Diagonal σ estimate for Z components from spectral
    blocks. Uses Whittle/Wishart large-sample approx:

        Var(Z_ij) ≈ (E_ii * G_jj) / M

    with ``G = inv(S_HH) @ inv(S_HH)ᴴ`` and
    ``E = S_EE - S_EH @ inv(S_HH) @ S_EHᴴ``.
    Returns a (2, 2) float array. Zeros if ``M`` is falsy.
    """
    if M is None or M <= 0:
        return np.zeros((2, 2), float)

    Ginv = _safe_inv2(S_HH, ridge=ridge)
    E = S_EE - S_EH @ Ginv @ np.conjugate(S_EH.T)

    gdiag = np.diag(Ginv @ np.conjugate(Ginv.T)).real
    ediag = np.diag(E).real

    out = np.zeros((2, 2), float)
    for i in range(2):
        for j in range(2):
            v = (ediag[i] * gdiag[j]) / float(M)
            out[i, j] = np.sqrt(max(0.0, float(v)))
    return out


def tipper_error_from_blocks(
    S_ZH: np.ndarray,
    S_ZZ: np.ndarray,
    S_HH: np.ndarray,
    *,
    M: float | None,
    ridge: float | None = None,
) -> np.ndarray:
    r"""
    Diagonal σ estimate for tipper components. Uses

        Var(T_j) ≈ (Ezz * G_jj) / M

    with ``Ezz = S_ZZ - S_ZH @ inv(S_HH) @ S_ZHᴴ`` and
    ``G = inv(S_HH) @ inv(S_HH)ᴴ``. Returns shape (1, 2).
    Zeros if ``M`` is falsy.
    """
    if M is None or M <= 0:
        return np.zeros((1, 2), float)

    Ginv = _safe_inv2(S_HH, ridge=ridge)
    # ensure scalar (1x1) for S_ZZ
    S_ZZ = np.asarray(S_ZZ, complex)
    if S_ZZ.shape == ():
        S_ZZ = S_ZZ.reshape(1, 1)

    Ezz = S_ZZ - S_ZH @ Ginv @ np.conjugate(S_ZH.T)
    ezz = float(np.real(Ezz.squeeze()))
    gdiag = np.diag(Ginv @ np.conjugate(Ginv.T)).real

    out = np.zeros((1, 2), float)
    for j in range(2):
        v = (ezz * gdiag[j]) / float(M)
        out[0, j] = np.sqrt(max(0.0, float(v)))
    return out


def compute_errors_from_S(
    S: np.ndarray,
    e_idx: tuple[int, int],
    h_idx: tuple[int, int],
    hz_idx: int | None = None,
    *,
    M: float | None,
    ridge: float | None = None,
) -> tuple[np.ndarray, np.ndarray | None]:
    r"""
    Estimate 1-sigma uncertainties for the impedance tensor
    and, optionally, the tipper at a single frequency.

    The function operates on one complex, Hermitian cross-
    spectral density matrix ``S`` of shape ``(nchan, nchan)``.
    It extracts the required sub-blocks using ``e_idx``,
    ``h_idx`` and (optionally) ``hz_idx``, applies an
    optional ridge term to stabilize the magnetic sub-block,
    and performs first-order error propagation under a
    complex-Wishart model with effective degrees of freedom
    ``M``.

    Parameters
    ----------
    S : ndarray
        Complex Hermitian spectra of shape ``(nchan, nchan)``.
    e_idx : tuple of int
        Indices of the two electric channels, ordered as
        ``(EX, EY)`` within ``S``.
    h_idx : tuple of int
        Indices of the two horizontal magnetic channels,
        ordered as ``(HX, HY)``.
    hz_idx : int, optional
        Index of the vertical magnetic channel ``HZ`` within
        ``S``.  If omitted, no tipper uncertainty is computed.
    M : float, optional
        Effective degrees of freedom used for the variance
        scaling.  Larger values reduce the estimated error as
        ``1 / sqrt(M)``.  If ``None``, no uncertainty is
        computed and NaNs are returned.
    ridge : float, optional
        Non-negative Tikhonov regularization added to the
        magnetic auto/cross block prior to inversion,
        ``S_HH + ridge * I``.

    Returns
    -------
    z_err : ndarray
        Real array of shape ``(2, 2)`` with the 1-sigma
        standard errors for ``Z`` components in the order
        ``[[zxx, zxy], [zyx, zyy]]``.
    tip_err : ndarray or None
        Real array of shape ``(1, 2)`` with the 1-sigma
        standard errors for ``(tx, ty)``.  ``None`` if
        ``hz_idx`` is not provided.

    Raises
    ------
    ValueError
        If shapes or indices are inconsistent with ``S``.
    LinAlgError
        If the stabilized magnetic block is singular.

    Notes
    -----
    Let ``Z = S_EH @ inv(S_HH)`` and
    ``T = S_ZH @ inv(S_HH)``.  Under a complex-Wishart noise
    model the covariance of the least-squares estimators can
    be approximated with standard first-order error
    propagation, leading to uncertainties that scale with
    ``1 / sqrt(M)``.  The returned errors are real magnitudes
    (per complex component) suitable for populating
    ``Z.z_err`` and tipper error arrays.

    Examples
    --------
    >>> z_e, t_e = compute_errors_from_S(
    ...     S,
    ...     e_idx=(2, 3),
    ...     h_idx=(0, 1),
    ...     hz_idx=4,
    ...     M=24.0,
    ...     ridge=1e-6,
    ... )
    >>> z_e.shape, (t_e is None) or t_e.shape
    ((2, 2), (1, 2))

    See Also
    --------
    effective_dof_from_meta
        Helper to infer ``M`` from metadata.
    Spectra.to_Z
        Uses this routine when ``estimate_error=True``.

    References
    ----------
    .. [1] Chave, A. D., & Jones, A. G. (2012). *The
           Magnetotelluric Method: Theory and Practice*.
           Cambridge University Press.
    .. [2] Bendat, J. S., & Piersol, A. G. (2011). *Random
           Data: Analysis and Measurement Procedures*. Wiley.
    """

    S = np.asarray(S, complex)

    ex, ey = e_idx
    hx, hy = h_idx

    H = S[np.ix_((hx, hy), (hx, hy))]
    EH = S[np.ix_((ex, ey), (hx, hy))]
    EE = S[np.ix_((ex, ey), (ex, ey))]

    z_err = z_error_from_blocks(
        EE,
        EH,
        H,
        M=M,
        ridge=ridge,
    )

    tip_err = None
    if hz_idx is not None:
        ZH = S[np.ix_((hz_idx,), (hx, hy))]
        ZZ = S[np.ix_((hz_idx,), (hz_idx,))]
        tip_err = tipper_error_from_blocks(
            ZH,
            ZZ,
            H,
            M=M,
            ridge=ridge,
        )

    return z_err, tip_err
