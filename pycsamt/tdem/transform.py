# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
TEM time-domain → frequency-domain transforms.

Three loop configurations are fully supported:

**Central / coincident loop** (``offset = 0``)
    Receiver at the centre of the transmitter loop.  Uses the Ward &
    Hohmann (1988) late-time formula directly:

    .. math::

        \\rho_a(t) = \\left(
            \\frac{M\\,\\mu_0^{5/2}}
                  {10\\sqrt{\\pi}\\,|\\partial B_z/\\partial t|\\,t^{5/2}}
        \\right)^{2/3}

**In-loop / large-loop** (``0 < offset < inner_radius``)
    Receiver is inside the transmitter loop but not at the centre.
    A Biot-Savart static-field correction factor
    :math:`\\eta = H_z^{\\rm static}(r_x, r_y)\\,/\\,H_z^{\\rm static}(0,0)`
    is applied to the effective transmitter moment:

    .. math::

        \\rho_a(t) = \\left(
            \\frac{M\\,\\eta\\,\\mu_0^{5/2}}
                  {10\\sqrt{\\pi}\\,|\\partial B_z/\\partial t|\\,t^{5/2}}
        \\right)^{2/3}

    The Biot-Savart integral is evaluated analytically for rectangular /
    square loops and numerically for circular loops.

**Offset / separated loop** (``offset ≥ inner_radius``)
    Transmitter and receiver are separated by horizontal distance *d*.
    Uses the Ward & Hohmann (1988) offset-dipole formula:

    .. math::

        \\rho_a(t) = \\left(
            \\frac{M\\,\\mu_0^{5/2}}
                  {20\\sqrt{\\pi}\\,|\\partial B_z/\\partial t|\\,d^3\\,t^{5/2}}
        \\right)^{2/3}

:class:`LateTimeTransform`
    Fast, approximate method — standard industry approach.

:class:`FourierTransform`
    Rigorous Fourier cosine transform with Kramers-Kronig reconstruction.
    Valid at all time gates; first-order waveform deconvolution for all
    waveform types (Fitterman & Stewart 1986).

:class:`TEMtoEDI`
    High-level dispatcher → :class:`~pycsamt.seg.collection.EDICollection`.
"""

from __future__ import annotations

import warnings
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from ._base import TEMSounding

__all__ = [
    "LateTimeTransform",
    "FourierTransform",
    "TEMtoEDI",
    # geometry helpers (importable for testing)
    "build_adjacency",
    "_rho_a_late_time",
    "_rho_a_in_loop",
    "_rho_a_offset_loop",
    "_biot_savart_rect_hz",
    "_biot_savart_circle_hz",
    "_biot_savart_rect_hz_segments",
    "_biot_savart_circle_hz_segments",
    "_in_loop_geometry_factor",
    "_in_loop_geometry_factor_td",
    "_pseudo_freq",
    "_phase_from_rho",
    "_build_z_array",
    "_z_error_from_rho_error",
    "_cosine_transform_1d",
    "_kramers_kronig_re",
    "_waveform_moments",
    "_apply_waveform_correction",
    "MU0",
]

MU0 = 4.0 * np.pi * 1e-7          # H/m


# ---------------------------------------------------------------------------
# Biot-Savart static field helpers
# ---------------------------------------------------------------------------

def _hz_horiz_seg(y0: float, x1: float, x2: float,
                  rx: float, ry: float) -> float:
    """Hz from a horizontal wire segment at y=y0, running x1→x2."""
    d = ry - y0
    if abs(d) < 1e-14:
        return 0.0
    u1, u2 = x1 - rx, x2 - rx
    return (1.0 / d) * (u2 / np.sqrt(d*d + u2*u2) - u1 / np.sqrt(d*d + u1*u1))


def _hz_vert_seg(x0: float, y1: float, y2: float,
                 rx: float, ry: float) -> float:
    """Hz from a vertical wire segment at x=x0, running y1→y2."""
    d = rx - x0
    if abs(d) < 1e-14:
        return 0.0
    u1, u2 = y1 - ry, y2 - ry
    return -(1.0 / d) * (u2 / np.sqrt(d*d + u2*u2) - u1 / np.sqrt(d*d + u1*u1))


def _biot_savart_rect_hz(
    rx: float,
    ry: float,
    a: float,
    b: float,
) -> float:
    r"""
    Static :math:`H_z` per unit current at ``(rx, ry)`` inside a
    rectangular transmitter loop of half-sides ``(a, b)`` centred at
    the origin (analytical Biot-Savart).

    The loop current runs counter-clockwise when viewed from above:
    bottom (+x) → right (+y) → top (−x) → left (−y).

    Parameters
    ----------
    rx, ry : float
        Receiver position in metres.  Must be strictly inside
        the loop (|rx| < a, |ry| < b) to avoid a singularity.
    a, b : float
        Half-sides of the rectangular loop (m).  For a square loop
        pass ``a = b = side / 2``.

    Returns
    -------
    float
        :math:`H_z` in A m⁻¹ per ampere of transmitter current.

    References
    ----------
    .. [1] Griffiths, D.J. (1999). *Introduction to Electrodynamics*,
       3rd ed., Prentice Hall — §5.4 (Biot-Savart law for line current).
    """
    hz = (
        _hz_horiz_seg(-b, -a, +a, rx, ry) +   # bottom: (−a,−b) → (+a,−b)
        _hz_vert_seg(+a,  -b, +b, rx, ry) +    # right:  (+a,−b) → (+a,+b)
        _hz_horiz_seg(+b, +a, -a, rx, ry) +    # top:    (+a,+b) → (−a,+b)
        _hz_vert_seg(-a,  +b, -b, rx, ry)      # left:   (−a,+b) → (−a,−b)
    )
    return hz / (4.0 * np.pi)


def _biot_savart_rect_hz_segments(
    rx: float,
    ry: float,
    a: float,
    b: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    r"""
    Per-segment :math:`H_z / (4\pi)` and perpendicular distances for a
    rectangular transmitter loop of half-sides ``(a, b)``.

    Returns four arrays of length 4 (segments: bottom, right, top, left):

    * ``hz_rx``   — Hz contribution at ``(rx, ry)`` per segment.
    * ``hz_00``   — Hz contribution at ``(0, 0)`` per segment.
    * ``dist_rx`` — Perpendicular distance from ``(rx, ry)`` to each
      segment's line.
    * ``dist_00`` — Perpendicular distance from ``(0, 0)`` to each
      segment's line.

    Used by :func:`_in_loop_geometry_factor_td` to weight each segment's
    contribution by the EM diffusion length at each time gate.
    """
    hz_rx   = np.array([
        _hz_horiz_seg(-b, -a, +a, rx,  ry),    # bottom
        _hz_vert_seg( +a, -b, +b, rx,  ry),    # right
        _hz_horiz_seg(+b, +a, -a, rx,  ry),    # top
        _hz_vert_seg( -a, +b, -b, rx,  ry),    # left
    ]) / (4.0 * np.pi)

    hz_00   = np.array([
        _hz_horiz_seg(-b, -a, +a, 0.0, 0.0),
        _hz_vert_seg( +a, -b, +b, 0.0, 0.0),
        _hz_horiz_seg(+b, +a, -a, 0.0, 0.0),
        _hz_vert_seg( -a, +b, -b, 0.0, 0.0),
    ]) / (4.0 * np.pi)

    # Perpendicular distances from (rx,ry) and (0,0) to each segment's line
    dist_rx = np.array([abs(ry + b), abs(a - rx), abs(b - ry), abs(rx + a)])
    dist_00 = np.array([b,           a,            b,           a          ])

    return hz_rx, hz_00, dist_rx, dist_00


def _biot_savart_circle_hz(
    rx: float,
    ry: float,
    radius: float,
    n_seg: int = 720,
) -> float:
    r"""
    Static :math:`H_z` per unit current at ``(rx, ry)`` inside a
    circular transmitter loop of given ``radius`` (numerical
    Biot-Savart via piecewise-linear integration).

    Parameters
    ----------
    rx, ry : float
        Receiver position (m).  Should be inside the loop (r < radius).
    radius : float
        Loop radius (m).
    n_seg : int
        Number of equal arc segments for the integration.  720 gives
        relative errors < 10⁻⁵ for rx/radius < 0.9.

    Returns
    -------
    float
        :math:`H_z` in A m⁻¹ per ampere of transmitter current.
    """
    theta = np.linspace(0.0, 2.0 * np.pi, n_seg, endpoint=False)
    dtheta = 2.0 * np.pi / n_seg

    xw = radius * np.cos(theta)
    yw = radius * np.sin(theta)
    dlx = -radius * np.sin(theta) * dtheta   # tangent × arc length
    dly = +radius * np.cos(theta) * dtheta

    dx = rx - xw
    dy = ry - yw
    r3 = (dx * dx + dy * dy) ** 1.5

    hz = np.sum((dlx * dy - dly * dx) / r3)
    return hz / (4.0 * np.pi)


def _biot_savart_circle_hz_segments(
    rx: float,
    ry: float,
    radius: float,
    n_seg: int = 360,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    r"""
    Per-element :math:`H_z / (4\pi)` and distances for a circular loop.

    The circle is discretised into ``n_seg`` equal arc elements.  For each
    element the Biot-Savart contribution at ``(rx, ry)`` and at ``(0, 0)``
    is returned together with the element-to-receiver distances.

    Returns
    -------
    hz_rx : ndarray (n_seg,)
    hz_00 : ndarray (n_seg,)
    dist_rx : ndarray (n_seg,)
        Distance from each element to ``(rx, ry)``.
    dist_00 : ndarray (n_seg,)
        Distance from each element to ``(0, 0)`` = ``radius`` (uniform).
    """
    theta  = np.linspace(0.0, 2.0 * np.pi, n_seg, endpoint=False)
    dt     = 2.0 * np.pi / n_seg

    xw  = radius * np.cos(theta)
    yw  = radius * np.sin(theta)
    dlx = -radius * np.sin(theta) * dt
    dly = +radius * np.cos(theta) * dt

    # At (rx, ry)
    dx_rx = rx - xw;  dy_rx = ry - yw
    r3_rx = (dx_rx ** 2 + dy_rx ** 2) ** 1.5
    hz_rx = (dlx * dy_rx - dly * dx_rx) / r3_rx / (4.0 * np.pi)

    # At (0, 0)  — all wire elements at the same distance = radius
    hz_00 = (dlx * (-yw) - dly * (-xw)) / radius ** 3 / (4.0 * np.pi)

    dist_rx = np.sqrt(dx_rx ** 2 + dy_rx ** 2)
    dist_00 = radius * np.ones(n_seg)

    return hz_rx, hz_00, dist_rx, dist_00


def _in_loop_geometry_factor(
    rx: float,
    ry: float,
    loop_shape: str,
    loop_dims: tuple,
) -> float:
    r"""
    Geometric correction factor
    :math:`\eta = H_z^{\rm static}(r_x, r_y)\,/\,H_z^{\rm static}(0, 0)`
    for a receiver inside the transmitter loop at position ``(rx, ry)``
    relative to the loop centre.

    Used by :func:`_rho_a_in_loop` to correct the effective transmitter
    moment for non-central in-loop receiver positions.

    Parameters
    ----------
    rx, ry : float
        Receiver offset from loop centre (m).
    loop_shape : str
        ``'square'``, ``'rectangular'``, or ``'circle'``.
    loop_dims : tuple of float
        Shape parameters:

        * ``'square'``      → ``(side_length,)``
        * ``'rectangular'`` → ``(length, width)``
        * ``'circle'``      → ``(radius,)``

    Returns
    -------
    float
        Correction factor η > 0.  Equal to 1.0 when (rx, ry) = (0, 0).

    Notes
    -----
    At late time (diffusion time >> loop size), η has no effect because
    the induced eddy current "smoke-ring" has expanded far beyond the
    loop and the field is spatially uniform.  The correction is most
    significant at early-to-mid time gates and for large offsets r/L.

    References
    ----------
    .. [1] Nabighian, M.N. & Macnae, J.C. (1991).  Time domain
       electromagnetic prospecting methods.  In *Electromagnetic Methods
       in Applied Geophysics*, Vol. 2, pp. 427–520. SEG.
    .. [2] Spies, B.R. & Frischknecht, F.C. (1991).  Electromagnetic
       sounding.  Ibid., pp. 285–425.
    """
    if abs(rx) < 1e-10 and abs(ry) < 1e-10:
        return 1.0

    if loop_shape in ("square", "rectangular"):
        if loop_shape == "square":
            a = b = float(loop_dims[0]) / 2.0
        else:
            a = float(loop_dims[0]) / 2.0
            b = float(loop_dims[1]) / 2.0
        hz_rx = _biot_savart_rect_hz(rx, ry, a, b)
        hz_00 = _biot_savart_rect_hz(0.0, 0.0, a, b)
    elif loop_shape == "circle":
        r = float(loop_dims[0])
        hz_rx = _biot_savart_circle_hz(rx, ry, r)
        hz_00 = _biot_savart_circle_hz(0.0, 0.0, r)
    else:
        # Unknown shape — treat as equivalent-area square
        a = b = np.sqrt(float(loop_dims[0])) / 2.0
        hz_rx = _biot_savart_rect_hz(rx, ry, a, b)
        hz_00 = _biot_savart_rect_hz(0.0, 0.0, a, b)

    if hz_00 < 1e-30:
        return 1.0
    return float(max(hz_rx / hz_00, 0.01))   # guard against unphysical values


def _in_loop_geometry_factor_td(
    rx: float,
    ry: float,
    loop_shape: str,
    loop_dims: tuple,
    t: np.ndarray,
    rho_a: np.ndarray,
) -> np.ndarray:
    r"""
    Time-dependent geometry factor :math:`\eta(t)` for in-loop receivers.

    Replaces the DC Biot-Savart approximation with a physically correct
    transition between the near-field (early-time) and far-field (late-time)
    regimes.  Each transmitter segment :math:`i` contributes to the field
    at the receiver with a weight

    .. math::

        w_i(t) = \operatorname{erf}\!\left(
                     \frac{\lambda(t)}{d_i}
                 \right), \quad
        \lambda(t) = \sqrt{\frac{4\,\rho_a(t)\,t}{\mu_0}}

    where :math:`d_i` is the perpendicular distance from the receiver to
    segment :math:`i`, and :math:`\lambda(t)` is the EM diffusion length.

    .. math::

        \eta(t) =
            \frac{\displaystyle\sum_i w_i^{(r_x)}(t)\,
                  \Delta H_z^{(i)}(r_x, r_y)}
                 {\displaystyle\sum_i w_i^{(0)}(t)\,
                  \Delta H_z^{(i)}(0, 0)}

    **Limits**:

    * Late time (:math:`\lambda \gg d_i`): :math:`w_i \to 1` for all
      segments → :math:`\eta \to \eta_\mathrm{static}` (Biot-Savart limit).
    * Early time (:math:`\lambda \ll d_i`): :math:`w_i \propto \lambda/d_i`
      → segments at small :math:`d_i` dominate (near-field limit).

    Parameters
    ----------
    rx, ry : float
        Receiver offset from loop centre (m).
    loop_shape : str
        ``'square'``, ``'rectangular'``, or ``'circle'``.
    loop_dims : tuple
        Loop geometry parameters.
    t : ndarray (n,)
        Gate centre times (s).
    rho_a : ndarray (n,)
        Current estimate of apparent resistivity (Ω m).

    Returns
    -------
    eta_t : ndarray (n,)
        Time-dependent geometry factor.

    References
    ----------
    .. [1] Nabighian, M.N. & Macnae, J.C. (1991). SEG, Vol. 2, eq. 4.34.
    .. [2] Christiansen, A.V. et al. (2009). *Geophysics*, 74, F35–F46.
    """
    from scipy.special import erf

    t   = np.asarray(t)
    n   = len(t)

    if abs(rx) < 1e-10 and abs(ry) < 1e-10:
        return np.ones(n)

    # --- get per-segment Hz and distances ---
    if loop_shape in ("square", "rectangular"):
        if loop_shape == "square":
            a = b = float(loop_dims[0]) / 2.0
        else:
            a = float(loop_dims[0]) / 2.0
            b = float(loop_dims[1]) / 2.0
        hz_rx, hz_00, dist_rx, dist_00 = _biot_savart_rect_hz_segments(rx, ry, a, b)
    elif loop_shape == "circle":
        r_loop = float(loop_dims[0])
        hz_rx, hz_00, dist_rx, dist_00 = _biot_savart_circle_hz_segments(rx, ry, r_loop)
    else:
        a = b = np.sqrt(float(loop_dims[0])) / 2.0
        hz_rx, hz_00, dist_rx, dist_00 = _biot_savart_rect_hz_segments(rx, ry, a, b)

    # --- diffusion length per gate ---
    rho_safe = np.maximum(np.asarray(rho_a), 1e-6)
    lam      = np.sqrt(4.0 * rho_safe * t / MU0)        # (n,)

    # --- per-segment weights (n_seg, n) ---
    d_rx = np.maximum(dist_rx[:, None], 1e-12)           # avoid /0
    d_00 = np.maximum(dist_00[:, None], 1e-12)
    w_rx = erf(lam[None, :] / d_rx)                      # (n_seg, n)
    w_00 = erf(lam[None, :] / d_00)

    # --- time-dependent weighted Hz ---
    hz_td_rx = (hz_rx[:, None] * w_rx).sum(axis=0)       # (n,)
    hz_td_00 = (hz_00[:, None] * w_00).sum(axis=0)

    with np.errstate(divide="ignore", invalid="ignore"):
        eta_t = np.where(np.abs(hz_td_00) > 1e-30,
                         hz_td_rx / hz_td_00,
                         np.ones(n))

    return eta_t


def _loop_inner_radius(loop_shape: str, loop_dims: tuple) -> float:
    """
    Maximum Tx–Rx offset for the receiver to still lie inside the loop.

    * ``'square'``      → half the side length
    * ``'rectangular'`` → half the shorter side
    * ``'circle'``      → the loop radius
    """
    if loop_shape == "circle":
        return float(loop_dims[0])
    if loop_shape == "square":
        return float(loop_dims[0]) / 2.0
    # rectangular: min(length, width) / 2
    return min(float(loop_dims[0]), float(loop_dims[1])) / 2.0


# ---------------------------------------------------------------------------
# Apparent-resistivity formulas
# ---------------------------------------------------------------------------

def _rho_a_late_time(
    dBdt: np.ndarray,
    t: np.ndarray,
    moment: float,
) -> np.ndarray:
    r"""
    Late-time apparent resistivity for a coincident / central-loop system.

    Derived from the late-time asymptotic response of a magnetic
    dipole source over a uniform half-space (Ward & Hohmann 1988,
    eq. 4.96; Nabighian & Macnae 1991):

    .. math::

        \rho_a(t) =
            \left(\frac{M\,\mu_0^{5/2}}
                  {10\,\sqrt{\pi}\,|\partial B_z/\partial t|\,t^{5/2}}
            \right)^{2/3}

    Parameters
    ----------
    dBdt : ndarray (n,)
        :math:`\partial B_z/\partial t` in T s⁻¹.
    t : ndarray (n,)
        Time gates in seconds.
    moment : float
        Transmitter magnetic moment :math:`M = I n A` in A m².

    Returns
    -------
    ndarray (n,)
        Apparent resistivity in Ω m.

    References
    ----------
    .. [1] Ward, S.H. & Hohmann, G.W. (1988).  Electromagnetic Theory
       for Geophysical Applications. Chapter 4, eq. 4.96.
    .. [2] Nabighian, M.N. & Macnae, J.C. (1991).  Time domain
       electromagnetic prospecting methods.  SEG, Vol. 2.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        arg = (moment * MU0 ** 2.5) / (
            10.0 * np.sqrt(np.pi) * np.abs(dBdt) * t ** 2.5
        )
        rho = arg ** (2.0 / 3.0)
    return np.where(np.isfinite(rho) & (rho > 0.0), rho, np.nan)


def _rho_a_in_loop(
    dBdt: np.ndarray,
    t: np.ndarray,
    moment: float,
    loop_shape: str,
    loop_dims: tuple,
    rx: float,
    ry: float = 0.0,
    n_iter: int = 3,
) -> np.ndarray:
    r"""
    Apparent resistivity for an in-loop receiver using a time-dependent
    geometry correction.

    Starting from the static Biot-Savart factor
    :math:`\eta_0 = H_z^\mathrm{static}(r_x, r_y) / H_z^\mathrm{static}(0,0)`,
    the algorithm iteratively refines the per-gate factor
    :math:`\eta(t)` via :func:`_in_loop_geometry_factor_td`:

    .. math::

        \rho_a^{(k+1)}(t) = \left(
            \frac{M\,\eta^{(k)}(t)\,\mu_0^{5/2}}
                  {10\,\sqrt{\pi}\,|\partial B_z/\partial t|\,t^{5/2}}
        \right)^{2/3}

    At late time :math:`\eta(t) \to \eta_0` (static limit); at early
    time the correction accounts for the proximity of the receiver to
    individual transmitter segments.  Convergence is typically reached
    in 2–3 iterations.

    Parameters
    ----------
    dBdt : ndarray (n,)
    t : ndarray (n,)
    moment : float
    loop_shape : str
        ``'square'``, ``'rectangular'``, or ``'circle'``.
    loop_dims : tuple
        Loop geometry parameters.
    rx : float
        Receiver x-offset from loop centre (m).
    ry : float
        Receiver y-offset from loop centre (m).  Default 0.
    n_iter : int
        Number of correction iterations (0 = static only).  Default 3.

    Returns
    -------
    ndarray (n,)
        Corrected apparent resistivity in Ω m.

    References
    ----------
    .. [1] Nabighian & Macnae (1991), Table 4.1.
    .. [2] Christiansen, A.V. et al. (2009). *Geophysics*, 74, F35–F46.
    """
    # Seed: static correction
    eta_0 = _in_loop_geometry_factor(rx, ry, loop_shape, loop_dims)
    rho_a = _rho_a_late_time(dBdt, t, moment * eta_0)

    if n_iter == 0 or (abs(rx) < 1e-10 and abs(ry) < 1e-10):
        return rho_a

    for _ in range(n_iter):
        # Fill NaN/non-positive ρ_a with median for stable η(t) computation
        finite = rho_a[np.isfinite(rho_a) & (rho_a > 0.0)]
        fill   = float(np.median(finite)) if len(finite) > 0 else 100.0
        rho_for_eta = np.where(np.isfinite(rho_a) & (rho_a > 0.0), rho_a, fill)

        eta_t = _in_loop_geometry_factor_td(rx, ry, loop_shape, loop_dims, t, rho_for_eta)
        rho_a = _rho_a_late_time(dBdt, t, moment * eta_t)

    return rho_a


def _rho_a_offset_loop(
    dBdt: np.ndarray,
    t: np.ndarray,
    moment: float,
    offset: float,
) -> np.ndarray:
    r"""
    Late-time apparent resistivity for a separated-loop system where
    the receiver is at horizontal distance ``offset`` from the Tx.

    Uses the Ward & Hohmann (1988) late-time formula for a vertical
    magnetic dipole measured at horizontal distance :math:`d`:

    .. math::

        \frac{\partial H_z}{\partial t}
            = -\frac{M\,\mu_0^{5/2}\,\sigma^{3/2}}
                    {20\,\sqrt{\pi}\,d^3\,t^{5/2}}

    Solving for apparent resistivity:

    .. math::

        \rho_a(t) = \left(
            \frac{M\,\mu_0^{5/2}}
                  {20\,\sqrt{\pi}\,|\partial B_z/\partial t|\,d^3\,t^{5/2}}
        \right)^{2/3}

    Parameters
    ----------
    dBdt : ndarray (n,)
        Measured :math:`\partial B_z/\partial t` at the Rx position (T s⁻¹).
    t : ndarray (n,)
        Time gates in seconds.
    moment : float
        Transmitter moment :math:`M = I n A` in A m².
    offset : float
        Horizontal Tx–Rx separation :math:`d` in metres.

    Returns
    -------
    ndarray (n,)
        Apparent resistivity in Ω m.

    Raises
    ------
    ValueError
        If ``offset <= 0``.

    Notes
    -----
    The factor 20 (vs 10 in the central-loop formula) and the
    :math:`d^3` denominator arise from the angular geometry of the
    off-axis dipole field.  The formula is valid when :math:`d` is
    large compared to the loop size and large compared to the skin
    depth (:math:`d\sqrt{\mu_0 \sigma / (4t)} \ll 1`).

    References
    ----------
    .. [1] Ward, S.H. & Hohmann, G.W. (1988).  Chapter 4, eq. 4.97.
    .. [2] Nabighian, M.N. (1979). Quasi-static transient response of
       a conducting half-space. *Geophysics*, 44, 1700–1705.
    """
    if offset <= 0.0:
        raise ValueError(
            f"offset must be > 0 for a separated-loop configuration "
            f"(got {offset:.3g} m).  Use offset=0 for central-loop."
        )
    with np.errstate(divide="ignore", invalid="ignore"):
        arg = (moment * MU0 ** 2.5) / (
            20.0 * np.sqrt(np.pi) * np.abs(dBdt) * (offset ** 3) * t ** 2.5
        )
        rho = arg ** (2.0 / 3.0)
    return np.where(np.isfinite(rho) & (rho > 0.0), rho, np.nan)


# ---------------------------------------------------------------------------
# Shared downstream helpers (unchanged from v1)
# ---------------------------------------------------------------------------

def _pseudo_freq(t: np.ndarray, convention: str = "skin_depth") -> np.ndarray:
    r"""
    Map time gates to equivalent pseudo-frequencies.

    ``"skin_depth"`` (default) :
        :math:`f = 1 / (2\pi t)`

    ``"diffusion"`` :
        :math:`f = 1 / (4\pi t)`
    """
    if convention == "skin_depth":
        return 1.0 / (2.0 * np.pi * t)
    if convention == "diffusion":
        return 1.0 / (4.0 * np.pi * t)
    raise ValueError(
        f"Unknown convention '{convention}'. Use 'skin_depth' or 'diffusion'."
    )


def _phase_from_rho(
    rho_a: np.ndarray,
    freq: np.ndarray,
    mode: str = "homogeneous",
) -> np.ndarray:
    r"""
    Estimate MT phase from the apparent resistivity curve.

    ``"homogeneous"`` — uniform half-space: :math:`\phi_{xy} = 45°`.

    ``"weidelt"`` — dispersion relation (Weidelt 1972, Boehl et al.
    1983):

    .. math::

        \phi(\omega) \approx \frac{\pi}{4}
            + \frac{1}{2}\frac{\partial\ln\rho_a}{\partial\ln\omega}
    """
    if mode == "homogeneous":
        return np.full_like(rho_a, 45.0)

    if mode == "weidelt":
        log_rho = np.log(rho_a)
        log_omega = np.log(2.0 * np.pi * freq)
        d_log_rho = np.gradient(log_rho, log_omega)
        phase_deg = np.degrees(np.pi / 4.0 + 0.5 * d_log_rho)
        phase_deg = np.clip(phase_deg, 0.0, 90.0)
        phase_deg[~np.isfinite(rho_a)] = np.nan
        return phase_deg

    raise ValueError(f"Unknown phase mode '{mode}'. Use 'homogeneous' or 'weidelt'.")


def _build_z_array(
    rho_a: np.ndarray,
    phase_xy_deg: np.ndarray,
    freq: np.ndarray,
) -> np.ndarray:
    r"""
    Construct complex impedance tensor array from apparent resistivity
    and phase (1-D earth assumption).

    .. math::

        |Z_{xy}(\omega)| = \sqrt{\rho_a(\omega)\,\omega\,\mu_0}
    """
    omega = 2.0 * np.pi * freq
    Z_mag = np.sqrt(rho_a * omega * MU0)
    phi_rad = np.radians(phase_xy_deg)
    Z_xy = Z_mag * (np.cos(phi_rad) + 1j * np.sin(phi_rad))
    n = len(freq)
    Z_arr = np.zeros((n, 2, 2), dtype=complex)
    Z_arr[:, 0, 1] = Z_xy
    Z_arr[:, 1, 0] = -Z_xy
    return Z_arr


def _z_error_from_rho_error(
    rho_a: np.ndarray,
    rho_err: np.ndarray,
    freq: np.ndarray,
) -> np.ndarray:
    r"""
    Propagate apparent-resistivity uncertainty to impedance uncertainty.

    .. math::

        \delta|Z| = \frac{\delta\rho_a}{2\rho_a} \cdot |Z|
    """
    omega = 2.0 * np.pi * freq
    Z_mag = np.sqrt(rho_a * omega * MU0)
    with np.errstate(divide="ignore", invalid="ignore"):
        dZ = np.where(rho_a > 0.0, rho_err / (2.0 * rho_a) * Z_mag, np.nan)
    n = len(freq)
    Z_err = np.zeros((n, 2, 2))
    Z_err[:, 0, 1] = dZ
    Z_err[:, 1, 0] = dZ
    return Z_err


# Expose build_adjacency from gcn for import convenience; harmless if AI not installed
try:
    from ..ai.nets.gcn import build_adjacency
except Exception:
    build_adjacency = None   # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Waveform deconvolution helpers
# ---------------------------------------------------------------------------

def _waveform_moments(waveform, n_samples: int = 2000) -> Tuple[float, float]:
    r"""
    Compute the effective time shift and gate-rejection window for any
    transmitter waveform.

    The turn-off impulse response is :math:`w(\tau) = -dI/d\tau \geq 0`.
    Two cases are handled:

    **Turn-off starts at** :math:`t = 0` (e.g. ``RampWaveform``):
    :math:`I(0) > 0`.  The normalised first moment is

    .. math::

        \tau_\mathrm{eff} = \int_0^{\tau_\mathrm{win}} \tau\,w(\tau)\,d\tau

    and :math:`\tau_\mathrm{window}` is the 99th-percentile of the
    cumulative turn-off.  Gates with :math:`t_\mathrm{corr} \leq
    \tau_\mathrm{window}` are discarded.

    **Turn-off completed before** :math:`t = 0` (e.g. ``HalfSineWaveform``):
    :math:`I(0) \approx 0`.  The last turn-off in the preceding half-period
    is identified and its centroid gives :math:`\tau_\mathrm{eff}` (positive,
    measured from centroid to :math:`t = 0`).  No gate rejection is applied
    (:math:`\tau_\mathrm{window} = 0`).

    For a :class:`~pycsamt.tdem.waveform.RampWaveform` the result is
    analytic: :math:`\tau_\mathrm{eff} = \tau_r/2`,
    :math:`\tau_\mathrm{window} = \tau_r`.

    Parameters
    ----------
    waveform :
        Any waveform object with a ``current_at(t)`` method and a
        ``base_frequency`` attribute, or ``None``.
    n_samples : int
        Number of integration points for numerical waveforms.

    Returns
    -------
    tau_eff : float
        Effective time shift (seconds).
    tau_window : float
        Gate-rejection threshold (seconds).
    """
    from .waveform import SquareWaveform, RampWaveform

    if waveform is None or isinstance(waveform, SquareWaveform):
        return 0.0, 0.0

    if isinstance(waveform, RampWaveform):
        t_ramp = float(waveform.ramp_off)
        return t_ramp / 2.0, t_ramp

    hp = 0.5 / float(waveform.base_frequency)   # half-period
    I0 = float(waveform.current_at(np.array([0.0]))[0])

    if I0 > 0.01:
        # Turn-off begins at (or just after) t=0 — sample [0, 2×hp]
        tau_grid = np.linspace(0.0, 2.0 * hp, n_samples + 1)
        I        = waveform.current_at(tau_grid)
        dI       = np.diff(I)
        tau_mid  = 0.5 * (tau_grid[:-1] + tau_grid[1:])
        w        = np.where(dI < 0.0, -dI, 0.0)

        total = float(w.sum())
        if total <= 0.0:
            return 0.0, 0.0

        w_norm   = w / total
        tau_eff  = float(np.dot(tau_mid, w_norm))

        cumulative = np.cumsum(w_norm)
        idx        = int(np.searchsorted(cumulative, 0.99))
        tau_window = float(tau_mid[min(idx, len(tau_mid) - 1)])

        return tau_eff, tau_window

    else:
        # Turn-off completed before t=0 — find the last turn-off in [-2hp, 0]
        tau_grid = np.linspace(-2.0 * hp, 0.0, n_samples + 1)
        I        = waveform.current_at(tau_grid)
        dI       = np.diff(I)
        tau_mid  = 0.5 * (tau_grid[:-1] + tau_grid[1:])
        w        = np.where(dI < 0.0, -dI, 0.0)

        # Isolate the last contiguous turn-off (nearest to t=0)
        inc_idx  = np.where(dI > 0.0)[0]
        start    = int(inc_idx[-1]) + 1 if len(inc_idx) > 0 else 0
        w_last   = np.zeros_like(w)
        w_last[start:] = w[start:]

        total = float(w_last.sum())
        if total <= 0.0:
            return 0.0, 0.0

        w_norm  = w_last / total
        # centroid is at some tau < 0; tau_eff is the positive lag to t=0
        tau_eff = -float(np.dot(tau_mid, w_norm))
        return tau_eff, 0.0   # turn-off already complete: no gate rejection


def _apply_waveform_correction(
    dBdt: np.ndarray,
    t: np.ndarray,
    waveform,
    error: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
    r"""
    First-order waveform deconvolution for any transmitter waveform.

    For a transmitter with a finite turn-off, the measured :math:`\partial
    B_z/\partial t` is the convolution of the ideal step-off response
    :math:`b(t)` with the normalised turn-off impulse
    :math:`w(\tau)=-dI/d\tau`:

    .. math::

        V_\mathrm{meas}(t) = \int_0^{\tau_\mathrm{max}} w(\tau)\,b(t-\tau)\,d\tau

    A first-order Taylor expansion (Fitterman & Stewart 1986) yields

    .. math::

        b(t) \approx V_\mathrm{meas}(t + \tau_\mathrm{eff})
                     \left(1 + \frac{\tau_\mathrm{eff}}{2\,t}\right)

    where :math:`\tau_\mathrm{eff}` is the first moment of :math:`w`.
    Gates measured within the turn-off window (:math:`t_\mathrm{corr}
    \leq \tau_\mathrm{window}`) are discarded.

    Waveform-type dispatch:

    * :class:`~pycsamt.tdem.waveform.SquareWaveform` — no correction (ideal).
    * :class:`~pycsamt.tdem.waveform.RampWaveform` — analytic
      :math:`\tau_\mathrm{eff}=\tau_r/2`, recovering the classic
      Fitterman & Stewart formula exactly.
    * :class:`~pycsamt.tdem.waveform.HalfSineWaveform` and
      :class:`~pycsamt.tdem.waveform.CustomWaveform` — :math:`\tau_\mathrm{eff}`
      computed numerically via :func:`_waveform_moments`.

    Parameters
    ----------
    dBdt : np.ndarray
        Measured :math:`\partial B_z/\partial t` gate values.
    t : np.ndarray
        Gate centre times in seconds.
    waveform :
        Waveform descriptor or ``None``.
    error : np.ndarray or None
        Per-gate measurement uncertainty; corrected identically to ``dBdt``.

    Returns
    -------
    dBdt_corr : np.ndarray
    t_corr : np.ndarray
    error_corr : np.ndarray or None
    """
    tau_eff, tau_window = _waveform_moments(waveform)
    if tau_eff <= 0.0:
        return dBdt, t, error

    t_corr   = t + tau_eff
    amp_corr = 1.0 + tau_eff / (2.0 * t_corr)
    valid    = t_corr > tau_window

    err_corr = error[valid] if error is not None else None
    return dBdt[valid] * amp_corr[valid], t_corr[valid], err_corr


# ---------------------------------------------------------------------------
# LateTimeTransform
# ---------------------------------------------------------------------------

_CONFIG_CENTRAL = "central"
_CONFIG_IN_LOOP = "in_loop"
_CONFIG_OFFSET  = "offset"


class LateTimeTransform:
    r"""
    Convert TEM soundings to frequency-domain apparent impedance using
    the late-time apparent-resistivity approximation.

    Three loop configurations are handled automatically based on the
    ``offset`` and ``loop_dims`` attributes of each
    :class:`~pycsamt.tdem.TEMSounding`:

    * **central / coincident loop** (``offset == 0``): standard
      Ward & Hohmann (1988) formula.
    * **in-loop off-centre receiver** (``0 < offset < inner_radius``):
      Biot-Savart geometric correction applied to the effective moment.
    * **separated / offset loop** (``offset ≥ inner_radius``): offset
      dipole formula with :math:`d^3` denominator.

    Parameters
    ----------
    freq_convention : str
        Time-to-pseudo-frequency mapping:
        ``"skin_depth"`` (default) or ``"diffusion"``.
    phase_mode : str
        MT phase estimate: ``"homogeneous"`` (default, 45°) or
        ``"weidelt"`` (dispersion relation).
    drop_nan : bool
        Remove gates with undefined ρ_a.  Default ``True``.
    loop_geometry_correction : bool
        Apply the Biot-Savart in-loop correction when the receiver is
        off-centre.  Set ``False`` to use the central-loop formula for
        all configurations (legacy behaviour).  Default ``True``.
    in_loop_n_iter : int
        Number of iterations for the time-dependent geometry correction
        (see :func:`_rho_a_in_loop`).  ``0`` uses the static Biot-Savart
        approximation only.  Default 3.
    waveform_correction : bool
        Apply first-order waveform deconvolution when
        ``sounding.waveform`` is set (Fitterman & Stewart 1986).
        Supports :class:`~pycsamt.tdem.waveform.RampWaveform`
        (analytic), :class:`~pycsamt.tdem.waveform.HalfSineWaveform`,
        and :class:`~pycsamt.tdem.waveform.CustomWaveform` (both
        numerical).  :class:`~pycsamt.tdem.waveform.SquareWaveform`
        is a no-op (ideal step-off).  Default ``True``.

    Examples
    --------
    Central-loop sounding:

    >>> import numpy as np
    >>> from pycsamt.tdem import TEMSounding
    >>> from pycsamt.tdem.transform import LateTimeTransform
    >>> t = np.logspace(-5, -2, 30)
    >>> M = 8.0 * 100.0 ** 2
    >>> MU0_loc = 4e-7 * np.pi
    >>> dBdt = M * MU0_loc**2.5 / (10*np.sqrt(np.pi)*100.**1.5*t**2.5)
    >>> snd = TEMSounding(t, dBdt, current=8.0, tx_area=100.**2)
    >>> tr = LateTimeTransform()
    >>> result = tr.transform(snd)
    >>> result["freq"].shape
    (30,)

    Offset-loop sounding:

    >>> snd_off = TEMSounding(t, dBdt, current=8.0, tx_area=100.**2,
    ...                       offset=500.0, loop_dims=(100.0,))
    >>> result_off = tr.transform(snd_off)
    """

    def __init__(
        self,
        freq_convention: str = "skin_depth",
        phase_mode: str = "homogeneous",
        drop_nan: bool = True,
        loop_geometry_correction: bool = True,
        in_loop_n_iter: int = 3,
        waveform_correction: bool = True,
    ) -> None:
        self.freq_convention = freq_convention
        self.phase_mode = phase_mode
        self.drop_nan = drop_nan
        self.loop_geometry_correction = loop_geometry_correction
        self.in_loop_n_iter = int(in_loop_n_iter)
        self.waveform_correction = bool(waveform_correction)

    # ── loop-config detection ────────────────────────────────────────────

    @staticmethod
    def _detect_config(sounding: TEMSounding) -> str:
        """Return ``'central'``, ``'in_loop'``, or ``'offset'``."""
        d = float(sounding.offset)
        if d <= 0.0:
            return _CONFIG_CENTRAL
        r_inner = _loop_inner_radius(sounding.loop_shape, sounding.loop_dims)
        return _CONFIG_IN_LOOP if d < r_inner else _CONFIG_OFFSET

    # ── apparent-resistivity dispatch ───────────────────────────────────

    def _compute_rho_a(
        self,
        sounding: TEMSounding,
        dBdt: np.ndarray,
        t: np.ndarray,
    ) -> np.ndarray:
        """Route to the correct ρ_a formula for the sounding geometry."""
        M   = sounding.moment
        cfg = self._detect_config(sounding)

        if cfg == _CONFIG_CENTRAL or not self.loop_geometry_correction:
            return _rho_a_late_time(dBdt, t, M)

        if cfg == _CONFIG_IN_LOOP:
            # Use rx_position if provided (2-D offset), else scalar offset along x
            if getattr(sounding, "rx_position", None) is not None:
                rx, ry = float(sounding.rx_position[0]), float(sounding.rx_position[1])
            else:
                rx, ry = float(sounding.offset), 0.0
            return _rho_a_in_loop(
                dBdt, t, M,
                sounding.loop_shape, sounding.loop_dims,
                rx, ry,
                n_iter=self.in_loop_n_iter,
            )

        # _CONFIG_OFFSET
        return _rho_a_offset_loop(dBdt, t, M, sounding.offset)

    # ── public transform ─────────────────────────────────────────────────

    def transform(self, sounding: TEMSounding) -> Dict:
        r"""
        Transform one :class:`~pycsamt.tdem.TEMSounding` to
        frequency-domain arrays.

        Parameters
        ----------
        sounding : TEMSounding

        Returns
        -------
        dict with keys:
            ``freq``, ``Z``, ``Z_err``, ``rho_a``, ``phase_xy``,
            ``station_name``, ``x``, ``y``, ``elevation``,
            ``loop_config``.
        """
        dBdt  = sounding.dBdt()
        t     = sounding.time_gates.copy()
        error = sounding.error.copy() if sounding.error is not None else None

        if self.waveform_correction:
            dBdt, t, error = _apply_waveform_correction(
                dBdt, t, sounding.waveform, error=error
            )
            if len(t) < 2:
                raise ValueError(
                    "After waveform correction fewer than 2 gates remain. "
                    "Reduce ramp_off or supply more early-time gates."
                )

        rho_a = self._compute_rho_a(sounding, dBdt, t)
        freq  = _pseudo_freq(t, self.freq_convention)

        if self.drop_nan:
            mask = np.isfinite(rho_a) & np.isfinite(freq) & (freq > 0.0)
            rho_a  = rho_a[mask]
            freq   = freq[mask]
            t      = t[mask]
            dBdt   = dBdt[mask]
            err_raw = error[mask] if error is not None else None
        else:
            err_raw = error

        # sort ascending frequency
        order  = np.argsort(freq)
        freq   = freq[order]
        rho_a  = rho_a[order]
        if err_raw is not None:
            err_raw = err_raw[order]

        # ρ_a error → Z_err
        if err_raw is not None:
            rho_err = (2.0 / 3.0) * rho_a * np.abs(
                err_raw / (dBdt[order] + 1e-300)
            )
        else:
            rho_err = None

        phase_xy = _phase_from_rho(rho_a, freq, mode=self.phase_mode)
        Z_arr    = _build_z_array(rho_a, phase_xy, freq)

        Z_err = (
            _z_error_from_rho_error(rho_a, rho_err, freq)
            if rho_err is not None
            else np.full_like(Z_arr, np.nan, dtype=float)
        )

        return {
            "freq":         freq,
            "Z":            Z_arr,
            "Z_err":        Z_err,
            "rho_a":        rho_a,
            "phase_xy":     phase_xy,
            "station_name": sounding.station_name,
            "x":            sounding.x,
            "y":            sounding.y,
            "elevation":    sounding.elevation,
            "loop_config":  self._detect_config(sounding),
        }

    def transform_many(
        self,
        soundings: Sequence[TEMSounding],
    ) -> List[Dict]:
        """Transform a list of soundings.  Returns a list of result dicts."""
        return [self.transform(s) for s in soundings]


# ---------------------------------------------------------------------------
# Fourier / Kramers-Kronig numerical helpers
# ---------------------------------------------------------------------------

def _cosine_transform_1d(
    g: np.ndarray,
    t: np.ndarray,
    omega_arr: np.ndarray,
    n_interp: int = 0,
) -> np.ndarray:
    r"""
    Half-range Fourier cosine transform via trapezoidal quadrature in
    log-time space.

    .. math::

        F_c[g](\omega) = \int_0^\infty g(t)\,\cos(\omega t)\,\mathrm{d}t
        \approx \int_{\ln t_0}^{\ln t_N}
                g(t)\,\cos(\omega t)\,t\,\mathrm{d}(\ln t)

    Parameters
    ----------
    g : ndarray (n,)
        Integrand at discrete time gates.
    t : ndarray (n,)
        Strictly positive, monotone increasing time gates in seconds.
    omega_arr : ndarray (m,)
        Angular frequencies [rad/s] at which to evaluate.
    n_interp : int
        If > ``len(t)``, interpolate *g* onto a denser log-spaced grid
        before integration.  Default 0 (no interpolation).

    Returns
    -------
    ndarray (m,)
    """
    if n_interp > len(t):
        t_fine = np.exp(
            np.linspace(np.log(t[0]), np.log(t[-1]), int(n_interp))
        )
        log_abs_g = np.log(np.maximum(np.abs(g), 1e-300))
        g_fine = np.sign(np.interp(np.log(t_fine), np.log(t), g)) * np.exp(
            np.interp(np.log(t_fine), np.log(t), log_abs_g)
        )
    else:
        t_fine = t
        g_fine = g

    log_t  = np.log(t_fine)
    result = np.empty(len(omega_arr))
    for i, w in enumerate(omega_arr):
        integrand = g_fine * np.cos(w * t_fine) * t_fine   # ×t : d(log t)→dt
        result[i] = np.trapezoid(integrand, log_t)
    return result


def _kramers_kronig_re(
    omega: np.ndarray,
    im_k: np.ndarray,
) -> np.ndarray:
    r"""
    Compute :math:`\mathrm{Re}[K(\omega)]` from
    :math:`\mathrm{Im}[K(\omega)]` via the one-sided K-K relation:

    .. math::

        \mathrm{Re}[K(\omega)] =
            \frac{2}{\pi}\,\mathrm{P.V.}
            \int_0^\infty
            \frac{\omega'\,\mathrm{Im}[K(\omega')]}{\omega'^2 - \omega^2}
            \,\mathrm{d}\omega'

    The Cauchy principal value is handled by interpolating across the
    singular point :math:`\omega' = \omega` in log-space.

    Parameters
    ----------
    omega : ndarray (n,)   angular frequencies, log-spaced recommended
    im_k  : ndarray (n,)   imaginary part of the kernel

    Returns
    -------
    ndarray (n,)   real part of the kernel
    """
    n         = len(omega)
    log_omega = np.log(omega)
    num_base  = omega ** 2 * im_k    # ω'^2 Im[K] — numerator in log-ω integrand
    re_k      = np.zeros(n)

    for i in range(n):
        w2       = omega[i] ** 2
        den      = omega ** 2 - w2   # zero only at j == i
        tol      = 1e-10 * max(w2, 1e-100)
        mask_ok  = np.abs(den) > tol

        integrand           = np.where(mask_ok, num_base / np.where(mask_ok, den, 1.0), 0.0)
        bad                 = np.where(~mask_ok)[0]
        for j in bad:
            left  = j - 1 if j > 0     else j + 1
            right = j + 1 if j < n - 1 else j - 1
            integrand[j] = 0.5 * (integrand[left] + integrand[right])

        re_k[i] = (2.0 / np.pi) * np.trapezoid(integrand, log_omega)

    return re_k


# ---------------------------------------------------------------------------
# FourierTransform
# ---------------------------------------------------------------------------

class FourierTransform:
    r"""
    Rigorous TDEM → MT impedance via numerical Fourier cosine transform
    and Kramers-Kronig reconstruction (Meju 1996; Christensen 1990).

    The 1-D step-off forward relation gives:

    .. math::

        \frac{\partial B_z}{\partial t}(t) =
            \frac{2M}{\pi} \int_0^\infty
            \omega\,\mathrm{Im}[K(\omega)]\,\cos(\omega t)\,\mathrm{d}\omega

    Inverting by the half-range cosine transform:

    .. math::

        \mathrm{Im}[K(\omega)] =
            \frac{1}{M\omega}
            \int_0^\infty \frac{\partial B_z}{\partial t}(t)
            \cos(\omega t)\,\mathrm{d}t

    The real part :math:`\mathrm{Re}[K(\omega)]` is recovered via the
    Kramers-Kronig relation, then:

    .. math::

        Z(\omega) = i\omega\mu_0 K(\omega),\qquad
        \rho_a(\omega) = \omega\mu_0\,|K(\omega)|^2

    Parameters
    ----------
    n_freq : int
        Number of output frequency points.  Default 50.
    freq_min : float or None
        Minimum output frequency [Hz] (``None`` → auto from time gates).
    freq_max : float or None
        Maximum output frequency [Hz] (``None`` → auto from time gates).
    n_aux : int
        Auxiliary frequency points for the K-K integral (≥ 4 × n_freq).
        Default 200.
    n_interp : int
        Interpolation grid for the cosine transform (0 = disabled).
        Default 400.
    waveform_correction : bool
        Apply first-order waveform deconvolution (Fitterman & Stewart
        1986) when ``sounding.waveform`` is set.  Supports all waveform
        types via :func:`_apply_waveform_correction`.  Default ``True``.
    drop_nan : bool
        Remove output frequencies with undefined ρ_a.  Default ``True``.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.tdem import TEMSounding
    >>> from pycsamt.tdem.transform import FourierTransform
    >>> t = np.logspace(-5, -2, 40)
    >>> MU0_v = 4e-7 * np.pi; rho0 = 100.0; M = 8.0 * 1e4
    >>> dBdt = M * MU0_v**2.5 / (10*np.sqrt(np.pi)*rho0**1.5*t**2.5)
    >>> snd = TEMSounding(t, dBdt, current=8.0, tx_area=1e4)
    >>> res = FourierTransform(n_freq=30).transform(snd)
    >>> np.isfinite(res["rho_a"]).all()
    True

    References
    ----------
    .. [1] Meju, M.A. (1996). *Geophysics*, 61, 1579-1590.
    .. [2] Christensen, N.B. (1990). *Geophysical Prospecting*, 38, 545-568.
    .. [3] Ward & Hohmann (1988). Chapter 4.
    """

    def __init__(
        self,
        n_freq: int = 50,
        freq_min: Optional[float] = None,
        freq_max: Optional[float] = None,
        n_aux: int = 200,
        n_interp: int = 400,
        waveform_correction: bool = True,
        drop_nan: bool = True,
    ) -> None:
        self.n_freq              = int(n_freq)
        self.freq_min            = freq_min
        self.freq_max            = freq_max
        self.n_aux               = int(max(n_aux, 4 * n_freq))
        self.n_interp            = int(n_interp)
        self.waveform_correction = bool(waveform_correction)
        self.drop_nan            = bool(drop_nan)

    # ── frequency grids ──────────────────────────────────────────────────

    def _output_freqs(self, t: np.ndarray) -> np.ndarray:
        f_lo = self.freq_min if self.freq_min is not None else 1.0 / (2.0 * np.pi * t[-1])
        f_hi = self.freq_max if self.freq_max is not None else 1.0 / (2.0 * np.pi * t[0])
        f_lo = max(f_lo, 1e-8)
        return np.logspace(np.log10(f_lo), np.log10(f_hi), self.n_freq)

    def _aux_freqs(self, t: np.ndarray) -> np.ndarray:
        f_lo = max(0.5 / (2.0 * np.pi * t[-1]), 1e-8)
        f_hi = 2.0  / (2.0 * np.pi * t[0])
        return np.logspace(np.log10(f_lo), np.log10(f_hi), self.n_aux)

    # ── core kernel computation ──────────────────────────────────────────

    def _compute_kernel(
        self,
        dBdt: np.ndarray,
        t: np.ndarray,
        moment: float,
        omega_aux: np.ndarray,
        omega_out: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Return (re_k, im_k) at omega_out."""
        fc       = _cosine_transform_1d(dBdt / moment, t, omega_aux, self.n_interp)
        im_k_aux = np.where(omega_aux > 0.0, fc / omega_aux, 0.0)
        re_k_aux = _kramers_kronig_re(omega_aux, im_k_aux)

        log_aux = np.log(omega_aux)
        log_out = np.log(omega_out)
        im_k = np.interp(log_out, log_aux, im_k_aux)
        re_k = np.interp(log_out, log_aux, re_k_aux)
        return re_k, im_k

    # ── public API ───────────────────────────────────────────────────────

    def transform(self, sounding: TEMSounding) -> Dict:
        r"""
        Transform one :class:`~pycsamt.tdem.TEMSounding` to
        frequency-domain arrays.

        Returns the same dict structure as
        :meth:`LateTimeTransform.transform` plus key ``"method":
        "fourier"``.
        """
        dBdt = sounding.dBdt()
        t    = sounding.time_gates.copy()

        if len(t) < 4:
            raise ValueError(
                f"FourierTransform needs ≥ 4 time gates; "
                f"sounding '{sounding.station_name}' has {len(t)}."
            )

        if self.waveform_correction:
            dBdt, t, _ = _apply_waveform_correction(
                dBdt, t, sounding.waveform
            )
            if len(t) < 4:
                raise ValueError(
                    "After waveform correction fewer than 4 gates remain. "
                    "Reduce ramp_off or supply more early-time gates."
                )

        freqs_out = self._output_freqs(t)
        omega_out = 2.0 * np.pi * freqs_out
        omega_aux = 2.0 * np.pi * self._aux_freqs(t)

        re_k, im_k = self._compute_kernel(
            dBdt, t, sounding.moment, omega_aux, omega_out
        )

        # Z = iωμ₀K = -ωμ₀ Im[K] + i ωμ₀ Re[K]
        Z_xy    = -omega_out * MU0 * im_k + 1j * omega_out * MU0 * re_k
        rho_a   = np.abs(Z_xy) ** 2 / (omega_out * MU0)
        phase_xy = np.clip(np.angle(Z_xy, deg=True), 0.0, 90.0)

        # error propagation
        if sounding.error is not None:
            fc_err   = _cosine_transform_1d(
                sounding.error, sounding.time_gates, omega_aux, self.n_interp
            )
            ik_err_a = np.where(omega_aux > 0.0,
                                fc_err / (sounding.moment * omega_aux), 0.0)
            ik_err   = np.interp(np.log(omega_out), np.log(omega_aux), ik_err_a)
            rho_err  = 2.0 * np.abs(ik_err) * rho_a / np.maximum(np.abs(im_k), 1e-300)
        else:
            rho_err = None

        n     = len(freqs_out)
        Z_arr = np.zeros((n, 2, 2), dtype=complex)
        Z_arr[:, 0, 1] =  Z_xy
        Z_arr[:, 1, 0] = -Z_xy
        Z_err = (
            _z_error_from_rho_error(rho_a, rho_err, freqs_out)
            if rho_err is not None
            else np.full_like(Z_arr, np.nan, dtype=float)
        )

        if self.drop_nan:
            mask      = np.isfinite(rho_a) & (rho_a > 0.0) & np.isfinite(phase_xy)
            freqs_out = freqs_out[mask]
            rho_a     = rho_a[mask]
            phase_xy  = phase_xy[mask]
            Z_arr     = Z_arr[mask]
            Z_err     = Z_err[mask]

        order = np.argsort(freqs_out)
        return {
            "freq":         freqs_out[order],
            "Z":            Z_arr[order],
            "Z_err":        Z_err[order],
            "rho_a":        rho_a[order],
            "phase_xy":     phase_xy[order],
            "station_name": sounding.station_name,
            "x":            sounding.x,
            "y":            sounding.y,
            "elevation":    sounding.elevation,
            "loop_config":  "fourier",
            "method":       "fourier",
        }

    def transform_many(
        self,
        soundings: Sequence[TEMSounding],
    ) -> List[Dict]:
        """Transform a list of soundings.  Returns a list of result dicts."""
        return [self.transform(s) for s in soundings]


# ---------------------------------------------------------------------------
# TEMtoEDI — high-level dispatcher
# ---------------------------------------------------------------------------

class TEMtoEDI:
    r"""
    Convert one or more TEM soundings to a
    :class:`~pycsamt.seg.collection.EDICollection`.

    Wraps either :class:`LateTimeTransform` (default) or
    :class:`FourierTransform` and writes one synthetic EDI file per
    sounding into ``out_dir``.

    All three loop configurations (central, in-loop, offset) are handled
    automatically via :class:`LateTimeTransform`.

    Parameters
    ----------
    method : str
        ``"late_time"`` (default) or ``"fourier"``.
    freq_convention : str
        ``"skin_depth"`` (default) or ``"diffusion"``.
    phase_mode : str
        ``"homogeneous"`` (default) or ``"weidelt"``.
    loop_geometry_correction : bool
        Forward to :class:`LateTimeTransform`.  Default ``True``.
    out_dir : str or Path
        EDI output directory.  Default ``"edi_out/tem"``.
    verbose : int
        Verbosity.  Default 0.

    Examples
    --------
    Central-loop:

    >>> from pycsamt.tdem import TEMSounding, TEMtoEDI
    >>> import numpy as np
    >>> t = np.logspace(-5, -2, 25)
    >>> dBdt = 5e-5 * t ** (-5.0 / 2.0)
    >>> snd = TEMSounding(t, dBdt, current=8.0, tx_area=1e4,
    ...                   station_name="S01")
    >>> conv = TEMtoEDI(method="late_time")
    >>> coll = conv.transform(snd)

    Offset (separated) loop — just set ``offset``:

    >>> snd_off = TEMSounding(t, dBdt, current=8.0, tx_area=1e4,
    ...                       offset=500.0, loop_dims=(100.0,),
    ...                       station_name="S01")
    >>> coll_off = conv.transform(snd_off)
    """

    def __init__(
        self,
        method: str = "late_time",
        freq_convention: str = "skin_depth",
        phase_mode: str = "homogeneous",
        loop_geometry_correction: bool = True,
        out_dir: Union[str, "Path"] = "edi_out/tem",
        verbose: int = 0,
    ) -> None:
        self.method  = method.lower()
        self.out_dir = out_dir
        self.verbose = verbose

        if self.method == "late_time":
            self._transformer = LateTimeTransform(
                freq_convention=freq_convention,
                phase_mode=phase_mode,
                loop_geometry_correction=loop_geometry_correction,
            )
        elif self.method == "fourier":
            self._transformer = FourierTransform()
        else:
            raise ValueError(
                f"Unknown method '{method}'. Use 'late_time' or 'fourier'."
            )

    def transform(self, sounding: TEMSounding):
        result = self._transformer.transform(sounding)
        edi    = self._result_to_edifile(result)
        return self._make_collection([edi])

    def transform_many(self, soundings: Sequence[TEMSounding]):
        results = self._transformer.transform_many(soundings)
        edis    = [self._result_to_edifile(r) for r in results]
        return self._make_collection(edis)

    def save(self, sounding_or_soundings, path=None):
        """Transform and write EDI files.  Returns list of written paths."""
        from pathlib import Path as _Path

        out = _Path(self.out_dir if path is None else path)
        out.mkdir(parents=True, exist_ok=True)

        soundings = (
            [sounding_or_soundings]
            if isinstance(sounding_or_soundings, TEMSounding)
            else list(sounding_or_soundings)
        )

        written = []
        for snd in soundings:
            result = self._transformer.transform(snd)
            edi    = self._result_to_edifile(result)
            name   = (result["station_name"] or "tem_site") + ".edi"
            fp     = str(out / name)
            edi.write(new_edifn=name, savepath=str(out))
            written.append(fp)
            if self.verbose:
                cfg = result.get("loop_config", "?")
                print(f"  Wrote {fp}  ({result['freq'].size} freq, cfg={cfg})")

        return written

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _result_to_edifile(self, result: Dict):
        from ..seg.edi import EDIFile

        edi = EDIFile()
        edi.Z.freq  = result["freq"]
        edi.Z.z     = result["Z"]
        z_err = self._finite_z_error(result.get("Z_err"))
        if z_err is not None:
            edi.Z.z_err = z_err
        self._set_head(edi, result)
        return edi

    @staticmethod
    def _finite_z_error(z_err):
        """Return finite non-negative impedance errors when present."""
        if z_err is None:
            return None
        arr = np.asarray(z_err, dtype=float)
        if not np.isfinite(arr).any():
            return None
        arr = np.where(np.isfinite(arr), arr, 0.0)
        arr = np.maximum(arr, 0.0)
        return arr

    @staticmethod
    def _set_head(edi, result: Dict) -> None:
        cfg = result.get("loop_config", "central")
        try:
            from ..seg.heads import Head
            head = Head()
            head.dataid = result["station_name"] or "TEM_SITE"
            head.lon    = float(result.get("x", 0.0))
            head.lat    = float(result.get("y", 0.0))
            head.elev   = float(result.get("elevation", 0.0))
            head.acqby  = "pyCSAMT.tdem"
            head.fileby = "pyCSAMT.tdem"
            edi.add_section("head", head)
        except Exception:
            pass

        try:
            from ..seg.heads import Info
            info = Info()
            info.info_list = [
                "Data origin: TEM (time-domain EM) transformed to MT equivalent",
                "Transform: pycsamt.tdem.transform.TEMtoEDI",
                f"Loop configuration: {cfg}",
                "Z_xx = Z_yy = 0 (1-D assumption); Z_yx = -Z_xy",
            ]
            edi.add_section("info", info)
        except Exception:
            pass

    @staticmethod
    def _make_collection(edis: list):
        from ..seg.collection import EDICollection
        return EDICollection(items=edis)
