# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
Coordinate frames of impedance data: sensor, geographic, profile.

Field EDIs are often in the **sensor frame** (x along the Ex dipole,
say 138 deg) with ``ZROT = 0``, while most programs read
``ZROT = 0`` as x = geographic north. This module rotates tensors,
errors, tippers *and* the sensor metadata consistently, so the
data frame is always

.. math::

    a_{frame} = a_{HX} + \mathrm{ZROT},

the HX-coil azimuth plus the recorded rotation (both clockwise from
north, degrees).

Rotation convention (as :meth:`pycsamt.z.z.Z.rotate`): rotating by
:math:`\theta` (clockwise positive) gives axes at
:math:`a + \theta`,

.. math::

    Z' = R Z R^{\mathsf T}, \qquad
    R = \begin{pmatrix}\cos\theta & \sin\theta \\
                       -\sin\theta & \cos\theta\end{pmatrix}.

Masked data
-----------
A rotation by an arbitrary angle mixes all four elements, so a
frequency with any missing (NaN) element becomes fully missing. For
multiples of 90 deg the rotation is an exact signed **permutation**
(e.g. :math:`\theta = -90^\circ` maps :math:`Z'_{xy} = -Z_{yx}`,
:math:`Z'_{yx} = -Z_{xy}`, :math:`Z'_{xx} = Z_{yy}`), applied by
indexing so that per-element masks follow their element. This is
the case needed to move profile-frame data into a strike frame
(x along strike) for 2-D codes such as Occam2D.
"""

from __future__ import annotations

from typing import Any

import numpy as np

__all__ = ["rotate_tensor", "rotate_meas", "frame_azimuth",
           "to_frame"]


def _rmat(theta: float) -> np.ndarray:
    t = np.radians(theta)
    return np.array([[np.cos(t), np.sin(t)], [-np.sin(t), np.cos(t)]])


def _quarter_turns(theta: float, tol: float = 1e-9) -> int | None:
    q = theta / 90.0
    k = round(q)
    return int(k) % 4 if abs(q - k) < tol else None


def rotate_tensor(
    z: np.ndarray,
    z_err: np.ndarray | None,
    theta: float,
) -> tuple[np.ndarray, np.ndarray | None]:
    """Rotate impedance tensors (and errors) by *theta* degrees.

    Parameters
    ----------
    z : ndarray of complex, shape (n, 2, 2)
        Impedances; NaN marks missing elements.
    z_err : ndarray, shape (n, 2, 2), optional
        Standard deviations of :math:`|Z_{ij}|`.
    theta : float
        Clockwise rotation (degrees); new axes at old + theta.

    Returns
    -------
    z_rot : ndarray of complex
    z_err_rot : ndarray or None
        For multiples of 90 deg the errors are permuted; otherwise
        they are propagated with
        :func:`pycsamt.utils.zmath.rotatematrix_incl_errors`
        (independent element errors).

    Notes
    -----
    See the module notes on masked data: multiples of 90 deg keep
    per-element masks, other angles mask whole frequencies that have
    any NaN element.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.emtools.frames import rotate_tensor
    >>> z = np.array([[[0, 1 + 1j], [np.nan, 0]]])
    >>> zr, _ = rotate_tensor(z, None, -90.0)
    >>> zr[0, 0, 1], zr[0, 1, 0]
    (np.complex128(nan+nanj), np.complex128(-1-1j))
    """
    z = np.asarray(z, dtype=complex)
    e = None if z_err is None else np.asarray(z_err, dtype=float)
    k = _quarter_turns(float(theta))
    if k is not None:
        r = np.rint(_rmat(90.0 * k)).astype(int)
        zr = np.empty_like(z)
        er = None if e is None else np.empty_like(e)
        src = {i: int(np.flatnonzero(r[i])[0]) for i in range(2)}
        for i in range(2):
            for j in range(2):
                a, b = src[i], src[j]
                sign = r[i, a] * r[j, b]
                zr[:, i, j] = sign * z[:, a, b]
                if er is not None:
                    er[:, i, j] = e[:, a, b]
        return zr, er

    from pycsamt.utils.zmath import rotatematrix_incl_errors

    zr = np.full_like(z, np.nan + 1j * np.nan)
    er = None if e is None else np.full_like(e, np.nan)
    for n in range(z.shape[0]):
        if not np.all(np.isfinite(z[n])):
            continue
        en = None
        if e is not None:
            en = np.where(np.isfinite(e[n]), e[n], 0.0)
        m, em = rotatematrix_incl_errors(z[n], float(theta), en)
        zr[n] = m
        if er is not None:
            er[n] = em
    return zr, er


def rotate_meas(definemeas: Any, theta: float) -> None:
    """Rotate sensor azimuths and dipole end-points in place.

    HMEAS azimuths increase by *theta*; EMEAS/HMEAS ``x``/``y``
    (north/east, m) are expressed in axes rotated by *theta*.
    Rotating by ``-a_HX`` thus describes the sensors as if laid
    north/east, which is what ``ZROT`` is relative to.
    """
    t = np.radians(-float(theta))
    c, s = np.cos(t), np.sin(t)

    def rot(x, y):
        return (round(x * c + y * s, 3), round(-x * s + y * c, 3))

    for h in getattr(definemeas, "hmeas", []) or []:
        if getattr(h, "azm", None) is not None:
            h.azm = round((float(h.azm) + float(theta)) % 360.0, 3)
        if getattr(h, "x", None) is not None:
            h.x, h.y = rot(float(h.x), float(h.y))
    for m in getattr(definemeas, "emeas", []) or []:
        if getattr(m, "x", None) is not None:
            m.x, m.y = rot(float(m.x), float(m.y))
            m.x2, m.y2 = rot(float(m.x2), float(m.y2))


def frame_azimuth(edi: Any) -> float:
    """Azimuth of the data x-axis: HX azimuth + ZROT (deg).

    Raises
    ------
    ValueError
        If ZROT varies with frequency (not a single frame).
    """
    dm = edi.get_section("definemeas")
    a_hx = 0.0
    for h in getattr(dm, "hmeas", []) or []:
        if str(getattr(h, "chtype", "")).upper() == "HX":
            a_hx = float(h.azm or 0.0)
            break
    zr = np.atleast_1d(np.asarray(
        getattr(edi.Z, "rotation_angle", 0.0), dtype=float))
    zr = zr[np.isfinite(zr)]
    if zr.size and np.ptp(zr) > 1e-6:
        raise ValueError("ZROT varies with frequency.")
    return float((a_hx + (zr[0] if zr.size else 0.0)) % 360.0)


def to_frame(edi: Any, azimuth: float) -> float:
    """Rotate an EDI (in place) so its data x-axis is at *azimuth*.

    Impedances, errors and tipper are rotated by
    ``azimuth - frame_azimuth(edi)``; the sensor metadata is rotated
    to north/east (HX at 0 deg) and ``ZROT`` set to *azimuth*, so
    ``frame_azimuth`` returns *azimuth* afterwards and any reader
    honouring ``ZROT`` sees a consistent file.

    Parameters
    ----------
    edi : EDIFile
        Station to update in place.
    azimuth : float
        Target azimuth of the x-axis (deg clockwise from north):
        0 for geographic, the line azimuth for a profile frame, the
        line azimuth - 90 for a strike frame perpendicular to it.

    Returns
    -------
    float
        Rotation applied (deg).
    """
    from pycsamt.z.z import Z

    a0 = frame_azimuth(edi)
    theta = float(azimuth) - a0
    theta = (theta + 180.0) % 360.0 - 180.0
    f = np.asarray(edi.Z.freq, float)
    zr, er = rotate_tensor(edi.Z.z, getattr(edi.Z, "z_err", None),
                           theta)
    edi.Z = Z(z_array=zr, z_err_array=er, freq=f)
    edi.Z.rotation_angle = np.full(f.size, float(azimuth) % 360.0)
    tip = getattr(getattr(edi, "Tip", None), "tipper", None)
    if tip is not None and np.asarray(tip).size and np.any(
            np.asarray(tip) != 0):
        r = _rmat(theta)
        t = np.asarray(tip, dtype=complex)
        edi.Tip.tipper = np.einsum("...ij,kj->...ik", t, r)
        edi.Tip.rotation_angle = np.full(f.size,
                                         float(azimuth) % 360.0)
    dm = edi.get_section("definemeas")
    if dm is not None:
        a_hx = next((float(h.azm) for h in dm.hmeas
                     if str(h.chtype).upper() == "HX"), 0.0)
        rotate_meas(dm, -a_hx)
    return theta
