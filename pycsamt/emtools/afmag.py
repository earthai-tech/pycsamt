# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

r"""AFMAG-specific processing, diagnostics, and plotting.

AFMAG (Audio-Frequency Magnetics) is a passive, magnetic-field-only
method: there is no electric-field channel, so the usual impedance
tensor ``Z`` is not part of its data model. What AFMAG *does* measure —
whether the historical two-coil comparator or a modern tensor system
such as ZTEM/AirMt — is, in modern MT terms, exactly the same object
:class:`~pycsamt.site.base.Site` already carries as :attr:`Site.tipper`
(the complex ``(Tzx, Tzy)`` pair, ``Hz = Tzx*Hx + Tzy*Hy``); pyCSAMT's
own ZTEM adapter documents this equivalence explicitly (see
:mod:`pycsamt.airborne.ztem`). This module is therefore built, like
:mod:`pycsamt.emtools.tf`, entirely on :attr:`Site.tipper` — it never
touches :attr:`Site.z`.

Two genuinely different things live here side by side:

* Sections "Tilt-angle diagnostics" and "Motion-coupling QC" are
  ordinary ``emtools`` functions: they take ``sites`` (coerced through
  :func:`~pycsamt.emtools._core.ensure_sites`, like every other
  ``emtools`` function) and either return a
  :class:`~pycsamt.site.base.Sites` (the one mutating function,
  :func:`flag_motion_susceptible_band`) or a tidy
  :class:`pandas.DataFrame` (matching the convention already used by
  :func:`~pycsamt.emtools.spectra.psd_table`,
  :func:`~pycsamt.emtools.spectra.coherence_table`, and friends — not
  every ``emtools`` "processing" function returns ``Sites``; only the
  data-mutating ones do).
* Section "Motion-coupling physics" implements the rotation-matrix /
  Euler-angle motion-induced-noise method of Liu et al. (2018) directly
  from their Eq. 1-14. These functions are deliberately **not**
  ``Sites``-based: they operate on raw attitude (yaw/pitch/roll) time
  series and geomagnetic field geometry, neither of which
  :class:`~pycsamt.site.base.Site` carries (it is strictly
  frequency-domain, with no time axis and no attitude field — the only
  attitude container in pyCSAMT today is
  :class:`~pycsamt.airborne.navigation.NavigationTrack`, which belongs
  to the unrelated :mod:`pycsamt.airborne` line/dataset model). They
  are included here, as pure reusable physics, because the paper's
  actual deliverable is exactly this math; forcing it onto ``Sites``
  would require either fabricating time-series storage that does not
  exist or silently doing nothing. The one bridge between the two
  halves is :func:`motion_susceptibility_table` /
  :func:`flag_motion_susceptible_band`, which use the Section-1 math
  together with a nominal attitude-amplitude envelope to flag which
  *stations* are most exposed to motion noise in their low-frequency
  band — the paper's own finding is that this noise concentrates at
  low frequency and depends on attitude amplitude and geomagnetic
  geometry, not on a raw time series.

No apparent-resistivity/conductivity formula is derived from tilt angle
here: unlike full MT, classical AFMAG tilt data has no simple closed
form for that (nothing in Ward 1959 or Liu et al. 2018 provides one),
so none is invented.

References
----------
.. [Ward1959] Ward, S. H. (1959). AFMAG -- Airborne and ground.
   Geophysics, 24(4), 761-787.
.. [Liu2018] Liu, F., Huang, L., Pang, Y., Shi, Z., Xiao, P., & Fang, G.
   (2018). Airborne AFMAG method motion-induced noise simulation and
   suppression. Journal of Applied Geophysics, 158, 129-138.
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..api.labels import LOG10_PERIOD_LABEL
from ..api.plot import add_colorbar, add_polar_colorbar
from ..api.station import PYCSAMT_STATION_RENDERING
from ..api.style import PYCSAMT_STYLE
from ._core import (
    _apply_each,
    _axes_list,
    _get_t_block,
    _iter_items,
    _name,
    ensure_sites,
    hide_polar_radius_labels,
)

__all__ = [
    # motion-coupling physics (Liu et al. 2018)
    "euler_rotation_matrix",
    "geomagnetic_field_direction",
    "coil_normal_direction",
    "motion_coupling_cosine",
    "motion_coupling_angle",
    "simulate_motion_induced_voltage",
    "correct_motion_induced_noise",
    # tilt-angle diagnostics
    "afmag_tilt_angles",
    # motion-coupling QC
    "motion_susceptibility_table",
    "flag_motion_susceptible_band",
    # plots
    "plot_afmag_tilt_profile",
    "plot_afmag_tilt_psection",
    "plot_afmag_tilt_polar",
    "plot_motion_coupling",
    "plot_motion_susceptibility_map",
    "plot_afmag_correction_comparison",
]

# Historical original-AFMAG comparator band; reused as the default
# "motion-vulnerable" reference range rather than inventing a new one.
# Mirrors ``pycsamt.airborne.afmag.constants
# .AFMAG_ORIGINAL_TYPICAL_FREQUENCIES_HZ`` without importing
# ``pycsamt.airborne`` here (a heavier, unrelated package graph).
_DEFAULT_SUSCEPTIBLE_BAND_HZ = (150.0, 510.0)


# ─────────────────────────────────────────────────────────────────────────
# Motion-coupling physics (Liu et al. 2018, Eq. 1-14 and appendix)
# ─────────────────────────────────────────────────────────────────────────


def euler_rotation_matrix(
    yaw: Any,
    pitch: Any,
    roll: Any,
    *,
    degrees: bool = True,
) -> np.ndarray:
    r"""Return the ZYX "airline convention" attitude rotation matrix.

    Implements Liu et al. (2018) Eq. 3-4 (equivalently their appendix
    A1-A6): the standard aerospace Tait-Bryan rotation
    :math:`R_{LS} = R_z(\psi)\,R_y(\theta)\,R_x(\phi)` from the
    sensor/body frame (S) to the local, fixed east-north-up frame (L),
    where :math:`\psi` = yaw (about :math:`Z_S`), :math:`\theta` =
    pitch (about the once-rotated :math:`Y_S'`), :math:`\phi` = roll
    (about the twice-rotated :math:`X_S''`).

    Parameters
    ----------
    yaw, pitch, roll : array-like
        Attitude angles, broadcastable to a common shape ``(...,)``.
        Degrees by default; radians when ``degrees=False``.
    degrees : bool, default True
        Whether *yaw*, *pitch*, *roll* are in degrees.

    Returns
    -------
    ndarray of shape (..., 3, 3)
        Rotation matrix ``R_LS`` for each broadcast attitude sample;
        ``x_L = R_LS @ x_S``.

    Notes
    -----
    ``R_LS`` is a proper rotation: it is always orthogonal
    (``R_LS.T == inv(R_LS)``) with determinant ``+1``, for any input
    angles. Identity input (``yaw=pitch=roll=0``) returns the identity
    matrix.

    Examples
    --------
    >>> from pycsamt.emtools.afmag import euler_rotation_matrix
    >>> R = euler_rotation_matrix(0.0, 0.0, 0.0)
    >>> np.allclose(R, np.eye(3))
    True
    """
    yaw_arr = np.asarray(yaw, dtype=float)
    pitch_arr = np.asarray(pitch, dtype=float)
    roll_arr = np.asarray(roll, dtype=float)
    yaw_arr, pitch_arr, roll_arr = np.broadcast_arrays(
        yaw_arr, pitch_arr, roll_arr
    )
    if degrees:
        yaw_arr = np.deg2rad(yaw_arr)
        pitch_arr = np.deg2rad(pitch_arr)
        roll_arr = np.deg2rad(roll_arr)

    cy, sy = np.cos(yaw_arr), np.sin(yaw_arr)
    cp, sp = np.cos(pitch_arr), np.sin(pitch_arr)
    cr, sr = np.cos(roll_arr), np.sin(roll_arr)

    r = np.empty(yaw_arr.shape + (3, 3), dtype=float)
    r[..., 0, 0] = cy * cp
    r[..., 0, 1] = cy * sp * sr - sy * cr
    r[..., 0, 2] = cy * sp * cr + sy * sr
    r[..., 1, 0] = sy * cp
    r[..., 1, 1] = sy * sp * sr + cy * cr
    r[..., 1, 2] = sy * sp * cr - cy * sr
    r[..., 2, 0] = -sp
    r[..., 2, 1] = cp * sr
    r[..., 2, 2] = cp * cr
    return r


def geomagnetic_field_direction(
    inclination: Any,
    declination: Any,
    *,
    degrees: bool = True,
) -> np.ndarray:
    r"""Return the local geomagnetic field unit vector :math:`B_E`.

    In the (east, north, up) local frame used throughout this module
    (Liu et al. 2018 Eq. 8-9), with *inclination* the dip angle below
    horizontal (positive downward) and *declination* the horizontal
    bearing from geographic north toward east:

    .. math::

        \hat{B}_E = (\cos I \sin D,\ \cos I \cos D,\ -\sin I)

    Parameters
    ----------
    inclination, declination : array-like
        Geomagnetic inclination and declination, broadcastable to a
        common shape ``(...,)``. Degrees by default.
    degrees : bool, default True
        Whether *inclination* and *declination* are in degrees.

    Returns
    -------
    ndarray of shape (..., 3)
        Unit vector(s) ``(east, north, up)``.

    Examples
    --------
    >>> from pycsamt.emtools.afmag import geomagnetic_field_direction
    >>> np.round(geomagnetic_field_direction(90.0, 0.0), 6)  # straight down
    array([ 0.,  0., -1.])
    """
    inc = np.asarray(inclination, dtype=float)
    dec = np.asarray(declination, dtype=float)
    inc, dec = np.broadcast_arrays(inc, dec)
    if degrees:
        inc = np.deg2rad(inc)
        dec = np.deg2rad(dec)
    east = np.cos(inc) * np.sin(dec)
    north = np.cos(inc) * np.cos(dec)
    up = -np.sin(inc)
    return np.stack([east, north, up], axis=-1)


def coil_normal_direction(
    yaw: Any,
    pitch: Any,
    roll: Any,
    *,
    degrees: bool = True,
) -> np.ndarray:
    r"""Return the coil-normal unit vector :math:`N_L` in the local frame.

    Implements Liu et al. (2018) Eq. 10-12 for a single z-axis coil
    (:math:`N_S = (0, 0, 1)^T`, "we use z-axis hereinafter for the sake
    of simplification"): :math:`N_L = R_{LS}\,N_S`, i.e. the third
    column of :func:`euler_rotation_matrix`.

    Parameters
    ----------
    yaw, pitch, roll : array-like
        Forwarded to :func:`euler_rotation_matrix`.
    degrees : bool, default True
        Forwarded to :func:`euler_rotation_matrix`.

    Returns
    -------
    ndarray of shape (..., 3)
        Unit vector(s) ``(east, north, up)``. Identity attitude
        (``yaw=pitch=roll=0``) returns straight up, ``(0, 0, 1)``.
    """
    r = euler_rotation_matrix(yaw, pitch, roll, degrees=degrees)
    return r[..., :, 2]


def motion_coupling_cosine(
    yaw: Any,
    pitch: Any,
    roll: Any,
    inclination: Any,
    declination: Any,
    *,
    degrees: bool = True,
) -> np.ndarray:
    r"""Return :math:`\cos\theta(t)`, the coil/field coupling factor.

    Implements Liu et al. (2018) Eq. 13-14:
    :math:`\cos\theta(t) = \hat{B}_E \cdot N_L(t)`, the cosine of the
    angle between the (fixed) geomagnetic field direction and the
    (time-varying) coil-normal direction. This is the quantity whose
    time derivative drives the motion-induced noise voltage; see
    :func:`simulate_motion_induced_voltage`.

    Parameters
    ----------
    yaw, pitch, roll : array-like
        Attitude time series, broadcastable to a common shape.
    inclination, declination : array-like
        Geomagnetic field geometry, broadcastable against the attitude
        shape (a single site typically supplies scalars here while
        attitude varies over time).
    degrees : bool, default True
        Whether all five angle inputs are in degrees.

    Returns
    -------
    ndarray
        :math:`\cos\theta(t)`, broadcast shape of the inputs, values
        in ``[-1, 1]``.

    Examples
    --------
    At identity attitude the coil normal points straight up; at 90
    degrees inclination (a magnetic pole) the field points straight
    down, so the two are exactly anti-aligned and
    :math:`\cos\theta = -1`:

    >>> from pycsamt.emtools.afmag import motion_coupling_cosine
    >>> motion_coupling_cosine(0.0, 0.0, 0.0, 90.0, 0.0)
    -1.0
    """
    n_l = coil_normal_direction(yaw, pitch, roll, degrees=degrees)
    b_e = geomagnetic_field_direction(
        inclination, declination, degrees=degrees
    )
    n_l, b_e = np.broadcast_arrays(n_l, b_e)
    cos_theta = np.sum(b_e * n_l, axis=-1)
    return np.clip(cos_theta, -1.0, 1.0)


def motion_coupling_angle(
    yaw: Any,
    pitch: Any,
    roll: Any,
    inclination: Any,
    declination: Any,
    *,
    degrees: bool = True,
) -> np.ndarray:
    r"""Return :math:`\theta(t)` in degrees.

    See :func:`motion_coupling_cosine` for the underlying computation.
    """
    cos_theta = motion_coupling_cosine(
        yaw, pitch, roll, inclination, declination, degrees=degrees
    )
    return np.rad2deg(np.arccos(cos_theta))


def simulate_motion_induced_voltage(
    cos_theta: Any,
    *,
    dt: float,
    gain: float = 1.0,
    axis: int = -1,
) -> np.ndarray:
    r"""Return the simulated motion-induced noise voltage :math:`V(t)`.

    Implements Liu et al. (2018) Eq. 2, :math:`V(t) \propto
    -\,d(\cos\theta(t))/dt`, via a centred finite difference
    (:func:`numpy.gradient`).

    Parameters
    ----------
    cos_theta : array-like
        :math:`\cos\theta(t)`, typically from
        :func:`motion_coupling_cosine`.
    dt : float
        Sample interval in seconds along *axis*.
    gain : float, default 1.0
        Bundles the instrument-specific scale factor :math:`S \cdot N
        \cdot |B_E|` (coil area x turns x field magnitude) from Eq. 2.
        The paper calibrates this from the instrument rather than
        deriving it analytically, so it is left as an explicit,
        user-supplied factor rather than a fabricated constant.
    axis : int, default -1
        Time axis of *cos_theta*.

    Returns
    -------
    ndarray
        Simulated noise voltage, same shape as *cos_theta*.
    """
    arr = np.asarray(cos_theta, dtype=float)
    if dt <= 0.0 or not np.isfinite(dt):
        raise ValueError("dt must be finite and positive")
    return -float(gain) * np.gradient(arr, dt, axis=axis)


def correct_motion_induced_noise(
    measured: Any,
    predicted_noise: Any,
) -> np.ndarray:
    """Return the motion-noise-corrected signal (Liu et al. 2018 Sec. 5).

    The corrected signal is the measured signal minus the predicted
    noise voltage from :func:`simulate_motion_induced_voltage`. This
    function exists to give that one-line step an obvious, documented
    name rather than leaving every caller to subtract the two arrays
    themselves.

    Parameters
    ----------
    measured : array-like
        Raw movement-system signal.
    predicted_noise : array-like
        Predicted motion-induced noise, same shape as *measured*.

    Returns
    -------
    ndarray
        ``measured - predicted_noise``.

    Raises
    ------
    ValueError
        If *measured* and *predicted_noise* have different shapes.
    """
    m = np.asarray(measured, dtype=float)
    n = np.asarray(predicted_noise, dtype=float)
    if m.shape != n.shape:
        raise ValueError(
            "measured and predicted_noise must have the same shape: "
            f"{m.shape} != {n.shape}"
        )
    return m - n


# ─────────────────────────────────────────────────────────────────────────
# Tilt-angle diagnostics (classical AFMAG readout <-> Site.tipper)
# ─────────────────────────────────────────────────────────────────────────


def _tilt_components(t: np.ndarray) -> dict[str, np.ndarray]:
    """Return real/imag magnitude, azimuth, and resultant for one tipper."""
    tx, ty = t[:, 0], t[:, 1]
    out = {}
    for label, tx_c, ty_c in (
        ("real", np.real(tx), np.real(ty)),
        ("imag", np.imag(tx), np.imag(ty)),
    ):
        out[f"tilt_{label}_deg"] = np.rad2deg(
            np.arctan(np.hypot(tx_c, ty_c))
        )
        out[f"tilt_{label}_azimuth_deg"] = np.rad2deg(
            np.arctan2(ty_c, tx_c)
        )
    out["tilt_resultant_deg"] = np.rad2deg(
        np.arctan(np.hypot(np.abs(tx), np.abs(ty)))
    )
    return out


def afmag_tilt_angles(
    sites: Any,
    *,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    r"""Return per-station, per-frequency AFMAG tilt angles as a table.

    The classical AFMAG comparator readout is, in modern tipper terms,
    :math:`\arctan|T|` decomposed into in-phase ("real") and
    quadrature ("imag") components (Ward 1959) — the same real/imag
    split :mod:`pycsamt.emtools.tf`'s induction-arrow functions already
    use (:func:`~pycsamt.emtools.tf.plot_induction_arrows`), computed
    directly from :attr:`~pycsamt.site.base.Site.tipper`. No electric
    field or apparent-resistivity quantity is involved or derived.

    Parameters
    ----------
    sites : Sites-like
        Anything accepted by
        :func:`~pycsamt.emtools._core.ensure_sites`.
    recursive, on_dup, strict, verbose
        Forwarded to :func:`~pycsamt.emtools._core.ensure_sites`.

    Returns
    -------
    pandas.DataFrame
        Columns: ``station``, ``freq``, ``period``,
        ``tilt_real_deg``, ``tilt_real_azimuth_deg``,
        ``tilt_imag_deg``, ``tilt_imag_azimuth_deg``,
        ``tilt_resultant_deg``. Stations with no tipper are omitted,
        not filled with a fabricated zero.
    """
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    rows: list[dict[str, Any]] = []
    for i, ed in enumerate(_iter_items(S)):
        station = _name(ed, i)
        _, t, fr = _get_t_block(ed)
        if t is None or fr is None:
            continue
        comp = _tilt_components(t)
        for k in range(fr.size):
            row = {
                "station": station,
                "freq": float(fr[k]),
                "period": float(1.0 / fr[k]) if fr[k] != 0 else np.nan,
            }
            row.update({key: float(val[k]) for key, val in comp.items()})
            rows.append(row)
    return pd.DataFrame.from_records(rows)


# ─────────────────────────────────────────────────────────────────────────
# Motion-coupling QC (bridges the physics and tipper-diagnostic halves)
# ─────────────────────────────────────────────────────────────────────────


def _station_coordinates(ed: Any) -> tuple[float, float] | None:
    coords = getattr(ed, "coords", None)
    if coords is None:
        return None
    try:
        lat, lon = float(coords[0]), float(coords[1])
    except (TypeError, ValueError, IndexError):
        return None
    if np.isfinite(lat) and np.isfinite(lon):
        return lat, lon
    return None


def _susceptibility_rows(
    S: Any,
    *,
    inclination: float,
    declination: float,
    roll_amplitude_deg: float,
    pitch_amplitude_deg: float,
    yaw_amplitude_deg: float,
    n_phase: int = 180,
) -> list[dict[str, Any]]:
    """Sweep a nominal attitude envelope and score the resulting swing."""
    phase = np.linspace(0.0, 2.0 * np.pi, int(n_phase), endpoint=False)
    yaw = yaw_amplitude_deg * np.sin(phase)
    pitch = pitch_amplitude_deg * np.sin(phase)
    roll = roll_amplitude_deg * np.cos(phase)
    cos_theta = motion_coupling_cosine(
        yaw, pitch, roll, inclination, declination
    )
    swing = float(np.max(cos_theta) - np.min(cos_theta))

    rows: list[dict[str, Any]] = []
    for i, ed in enumerate(_iter_items(S)):
        rows.append(
            {
                "station": _name(ed, i),
                "inclination_deg": float(inclination),
                "declination_deg": float(declination),
                "cos_theta_min": float(np.min(cos_theta)),
                "cos_theta_max": float(np.max(cos_theta)),
                "susceptibility_score": swing,
            }
        )
    return rows


def motion_susceptibility_table(
    sites: Any,
    *,
    inclination: float,
    declination: float,
    roll_amplitude_deg: float,
    pitch_amplitude_deg: float,
    yaw_amplitude_deg: float = 0.0,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    r"""Score each station's exposure to motion-induced noise.

    Sweeps a nominal sinusoidal attitude envelope (amplitude only, no
    assumed platform frequency) through
    :func:`motion_coupling_cosine` at the survey's geomagnetic
    geometry, and reports the resulting :math:`\cos\theta` swing —
    the geometric factor that Liu et al. (2018) show drives
    motion-induced noise. This is independent of EM frequency; the
    paper's finding that the noise concentrates at low EM frequency
    comes from the platform's own physical oscillation spectrum
    (typically a few Hz), not from a frequency-dependence of this
    coupling geometry itself.

    Parameters
    ----------
    sites : Sites-like
        Anything accepted by
        :func:`~pycsamt.emtools._core.ensure_sites`.
    inclination, declination : float
        Geomagnetic inclination/declination in degrees, assumed
        common to the survey (a local, small-region approximation
        matching the paper's own "geomagnetic field is a constant
        vector in a local region" assumption).
    roll_amplitude_deg, pitch_amplitude_deg : float
        Nominal peak platform roll/pitch amplitude in degrees.
    yaw_amplitude_deg : float, default 0.0
        Nominal peak yaw amplitude; the paper finds yaw has no effect
        for a z-axis coil (rotation about the coil's own normal), so
        this defaults to zero rather than an assumed value.
    recursive, on_dup, strict, verbose
        Forwarded to :func:`~pycsamt.emtools._core.ensure_sites`.

    Returns
    -------
    pandas.DataFrame
        Columns: ``station``, ``inclination_deg``, ``declination_deg``,
        ``cos_theta_min``, ``cos_theta_max``, ``susceptibility_score``
        (``cos_theta_max - cos_theta_min``; larger means more exposed).
    """
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    rows = _susceptibility_rows(
        S,
        inclination=inclination,
        declination=declination,
        roll_amplitude_deg=roll_amplitude_deg,
        pitch_amplitude_deg=pitch_amplitude_deg,
        yaw_amplitude_deg=yaw_amplitude_deg,
    )
    return pd.DataFrame.from_records(rows)


def flag_motion_susceptible_band(
    sites: Any,
    *,
    inclination: float,
    declination: float,
    roll_amplitude_deg: float,
    pitch_amplitude_deg: float,
    yaw_amplitude_deg: float = 0.0,
    band_hz: tuple[float, float] = _DEFAULT_SUSCEPTIBLE_BAND_HZ,
    threshold: float = 0.05,
    action: str = "mask",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Any:
    r"""Flag/mask a station's low-frequency tipper band as motion-suspect.

    For stations whose :func:`motion_susceptibility_table` score
    exceeds *threshold*, the tipper values within *band_hz* are either
    masked (set to ``nan``, distinguishing "suspect" from "measured
    zero") or dropped entirely, mirroring the same
    :func:`~pycsamt.emtools._core.ensure_sites` /
    :func:`~pycsamt.emtools._core._apply_each` mutation contract used
    by :func:`~pycsamt.emtools.remove_noise.notch_powerline` and
    :func:`~pycsamt.emtools.remove_noise.drop_freqs_manual`.

    This does not attempt the paper's literal time-domain subtraction
    (see the module docstring for why: ``Sites`` carries no raw time
    series to subtract from). It is a QC gate, not a noise-removal
    algorithm — flagging a plausibly contaminated band so it can be
    reviewed or excluded, the same way the rest of ``emtools`` treats
    other unreliable-frequency conditions.

    Parameters
    ----------
    sites : Sites-like
        Anything accepted by
        :func:`~pycsamt.emtools._core.ensure_sites`.
    inclination, declination, roll_amplitude_deg, pitch_amplitude_deg,
    yaw_amplitude_deg
        Forwarded to :func:`motion_susceptibility_table`.
    band_hz : (float, float), default (150.0, 510.0)
        Frequency band to flag on a susceptible station, in Hz.
        Defaults to the historical original-AFMAG comparator band.
    threshold : float, default 0.05
        Minimum ``susceptibility_score`` (see
        :func:`motion_susceptibility_table`) for a station to be
        flagged at all.
    action : {"mask", "drop"}, default "mask"
        ``"mask"`` sets flagged tipper values to ``nan`` in place;
        ``"drop"`` removes the corresponding frequency rows entirely.
    inplace, recursive, on_dup, strict, verbose
        Standard ``emtools`` processing-function tail; see
        :func:`~pycsamt.emtools.remove_noise.notch_powerline` for the
        established convention this mirrors.

    Returns
    -------
    Sites
        The (optionally new) sites collection with flagged bands
        masked or dropped.

    Raises
    ------
    ValueError
        If *action* is not ``"mask"`` or ``"drop"``.
    """
    action = str(action).strip().lower()
    if action not in {"mask", "drop"}:
        raise ValueError("action must be 'mask' or 'drop'")
    lo, hi = float(band_hz[0]), float(band_hz[1])

    def _one(Si: Any) -> Any:
        rows = _susceptibility_rows(
            Si,
            inclination=inclination,
            declination=declination,
            roll_amplitude_deg=roll_amplitude_deg,
            pitch_amplitude_deg=pitch_amplitude_deg,
            yaw_amplitude_deg=yaw_amplitude_deg,
        )
        if not rows or rows[0]["susceptibility_score"] < threshold:
            return Si
        for ed in _iter_items(Si):
            T, t, fr = _get_t_block(ed)
            if T is None or t is None:
                continue
            in_band = (fr >= lo) & (fr <= hi)
            if not np.any(in_band):
                continue
            if action == "mask":
                t[in_band, :] = np.nan + 1j * np.nan
                if hasattr(T, "tipper"):
                    T.tipper = t
                elif hasattr(T, "T"):
                    T.T = t
            else:
                keep = ~in_band
                new_t = t[keep, :]
                new_fr = fr[keep]
                for name in ("tipper", "T"):
                    if hasattr(T, name):
                        setattr(T, name, new_t)
                        break
                if hasattr(T, "freq"):
                    T.freq = new_fr
                elif hasattr(ed, "freq"):
                    ed.freq = new_fr
        return Si

    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# ─────────────────────────────────────────────────────────────────────────
# Plots
# ─────────────────────────────────────────────────────────────────────────


def plot_afmag_tilt_profile(
    sites: Any,
    *,
    component: str = "real",
    frequency_hz: float | None = None,
    period_s: float | None = None,
    figsize: tuple[float, float] = (9.5, 4.0),
    station_label_step: int | None = 1,
    station_preset: str = "pseudosection",
    station_style: Any | None = None,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    r"""Plot the classic AFMAG flight-line tilt-angle profile.

    One value per station at a single reference frequency/period — the
    historical field presentation of AFMAG data, analogous to
    :func:`~pycsamt.emtools.remove_noise.plot_emap_filter_profile`'s
    station-profile layout.

    Parameters
    ----------
    sites : Sites-like
        Anything accepted by
        :func:`~pycsamt.emtools._core.ensure_sites`.
    component : {"real", "imag", "resultant"}, default "real"
        Which :func:`afmag_tilt_angles` column to plot.
    frequency_hz, period_s : float, optional
        Reference frequency/period; nearest available value is used
        per station. Exactly one may be given; the median frequency
        across all stations is used when neither is given.
    figsize : (float, float), default (9.5, 4.0)
        Used only when *ax* is not supplied.
    station_label_step, station_preset, station_style
        Forwarded to
        :meth:`~pycsamt.api.station.PyCSAMTStationRendering.style_for`
        via :func:`_apply_station_rendering`, matching the top-of-
        section station convention used throughout pyCSAMT.
    ax : matplotlib.axes.Axes, optional
        Existing axes to draw on.

    Returns
    -------
    matplotlib.axes.Axes
    """
    table = afmag_tilt_angles(sites)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    if table.empty:
        ax.text(0.5, 0.5, "no tipper data", ha="center", va="center")
        return ax

    col = {
        "real": "tilt_real_deg",
        "imag": "tilt_imag_deg",
        "resultant": "tilt_resultant_deg",
    }.get(component, "tilt_real_deg")

    if frequency_hz is None and period_s is not None:
        frequency_hz = 1.0 / max(float(period_s), 1e-24)
    if frequency_hz is None:
        frequency_hz = float(table["freq"].median())

    stations = list(dict.fromkeys(table["station"]))
    values = []
    used_freqs = []
    for station in stations:
        sub = table[table["station"] == station]
        idx = (sub["freq"] - frequency_hz).abs().idxmin()
        values.append(float(sub.loc[idx, col]))
        used_freqs.append(float(sub.loc[idx, "freq"]))

    x = np.arange(len(stations))
    _ml = PYCSAMT_STYLE.multiline
    ax.plot(x, values, lw=_ml.lw, alpha=_ml.alpha, marker="o")
    _apply_station_rendering(
        ax,
        stations,
        station_label_step=station_label_step,
        station_preset=station_preset,
        station_style=station_style,
    )
    ref = float(np.nanmedian(used_freqs))
    ax.set_ylabel(f"AFMAG tilt angle ({component}), deg")
    ax.set_title(f"AFMAG tilt profile at {ref:.4g} Hz", fontsize=10)
    ax.grid(True, ls=":", alpha=0.35)
    return ax


def plot_afmag_tilt_psection(
    sites: Any,
    *,
    component: str = "resultant",
    cmap: str = "RdYlBu_r",
    clim: tuple[float, float] | None = None,
    clim_pct: tuple[float, float] = (2.0, 98.0),
    figsize: tuple[float, float] = (9.0, 5.0),
    station_label_step: int | None = 1,
    station_preset: str = "pseudosection",
    station_style: Any | None = None,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    r"""Plot an AFMAG tilt-angle pseudosection (station x log-period).

    Parameters
    ----------
    sites : Sites-like
        Anything accepted by
        :func:`~pycsamt.emtools._core.ensure_sites`.
    component : {"real", "imag", "resultant"}, default "resultant"
        Which :func:`afmag_tilt_angles` column to image.
    cmap : str, default "RdYlBu_r"
        Colormap name.
    clim : (float, float), optional
        Explicit color limits; overrides *clim_pct*.
    clim_pct : (float, float), default (2.0, 98.0)
        Percentile color limits when *clim* is not given.
    figsize : (float, float), default (9.0, 5.0)
        Used only when *ax* is not supplied.
    station_label_step, station_preset, station_style
        See :func:`plot_afmag_tilt_profile`.
    ax : matplotlib.axes.Axes, optional
        Existing axes to draw on.

    Returns
    -------
    matplotlib.axes.Axes
    """
    table = afmag_tilt_angles(sites)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    if table.empty:
        ax.text(0.5, 0.5, "no tipper data", ha="center", va="center")
        return ax

    col = {
        "real": "tilt_real_deg",
        "imag": "tilt_imag_deg",
        "resultant": "tilt_resultant_deg",
    }.get(component, "tilt_resultant_deg")

    stations = list(dict.fromkeys(table["station"]))
    freqs = np.sort(table["freq"].unique())[::-1]
    grid = np.full((freqs.size, len(stations)), np.nan, dtype=float)
    for j, station in enumerate(stations):
        sub = table[table["station"] == station]
        for _, row in sub.iterrows():
            i = int(np.argmin(np.abs(freqs - row["freq"])))
            grid[i, j] = row[col]

    if clim is None:
        finite = grid[np.isfinite(grid)]
        vmin, vmax = (
            np.percentile(finite, clim_pct) if finite.size else (0.0, 1.0)
        )
    else:
        vmin, vmax = float(clim[0]), float(clim[1])

    extent = (
        -0.5,
        len(stations) - 0.5,
        np.log10(1.0 / freqs[0]),
        np.log10(1.0 / freqs[-1]),
    )
    im = ax.imshow(
        grid,
        aspect="auto",
        origin="upper",
        interpolation="nearest",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        extent=extent,
    )
    ax.set_ylabel(LOG10_PERIOD_LABEL)
    ax.set_title(f"AFMAG tilt pseudosection [{component}]", fontsize=10)
    _apply_station_rendering(
        ax,
        stations,
        station_label_step=station_label_step,
        station_preset=station_preset,
        station_style=station_style,
    )
    add_colorbar(im, ax, label="tilt angle (deg)")
    return ax


def plot_afmag_tilt_polar(
    sites: Any,
    *,
    station: str | None = None,
    component: str = "real",
    cmap: str = "viridis",
    figsize: tuple[float, float] = (5.5, 5.5),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    r"""Polar view of AFMAG tilt: azimuth and magnitude vs. period.

    Direct AFMAG-labelled counterpart of
    :func:`~pycsamt.emtools.tf.plot_tipper_polar`: each frequency is
    one scatter point, colour encodes :math:`\log_{10}(\text{period})`.

    Parameters
    ----------
    sites : Sites-like
    station : str, optional
        Station to plot; defaults to the first (sorted) station.
    component : {"real", "imag"}, default "real"
    cmap : str, default "viridis"
    figsize : (float, float), default (5.5, 5.5)
    recursive, on_dup, strict, verbose
        Forwarded to :func:`~pycsamt.emtools._core.ensure_sites`.
    ax : matplotlib.axes.Axes, optional
        Existing polar axes to draw on.

    Returns
    -------
    matplotlib.axes.Axes
    """
    import matplotlib.colors as mcolors

    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    pool = {_name(ed, i): ed for i, ed in enumerate(_iter_items(S))}
    if not pool:
        raise RuntimeError("no sites")
    st = station if station is not None else sorted(pool)[0]
    ed = pool.get(st)
    if ed is None:
        raise RuntimeError(f"station not found: {st!r}")

    _, t, fr = _get_t_block(ed)
    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, polar=True)
    if t is None or fr is None:
        hide_polar_radius_labels(ax)
        ax.set_title("no tipper")
        return ax

    tx, ty = t[:, 0], t[:, 1]
    u, v = (
        (np.real(tx), np.real(ty))
        if component == "real"
        else (np.imag(tx), np.imag(ty))
    )
    azimuth = np.arctan2(v, u)
    magnitude = np.rad2deg(np.arctan(np.hypot(u, v)))
    per = 1.0 / fr
    log_per = np.log10(np.maximum(per, 1e-9))
    norm = mcolors.Normalize(vmin=log_per.min(), vmax=log_per.max())

    sc = ax.scatter(
        azimuth,
        magnitude,
        c=log_per,
        cmap=cmap,
        norm=norm,
        s=30,
        edgecolors="none",
        zorder=4,
    )
    order = np.argsort(per)
    ax.plot(
        azimuth[order], magnitude[order], lw=1.0, alpha=0.4, color="0.6"
    )
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.grid(True, alpha=0.3)
    hide_polar_radius_labels(ax)
    ax.set_title(f"{st} — AFMAG tilt polar [{component}]", fontsize=10)
    cbar = add_polar_colorbar(sc, ax, label=r"$\log_{10}T$ (s)")
    cbar.ax.tick_params(labelsize=7)
    return ax


def plot_motion_coupling(
    yaw: Any,
    pitch: Any,
    roll: Any,
    inclination: float,
    declination: float,
    *,
    x: Any | None = None,
    xlabel: str = "sample",
    degrees: bool = True,
    ax: plt.Axes | None = None,
    figsize: tuple[float, float] = (8.0, 4.0),
) -> plt.Axes:
    r"""Plot :math:`\theta(t)` and :math:`\cos\theta(t)` for attitude data.

    Direct visualization of the paper's own Fig. 2-5 style curves,
    from raw attitude arrays -- not ``Sites``-based (see the module
    docstring for why).

    Parameters
    ----------
    yaw, pitch, roll : array-like
        Attitude time series (or a swept parameter), same shape.
    inclination, declination : float
        Geomagnetic field geometry in degrees.
    x : array-like, optional
        X-axis values (time or sample index); defaults to a plain
        integer index.
    xlabel : str, default "sample"
        X-axis label used when *x* is not supplied a label of its own.
    degrees : bool, default True
        Whether *yaw*/*pitch*/*roll* are in degrees.
    ax : matplotlib.axes.Axes, optional
        Existing axes to draw on; a right-hand twin axes is added for
        :math:`\cos\theta`.
    figsize : (float, float), default (8.0, 4.0)
        Used only when *ax* is not supplied.

    Returns
    -------
    matplotlib.axes.Axes
        The primary (:math:`\theta`) axes; the twin axes is reachable
        via ``ax.figure.axes[-1]``.
    """
    theta = motion_coupling_angle(
        yaw, pitch, roll, inclination, declination, degrees=degrees
    )
    cos_theta = motion_coupling_cosine(
        yaw, pitch, roll, inclination, declination, degrees=degrees
    )
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    xs = np.arange(theta.size) if x is None else np.asarray(x)

    ax.plot(xs, theta, color="tab:blue", label=r"$\theta(t)$")
    ax.set_ylabel(r"$\theta(t)$ (deg)", color="tab:blue")
    ax.tick_params(axis="y", labelcolor="tab:blue")
    ax.set_xlabel(xlabel if x is None else "")

    twin = ax.twinx()
    twin.plot(xs, cos_theta, color="tab:red", label=r"$\cos\theta(t)$")
    twin.set_ylabel(r"$\cos\theta(t)$", color="tab:red")
    twin.tick_params(axis="y", labelcolor="tab:red")

    ax.set_title("Motion-coupling angle (Liu et al. 2018)", fontsize=10)
    ax.grid(True, ls=":", alpha=0.35)
    return ax


def plot_motion_susceptibility_map(
    sites: Any,
    *,
    inclination: float,
    declination: float,
    roll_amplitude_deg: float,
    pitch_amplitude_deg: float,
    yaw_amplitude_deg: float = 0.0,
    cmap: str = "magma_r",
    figsize: tuple[float, float] = (7.0, 5.5),
    ax: plt.Axes | None = None,
) -> plt.Axes:
    r"""Map each station's motion-noise susceptibility score.

    Uses station coordinates when available
    (:attr:`~pycsamt.site.base.Site.coords`), otherwise falls back to
    a station-index profile -- the same fallback strategy
    :mod:`pycsamt.emtools.tf` uses for stations with no geo-referencing.

    Parameters
    ----------
    sites : Sites-like
    inclination, declination, roll_amplitude_deg, pitch_amplitude_deg,
    yaw_amplitude_deg
        Forwarded to :func:`motion_susceptibility_table`.
    cmap : str, default "magma_r"
    figsize : (float, float), default (7.0, 5.5)
        Used only when *ax* is not supplied.
    ax : matplotlib.axes.Axes, optional
        Existing axes to draw on.

    Returns
    -------
    matplotlib.axes.Axes
    """
    S = ensure_sites(sites, recursive=True, strict=False)
    table = motion_susceptibility_table(
        S,
        inclination=inclination,
        declination=declination,
        roll_amplitude_deg=roll_amplitude_deg,
        pitch_amplitude_deg=pitch_amplitude_deg,
        yaw_amplitude_deg=yaw_amplitude_deg,
    )
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    if table.empty:
        ax.text(0.5, 0.5, "no sites", ha="center", va="center")
        return ax

    xs, ys = [], []
    has_coords = True
    for i, ed in enumerate(_iter_items(S)):
        coord = _station_coordinates(ed)
        if coord is None:
            has_coords = False
            xs.append(float(i))
            ys.append(0.0)
        else:
            lat, lon = coord
            xs.append(lon)
            ys.append(lat)

    sc = ax.scatter(
        xs,
        ys,
        c=table["susceptibility_score"],
        cmap=cmap,
        s=80,
        edgecolors="k",
        linewidths=0.5,
    )
    ax.set_xlabel("Longitude" if has_coords else "Station index")
    ax.set_ylabel("Latitude" if has_coords else "")
    ax.set_title("AFMAG motion-noise susceptibility", fontsize=10)
    add_colorbar(sc, ax, label=r"$\cos\theta$ swing (susceptibility)")
    return ax


def plot_afmag_correction_comparison(
    before_sites: Any,
    after_sites: Any | None = None,
    *,
    component: str = "resultant",
    inclination: float | None = None,
    declination: float | None = None,
    roll_amplitude_deg: float | None = None,
    pitch_amplitude_deg: float | None = None,
    yaw_amplitude_deg: float = 0.0,
    cmap: str = "RdYlBu_r",
    delta_cmap: str = "RdBu_r",
    figsize: tuple[float, float] = (11.0, 8.2),
    axes: Any | None = None,
    **flag_kws: Any,
) -> plt.Figure:
    r"""Plot before/after/delta AFMAG tilt pseudosections.

    Structurally identical to
    :func:`~pycsamt.emtools.remove_noise.plot_emap_filter_psection`'s
    triptych: before, after, and :math:`\Delta` (after - before)
    panels. If *after_sites* is omitted, it is computed by calling
    :func:`flag_motion_susceptible_band` on *before_sites* with
    *inclination*, *declination*, *roll_amplitude_deg*, and
    *pitch_amplitude_deg* (all then required).

    Parameters
    ----------
    before_sites : Sites-like
    after_sites : Sites-like, optional
        Already-corrected sites; computed from *before_sites* when
        omitted.
    component : {"real", "imag", "resultant"}, default "resultant"
    inclination, declination, roll_amplitude_deg, pitch_amplitude_deg,
    yaw_amplitude_deg
        Forwarded to :func:`flag_motion_susceptible_band` when
        *after_sites* is omitted; *inclination*, *declination*,
        *roll_amplitude_deg*, and *pitch_amplitude_deg* are then
        required.
    cmap, delta_cmap : str
        Colormaps for the before/after and delta panels.
    figsize : (float, float), default (11.0, 8.2)
        Used only when *axes* is not supplied.
    axes : sequence of 3 Axes, optional
        Existing axes (before, after, delta) to draw on.
    **flag_kws
        Forwarded to :func:`flag_motion_susceptible_band` (for example
        ``band_hz``, ``threshold``, ``action``).

    Returns
    -------
    matplotlib.Figure

    Raises
    ------
    ValueError
        If *after_sites* is omitted and *inclination*, *declination*,
        *roll_amplitude_deg*, or *pitch_amplitude_deg* is not given.
    """
    if after_sites is None:
        required = {
            "inclination": inclination,
            "declination": declination,
            "roll_amplitude_deg": roll_amplitude_deg,
            "pitch_amplitude_deg": pitch_amplitude_deg,
        }
        missing = [k for k, v in required.items() if v is None]
        if missing:
            raise ValueError(
                "after_sites was not given, so "
                f"{', '.join(missing)} must be; got None"
            )
        after_sites = flag_motion_susceptible_band(
            before_sites,
            inclination=inclination,
            declination=declination,
            roll_amplitude_deg=roll_amplitude_deg,
            pitch_amplitude_deg=pitch_amplitude_deg,
            yaw_amplitude_deg=yaw_amplitude_deg,
            inplace=False,
            **flag_kws,
        )

    col = {
        "real": "tilt_real_deg",
        "imag": "tilt_imag_deg",
        "resultant": "tilt_resultant_deg",
    }.get(component, "tilt_resultant_deg")

    before = afmag_tilt_angles(before_sites)
    after = afmag_tilt_angles(after_sites)

    axes_given = _axes_list(axes, 3) if axes is not None else None
    if axes_given is None:
        fig, axes_arr = plt.subplots(
            3, 1, figsize=figsize, sharex=True,
            gridspec_kw={"hspace": 0.24},
        )
        axes_arr = np.asarray(axes_arr, dtype=object).ravel()
    else:
        axes_arr = np.asarray(axes_given, dtype=object)
        fig = axes_arr[0].figure

    if before.empty or after.empty:
        axes_arr[1].text(
            0.5, 0.5, "no paired data", ha="center", va="center"
        )
        return fig

    stations = [
        s for s in dict.fromkeys(before["station"])
        if s in set(after["station"])
    ]
    freqs = np.sort(before["freq"].unique())[::-1]

    def _grid(table: pd.DataFrame) -> np.ndarray:
        g = np.full((freqs.size, len(stations)), np.nan, dtype=float)
        for j, station in enumerate(stations):
            sub = table[table["station"] == station]
            for _, row in sub.iterrows():
                i = int(np.argmin(np.abs(freqs - row["freq"])))
                g[i, j] = row[col]
        return g

    before_grid = _grid(before)
    after_grid = _grid(after)
    delta_grid = after_grid - before_grid

    finite = np.concatenate([before_grid.ravel(), after_grid.ravel()])
    finite = finite[np.isfinite(finite)]
    vmin, vmax = (
        np.percentile(finite, (2.0, 98.0)) if finite.size else (0.0, 1.0)
    )
    dvals = np.abs(delta_grid[np.isfinite(delta_grid)])
    dvlim = max(
        float(np.nanpercentile(dvals, 95.0)) if dvals.size else 1e-6,
        1e-6,
    )

    extent = (
        -0.5,
        len(stations) - 0.5,
        np.log10(1.0 / freqs[0]),
        np.log10(1.0 / freqs[-1]),
    )
    panels = [
        (before_grid, axes_arr[0], cmap, vmin, vmax, "Before"),
        (after_grid, axes_arr[1], cmap, vmin, vmax, "After"),
        (delta_grid, axes_arr[2], delta_cmap, -dvlim, dvlim, r"$\Delta$"),
    ]
    images = []
    for data, ax, cm, lo, hi, title in panels:
        im = ax.imshow(
            data,
            aspect="auto",
            origin="upper",
            interpolation="nearest",
            cmap=cm,
            vmin=lo,
            vmax=hi,
            extent=extent,
        )
        images.append(im)
        ax.set_ylabel(LOG10_PERIOD_LABEL)
        ax.set_title(title, fontsize=9)

    _apply_station_rendering(axes_arr[0], stations)
    for ax in axes_arr[1:]:
        ax.tick_params(
            axis="x", labeltop=False, labelbottom=False, top=False,
            bottom=False,
        )
    cb_main = fig.colorbar(
        images[0], ax=[axes_arr[0], axes_arr[1]], fraction=0.018, pad=0.012,
    )
    cb_main.set_label(f"AFMAG tilt ({component}), deg")
    cb_delta = fig.colorbar(
        images[2], ax=axes_arr[2], fraction=0.018, pad=0.012,
    )
    cb_delta.set_label(r"$\Delta$ tilt (deg)")
    return fig


def _apply_station_rendering(
    ax: plt.Axes,
    labels: list[str],
    *,
    station_label_step: int | None = 1,
    station_preset: str = "pseudosection",
    station_style: Any | None = None,
) -> None:
    """Apply the shared top-of-section station rendering convention."""
    import copy

    style = station_style or PYCSAMT_STATION_RENDERING.style_for(
        station_preset
    )
    style = copy.copy(style)
    style.side = "top"
    style.max_labels = max(int(style.max_labels), len(labels))
    style.every = (
        1 if station_label_step is None else int(station_label_step)
    )
    x = np.arange(len(labels), dtype=float)
    style.apply(ax, x, labels, xlim=(-0.5, len(labels) - 0.5))
