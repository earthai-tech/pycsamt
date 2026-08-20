# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

r"""ZTEM-specific spatial-domain processing, diagnostics, and plotting.

ZTEM (Z-axis Tipper Electromagnetic; Lo and Zang 2008) is an airborne
AFMAG technique: a helicopter-towed vertical-field (:math:`H_z`) coil
is referenced to fixed-ground horizontal-field (:math:`H_x`,
:math:`H_y`) base-station coils, giving exactly the same
:math:`H_z = T_{zx}H_x + T_{zy}H_y` tipper relationship that
:class:`~pycsamt.site.base.Site` already carries as
:attr:`~pycsamt.site.base.Site.tipper` -- the equivalence
:mod:`pycsamt.emtools.afmag` documents and that
:mod:`pycsamt.airborne.ztem` encodes at the file-format/adapter layer.
This module is therefore built, like :mod:`pycsamt.emtools.afmag` and
:mod:`pycsamt.emtools.tf`, entirely on ``Site.tipper``.

It does not repeat what those two modules already do:

* :mod:`pycsamt.emtools.afmag` gives per-station, per-frequency tilt-
  angle diagnostics (:func:`~pycsamt.emtools.afmag.afmag_tilt_angles`)
  and motion-coupling QC -- purely *local* quantities, no derivative
  across stations.
* :mod:`pycsamt.emtools.tf` gives generic tipper visualization
  (induction arrows, polar plots, period sections) that applies
  equally to ground MT tipper and airborne ZTEM tipper alike.

What is missing from both, and what ZTEM survey practice specifically
relies on for a first-pass, pre-inversion read of a flight line
(Sattel and Witherly 2012), is *along-profile spatial-derivative*
processing -- comparing a station's tipper to its neighbours' rather
than looking at one station in isolation:

* :func:`total_divergence_table` -- the horizontal derivative of the
  in-line tipper component along the flight line. Lo and Zang (2008)
  define the (map-grid) "Total Divergence" as
  :math:`DT = \partial T_{zx}/\partial x + \partial T_{zy}/\partial y`
  (the literal equation is reproduced, e.g., as eq. 6 of
  wang2025). Sattel and Witherly (2012) note that "for the 2D case,
  the horizontal derivative is equivalent to the Peaker (Pedersen
  et al., 1994) and the total divergence" -- i.e. along one flight
  line the two published image products, DT and the VLF-style
  "Peaker", coincide with a single along-line derivative
  :math:`\partial T/\partial x`. Reproducing the full map-grid DT
  (needing a true 2-D grid across multiple flight lines, not a
  single :class:`~pycsamt.site.base.Sites` profile) is out of scope
  here; this function computes exactly the along-profile quantity
  the source papers show is equivalent for the 2D/profile case, and
  says so in its own docstring rather than silently overclaiming the
  3-D grid product.
* :func:`phase_rotate_table` -- the "phase-rotated response" image
  product Sattel and Witherly (2012, their Fig. 2) show converts a
  tipper crossover anomaly (over a lateral contact) into a peak
  anomaly, "obtained by reduction-to-pole (RTP) filtering and Hilbert
  transformation", noting "little difference between the RTP-filtered
  and the Hilbert-transformed profiles". This function implements the
  Hilbert-transform half of that pair directly
  (:func:`scipy.signal.hilbert`), a standard, unambiguous spatial
  analytic-signal transform.
* :func:`mask_outside_ztem_band` -- a ``Sites``-in/``Sites``-out QC
  gate (mirroring
  :func:`~pycsamt.emtools.afmag.flag_motion_susceptible_band`'s
  mask/drop contract) that reuses the survey's own published usable
  bandwidth, :attr:`pycsamt.airborne.ztem.ZTEMSystemSpec
  .practical_frequency_range_hz`, rather than inventing a new band
  definition -- the one function in this module meant to sit *inside*
  a processing pipeline rather than only produce a diagnostic table.

Karous-Hjelt (1983) pseudo-depth current-density sections and the
Becken and Pedersen (2003) tipper-gradient apparent-resistivity/phase
transform are both mentioned by Sattel and Witherly (2012) as further
image products, but neither closed-form is reproduced in the papers
available locally (``data/ZTEM/``), so neither is implemented here --
consistent with not silently converting an unverified formula into a
"working" function.

References
----------
.. [Lo2008] Lo, B., and Zang, M. (2008). Numerical modeling of Z-TEM
   (airborne AFMAG) responses to guide exploration strategies. SEG
   Expanded Abstracts, 27, 1098-1101.
.. [Legault2012] Legault, J. M., Zhao, S., and Fitch, R. (2012). ZTEM
   airborne AFMAG survey results over low sulphidation epithermal
   gold-silver vein systems at Gold Springs, south eastern Nevada.
   22nd International Geophysical Conference and Exhibition (ASEG),
   Brisbane.
.. [Sattel2012] Sattel, D., and Witherly, K. (2012). An overview of
   ZTEM data interpretation tools. 2012 NFEM Forum.
.. [Pedersen1994] Pedersen, L. B., Qian, W., Dynesius, L., and Zhang,
   P. (1994). An airborne tensor VLF system. From concept to
   realization. Geophysical Prospecting, 42, 863-883.
.. [wang2025] Wang, Y., Qu, J., Chen, T., Zhou, S., and Li, Y. (2025).
   Studies of three dimensional staggered-grid finite difference for
   Z-axis tipper electromagnetic numerical simulation. Frontiers in
   Earth Science, 13:1496312.
"""

from __future__ import annotations

import os.path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..api.station import PYCSAMT_STATION_RENDERING
from ..api.style import PYCSAMT_STYLE
from ._core import (
    _apply_each,
    _axes_list,
    _get_t_block,
    _iter_items,
    _name,
    _station_positions,
    _unwrap,
    ensure_any_sites,
)

__all__ = [
    # raw crossover diagnostics
    "ztem_crossover_diagnostics",
    # spatial-domain processing (tables)
    "total_divergence_table",
    "phase_rotate_table",
    # pipeline QC (Sites in, Sites out)
    "mask_outside_ztem_band",
    # plots
    "plot_ztem_tipper_profile",
    "plot_ztem_divergence_profile",
    "plot_ztem_divergence_psection",
    "plot_ztem_divergence_psection_grid",
    "plot_ztem_phase_rotation_profile",
    "plot_ztem_band_mask_psection",
    "plot_ztem_flight_lines",
    "plot_ztem_map",
]

_DIVERGENCE_COLS = [
    "station_a",
    "station_b",
    "x_m",
    "dx_m",
    "freq_hz",
    "period_s",
    "divergence_real",
    "divergence_imag",
    "divergence_abs",
]

_PHASE_ROTATE_COLS = [
    "x_m",
    "nearest_station",
    "freq_hz",
    "period_s",
    "raw",
    "rotated",
    "envelope",
]


# ─────────────────────────────────────────────────────────────────────────
# Private helpers
# ─────────────────────────────────────────────────────────────────────────


def _build_tipper_grid(
    sites: Any,
    spacing_m: float,
    recursive: bool,
    on_dup: str,
    strict: bool,
    verbose: int,
) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return station-ordered ``Tzx``/``Tzy`` grids for a ZTEM profile.

    Returns
    -------
    names : list of str
        Station names, ordered by position along the profile.
    x_pos : ndarray, shape (N,)
        Station chainage [m].
    freqs : ndarray, shape (F,)
        Sorted unique frequencies [Hz] across the profile.
    tzx_grid, tzy_grid : ndarray of complex, shape (N, F)
        ``NaN`` where a station has no value at that frequency.
    """
    S = ensure_any_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    items = list(_iter_items(S))
    if not items:
        empty = np.zeros((0, 0), dtype=complex)
        return [], np.array([]), np.array([]), empty, empty

    names: list[str] = []
    raw_items: list[Any] = []
    freq_t: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for i, ed in enumerate(items):
        T, t, fr = _get_t_block(ed)
        if T is None or t is None or fr is None or fr.size == 0:
            continue
        nm = _name(ed, i)
        names.append(nm)
        # _station_positions reads east/north (or lat/lon) straight off
        # each object; a Sites-wrapped Site hides those behind
        # ``.edi``, so unwrap here the same way
        # gradient_imaging._build_rho_grid does for the impedance case.
        raw_items.append(_unwrap(ed))
        freq_t[nm] = (fr, t)

    if not names:
        empty = np.zeros((0, 0), dtype=complex)
        return [], np.array([]), np.array([]), empty, empty

    x_pos = _station_positions(raw_items, spacing_m=spacing_m)
    order = np.argsort(x_pos)
    names = [names[k] for k in order]
    x_pos = x_pos[order]

    all_f: set = set()
    for fr, _ in freq_t.values():
        all_f.update(fr.tolist())
    freqs = np.array(sorted(all_f), dtype=float)
    if freqs.size == 0:
        empty = np.full((len(names), 0), np.nan, dtype=complex)
        return names, x_pos, freqs, empty, empty

    f_idx = {float(f): k for k, f in enumerate(freqs)}
    tzx_grid = np.full((len(names), freqs.size), np.nan, dtype=complex)
    tzy_grid = np.full((len(names), freqs.size), np.nan, dtype=complex)
    for i, nm in enumerate(names):
        fr, t = freq_t[nm]
        for j in range(fr.size):
            k = f_idx.get(float(fr[j]))
            if k is not None:
                tzx_grid[i, k] = t[j, 0]
                tzy_grid[i, k] = t[j, 1]
    return names, x_pos, freqs, tzx_grid, tzy_grid


def _resolve_target_frequency(
    freqs: np.ndarray,
    frequency_hz: float | None,
    period_s: float | None,
) -> tuple[int, float]:
    """Return ``(index, value)`` of the frequency nearest a target."""
    if frequency_hz is not None and period_s is not None:
        raise ValueError(
            "give at most one of frequency_hz or period_s, not both"
        )
    if frequency_hz is not None:
        target = float(frequency_hz)
    elif period_s is not None:
        target = 1.0 / max(float(period_s), 1e-24)
    else:
        target = float(np.median(freqs))
    k = int(np.argmin(np.abs(freqs - target)))
    return k, float(freqs[k])


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


# ─────────────────────────────────────────────────────────────────────────
# Total divergence / Peaker (Lo and Zang 2008; Sattel and Witherly 2012)
# ─────────────────────────────────────────────────────────────────────────


def total_divergence_table(
    sites: Any,
    spacing_m: float = 200.0,
    *,
    component: str = "tzx",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    r"""Along-profile ZTEM total-divergence / Peaker table.

    Computes the along-line horizontal derivative of the selected
    tipper component by first-order central-in-space finite
    differences between adjacent stations, at every frequency:

    .. math::

        DT(j,\,f) \approx
            \frac{T(j+1,\,f) - T(j,\,f)}{x(j+1) - x(j)}

    where stations are ordered by chainage along the profile (see
    :func:`~pycsamt.emtools._core._station_positions`). Per Sattel and
    Witherly (2012), this single along-line derivative *is* both the
    "Total Divergence" (Lo and Zang 2008) and the VLF-style "Peaker"
    (Pedersen et al. 1994) in the 2D/profile case -- the full 3-D
    map-grid divergence
    (:math:`\partial T_{zx}/\partial x + \partial T_{zy}/\partial y`)
    would additionally require a genuine cross-line (``y``) sampling
    that a single :class:`~pycsamt.site.base.Sites` profile does not
    carry, and is not attempted here.

    .. warning::

       *sites* is assumed to be **one flight line**. Chainage comes
       from :func:`~pycsamt.emtools._core._station_positions`, which
       projects every station onto a single bearing; passing a
       multi-line survey directly differentiates across line
       boundaries too, producing a physically meaningless value at
       every line-to-line join. Pre-filter to one line first (e.g.
       :meth:`~pycsamt.airborne.site.AirborneSites.select` on a
       per-line predicate) before calling this function on a
       multi-line dataset -- see
       :func:`~pycsamt.emtools.ztem.plot_ztem_map`'s own
       ``quantity="divergence"`` branch for a worked example that
       does this per-line grouping automatically.

    Parameters
    ----------
    sites : Sites-like or AirborneSites-like
        Anything accepted by
        :func:`~pycsamt.emtools._core.ensure_any_sites`: ground
        :class:`~pycsamt.site.base.Sites` when the data carries an
        impedance channel too, or
        :class:`~pycsamt.airborne.site.AirborneSites` for genuine
        tipper-only ZTEM/AFMAG data (a path/directory of EMTF-XML is
        routed automatically based on what it contains).
    spacing_m : float, default 200.0
        Fall-back inter-station spacing [m] used only when no station
        coordinates are available; see
        :func:`~pycsamt.emtools._core._station_positions`.
    component : {"tzx", "tzy"}, default "tzx"
        Tipper component to differentiate. ``"tzx"`` is the classical
        in-line (flight-direction) choice (Legault et al. 2012);
        ``"tzy"`` highlights structures striking across the line
        (Sattel and Witherly 2012).
    recursive, on_dup, strict, verbose
        Forwarded to :func:`~pycsamt.emtools._core.ensure_any_sites`.

    Returns
    -------
    pandas.DataFrame
        One row per (adjacent-station pair, frequency). Columns:
        ``station_a``, ``station_b`` (left/right station of the
        pair), ``x_m`` (pair midpoint chainage), ``dx_m`` (station
        spacing), ``freq_hz``, ``period_s``, ``divergence_real``,
        ``divergence_imag`` [each in units of tipper per metre],
        ``divergence_abs``. Pairs/frequencies with a missing tipper
        value on either side are omitted, not filled with zero.

    Raises
    ------
    ValueError
        If *component* is not ``"tzx"`` or ``"tzy"``.
    """
    component = str(component).strip().lower()
    if component not in ("tzx", "tzy"):
        raise ValueError(
            f"component must be 'tzx' or 'tzy'; got {component!r}"
        )

    names, x_pos, freqs, tzx_grid, tzy_grid = _build_tipper_grid(
        sites, spacing_m, recursive, on_dup, strict, verbose
    )
    grid = tzx_grid if component == "tzx" else tzy_grid

    n_st, n_f = grid.shape
    if n_st < 2 or n_f == 0:
        return pd.DataFrame(columns=_DIVERGENCE_COLS)

    dx = x_pos[1:] - x_pos[:-1]
    x_mid = 0.5 * (x_pos[1:] + x_pos[:-1])
    with np.errstate(invalid="ignore", divide="ignore"):
        ddx = (grid[1:, :] - grid[:-1, :]) / dx[:, None]

    rows: list[dict[str, Any]] = []
    for j in range(n_st - 1):
        for k in range(n_f):
            val = ddx[j, k]
            if not (np.isfinite(val.real) and np.isfinite(val.imag)):
                continue
            rows.append(
                {
                    "station_a": names[j],
                    "station_b": names[j + 1],
                    "x_m": float(x_mid[j]),
                    "dx_m": float(dx[j]),
                    "freq_hz": float(freqs[k]),
                    "period_s": 1.0 / max(float(freqs[k]), 1e-12),
                    "divergence_real": float(val.real),
                    "divergence_imag": float(val.imag),
                    "divergence_abs": float(np.abs(val)),
                }
            )
    if not rows:
        return pd.DataFrame(columns=_DIVERGENCE_COLS)
    return pd.DataFrame(rows, columns=_DIVERGENCE_COLS)


# ─────────────────────────────────────────────────────────────────────────
# Phase rotation (Sattel and Witherly 2012, Fig. 2: crossover -> peak)
# ─────────────────────────────────────────────────────────────────────────


def phase_rotate_table(
    sites: Any,
    *,
    frequency_hz: float | None = None,
    period_s: float | None = None,
    component: str = "tzx",
    part: str = "real",
    spacing_m: float = 200.0,
    n_resample: int | None = None,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    r"""Hilbert-transform "phase-rotated" ZTEM profile at one frequency.

    Reproduces the Hilbert-transform half of the phase-rotation image
    product described by Sattel and Witherly (2012, Fig. 2): the
    tipper component's along-profile crossover anomaly (odd about the
    causative contact) is converted into a peak anomaly (even about
    it) by taking the imaginary part of its spatial analytic signal
    (:func:`scipy.signal.hilbert`). Because the Hilbert transform
    assumes uniform sampling, the selected component is first linearly
    interpolated onto a uniform grid along chainage; the returned
    table is indexed by that uniform grid (with the nearest real
    station attached for reference), not by the original, generally
    unevenly spaced, station positions.

    .. warning::

       Like :func:`total_divergence_table`, *sites* is assumed to be
       one flight line: chainage is a single-bearing projection of
       every station, so a multi-line survey passed directly gets
       interpolated across line boundaries too. Pre-filter to one
       line first for a multi-line dataset.

    Parameters
    ----------
    sites : Sites-like or AirborneSites-like
        Anything accepted by
        :func:`~pycsamt.emtools._core.ensure_any_sites`: ground
        :class:`~pycsamt.site.base.Sites` when the data carries an
        impedance channel too, or
        :class:`~pycsamt.airborne.site.AirborneSites` for genuine
        tipper-only ZTEM/AFMAG data (a path/directory of EMTF-XML is
        routed automatically based on what it contains).
    frequency_hz, period_s : float, optional
        Target frequency/period; the nearest available frequency is
        used. At most one may be given; the median frequency across
        the profile is used when neither is given.
    component : {"tzx", "tzy"}, default "tzx"
    part : {"real", "imag"}, default "real"
        Which part of the complex tipper to phase-rotate. ``"real"``
        (in-phase) matches the classical VLF/ZTEM crossover-to-peak
        image product.
    spacing_m : float, default 200.0
        Fall-back inter-station spacing [m]; see
        :func:`~pycsamt.emtools._core._station_positions`.
    n_resample : int, optional
        Number of points on the uniform resampling grid. Defaults to
        the number of stations with a valid value (minimum 64).
    recursive, on_dup, strict, verbose
        Forwarded to :func:`~pycsamt.emtools._core.ensure_any_sites`.

    Returns
    -------
    pandas.DataFrame
        Columns: ``x_m`` (uniform-grid chainage), ``nearest_station``,
        ``freq_hz``, ``period_s``, ``raw`` (the interpolated, un-
        rotated component), ``rotated`` (its Hilbert transform),
        ``envelope`` (the analytic-signal magnitude
        :math:`\sqrt{raw^2 + rotated^2}`).

    Raises
    ------
    ValueError
        If *component* is not ``"tzx"``/``"tzy"``, if *part* is not
        ``"real"``/``"imag"``, or if both *frequency_hz* and
        *period_s* are given.
    """
    from scipy.signal import hilbert

    component = str(component).strip().lower()
    if component not in ("tzx", "tzy"):
        raise ValueError(
            f"component must be 'tzx' or 'tzy'; got {component!r}"
        )
    part = str(part).strip().lower()
    if part not in ("real", "imag"):
        raise ValueError(f"part must be 'real' or 'imag'; got {part!r}")

    names, x_pos, freqs, tzx_grid, tzy_grid = _build_tipper_grid(
        sites, spacing_m, recursive, on_dup, strict, verbose
    )
    grid = tzx_grid if component == "tzx" else tzy_grid
    if len(names) < 4 or freqs.size == 0:
        return pd.DataFrame(columns=_PHASE_ROTATE_COLS)

    k, f0 = _resolve_target_frequency(freqs, frequency_hz, period_s)
    col = grid[:, k]
    valid = np.isfinite(col.real) & np.isfinite(col.imag)
    if int(valid.sum()) < 4:
        return pd.DataFrame(columns=_PHASE_ROTATE_COLS)

    names_arr = np.array(names)
    xv = x_pos[valid]
    order = np.argsort(xv)
    xv = xv[order]
    yv = (col[valid].real if part == "real" else col[valid].imag)[order]
    nv = names_arr[valid][order]

    n = int(n_resample) if n_resample else max(len(xv), 64)
    x_uniform = np.linspace(float(xv[0]), float(xv[-1]), n)
    y_uniform = np.interp(x_uniform, xv, yv)

    analytic = hilbert(y_uniform)
    rotated = np.imag(analytic)
    envelope = np.abs(analytic)

    near = np.clip(np.searchsorted(xv, x_uniform), 1, xv.size - 1)
    prev = near - 1
    use_prev = np.abs(x_uniform - xv[prev]) <= np.abs(x_uniform - xv[near])
    near = np.where(use_prev, prev, near)
    nearest_station = nv[near]

    period_val = 1.0 / max(f0, 1e-12)
    return pd.DataFrame(
        {
            "x_m": x_uniform,
            "nearest_station": nearest_station,
            "freq_hz": np.full(n, f0),
            "period_s": np.full(n, period_val),
            "raw": y_uniform,
            "rotated": rotated,
            "envelope": envelope,
        },
        columns=_PHASE_ROTATE_COLS,
    )


# ─────────────────────────────────────────────────────────────────────────
# Pipeline QC: mask/drop frequencies outside the published ZTEM band
# ─────────────────────────────────────────────────────────────────────────


def mask_outside_ztem_band(
    sites: Any,
    *,
    band_hz: tuple[float, float] | None = None,
    system_spec: Any | None = None,
    action: str = "mask",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Any:
    r"""Mask/drop tipper frequencies outside the usable ZTEM band.

    Reuses the published usable bandwidth already carried by
    :class:`pycsamt.airborne.ztem.ZTEMSystemSpec` (default
    ``practical_frequency_range_hz`` of 22-720 Hz) rather than
    inventing a new band definition. For ground ``Sites`` input, this
    mirrors the same ``ensure_sites`` -> ``_apply_each`` mutation
    contract used by
    :func:`~pycsamt.emtools.afmag.flag_motion_susceptible_band` and
    :func:`~pycsamt.emtools.remove_noise.notch_powerline`: this is the
    one function in this module meant to sit inside a processing
    pipeline (container in, container out) rather than only produce a
    diagnostic table. For :class:`~pycsamt.airborne.site.AirborneSites`
    input, ``action="drop"`` is refused (see Raises) because it would
    leave the ``EMTF`` document's shared period axis inconsistent
    with the tipper transfer function's own periods; only
    ``action="mask"`` is offered there, matching
    :func:`~pycsamt.emtools.mobilemt.mask_outside_mobilemt_band`'s
    identical restriction for the identical reason.

    Parameters
    ----------
    sites : Sites-like or AirborneSites-like
        Anything accepted by
        :func:`~pycsamt.emtools._core.ensure_any_sites`: ground
        :class:`~pycsamt.site.base.Sites` when the data carries an
        impedance channel too, or
        :class:`~pycsamt.airborne.site.AirborneSites` for genuine
        tipper-only ZTEM/AFMAG data (a path/directory of EMTF-XML is
        routed automatically based on what it contains).
    band_hz : (float, float), optional
        Explicit ``(low, high)`` band in Hz. Mutually exclusive with
        *system_spec*; when neither is given, a default
        :class:`~pycsamt.airborne.ztem.ZTEMSystemSpec`'s
        ``practical_frequency_range_hz`` is used.
    system_spec : pycsamt.airborne.ztem.ZTEMSystemSpec, optional
        Survey-specific system metadata to read the band from.
    action : {"mask", "drop"}, default "mask"
        ``"mask"`` sets out-of-band tipper values to ``nan`` in place;
        ``"drop"`` removes the corresponding frequency rows entirely.
    inplace, recursive, on_dup, strict, verbose
        Standard ``emtools`` processing-function tail; see
        :func:`~pycsamt.emtools.afmag.flag_motion_susceptible_band`
        for the established convention this mirrors.

    Returns
    -------
    Sites
        The (optionally new) sites collection with out-of-band tipper
        frequencies masked or dropped.

    Raises
    ------
    ValueError
        If *action* is not ``"mask"`` or ``"drop"``; if both
        *band_hz* and *system_spec* are given; or if *action* is
        ``"drop"`` and *sites* resolves to
        :class:`~pycsamt.airborne.site.AirborneSites`.
    TypeError
        If *system_spec* is given and is not a
        :class:`~pycsamt.airborne.ztem.ZTEMSystemSpec`.
    """
    from ..airborne.site import AirborneSites
    from ..airborne.ztem import ZTEMSystemSpec

    action = str(action).strip().lower()
    if action not in {"mask", "drop"}:
        raise ValueError("action must be 'mask' or 'drop'")
    if band_hz is not None and system_spec is not None:
        raise ValueError(
            "give at most one of band_hz or system_spec, not both"
        )
    if system_spec is not None and not isinstance(
        system_spec, ZTEMSystemSpec
    ):
        raise TypeError("system_spec must be a ZTEMSystemSpec or None")

    if band_hz is not None:
        lo, hi = float(band_hz[0]), float(band_hz[1])
    else:
        spec = system_spec or ZTEMSystemSpec()
        lo, hi = spec.practical_frequency_range_hz

    def _one(Si: Any) -> Any:
        for ed in _iter_items(Si):
            T, t, fr = _get_t_block(ed)
            if T is None or t is None:
                continue
            keep = (fr >= lo) & (fr <= hi)
            if keep.all():
                continue
            if action == "mask":
                t[~keep, :] = np.nan + 1j * np.nan
                if hasattr(T, "tipper"):
                    T.tipper = t
                elif hasattr(T, "T"):
                    T.T = t
            else:
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

    S = ensure_any_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if isinstance(S, AirborneSites):
        if action == "drop":
            raise ValueError(
                "action='drop' is not supported for AirborneSites "
                "(it would leave the EMTF document's shared period "
                "axis inconsistent with the tipper transfer "
                "function's own periods); use action='mask' instead"
            )
        # AirborneSite/EMTF are plain deepcopy-able objects with no
        # EDICollection-style container to rebuild, so the ground
        # _apply_each machinery (which does exactly that rebuilding)
        # is neither needed nor applicable here.
        if not inplace:
            import copy

            S = copy.deepcopy(S)
        _one(S)
        return S
    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# ─────────────────────────────────────────────────────────────────────────
# Raw crossover diagnostics (Legault et al. 2012, Fig. 6)
# ─────────────────────────────────────────────────────────────────────────


def ztem_crossover_diagnostics(
    sites: Any,
    *,
    frequency_hz: float | None = None,
    period_s: float | None = None,
    component: str = "tzx",
    spacing_m: float = 200.0,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> dict[str, Any]:
    r"""Legault et al. (2012, Fig. 6)-style raw in-phase/quadrature crossover.

    Their own synthetic forward-model example over a mushroom-shaped
    epithermal target reads a negative-to-positive **in-phase**
    crossover directly above the target, generally accompanied by a
    (noisier) negative-to-positive **quadrature** crossover at every
    frequency -- the qualitative, single-frequency, pre-processing
    read of a flight line before any derivative or transform is
    applied (contrast :func:`total_divergence_table`/
    :func:`phase_rotate_table`, both of which operate on this same raw
    tipper but convert the crossover into a peak). This function finds
    those two crossovers and reports the peak-to-peak swing of each,
    the same crossover/amplitude measurements
    :func:`~pycsamt.emtools.afmag.original_afmag_conductor_diagnostics`
    reports for the AFMAG comparator, applied here to the real/
    imaginary parts of one tipper component instead of two hardware
    frequencies.

    .. warning::

       Like :func:`total_divergence_table`, *sites* is assumed to be
       one flight line; pre-filter a multi-line survey to one line
       first (see that function's own warning for why).

    Parameters
    ----------
    sites : Sites-like or AirborneSites-like
        Anything accepted by
        :func:`~pycsamt.emtools._core.ensure_any_sites`.
    frequency_hz, period_s : float, optional
        Reference frequency/period; nearest available value is used.
        At most one may be given; the median frequency across the
        profile is used when neither is given.
    component : {"tzx", "tzy"}, default "tzx"
        ``"tzx"`` is the classical in-line choice (Legault et al.
        2012); ``"tzy"`` highlights cross-line structure (Sattel and
        Witherly 2012).
    spacing_m : float, default 200.0
        Forwarded to :func:`~pycsamt.emtools._core._station_positions`.
    recursive, on_dup, strict, verbose
        Forwarded to :func:`~pycsamt.emtools._core.ensure_any_sites`.

    Returns
    -------
    dict
        Keys: ``freq_hz``, ``crossover_real_m``, ``crossover_imag_m``
        (along-profile position, ``nan`` if that part never changes
        sign between its own max and min), ``peak_to_peak_real``,
        ``peak_to_peak_imag``, and ``profile`` -- a
        :class:`pandas.DataFrame` with columns ``station``,
        ``position_m``, ``real``, ``imag`` (the raw, dimensionless
        tipper values; multiply by 100 for the percent convention
        Legault et al. 2012 and Sattel and Witherly 2012 both plot).

    Raises
    ------
    ValueError
        If *component* is not ``"tzx"``/``"tzy"``, or fewer than 2
        stations have a usable value at the resolved frequency.
    """
    component = str(component).strip().lower()
    if component not in ("tzx", "tzy"):
        raise ValueError(
            f"component must be 'tzx' or 'tzy'; got {component!r}"
        )

    names, x_pos, freqs, tzx_grid, tzy_grid = _build_tipper_grid(
        sites, spacing_m, recursive, on_dup, strict, verbose
    )
    grid = tzx_grid if component == "tzx" else tzy_grid
    if len(names) < 2 or freqs.size == 0:
        raise ValueError(
            "ztem_crossover_diagnostics requires at least 2 stations "
            "with tipper data"
        )

    k, f0 = _resolve_target_frequency(freqs, frequency_hz, period_s)
    col = grid[:, k]
    valid = np.isfinite(col.real) & np.isfinite(col.imag)
    if int(valid.sum()) < 2:
        raise ValueError(
            "ztem_crossover_diagnostics: fewer than 2 stations have a "
            f"usable value at {f0:.4g} Hz"
        )

    names_arr = np.array(names)[valid]
    xv = x_pos[valid]
    order = np.argsort(xv)
    xv = xv[order]
    names_arr = names_arr[order]
    re = col[valid].real[order]
    im = col[valid].imag[order]

    profile = pd.DataFrame(
        {
            "station": names_arr,
            "position_m": xv,
            "real": re,
            "imag": im,
        }
    )

    def _crossover(x: np.ndarray, y: np.ndarray) -> float:
        if x.size < 2:
            return float("nan")
        i_max, i_min = int(np.argmax(y)), int(np.argmin(y))
        lo, hi = (i_max, i_min) if i_max < i_min else (i_min, i_max)
        for i in range(lo, hi):
            y0, y1 = y[i], y[i + 1]
            if y0 == 0.0:
                return float(x[i])
            if np.sign(y0) != np.sign(y1):
                t = y0 / (y0 - y1)
                return float(x[i] + t * (x[i + 1] - x[i]))
        return float("nan")

    return {
        "freq_hz": f0,
        "crossover_real_m": _crossover(xv, re),
        "crossover_imag_m": _crossover(xv, im),
        "peak_to_peak_real": float(np.max(re) - np.min(re)),
        "peak_to_peak_imag": float(np.max(im) - np.min(im)),
        "profile": profile,
    }


# ─────────────────────────────────────────────────────────────────────────
# Plots
# ─────────────────────────────────────────────────────────────────────────


def plot_ztem_tipper_profile(
    sites: Any,
    *,
    frequency_hz: float | None = None,
    period_s: float | None = None,
    component: str = "tzx",
    as_percent: bool = True,
    figsize: tuple[float, float] = (9.5, 4.2),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    r"""Legault et al. (2012, Fig. 6)-style raw in-phase/quadrature profile.

    The classic ZTEM field-presentation figure: real (in-phase) and
    imaginary (quadrature) tipper plotted together at one frequency,
    in percent, along real flight-line chainage -- Fig. 6's own
    "METERS" x-axis, not a discrete station index -- with each part's
    crossover marked; see :func:`ztem_crossover_diagnostics`.

    Parameters
    ----------
    sites : Sites-like or AirborneSites-like
        Anything accepted by
        :func:`~pycsamt.emtools._core.ensure_any_sites`.
    frequency_hz, period_s, component
        Forwarded to :func:`ztem_crossover_diagnostics`.
    as_percent : bool, default True
        Multiply the plotted tipper values by 100, matching the
        percent convention both cited papers use. Set ``False`` to
        plot the raw dimensionless tipper instead.
    figsize : (float, float), default (9.5, 4.2)
        Used only when *ax* is not supplied.
    recursive, on_dup, strict, verbose
        Forwarded to :func:`ztem_crossover_diagnostics`.
    ax : matplotlib.axes.Axes, optional
        Existing axes to draw on.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    try:
        diag = ztem_crossover_diagnostics(
            sites,
            frequency_hz=frequency_hz,
            period_s=period_s,
            component=component,
            recursive=recursive,
            on_dup=on_dup,
            strict=strict,
            verbose=verbose,
        )
    except ValueError:
        ax.text(0.5, 0.5, "no tipper data", ha="center", va="center")
        return ax

    profile = diag["profile"]
    scale = 100.0 if as_percent else 1.0
    x = profile["position_m"].to_numpy()

    ax.axhline(0.0, color="0.8", lw=0.8)
    ax.plot(
        x, profile["real"] * scale, "-", color="tab:blue", marker="o",
        ms=4, label="in-phase",
    )
    ax.plot(
        x, profile["imag"] * scale, "--", color="tab:red", marker="x",
        ms=5, label="quadrature",
    )
    for cross, color in (
        (diag["crossover_real_m"], "tab:blue"),
        (diag["crossover_imag_m"], "tab:red"),
    ):
        if np.isfinite(cross):
            ax.axvline(cross, color=color, lw=1.0, ls=":", alpha=0.7)

    ax.set_xlabel("Position along flight line (m)")
    tlabel = r"$T_{zx}$" if component == "tzx" else r"$T_{zy}$"
    ax.set_ylabel(f"{tlabel} ({'%' if as_percent else 'dimensionless'})")
    ax.set_title(
        f"ZTEM raw in-phase/quadrature profile at {diag['freq_hz']:.4g} Hz",
        fontsize=10,
    )
    ax.legend(fontsize=8, loc="best")
    ax.grid(True, ls=":", alpha=0.35)
    return ax


def plot_ztem_divergence_profile(
    sites: Any,
    *,
    component: str = "tzx",
    part: str = "real",
    frequency_hz: float | None = None,
    period_s: float | None = None,
    spacing_m: float = 200.0,
    figsize: tuple[float, float] = (9.5, 4.0),
    station_label_step: int | None = 1,
    station_preset: str = "pseudosection",
    station_style: Any | None = None,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    r"""Plot the ZTEM total-divergence / Peaker flight-line profile.

    One value per adjacent-station pair at a single reference
    frequency/period -- the along-line, pre-gridding form of the
    "Total Divergence" / "Peaker" image product (see
    :func:`total_divergence_table`).

    Parameters
    ----------
    sites : Sites-like or AirborneSites-like
        Anything accepted by
        :func:`~pycsamt.emtools._core.ensure_any_sites`: ground
        :class:`~pycsamt.site.base.Sites` when the data carries an
        impedance channel too, or
        :class:`~pycsamt.airborne.site.AirborneSites` for genuine
        tipper-only ZTEM/AFMAG data (a path/directory of EMTF-XML is
        routed automatically based on what it contains).
    component : {"tzx", "tzy"}, default "tzx"
    part : {"real", "imag"}, default "real"
    frequency_hz, period_s : float, optional
        Reference frequency/period; nearest available value is used
        per station pair. At most one may be given; the median
        frequency across all pairs is used when neither is given.
    spacing_m : float, default 200.0
        Forwarded to :func:`total_divergence_table`.
    figsize : (float, float), default (9.5, 4.0)
        Used only when *ax* is not supplied.
    station_label_step, station_preset, station_style
        Forwarded to
        :meth:`~pycsamt.api.station.PyCSAMTStationRendering.style_for`
        via the shared top-of-section station convention.
    ax : matplotlib.axes.Axes, optional
        Existing axes to draw on.

    Returns
    -------
    matplotlib.axes.Axes
    """
    df = total_divergence_table(sites, spacing_m, component=component)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    if df.empty:
        ax.text(0.5, 0.5, "no tipper data", ha="center", va="center")
        return ax

    col = "divergence_real" if part == "real" else "divergence_imag"

    if frequency_hz is None and period_s is not None:
        frequency_hz = 1.0 / max(float(period_s), 1e-24)
    if frequency_hz is None:
        frequency_hz = float(df["freq_hz"].median())

    pairs = list(dict.fromkeys(zip(df["station_a"], df["station_b"])))
    values: list[float] = []
    used_freqs: list[float] = []
    for a, b in pairs:
        sub = df[(df["station_a"] == a) & (df["station_b"] == b)]
        idx = (sub["freq_hz"] - frequency_hz).abs().idxmin()
        values.append(float(sub.loc[idx, col]))
        used_freqs.append(float(sub.loc[idx, "freq_hz"]))

    labels = [a for a, _ in pairs]
    x = np.arange(len(pairs))
    _ml = PYCSAMT_STYLE.multiline
    ax.axhline(0.0, color="0.8", lw=0.8)
    ax.plot(
        x, values, lw=_ml.lw, alpha=_ml.alpha, marker="o", color="tab:blue"
    )
    _apply_station_rendering(
        ax,
        labels,
        station_label_step=station_label_step,
        station_preset=station_preset,
        station_style=station_style,
    )
    ref = float(np.nanmedian(used_freqs))
    tlabel = r"$T_{zx}$" if component == "tzx" else r"$T_{zy}$"
    ax.set_ylabel(rf"$\partial${tlabel}$/\partial x$  [{part}]  (m$^{{-1}}$)")
    ax.set_title(
        f"ZTEM total divergence profile at {ref:.4g} Hz", fontsize=10
    )
    ax.grid(True, ls=":", alpha=0.35)
    return ax


def _divergence_grid(
    line_items: Any, component: str, part: str, spacing_m: float,
) -> tuple[np.ndarray, list[str], np.ndarray]:
    """Return ``(grid, pair_labels, freqs)`` for one flight line's
    total-divergence pseudosection, shared by
    :func:`plot_ztem_divergence_psection` and
    :func:`plot_ztem_divergence_psection_grid`."""
    df = total_divergence_table(line_items, spacing_m, component=component)
    if df.empty:
        return np.zeros((0, 0)), [], np.zeros(0)
    col = "divergence_real" if part == "real" else "divergence_imag"
    pairs = list(dict.fromkeys(zip(df["station_a"], df["station_b"])))
    labels = [a for a, _ in pairs]
    freqs = np.sort(df["freq_hz"].unique())[::-1]
    grid = np.full((freqs.size, len(pairs)), np.nan, dtype=float)
    pair_idx = {p: j for j, p in enumerate(pairs)}
    freq_idx = {float(f): i for i, f in enumerate(freqs)}
    for _, row in df.iterrows():
        j = pair_idx.get((row["station_a"], row["station_b"]))
        i = freq_idx.get(float(row["freq_hz"]))
        if i is not None and j is not None:
            grid[i, j] = row[col]
    return grid, labels, freqs


def _draw_divergence_psection(
    ax: plt.Axes,
    grid: np.ndarray,
    labels: list[str],
    freqs: np.ndarray,
    *,
    cmap: str,
    vmin: float,
    vmax: float,
    show_grid: bool,
    show_contour: bool,
    n_contour_levels: int,
    station_label_step: int | None,
    station_preset: str,
    station_style: Any | None,
) -> Any:
    """Draw one divergence pseudosection panel (no colorbar of its own);
    returns the ``imshow`` mappable so callers can attach a shared or
    per-panel colorbar. Factored out of
    :func:`plot_ztem_divergence_psection` so
    :func:`plot_ztem_divergence_psection_grid` can reuse it verbatim
    for every panel with one shared colour scale."""
    from ..api.labels import LOG10_PERIOD_LABEL

    extent = (
        -0.5,
        len(labels) - 0.5,
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

    nf, ns = grid.shape
    y_top, y_bottom = extent[3], extent[2]
    dy = (y_top - y_bottom) / nf
    y_centers = y_top - (np.arange(nf) + 0.5) * dy

    if show_grid:
        ax.set_xticks(np.arange(-0.5, ns, 1.0), minor=True)
        ax.set_yticks(y_top - np.arange(nf + 1) * dy, minor=True)
        ax.grid(which="minor", color="0.25", linewidth=0.6, alpha=0.5)
        ax.tick_params(which="minor", length=0)

    if show_contour and n_contour_levels >= 3 and nf >= 2 and ns >= 2:
        finite = grid[np.isfinite(grid)]
        if finite.size:
            levels = np.linspace(vmin, vmax, n_contour_levels)[1:-1]
            if levels.size:
                xx, yy = np.meshgrid(np.arange(ns), y_centers)
                masked = np.ma.masked_invalid(grid)
                cs = ax.contour(
                    xx, yy, masked, levels=levels, colors="black",
                    linewidths=0.8,
                )
                ax.clabel(cs, fmt="%.1g", fontsize=6, inline=True)

    _apply_station_rendering(
        ax,
        labels,
        station_label_step=station_label_step,
        station_preset=station_preset,
        station_style=station_style,
    )
    return im


def plot_ztem_divergence_psection(
    sites: Any,
    *,
    component: str = "tzx",
    part: str = "real",
    spacing_m: float = 200.0,
    cmap: str = "RdBu_r",
    clim: tuple[float, float] | None = None,
    clim_pct: float = 95.0,
    show_grid: bool = True,
    show_contour: bool = True,
    n_contour_levels: int = 3,
    figsize: tuple[float, float] = (9.0, 5.0),
    station_label_step: int | None = 1,
    station_preset: str = "pseudosection",
    station_style: Any | None = None,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    r"""Plot a ZTEM total-divergence pseudosection (station x log-period).

    A diverging, zero-centred colour scale is used, matching the
    physical sign convention of a spatial derivative (positive on one
    side of an anomaly, negative on the other -- see
    :func:`total_divergence_table`). Optional cell-boundary gridlines
    and a contour overlay (default ``n_contour_levels=3``, one
    interior level -- here the physically meaningful zero-divergence
    line itself) match the same ``imshow``/``contour`` convention used
    by :func:`~pycsamt.emtools.afmag.plot_airmt_tilt_psection` and
    :mod:`pycsamt.emtools.fieldzone`'s own pseudosections. For several
    flight lines compared side by side on one shared colour scale, see
    :func:`plot_ztem_divergence_psection_grid`.

    Parameters
    ----------
    sites : Sites-like or AirborneSites-like
        Anything accepted by
        :func:`~pycsamt.emtools._core.ensure_any_sites`: ground
        :class:`~pycsamt.site.base.Sites` when the data carries an
        impedance channel too, or
        :class:`~pycsamt.airborne.site.AirborneSites` for genuine
        tipper-only ZTEM/AFMAG data (a path/directory of EMTF-XML is
        routed automatically based on what it contains).
    component : {"tzx", "tzy"}, default "tzx"
    part : {"real", "imag"}, default "real"
    spacing_m : float, default 200.0
        Forwarded to :func:`total_divergence_table`.
    cmap : str, default "RdBu_r"
        Diverging colormap name.
    clim : (float, float), optional
        Explicit, zero-centred color limits; overrides *clim_pct*.
    clim_pct : float, default 95.0
        Percentile of ``|divergence|`` used to size a symmetric colour
        range when *clim* is not given.
    show_grid : bool, default True
        Draw thin gridlines at every station-pair/period cell
        boundary.
    show_contour : bool, default True
        Overlay *n_contour_levels* - 2 evenly-spaced contour lines
        with inline labels; with the default zero-centred colour
        scale and 3 levels, the single interior level drawn is the
        zero-divergence contour itself, i.e. the crossover/conductor
        axis at every period simultaneously.
    n_contour_levels : int, default 3
        Number of evenly-spaced levels spanning *clim* before dropping
        the two outermost; must be at least 3 for any line to be
        drawn. Kept low deliberately -- a coarse station/period grid
        does not support many contour levels without the lines
        tangling into visual noise (see
        :func:`~pycsamt.emtools.afmag.plot_airmt_tilt_psection`'s
        docstring for the same reasoning).
    figsize : (float, float), default (9.0, 5.0)
        Used only when *ax* is not supplied.
    station_label_step, station_preset, station_style
        See :func:`plot_ztem_divergence_profile`.
    ax : matplotlib.axes.Axes, optional
        Existing axes to draw on.

    Returns
    -------
    matplotlib.axes.Axes
    """
    from ..api.plot import add_colorbar

    grid, labels, freqs = _divergence_grid(sites, component, part, spacing_m)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    if grid.size == 0:
        ax.text(0.5, 0.5, "no tipper data", ha="center", va="center")
        return ax

    if clim is None:
        finite = np.abs(grid[np.isfinite(grid)])
        vabs = float(np.percentile(finite, clim_pct)) if finite.size else 1.0
        vabs = max(vabs, 1e-12)
        clim = (-vabs, vabs)
    vmin, vmax = clim

    im = _draw_divergence_psection(
        ax, grid, labels, freqs, cmap=cmap, vmin=vmin, vmax=vmax,
        show_grid=show_grid, show_contour=show_contour,
        n_contour_levels=n_contour_levels,
        station_label_step=station_label_step,
        station_preset=station_preset, station_style=station_style,
    )
    tlabel = "Tzx" if component == "tzx" else "Tzy"
    ax.set_title(
        f"ZTEM total divergence pseudosection [{tlabel}, {part}]",
        fontsize=10,
    )
    add_colorbar(im, ax, label=r"$\partial T/\partial x$ (m$^{-1}$)")
    return ax


def plot_ztem_divergence_psection_grid(
    sites: Any,
    *,
    component: str = "tzx",
    part: str = "real",
    spacing_m: float = 200.0,
    max_lines: int = 6,
    n_cols: int = 3,
    cmap: str = "seismic",
    clim: tuple[float, float] | None = None,
    clim_pct: float = 95.0,
    show_grid: bool = True,
    show_contour: bool = True,
    n_contour_levels: int = 3,
    panel_size: tuple[float, float] = (4.3, 3.4),
    station_label_step: int | None = 2,
    station_preset: str = "pseudosection",
    station_style: Any | None = None,
    axes: Any | None = None,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> plt.Figure:
    r"""Compare several flight lines' divergence pseudosections at once.

    Every panel shares one colour scale (unlike calling
    :func:`plot_ztem_divergence_psection` once per line, where each
    panel would size its own), so colour differences between lines are
    directly comparable -- the multi-line counterpart of a single
    :func:`plot_ztem_divergence_psection` call, laid out on a grid
    with *n_cols* columns. Flight lines are detected from station
    geometry (see :func:`_detect_line_groups`); when more than
    *max_lines* are found, a spatially even subset is kept rather than
    just the first *max_lines* encountered.

    Parameters
    ----------
    sites : Sites-like or AirborneSites-like
        Anything accepted by
        :func:`~pycsamt.emtools._core.ensure_any_sites`, spanning
        several flight lines.
    component, part, spacing_m
        Forwarded to :func:`total_divergence_table` for every line.
    max_lines : int, default 6
        Maximum number of lines to draw.
    n_cols : int, default 3
        Number of grid columns; rows are added as needed.
    cmap : str, default "seismic"
        Diverging colormap name.
    clim : (float, float), optional
        Explicit, zero-centred color limits shared by every panel;
        overrides *clim_pct*.
    clim_pct : float, default 95.0
        Percentile of ``|divergence|``, pooled across every drawn
        line, used to size the shared symmetric colour range when
        *clim* is not given.
    show_grid, show_contour, n_contour_levels
        See :func:`plot_ztem_divergence_psection`.
    panel_size : (float, float), default (4.3, 3.4)
        Per-panel figure size in inches; the full figure scales with
        the number of rows/columns actually used. Ignored when *axes*
        is supplied.
    station_label_step, station_preset, station_style
        See :func:`plot_ztem_divergence_profile`.
    axes : sequence of Axes, optional
        Existing axes to draw the (up to *max_lines*) panels on,
        flattened in the same row-major order the auto-created grid
        would use; must provide at least as many axes as lines are
        actually drawn. When not given, a new figure and grid of axes
        is created.
    recursive, on_dup, strict, verbose
        Forwarded to :func:`~pycsamt.emtools._core.ensure_any_sites`.

    Returns
    -------
    matplotlib.Figure
    """
    from matplotlib.ticker import ScalarFormatter

    S = ensure_any_sites(
        sites, recursive=recursive, on_dup=on_dup, strict=strict,
        verbose=verbose,
    )
    items, lon, lat = _station_lonlat(S)
    if not items:
        fig, ax = plt.subplots(figsize=panel_size)
        ax.text(0.5, 0.5, "no station coordinates", ha="center", va="center")
        return fig

    groups = _detect_line_groups(lon, lat)
    group_ids = sorted(set(groups.tolist()))
    if len(group_ids) > max_lines:
        pick = np.linspace(0, len(group_ids) - 1, max_lines)
        pick = sorted(set(int(round(p)) for p in pick))
        group_ids = [group_ids[i] for i in pick]

    panels = []
    for g in group_ids:
        idx = np.where(groups == g)[0]
        if idx.size < 2:
            continue
        line_items = [items[i] for i in idx]
        grid, labels, freqs = _divergence_grid(
            line_items, component, part, spacing_m,
        )
        if grid.size == 0:
            continue
        doc = getattr(line_items[0], "emtf", None)
        line_name = None
        if doc is not None:
            line_name = (
                doc.metadata.get("notes", {})
                .get("ZTEM", {})
                .get("LineId")
            )
        panels.append((line_name or f"Line {g}", grid, labels, freqs))

    if not panels:
        fig, ax = plt.subplots(figsize=panel_size)
        ax.text(0.5, 0.5, "no tipper data", ha="center", va="center")
        return fig

    if clim is None:
        pooled = np.concatenate(
            [np.abs(g[np.isfinite(g)]).ravel() for _, g, _, _ in panels]
        )
        vabs = float(np.percentile(pooled, clim_pct)) if pooled.size else 1.0
        vabs = max(vabs, 1e-12)
        clim = (-vabs, vabs)
    vmin, vmax = clim

    n_panels = len(panels)
    n_cols_eff = max(1, min(n_cols, n_panels))
    n_rows = int(np.ceil(n_panels / n_cols_eff))
    axes_given = _axes_list(axes, n_panels) if axes is not None else None
    if axes_given is None:
        fig, axes_arr = plt.subplots(
            n_rows, n_cols_eff,
            figsize=(panel_size[0] * n_cols_eff, panel_size[1] * n_rows),
            sharey=True,
            squeeze=False,
            gridspec_kw={"hspace": 1.05, "wspace": 0.25},
        )
        axes_flat = axes_arr.ravel()
    else:
        axes_flat = np.asarray(axes_given, dtype=object)
        fig = axes_flat[0].figure

    im = None
    for k, (name, grid, labels, freqs) in enumerate(panels):
        ax = axes_flat[k]
        im = _draw_divergence_psection(
            ax, grid, labels, freqs, cmap=cmap, vmin=vmin, vmax=vmax,
            show_grid=show_grid, show_contour=show_contour,
            n_contour_levels=n_contour_levels,
            station_label_step=station_label_step,
            station_preset=station_preset, station_style=station_style,
        )
        # extra pad pushes the title clear of the top-side station
        # ticks/labels _apply_station_rendering already draws there.
        ax.set_title(str(name), fontsize=9, pad=26)
        if k % n_cols_eff != 0:
            ax.set_ylabel("")

    for k in range(n_panels, axes_flat.size):
        axes_flat[k].axis("off")

    tlabel = "Tzx" if component == "tzx" else "Tzy"
    fig.text(
        0.5, 1.0 + 0.10 / n_rows,
        f"ZTEM total divergence pseudosections [{tlabel}, {part}]"
        f" -- {n_panels} lines",
        fontsize=11, ha="center",
    )
    fmt = ScalarFormatter(useMathText=True)
    fmt.set_powerlimits((-1, 1))
    cbar = fig.colorbar(
        im, ax=axes_flat[:n_panels].tolist(), format=fmt,
        shrink=0.85, pad=0.02,
    )
    cbar.set_label(r"$\partial T/\partial x$ (m$^{-1}$)")
    return fig


def plot_ztem_phase_rotation_profile(
    sites: Any,
    *,
    component: str = "tzx",
    part: str = "real",
    frequency_hz: float | None = None,
    period_s: float | None = None,
    spacing_m: float = 200.0,
    n_resample: int | None = None,
    figsize: tuple[float, float] = (9.5, 4.2),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    r"""Plot raw vs. Hilbert-phase-rotated ZTEM response at one frequency.

    Direct reproduction of the crossover-to-peak comparison in Sattel
    and Witherly (2012, Fig. 2); see :func:`phase_rotate_table`.

    Parameters
    ----------
    sites : Sites-like or AirborneSites-like
        Anything accepted by
        :func:`~pycsamt.emtools._core.ensure_any_sites`: ground
        :class:`~pycsamt.site.base.Sites` when the data carries an
        impedance channel too, or
        :class:`~pycsamt.airborne.site.AirborneSites` for genuine
        tipper-only ZTEM/AFMAG data (a path/directory of EMTF-XML is
        routed automatically based on what it contains).
    component, part, frequency_hz, period_s, spacing_m, n_resample
        Forwarded to :func:`phase_rotate_table`.
    figsize : (float, float), default (9.5, 4.2)
        Used only when *ax* is not supplied.
    recursive, on_dup, strict, verbose
        Forwarded to :func:`phase_rotate_table`.
    ax : matplotlib.axes.Axes, optional
        Existing axes to draw on.

    Returns
    -------
    matplotlib.axes.Axes
    """
    df = phase_rotate_table(
        sites,
        frequency_hz=frequency_hz,
        period_s=period_s,
        component=component,
        part=part,
        spacing_m=spacing_m,
        n_resample=n_resample,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    if df.empty:
        ax.text(0.5, 0.5, "no tipper data", ha="center", va="center")
        return ax

    ax.axhline(0.0, color="0.8", lw=0.8)
    ax.plot(
        df["x_m"], df["raw"], lw=1.6, color="tab:blue",
        label="raw (crossover)",
    )
    ax.plot(
        df["x_m"], df["rotated"], lw=1.6, color="tab:red", ls="--",
        label="Hilbert-rotated (peak)",
    )
    ax.plot(
        df["x_m"], df["envelope"], lw=1.2, color="0.35", ls=":",
        label="envelope",
    )
    ax.set_xlabel("Position along profile (m)")
    tlabel = r"$T_{zx}$" if component == "tzx" else r"$T_{zy}$"
    ax.set_ylabel(f"{tlabel}  [{part}]")
    f0 = float(df["freq_hz"].iloc[0])
    ax.set_title(
        f"ZTEM phase-rotated response at {f0:.4g} Hz", fontsize=10
    )
    ax.legend(fontsize=8, framealpha=0.85)
    ax.grid(True, ls=":", alpha=0.3)
    return ax


def plot_ztem_band_mask_psection(
    sites: Any,
    *,
    band_hz: tuple[float, float] | None = None,
    system_spec: Any | None = None,
    component: str = "abs",
    cmap: str = "RdBu_r",
    figsize: tuple[float, float] = (9.5, 8.0),
    axes: Any | None = None,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> plt.Figure:
    r"""Plot before/after :math:`|T|` pseudosections around the ZTEM
    usable band.

    Reuses :func:`~pycsamt.emtools.tf.plot_induction_section` for both
    panels rather than re-implementing pseudosection gridding, and
    :func:`mask_outside_ztem_band` to compute the "after" sites.

    Parameters
    ----------
    sites : Sites-like
        Anything accepted by
        :func:`~pycsamt.emtools._core.ensure_any_sites`: ground
        :class:`~pycsamt.site.base.Sites` or tipper-only
        :class:`~pycsamt.airborne.site.AirborneSites` transparently,
        since :func:`~pycsamt.emtools.tf.plot_induction_section` now
        accepts both.
    band_hz, system_spec
        Forwarded to :func:`mask_outside_ztem_band`.
    component : {"real", "imag", "abs"}, default "abs"
        Forwarded to
        :func:`~pycsamt.emtools.tf.plot_induction_section`.
    cmap : str, default "RdBu_r"
    figsize : (float, float), default (9.5, 8.0)
        Used only when *axes* is not supplied.
    axes : sequence of 2 Axes, optional
        Existing axes (before, after) to draw on.
    recursive, on_dup, strict, verbose
        Forwarded to :func:`~pycsamt.emtools._core.ensure_any_sites`
        and :func:`mask_outside_ztem_band`.

    Returns
    -------
    matplotlib.Figure
    """
    from .tf import plot_induction_section

    before = ensure_any_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    after = mask_outside_ztem_band(
        before,
        band_hz=band_hz,
        system_spec=system_spec,
        action="mask",
        inplace=False,
    )

    axes_given = _axes_list(axes, 2) if axes is not None else None
    if axes_given is None:
        fig, axes_arr = plt.subplots(
            2, 1, figsize=figsize, sharex=True,
            gridspec_kw={"hspace": 0.3},
        )
    else:
        axes_arr = np.asarray(axes_given, dtype=object)
        fig = axes_arr[0].figure

    plot_induction_section(
        before, component=component, cmap=cmap, ax=axes_arr[0],
        title="Before ZTEM band mask",
    )
    plot_induction_section(
        after, component=component, cmap=cmap, ax=axes_arr[1],
        title="After ZTEM band mask",
    )
    return fig


# ─────────────────────────────────────────────────────────────────────────
# Map view (Legault et al. 2012, Fig. 7; Sattel and Witherly 2012, Fig. 7-11)
# ─────────────────────────────────────────────────────────────────────────


def _station_lonlat(sites: Any) -> tuple[list[Any], np.ndarray, np.ndarray]:
    """Return ``(items, lon, lat)`` for every site with finite coords."""
    items = []
    lon: list[float] = []
    lat: list[float] = []
    for ed in _iter_items(ensure_any_sites(sites)):
        coords = getattr(ed, "coords", None)
        if coords is None:
            continue
        try:
            la, lo = float(coords[0]), float(coords[1])
        except (TypeError, ValueError, IndexError):
            continue
        if np.isfinite(la) and np.isfinite(lo):
            items.append(ed)
            lon.append(lo)
            lat.append(la)
    return items, np.asarray(lon), np.asarray(lat)


def _detect_line_groups(
    lon: np.ndarray, lat: np.ndarray, n_bins: int = 50
) -> np.ndarray:
    r"""Group stations into flight lines from their (lon, lat) alone.

    A parallel, axis-aligned block survey has one coordinate that
    stays nearly constant *within* a line and jumps *between* lines
    (the cross-line axis) and another that varies continuously along
    each line (the along-line axis). This picks whichever of
    longitude/latitude clusters into fewer, coarser groups after
    binning at roughly ``span / n_bins`` and uses that as the
    cross-line axis. It is a real, general heuristic for this common
    survey geometry, not a full arbitrary-azimuth line detector.

    Returns
    -------
    numpy.ndarray of int
        One line-group index per input station, in encounter order
        within each detected group.
    """
    if lon.size == 0:
        return np.zeros(0, dtype=int)

    def _groups(coord: np.ndarray) -> np.ndarray:
        span = float(np.ptp(coord)) if coord.size > 1 else 0.0
        tol = max(span / max(n_bins, 1), 1e-9)
        rounded = np.round(coord / tol) * tol
        uniq = sorted(set(rounded.tolist()))
        idx = {v: i for i, v in enumerate(uniq)}
        return np.array([idx[v] for v in rounded.tolist()])

    g_lat = _groups(lat)
    g_lon = _groups(lon)
    return g_lat if g_lat.max(initial=0) <= g_lon.max(initial=0) else g_lon


def plot_ztem_flight_lines(
    sites: Any,
    *,
    figsize: tuple[float, float] = (7.0, 6.0),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    r"""Plot a Sattel and Witherly (2012, Fig. 7)-style flight-line map.

    Every detected flight line (see :func:`_detect_line_groups`) is
    drawn as its own connected navigation trace, coloured distinctly
    (a ``viridis`` sample per line) and labelled near its first
    station -- with the real flight-line identifier when a
    technology note carries one (e.g. ZTEM's own
    ``metadata["notes"]["ZTEM"]["LineId"]``), or else a generic
    ``L1``, ``L2``, ... in detected-group order, which need not match
    any real line numbering -- with station markers. The plan-view
    counterpart of every other function in this module, which reads a
    single profile at a time. Longitude tick labels are rotated 45
    degrees to avoid overlapping.

    Parameters
    ----------
    sites : Sites-like or AirborneSites-like
        Anything accepted by
        :func:`~pycsamt.emtools._core.ensure_any_sites`.
    figsize : (float, float), default (7.0, 6.0)
        Used only when *ax* is not supplied.
    recursive, on_dup, strict, verbose
        Forwarded to :func:`~pycsamt.emtools._core.ensure_any_sites`.
    ax : matplotlib.axes.Axes, optional
        Existing axes to draw on.

    Returns
    -------
    matplotlib.axes.Axes
    """
    S = ensure_any_sites(
        sites, recursive=recursive, on_dup=on_dup, strict=strict,
        verbose=verbose,
    )
    items, lon, lat = _station_lonlat(S)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    if lon.size == 0:
        ax.text(0.5, 0.5, "no station coordinates", ha="center", va="center")
        return ax

    groups = _detect_line_groups(lon, lat)
    group_ids = sorted(set(groups.tolist()))
    line_cmap = plt.get_cmap("viridis", max(len(group_ids), 1))

    # Resolve every line's display label first, in a separate pass,
    # so a common prefix shared by every real flight-line identifier
    # (e.g. "gold_springs_nv_L" in "gold_springs_nv_L1" ..
    # "gold_springs_nv_L7") can be stripped before drawing -- the
    # full identifier is unique but far too wide to label a compact
    # map with, and every line sharing the same survey prefix carries
    # no distinguishing information anyway.
    orders: list[np.ndarray] = []
    labels: list[str] = []
    for g in group_ids:
        idx = np.where(groups == g)[0]
        order = (
            idx[np.argsort(lon[idx])]
            if np.ptp(lon[idx]) >= np.ptp(lat[idx])
            else idx[np.argsort(lat[idx])]
        )
        orders.append(order)
        doc = getattr(items[idx[0]], "emtf", None)
        line_name = (
            doc.metadata.get("notes", {}).get("ZTEM", {}).get("LineId")
            if doc is not None
            else None
        )
        labels.append(line_name or f"L{len(labels) + 1}")

    if len(labels) > 1:
        prefix = os.path.commonprefix(labels)
        # only strip up to the last safe boundary, and only if doing
        # so still leaves every label non-empty and distinguishable
        cut = max(prefix.rfind("_"), prefix.rfind("-")) + 1
        if cut > 0 and all(len(label) > cut for label in labels):
            labels = [label[cut:] for label in labels]

    for k, order in enumerate(orders):
        gx, gy = lon[order], lat[order]
        color = line_cmap(k)
        ax.plot(gx, gy, "-", color=color, lw=1.4, zorder=1)
        ax.scatter(
            gx, gy, s=26, color=color, edgecolors="white", linewidths=0.4,
            zorder=2,
        )
        # label each line once, just past its first station, alternating
        # sides so consecutive close-together lines don't collide
        dx = (float(np.ptp(lon)) or 1.0) * 0.015
        ha = "left" if k % 2 == 0 else "right"
        label_x = gx[0] + dx if ha == "left" else gx[0] - dx
        ax.text(
            label_x, gy[0], labels[k],
            ha=ha, va="center", fontsize=8, fontweight="bold", color=color,
        )

    ax.margins(x=0.10)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(
        f"ZTEM flight lines ({len(group_ids)} lines, {len(items)} stations)",
        fontsize=10,
    )
    ax.ticklabel_format(useOffset=False, style="plain")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    ax.grid(True, ls=":", alpha=0.3)
    return ax


def plot_ztem_map(
    sites: Any,
    *,
    quantity: str = "tipper",
    part: str = "real",
    component: str = "tzx",
    frequency_hz: float | None = None,
    period_s: float | None = None,
    n_grid: int = 120,
    cmap: str = "RdBu_r",
    clim: tuple[float, float] | None = None,
    clim_pct: float = 95.0,
    show_stations: bool = True,
    figsize: tuple[float, float] = (8.0, 6.5),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    r"""Legault et al. (2012, Fig. 7) and Sattel and Witherly (2012,
    Fig. 8-11)-style map-view grid.

    Interpolates one scalar field at one frequency across every
    flight line in *sites* onto a regular map grid
    (:func:`scipy.interpolate.griddata`, linear inside the convex
    hull of the stations, unfilled -- left ``nan`` -- outside it
    rather than extrapolated) and images it with a diverging,
    zero-centred colour scale -- the genuine multi-line map product
    both papers show (their "DT map"/"XIP grid"/"phase-rotated grid"),
    as opposed to every other function in this module, which reads
    one flight line as a profile or pseudosection.

    Parameters
    ----------
    sites : Sites-like or AirborneSites-like
        Anything accepted by
        :func:`~pycsamt.emtools._core.ensure_any_sites`. A genuine
        map needs several roughly-parallel flight lines; a single
        line still renders, as a thin interpolated strip along it.
    quantity : {"tipper", "divergence"}, default "tipper"
        ``"tipper"`` images the raw, un-processed tipper component
        (Legault et al. 2012, Fig. 7's own "In-Phase" map);
        ``"divergence"`` images the along-line total-divergence /
        Peaker value (:func:`total_divergence_table`) at each
        station's own flight line, matching Sattel and Witherly
        (2012)'s "DT" grid.
    part : {"real", "imag"}, default "real"
    component : {"tzx", "tzy"}, default "tzx"
    frequency_hz, period_s : float, optional
        Reference frequency/period; nearest available value is used
        per station. At most one may be given; the median frequency
        is used when neither is given.
    n_grid : int, default 120
        Number of grid points along the longer map axis; the shorter
        axis is scaled to preserve the survey's aspect ratio.
    cmap : str, default "RdBu_r"
    clim : (float, float), optional
        Explicit, zero-centred color limits; overrides *clim_pct*.
    clim_pct : float, default 95.0
        Percentile of ``|value|`` used to size a symmetric colour
        range when *clim* is not given.
    show_stations : bool, default True
        Overlay the actual station positions as small markers.
    figsize : (float, float), default (8.0, 6.5)
        Used only when *ax* is not supplied.
    recursive, on_dup, strict, verbose
        Forwarded to :func:`~pycsamt.emtools._core.ensure_any_sites`.
    ax : matplotlib.axes.Axes, optional
        Existing axes to draw on.

    Returns
    -------
    matplotlib.axes.Axes

    Raises
    ------
    ValueError
        If *quantity* is not ``"tipper"``/``"divergence"``, or
        *part*/*component* is invalid.
    """
    from scipy.interpolate import griddata

    quantity = str(quantity).strip().lower()
    if quantity not in ("tipper", "divergence"):
        raise ValueError(
            f"quantity must be 'tipper' or 'divergence'; got {quantity!r}"
        )
    component = str(component).strip().lower()
    if component not in ("tzx", "tzy"):
        raise ValueError(
            f"component must be 'tzx' or 'tzy'; got {component!r}"
        )
    part = str(part).strip().lower()
    if part not in ("real", "imag"):
        raise ValueError(f"part must be 'real' or 'imag'; got {part!r}")

    S = ensure_any_sites(
        sites, recursive=recursive, on_dup=on_dup, strict=strict,
        verbose=verbose,
    )
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    items, lon, lat = _station_lonlat(S)
    if items:
        groups = _detect_line_groups(lon, lat)
    else:
        groups = np.zeros(0, dtype=int)

    if quantity == "tipper":
        values_l: list[float] = []
        keep: list[int] = []
        for i, ed in enumerate(items):
            T, t, fr = _get_t_block(ed)
            if T is None or t is None or fr is None or fr.size == 0:
                continue
            k, _ = _resolve_target_frequency(fr, frequency_hz, period_s)
            v = t[k, 0] if component == "tzx" else t[k, 1]
            values_l.append(v.real if part == "real" else v.imag)
            keep.append(i)
        lon, lat = lon[keep], lat[keep]
        values = np.asarray(values_l, dtype=float)
        tlabel = r"$T_{zx}$" if component == "tzx" else r"$T_{zy}$"
        label = f"{tlabel} [{part}]"
    else:
        # total_divergence_table's along-profile chainage assumes one
        # flight line: differentiating across a line boundary would
        # produce a physically meaningless value. Each detected line
        # is therefore differentiated separately, then the resulting
        # (lon, lat, divergence) triples from every line are pooled
        # for one combined map.
        col = "divergence_real" if part == "real" else "divergence_imag"
        lon_l: list[float] = []
        lat_l: list[float] = []
        values_l = []
        for g in sorted(set(groups.tolist())):
            idx = np.where(groups == g)[0]
            if idx.size < 2:
                continue
            line_items = [items[i] for i in idx]
            df = total_divergence_table(line_items, component=component)
            if df.empty:
                continue
            freqs_here = np.sort(df["freq_hz"].unique())
            _, f0 = _resolve_target_frequency(
                freqs_here, frequency_hz, period_s,
            )
            sub = df[np.isclose(df["freq_hz"], f0)]
            name_to_ll = {
                _name(ed, j): (lon[idx[j]], lat[idx[j]])
                for j, ed in enumerate(line_items)
            }
            for _, row in sub.iterrows():
                ll = name_to_ll.get(row["station_a"])
                if ll is None:
                    continue
                lon_l.append(ll[0])
                lat_l.append(ll[1])
                values_l.append(row[col])
        lon = np.asarray(lon_l)
        lat = np.asarray(lat_l)
        values = np.asarray(values_l, dtype=float)
        label = rf"$\partial T/\partial x$ [{part}]  (m$^{{-1}}$)"
    if lon.size < 3 or values.size < 3:
        ax.text(
            0.5, 0.5, "not enough stations to map", ha="center",
            va="center",
        )
        return ax

    lon_span = max(float(np.ptp(lon)), 1e-9)
    lat_span = max(float(np.ptp(lat)), 1e-9)
    if lon_span >= lat_span:
        nx = n_grid
        ny = max(int(round(n_grid * lat_span / lon_span)), 8)
    else:
        ny = n_grid
        nx = max(int(round(n_grid * lon_span / lat_span)), 8)
    grid_lon = np.linspace(lon.min(), lon.max(), nx)
    grid_lat = np.linspace(lat.min(), lat.max(), ny)
    gx, gy = np.meshgrid(grid_lon, grid_lat)
    grid_z = griddata(
        (lon, lat), values, (gx, gy), method="linear",
    )

    if clim is None:
        finite = np.abs(values[np.isfinite(values)])
        vabs = float(np.percentile(finite, clim_pct)) if finite.size else 1.0
        vabs = max(vabs, 1e-12)
        clim = (-vabs, vabs)
    vmin, vmax = clim

    im = ax.pcolormesh(
        gx, gy, grid_z, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto",
    )
    if show_stations:
        ax.scatter(
            lon, lat, s=10, color="black", alpha=0.5, zorder=3,
        )
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.ticklabel_format(useOffset=False, style="plain")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    q_label = "In-Phase" if part == "real" else "Quadrature"
    title_q = "raw tipper" if quantity == "tipper" else "total divergence"
    ax.set_title(f"ZTEM {q_label} {title_q} map", fontsize=10)
    from ..api.plot import add_colorbar

    add_colorbar(im, ax, label=label)
    return ax
