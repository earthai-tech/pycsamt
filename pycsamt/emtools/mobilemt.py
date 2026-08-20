# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

r"""MobileMT-specific processing, diagnostics, and plotting.

MobileMT (Prikhodko et al. 2022) is a natural-source airborne EM
technology built on a fundamentally different scientific object than
ZTEM/AFMAG: three orthogonal airborne magnetic coils
(:math:`H_x, H_y, H_z`) referenced to a fixed-ground horizontal
electric dipole pair (:math:`E_x, E_y`), giving a complex
**admittance** tensor

.. math::

    \begin{pmatrix} H_x \\ H_y \\ H_z \end{pmatrix} =
    \begin{pmatrix}
        Y_{xx} & Y_{xy} \\
        Y_{yx} & Y_{yy} \\
        Y_{hzx} & Y_{hzy}
    \end{pmatrix}
    \begin{pmatrix} E_x \\ E_y \end{pmatrix}

-- the reciprocal relation of the classical MT impedance tensor
:math:`Z` (which gives :math:`E` from :math:`H`), not a tipper. This
is why :mod:`pycsamt.airborne.mobilemt` maps it onto a dedicated
``mobilemt_admittance`` :class:`~pycsamt.emtf.TransferFunction` of
shape ``(nf, 3, 2)`` rather than onto
:attr:`~pycsamt.site.base.Site.z` or
:attr:`~pycsamt.site.base.Site.tipper`, and why converting it to an
EDI/:class:`~pycsamt.site.base.Site` is refused outright
(``EMTF.to_edi`` raises ``EMTFEDIConversionError`` for this
transfer-function type -- see
``pycsamt.airborne.mobilemt.tests.test_mobilemt_interop``). Every
``emtools`` module before this one accepts and returns
:class:`~pycsamt.site.base.Sites`; this module cannot honestly do
that, so it works instead on the container hierarchy that already
carries MobileMT's real scientific content:
:class:`~pycsamt.airborne.AirborneEMDataset` ->
:class:`~pycsamt.airborne.AirborneEMLine` ->
:class:`~pycsamt.airborne.AirborneEMRecord`. The public functions
below accept a dataset (or one line, normalized through
:func:`ensure_mobilemt_dataset`, this module's ``ensure_sites``
counterpart) and return either a tidy table or a
:class:`~pycsamt.airborne.AirborneEMDataset` -- the closest honest
analogue of the rest of ``emtools``'s "sites in, sites/table out"
contract.

:func:`ensure_mobilemt_dataset` also accepts anything
:func:`~pycsamt.airborne.site.ensure_asites` does -- an
:class:`~pycsamt.airborne.site.AirborneSites`/
:class:`~pycsamt.airborne.site.AirborneSite`, or a bare path/directory
of EMTF-XML files -- regrouping it into an
:class:`~pycsamt.airborne.AirborneEMDataset` by each site's
:attr:`~pycsamt.airborne.site.AirborneSite.line_id` (a fresh,
single-line grouping when that is unset) so a
:class:`~pycsamt.airborne.site.AirborneSites` produced elsewhere in
``emtools`` -- or a directory of raw MobileMT EMTF-XML -- can flow
straight into this module without an intermediate
:class:`~pycsamt.airborne.AirborneEMDataset` construction step. This
is deliberately one-directional, unlike ZTEM/AFMAG's
``ensure_any_sites``: the module still only ever *returns* a dataset
(never :class:`~pycsamt.airborne.site.AirborneSites`), since its
functions are organized around flight lines
(:func:`plot_mobilemt_admittance_profile` and friends plot one line's
along-line chainage), not a flat station list.

Two kinds of quantity are computed here, and the distinction is kept
visible in every column name:

* **Scale-invariant tensor diagnostics** -- :func:`admittance_skew_table`
  (a Swift 1967-style skew ratio, :math:`|Y_{xx}+Y_{yy}|
  /|Y_{xy}-Y_{yx}|`, applied to the horizontal 2x2 admittance
  submatrix by direct algebraic analogy to the same ratio already
  used for the impedance tensor). Being a ratio of magnitudes, it
  needs no absolute physical constant and is safe to compute directly
  from any admittance tensor.
* **Theoretical apparent conductivity/phase** --
  :func:`admittance_determinant_table`'s ``theoretical_*`` columns.
  In the co-located-sensor limit, the MobileMT admittance tensor
  equals the classical MT admittance :math:`Z^{-1}` (stated
  explicitly by Zhdanov et al. 2024 and Sattel et al. 2019). Applying
  that identity to pyCSAMT's own, already-shipped and tested
  Berdichevsky-determinant convention for :math:`Z`
  (:class:`pycsamt.z.resphase.ResPhase`'s ``res_det``/``phase_det``,
  :math:`\rho_a = 0.2\,|\det Z|/f`, :math:`\varphi = \arg\sqrt{\det Z}`)
  gives, by direct algebraic substitution (:math:`\det Y = 1/\det Z`):

  .. math::

      \sigma_a = 5\,f\,|\det Y|, \qquad
      \varphi_a = -\arg\sqrt{\det Y}

  This is a derived theoretical quantity, not a reproduction of
  MobileMT's proprietary processed apparent-conductivity product --
  :func:`~pycsamt.airborne.mobilemt.build_mobilemt_record`'s own
  docstring explicitly declines to derive that vendor quantity
  ("a verified delivery schema is required before pyCSAMT should
  codify its exact exported representation"), and this module
  respects that same restraint by never presenting the derived
  ``theoretical_*`` columns as the vendor product. When a record
  already carries the vendor-delivered
  ``fields["apparent_conductivity"]``, every table here reports it
  alongside as ``apparent_conductivity_native_Sm`` for direct
  comparison, and :func:`plot_mobilemt_conductivity_psection` can
  plot either one explicitly via its ``source`` argument.

References
----------
.. [Prikhodko2022] Prikhodko, A., Bagrianski, A., Kuzmin, P., and
   Sirohey, A. (2022). Natural field airborne electromagnetics --
   history of development and current exploration capabilities.
   Minerals, 12(5), 583.
.. [Sattel2019] Sattel, D., Witherly, K., and Kaminski, V. (2019). A
   brief analysis of MobileMT data. SEG International Exposition and
   Annual Meeting, D043S102R007.
.. [Zhdanov2024] Zhdanov, M. S., Gribenko, A., Prikhodko, A., Sabra,
   H. E., Jorgensen, M., and Cox, L. H. (2024). Three-dimensional
   MobileMT and TMI data inversions for mineral exploration. 1st ASEG
   DISCOVER Symposium.
.. [Swift1967] Swift, C. M. (1967). A magnetotelluric investigation
   of an electrical conductivity anomaly in the southwestern United
   States. PhD thesis, MIT.
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..airborne import (
    AirborneEMDataset,
    AirborneEMLine,
    AirborneEMRecord,
    NavigationTrack,
)
from ..airborne.mobilemt import (
    MOBILEMT_ADMITTANCE_TAG,
    MOBILEMT_APPARENT_CONDUCTIVITY_FIELD,
    MobileMTSystemSpec,
)
from ..api.labels import LOG10_PERIOD_LABEL
from ..api.plot import add_colorbar
from ..api.style import PYCSAMT_STYLE
from ._core import _WGS84_ELLIPSOID_ID

__all__ = [
    "ensure_mobilemt_dataset",
    # tables
    "admittance_table",
    "admittance_determinant_table",
    "admittance_skew_table",
    # pipeline QC (dataset in, dataset out)
    "mask_outside_mobilemt_band",
    # plots
    "plot_mobilemt_admittance_profile",
    "plot_mobilemt_conductivity_psection",
    "plot_mobilemt_skew_profile",
]

_COMPONENT_INDEX = {
    "xx": (0, 0),
    "xy": (0, 1),
    "yx": (1, 0),
    "yy": (1, 1),
    "hzx": (2, 0),
    "hzy": (2, 1),
}

_ADMITTANCE_COLS = [
    "line_id",
    "sample_id",
    "x_m",
    "freq_hz",
    "period_s",
    "Yxx_real",
    "Yxx_imag",
    "Yxy_real",
    "Yxy_imag",
    "Yyx_real",
    "Yyx_imag",
    "Yyy_real",
    "Yyy_imag",
    "Yhzx_real",
    "Yhzx_imag",
    "Yhzy_real",
    "Yhzy_imag",
    "apparent_conductivity_native_Sm",
]

_DETERMINANT_COLS = [
    "line_id",
    "sample_id",
    "x_m",
    "freq_hz",
    "period_s",
    "det_abs",
    "theoretical_sigma_a_Sm",
    "theoretical_rho_a_ohm_m",
    "theoretical_phase_deg",
    "apparent_conductivity_native_Sm",
]

_SKEW_COLS = [
    "line_id",
    "sample_id",
    "x_m",
    "freq_hz",
    "period_s",
    "skew",
]


# ─────────────────────────────────────────────────────────────────────────
# Private helpers
# ─────────────────────────────────────────────────────────────────────────


def _line_chainage(
    line: AirborneEMLine, spacing_m: float = 50.0
) -> np.ndarray:
    """Return per-sample chainage [m] along one flight line.

    Mirrors :func:`~pycsamt.emtools._core._station_positions`, adapted
    to a :class:`~pycsamt.airborne.NavigationTrack`: prefers real
    easting/northing, falls back to a local-metre projection of
    latitude/longitude, and falls back further to a uniform
    ``spacing_m`` step when neither is usable.
    """
    nav = line.navigation
    n = nav.n_samples
    east = north = None
    if nav.has_projected_coordinates:
        east = np.asarray(nav.easting, dtype=float)
        north = np.asarray(nav.northing, dtype=float)
    elif nav.has_geographic_coordinates:
        from ..gis.utils import ll_to_utm

        lat = np.asarray(nav.latitude, dtype=float)
        lon = np.asarray(nav.longitude, dtype=float)
        east = np.full(n, np.nan)
        north = np.full(n, np.nan)
        zones: set[str] = set()
        for i in range(n):
            if not (np.isfinite(lat[i]) and np.isfinite(lon[i])):
                continue
            try:
                zone, e, no = ll_to_utm(_WGS84_ELLIPSOID_ID, lat[i], lon[i])
            except Exception:
                continue
            east[i], north[i] = e, no
            zones.add(zone)
        if len(zones) > 1:
            east = north = None

    if east is None or north is None:
        return np.arange(n, dtype=float) * spacing_m

    valid = np.isfinite(east) & np.isfinite(north)
    idx = np.where(valid)[0]
    if idx.size < 2:
        return np.arange(n, dtype=float) * spacing_m

    origin = np.array([east[idx[0]], north[idx[0]]])
    direction = np.array([east[idx[-1]], north[idx[-1]]]) - origin
    norm = float(np.linalg.norm(direction))
    if norm < 1.0:
        return np.arange(n, dtype=float) * spacing_m
    direction /= norm

    pos = np.arange(n, dtype=float) * spacing_m
    for i in idx:
        pos[i] = float(
            np.dot([east[i] - origin[0], north[i] - origin[1]], direction)
        )
    return pos


def _dataset_records(
    ds: AirborneEMDataset,
) -> list[tuple[str, str, float, AirborneEMRecord]]:
    """Return ``(line_id, sample_id, x_m, record)`` for every record."""
    out: list[tuple[str, str, float, AirborneEMRecord]] = []
    for line in ds.iter_lines():
        pos = _line_chainage(line)
        for idx, sample_id in enumerate(line.navigation.sample_ids):
            record = line.records.get(sample_id)
            if record is None:
                continue
            out.append((line.line_id, sample_id, float(pos[idx]), record))
    return out


def _get_admittance(
    record: AirborneEMRecord,
) -> tuple[Any, np.ndarray | None]:
    """Return ``(TransferFunction, freq)`` for one MobileMT record."""
    if record.emtf is None:
        return None, None
    tf = record.emtf.get_transfer_function(MOBILEMT_ADMITTANCE_TAG)
    if tf is None:
        return None, None
    freq = record.emtf.frequency
    if freq is None:
        return None, None
    return tf, np.asarray(freq, dtype=float)


def _native_sigma(record: AirborneEMRecord) -> np.ndarray | None:
    value = record.fields.get(MOBILEMT_APPARENT_CONDUCTIVITY_FIELD)
    return None if value is None else np.asarray(value, dtype=float)


def _component_series(Y: np.ndarray, component: str) -> np.ndarray:
    """Return one complex admittance component series, or its determinant."""
    component = str(component).strip().lower()
    if component == "det":
        y2 = Y[:, :2, :]
        return y2[:, 0, 0] * y2[:, 1, 1] - y2[:, 0, 1] * y2[:, 1, 0]
    idx = _COMPONENT_INDEX.get(component)
    if idx is None:
        raise ValueError(
            "component must be one of "
            f"{sorted(_COMPONENT_INDEX) + ['det']}; got {component!r}"
        )
    return Y[:, idx[0], idx[1]]


def _resolve_target_frequency(
    freqs: np.ndarray,
    frequency_hz: float | None,
    period_s: float | None,
) -> tuple[int, float]:
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


def _recover_native_apparent_conductivity(doc: Any) -> np.ndarray | None:
    r"""Recover a native ``apparent_conductivity`` vector lost on read.

    :func:`~pycsamt.airborne.site.AirborneSite.from_xml` builds a bare
    ``AirborneEMRecord(sample_id=..., emtf=doc)``: it never carries
    ``AirborneEMRecord.fields`` (an in-memory-only attribute with no
    EMTF-XML representation of its own). The synthetic-data generator
    already works around this on write by stashing the vector as
    comma-separated text in ``doc.metadata["notes"]["MobileMT"]
    ["ApparentConductivitySm"]`` -- the one part of an EMTF-XML
    document that does round-trip losslessly -- but nothing previously
    read it back out again, so every MobileMT survey loaded through
    :func:`ensure_mobilemt_dataset` (path/:class:`AirborneSites`
    input) silently reported ``apparent_conductivity_native_Sm`` as
    all-``nan``, even for files that do carry the vendor value. This
    reads that same note back into a real array whenever present.
    """
    if doc is None:
        return None
    raw = doc.metadata.get("notes", {}).get("MobileMT", {}).get(
        "ApparentConductivitySm"
    )
    if not raw:
        return None
    try:
        return np.array(
            [float(v) for v in str(raw).split(",") if v.strip()],
            dtype=float,
        )
    except ValueError:
        return None


def _dataset_from_asites(asites: Any) -> AirborneEMDataset:
    """Regroup an :class:`~pycsamt.airborne.site.AirborneSites` into a
    dataset, one flight line per distinct ``line_id`` (unset ids all
    fall into one ``"L001"`` line).

    Each new :class:`~pycsamt.airborne.NavigationTrack` carries
    latitude/longitude/elevation straight from
    :attr:`~pycsamt.airborne.site.AirborneSite.coords` whenever at
    least one site in the group has a finite value, so
    :func:`_line_chainage` can still resolve real along-line
    positions afterward. Every underlying
    :class:`~pycsamt.airborne.AirborneEMRecord` is reused as-is (no
    copy) *unless* its native apparent conductivity can be recovered
    via :func:`_recover_native_apparent_conductivity` and is not
    already present, in which case a new record carrying the
    recovered field replaces it -- the original object the caller
    holds is never mutated in place.
    """
    groups: dict[str, list] = {}
    for s in asites:
        groups.setdefault(s.line_id or "L001", []).append(s)

    lines: dict[str, AirborneEMLine] = {}
    for line_id, items in groups.items():
        sample_ids = tuple(s.sample_id for s in items)
        lat = np.array([s.coords[0] for s in items], dtype=float)
        lon = np.array([s.coords[1] for s in items], dtype=float)
        elev = np.array([s.coords[2] for s in items], dtype=float)
        nav_kwargs: dict[str, Any] = {"sample_ids": sample_ids}
        if np.any(np.isfinite(lat)) and np.any(np.isfinite(lon)):
            nav_kwargs["latitude"] = lat
            nav_kwargs["longitude"] = lon
        if np.any(np.isfinite(elev)):
            nav_kwargs["terrain_elevation"] = elev
        line = AirborneEMLine(
            line_id=line_id, navigation=NavigationTrack(**nav_kwargs),
        )
        for s in items:
            record = s.record
            if MOBILEMT_APPARENT_CONDUCTIVITY_FIELD not in record.fields:
                recovered = _recover_native_apparent_conductivity(s.emtf)
                if recovered is not None:
                    record = AirborneEMRecord(
                        sample_id=record.sample_id,
                        emtf=record.emtf,
                        fields={
                            **record.fields,
                            MOBILEMT_APPARENT_CONDUCTIVITY_FIELD: recovered,
                        },
                        quality=record.quality,
                        attrs=record.attrs,
                    )
            line.add_record(record)
        lines[line.line_id] = line

    return AirborneEMDataset(name="mobilemt", lines=lines)


def _pick_line(ds: AirborneEMDataset, line_id: str | None) -> AirborneEMLine:
    if line_id is not None:
        line = ds.get_line(line_id)
        if line is None:
            raise ValueError(f"line not found: {line_id!r}")
        return line
    for line in ds.iter_lines():
        return line
    raise ValueError("dataset has no flight lines")


# ─────────────────────────────────────────────────────────────────────────
# Input normalization (this module's ``ensure_sites`` counterpart)
# ─────────────────────────────────────────────────────────────────────────


def ensure_mobilemt_dataset(obj: Any) -> AirborneEMDataset:
    """Normalize a dataset or single line to an :class:`AirborneEMDataset`.

    The single entry-point validator for every public function in
    this module, mirroring the role
    :func:`~pycsamt.emtools._core.ensure_sites` plays for the rest of
    ``emtools``.

    Parameters
    ----------
    obj : AirborneEMDataset or AirborneEMLine or AirborneSites or \
AirborneSite or str or pathlib.Path
        Accepted as-is when already a dataset; a single line is
        wrapped in a new one-line dataset. An
        :class:`~pycsamt.airborne.site.AirborneSites`/
        :class:`~pycsamt.airborne.site.AirborneSite`, or a path to a
        single EMTF-XML file or a directory of them, is first coerced
        via :func:`~pycsamt.airborne.site.ensure_asites` and then
        regrouped into flight lines by
        :attr:`~pycsamt.airborne.site.AirborneSite.line_id` (see
        :func:`_dataset_from_asites`).

    Returns
    -------
    AirborneEMDataset

    Raises
    ------
    TypeError
        If *obj* is none of the accepted types.
    """
    if isinstance(obj, AirborneEMDataset):
        return obj
    if isinstance(obj, AirborneEMLine):
        return AirborneEMDataset(name=obj.line_id, lines={obj.line_id: obj})

    from ..airborne.site import AirborneSite, AirborneSites, ensure_asites

    if isinstance(obj, AirborneSite):
        return _dataset_from_asites(AirborneSites([obj]))
    if isinstance(obj, (AirborneSites, str, Path)):
        return _dataset_from_asites(ensure_asites(obj))
    raise TypeError(
        "ensure_mobilemt_dataset: expected an AirborneEMDataset, "
        "AirborneEMLine, AirborneSites, AirborneSite, or a path to "
        f"EMTF-XML; got {type(obj).__name__!r}"
    )


# ─────────────────────────────────────────────────────────────────────────
# Tables
# ─────────────────────────────────────────────────────────────────────────


def admittance_table(dataset: Any) -> pd.DataFrame:
    r"""Return a tidy per-(line, sample, frequency) admittance table.

    Parameters
    ----------
    dataset : AirborneEMDataset or AirborneEMLine
        Anything accepted by :func:`ensure_mobilemt_dataset`.

    Returns
    -------
    pandas.DataFrame
        Columns: ``line_id``, ``sample_id``, ``x_m`` (chainage along
        the flight line, see :func:`~pycsamt.emtools._core
        ._station_positions`'s along-profile convention),
        ``freq_hz``, ``period_s``, the real/imaginary parts of every
        entry of the horizontal 2x2 admittance
        (``Yxx``, ``Yxy``, ``Yyx``, ``Yyy``) and of the vertical-field
        row (``Yhzx``, ``Yhzy``), and
        ``apparent_conductivity_native_Sm`` -- the vendor-delivered
        processed field (:data:`~pycsamt.airborne.mobilemt
        .MOBILEMT_APPARENT_CONDUCTIVITY_FIELD`) when present, ``NaN``
        otherwise.
    """
    ds = ensure_mobilemt_dataset(dataset)
    rows: list[dict[str, Any]] = []
    for line_id, sample_id, x_m, record in _dataset_records(ds):
        tf, freq = _get_admittance(record)
        if tf is None or freq is None:
            continue
        Y = np.asarray(tf.data, dtype=complex)
        sigma_native = _native_sigma(record)
        for k in range(freq.size):
            y = Y[k]
            rows.append(
                {
                    "line_id": line_id,
                    "sample_id": sample_id,
                    "x_m": x_m,
                    "freq_hz": float(freq[k]),
                    "period_s": 1.0 / max(float(freq[k]), 1e-12),
                    "Yxx_real": float(y[0, 0].real),
                    "Yxx_imag": float(y[0, 0].imag),
                    "Yxy_real": float(y[0, 1].real),
                    "Yxy_imag": float(y[0, 1].imag),
                    "Yyx_real": float(y[1, 0].real),
                    "Yyx_imag": float(y[1, 0].imag),
                    "Yyy_real": float(y[1, 1].real),
                    "Yyy_imag": float(y[1, 1].imag),
                    "Yhzx_real": float(y[2, 0].real),
                    "Yhzx_imag": float(y[2, 0].imag),
                    "Yhzy_real": float(y[2, 1].real),
                    "Yhzy_imag": float(y[2, 1].imag),
                    "apparent_conductivity_native_Sm": (
                        float(sigma_native[k])
                        if sigma_native is not None
                        and k < sigma_native.size
                        and np.isfinite(sigma_native[k])
                        else np.nan
                    ),
                }
            )
    if not rows:
        return pd.DataFrame(columns=_ADMITTANCE_COLS)
    return pd.DataFrame(rows, columns=_ADMITTANCE_COLS)


def admittance_determinant_table(dataset: Any) -> pd.DataFrame:
    r"""Return the theoretical Berdichevsky-determinant admittance table.

    See the module docstring for the full derivation. In brief, using
    the horizontal 2x2 admittance submatrix
    :math:`Y = \begin{pmatrix}Y_{xx}&Y_{xy}\\Y_{yx}&Y_{yy}\end{pmatrix}`
    and the co-located-sensor identity :math:`Y=Z^{-1}` (Zhdanov et
    al. 2024; Sattel et al. 2019), applying pyCSAMT's own
    :math:`Z`-determinant convention
    (:class:`pycsamt.z.resphase.ResPhase`) by substitution gives:

    .. math::

        Y_{\mathrm{eff}} = \sqrt{\det Y}, \qquad
        \sigma_a = 5\,f\,|Y_{\mathrm{eff}}|^2, \qquad
        \varphi_a = -\arg(Y_{\mathrm{eff}})

    Parameters
    ----------
    dataset : AirborneEMDataset or AirborneEMLine
        Anything accepted by :func:`ensure_mobilemt_dataset`.

    Returns
    -------
    pandas.DataFrame
        Columns: ``line_id``, ``sample_id``, ``x_m``, ``freq_hz``,
        ``period_s``, ``det_abs`` (:math:`|\det Y|`),
        ``theoretical_sigma_a_Sm``, ``theoretical_rho_a_ohm_m``
        (:math:`1/\sigma_a`), ``theoretical_phase_deg``, and
        ``apparent_conductivity_native_Sm`` (the vendor-delivered
        field, for direct comparison, ``NaN`` when absent). Samples
        with a non-finite determinant are omitted.

    Notes
    -----
    The ``theoretical_*`` columns are a derived quantity assuming
    ideal co-located sensors; they are **not** a reproduction of
    MobileMT's proprietary processed apparent-conductivity output.
    Prefer ``apparent_conductivity_native_Sm`` whenever it is present.
    """
    ds = ensure_mobilemt_dataset(dataset)
    rows: list[dict[str, Any]] = []
    for line_id, sample_id, x_m, record in _dataset_records(ds):
        tf, freq = _get_admittance(record)
        if tf is None or freq is None:
            continue
        Y = np.asarray(tf.data, dtype=complex)
        y2 = Y[:, :2, :]
        det = y2[:, 0, 0] * y2[:, 1, 1] - y2[:, 0, 1] * y2[:, 1, 0]
        with np.errstate(invalid="ignore"):
            y_eff = det**0.5
        sigma_a = 5.0 * freq * np.abs(y_eff) ** 2
        with np.errstate(divide="ignore"):
            rho_a = 1.0 / sigma_a
        phase_a = -np.degrees(np.angle(y_eff))
        sigma_native = _native_sigma(record)

        for k in range(freq.size):
            if not np.isfinite(det[k].real) or not np.isfinite(det[k].imag):
                continue
            rows.append(
                {
                    "line_id": line_id,
                    "sample_id": sample_id,
                    "x_m": x_m,
                    "freq_hz": float(freq[k]),
                    "period_s": 1.0 / max(float(freq[k]), 1e-12),
                    "det_abs": float(np.abs(det[k])),
                    "theoretical_sigma_a_Sm": float(sigma_a[k]),
                    "theoretical_rho_a_ohm_m": float(rho_a[k]),
                    "theoretical_phase_deg": float(phase_a[k]),
                    "apparent_conductivity_native_Sm": (
                        float(sigma_native[k])
                        if sigma_native is not None
                        and k < sigma_native.size
                        and np.isfinite(sigma_native[k])
                        else np.nan
                    ),
                }
            )
    if not rows:
        return pd.DataFrame(columns=_DETERMINANT_COLS)
    return pd.DataFrame(rows, columns=_DETERMINANT_COLS)


def admittance_skew_table(dataset: Any) -> pd.DataFrame:
    r"""Return a Swift (1967)-style skew table for the admittance tensor.

    .. math::

        \mathrm{skew} = \frac{|Y_{xx} + Y_{yy}|}{|Y_{xy} - Y_{yx}|}

    applied to the horizontal 2x2 admittance submatrix by direct
    algebraic analogy to the identical ratio already used for the
    impedance tensor elsewhere in pyCSAMT. Being a ratio of
    magnitudes, it needs no absolute physical constant and is safe to
    compute directly, unlike :func:`admittance_determinant_table`'s
    ``theoretical_*`` columns. Large values flag departures from an
    ideal 1D/2D-consistent admittance tensor (instrument coupling,
    cultural noise, genuinely 3-D structure).

    Parameters
    ----------
    dataset : AirborneEMDataset or AirborneEMLine
        Anything accepted by :func:`ensure_mobilemt_dataset`.

    Returns
    -------
    pandas.DataFrame
        Columns: ``line_id``, ``sample_id``, ``x_m``, ``freq_hz``,
        ``period_s``, ``skew``. Non-finite values are omitted.
    """
    ds = ensure_mobilemt_dataset(dataset)
    rows: list[dict[str, Any]] = []
    for line_id, sample_id, x_m, record in _dataset_records(ds):
        tf, freq = _get_admittance(record)
        if tf is None or freq is None:
            continue
        Y = np.asarray(tf.data, dtype=complex)
        yxx, yxy = Y[:, 0, 0], Y[:, 0, 1]
        yyx, yyy = Y[:, 1, 0], Y[:, 1, 1]
        with np.errstate(divide="ignore", invalid="ignore"):
            skew = np.abs(yxx + yyy) / np.abs(yxy - yyx)
        for k in range(freq.size):
            if not np.isfinite(skew[k]):
                continue
            rows.append(
                {
                    "line_id": line_id,
                    "sample_id": sample_id,
                    "x_m": x_m,
                    "freq_hz": float(freq[k]),
                    "period_s": 1.0 / max(float(freq[k]), 1e-12),
                    "skew": float(skew[k]),
                }
            )
    if not rows:
        return pd.DataFrame(columns=_SKEW_COLS)
    return pd.DataFrame(rows, columns=_SKEW_COLS)


# ─────────────────────────────────────────────────────────────────────────
# Pipeline QC: mask frequencies outside the published MobileMT band
# ─────────────────────────────────────────────────────────────────────────


def mask_outside_mobilemt_band(
    dataset: Any,
    *,
    band_hz: tuple[float, float] | None = None,
    system_spec: MobileMTSystemSpec | None = None,
    inplace: bool = False,
) -> AirborneEMDataset:
    r"""Mask admittance/conductivity outside the usable MobileMT band.

    Reuses the published usable bandwidth already carried by
    :class:`~pycsamt.airborne.mobilemt.MobileMTSystemSpec` (default
    ``nominal_frequency_range_hz`` of 19-26,000 Hz) rather than
    inventing a new band definition. This is the one function in this
    module meant to sit inside a processing pipeline (dataset in,
    dataset out) rather than only produce a diagnostic table -- the
    closest analogue here to
    :func:`~pycsamt.emtools.afmag.flag_motion_susceptible_band` and
    :func:`~pycsamt.emtools.ztem.mask_outside_ztem_band`.

    Unlike those two functions, only masking is offered (no
    ``action="drop"``): each :class:`~pycsamt.airborne.AirborneEMRecord`
    packages its admittance transfer function and any auxiliary
    per-frequency fields (variance, covariances, native apparent
    conductivity) around one shared period axis, and safely dropping
    frequencies would require rebuilding all of them consistently.
    Masking with ``nan`` needs no such reconstruction and never
    confuses "known bad" with a physical zero.

    Parameters
    ----------
    dataset : AirborneEMDataset or AirborneEMLine
        Anything accepted by :func:`ensure_mobilemt_dataset`.
    band_hz : (float, float), optional
        Explicit ``(low, high)`` band in Hz. Mutually exclusive with
        *system_spec*; when neither is given, a default
        :class:`~pycsamt.airborne.mobilemt.MobileMTSystemSpec`'s
        ``nominal_frequency_range_hz`` is used.
    system_spec : MobileMTSystemSpec, optional
        Survey-specific system metadata to read the band from.
    inplace : bool, default False
        When ``False`` (default), a deep copy of *dataset* is masked
        and returned, leaving the input untouched.

    Returns
    -------
    AirborneEMDataset
        The (optionally new) dataset with out-of-band admittance
        values and native apparent-conductivity samples set to
        ``nan``.

    Raises
    ------
    ValueError
        If both *band_hz* and *system_spec* are given.
    TypeError
        If *system_spec* is given and is not a
        :class:`~pycsamt.airborne.mobilemt.MobileMTSystemSpec`.
    """
    if band_hz is not None and system_spec is not None:
        raise ValueError(
            "give at most one of band_hz or system_spec, not both"
        )
    if system_spec is not None and not isinstance(
        system_spec, MobileMTSystemSpec
    ):
        raise TypeError("system_spec must be a MobileMTSystemSpec or None")

    if band_hz is not None:
        lo, hi = float(band_hz[0]), float(band_hz[1])
    else:
        spec = system_spec or MobileMTSystemSpec()
        lo, hi = spec.nominal_frequency_range_hz

    ds = ensure_mobilemt_dataset(dataset)
    if not inplace:
        ds = copy.deepcopy(ds)

    for line in ds.iter_lines():
        for record in line.iter_records():
            tf, freq = _get_admittance(record)
            if tf is None or freq is None:
                continue
            keep = (freq >= lo) & (freq <= hi)
            if keep.all():
                continue
            tf.data[~keep, :, :] = np.nan + 1j * np.nan
            sigma_native = _native_sigma(record)
            if sigma_native is not None:
                sigma_native = sigma_native.copy()
                sigma_native[~keep[: sigma_native.size]] = np.nan
                record.fields[MOBILEMT_APPARENT_CONDUCTIVITY_FIELD] = (
                    sigma_native
                )
    return ds


# ─────────────────────────────────────────────────────────────────────────
# Plots
# ─────────────────────────────────────────────────────────────────────────


def plot_mobilemt_admittance_profile(
    dataset: Any,
    *,
    line_id: str | None = None,
    component: str = "det",
    part: str = "abs",
    frequency_hz: float | None = None,
    period_s: float | None = None,
    figsize: tuple[float, float] = (9.5, 4.0),
    ax: plt.Axes | None = None,
) -> plt.Axes:
    r"""Plot one admittance component along one flight line.

    Parameters
    ----------
    dataset : AirborneEMDataset or AirborneEMLine
        Anything accepted by :func:`ensure_mobilemt_dataset`.
    line_id : str, optional
        Flight line to plot; defaults to the first line in *dataset*.
    component : {"xx", "xy", "yx", "yy", "hzx", "hzy", "det"}, default "det"
        Admittance entry to plot, or ``"det"`` for the horizontal
        2x2 determinant (see :func:`admittance_determinant_table`).
    part : {"real", "imag", "abs"}, default "abs"
    frequency_hz, period_s : float, optional
        Reference frequency/period; nearest available value is used
        per sample. At most one may be given; the median frequency
        is used when neither is given.
    figsize : (float, float), default (9.5, 4.0)
        Used only when *ax* is not supplied.
    ax : matplotlib.axes.Axes, optional
        Existing axes to draw on.

    Returns
    -------
    matplotlib.axes.Axes
    """
    ds = ensure_mobilemt_dataset(dataset)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    try:
        line = _pick_line(ds, line_id)
    except ValueError:
        ax.text(0.5, 0.5, "no flight lines", ha="center", va="center")
        return ax

    pos = _line_chainage(line)
    xs: list[float] = []
    ys: list[float] = []
    freqs_used: list[float] = []
    for idx, sample_id in enumerate(line.navigation.sample_ids):
        record = line.records.get(sample_id)
        if record is None:
            continue
        tf, freq = _get_admittance(record)
        if tf is None or freq is None:
            continue
        k, f0 = _resolve_target_frequency(freq, frequency_hz, period_s)
        val = _component_series(
            np.asarray(tf.data, dtype=complex)[k : k + 1], component
        )[0]
        if not np.isfinite(val.real) or not np.isfinite(val.imag):
            continue
        v = {"real": val.real, "imag": val.imag, "abs": abs(val)}[part]
        xs.append(pos[idx])
        ys.append(float(v))
        freqs_used.append(f0)

    if not xs:
        ax.text(0.5, 0.5, "no admittance data", ha="center", va="center")
        return ax

    order = np.argsort(xs)
    xs_arr = np.asarray(xs)[order]
    ys_arr = np.asarray(ys)[order]
    ref = float(np.nanmedian(freqs_used))

    _ml = PYCSAMT_STYLE.multiline
    ax.axhline(0.0, color="0.85", lw=0.8)
    ax.plot(xs_arr, ys_arr, lw=_ml.lw, alpha=_ml.alpha, color="tab:blue")
    ax.set_xlabel("Position along flight line (m)")
    ax.set_ylabel(f"Y[{component}]  [{part}]")
    ax.set_title(
        f"MobileMT admittance profile — {line.line_id} @ {ref:.4g} Hz",
        fontsize=10,
    )
    ax.grid(True, ls=":", alpha=0.3)
    return ax


def plot_mobilemt_conductivity_psection(
    dataset: Any,
    *,
    line_id: str | None = None,
    source: str = "theoretical",
    cmap: str = "viridis",
    clim: tuple[float, float] | None = None,
    clim_pct: tuple[float, float] = (2.0, 98.0),
    figsize: tuple[float, float] = (9.0, 5.0),
    ax: plt.Axes | None = None,
) -> plt.Axes:
    r"""Plot an apparent-conductivity pseudosection for one flight line.

    Parameters
    ----------
    dataset : AirborneEMDataset or AirborneEMLine
        Anything accepted by :func:`ensure_mobilemt_dataset`.
    line_id : str, optional
        Flight line to plot; defaults to the first line in *dataset*.
    source : {"theoretical", "native"}, default "theoretical"
        ``"theoretical"`` plots
        :func:`admittance_determinant_table`'s derived
        ``theoretical_sigma_a_Sm`` (see the module docstring for the
        caveat); ``"native"`` plots the vendor-delivered
        ``apparent_conductivity_native_Sm`` field, when present.
    cmap : str, default "viridis"
    clim : (float, float), optional
        Explicit color limits; overrides *clim_pct*.
    clim_pct : (float, float), default (2.0, 98.0)
        Percentile color limits when *clim* is not given.
    figsize : (float, float), default (9.0, 5.0)
        Used only when *ax* is not supplied.
    ax : matplotlib.axes.Axes, optional
        Existing axes to draw on.

    Returns
    -------
    matplotlib.axes.Axes

    Raises
    ------
    ValueError
        If *source* is not ``"theoretical"`` or ``"native"``.
    """
    source = str(source).strip().lower()
    if source not in {"theoretical", "native"}:
        raise ValueError(
            f"source must be 'theoretical' or 'native'; got {source!r}"
        )

    ds = ensure_mobilemt_dataset(dataset)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    try:
        line = _pick_line(ds, line_id)
    except ValueError:
        ax.text(0.5, 0.5, "no flight lines", ha="center", va="center")
        return ax

    one_line = AirborneEMDataset(name=line.line_id, lines={line.line_id: line})
    if source == "theoretical":
        df = admittance_determinant_table(one_line)
        col = "theoretical_sigma_a_Sm"
    else:
        df = admittance_table(one_line)
        col = "apparent_conductivity_native_Sm"

    if df.empty or df[col].isna().all():
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        return ax

    order_df = df.drop_duplicates("sample_id").sort_values("x_m")
    samples = order_df["sample_id"].tolist()
    x_sorted = order_df["x_m"].to_numpy()
    freqs = np.sort(df["freq_hz"].unique())[::-1]

    grid = np.full((freqs.size, len(samples)), np.nan)
    s_idx = {s: j for j, s in enumerate(samples)}
    f_idx = {float(f): i for i, f in enumerate(freqs)}
    for _, row in df.iterrows():
        j = s_idx.get(row["sample_id"])
        i = f_idx.get(float(row["freq_hz"]))
        if i is not None and j is not None:
            grid[i, j] = row[col]

    finite = grid[np.isfinite(grid)]
    if clim is None:
        vmin, vmax = (
            np.percentile(finite, clim_pct) if finite.size else (0.0, 1.0)
        )
    else:
        vmin, vmax = clim

    extent = (
        -0.5,
        len(samples) - 0.5,
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
    n_xtick = min(8, len(samples))
    step = max(1, len(samples) // n_xtick)
    xt_idx = np.arange(0, len(samples), step)
    ax.set_xticks(xt_idx)
    ax.set_xticklabels(
        [f"{x_sorted[k]:.0f}" for k in xt_idx],
        rotation=45,
        ha="right",
        fontsize=8,
    )
    ax.set_xlabel("Position along flight line (m)")
    ax.set_ylabel(LOG10_PERIOD_LABEL)
    label = "theoretical" if source == "theoretical" else "native"
    ax.set_title(
        f"MobileMT apparent conductivity [{label}] — {line.line_id}",
        fontsize=10,
    )
    add_colorbar(im, ax, label=r"$\sigma_a$ (S/m)")
    return ax


def plot_mobilemt_skew_profile(
    dataset: Any,
    *,
    line_id: str | None = None,
    frequency_hz: float | None = None,
    period_s: float | None = None,
    figsize: tuple[float, float] = (9.5, 4.0),
    ax: plt.Axes | None = None,
) -> plt.Axes:
    r"""Plot the admittance skew profile along one flight line.

    See :func:`admittance_skew_table` for the underlying formula.

    Parameters
    ----------
    dataset : AirborneEMDataset or AirborneEMLine
        Anything accepted by :func:`ensure_mobilemt_dataset`.
    line_id : str, optional
        Flight line to plot; defaults to the first line in *dataset*.
    frequency_hz, period_s : float, optional
        Reference frequency/period; nearest available value is used
        per sample. At most one may be given; the median frequency
        is used when neither is given.
    figsize : (float, float), default (9.5, 4.0)
        Used only when *ax* is not supplied.
    ax : matplotlib.axes.Axes, optional
        Existing axes to draw on.

    Returns
    -------
    matplotlib.axes.Axes
    """
    ds = ensure_mobilemt_dataset(dataset)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    try:
        line = _pick_line(ds, line_id)
    except ValueError:
        ax.text(0.5, 0.5, "no flight lines", ha="center", va="center")
        return ax

    one_line = AirborneEMDataset(name=line.line_id, lines={line.line_id: line})
    df = admittance_skew_table(one_line)
    if df.empty:
        ax.text(0.5, 0.5, "no admittance data", ha="center", va="center")
        return ax

    freqs = df["freq_hz"].to_numpy()
    _, target = _resolve_target_frequency(
        np.sort(np.unique(freqs)), frequency_hz, period_s
    )
    xs: list[float] = []
    ys: list[float] = []
    for _sample_id, sub in df.groupby("sample_id", sort=False):
        idx = (sub["freq_hz"] - target).abs().idxmin()
        xs.append(float(sub.loc[idx, "x_m"]))
        ys.append(float(sub.loc[idx, "skew"]))

    order = np.argsort(xs)
    xs_arr = np.asarray(xs)[order]
    ys_arr = np.asarray(ys)[order]

    _ml = PYCSAMT_STYLE.multiline
    ax.plot(xs_arr, ys_arr, lw=_ml.lw, alpha=_ml.alpha, color="tab:purple")
    ax.set_xlabel("Position along flight line (m)")
    ax.set_ylabel("Admittance skew")
    ax.set_title(
        f"MobileMT admittance skew — {line.line_id} @ {target:.4g} Hz",
        fontsize=10,
    )
    ax.grid(True, ls=":", alpha=0.3)
    return ax
