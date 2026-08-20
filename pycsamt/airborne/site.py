# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

r"""Airborne station containers: the non-EDI counterpart of ``pycsamt.site``.

:class:`~pycsamt.site.base.Site`/:class:`~pycsamt.site.base.Sites`
are, by design, EDI-shaped: a 2x2 impedance :math:`Z`, a 1x2 tipper,
one flat station per file. Every numeric accessor on ``Site``
(``.tipper``, ``.z``, ``.rho``, ``.phase``) routes through
:attr:`~pycsamt.site.base.Site.edi`, and
:func:`pycsamt.emtf.converters.edi.emtf_to_edi` deliberately refuses
to build an EDI from a tipper-only document (no impedance transfer
function) -- ``pycsamt.airborne.ztem.tests.test_ztem_interop
.test_ztem_is_not_silently_coerced_to_standard_edi`` enforces exactly
this, and ``pycsamt.airborne.mobilemt.tests.test_mobilemt_interop``
gives MobileMT's admittance the identical guarantee. So a genuine
ZTEM/AFMAG measurement -- tipper-only, no electric channel at all --
or a MobileMT measurement -- a 3x2 admittance tensor, not a tipper --
cannot reach ``Site``/``Sites`` today, and forcing either into that
2x2-impedance-shaped model would mean either bending ``Site`` until
it stops meaning "ground EDI station," or leaving MobileMT
permanently unrepresentable regardless.

This module is the parallel container family built instead, directly
on the airborne scientific model that already exists
(:class:`~pycsamt.airborne.AirborneEMRecord`/
:class:`~pycsamt.airborne.AirborneEMLine`/
:class:`~pycsamt.airborne.AirborneEMDataset`, holding a real
:class:`~pycsamt.emtf.document.EMTF` document per sample -- no EDI
bridge, no impedance requirement):

* :class:`AirborneSite` -- one flight-line sample. Read directly
  from its :class:`~pycsamt.emtf.document.EMTF` document:
  ``.tipper``, ``.freq``, and, unlike ``Site``, ``.admittance`` (for
  MobileMT) alongside the typed metadata surface ``Site`` already
  exposes (``.site_meta``, ``.processing``, ``.provenance``, ...),
  since ``AirborneEMRecord.emtf`` *is* an ``EMTF`` document too.
* :class:`AirborneSites` -- an ordered collection, mirroring
  :class:`~pycsamt.site.base.Sites`'s indexing/lookup/select/map/
  closest/write_xml surface where the concepts genuinely carry over.
  Bulk-editing helpers (``edit_all``, ``with_topography``,
  ``to_profile``) are deliberately not mirrored yet -- they are
  real ``Sites`` features this module does not need for the
  ``emtools`` integration it exists to unblock, so they are left out
  rather than stubbed.
* :func:`ensure_asites` -- the ``ensure_sites`` counterpart: the
  single entry point every future airborne-aware ``emtools``
  function should coerce its input through.

Reads EMTF-XML today (:meth:`AirborneSites.from_xml_dir`,
:meth:`AirborneSite.from_xml`); a native H5 reader is a real,
separate task for whenever a verified delivery schema exists (see
:mod:`pycsamt.airborne.io`'s own module docstring for why pyCSAMT
does not guess native formats), not something this module invents
now.

Ordering
--------
:class:`~pycsamt.site.base.Sites` has to *infer* a spatial line
order from station coordinates, because ground stations can arrive
in arbitrary order. Airborne samples do not have that problem --
:class:`~pycsamt.airborne.NavigationTrack` already records the
definitive flight-line order -- so :class:`AirborneSites` has no
``ordered()``/``ordering`` analogue: :meth:`AirborneSites.from_line`
and :meth:`AirborneSites.from_dataset` simply preserve navigation
order, which is already correct.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd

from ..constants import _EARTH_R
from ..core.base import CoreObject
from .base import AirborneEMDataset, AirborneEMLine, AirborneEMRecord
from .navigation import NavigationTrack
from .validation import clean_identifier, emtf_class

if TYPE_CHECKING:
    from ..emtf.document import EMTF

__all__ = ["AirborneSite", "AirborneSites", "ensure_asites"]

_MOBILEMT_ADMITTANCE_TAG = "mobilemt_admittance"
_AFMAG_INTERSTATION_TAG = "interstation_transfer_functions"
_AFMAG_TILT_TAG = "afmag_tilt"
_AFMAG_AP_TAG = "airmt_amplification_parameter"
_VALID_ON_DUP = frozenset({"replace", "keep_first", "keep_last", "raise"})


def _nan3() -> tuple[float, float, float]:
    return (float("nan"), float("nan"), float("nan"))


def _navigation_coords(
    nav: NavigationTrack, index: int
) -> tuple[float, float, float]:
    """Return ``(lat, lon, elev)`` for one navigation sample index."""
    lat = (
        float(nav.latitude[index])
        if nav.latitude is not None
        else float("nan")
    )
    lon = (
        float(nav.longitude[index])
        if nav.longitude is not None
        else float("nan")
    )
    if nav.platform_elevation is not None:
        elev = float(nav.platform_elevation[index])
    elif nav.terrain_elevation is not None:
        elev = float(nav.terrain_elevation[index])
    else:
        elev = float("nan")
    return (lat, lon, elev)


def _haversine_m(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Great-circle distance [m] between two lat/lon points.

    Self-contained on purpose: :func:`pycsamt.site.location.distance`
    would work too, but importing it pulls in ``pycsamt.site``'s own
    EDI-oriented helpers, which this module deliberately stays free
    of (see the module docstring).
    """
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dphi = np.radians(lat2 - lat1)
    dlmb = np.radians(lon2 - lon1)
    a = (
        np.sin(dphi / 2.0) ** 2
        + np.cos(p1) * np.cos(p2) * np.sin(dlmb / 2.0) ** 2
    )
    return float(2.0 * _EARTH_R * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0))))


def _stem_from_source(source: Any) -> str | None:
    if isinstance(source, (str, Path)):
        stem = Path(source).stem
        return stem or None
    return None


class AirborneSite(CoreObject):
    r"""One airborne flight-line sample, read directly from its EMTF.

    Unlike :class:`~pycsamt.site.base.Site`, every numeric accessor
    here reads straight from the wrapped
    :class:`~pycsamt.airborne.AirborneEMRecord`'s
    :class:`~pycsamt.emtf.document.EMTF` document -- there is no EDI
    bridge to materialize and no impedance requirement to satisfy.

    Parameters
    ----------
    record : AirborneEMRecord
        The sample this site wraps.
    line_id : str, optional
        Owning flight-line identifier, for provenance/grouping.
        ``None`` when unknown (e.g. a bare, line-less record).
    technology : str, optional
        Technology label (``"ztem"``, ``"mobilemt"``, ``"afmag"``,
        ...), the explicit differentiator against ground MT. Falls
        back to the record's ``EMTF.subtype`` when not given
        explicitly; see :attr:`technology`.
    coords : (float, float, float), optional
        Explicit ``(lat, lon, elev)`` override, normally supplied by
        the constructing classmethod
        (:meth:`AirborneSites.from_line`/``from_dataset``) from the
        parent line's :class:`~pycsamt.airborne.NavigationTrack`.
        When omitted, :attr:`coords` falls back to the record's own
        ``EMTF.site.location``, then to ``(nan, nan, nan)``.

    Raises
    ------
    TypeError
        If *record* is not an :class:`AirborneEMRecord`.

    See Also
    --------
    AirborneSites : Ordered collection of these.
    pycsamt.site.base.Site : The ground-MT, EDI-shaped counterpart.
    """

    def __init__(
        self,
        record: AirborneEMRecord,
        *,
        line_id: str | None = None,
        technology: str | None = None,
        coords: tuple[float, float, float] | None = None,
    ) -> None:
        if not isinstance(record, AirborneEMRecord):
            raise TypeError("record must be an AirborneEMRecord")
        self._record = record
        self._line_id = (
            None if line_id is None else clean_identifier(
                line_id, name="line_id",
            )
        )
        self._technology = (
            str(technology).strip().lower() if technology else None
        )
        self._coords = coords

    @classmethod
    def from_xml(
        cls,
        source: Any,
        *,
        line_id: str | None = None,
        technology: str | None = None,
    ) -> AirborneSite:
        r"""Build one site directly from an EMTF-XML file or document.

        Parameters
        ----------
        source : pycsamt.emtf.document.EMTF or str or pathlib.Path
            Parsed document, or a path read via
            :meth:`~pycsamt.emtf.document.EMTF.from_xml`.
        line_id, technology : str, optional
            Forwarded to the constructor.

        Returns
        -------
        AirborneSite

        Notes
        -----
        The sample identifier is resolved from
        ``document.site.site_id``, then ``document.station``, then
        the source file's stem, in that order.
        """
        EMTF = emtf_class()
        doc = source if isinstance(source, EMTF) else EMTF.from_xml(source)
        sample_id = None
        if doc.site is not None and doc.site.site_id:
            sample_id = str(doc.site.site_id)
        elif doc.station:
            sample_id = str(doc.station)
        else:
            sample_id = _stem_from_source(source)
        record = AirborneEMRecord(
            sample_id=sample_id or "site", emtf=doc,
        )
        return cls(record, line_id=line_id, technology=technology)

    @property
    def record(self) -> AirborneEMRecord:
        """The wrapped :class:`~pycsamt.airborne.AirborneEMRecord`."""
        return self._record

    @property
    def emtf(self) -> EMTF | None:
        """The underlying :class:`~pycsamt.emtf.document.EMTF`."""
        return self._record.emtf

    @property
    def tf(self) -> EMTF | None:
        """Alias for :attr:`emtf`, matching ``Site.tf``'s naming."""
        return self._record.emtf

    @property
    def sample_id(self) -> str:
        """Navigation sample identifier (stable; independent of
        any richer :attr:`name` resolved from metadata)."""
        return self._record.sample_id

    @property
    def line_id(self) -> str | None:
        """Owning flight-line identifier, or ``None`` if unknown."""
        return self._line_id

    @property
    def technology(self) -> str | None:
        """Technology label (``"ztem"``, ``"mobilemt"``, ...).

        Explicit at construction, else the record's ``EMTF.subtype``,
        else ``None``.
        """
        if self._technology:
            return self._technology
        doc = self._record.emtf
        if doc is not None and doc.subtype:
            return str(doc.subtype).strip().lower()
        return None

    @property
    def name(self) -> str:
        r"""Station identifier resolved from metadata, else
        :attr:`sample_id`.

        Resolution order: ``EMTF.site.site_id``, ``EMTF.station``,
        :attr:`sample_id`.
        """
        doc = self._record.emtf
        if doc is not None:
            if doc.site is not None and doc.site.site_id:
                return str(doc.site.site_id)
            if doc.station:
                return str(doc.station)
        return self._record.sample_id

    @property
    def station(self) -> str:
        """Alias for :attr:`name`, matching ``Site.name``'s role."""
        return self.name

    @property
    def coords(self) -> tuple[float, float, float]:
        r"""``(lat, lon, elev)`` in decimal degrees and metres.

        Explicit constructor override first, then the record's own
        ``EMTF.site.location``, then ``(nan, nan, nan)`` -- never a
        fabricated position.
        """
        if self._coords is not None:
            return self._coords
        doc = self._record.emtf
        if doc is not None and doc.site is not None:
            loc = doc.site.location
            if loc is not None:
                lat = loc.latitude if loc.latitude is not None else np.nan
                lon = loc.longitude if loc.longitude is not None else np.nan
                elev = (
                    loc.elevation if loc.elevation is not None else np.nan
                )
                return (float(lat), float(lon), float(elev))
        return _nan3()

    @property
    def freq(self) -> np.ndarray | None:
        """Frequency vector [Hz], or ``None`` if unknown."""
        doc = self._record.emtf
        return None if doc is None else doc.frequency

    @property
    def tipper(self) -> np.ndarray | None:
        """Tipper array, shape ``(nf, 1, 2)``, or ``None`` if absent."""
        doc = self._record.emtf
        return None if doc is None else doc.tipper

    @tipper.setter
    def tipper(self, value: Any) -> None:
        r"""Write back into the underlying tipper transfer function.

        Accepts either the canonical ``(nf, 1, 2)`` shape or the
        ``(nf, 2)`` form (broadcast to the former). This exists so
        in-place QC edits (masking a frequency band, for example)
        write straight through to :attr:`emtf` -- mirroring what
        ground EDI-backed ``Site`` objects already allow via their
        mutable ``Tip.tipper`` attribute.

        Raises
        ------
        ValueError
            If this site has no ``EMTF`` document, or the document
            has no tipper transfer function to write into.
        """
        doc = self._record.emtf
        if doc is None:
            raise ValueError(f"{self.name!r} has no EMTF document")
        tf = doc.tipper_tf
        if tf is None:
            raise ValueError(
                f"{self.name!r} has no tipper transfer function"
            )
        arr = np.asarray(value, dtype=complex)
        if arr.ndim == 2 and arr.shape[-1] == 2:
            arr = arr[:, None, :]
        tf.data = arr

    @property
    def z(self) -> np.ndarray | None:
        """Impedance array, shape ``(nf, 2, 2)``, or ``None`` if
        absent (the normal case for ZTEM/AFMAG/MobileMT)."""
        doc = self._record.emtf
        return None if doc is None else doc.z

    @property
    def admittance(self) -> np.ndarray | None:
        r"""MobileMT admittance array, shape ``(nf, 3, 2)``.

        ``None`` for any record without an attached
        ``mobilemt_admittance`` transfer function -- there is no
        analogue of this accessor on ``Site``, since a 3x2 admittance
        cannot be represented by ``Site.z``'s 2x2 shape.
        """
        doc = self._record.emtf
        if doc is None:
            return None
        tf = doc.get_transfer_function(_MOBILEMT_ADMITTANCE_TAG)
        return None if tf is None else tf.data

    @property
    def interstation_tensor(self) -> np.ndarray | None:
        r"""Tensor AFMAG/AirMt interstation magnetic TF, shape
        ``(nf, 3, 2)``.

        ``None`` for any record without an attached
        ``interstation_transfer_functions`` transfer function.
        Deliberately a separate accessor from :attr:`admittance`
        (also ``(nf, 3, 2)``-shaped): AirMt's tensor relates
        ground-reference *magnetic* fields to airborne magnetic
        fields (Hx,Hy -> Hx,Hy,Hz), while MobileMT's admittance
        relates ground *electric* fields to airborne magnetic fields
        (Ex,Ey -> Hx,Hy,Hz) -- physically different responses that
        happen to share a matrix shape, so this module keeps them
        under different names rather than one generic "3x2 tensor"
        accessor that would blur which is which. See
        :mod:`pycsamt.airborne.afmag`'s module docstring for why the
        original-comparator generation (:attr:`afmag_tilt_deg`) is
        kept separate again from this one.
        """
        doc = self._record.emtf
        if doc is None:
            return None
        tf = doc.get_transfer_function(_AFMAG_INTERSTATION_TAG)
        return None if tf is None else tf.data

    @property
    def afmag_tilt_deg(self) -> np.ndarray | None:
        r"""Original comparator AFMAG scalar tilt/deflection, shape
        ``(nf,)``.

        ``None`` for any record without an attached ``afmag_tilt``
        transfer function. Unlike :attr:`interstation_tensor` or a
        ground tipper, this is a single real number per frequency --
        the historical comparator has no polarization-ellipse
        decomposition available at all, only a line-direction
        deflection (see :mod:`pycsamt.airborne.afmag`'s module
        docstring).
        """
        doc = self._record.emtf
        if doc is None:
            return None
        tf = doc.get_transfer_function(_AFMAG_TILT_TAG)
        if tf is None:
            return None
        return np.asarray(tf.data, dtype=complex).real.reshape(-1)

    @property
    def afmag_amplification_parameter(self) -> np.ndarray | None:
        r"""AirMt rotation-invariant amplification parameter, shape
        ``(nf,)``.

        ``None`` for any record without an attached
        ``airmt_amplification_parameter`` transfer function; see
        :func:`pycsamt.airborne.afmag.compute_airmt_amplification_parameter`
        for the formula that derives it from
        :attr:`interstation_tensor`.
        """
        doc = self._record.emtf
        if doc is None:
            return None
        tf = doc.get_transfer_function(_AFMAG_AP_TAG)
        if tf is None:
            return None
        return np.asarray(tf.data, dtype=complex).reshape(-1)

    @property
    def quality(self) -> dict[str, Any]:
        """Sample-level quality flags/scores (technology-defined)."""
        return dict(self._record.quality)

    @property
    def fields(self) -> dict[str, Any]:
        r"""Auxiliary decoded fields (technology-defined).

        For MobileMT, this is where a vendor-delivered native
        ``apparent_conductivity`` vector lives when present; see
        :data:`pycsamt.airborne.mobilemt.MOBILEMT_APPARENT_CONDUCTIVITY_FIELD`.
        """
        return dict(self._record.fields)

    @property
    def site_meta(self) -> Any:
        """``pycsamt.metadata.SiteMeta`` (via :attr:`emtf`)."""
        doc = self._record.emtf
        return None if doc is None else doc.site

    @property
    def site_layout(self) -> Any:
        """``pycsamt.metadata.SiteLayout`` channel geometry."""
        doc = self._record.emtf
        return None if doc is None else doc.site_layout

    @property
    def provenance(self) -> Any:
        """``pycsamt.metadata.ProvenanceMeta`` creator/submitter info."""
        doc = self._record.emtf
        return None if doc is None else doc.provenance

    @property
    def processing(self) -> Any:
        """``pycsamt.metadata.ProcessingMeta`` processing/software info."""
        doc = self._record.emtf
        return None if doc is None else doc.processing

    @property
    def copyright(self) -> Any:
        """``pycsamt.metadata.CopyrightInfo`` release/citation info."""
        doc = self._record.emtf
        return None if doc is None else doc.copyright

    @property
    def quality_meta(self) -> Any:
        """``pycsamt.metadata.TransferFunctionQuality`` QC rating."""
        doc = self._record.emtf
        return None if doc is None else doc.quality

    def has_component(self, comp: str) -> bool:
        r"""Whether *comp* exists and has at least one finite value.

        Parameters
        ----------
        comp : str
            ``"tip"``/``"tx"``/``"ty"``/``"tipper"`` for the tipper;
            ``"admittance"``/``"y"`` for the MobileMT admittance;
            ``"interstation_tensor"``/``"ti"`` for the AirMt tensor;
            ``"afmag_tilt"``/``"tilt"`` for the original-comparator
            AFMAG scalar; ``"amplification_parameter"``/``"ap"`` for
            the AirMt derived parameter; anything else is looked up
            against :attr:`z`.

        Returns
        -------
        bool
        """
        c = str(comp).strip().lower()
        if c in ("tip", "tx", "ty", "tipper"):
            arr = self.tipper
        elif c in ("admittance", "y"):
            arr = self.admittance
        elif c in ("interstation_tensor", "ti"):
            arr = self.interstation_tensor
        elif c in ("afmag_tilt", "tilt"):
            arr = self.afmag_tilt_deg
        elif c in ("amplification_parameter", "ap"):
            arr = self.afmag_amplification_parameter
        else:
            arr = self.z
        if arr is None:
            return False
        a = np.asarray(arr)
        return a.size > 0 and bool(np.any(np.isfinite(a)))

    def to_dataframe(self, kind: str = "tipper") -> pd.DataFrame:
        r"""Export this site's data to a tidy :class:`pandas.DataFrame`.

        Parameters
        ----------
        kind : {"tipper", "admittance", "z"}, default "tipper"

        Returns
        -------
        pandas.DataFrame
            Indexed by frequency (name ``"f"``). Columns depend on
            *kind*: ``Tx, Ty``; ``Yxx, Yxy, Yyx, Yyy, Yhzx, Yhzy``;
            or ``Zxx, Zxy, Zyx, Zyy``.

        Raises
        ------
        ValueError
            If *kind* is not recognized.
        """
        k = str(kind).strip().lower()
        if k in ("tip", "tipper", "t"):
            arr, cols = self.tipper, ("Tx", "Ty")
        elif k in ("admittance", "y"):
            arr = self.admittance
            cols = ("Yxx", "Yxy", "Yyx", "Yyy", "Yhzx", "Yhzy")
        elif k in ("z", "imp", "impedance"):
            arr, cols = self.z, ("Zxx", "Zxy", "Zyx", "Zyy")
        else:
            raise ValueError(f"Unknown kind: {kind!r}")

        freq = self.freq
        idx = pd.Index(
            np.asarray(freq, dtype=float)
            if freq is not None
            else np.array([], dtype=float),
            name="f",
        )
        if arr is None:
            return pd.DataFrame(index=idx)
        flat = np.asarray(arr).reshape(np.asarray(arr).shape[0], -1)
        data: dict[str, Any] = {}
        for i, c in enumerate(cols):
            data[c] = (
                flat[:, i]
                if i < flat.shape[1]
                else np.full(idx.size, np.nan)
            )
        return pd.DataFrame(data, index=idx)

    def summary(self) -> dict[str, Any]:
        r"""Summarize identity, geometry, and data coverage.

        Returns
        -------
        dict
            Keys: ``name``, ``line_id``, ``sample_id``,
            ``technology``, ``nfreq``, ``lat``, ``lon``, ``elev``,
            ``tipper``, ``admittance`` (booleans).
        """
        lat, lon, elev = self.coords
        freq = self.freq
        nfreq = 0 if freq is None else int(np.asarray(freq).size)
        return {
            "name": self.name,
            "line_id": self.line_id,
            "sample_id": self.sample_id,
            "technology": self.technology,
            "nfreq": nfreq,
            "lat": lat,
            "lon": lon,
            "elev": elev,
            "tipper": self.has_component("tipper"),
            "admittance": self.has_component("admittance"),
        }

    def to_xml(self, target: Any = None, **kwargs: Any) -> Any:
        r"""Serialize this site's document to EMTF XML.

        Parameters
        ----------
        target : str or pathlib.Path, optional
            Destination path. If ``None``, the XML is returned as a
            string.
        **kwargs
            Forwarded to
            :meth:`~pycsamt.emtf.document.EMTF.write_xml`/``to_xml``.

        Returns
        -------
        str or Any

        Raises
        ------
        ValueError
            If this site has no attached ``EMTF`` document.
        """
        doc = self._record.emtf
        if doc is None:
            raise ValueError(f"{self.name!r} has no EMTF document")
        if target is None:
            return doc.to_xml(**kwargs)
        return doc.write_xml(target, **kwargs)

    def __repr__(self) -> str:
        s = self.summary()
        return (
            f"AirborneSite(name={s['name']!r}, "
            f"technology={s['technology']!r}, nfreq={s['nfreq']}, "
            f"coords=({s['lat']:.5f},{s['lon']:.5f},{s['elev']:.1f}))"
        )


def _coerce_one(
    item: Any, *, line_id: str | None = None
) -> AirborneSite:
    """Coerce one heterogeneous item to an :class:`AirborneSite`."""
    if isinstance(item, AirborneSite):
        return item
    if isinstance(item, AirborneEMRecord):
        return AirborneSite(item, line_id=line_id)
    if isinstance(item, (str, Path)):
        return AirborneSite.from_xml(item, line_id=line_id)
    EMTF = emtf_class()
    if isinstance(item, EMTF):
        return AirborneSite.from_xml(item, line_id=line_id)
    raise TypeError(
        "cannot coerce to AirborneSite: "
        f"{type(item).__name__}; pass an AirborneSite, "
        "AirborneEMRecord, EMTF, or a path to an EMTF-XML file"
    )


def _is_single_item(x: Any) -> bool:
    """Whether *x* is one item rather than an iterable of items."""
    if isinstance(x, (AirborneSite, AirborneEMRecord, str, Path)):
        return True
    return not hasattr(x, "__iter__")


def _dedupe(sites: AirborneSites, *, on_dup: str) -> AirborneSites:
    r"""Apply a duplicate-name policy over an :class:`AirborneSites`.

    Mirrors :func:`pycsamt.site.base.to_sites`'s ``on_dup`` contract.
    """
    policy = str(on_dup).strip().lower()
    if policy not in _VALID_ON_DUP:
        raise ValueError(f"on_dup must be one of {sorted(_VALID_ON_DUP)}")
    seen: dict[str, AirborneSite] = {}
    order: list[str] = []
    for s in sites:
        key = s.name.lower()
        if key not in seen:
            order.append(key)
            seen[key] = s
            continue
        if policy == "raise":
            raise ValueError(f"duplicate airborne site name: {s.name!r}")
        if policy == "keep_first":
            continue
        # "keep_last" and "replace" both take the newest occurrence.
        seen[key] = s
    return AirborneSites([seen[k] for k in order])


class AirborneSites(CoreObject):
    r"""Ordered collection of :class:`AirborneSite` objects.

    The airborne counterpart of :class:`~pycsamt.site.base.Sites`.
    See the module docstring for why order here is simply preserved
    navigation order rather than something to infer.

    Parameters
    ----------
    items : AirborneSite, AirborneEMRecord, EMTF, str, Path, or iterable
        A single item, or an iterable of them. Raw
        :class:`~pycsamt.airborne.AirborneEMRecord`/``EMTF``/path
        items are coerced via :meth:`AirborneSite.from_xml`-style
        construction with no line context; prefer
        :meth:`from_line`/:meth:`from_dataset` when that context is
        available.

    Raises
    ------
    TypeError
        If an item cannot be coerced to :class:`AirborneSite`.

    See Also
    --------
    ensure_asites : Flexible entry-point coercion, including
        directories of EMTF-XML files and duplicate-name policy.
    pycsamt.site.base.Sites : The ground-MT counterpart.
    """

    def __init__(self, items: Any) -> None:
        if _is_single_item(items):
            items = [items]
        self._items: list[AirborneSite] = [
            _coerce_one(it) for it in list(items)
        ]

    # -- construction -----------------------------------------------

    @classmethod
    def from_xml_dir(
        cls,
        path: str | Path,
        *,
        recursive: bool = True,
        pattern: str = "*.xml",
        line_id: str | None = None,
        strict: bool = False,
    ) -> AirborneSites:
        r"""Read every EMTF-XML file under *path* into one container.

        Parameters
        ----------
        path : str or pathlib.Path
            A single EMTF-XML file, or a directory to search.
        recursive : bool, default True
            Search subdirectories too (``Path.rglob``) instead of
            only the top level (``Path.glob``).
        pattern : str, default "\*.xml"
            Glob pattern used when *path* is a directory.
        line_id : str, optional
            Forwarded to every :meth:`AirborneSite.from_xml` call.
        strict : bool, default False
            If ``True``, a file that fails to parse raises instead
            of being skipped, and an empty result raises too.

        Returns
        -------
        AirborneSites
            Sites in sorted-filename order.

        Raises
        ------
        ValueError
            If *strict* and nothing could be read.
        """
        p = Path(path)
        files = (
            [p]
            if p.is_file()
            else sorted(p.rglob(pattern) if recursive else p.glob(pattern))
        )
        items: list[AirborneSite] = []
        for f in files:
            try:
                items.append(AirborneSite.from_xml(f, line_id=line_id))
            except Exception:
                if strict:
                    raise
                continue
        if strict and not items:
            raise ValueError(
                f"no EMTF-XML files could be read from {path!r}"
            )
        return cls(items)

    @classmethod
    def from_line(
        cls,
        line: AirborneEMLine,
        *,
        technology: str | None = None,
    ) -> AirborneSites:
        r"""Build a container from one already-constructed flight line.

        Parameters
        ----------
        line : AirborneEMLine
            Records are visited via
            :meth:`~pycsamt.airborne.AirborneEMLine.iter_records`,
            i.e. in navigation order, skipping samples with no
            attached record.
        technology : str, optional
            Forwarded to every :class:`AirborneSite`; falls back to
            ``line.attrs["technology"]`` when not given.

        Returns
        -------
        AirborneSites

        Raises
        ------
        TypeError
            If *line* is not an :class:`AirborneEMLine`.
        """
        if not isinstance(line, AirborneEMLine):
            raise TypeError("line must be an AirborneEMLine")
        tech = technology or line.attrs.get("technology")
        nav = line.navigation
        items = [
            AirborneSite(
                record,
                line_id=line.line_id,
                technology=tech,
                coords=_navigation_coords(
                    nav, nav.index_of(record.sample_id)
                ),
            )
            for record in line.iter_records()
        ]
        return cls(items)

    @classmethod
    def from_dataset(
        cls,
        dataset: AirborneEMDataset,
        *,
        technology: str | None = None,
    ) -> AirborneSites:
        r"""Flatten every line of a dataset into one container.

        Parameters
        ----------
        dataset : AirborneEMDataset
            Lines are visited via
            :meth:`~pycsamt.airborne.AirborneEMDataset.iter_lines`
            in insertion order; see :meth:`from_line` for the
            per-line ordering.
        technology : str, optional
            Forwarded to :meth:`from_line` for every line; falls
            back to ``dataset.attrs["technology"]`` when not given.

        Returns
        -------
        AirborneSites

        Raises
        ------
        TypeError
            If *dataset* is not an :class:`AirborneEMDataset`.
        """
        if not isinstance(dataset, AirborneEMDataset):
            raise TypeError("dataset must be an AirborneEMDataset")
        tech = technology or dataset.attrs.get("technology")
        items: list[AirborneSite] = []
        for line in dataset.iter_lines():
            items.extend(cls.from_line(line, technology=tech).as_list())
        return cls(items)

    # -- container protocol ------------------------------------------

    def __len__(self) -> int:
        return len(self._items)

    def __iter__(self):
        return iter(self._items)

    def __getitem__(self, key: int | str) -> AirborneSite:
        r"""Retrieve by zero-based index or case-insensitive name.

        Raises
        ------
        KeyError
            If *key* is a name and no site matches.
        """
        if isinstance(key, int):
            return self._items[key]
        nm = str(key).lower()
        for s in self._items:
            if s.name.lower() == nm:
                return s
        raise KeyError(key)

    def by_index(self, i: int) -> AirborneSite:
        """Retrieve by zero-based index."""
        return self._items[i]

    def get(self, name: str) -> AirborneSite | None:
        """Safe lookup by case-insensitive name; ``None`` if absent."""
        try:
            return self[name]
        except Exception:
            return None

    def as_list(self) -> list[AirborneSite]:
        """The underlying list of :class:`AirborneSite` objects."""
        return list(self._items)

    def to_emtf_list(self) -> list[Any]:
        """The underlying :class:`~pycsamt.emtf.document.EMTF`
        documents, in site order."""
        return [s.emtf for s in self._items]

    @property
    def technologies(self) -> tuple[str, ...]:
        """Sorted, deduplicated :attr:`AirborneSite.technology` values
        present in this container."""
        return tuple(
            sorted({s.technology for s in self._items if s.technology})
        )

    @property
    def line_ids(self) -> tuple[str, ...]:
        """Sorted, deduplicated :attr:`AirborneSite.line_id` values
        present in this container."""
        return tuple(sorted({s.line_id for s in self._items if s.line_id}))

    # -- selection / mapping ------------------------------------------

    def select(
        self,
        names: Any = None,
        predicate: Any = None,
    ) -> AirborneSites:
        r"""Filter by explicit names or by a boolean predicate.

        Parameters
        ----------
        names : sequence of str, optional
            Case-insensitive names to retain; takes precedence over
            *predicate*.
        predicate : callable, optional
            ``predicate(site) -> bool``.

        Returns
        -------
        AirborneSites
            A new container; a shallow copy when neither argument is
            given.
        """
        if names:
            wanted = {str(n).lower() for n in names}
            out = [s for s in self._items if s.name.lower() in wanted]
        elif predicate:
            out = [s for s in self._items if predicate(s)]
        else:
            out = list(self._items)
        return AirborneSites(out)

    def map(self, fn: Any) -> list[Any]:
        """Apply ``fn(site) -> Any`` to every site; collect results."""
        return [fn(s) for s in self._items]

    def closest(
        self,
        lat: float,
        lon: float,
        tol: float | None = None,
    ) -> AirborneSite | None:
        r"""Nearest site to a target coordinate (great-circle distance).

        Parameters
        ----------
        lat, lon : float
            Target location in decimal degrees.
        tol : float, optional
            Maximum allowed distance in metres; farther than that
            returns ``None``.

        Returns
        -------
        AirborneSite or None
            ``None`` if every site lacks finite coordinates, or the
            nearest is farther than *tol*.
        """
        best: AirborneSite | None = None
        best_d = float("inf")
        for s in self._items:
            sla, slo, _ = s.coords
            if not (np.isfinite(sla) and np.isfinite(slo)):
                continue
            d = _haversine_m(lat, lon, sla, slo)
            if d < best_d:
                best, best_d = s, d
        if best is None:
            return None
        if tol is not None and best_d > tol:
            return None
        return best

    # -- persistence ---------------------------------------------------

    def write_xml(self, outdir: str | Path, **kwargs: Any) -> list[Path]:
        r"""Write one EMTF-XML file per site into a directory.

        Parameters
        ----------
        outdir : str or pathlib.Path
            Destination directory; created if missing.
        **kwargs
            Forwarded to :meth:`AirborneSite.to_xml` for each site.

        Returns
        -------
        list of pathlib.Path
            Paths written, named ``"{name}.xml"``.
        """
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        paths: list[Path] = []
        for s in self._items:
            p = outdir / f"{s.name}.xml"
            s.to_xml(p, **kwargs)
            paths.append(p)
        return paths

    def __repr__(self) -> str:
        techs = ",".join(self.technologies) or "?"
        return f"AirborneSites(n={len(self)}, technologies={techs!r})"


def ensure_asites(
    obj: Any,
    *,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> AirborneSites:
    r"""Normalize arbitrary airborne input to an :class:`AirborneSites`.

    The single entry-point coercion for airborne-aware ``emtools``
    functions, mirroring the role
    :func:`~pycsamt.emtools._core.ensure_sites` plays for the rest of
    ``emtools``.

    Parameters
    ----------
    obj : Any
        Accepts: an existing :class:`AirborneSites`; an
        :class:`~pycsamt.airborne.AirborneEMDataset`
        (:meth:`AirborneSites.from_dataset`); an
        :class:`~pycsamt.airborne.AirborneEMLine`
        (:meth:`AirborneSites.from_line`); a path to a single
        EMTF-XML file or a directory of them
        (:meth:`AirborneSites.from_xml_dir`); or an iterable mixing
        any of the above with bare
        :class:`~pycsamt.airborne.AirborneEMRecord`/``EMTF``/path
        items (each coerced with no line context).
    recursive : bool, default True
        Forwarded to :meth:`AirborneSites.from_xml_dir` for any
        directory encountered.
    on_dup : {"replace", "keep_first", "keep_last", "raise"}, default "replace"
        Duplicate-name policy; see :func:`pycsamt.site.base.to_sites`
        for the identical semantics on the ground side.
    strict : bool, default False
        If ``True``, raise when nothing can be resolved (or, for a
        directory, when no file parses).
    verbose : int, default 0
        ``>0`` warns when the result is empty and *strict* is
        ``False``.

    Returns
    -------
    AirborneSites

    Raises
    ------
    ValueError
        If *obj* is ``None``; if *on_dup* is invalid; or, in
        *strict* mode, if nothing could be resolved.
    """
    if obj is None:
        raise ValueError(
            "ensure_asites: got None. Provide a path, AirborneEMLine, "
            "AirborneEMDataset, AirborneSites, or an iterable of "
            "AirborneSite-like items."
        )

    if isinstance(obj, AirborneSites):
        # Matches to_sites()'s identity contract: an already-coerced
        # container is returned unchanged, not silently rebuilt.
        if str(on_dup).strip().lower() not in _VALID_ON_DUP:
            raise ValueError(f"on_dup must be one of {sorted(_VALID_ON_DUP)}")
        if len(obj) == 0 and strict:
            raise ValueError(
                "ensure_asites(strict=True): no airborne sites were "
                "resolved from the given input."
            )
        return obj
    elif isinstance(obj, AirborneEMDataset):
        result = AirborneSites.from_dataset(obj)
    elif isinstance(obj, AirborneEMLine):
        result = AirborneSites.from_line(obj)
    elif isinstance(obj, (str, Path)):
        result = AirborneSites.from_xml_dir(
            obj, recursive=recursive, strict=strict,
        )
    else:
        items: list[AirborneSite] = []
        candidates = list(obj) if hasattr(obj, "__iter__") else [obj]
        for it in candidates:
            if isinstance(it, AirborneEMDataset):
                items.extend(AirborneSites.from_dataset(it).as_list())
            elif isinstance(it, AirborneEMLine):
                items.extend(AirborneSites.from_line(it).as_list())
            elif isinstance(it, (str, Path)) and Path(it).is_dir():
                items.extend(
                    AirborneSites.from_xml_dir(
                        it, recursive=recursive, strict=False,
                    ).as_list()
                )
            else:
                try:
                    items.append(_coerce_one(it))
                except TypeError:
                    if strict:
                        raise
                    continue
        result = AirborneSites(items)

    result = _dedupe(result, on_dup=on_dup)

    if len(result) == 0:
        if strict:
            raise ValueError(
                "ensure_asites(strict=True): no airborne sites were "
                "resolved from the given input."
            )
        if verbose > 0:
            import warnings

            warnings.warn(
                "ensure_asites: no airborne sites were resolved from "
                "the given input.",
                RuntimeWarning,
                stacklevel=2,
            )
    return result
