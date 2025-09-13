# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

from typing import Any, Optional
from pathlib import Path
import numpy as np

from ..seg.collection import EDICollection
from ..jones.collection import JCollection
from ..seg.edi import EDIFile
from ..jones.j import JFile
from ..zonge.avg import AVG

from ._transformers import TransformerMixin
from .base import TFBundle, ensure_station
from .config import get_config


__all__ = [
    "AVGtoEDI",
    "JtoEDI", 
]


class AVGtoEDI(TransformerMixin):
    r"""
    Convert a Zonge ``AVG`` source to an ``EDICollection``.
    
    This transformer orchestrates extraction of transfer
    functions from a Zonge :class:`~pycsamt.zonge.avg.AVG`
    object (or file path), finalizes the neutral payload
    (:class:`~pycsamt.core.base.TFBundle`), and emits EDI
    objects. The result is always an
    :class:`~pycsamt.seg.collection.EDICollection`, even for
    single-site inputs.
    
    The pipeline is:
    
    1. Parse TFs from the AVG using ``to_tensor``. When ``Z`` is
       not present, fall back to ``(rho, phase)``.
    2. Finalize each bundle using
       :class:`~pycsamt.core._transformers.TransformerMixin`:
       validate the station name, order and de-duplicate
       frequencies, and fill missing parts if enabled by the
       global config.
    3. Materialize :class:`~pycsamt.seg.edi.EDIFile` objects and
       attach site metadata (e.g., coordinates) if available
       from ``avg.topo.frame``.
    
    Notes
    -----
    * Station naming obeys the global policy from
      :mod:`pycsamt.core.config`.
    * Frequency sorting and de-duplication follow ``freq_order``
      and ``freq_tol``.
    * If only ``(rho, phase)`` are present, ``Z`` can be
      reconstructed when ``compute_z_from_res`` is enabled.
    * If ``avg.topo`` exposes a DataFrame-like ``frame`` with
      columns ``station``, ``latitude``, ``longitude`` and
      optionally ``elevation``, a ``>HEAD`` section is updated.
    
    See Also
    --------
    pycsamt.core._transformers.TransformerMixin
        Provides finalization steps and overridable hooks.
    pycsamt.core.base.TFBundle
        Neutral, backend-agnostic payload.
    pycsamt.seg.edi.EDIFile
        Concrete EDI object.
    pycsamt.seg.collection.EDICollection
        Collection wrapper for multiple EDI files.
    
    Examples
    --------
    >>> from pycsamt.core.transformers import AVGtoEDI
    >>> # path or AVG object are both accepted
    >>> # out = AVGtoEDI().transform(\"/path/to/data.avg\")  # doctest:+SKIP
    """

    class _HeadStub:
        r"""
        Private header stub used to seed the EDI ``>HEAD`` section.
        
        Parameters
        ----------
        dataid : str
            EDI data identifier (often the station name).
        lat, lon, elev : float or None
            Optional site coordinates and elevation.
        empty : float or None
            Sentinel for empty numeric values.
        
        Notes
        -----
        Small convenience container; semantics are established by
        the :mod:`pycsamt.seg.edi` implementation.
        """

        def __init__(
            self,
            dataid: str,
            *,
            lat: float | None = None,
            lon: float | None = None,
            elev: float | None = None,
            empty: float | None = None,
        ) -> None:
            self.dataid = dataid
            self.lat = lat
            self.long = lon
            self.elev = elev
            if empty is not None:
                self.empty = float(empty)

    def _as_avg(self, src: Any) -> Any:

        if isinstance(src, AVG):
            return src
        if isinstance(src, (str, Path)):
            return AVG.from_file(src)
        raise TypeError("source must be AVG or path")

    def _empty(self) -> float:
        return float(get_config().empty)

    def compute_z_from_res(self, b: TFBundle) -> TFBundle:
        r"""
        Reconstruct ``Z`` from apparent resistivity and phase.
        
        If a bundle lacks ``Z`` but holds ``(rho, phase)`` and the
        global config enables ``compute_z_from_res``, this method
        builds a complex tensor whose magnitude and angle satisfy
        the standard MT relations:
        
        * ``|Z| = sqrt( μ0 * ω * ρa )``
        * ``phase`` is interpreted as milliradians and converted to
          radians before forming ``Z = |Z|(cos φ + i sin φ)``.
        
        Parameters
        ----------
        b : TFBundle
            Input bundle with ``rho`` and ``phase`` set.
        
        Returns
        -------
        TFBundle
            The same bundle with ``z`` populated.
        
        Notes
        -----
        Ensure phase units are consistent. Zonge workflows often
        store milliradians; the implementation converts by dividing
        by ``1000``. If units differ, override this method.
        
        References
        ----------
        .. [1] Simpson, F. & Bahr, K. (2005). *Practical MT*.
        .. [2] Egbert, G. D. (1997). Robust MT processing.
        """

        if b.freq is None or b.rho is None or b.phase is None:
            return b
        f = np.asarray(b.freq, float)
        wmu = 2.0 * float(np.pi) * self.MU0 * f
        amp = np.sqrt(np.asarray(b.rho) * wmu)
        # assume milliradians → radians
        phi = np.asarray(b.phase, float) / 1000.0
        c = np.cos(phi)
        s = np.sin(phi)
        z = amp[..., None, None] * (c + 1j * s)
        b.z = z  # type: ignore[assignment]
        return b

    def _iter_bundles(self, avg: Any) -> list[TFBundle]:
        r"""
        Extract per-station :class:`TFBundle` objects from an AVG.
        
        Returns
        -------
        list of TFBundle
            One bundle per station. May be empty if no TFs exist.
        """

        z_t, f, st = None, None, None
        try:
            z_t, f, st = avg.to_tensor(
                var="z",
                station=None,
                sort_freq=True,
                align="union",
            )
        except Exception:
            pass
        z_e = None
        try:
            z_e, _, _ = avg.to_tensor(var="z_err")
        except Exception:
            pass
        r_t = p_t = None
        if z_t is None:
            try:
                r_t, f, st = avg.to_tensor(var="rho")
                p_t, _, _ = avg.to_tensor(var="phase")
            except Exception:
                pass
        out: list[TFBundle] = []
        if z_t is None and r_t is None:
            return out

        def _as4(a: Any) -> Any:
            if a is None:
                return None
            if a.ndim == 3:
                return a[None, ...]
            return a

        z_t = _as4(z_t)
        z_e = _as4(z_e)
        r_t = _as4(r_t)
        p_t = _as4(p_t)
        n_site = int(len(st)) if st is not None else 1
        for i in range(n_site):
            z = None if z_t is None else z_t[i]
            ze = None if z_e is None else z_e[i]
            r = None if r_t is None else r_t[i]
            p = None if p_t is None else p_t[i]
            sid = None if st is None else st[i]
            name = None
            b = TFBundle(
                station=name,
                station_id=sid,
                freq=f,
                z=z,
                z_err=ze,
                rho=r,
                phase=p,
                tipper=None,
                tipper_err=None,
            )
            out.append(b)
        return out

    def extract(self, source: Any) -> TFBundle:
        r"""
        Pull the first :class:`TFBundle` from an AVG object or path.
        
        Parameters
        ----------
        source : Any
            :class:`AVG` instance or path to an ``.avg`` file.
        
        Returns
        -------
        TFBundle
            The first available bundle (typically first station).
        
        Raises
        ------
        ValueError
            If no transfer functions could be obtained.
        
        Notes
        -----
        Multi-station inputs are handled by :meth:`transform`,
        which iterates and converts all sites to an EDI collection.
        """
        avg = self._as_avg(source)
        bs = self._iter_bundles(avg)
        if not bs:
            raise ValueError("no transfer functions in AVG")
        return bs[0]

    def emit_edi(self, bundle: TFBundle) -> Any:
        r"""
        Materialize an :class:`EDIFile` from a finalized bundle.
        
        Populates ``Z`` fields (``_freq``, ``_z``, ``_z_err``). If
        ``Z`` is missing but ``(rho, phase)`` are present, uses the
        backend method to set resistivity and phase. If tipper data
        exist, populates tipper arrays and derived attributes.
        
        Parameters
        ----------
        bundle : TFBundle
            Finalized transfer-function payload.
        
        Returns
        -------
        EDIFile
            A concrete EDI object.
        
        Notes
        -----
        Errors during optional computations are silenced to keep
        the transformation robust on imperfect field data.
        """

        ed = EDIFile(verbose=0)
        stub = self._HeadStub(
            bundle.station or "",
            empty=self._empty(),
        )
        ed.add_section("head", stub)
        if bundle.freq is not None:
            ed.Z._freq = np.asarray(bundle.freq, float)
        if bundle.z is not None:
            ed.Z._z = np.asarray(bundle.z, complex)
        if bundle.z_err is not None:
            ed.Z._z_err = np.asarray(bundle.z_err, float)
        if ed.Z._z is None and bundle.rho is not None:
            try:
                ed.Z.set_res_phase(
                    np.asarray(bundle.rho, float),
                    np.asarray(bundle.phase, float),
                    ed.Z._freq,
                )
            except Exception:
                pass
        try:
            ed.Z.compute_resistivity_phase()
        except Exception:
            pass
        if bundle.tipper is not None:
            ed.Tip._freq = ed.Z._freq
            ed.Tip._tipper = np.asarray(bundle.tipper, complex)
            if bundle.tipper_err is not None:
                ed.Tip._tipper_err = np.asarray(
                    bundle.tipper_err, float
                )
            try:
                ed.Tip.compute_amp_phase()
                ed.Tip.compute_mag_direction()
            except Exception:
                pass
        return ed

    def post_emit(
        self,
        edi_obj: Any,
        source: Any,
        bundle: TFBundle,
    ) -> Any:
        r"""
        Attach station metadata and optional location to the EDI.
        
        This step applies the global station naming policy and, if
        ``avg.topo.frame`` is present, injects coordinates into the
        ``>HEAD`` section.
        
        Parameters
        ----------
        edi_obj : EDIFile
            Newly created EDI object.
        source : AVG
            Original AVG used to derive the bundle.
        bundle : TFBundle
            Finalized bundle for this station.
        
        Returns
        -------
        EDIFile
            The updated EDI object.
        
        Notes
        -----
        When a matching row is found in ``topo.frame`` (based on the
        ``station`` column), latitude, longitude and elevation are
        copied if available.
        """

        try:
            nm = ensure_station(
                bundle.station,
                bundle.station_id,
            )
            edi_obj.station = nm
        except Exception:
            pass
        try:
            topo = getattr(source, "topo", None)
            if topo is None or getattr(topo, "frame", None) is None:
                return edi_obj
            tb = topo.frame
            col = "station"
            if col in tb.columns:
                sid = bundle.station_id
                row = tb[tb[col] == sid]
                if not row.empty:
                    lat = float(row["latitude"].iloc[0])
                    lon = float(row["longitude"].iloc[0])
                    if "elevation" in row.columns:
                        elv = float(row["elevation"].iloc[0])
                    else:
                        elv = float("nan")
                    h = self._HeadStub(
                        edi_obj.station or "",
                        lat=lat,
                        lon=lon,
                        elev=elv,
                        empty=self._empty(),
                    )
                    edi_obj.add_section("head", h)
        except Exception:
            pass
        return edi_obj

    def transform(
        self,
        source: Any,
        *,
        name: Optional[str] = None,
        station_id: Optional[str | int] = None,
    ) -> Any:
        r"""
        Convert an AVG source to an :class:`EDICollection`.
        
        Handles both objects and file paths. All stations found in
        the source are extracted, finalized, converted to EDI, and
        collected.
        
        Parameters
        ----------
        source : Any
            :class:`AVG` instance or path-like.
        name : str, optional
            Preferred station name override, applied during
            finalization when present.
        station_id : str or int, optional
            Identifier used by the naming policy if a name must be
            synthesized.
        
        Returns
        -------
        EDICollection
            A collection of :class:`EDIFile` objects, one per site.
        
        Examples
        --------
        >>> from pycsamt.core.transformers import AVGtoEDI
        >>> # out = AVGtoEDI().transform(\"/path/site.avg\")  # doctest:+SKIP
        """

        avg = self._as_avg(source)
        bundles = self._iter_bundles(avg)
        edis: list[Any] = []
        for b in bundles:
            b = self._finalize(
                b,
                name=name,
                station_id=station_id,
            )
            ed = self.emit_edi(b)
            ed = self.post_emit(ed, avg, b)
            edis.append(ed)
        
        return EDICollection(items=edis, verbose=0)

class JtoEDI(TransformerMixin):
    r"""
    Convert a Jones ``J`` source to EDI or an EDI collection.
    
    This transformer ingests a :class:`~pycsamt.jones.j.JFile`
    instance or a path to a ``.j`` file, extracts transfer
    functions (``Z``, tipper, or fallback ``rho/phase``),
    finalizes a neutral payload
    (:class:`~pycsamt.core.base.TFBundle`), and emits concrete
    :class:`~pycsamt.seg.edi.EDIFile` objects. When given a
    :class:`~pycsamt.jones.collection.JCollection`, it produces
    an :class:`~pycsamt.seg.collection.EDICollection`.
    
    The pipeline mirrors the AVG workflow:
    
    1. Extract TF arrays and metadata from the J structure,
       tolerating variant attribute names often found in legacy
       files.
    2. Finalize each bundle using the mixin logic: ensure a
       valid station name, order and de-duplicate frequencies,
       and fill missing TF parts if configured.
    3. Emit EDI and optionally attach site coordinates if a
       header-like structure is present on the J side.
    
    Notes
    -----
    * Naming, frequency sorting, and de-duplication follow the
      global configuration in :mod:`pycsamt.core.config`.
    * If only ``(rho, phase)`` are present, ``Z`` can be
      reconstructed when allowed by
      ``compute_z_from_res`` (via the mixin hook).
    * A ``Head`` or ``head`` object on ``JFile`` with attributes
      ``lat``, ``long``, and ``elev`` is used to populate the EDI
      ``>HEAD`` section when available.
    
    See Also
    --------
    pycsamt.core._transformers.TransformerMixin
        Finalization steps and overridable hooks.
    pycsamt.core.transformers.AVGtoEDI
        Companion converter for Zonge AVG.
    pycsamt.core.base.TFBundle
        Neutral payload used across backends.
    
    Examples
    --------
    >>> from pycsamt.core.transformers import JtoEDI
    >>> # Single file → EDIFile
    >>> # edi = JtoEDI().transform(\"/path/site.j\")  # doctest:+SKIP
    >>> # Collection → EDICollection
    >>> # out = JtoEDI().transform(j_collection)     # doctest:+SKIP
    
    References
    ----------
    .. [1] Simpson, F. & Bahr, K. (2005). *Practical MT*.
    .. [2] Egbert, G. D. (1997). Robust MT processing.
    """

    class _HeadStub:
        r"""
        Private header stub for seeding the EDI ``>HEAD`` section.
        
        Parameters
        ----------
        dataid : str
            Identifier for the dataset, usually station name.
        lat, lon, elev : float or None
            Optional geographic location and elevation.
        empty : float or None
            Empty sentinel copied into the header.
        """

        def __init__(
            self,
            dataid: str,
            *,
            lat: float | None = None,
            lon: float | None = None,
            elev: float | None = None,
            empty: float | None = None,
        ) -> None:
            self.dataid = dataid
            self.lat = lat
            self.long = lon
            self.elev = elev
            if empty is not None:
                self.empty = float(empty)

    def _as_jfile(self, src: Any) -> JFile:
        
        if isinstance(src, JFile):
            return src
        if isinstance(src, (str, Path)):
            return JFile.from_file(src)
        raise TypeError("source must be JFile or path")

    def _empty(self) -> float:
        return float(get_config().empty)

    def _bundle_from_j(self, jf: JFile) -> TFBundle:
        r"""
        Build a :class:`TFBundle` by probing known J attributes.

        Handles common aliases, e.g. ``Z``/``z``, ``ResPhase``/``RP``,
        and tipper fields. Missing parts are left as ``None`` and can
        be filled during finalization if configured.
        """

        def g(o: Any, *ns: str) -> Any:
            for n in ns:
                if hasattr(o, n):
                    return getattr(o, n)
            return None

        z = g(jf, "Z", "z")
        tip = g(jf, "Tipper", "Tip", "tip")
        rp = g(jf, "ResPhase", "resphase", "RP")

        freq = None
        z_arr = z.z if z and hasattr(z, "z") else g(z, "Z")
        z_err = g(z, "z_err", "Z_err", "z_error")
        if z and hasattr(z, "freq"):
            freq = z.freq
        elif hasattr(jf, "freq"):
            freq = jf.freq

        rho = phase = None
        if rp is not None:
            rho = g(rp, "rho", "resistivity")
            phase = g(rp, "phase", "phi")
        else:
            rho = g(jf, "rho", "resistivity")
            phase = g(jf, "phase", "phi")

        tip_arr = None
        tip_err = None
        if tip is not None:
            tip_arr = g(tip, "tipper", "T", "t")
            tip_err = g(tip, "tipper_err", "tip_err")

        name = g(jf, "station", "site", "name")
        lat = g(jf, "lat", "latitude")
        lon = g(jf, "lon", "longitude")
        elev = g(jf, "elev", "elevation")

        return TFBundle(
            freq=freq,
            z=z_arr,
            z_err=z_err,
            tipper=tip_arr,
            tipper_err=tip_err,
            rho=rho,
            phase=phase,
            station=name,
            lat=lat,
            lon=lon,
            elev=elev,
        )

    def extract(self, source: Any) -> TFBundle:
        r"""
        Extract a :class:`TFBundle` from a J file or object.
        
        Parameters
        ----------
        source : Any
            :class:`JFile` instance or a path to a ``.j`` file.
        
        Returns
        -------
        TFBundle
            Bundle with TF arrays and basic site info.
        
        Raises
        ------
        ValueError
            If no usable transfer functions are present.
        
        Notes
        -----
        This method does not reorder or de-duplicate; those steps
        are applied in the finalization phase.
        """

        jf = self._as_jfile(source)
        b = self._bundle_from_j(jf)
        if b.is_empty():
            raise ValueError("no transfer functions in J file")
        return b

    def emit_edi(self, bundle: TFBundle) -> Any:
        r"""
        Create an :class:`EDIFile` from a finalized bundle.
        
        Populates impedance arrays when present. If ``Z`` is absent
        but ``(rho, phase)`` are available, the backend call is used
        to set resistivity/phase onto the EDI structure. Tipper data
        are also populated when provided.
        
        Parameters
        ----------
        bundle : TFBundle
            Finalized transfer-function payload.
        
        Returns
        -------
        EDIFile
            Concrete EDI object ready for downstream use.
        
        Notes
        -----
        Optional computations are guarded; failures are ignored to
        remain robust on imperfect field data.
        """

        ed = EDIFile(verbose=0)
        stub = self._HeadStub(
            bundle.station or "",
            empty=self._empty(),
        )
        ed.add_section("head", stub)
        if bundle.freq is not None:
            ed.Z._freq = np.asarray(bundle.freq, float)
        if bundle.z is not None:
            ed.Z._z = np.asarray(bundle.z, complex)
        if bundle.z_err is not None:
            ed.Z._z_err = np.asarray(bundle.z_err, float)
        if ed.Z._z is None and bundle.rho is not None:
            try:
                ed.Z.set_res_phase(
                    np.asarray(bundle.rho, float),
                    np.asarray(bundle.phase, float),
                    ed.Z._freq,
                )
            except Exception:
                pass
        try:
            ed.Z.compute_resistivity_phase()
        except Exception:
            pass
        if bundle.tipper is not None:
            ed.Tip._freq = ed.Z._freq
            ed.Tip._tipper = np.asarray(bundle.tipper, complex)
            if bundle.tipper_err is not None:
                ed.Tip._tipper_err = np.asarray(
                    bundle.tipper_err, float
                )
            try:
                ed.Tip.compute_amp_phase()
                ed.Tip.compute_mag_direction()
            except Exception:
                pass
        return ed

    def post_emit(
        self,
        edi_obj: Any,
        source: Any,
        bundle: TFBundle,
    ) -> Any:
        r"""
        Attach naming and optional coordinates to the EDI.
        
        Applies the global station naming policy and, when the J
        object exposes a ``Head``/``head`` with ``lat``, ``long``,
        and ``elev``, updates the EDI header accordingly.
        
        Parameters
        ----------
        edi_obj : EDIFile
            Newly created EDI object.
        source : JFile
            Original J source used for extraction.
        bundle : TFBundle
            Finalized bundle for this site.
        
        Returns
        -------
        EDIFile
            The same EDI object, augmented in place.
        """

        try:
            nm = ensure_station(
                bundle.station,
                bundle.station_id,
            )
            edi_obj.station = nm
        except Exception:
            pass
        try:
            head = getattr(source, "Head", None)
            if head is None:
                head = getattr(source, "head", None)
            if head is not None:
                lat = getattr(head, "lat", None)
                lon = getattr(head, "long", None)
                elev = getattr(head, "elev", None)
                h = self._HeadStub(
                    edi_obj.station or "",
                    lat=lat,
                    lon=lon,
                    elev=elev,
                    empty=self._empty(),
                )
                edi_obj.add_section("head", h)
        except Exception:
            pass
        return edi_obj

    def transform(
        self,
        source: Any,
        *,
        name: Optional[str] = None,
        station_id: Optional[str | int] = None,
    ) -> Any:
        r"""
        Convert a J source to EDI or an EDI collection.
        
        A :class:`JCollection` yields an
        :class:`EDICollection` with one EDI per member. A single
        ``JFile`` (or path) yields a single :class:`EDIFile`.
        
        Parameters
        ----------
        source : Any
            :class:`JFile` instance, path-like, or
            :class:`JCollection`.
        name : str, optional
            Preferred station name override for finalization.
        station_id : str or int, optional
            Identifier used when a name must be synthesized.
        
        Returns
        -------
        EDIFile or EDICollection
            EDI materialization corresponding to the input form.
        
        Examples
        --------
        >>> from pycsamt.core.transformers import JtoEDI
        >>> # Single file
        >>> # edi = JtoEDI().transform(\"/data/site.j\")    # doctest:+SKIP
        >>> # Collection
        >>> # out = JtoEDI().transform(j_collection)       # doctest:+SKIP
        """

        if isinstance(source, JCollection):
            edis: list[Any] = []
            for jf in source:
                b = self._bundle_from_j(jf)
                b = self._finalize(
                    b,
                    name=name,
                    station_id=station_id,
                )
                ed = self.emit_edi(b)
                ed = self.post_emit(ed, jf, b)
                edis.append(ed)

            return EDICollection(items=edis, verbose=0)
        
        jf = self._as_jfile(source)
        b = self._bundle_from_j(jf)
        b = self._finalize(
            b,
            name=name,
            station_id=station_id,
        )
        ed = self.emit_edi(b)
        
        return self.post_emit(ed, jf, b)
