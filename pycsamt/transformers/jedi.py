# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from ..core.base import TFBundle, ensure_station
from ..core.config import get_config
from ..jones.collection import JCollection
from ..jones.j import JFile
from ..seg.collection import EDICollection
from ..seg.edi import EDIFile
from ..seg.heads import Head, Info
from ..seg.meas import (
    DefineMeas,
    Emeasurement,
    Hmeasurement,
)
from ..seg.mtemap import MTEMAP
from ..zonge.avg import AVG, BaseAVG
from ._base import TransformerMixin

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
            **kws,
        ) -> None:
            self.dataid = dataid
            self.lat = lat
            self.long = lon
            self.elev = elev
            if empty is not None:
                self.empty = float(empty)
            for key in list(kws.keys()):
                setattr(self, key, kws[key])

    def _as_avg(self, src: Any) -> Any:

        # BaseAVG is the documented extension point; AVG is kept
        # in the tuple because tests may rebind it to a stub that
        # does not inherit BaseAVG.
        if isinstance(src, (AVG, BaseAVG)):
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

        attrs = {
            "has_any_magnetic": self._has_any_magnetic(avg),
            "comp": self._unique_components_from_avg(avg),
        }

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
                attrs=attrs,
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

    def _unique_components_from_avg(self, avg) -> set[str]:
        """
        Return the distinct component labels present in
        the AVG measurement table.
        Works with CompMeas exposed via `avg.info.comp.frame`
        and its 'comp' column.
        """
        try:
            comp_col = avg.info.comp.frame.get("comp")
            if comp_col is None:
                return set()
            # Coerce to plain strings (handles categorical dtype).
            return set(map(str, comp_col.unique()))
        except Exception:
            return set()

    def _has_any_magnetic(self, avg) -> bool:
        """
        True if the AVG has any magnetic channel (Hx, Hy, Hz)
        or Bx/By/Bz implied
        in the component names (e.g., 'ExHy', 'EyHx', 'Zxy', 'Zyx', ...).
        """
        comps = self._unique_components_from_avg(avg)
        if not comps:
            return False
        # Look for 'H' or 'B' tokens anywhere in the component label.
        # This covers 'ExHy', 'EyHx', 'Bz', etc. Case-insensitive.
        return any(("H" in c.upper()) or ("B" in c.upper()) for c in comps)

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

        # ---------- seed HEAD + INFO (coordinates updated later)
        _head = Head()
        _head.dataid = (bundle.station or "").strip()
        _head.stdvers = "SEG 1.0"
        _head.progvers = "PYCSAMT"
        _head.empty = float(self._empty())
        # default zeros when topography is unknown (updated later)
        _head.lat = 0.0
        _head.long = 0.0
        _head.elev = 0.0
        ed.add_section("head", _head)

        _info = Info()
        # seed INFO (no .add; use update/info_text)
        _info.update(
            info_text=[
                f"SURVEY ID:{_head.dataid or ''}",
                "ROTATION=FIX",
            ]
        )
        ed.add_section("info", _info)

        # ---------- >=DEFINEMEAS: build ids for channels present
        attrs = bundle.attrs or {}
        comps = set(map(str, attrs.get("comp", [])))
        has_mag = bool(attrs.get("has_any_magnetic"))

        def _has(tok: str) -> bool:
            t = tok.upper()
            return any(t in c.upper() for c in comps)

        want_ex = True if not comps else _has("EX")
        want_ey = True if not comps else _has("EY")
        want_hx = has_mag and (_has("HX") or _has("BX"))
        want_hy = has_mag and (_has("HY") or _has("BY"))
        want_hz = has_mag and (_has("HZ") or _has("BZ"))

        dm = DefineMeas()
        dm.units = "M"
        dm.reftype = "CART"

        # seed origin from HEAD (updated later in post_emit)
        dm.reflat = _head.lat
        dm.reflong = _head.long
        dm.refelev = _head.elev

        _ids: dict[str, str] = {}

        def _hid(i: int) -> str:
            return f"{100 + i}.001"

        k = 0
        if want_hx:
            _ids["HX"] = _hid(k)
            k += 1
            dm.hmeas.append(
                Hmeasurement(
                    id=_ids["HX"],
                    chtype="HX",
                    x=0,
                    y=0,
                    z=0,
                    azm=0,
                )
            )
        if want_hy:
            _ids["HY"] = _hid(k)
            k += 1
            dm.hmeas.append(
                Hmeasurement(
                    id=_ids["HY"],
                    chtype="HY",
                    x=0,
                    y=0,
                    z=0,
                    azm=90,
                )
            )
        if want_hz:
            _ids["HZ"] = _hid(k)
            k += 1
            dm.hmeas.append(
                Hmeasurement(
                    id=_ids["HZ"],
                    chtype="HZ",
                    x=0,
                    y=0,
                    z=0,
                    azm=0,
                )
            )
        if want_ex:
            _ids["EX"] = _hid(k)
            k += 1
            dm.emeas.append(
                Emeasurement(
                    id=_ids["EX"],
                    chtype="EX",
                    x=0,
                    y=0,
                    z=0,
                    x2=0,
                    y2=0,
                    z2=0,
                )
            )
        if want_ey:
            _ids["EY"] = _hid(k)
            k += 1
            dm.emeas.append(
                Emeasurement(
                    id=_ids["EY"],
                    chtype="EY",
                    x=0,
                    y=0,
                    z=0,
                    x2=0,
                    y2=0,
                    z2=0,
                )
            )

        ed.add_section("definemeas", dm)

        # ---------- >=MTSECT or >=EMAPSECT
        mt = MTEMAP()
        mt.sectid = _head.dataid or ""
        freq = bundle.freq
        mt.nfreq = int(len(freq) if freq is not None else [])

        if has_mag:
            mt.hx = _ids.get("HX")
            mt.hy = _ids.get("HY")
            mt.hz = _ids.get("HZ")
            mt.ex = _ids.get("EX")
            mt.ey = _ids.get("EY")
            ed.add_section("mtsect", mt)
            ed.emap = False
        else:
            mt.ex = _ids.get("EX")
            mt.ey = _ids.get("EY")
            ed.add_section("mtsect", mt)
            ed.emap = True

        if bundle.freq is not None:
            ed.Z._freq = np.asarray(bundle.freq, dtype=float)
        if bundle.z is not None:
            ed.Z._z = np.asarray(bundle.z, dtype=complex)
        if bundle.z_err is not None:
            ed.Z._z_err = np.asarray(bundle.z_err, dtype=float)

        if bundle.phase is not None:
            ed.Z._phase = np.asarray(bundle.phase, dtype=float)

        if bundle.z is not None:
            # --- after you set ed.Z._freq and ed.Z._z ---
            z = np.asarray(ed.Z._z, dtype=complex)
            f = np.asarray(ed.Z._freq, dtype=float)

            # 1) keep only rows that have at least one finite component
            keep = np.isfinite(z.real).any(axis=(1, 2)) | np.isfinite(
                z.imag
            ).any(axis=(1, 2))
            if keep.ndim:  # guard for scalar shapes
                z = z[keep]
                f = f[keep]

            # 2) replace remaining NaN/Inf with zeros
            z = np.nan_to_num(z, nan=0.0, posinf=0.0, neginf=0.0)

            # and keep error array consistent if present
            ze = ed.Z._z_err
            if ze is not None:
                ze = np.asarray(ze, dtype=float)
                if keep.ndim:
                    ze = ze[keep]
                ze = np.nan_to_num(ze, nan=0.0, posinf=0.0, neginf=0.0)

            ed.Z._z = z
            ed.Z._freq = f
            ed.Z._z_err = ze

            # now safe:
            ed.Z.compute_resistivity_phase()

        elif (bundle.rho is not None) and (bundle.phase is not None):
            ed.Z._resistivity = np.asarray(bundle.rho, dtype=float)
            freq = np.asarray(ed.Z._freq, dtype=float)
            rho = np.asarray(ed.Z._resistivity, dtype=float)
            # mrad → degrees
            # mrad -> deg
            # phase in degrees (AVG metadata shows Unit.Phase: 'mrad').
            # Convert mrad → deg before calling the backend setter:
            phi_deg = ed.Z._phase * (180.0 / (1000.0 * np.pi))

            ed.Z.set_res_phase(rho, phi_deg, freq)
        else:
            # no Z and no (ρ,φ): Never rich here.
            raise ValueError("Neither Z nor (rho, phase) provided.")

        if bundle.tipper is not None:
            ed.Tip._freq = ed.Z._freq
            ed.Tip._tipper = np.asarray(bundle.tipper, dtype=complex)
            if bundle.tipper_err is not None:
                ed.Tip._tipper_err = np.asarray(
                    bundle.tipper_err, dtype=float
                )
            try:
                ed.Tip.compute_amp_phase()
                ed.Tip.compute_mag_direction()
            except Exception:
                pass
        return ed

    def _apply_topo_coords(
        self,
        ed,
        topo_frame,
        station_id,
    ):
        if topo_frame is None or getattr(topo_frame, "empty", True):
            return

        def _pick(names):
            for n in names:
                if n in topo_frame.columns:
                    return n
            return None

        sid_col = _pick(("station", "site", "id", "SITE", "STATION"))
        if sid_col is None:
            return

        row = topo_frame[topo_frame[sid_col] == station_id]

        if row.empty:
            row = topo_frame[topo_frame[sid_col] == str(station_id)]
            if row.empty:
                return

        def _num(v, d=0.0):
            try:
                return float(v)
            except Exception:
                return d

        lat_col = _pick(("latitude", "lat", "LAT"))
        lon_col = _pick(("longitude", "lon", "long", "LON", "LONG"))
        elv_col = _pick(("elevation", "elev", "alt", "ALT"))

        has_latlon = lat_col is not None and lon_col is not None
        lat = _num(row[lat_col].iloc[0]) if has_latlon else None
        lon = _num(row[lon_col].iloc[0]) if has_latlon else None
        elv = _num(row[elv_col].iloc[0]) if elv_col else 0.0

        h = ed.get_section("head")
        if h is None:
            h = self._ensure_head(
                ed,
                station=str(station_id),
                empty=self._empty(),
            )

        if has_latlon:
            h.lat = lat
            h.long = lon
        h.elev = elv

        dm = ed.get_section("definemeas")
        if dm is not None:
            if has_latlon:
                dm.reflat = lat
                dm.reflong = lon
            dm.refelev = elv

        # stub = self._HeadStub(
        #     ed.station or "",
        #     lat=lat,
        #     lon=lon,
        #     elev=elv,
        #     empty=self._empty(),
        # )
        # h = self._ensure_head(
        #     ed,
        #     station=stub.dataid,
        #     empty=stub.empty,
        # )
        # # copy stub attrs onto Head, but don't add_section again
        # for k in ("lat", "long", "elev"):
        #     v = getattr(stub, k, None)
        #     if v is not None:
        #         setattr(h, k, v)
        # # copy any extra kws carried by stub
        # for k, v in stub.__dict__.items():
        #     if k not in {"dataid", "lat", "long", "elev", "empty"}:
        #         setattr(h, k, v)

    def _enrich_from_avg_info(self, ed, source):
        # Accept flexible sources:
        # - mapping-like: source.info (dict-like)
        # - obj with attrs: source.info.<field>
        meta = getattr(source, "info", None)
        if meta is None:
            return

        def _get(k: str, default: str = ""):
            try:
                if hasattr(meta, "get"):
                    return meta.get(k, default)
                return getattr(meta, k, default)
            except Exception:
                return default

        def _infoln(s: str) -> str:
            return s if s.startswith("  ") else "  " + s

        # HEAD enrich
        h = ed.get_section("head")
        if h is None:
            h = self._ensure_head(
                ed,
                station=getattr(ed, "station", "") or "",
                empty=self._empty(),
            )

        val = _get("stdvers", None)
        if val:
            h.stdvers = str(val)
        val = _get("progvers", None)
        if val:
            h.progvers = str(val)
        val = _get("progdate", None)
        if val:
            h.progdate = str(val)
        val = _get("acqdate", None)
        if val:
            h.acqdate = str(val)
        val = _get("filedate", None)
        if val:
            h.filedate = str(val)
        val = _get("acqby", None)
        if val:
            h.acqby = str(val)
        val = _get("fileby", None)
        if val:
            h.fileby = str(val)
        val = _get("prospect", None)
        if val:
            h.prospect = str(val)
        val = _get("loc", None)
        if val:
            h.loc = str(val)
        val = _get("maxsect", None)
        if val is not None:
            try:
                h.maxsect = int(val)
            except Exception:
                pass
        val = _get("empty", None)
        if val is not None:
            try:
                h.empty = float(val)
            except Exception:
                pass

        # INFO enrich
        info = ed.get_section("info")
        if info is None:
            info = self._ensure_info(
                ed,
                survey_id=getattr(ed, "station", "") or "",
            )

        # extend free-text lines (preserved verbatim)
        txt = list(getattr(info, "info_text", []))
        add = []
        v = _get("survey_co", None)
        if v:
            add.append(f"SURVEY CO:{v}")
        v = _get("client_co", None)
        if v:
            add.append(f"CLIENT CO:{v}")
        v = _get("area", None)
        if v:
            add.append(f"AREA:{v}")
        # keep ROTATION and SURVEY ID if already present
        have_sid = any(
            s.strip().upper().startswith("SURVEY ID:") for s in txt
        )
        if not have_sid:
            sid = getattr(ed, "station", "") or ""
            add.insert(0, _infoln(f"SURVEY ID:{sid}"))
        if not any("ROTATION=" in s or "ROTATION:" in s for s in txt):
            add.append(_infoln("ROTATION=FIX"))

        if add:
            info.update(info_text=txt + add)

    def post_emit(self, edi_obj, source, bundle):
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
            # nm = (bundle.station or "").strip()
            nm = ensure_station(
                bundle.station,
                bundle.station_id,
            )
            edi_obj.station = nm

            # if nm:
            #     edi_obj.station = nm
        except Exception:
            pass

        # ensure HEAD/INFO exist early
        try:
            self._ensure_head(
                edi_obj,
                station=getattr(edi_obj, "station", "") or "",
                empty=self._empty(),
            )
            self._ensure_info(
                edi_obj,
                survey_id=getattr(edi_obj, "station", "") or "",
            )
        except Exception:
            pass

        # coords from topo if available
        try:
            topo = getattr(source, "topo", None)
            fr = getattr(topo, "frame", None)
            self._apply_topo_coords(
                edi_obj,
                fr,
                bundle.station_id,
            )
        except Exception:
            pass

        # enrich from AVG.info (if any)
        try:
            self._enrich_from_avg_info(edi_obj, source)
        except Exception:
            pass

        return edi_obj

    def transform(
        self,
        source: Any,
        *,
        name: str | None = None,
        station_id: str | int | None = None,
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
        # real HEAD/INFO sections up front
        self._ensure_head(
            ed,
            station=(bundle.station or ""),
            empty=self._empty(),
        )
        self._ensure_info(
            ed,
            survey_id=(bundle.station or ""),
        )
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
                ed.Tip._tipper_err = np.asarray(bundle.tipper_err, float)
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

        # 2) ensure HEAD + INFO exist, then fill defaults
        try:
            self._ensure_head(
                edi_obj,
                station=getattr(edi_obj, "station", "") or "",
                empty=self._empty(),
            )
            self._ensure_info(
                edi_obj,
                survey_id=getattr(edi_obj, "station", "") or "",
            )
        except Exception:
            pass

        # 3) copy coordinates from J headers (if present)
        try:
            self._apply_j_coords(edi_obj, source)
        except Exception:
            pass

        # 4) seed DefineMeas + MTSECT from actual TF content
        try:
            self._seed_meas_and_mtsect(edi_obj, source)
        except Exception:
            pass

        # 5) optional enrichment from J info (lightweight)
        try:
            info = getattr(source, "heads", None)
            if info is not None and getattr(info, "info", None):
                _info = edi_obj.get_section("info")
                # keep it small and predictable; don’t overwrite
                _info.setdefault(
                    "PROCESSEDBY",
                    f"pyCSAMT v{get_config().version}",
                )
                _info.setdefault("SIGNCONVENTION", "EXP(+I Ω T)")
        except Exception:
            pass

        return edi_obj

    def transform(
        self,
        source: Any,
        *,
        name: str | None = None,
        station_id: str | int | None = None,
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

    def _apply_j_coords(self, ed, jf) -> None:
        lat = getattr(jf, "lat", None)
        lon = getattr(jf, "lon", None)
        elv = getattr(jf, "elev", None)
        if lat is None and lon is None and elv is None:
            return
        h = self._ensure_head(
            ed,
            station=getattr(ed, "station", "") or "",
            empty=self._empty(),
        )
        if lat is not None:
            h.lat = float(lat)
        if lon is not None:
            h.long = float(lon)
        if elv is not None:
            h.elev = float(elv)

        dm = ed.get_section("definemeas")
        if dm is not None:
            dm.reflat = getattr(h, "lat", 0.0)
            dm.reflong = getattr(h, "long", 0.0)
            dm.refelev = getattr(h, "elev", 0.0)

    def _seed_meas_and_mtsect(self, ed, jf) -> None:

        z = getattr(ed.Z, "z", None)
        f = getattr(ed.Z, "freq", None)
        nfreq = int(f.size) if f is not None else 0

        tip = getattr(ed, "Tip", None)
        has_tip = bool(
            tip is not None and getattr(tip, "tipper", None) is not None
        )

        have_xy = have_yx = False
        if z is not None and np.size(z) > 0:
            z = np.asarray(z)
            if z.ndim == 3 and z.shape[1:] == (2, 2):
                have_xy = np.any(np.abs(z[:, 0, 1]) > 0)
                have_yx = np.any(np.abs(z[:, 1, 0]) > 0)

        ids = {
            "HX": "100.001",
            "HY": "100.002",
            "HZ": "100.003",
            "EX": "101.001",
            "EY": "101.002",
        }

        dm = self._ensure_definemeas(ed, units="M", reftype="CART")

        # ensure origin mirrors HEAD
        h = ed.get_section("head")
        if h is not None:
            dm.reflat = getattr(h, "lat", 0.0)
            dm.reflong = getattr(h, "long", 0.0)
            dm.refelev = getattr(h, "elev", 0.0)

        az = getattr(jf, "azimuth", None)
        azx = float(az) if az is not None else 0.0
        azy = (azx + 90.0) % 360.0

        def _has_h(id_):
            return any(getattr(m, "id", None) == id_ for m in dm.hmeas)

        def _has_e(id_):
            return any(getattr(m, "id", None) == id_ for m in dm.emeas)

        # H channels referenced later
        if have_yx and not _has_h(ids["HX"]):
            dm.hmeas.append(
                Hmeasurement(
                    id=ids["HX"],
                    chtype="HX",
                    x=0.0,
                    y=0.0,
                    z=0.0,
                    azm=azx,
                )
            )
        if have_xy and not _has_h(ids["HY"]):
            dm.hmeas.append(
                Hmeasurement(
                    id=ids["HY"],
                    chtype="HY",
                    x=0.0,
                    y=0.0,
                    z=0.0,
                    azm=azy,
                )
            )
        if has_tip and not _has_h(ids["HZ"]):
            dm.hmeas.append(
                Hmeasurement(
                    id=ids["HZ"],
                    chtype="HZ",
                    x=0.0,
                    y=0.0,
                    z=0.0,
                    azm=0.0,
                )
            )

        # E channels that pair with Z blocks
        if have_xy and not _has_e(ids["EX"]):
            dm.emeas.append(
                Emeasurement(
                    id=ids["EX"],
                    chtype="EX",
                    x=0.0,
                    y=0.0,
                    z=0.0,
                    x2=0.0,
                    y2=0.0,
                    z2=0.0,
                )
            )
        if have_yx and not _has_e(ids["EY"]):
            dm.emeas.append(
                Emeasurement(
                    id=ids["EY"],
                    chtype="EY",
                    x=0.0,
                    y=0.0,
                    z=0.0,
                    x2=0.0,
                    y2=0.0,
                    z2=0.0,
                )
            )

        mt = self._ensure_mtsect(
            ed,
            sectid=(getattr(ed, "station", "") or ""),
            nfreq=nfreq,
        )
        if have_xy:
            mt.hy = ids["HY"]
            mt.ex = ids["EX"]
        if have_yx:
            mt.hx = ids["HX"]
            mt.ey = ids["EY"]
        if has_tip:
            mt.hz = ids["HZ"]
