# -*- coding: utf-8 -*-
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
    """AVG → EDICollection transformer."""

    class _HeadStub:
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

    def _mu0(self) -> float:
        return 4.0e-7 * float(np.pi)

    # optional math hooks 
    def compute_z_from_res(self, b: TFBundle) -> TFBundle:
        if b.freq is None or b.rho is None or b.phase is None:
            return b
        f = np.asarray(b.freq, float)
        wmu = 2.0 * float(np.pi) * self._mu0() * f
        amp = np.sqrt(np.asarray(b.rho) * wmu)
        # assume milliradians → radians
        phi = np.asarray(b.phase, float) / 1000.0
        c = np.cos(phi)
        s = np.sin(phi)
        z = amp[..., None, None] * (c + 1j * s)
        b.z = z  # type: ignore[assignment]
        return b

    # core extraction 
    def _iter_bundles(self, avg: Any) -> list[TFBundle]:
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
        avg = self._as_avg(source)
        bs = self._iter_bundles(avg)
        if not bs:
            raise ValueError("no transfer functions in AVG")
        return bs[0]

    def emit_edi(self, bundle: TFBundle) -> Any:
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
    """J → EDI/EDICollection transformer."""

    class _HeadStub:
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
        jf = self._as_jfile(source)
        b = self._bundle_from_j(jf)
        if b.is_empty():
            raise ValueError("no transfer functions in J file")
        return b

    def emit_edi(self, bundle: TFBundle) -> Any:
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
