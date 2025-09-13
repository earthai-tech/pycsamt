# -*- coding: utf-8 -*-
from __future__ import annotations

from typing import Any, Optional
import numpy as np 

from .base import TFBundle, ensure_station
from .config import get_config

__all__ = [
    "TransformerMixin",
]


class TransformerMixin:
    """Lightweight template for AVG/J → EDI transformers.

    Subclasses implement ``extract`` and ``emit_edi``. This mixin
    finalizes a :class:`TFBundle`: validates station naming, orders and
    de‑duplicates frequencies, and (optionally) fills missing TF parts
    via small hooks that subclasses may override.
    """

    def extract(self, source: Any) -> TFBundle:  # noqa: D401
        raise NotImplementedError

    def emit_edi(self, bundle: TFBundle) -> Any:  # noqa: D401
        raise NotImplementedError

 
    def post_emit(
        self,
        edi_obj: Any,
        source: Any,
        bundle: TFBundle,
    ) -> Any:
        return edi_obj

    def compute_res_from_z(self, b: TFBundle) -> TFBundle:
        return b

    def compute_z_from_res(self, b: TFBundle) -> TFBundle:
        return b


    def _order_freq(self, b: TFBundle) -> TFBundle:
        if b.freq is None:
            return b
        order = get_config().freq_order
        if np is not None:
            arr = np.asarray(b.freq)
            idx = np.argsort(arr)
            if order == "desc":
                idx = idx[::-1]
            b.freq = arr[idx]
            if b.z is not None:
                b.z = b.z[idx]
            if b.z_err is not None:
                b.z_err = b.z_err[idx]
            if b.tipper is not None:
                b.tipper = b.tipper[idx]
            if b.tipper_err is not None:
                b.tipper_err = b.tipper_err[idx]
            if b.rho is not None:
                b.rho = b.rho[idx]
            if b.phase is not None:
                b.phase = b.phase[idx]
            return b
        seq = list(zip(b.freq, range(len(b.freq))))
        seq.sort(key=lambda x: x[0], reverse=(order == "desc"))
        take = [i for _, i in seq]
        b.freq = [b.freq[i] for i in take]
        if b.z is not None:
            b.z = [b.z[i] for i in take]
        if b.z_err is not None:
            b.z_err = [b.z_err[i] for i in take]
        if b.tipper is not None:
            b.tipper = [b.tipper[i] for i in take]
        if b.tipper_err is not None:
            b.tipper_err = [b.tipper_err[i] for i in take]
        if b.rho is not None:
            b.rho = [b.rho[i] for i in take]
        if b.phase is not None:
            b.phase = [b.phase[i] for i in take]
        return b

    def _dedup_freq(self, b: TFBundle) -> TFBundle:
        if b.freq is None:
            return b
        tol = get_config().freq_tol
        if np is not None:
            arr = np.asarray(b.freq, dtype=float)
            keep = [0]
            for i in range(1, arr.size):
                prev = arr[keep[-1]]
                if abs(arr[i] - prev) / max(prev, 1.0) > tol:
                    keep.append(i)
            idx = np.asarray(keep)
            b.freq = arr[idx]
            if b.z is not None:
                b.z = b.z[idx]
            if b.z_err is not None:
                b.z_err = b.z_err[idx]
            if b.tipper is not None:
                b.tipper = b.tipper[idx]
            if b.tipper_err is not None:
                b.tipper_err = b.tipper_err[idx]
            if b.rho is not None:
                b.rho = b.rho[idx]
            if b.phase is not None:
                b.phase = b.phase[idx]
            return b
        # fallback list mode
        out_f = []
        take = []
        for i, f in enumerate(b.freq):
            if not out_f:
                out_f.append(f)
                take.append(i)
                continue
            prev = out_f[-1]
            denom = prev if prev else 1.0
            if abs(f - prev) / denom > tol:
                out_f.append(f)
                take.append(i)
        b.freq = out_f
        if b.z is not None:
            b.z = [b.z[i] for i in take]
        if b.z_err is not None:
            b.z_err = [b.z_err[i] for i in take]
        if b.tipper is not None:
            b.tipper = [b.tipper[i] for i in take]
        if b.tipper_err is not None:
            b.tipper_err = [b.tipper_err[i] for i in take]
        if b.rho is not None:
            b.rho = [b.rho[i] for i in take]
        if b.phase is not None:
            b.phase = [b.phase[i] for i in take]
        return b

    def _fill_missing(self, b: TFBundle) -> TFBundle:
        cfg = get_config()
        have_z = b.z is not None
        have_rp = (b.rho is not None) and (b.phase is not None)
        if cfg.compute_res_from_z and have_z and not have_rp:
            b = self.compute_res_from_z(b)
        if cfg.compute_z_from_res and not have_z and have_rp:
            b = self.compute_z_from_res(b)
        return b

    def _finalize(
        self,
        b: TFBundle,
        *,
        name: Optional[str] = None,
        station_id: Optional[str | int] = None,
    ) -> TFBundle:
        sid = b.station_id if station_id is None else station_id
        b.station = ensure_station(name or b.station, sid)
        b = self._fill_missing(b)
        b = self._order_freq(b)
        b = self._dedup_freq(b)
        return b

    def transform(
        self,
        source: Any,
        *,
        name: Optional[str] = None,
        station_id: Optional[str | int] = None,
    ) -> Any:
        b = self.extract(source)
        b = self._finalize(
            b,
            name=name,
            station_id=station_id,
        )
        edi = self.emit_edi(b)
        return self.post_emit(edi, source, b)
