# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

from datetime import datetime
from pathlib import Path
import re
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

from ..gis.utils import dms_to_decimal, decimal_to_dms
from ..exceptions import (
    FileHandlingError,
    EdIDataError,
    HeaderError,
)
from .base import EDIComponentBase
from .validation import ( 
    _to_float_or_none, 
    _is_tag, 
    _extract_tag, 
    _norm_str, 
    IsEdi, 
)

__all__= [ 
    "Hmeasurement", 
    "Emeasurement", 
    "DefineMeas", 
    "MeasMixin", 
    "EMeasMixin", 
    "DefineMeasMixin", 
    "HMeasMixin"
    ]


def _slice_section(
    lines: Sequence[str],
    start_tag: str,
    after_tags: Optional[Sequence[str]] = None,
) -> Tuple[List[str], int, int]:
    """
    Return payload between `start_tag` and the next section begin
    (a line that starts with '>=' and is not `start_tag`).
    """
    after_tags = list(after_tags or [])
    start = None
    for i, ln in enumerate(lines):
        if _is_tag(ln, start_tag):
            start = i
            break
    if start is None:
        raise EdIDataError(f"Missing section: {start_tag!r}")

    j = start + 1
    stop = None
    while j < len(lines):
        t = _extract_tag(lines[j])
        if t and t.startswith(">=") and not t.startswith(start_tag):
            stop = j
            break
        # allow measurement lines (>HMEAS / >EMEAS) and comments (>!)
        j += 1
    if stop is None:
        stop = len(lines)

    # if after_tags are provided, prefer them as hard stops
    if after_tags:
        for k in range(start + 1, len(lines)):
            for a in after_tags:
                if _is_tag(lines[k], a):
                    stop = min(stop, k)
                    break

    payload = [ln.rstrip("\n") for ln in lines[start + 1 : stop]]
    return payload, start, stop


# Measurement line parsing

# TOKEN =  key=value   (value may contain +/-, dot, colon or text)
_TOKEN_RE = re.compile(
    r"(?P<k>[A-Za-z][A-Za-z0-9_]*)\s*=\s*(?P<v>[^\s]+)"
)


def _kv_tokens_from_line(line: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for m in _TOKEN_RE.finditer(line):
        k = m.group("k").strip().lower()
        v = m.group("v").strip().strip('"').strip("'")
        out[k] = v
    return out



class Hmeasurement(EDIComponentBase):
    # Magnetic field measurement (HMEAS)
    hmeasurementkey = [
        "id",
        "chtype",
        "x",
        "y",
        "z",
        "azm",
        "dip",
        "acqchan",
        "filter",
        "sensor",
        "gain",
        "measdate",
    ]
    _ALLOWED = {"HX", "HY", "HZ", "RHX", "RHY"}

    def __init__(self, **kws: Any):
        self.id: Optional[str] = None
        self.chtype: Optional[str] = None
        self.x: Optional[float] = 0.0
        self.y: Optional[float] = 0.0
        self.z: Optional[float] = 0.0
        self.azm: Optional[float] = 0.0
        self.dip: Optional[float] = None
        self.acqchan: Optional[str] = None
        self.filter: Optional[str] = None
        self.sensor: Optional[str] = None
        self.gain: Optional[float] = None
        self.measdate: Optional[str] = datetime.utcnow().strftime(
            "%Y-%m-%d %H:%M:%S"
        )
        if kws:
            self.update_from_dict(kws)

    def update_from_dict(self, d: Dict[str, Any]) -> "Hmeasurement":
        for k, v in d.items():
            kk = str(k).lower()
            if kk == "id":
                self.id = _norm_str(v)
            elif kk == "chtype":
                vv = _norm_str(v)
                self.chtype = vv.upper() if vv else None
            elif kk in ("x", "y", "z", "azm", "dip", "gain"):
                setattr(self, kk, _to_float_or_none(v))
            elif kk in ("acqchan", "filter", "sensor", "measdate"):
                setattr(self, kk, _norm_str(v))
        if self.chtype and self.chtype.upper() not in self._ALLOWED:
            self.chtype = self.chtype.upper()
        return self

    @classmethod
    def from_line(cls, line: str) -> "Hmeasurement":
        kv = _kv_tokens_from_line(line)
        return cls(**kv)

    def to_dict(self) -> Dict[str, Any]:
        return {k: getattr(self, k) for k in self.hmeasurementkey}

    def to_line(self) -> str:
        parts: List[str] = []
        for k in self.hmeasurementkey:
            v = getattr(self, k, None)
            if v in (None, "", "None"):
                continue
            parts.append(f"{k.upper()}={v}")
        return ">HMEAS " + " ".join(parts)


class Emeasurement(EDIComponentBase):
    # Electric field measurement (EMEAS)
    emeasurementkey = [
        "id",
        "chtype",
        "x",
        "y",
        "z",
        "x2",
        "y2",
        "z2",
        "acqchan",
        "filter",
        "sensor",
        "gain",
        "measdate",
    ]
    _ALLOWED = {"EX", "EY"}

    def __init__(self, **kws: Any):
        self.id: Optional[str] = None
        self.chtype: Optional[str] = None
        self.x: Optional[float] = 0.0
        self.y: Optional[float] = 0.0
        self.z: Optional[float] = 0.0
        self.x2: Optional[float] = 0.0
        self.y2: Optional[float] = 0.0
        self.z2: Optional[float] = 0.0
        self.acqchan: Optional[str] = None
        self.filter: Optional[str] = None
        self.sensor: Optional[str] = None
        self.gain: Optional[float] = None
        self.measdate: Optional[str] = datetime.utcnow().strftime(
            "%Y-%m-%d %H:%M:%S"
        )
        if kws:
            self.update_from_dict(kws)

    def update_from_dict(self, d: Dict[str, Any]) -> "Emeasurement":
        for k, v in d.items():
            kk = str(k).lower()
            if kk == "id":
                self.id = _norm_str(v)
            elif kk == "chtype":
                vv = _norm_str(v)
                self.chtype = vv.upper() if vv else None
            elif kk in ("x", "y", "z", "x2", "y2", "z2", "gain"):
                setattr(self, kk, _to_float_or_none(v))
            elif kk in ("acqchan", "filter", "sensor", "measdate"):
                setattr(self, kk, _norm_str(v))
        if self.chtype and self.chtype.upper() not in self._ALLOWED:
            self.chtype = self.chtype.upper()
        return self

    @classmethod
    def from_line(cls, line: str) -> "Emeasurement":
        kv = _kv_tokens_from_line(line)
        return cls(**kv)

    def to_dict(self) -> Dict[str, Any]:
        return {k: getattr(self, k) for k in self.emeasurementkey}

    def to_line(self) -> str:
        parts: List[str] = []
        for k in self.emeasurementkey:
            v = getattr(self, k, None)
            if v in (None, "", "None"):
                continue
            parts.append(f"{k.upper()}={v}")
        return ">EMEAS " + " ".join(parts)


class DefineMeas(EDIComponentBase):
    """
    Container for >=DEFINEMEAS payload including reference origin
    and lists of HMEAS / EMEAS.
    """

    def __init__(self, **kws: Any):
        self.maxchan: Optional[int] = None
        self.maxrun: Optional[int] = None
        self.maxmeas: Optional[int] = None
        self.units: Optional[str] = None
        self.reftype: Optional[str] = None

        self.reflat: Optional[float] = None
        self.reflong: Optional[float] = None
        self.refelev: Optional[float] = None

        self.hmeas: List[Hmeasurement] = []
        self.emeas: List[Emeasurement] = []

        if kws:
            for k, v in kws.items():
                setattr(self, str(k).lower(), v)

    @classmethod
    def from_file(cls, edi_fn: Union[str, Path]) -> "DefineMeas":
        if edi_fn is None:
            raise FileHandlingError("No EDI path provided.")
        p = Path(edi_fn)
        IsEdi._assert_edi(p, deep=True)

        lines = p.read_text(
            encoding="utf-8-sig", errors="replace"
        ).splitlines()

        payload, _, _ = _slice_section(
            lines, ">=DEFINEMEAS", after_tags=[">=MTSECT", ">=SPECTRASECT"]
        )
        obj = cls()
        return obj.read(payload)

    def read(self, lines: Sequence[str]) -> "DefineMeas":
        if not lines:
            raise HeaderError("Empty DEFINEMEAS payload.")

        for raw in lines:
            s = raw.strip()
            if not s:
                continue

            if _is_tag(s, ">HMEAS"):
                self.hmeas.append(Hmeasurement.from_line(s))
                continue

            if _is_tag(s, ">EMEAS"):
                self.emeas.append(Emeasurement.from_line(s))
                continue

            # plain KV lines (MAX*, UNITS, REF*, etc.)
            m = _TOKEN_RE.match(s)
            if not m:
                # allow comments like >!****
                continue

            kv = _kv_tokens_from_line(s)
            for k, v in kv.items():
                kk = k.lower()
                if kk in {"maxchan", "maxrun", "maxmeas"}:
                    try:
                        setattr(self, kk, int(float(v)))
                    except Exception:
                        setattr(self, kk, None)
                elif kk in {"units", "reftype"}:
                    setattr(self, kk, v.upper())
                elif kk == "reflat":
                    self.reflat = dms_to_decimal(v)
                elif kk == "reflong":
                    self.reflong = dms_to_decimal(v)
                elif kk == "refelev":
                    try:
                        self.refelev = float(v)
                    except Exception:
                        self.refelev = None
                # ignore unknown keys silently

        return self

    def write(self) -> List[str]:
        out: List[str] = [">=DEFINEMEAS\n"]
        if self.maxchan is not None:
            out.append(f"  MAXCHAN={self.maxchan}\n")
        if self.maxrun is not None:
            out.append(f"  MAXRUN={self.maxrun}\n")
        if self.maxmeas is not None:
            out.append(f"  MAXMEAS={self.maxmeas}\n")
        if self.units:
            out.append(f"  UNITS={self.units}\n")
        if self.reftype:
            out.append(f"  REFTYPE={self.reftype}\n")

        # optional reference origin (write DMS if numeric)
        if self.reflat is not None:
            out.append(f"  REFLAT={decimal_to_dms(float(self.reflat))}\n")
        if self.reflong is not None:
            out.append(f"  REFLONG={decimal_to_dms(float(self.reflong))}\n")
        if self.refelev is not None:
            out.append(f"  REFELEV={self.refelev}\n")

        # measurement lines
        for hm in self.hmeas:
            out.append("  " + hm.to_line() + "\n")
        for em in self.emeas:
            out.append("  " + em.to_line() + "\n")

        out.append("\n")
        return out

    # convenience
    def all_meas(self) -> List[Union[Hmeasurement, Emeasurement]]:
        return list(self.hmeas) + list(self.emeas)


class MeasMixin:
    """
    Facade: provide one-call access to parsed DefineMeas.
    """

    @classmethod
    def from_file(cls, edi_fn: Union[str, Path]) -> DefineMeas:
        return DefineMeas.from_file(edi_fn)


class DefineMeasMixin:
    """
    Facade for >=DEFINEMEAS container only.
    """

    @classmethod
    def from_file(cls, edi_fn: Union[str, Path]) -> DefineMeas:
        return DefineMeas.from_file(edi_fn)


class EMeasMixin:
    """
    Facade that returns only EMEAS list from file.
    """

    @classmethod
    def from_file(
        cls, edi_fn: Union[str, Path]
    ) -> List[Emeasurement]:
        dm = DefineMeas.from_file(edi_fn)
        return dm.emeas


class HMeasMixin:
    """
    Facade that returns only HMEAS list from file.
    """

    @classmethod
    def from_file(
        cls, edi_fn: Union[str, Path]
    ) -> List[Hmeasurement]:
        dm = DefineMeas.from_file(edi_fn)
        return dm.hmeas
