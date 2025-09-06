# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

from typing import Dict, Any, Optional
from datetime import datetime

from .base import EdiComponentBase


def _to_float_or_none(value: Any) -> Optional[float]:
    """Best-effort float cast; returns None on blanks/None/invalid."""
    if value in (None, "", "None"):
        return None
    try:
        return float(value)
    except Exception:
        return None


def _norm_str(value: Any) -> Optional[str]:
    """Normalize strings: strip quotes/whitespace; None stays None."""
    if value is None:
        return None
    s = str(value).strip().strip('"').strip("'")
    return s if s != "" else None


class Hmeasurement(EdiComponentBase):
    """
    Magnetic field measurement metadata (HMEAS).

    Attributes (ordered for write-out)
    ----------------------------------
    id : str | None
    chtype : str | None             # HX | HY | HZ | RHX | RHY
    x, y, z : float | None
    azm : float | None
    dip : float | None
    acqchan : str | None
    filter : str | None
    sensor : str | None
    gain : float | None
    measdate : str | None           # UTC timestamp string
    """

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
        # Defaults
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
        self.measdate: Optional[str] = datetime.utcnow(
            ).strftime("%Y-%m-%d %H:%M:%S")

        if kws:
            self.update_from_dict(kws)

    # -----------------
    # Public utilities
    # -----------------
    def update_from_dict(self, d: Dict[str, Any]) -> "Hmeasurement":
        """Update fields from a (possibly stringy) dict."""
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

        # Validate chtype if present
        if self.chtype and self.chtype.upper() not in self._ALLOWED:
            # Keep the value but normalize; writer can still output it.
            self.chtype = self.chtype.upper()

        return self

    def to_dict(self) -> Dict[str, Any]:
        """Ordered dict per EDI write order."""
        return {k: getattr(self, k) for k in self.hmeasurementkey}


class Emeasurement(EdiComponentBase):
    """
    Electric field measurement metadata (EMEAS).

    Attributes (ordered for write-out)
    ----------------------------------
    id : str | None
    chtype : str | None             # EX | EY
    x, y, z : float | None          # first electrode reference
    x2, y2, z2 : float | None       # second electrode reference
    acqchan : str | None
    filter : str | None
    sensor : str | None
    gain : float | None
    measdate : str | None           # UTC timestamp string
    """

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
        # Defaults
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
        self.measdate: Optional[str] = datetime.utcnow(
            ).strftime("%Y-%m-%d %H:%M:%S")

        if kws:
            self.update_from_dict(kws)

    # -----------------
    # Public utilities
    # -----------------
    def update_from_dict(self, d: Dict[str, Any]) -> "Emeasurement":
        """Update fields from a (possibly stringy) dict."""
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

        # Validate chtype if present
        if self.chtype and self.chtype.upper() not in self._ALLOWED:
            # Keep the value but normalize; writer can still output it.
            self.chtype = self.chtype.upper()

        return self

    def to_dict(self) -> Dict[str, Any]:
        """Ordered dict per EDI write order."""
        return {k: getattr(self, k) for k in self.emeasurementkey}

