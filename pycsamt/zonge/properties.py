# -*- coding: utf-8 -*-
#       Author: LKouadio <etanoyau@gmail.com>
#       License: LGPL-3.0-or-later
"""pycsamt.zonge.properties
Reusable, lightweight property containers for survey‑level metadata
that are *not* tied to a particular AVG record.  These classes expose a
simple set/get API so downstream objects can update information in a
uniform way.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from datetime    import datetime, timezone
from pathlib     import Path
from typing      import Dict, Any, Optional
import json

__all__ = [
    "SkipFlag", "Hardware", "Receiver", "Transmitter",
    "SurveyConfiguration", "SurveyAnnotation"
]


class SkipFlag:
    """Encapsulate Zonge *skip‑flag* quality codes.

    Codes
    -----
    * 2 – good quality
    * 1 – kept but flagged as *skip*
    * 0 – rejected / bad data
    * * – no data
    """

    _flag_map: Dict[str, str] = {
        "2": "good", "1": "skip", "0": "reject", "*": "nodata"
    }

    def __init__(self, value: str | int | None = "2") -> None:
        self._value = "2"
        self.set(value)

    # set / get follow the same naming in all property classes
    def set(self, value: str | int | None = None) -> None:
        if value is None:
            return
        val = str(value)
        if val not in self._flag_map:
            valid = ", ".join(self._flag_map)
            raise ValueError(f"skip‑flag must be one of: {valid}")
        self._value = val

    def get(self) -> str:
        return self._flag_map[self._value]

    @property
    def code(self) -> str:
        return self._value

    def __str__(self) -> str:
        return f"SkipFlag({self.code} → {self.get()})"

    __repr__ = __str__


@dataclass(slots=True)
class Hardware:
    """Minimal information block taken from the AVG header."""

    version:        str            = "7.76"
    source_file:    Optional[Path]  = None
    dated:          Optional[str]   = None
    processed:      Optional[str]   = None
    astatic_ver:    str            = "v3.60"
    updated:        Optional[str]   = None
    tma_points:     Optional[int]   = None
    tma_freq:       Optional[float] = None
    _extra:         Dict[str, Any]  = field(default_factory=dict, repr=False)

    def set(self, **kwargs) -> None:
        for key, val in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, val)
            else:
                self._extra[key] = val

    def get(self, key: str, default: Any = None) -> Any:
        return getattr(self, key, self._extra.get(key, default))

    def __str__(self) -> str:
        src = f", file='{self.source_file.name}'" if self.source_file else ""
        return f"Hardware(ver={self.version}{src})"

    def to_json(self, *, indent: int = 0) -> str:
        base = asdict(self)
        base.update(self._extra)
        base["source_file"] = str(base["source_file"]) if base["source_file"] else None
        return json.dumps(base, indent=indent, default=str)

@dataclass(slots=True)
class Receiver:
    """Receiver electrode / coil (Rx) metadata."""

    station:         int   | None=None                    # client station number
    length_m:        float  | None =None                   # dipole length or loop width
    gdp_station:     int    | None = None      # Rx.GdpStn
    hpr:             tuple[float, float,
                           float] | None = None  # heading-pitch-roll
    comps:           str   | None = "ExHy"     # component label
    azimuth_deg:     float | None = None       # alias for heading
    latitude:        float | None = None
    longitude:       float | None = None
    elevation:       float | None = None
    unit:            str   | None = "m"        # length unit
    notes:           str   | None = None

    def set(self, **kwargs) -> None:
        for k, v in kwargs.items():
            setattr(self, k, v)

    def get(self, key: str, default: Any = None) -> Any:
        return getattr(self, key, default)

    def __str__(self) -> str:
        return (f"Receiver(stn={self.station}, len={self.length_m:g} "
                f"{self.unit or 'm'}, cmp={self.comps})")


@dataclass(slots=True)
class Transmitter:
    """Transmitter (Tx) loop / bipole metadata."""

    station:         int    | None = None  
    length_m:        float | None = None       # bipole length / loop width
    gdp_station:     int    | None = None      # Tx.GdpStn
    tx_type:         str   | None = None       # Natural, Bipole, Loop …
    center:          tuple[float, float,
                           float] | None = None  # (E, N, elev)
    hpr:             tuple[float, float,
                           float] | None = None  # heading-pitch-roll
    current_a:       float | None = None
    frequency_hz:    float | None = None
    latitude:        float | None = None
    longitude:       float | None = None
    notes:           str   | None = None

    def set(self, **kwargs) -> None:
        for k, v in kwargs.items():
            setattr(self, k, v)

    def get(self, key: str, default: Any = None) -> Any:
        return getattr(self, key, default)

    def __str__(self) -> str:
        cur = f", I={self.current_a:g} A" if self.current_a else ""
        typ = f", type={self.tx_type}"   if self.tx_type   else ""
        return f"Transmitter(stn={self.station}{cur}{typ})"


@dataclass(slots=True)
class SurveyConfiguration:
    """Survey-level configuration taken from AVG headers."""

    # core identifiers
    survey_type:   str            = "CSAMT"
    array_type:    str | None     = None

    # line information
    line_name:     str | None     = None
    line_number:   float | None   = None
    line_azim_deg: float | None   = None

    # GDP / client station numbering
    stn_gdp_beg:   float | None   = None
    stn_gdp_inc:   float | None   = None
    stn_beg:       float | None   = None
    stn_inc:       float | None   = None
    stn_left:      float | None   = None
    stn_right:     float | None   = None

    # units
    unit_length:   str            = "m"
    unit_emag:     str            = "nV/m"
    unit_hfield:   str            = "pT"
    unit_phase:    str            = "mrad"

    # misc
    utm_zone:      int  | None    = None
    created:       str            = field(
        default_factory=lambda:
            datetime.now().isoformat(timespec="seconds")
    )
    _extra:        Dict[str, Any] = field(default_factory=dict, repr=False)


    def set(self, **kwargs) -> None:
        for k, v in kwargs.items():
            if hasattr(self, k):
                setattr(self, k, v)
            else:
                self._extra[k] = v

    def get(self, key: str, default: Any = None) -> Any:
        return getattr(self, key, self._extra.get(key, default))


    def __str__(self) -> str:
        line = self.line_name or "-"
        return f"Survey(type={self.survey_type}, line={line})"

    def to_json(self, *, indent: int = 0) -> str:
        data = asdict(self)
        data.update(self._extra)
        return json.dumps(data, indent=indent)


@dataclass(slots=True)
class SurveyAnnotation:
    """Project-level annotation block (client, area, dates, …)."""

    project_name:    str            = "CSAMT_survey"
    project_area:    str | None     = None
    customer_name:   str            = "Zonge Engineering"
    contractor_name: str            = "Zonge"
    project_label:   str            = "pyCSAMT"
    acq_date:        str            = field(
        default_factory=lambda:
            datetime.now(tz=timezone.utc).isoformat(timespec="seconds")
    )
    _extra:          Dict[str, Any] = field(default_factory=dict,
                                           repr=False)

    def set(self, **kwargs) -> None:
        for k, v in kwargs.items():
            if hasattr(self, k):
                setattr(self, k, v)
            else:
                self._extra[k] = v

    def get(self, key: str, default: Any = None) -> Any:
        return getattr(self, key, self._extra.get(key, default))


    def __str__(self) -> str:
        area = self.project_area or "-"
        return (f"Annotation(project={self.project_name}, "
                f"area={area})")

    def to_json(self, *, indent: int = 0) -> str:
        data = asdict(self)
        data.update(self._extra)
        return json.dumps(data, indent=indent)
