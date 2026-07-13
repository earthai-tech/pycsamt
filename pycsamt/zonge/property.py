# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.zonge.property
----------------------

Reusable, lightweight property containers for survey-level
metadata that are *not* tied to a particular AVG record.
All classes expose a simple ``set()/get()`` API, plus helpers
to read/write CSAVGW/ASTATIC ``$keyword=value`` style headers.

Design goals
~~~~~~~~~~~~
- Small, explicit data holders (dataclasses).
- Robust set/get with gentle validation.
- Round-trip helpers:
    • ``update_from_keywords(meta)``  ← dict from parser
    • ``to_keywords()``               → dict for writer
- Tolerate legacy aliases (e.g. ``XMTR`` ↔ ``Tx.GdpStn``).
- Keep lines short for readability in code reviews.
"""

from __future__ import annotations

import json
from dataclasses import asdict, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from ..compat.python import dc

__all__ = [
    "SkipFlag",
    "Hardware",
    "Receiver",
    "Transmitter",
    "SurveyConfiguration",
    "SurveyAnnotation",
]


class SkipFlag:
    """
    Encapsulate Zonge *skip-flag* quality codes.

    Codes
    -----
    2 → good quality
    1 → kept but *skip* on plots
    0 → rejected / bad data
    * → no data (placeholder)
    """

    _flag_map: dict[str, str] = {
        "2": "good",
        "1": "skip",
        "0": "reject",
        "*": "nodata",
    }

    def __init__(self, value: str | int | None = "2") -> None:
        self._value = "2"
        self.set(value)

    def set(self, value: str | int | None = None) -> None:
        """Set the flag using {0,1,2,'*'}."""
        if value is None:
            return
        code = str(value)
        if code not in self._flag_map:
            valid = ", ".join(self._flag_map)
            raise ValueError(f"skip-flag must be one of: {valid}")
        self._value = code

    def get(self) -> str:
        """Return the *label* ('good', 'skip', ...)."""
        return self._flag_map[self._value]

    @property
    def code(self) -> str:
        """Return the *raw* code as string."""
        return self._value

    def __str__(self) -> str:
        return f"SkipFlag({self.code} → {self.get()})"

    __repr__ = __str__


@dc(slots=True)
class Hardware:
    """
    Minimal provenance captured from banner / comment lines.

    Many of these are informational and may not appear as $keywords.
    We keep ``to_keywords`` minimal to avoid inventing conventions.
    """

    version: str = "7.76"
    source_file: Path | None = None
    dated: str | None = None
    processed: str | None = None
    astatic_ver: str = "v3.60"
    updated: str | None = None
    tma_points: int | None = None
    tma_freq: float | None = None
    _extra: dict[str, Any] = field(default_factory=dict, repr=False)

    # No formal CSAVGW keymap for banner lines; keep empty by design.
    KEYMAP: dict[str, str] = field(default_factory=dict, init=False)

    def set(self, **kwargs) -> None:
        """Set known fields; unknowns land in ``_extra``."""
        for key, val in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, val)
            else:
                self._extra[key] = val

    def get(self, key: str, default: Any = None) -> Any:
        """Get known field or fallback to ``_extra``."""
        return getattr(self, key, self._extra.get(key, default))

    def __str__(self) -> str:
        src = f", file='{self.source_file.name}'" if self.source_file else ""
        return f"Hardware(ver={self.version}{src})"

    def to_json(self, *, indent: int = 0) -> str:
        base = asdict(self)
        base.update(self._extra)
        if base["source_file"] is not None:
            base["source_file"] = str(base["source_file"])
        return json.dumps(base, indent=indent, default=str)

    def update_from_keywords(self, meta: dict[str, Any]) -> None:
        """
        Accept a plain dict (keys may or may not start with '$') and
        coerce a few well-known banner values to proper types.
        """
        m = {_norm_key(k): v for k, v in meta.items()}

        # source_file handling so tests can round-trip it
        if "source_file" in m:
            try:
                self.source_file = Path(str(m["source_file"]).strip())
            except Exception:
                # Fall back to a plain string if Path chokes
                self.source_file = Path(str(m["source_file"]))

        if "version" in m:
            self.version = str(m["version"]).strip()
        if "dated" in m:
            self.dated = str(m["dated"]).strip()
        if "processed" in m:
            self.processed = str(m["processed"]).strip()
        if "astatic_ver" in m:
            self.astatic_ver = str(m["astatic_ver"]).strip()
        if "updated" in m:
            self.updated = str(m["updated"]).strip()

        if "tma_points" in m:
            tp = _to_number(m["tma_points"])
            self.tma_points = (
                int(tp) if isinstance(tp, (int, float)) else None
            )
        if "tma_freq" in m:
            tf = _to_number(m["tma_freq"])
            self.tma_freq = (
                float(tf) if isinstance(tf, (int, float)) else None
            )

    def to_keywords(self) -> dict[str, Any]:
        """
        Provide a small, stable set for round-trip / tests.
        """
        out: dict[str, Any] = {}
        if self.version is not None:
            out["version"] = self.version
        if self.dated is not None:
            out["dated"] = self.dated
        if self.processed is not None:
            out["processed"] = self.processed
        if self.astatic_ver is not None:
            out["astatic_ver"] = self.astatic_ver
        if self.updated is not None:
            out["updated"] = self.updated
        if self.tma_points is not None:
            out["tma_points"] = int(self.tma_points)
        if self.tma_freq is not None:
            out["tma_freq"] = float(self.tma_freq)

        out["source_file"] = (
            str(self.source_file) if self.source_file is not None else None
        )

        return out


@dc(slots=True)
class Receiver:
    """Receiver electrode / coil metadata (``$Rx.*``)."""

    station: int | None = None  # $Rx.Stn  (client)
    length_m: float | None = None  # $Rx.Length (units)
    gdp_station: int | None = None  # $Rx.GdpStn
    hpr: tuple[float, float, float] | None = None  # $Rx.HPR
    comps: str | None = "ExHy"  # $Rx.Cmp
    azimuth_deg: float | None = None  # convenience (alias of HPR[0])
    latitude: float | None = None  # $GPS.Lat
    longitude: float | None = None  # $GPS.Lon
    elevation: float | None = None  # (rare; when present)
    unit: str | None = "m"  # $Unit.Length
    notes: str | None = None

    KEYMAP: dict[str, str] = field(
        default_factory=lambda: {
            "station": "Rx.Stn",
            "gdp_station": "Rx.GdpStn",
            "length_m": "Rx.Length",
            "hpr": "Rx.HPR",
            "comps": "Rx.Cmp",
            "latitude": "GPS.Lat",
            "longitude": "GPS.Lon",
        },
        init=False,
    )
    ALIASES: dict[str, str] = field(default_factory=dict, init=False)

    def set(self, **kwargs) -> None:
        """Set any attribute; no hard validation here."""
        for k, v in kwargs.items():
            setattr(self, k, v)

    def get(self, key: str, default: Any = None) -> Any:
        """Get attribute by name."""
        return getattr(self, key, default)

    def update_from_keywords(self, meta: dict[str, Any]) -> None:
        """Update from parsed $Rx.* / $GPS.* keys with typing."""
        # 1) normalize keys once (“$Rx.Stn” → “Rx.Stn”)
        norm = {_norm_key(k): v for k, v in meta.items()}

        # 2) apply mapping using normalized keys
        _apply_keywords(self, self.KEYMAP, norm, aliases=self.ALIASES)

        # 3) numeric coercions from normalized keys
        stn = _to_number(norm.get("Rx.Stn"))
        if isinstance(stn, (int, float)):
            self.station = int(stn)

        gdp = _to_number(norm.get("Rx.GdpStn"))
        if isinstance(gdp, (int, float)):
            self.gdp_station = int(gdp)

        length, unit = _parse_length(norm.get("Rx.Length"), self.unit or "m")
        if length is not None:
            self.length_m = float(length)
        self.unit = unit or self.unit

        hpr_raw = norm.get("Rx.HPR")
        if hpr_raw is not None:
            if isinstance(hpr_raw, (tuple, list)) and len(hpr_raw) == 3:
                self.hpr = tuple(float(_to_number(v) or 0.0) for v in hpr_raw)  # type: ignore
            else:
                parts = [
                    p.strip()
                    for p in str(hpr_raw).replace(";", ",").split(",")
                ]
                if len(parts) >= 3:
                    self.hpr = (
                        float(_to_number(parts[0]) or 0.0),
                        float(_to_number(parts[1]) or 0.0),
                        float(_to_number(parts[2]) or 0.0),
                    )

        if self.hpr and self.azimuth_deg is None:
            self.azimuth_deg = float(self.hpr[0])

        lat = _to_number(norm.get("GPS.Lat"))
        lon = _to_number(norm.get("GPS.Lon"))
        if isinstance(lat, (int, float)):
            self.latitude = float(lat)
        if isinstance(lon, (int, float)):
            self.longitude = float(lon)

        if "Rx.Cmp" in norm and norm["Rx.Cmp"] is not None:
            self.comps = str(norm["Rx.Cmp"])

    def to_keywords(self) -> dict[str, Any]:
        """Export standardized $Rx.* keys with units preserved."""
        out: dict[str, Any] = {}
        if self.gdp_station is not None:
            out["Rx.GdpStn"] = int(self.gdp_station)
        if self.station is not None:
            out["Rx.Stn"] = int(self.station)
        if self.length_m is not None:
            out["Rx.Length"] = f"{self.length_m:g} {self.unit or 'm'}"
        if self.hpr is not None:
            h, p, r = self.hpr
            out["Rx.HPR"] = f"{h:g},{p:g},{r:g}"
        if self.comps:
            out["Rx.Cmp"] = self.comps
        if self.latitude is not None:
            out["GPS.Lat"] = float(self.latitude)
        if self.longitude is not None:
            out["GPS.Lon"] = float(self.longitude)
        return out

    def __str__(self) -> str:
        unit = self.unit or "m"
        st = f"{self.station}" if self.station is not None else "?"
        ln = f"{self.length_m:g} {unit}" if self.length_m else "?"
        cmp_ = self.comps or "?"
        return f"Receiver(stn={st}, len={ln}, cmp={cmp_})"


@dc(slots=True)
class Transmitter:
    """Transmitter loop / bipole metadata (``$Tx.*``)."""

    station: int | None = None  # $Tx.Stn  (client)
    length_m: float | None = None  # $Tx.Length
    gdp_station: int | None = None  # $Tx.GdpStn  (alias XMTR)
    tx_type: str | None = None  # $Tx.Type
    center: tuple[float, float, float] | None = None  # $Tx.Center
    hpr: tuple[float, float, float] | None = None  # $Tx.HPR
    current_a: float | None = None
    frequency_hz: float | None = None
    latitude: float | None = None  # seldom in headers
    longitude: float | None = None
    notes: str | None = None

    KEYMAP: dict[str, str] = field(
        default_factory=lambda: {
            "station": "Tx.Stn",
            "gdp_station": "Tx.GdpStn",
            "length_m": "Tx.Length",
            "tx_type": "Tx.Type",
            "center": "Tx.Center",
            "hpr": "Tx.HPR",
        },
        init=False,
    )
    ALIASES: dict[str, str] = field(
        default_factory=lambda: {
            # legacy alias encountered in AMTAVG headers
            "XMTR": "gdp_station",
        },
        init=False,
    )

    def set(self, **kwargs) -> None:
        """Set any attribute; no hard validation here."""
        for k, v in kwargs.items():
            setattr(self, k, v)

    def get(self, key: str, default: Any = None) -> Any:
        """Get attribute by name."""
        return getattr(self, key, default)

    def update_from_keywords(self, meta: dict[str, Any]) -> None:
        """Update from $Tx.* (and legacy XMTR) with typing."""
        # 1) normalize keys once
        norm = {_norm_key(k): v for k, v in meta.items()}

        # 2) apply mapping with normalized keys
        _apply_keywords(self, self.KEYMAP, norm)

        # 3) numeric coercions from normalized keys
        stn = _to_number(norm.get("Tx.Stn"))
        if isinstance(stn, (int, float)):
            self.station = int(stn)

        gdp = _to_number(norm.get("Tx.GdpStn") or norm.get("XMTR"))
        if isinstance(gdp, (int, float)):
            self.gdp_station = int(gdp)
        elif "XMTR" in norm and self.gdp_station is None:
            try:
                self.gdp_station = int(norm["XMTR"])  # type: ignore[arg-type]
            except Exception:
                self.gdp_station = norm[
                    "XMTR"
                ]  # leave as string if truly odd

        length, _ = _parse_length(norm.get("Tx.Length"), "m")
        if length is not None:
            self.length_m = float(length)

        if "Tx.Type" in norm and norm["Tx.Type"] is not None:
            self.tx_type = str(norm["Tx.Type"])

        cen = norm.get("Tx.Center")
        if cen is not None:
            if isinstance(cen, (tuple, list)) and len(cen) == 3:
                self.center = (
                    float(_to_number(cen[0]) or 0.0),
                    float(_to_number(cen[1]) or 0.0),
                    float(_to_number(cen[2]) or 0.0),
                )
            else:
                parts = [
                    p.strip() for p in str(cen).replace(";", ",").split(",")
                ]
                if len(parts) >= 3:
                    self.center = (
                        float(_to_number(parts[0]) or 0.0),
                        float(_to_number(parts[1]) or 0.0),
                        float(_to_number(parts[2]) or 0.0),
                    )

        hpr_raw = norm.get("Tx.HPR")
        if hpr_raw is not None:
            if isinstance(hpr_raw, (tuple, list)) and len(hpr_raw) == 3:
                self.hpr = tuple(float(_to_number(v) or 0.0) for v in hpr_raw)  # type: ignore
            else:
                parts = [
                    p.strip()
                    for p in str(hpr_raw).replace(";", ",").split(",")
                ]
                if len(parts) >= 3:
                    self.hpr = (
                        float(_to_number(parts[0]) or 0.0),
                        float(_to_number(parts[1]) or 0.0),
                        float(_to_number(parts[2]) or 0.0),
                    )

    def to_keywords(self) -> dict[str, Any]:
        """Export standardized $Tx.* keys."""
        out: dict[str, Any] = {}
        if self.gdp_station is not None:
            out["Tx.GdpStn"] = int(self.gdp_station)
        if self.station is not None:
            out["Tx.Stn"] = int(self.station)
        if self.length_m is not None:
            out["Tx.Length"] = f"{self.length_m:g} m"
        if self.tx_type:
            out["Tx.Type"] = self.tx_type
        if self.center is not None:
            e, n, z = self.center
            out["Tx.Center"] = f"{e:g},{n:g},{z:g}"
        if self.hpr is not None:
            h, p, r = self.hpr
            out["Tx.HPR"] = f"{h:g},{p:g},{r:g}"
        return out

    def __str__(self) -> str:
        cur = f", I={self.current_a:g} A" if self.current_a else ""
        typ = f", type={self.tx_type}" if self.tx_type else ""
        stn = "?" if self.station is None else f"{self.station}"
        return f"Transmitter(stn={stn}{cur}{typ})"


@dc(slots=True)
class SurveyConfiguration:
    """Survey-level configuration taken from AVG headers."""

    # core identifiers
    survey_type: str = "CSAMT"  # $Survey.Type
    array_type: str | None = None  # $Survey.Array

    # line information
    line_name: str | None = None  # $Line.Name
    line_number: float | None = None  # $Line.Number
    line_azim_deg: float | None = None  # $Line.Azimuth

    # GDP / client station numbering
    stn_gdp_beg: float | None = None  # $Stn.GdpBeg
    stn_gdp_inc: float | None = None  # $Stn.GdpInc
    stn_beg: float | None = None  # $Stn.Beg
    stn_inc: float | None = None  # $Stn.Inc
    stn_left: float | None = None  # $Stn.Left
    stn_right: float | None = None  # $Stn.Right

    # units
    unit_length: str = "m"  # $Unit.Length
    unit_emag: str = "nV/Am"  # $Unit.E
    unit_hfield: str = "pT/A"  # $Unit.B
    unit_phase: str = "mrad"  # $Unit.Phase

    # misc
    utm_zone: int | None = None
    created: str = field(
        default_factory=lambda: datetime.now().isoformat(timespec="seconds")
    )
    _extra: dict[str, Any] = field(default_factory=dict, repr=False)

    KEYMAP: dict[str, str] = field(
        default_factory=lambda: {
            "survey_type": "Survey.Type",
            "array_type": "Survey.Array",
            "line_name": "Line.Name",
            "line_number": "Line.Number",
            "line_azim_deg": "Line.Azimuth",
            "stn_gdp_beg": "Stn.GdpBeg",
            "stn_gdp_inc": "Stn.GdpInc",
            "stn_beg": "Stn.Beg",
            "stn_inc": "Stn.Inc",
            "stn_left": "Stn.Left",
            "stn_right": "Stn.Right",
            "unit_length": "Unit.Length",
            "unit_emag": "Unit.E",
            "unit_hfield": "Unit.B",
            "unit_phase": "Unit.Phase",
        },
        init=False,
    )
    ALIASES: dict[str, str] = field(
        default_factory=lambda: {
            # common synonyms observed in real files
            "Units": "Unit.Length",
            "JobLine": "Line.Name",
            "BrgLine": "Line.Azimuth",
            "StnBeg": "Stn.GdpBeg",
            "StnDelt": "Stn.GdpInc",
            "LblFrst": "Stn.Left",
            "LblDelt": "Stn.Inc",
        },
        init=False,
    )

    def set(self, **kwargs) -> None:
        """
        Set known fields; unknowns are stored in ``_extra`` so
        nothing is lost on round-trip.
        """
        for k, v in kwargs.items():
            if hasattr(self, k):
                setattr(self, k, v)
            else:
                self._extra[k] = v

    def get(self, key: str, default: Any = None) -> Any:
        """Get known field or fallback to ``_extra``."""
        return getattr(self, key, self._extra.get(key, default))

    def update_from_keywords(self, meta: dict[str, Any]) -> None:
        """
        Populate from $... header dict (keys without '$').
        Coerce numeric fields so tests don't see strings.
        """
        # Merge aliases → external canonical keys
        unified: dict[str, Any] = {}
        for k, v in meta.items():
            key = _norm_key(k)
            key = self.ALIASES.get(key, key)
            unified[key] = v

        # Map into attributes
        _apply_keywords(self, self.KEYMAP, unified)

        # Coerce numerics
        for attr in (
            "line_number",
            "line_azim_deg",
            "stn_gdp_beg",
            "stn_gdp_inc",
            "stn_beg",
            "stn_inc",
            "stn_left",
            "stn_right",
        ):
            val = getattr(self, attr)
            num = _to_number(val)
            if isinstance(num, (int, float)):
                setattr(self, attr, float(num))

    def to_keywords(self) -> dict[str, Any]:
        """Export standardized $Survey.*, $Line.*, $Stn.*."""
        base = asdict(self)
        base.update(self._extra)
        # prune internals
        base.pop("_extra", None)
        base.pop("utm_zone", None)
        base.pop("created", None)
        # leave floats as floats (no stringification)
        return _kv_roundtrip(base, self.KEYMAP)

    def __str__(self) -> str:
        line = self.line_name or "-"
        return f"Survey(type={self.survey_type}, line={line})"

    def to_json(self, *, indent: int = 0) -> str:
        data = asdict(self)
        data.update(self._extra)
        return json.dumps(data, indent=indent)


@dc(slots=True)
class SurveyAnnotation:
    """Project-level annotation block (``$Job.*``)."""

    project_name: str = "CSAMTSurvey"  # $Job.Name
    project_area: str | None = None  # $Job.Area
    customer_name: str = "Zonge Engineering"  # $Job.For
    contractor_name: str = "Zonge"  # $Job.By
    project_label: str = "pyCSAMT"  # $Job.Number
    acq_date: str = field(  # $Job.Date
        default_factory=lambda: datetime.now(tz=timezone.utc).isoformat(
            timespec="seconds"
        )
    )
    _extra: dict[str, Any] = field(default_factory=dict, repr=False)

    KEYMAP: dict[str, str] = field(
        default_factory=lambda: {
            "project_name": "Job.Name",
            "project_area": "Job.Area",
            "customer_name": "Job.For",
            "contractor_name": "Job.By",
            "project_label": "Job.Number",
            "acq_date": "Job.Date",
        },
        init=False,
    )
    ALIASES: dict[str, str] = field(
        default_factory=lambda: {
            "Project": "Job.Name",
            "client": "Job.For",
            "company": "Job.By",
            "JobNumb": "Job.Number",
            "JobDate": "Job.Date",
        },
        init=False,
    )

    def set(self, **kwargs) -> None:
        """Set known fields; unknowns land in ``_extra``."""
        for k, v in kwargs.items():
            if hasattr(self, k):
                setattr(self, k, v)
            else:
                self._extra[k] = v

    def get(self, key: str, default: Any = None) -> Any:
        """Get known field or fallback to ``_extra``."""
        return getattr(self, key, self._extra.get(key, default))

    def update_from_keywords(self, meta: dict[str, Any]) -> None:
        """Update from parsed ``$Job.*`` keys (with alias support)."""
        unified: dict[str, Any] = {}
        for k, v in meta.items():
            key = _norm_key(k)
            key = self.ALIASES.get(key, key)
            unified[key] = v
        _apply_keywords(self, self.KEYMAP, unified)

    def to_keywords(self) -> dict[str, Any]:
        """Export standardized ``$Job.*`` keys."""
        base = asdict(self)
        base.update(self._extra)
        base.pop("_extra", None)
        return _kv_roundtrip(base, self.KEYMAP)

    def __str__(self) -> str:
        area = self.project_area or "-"
        return f"Annotation(project={self.project_name}, area={area})"

    def to_json(self, *, indent: int = 0) -> str:
        data = asdict(self)
        data.update(self._extra)
        return json.dumps(data, indent=indent)


def _norm_key(key: str) -> str:
    """
    Normalize header keys:

    - strips leading '$'
    - preserves section (e.g., 'Rx.Length')
    - trims surrounding whitespace
    """
    key = key.lstrip("$").strip()
    return key


def _kv_roundtrip(
    data: dict[str, Any],
    keymap: dict[str, str],
    *,
    drop_none: bool = True,
) -> dict[str, Any]:
    """
    Map internal attribute names → CSAVGW keys using *keymap*.

    Parameters
    ----------
    data
        Internal ``{attr: value}`` mapping (e.g., from dataclass).
    keymap
        Mapping of internal attr → "Section.Key" strings.
    drop_none
        When *True*, omit entries with ``None`` values.

    Returns
    -------
    dict
        ``{"Section.Key": value}`` dictionary suitable for writing
        back to ``.avg`` headers.
    """
    out: dict[str, Any] = {}
    for attr, ext in keymap.items():
        if attr not in data:
            continue
        val = data[attr]
        if drop_none and val is None:
            continue
        out[ext] = val
    return out


def _apply_keywords(
    obj: Any,
    keymap: dict[str, str],
    meta: dict[str, Any],
    *,
    aliases: dict[str, str] | None = None,
) -> None:
    """
    Update *obj* attributes from a parsed metadata dict.

    Notes
    -----
    - Keys may come with or without leading '$'.
    - *aliases* supports legacy synonyms (e.g., 'XMTR').
    """
    inv = {v: k for k, v in keymap.items()}
    if aliases:
        inv.update({_norm_key(k): v for k, v in aliases.items()})

    for k, v in meta.items():
        key = _norm_key(k)
        attr = inv.get(key)
        if attr is None:
            continue
        setattr(obj, attr, v)


def _to_number(x):
    """int if integral, else float, else original."""
    if x is None:
        return None
    s = str(x).strip()
    # empty/asterisk/null → None
    if s in {"", "*", "nan", "NaN", "None", "null"}:
        return None
    try:
        f = float(s)
        return int(f) if f.is_integer() else f
    except Exception:
        return x  # leave as string if not numeric


def _parse_length(s, default_unit="m"):
    """
    Parse '200 ft' or '50 m' → (value: float, unit: str).
    If only '200' provided, keep default_unit.
    """
    if s is None:
        return None, default_unit
    t = str(s).strip()
    parts = t.split()
    if not parts:
        return None, default_unit
    value = _to_number(parts[0])
    unit = parts[1] if len(parts) > 1 else default_unit
    # coerce unit spacing/format just a bit
    unit = unit.replace(".", "")
    return (float(value) if isinstance(value, (int, float)) else None, unit)
