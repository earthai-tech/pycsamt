# -*- coding: utf-8 -*-
"""
Site‑level properties parsed from the J‑format information block.

This module defines a lightweight container for the metadata
exposed in lines beginning with ``>KEY=VALUE`` in A.G. Jones J files.
It focuses on robust parsing and normalization (e.g., azimuth
wrapping, latitude/longitude hemispheres, and missing sentinels).

This revision leverages utilities from :mod:`pycsamt.gis.utils`
for robust coordinate parsing/validation (DMS handling,
hemisphere suffixes, and range checks).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Mapping, Tuple
import io
import math
import warnings

from ..gis.utils import (
    convert_position_str2float as _to_deg,  
    assert_lat_value as _assert_lat,        
    assert_lon_value as _assert_lon,        
)

from .config import (
    ENCODING_DEFAULT,
    MISSING_FLOAT,
)
from .utils import iter_lines, iter_info, JParseError


def _to_float_safe(s: str) -> float:
    """Best‑effort float parser for simple numeric fields (e.g., azimuth).

    Unlike :func:`_to_deg`, this does *not* parse DMS — it's meant for
    scalar numeric tokens where non‑numeric input should surface quickly.
    """
    s = str(s).strip()
    if s == "" or s == "None":
        raise ValueError("empty string")
    return float(s)


def _norm_azimuth(a: float) -> float:
    """Wrap azimuth to ``[0, 360)`` degrees."""
    a = a % 360.0
    return a if a >= 0.0 else a + 360.0


def _coerce_lat(v: float) -> float:
    return 90.0 if v > 90.0 else (-90.0 if v < -90.0 else v)


def _coerce_lon(v: float) -> float:
    # normalize to [-180, 180)
    return ((v + 180.0) % 360.0) - 180.0


@dataclass
class JSiteProperty:
    """Container for site properties in a J file.

    Parameters
    ----------
    azimuth : float | None, default=None
        Site X‑axis azimuth in degrees (true north reference).
    latitude : float | None, default=None
        Latitude in decimal degrees.
    longitude : float | None, default=None
        Longitude in decimal degrees.
    elevation : float | None, default=None
        Elevation in metres.
    extra : dict[str, str], optional
        Any unrecognized ``>KEY=VALUE`` pairs preserved verbatim.
    verbose : int | bool, default=0
        Verbosity level used for warnings during parsing.
    strict : bool, default=False
        If ``True``, out‑of‑range values raise errors instead of
        being coerced.
    """

    azimuth: float | None = None
    latitude: float | None = None
    longitude: float | None = None
    elevation: float | None = None
    extra: Dict[str, str] = field(default_factory=dict)
    verbose: int | bool = 0
    strict: bool = False

    @classmethod
    def from_file(
        cls,
        obj: io.TextIOBase | str | Path,
        *,
        encoding: str = ENCODING_DEFAULT,
        strict: bool = False,
        verbose: int | bool = 0,
    ) -> "JSiteProperty":
        """Create from a J file by reading the info block.

        Only the header information lines (``>KEY=VALUE``) are
        consumed. Later blocks (stations/data) are ignored here.
        """
        lines = iter_lines(obj, encoding=encoding)
        return cls.from_lines(lines, strict=strict, verbose=verbose)

    @classmethod
    def from_lines(
        cls,
        lines: Iterable[str],
        *,
        strict: bool = False,
        verbose: int | bool = 0,
    ) -> "JSiteProperty":
        """Create from an iterable of lines.

        Stops at the first non‑info, non‑comment, non‑blank line.
        """
        props: Dict[str, float | None] = {}
        extra: Dict[str, str] = {}

        for key, val in iter_info(lines):
            key_up = key.upper()
            try:
                if key_up == "AZIMUTH":
                    az = _to_float_safe(val)
                    props["azimuth"] = _norm_azimuth(az)
                elif key_up == "LATITUDE":
                    lat = _parse_latitude(val, strict=strict, verbose=verbose)
                    props["latitude"] = lat
                elif key_up == "LONGITUDE":
                    lon = _parse_longitude(val, strict=strict, verbose=verbose)
                    props["longitude"] = lon
                elif key_up == "ELEVATION":
                    el = _maybe_missing_float(val)
                    props["elevation"] = float(el) if el is not None else None
                else:
                    # preserve unknown keys verbatim
                    extra[key_up] = str(val).strip()
            except ValueError as exc:
                msg = f"Failed to parse {key_up}='{val}': {exc}"
                if strict:
                    raise JParseError(msg) from exc
                _vwarn(msg, verbose)

        return cls(
            azimuth=props.get("azimuth"),
            latitude=props.get("latitude"),
            longitude=props.get("longitude"),
            elevation=props.get("elevation"),
            extra=extra,
            verbose=verbose,
            strict=strict,
        )

    @classmethod
    def from_mapping(
        cls,
        mapping: Mapping[str, str | float | int | None],
        *,
        strict: bool = False,
        verbose: int | bool = 0,
    ) -> "JSiteProperty":
        """Create from a mapping of key→value pairs.

        Keys may be any case; values can be strings or numbers.
        """
        lines = (f">{k}={v}" for k, v in mapping.items())
        return cls.from_lines(lines, strict=strict, verbose=verbose)

    def asdict(self) -> Dict[str, float | None]:
        """Return a dictionary view of canonical fields only."""
        return {
            "azimuth": self.azimuth,
            "latitude": self.latitude,
            "longitude": self.longitude,
            "elevation": self.elevation,
        }

    @property
    def location(self) -> Tuple[float, float] | None:
        """Return ``(lat, lon)`` if both are defined, else ``None``."""
        if self.latitude is None or self.longitude is None:
            return None
        return (self.latitude, self.longitude)

    @property
    def azimuth_rad(self) -> float | None:
        """Azimuth in radians if defined."""
        return None if self.azimuth is None else math.radians(self.azimuth)

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        fields = ", ".join(
            f"{k}={v!r}" for k, v in self.asdict().items() if v is not None
        )
        if self.extra:
            fields += f", extra={len(self.extra)} keys"
        return f"JSiteProperty({fields})"


def _maybe_missing_float(value: str) -> float | None:
    try:
        v = _to_float_safe(value)
    except Exception:
        return None
    return None if math.isclose(v, MISSING_FLOAT, rel_tol=0.0, abs_tol=0.0) else v


def _parse_latitude(value: str, *, strict: bool, verbose: int | bool) -> float | None:
    """Parse latitude using GIS utilities.

    - In ``strict`` mode, delegate to :func:`assert_lat_value` which
      raises for out‑of‑range inputs.
    - Otherwise, parse with :func:`convert_position_str2float` and
      coerce into ``[-90, 90]`` with a warning if needed.
    """
    if strict:
        # may return None when input is explicitly missing
        return _assert_lat(value)

    lat = _to_deg(value)
    if lat is None:
        return None
    if not (-90.0 <= lat <= 90.0):
        _vwarn(f"Latitude {lat:.6f} coerced to range", verbose)
        lat = _coerce_lat(lat)
    return lat


def _parse_longitude(value: str, *, strict: bool, verbose: int | bool) -> float | None:
    """Parse longitude using GIS utilities.

    - In ``strict`` mode, delegate to :func:`assert_lon_value` which
      raises for out‑of‑range inputs.
    - Otherwise, parse with :func:`convert_position_str2float` and
      normalize to ``[-180, 180)``; warn if a wrap was applied.
    """
    if strict:
        return _assert_lon(value)

    lon = _to_deg(value)
    if lon is None:
        return None
    fixed = _coerce_lon(lon)
    if abs(fixed - lon) > 1e-12:
        _vwarn(f"Longitude {lon:.6f} normalized to {fixed:.6f}", verbose)
    return fixed


def _vwarn(msg: str, verbose: int | bool) -> None:
    if verbose:
        warnings.warn(msg, RuntimeWarning, stacklevel=2)


__all__ = [
    "JSiteProperty",
]
