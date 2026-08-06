"""
Site‑level properties parsed from the J‑format information block.

This module defines a lightweight container for the metadata
exposed in lines beginning with ``>KEY=VALUE`` in A.G. Jones J files.
It focuses on robust parsing and normalization (e.g., azimuth
wrapping, latitude/longitude hemispheres, and missing sentinels).
"""

from __future__ import annotations

import io
import math
import warnings
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path

from ..exceptions import JParseError
from ..gis.utils import (
    assert_lat_value as _assert_lat,
)
from ..gis.utils import (
    assert_lon_value as _assert_lon,
)
from ..gis.utils import (
    convert_position_str2float as _to_deg,
)
from .config import ENCODING_DEFAULT, MISSING_FLOAT
from .utils import iter_info, iter_lines

__all__ = ["JSiteProperty"]


@dataclass
class JSiteProperty:
    r"""
    Container for site properties parsed from J-format info.

    This class stores site-level metadata such as latitude,
    longitude, azimuth, and elevation that are declared in the
    J-format *information block* as ``>KEY=VALUE`` lines.

    Parameters
    ----------
    azimuth : float or None, optional
        Site X-axis azimuth in degrees, referenced to true
        north. Values are wrapped to ``[0, 360)`` when parsed.
    latitude : float or None, optional
        Latitude in decimal degrees. Inputs may include DMS or
        hemisphere suffixes and are normalized to ``[-90, 90]``.
    longitude : float or None, optional
        Longitude in decimal degrees. Inputs may include DMS or
        hemisphere suffixes and are normalized to ``[-180, 180)``.
    elevation : float or None, optional
        Elevation in metres. Missing sentinels are mapped to
        ``None``.
    extra : dict of (str -> str), optional
        Unrecognized ``>KEY=VALUE`` pairs preserved verbatim.
    verbose : int or bool, default 0
        Verbosity for non-fatal parsing warnings.
    strict : bool, default False
        If ``True``, out-of-range coordinates raise errors
        instead of being coerced to valid bounds.

    Attributes
    ----------
    location : tuple of (float, float) or None
        Convenience property returning ``(lat, lon)`` when both
        are available, else ``None``.
    azimuth_rad : float or None
        Azimuth converted to radians when available.

    Notes
    -----
    Parsing leverages robust coordinate helpers that accept DMS,
    hemisphere letters, and free whitespace. Longitudes are
    normalized into the half-open interval ``[-180, 180)``.

    Examples
    --------
    >>> lines = [\">LATITUDE = 41.9782\", \">LONGITUDE = 140.8958\"]\n
    >>> prop = JSiteProperty.from_lines(lines)\n
    >>> prop.location\n
    (41.9782, 140.8958)

    See Also
    --------
    Info
        High-level container of the information block which can
        be used to build a :class:`JSiteProperty`.

    References
    ----------
    .. [JSiteProperty-1] A. G. Jones (1994). Magnetotelluric data file J-format,
       version 2.0.
    """

    azimuth: float | None = None
    latitude: float | None = None
    longitude: float | None = None
    elevation: float | None = None
    extra: dict[str, str] = field(default_factory=dict)
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
    ) -> JSiteProperty:
        r"""
        Build a :class:`JSiteProperty` from a J file path or
        file-like object.

        Only the initial *information block* is consumed. Later
        station or data sections are ignored by this method.

        Parameters
        ----------
        obj : str or pathlib.Path or TextIO
            Path to a J file or an open text stream positioned at
            the beginning of the file.
        encoding : str, default 'utf-8'
            Text encoding used when reading from a file path.
        strict : bool, default False
            If ``True``, invalid coordinates raise. If ``False``,
            out-of-range values are coerced with a warning.
        verbose : int or bool, default 0
            Verbosity for non-fatal parsing warnings.

        Returns
        -------
        JSiteProperty
            Parsed site properties.

        Notes
        -----
        This method stops at the first line that is not a comment,
        not blank, and not a ``>KEY=VALUE`` record.

        Examples
        --------
        >>> prop = JSiteProperty.from_file(\"data/j/kb0-s001.txt\")\n
        >>> prop.azimuth is None or 0.0 <= prop.azimuth < 360.0\n
        True

        See Also
        --------
        from_lines
            Parse from an in-memory sequence of lines.
        from_mapping
            Quick constructor from a dict of keys and values.

        References
        ----------
        .. [JSiteProperty-from-file-1] A. G. Jones (1994). Magnetotelluric data file
           J-format, version 2.0.
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
    ) -> JSiteProperty:
        r"""
        Build a :class:`JSiteProperty` from an iterable of lines.

        The parser reads sequentially and collects ``>KEY=VALUE``
        pairs as long as lines are comments, blanks, or info
        records. It stops at the first non-conforming line.

        Parameters
        ----------
        lines : iterable of str
            Lines that include the information block content.
        strict : bool, default False
            If ``True``, invalid coordinates raise. If ``False``,
            out-of-range values are coerced with a warning.
        verbose : int or bool, default 0
            Verbosity for non-fatal parsing warnings.

        Returns
        -------
        JSiteProperty
            Parsed site properties.

        Notes
        -----
        Latitude and longitude accept flexible text including DMS
        and hemisphere tokens. Longitudes are normalized into
        ``[-180, 180)``.

        Examples
        --------
        >>> from pycsamt.
        >>> lines = [\n
        ...   \">AZIMUTH = -30.0\", \">LATITUDE = 41.9782\",\n
        ...   \">LONGITUDE = 140.8958\", \">ELEVATION = 0.0\",\n
        ... ]\n
        >>> prop = JSiteProperty.from_lines(lines)\n
        >>> prop.location[0] == 41.9782\n
        True

        See Also
        --------
        from_file
            Parse from a file path or open stream.
        from_mapping
            Quick constructor from a dict of keys and values.

        References
        ----------
        .. [JSiteProperty-from-lines-1] A. G. Jones (1994). Magnetotelluric data file
           J-format, version 2.0.
        """

        props: dict[str, float | None] = {}
        extra: dict[str, str] = {}
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
    ) -> JSiteProperty:
        r"""
        Build a :class:`JSiteProperty` from a dict of key-value
        pairs.

        Keys can be any case and will be normalized to the upper
        case tokens used by J files. Values can be strings or
        numbers. Unknown keys are preserved in ``extra``.

        Parameters
        ----------
        mapping : Mapping[str, Any]
            Pairs like ``{\"latitude\": 41.9, \"longitude\": 140.8}``.
        strict : bool, default False
            If ``True``, invalid coordinates raise. If ``False``,
            out-of-range values are coerced with a warning.
        verbose : int or bool, default 0
            Verbosity for non-fatal parsing warnings.

        Returns
        -------
        JSiteProperty
            Parsed site properties.

        Notes
        -----
        This is a convenience front-end that converts the mapping
        to synthetic ``>KEY=VALUE`` lines and forwards the work to
        :meth:`from_lines`.

        Examples
        --------
        >>> prop = JSiteProperty.from_mapping({\n
        ...   'azimuth': 370, 'latitude': '41:58N',\n
        ...   'longitude': '140:54E', 'elevation': 0,\n
        ... })\n
        >>> 0.0 <= (prop.azimuth or 0.0) < 360.0\n
        True

        See Also
        --------
        from_lines
            Parse from an in-memory sequence of lines.
        from_file
            Parse from a file path or open stream.

        References
        ----------
        .. [JSiteProperty-from-mapping-1] A. G. Jones (1994). Magnetotelluric data file
           J-format, version 2.0.
        """
        lines = (f">{k}={v}" for k, v in mapping.items())
        return cls.from_lines(lines, strict=strict, verbose=verbose)

    def asdict(self) -> dict[str, float | None]:
        return {
            "azimuth": self.azimuth,
            "latitude": self.latitude,
            "longitude": self.longitude,
            "elevation": self.elevation,
        }

    @property
    def location(self) -> tuple[float, float] | None:
        if self.latitude is None or self.longitude is None:
            return None
        return (self.latitude, self.longitude)

    @property
    def azimuth_rad(self) -> float | None:
        if self.azimuth is None:
            return None
        return math.radians(self.azimuth)

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        parts = []
        for k, v in self.asdict().items():
            if v is not None:
                parts.append(f"{k}={v!r}")
        s = ", ".join(parts)
        if self.extra:
            s = (
                f"{s}, extra={len(self.extra)} keys"
                if s
                else (f"extra={len(self.extra)} keys")
            )
        return f"JSiteProperty({s})"


def _maybe_missing_float(value: str) -> float | None:
    try:
        v = _to_float_safe(value)
    except Exception:
        return None
    if math.isclose(v, MISSING_FLOAT, rel_tol=0.0, abs_tol=0.0):
        return None
    return v


def _parse_latitude(
    value: str, *, strict: bool, verbose: int | bool
) -> float | None:
    """Parse latitude using GIS utilities.

    - In ``strict`` mode, delegate to :func:`assert_lat_value` which
      raises for out‑of‑range inputs.
    - Otherwise, parse with :func:`convert_position_str2float` and
      coerce into ``[-90, 90]`` with a warning if needed.
    """

    if strict:
        return _assert_lat(value)
    lat = _to_deg(value)
    if lat is None:
        # plain decimal fallback
        try:
            lat = float(str(value).strip())
        except Exception:
            return None
    if not (-90.0 <= lat <= 90.0):
        _vwarn(f"Latitude {lat:.6f} coerced to range", verbose)
        lat = _coerce_lat(lat)
    return lat


def _parse_longitude(
    value: str, *, strict: bool, verbose: int | bool
) -> float | None:
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
        # plain decimal fallback
        try:
            lon = float(str(value).strip())
        except Exception:
            return None

    fixed = _coerce_lon(lon)
    if abs(fixed - lon) > 1e-12:
        _vwarn(f"Longitude {lon:.6f} normalized to {fixed:.6f}", verbose)
    return fixed


def _vwarn(msg: str, verbose: int | bool) -> None:
    if verbose:
        warnings.warn(msg, RuntimeWarning, stacklevel=2)


def _to_float_safe(s: str) -> float:
    """Best‑effort float parser for simple
    numeric fields (e.g., azimuth).

    Unlike :func:`_to_deg`, this does *not*
    parse DMS — it's meant for scalar numeric
    tokens where non‑numeric input should
    surface quickly.
    """
    s = str(s).strip()
    if s == "" or s == "None":
        raise ValueError("empty string")
    return float(s)


def _norm_azimuth(a: float) -> float:
    a = a % 360.0
    return a if a >= 0.0 else a + 360.0


def _coerce_lat(v: float) -> float:
    return 90.0 if v > 90.0 else (-90.0 if v < -90.0 else v)


def _coerce_lon(v: float) -> float:
    return ((v + 180.0) % 360.0) - 180.0
