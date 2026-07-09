"""Shared validation and coercion helpers for the IoT subpackage.

These helpers centralise the small, defensive value-normalisation logic
used across the IoT modules so that packets, configs, and telemetry
schemas validate field values consistently. They intentionally have no
dependency on the rest of pyCSAMT beyond the standard library, numpy,
and the local enums.
"""

from __future__ import annotations

import math
import time
from enum import Enum
from typing import Any, Iterable, List, Optional

__all__ = [
    "as_nonempty_str",
    "as_optional_str",
    "as_bool",
    "as_finite_float",
    "as_optional_finite_float",
    "as_positive",
    "as_nonnegative",
    "as_optional_positive",
    "as_probability",
    "as_timestamp",
    "as_qos",
    "normalise_enum",
    "as_latitude",
    "as_longitude",
    "as_elevation",
    "utc_now",
]


_TRUE_STRINGS = frozenset({"true", "1", "yes", "y", "on", "t"})
_FALSE_STRINGS = frozenset({"false", "0", "no", "n", "off", "f"})


def as_nonempty_str(value: Any, name: str) -> str:
    """Return *value* as a stripped, non-empty string or raise."""
    text = str(value).strip()
    if not text:
        raise ValueError(f"{name} cannot be empty.")
    return text


def as_optional_str(value: Any, name: str) -> Optional[str]:
    """Return a stripped string, or ``None`` when *value* is empty/None."""
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def as_bool(value: Any, *, default: Optional[bool] = None) -> bool:
    """Interpret *value* as a boolean, parsing strings explicitly.

    Unlike :func:`bool`, string values such as ``"false"``, ``"0"``,
    ``"no"``, and ``"off"`` are treated as ``False``. ``NaN`` is treated
    as ``False``. Unrecognised strings raise :class:`ValueError` so that
    malformed telemetry is surfaced rather than silently accepted.
    """
    if isinstance(value, bool):
        return value
    if value is None:
        return bool(default) if default is not None else False
    if isinstance(value, (int, float)):
        if isinstance(value, float) and math.isnan(value):
            return False
        return bool(value)
    text = str(value).strip().lower()
    if text in _TRUE_STRINGS:
        return True
    if text in _FALSE_STRINGS or text in {"", "none", "null"}:
        return False
    raise ValueError(f"Cannot interpret {value!r} as a boolean.")


def as_finite_float(value: Any, name: str) -> float:
    """Return *value* as a finite float or raise on ``NaN``/``inf``."""
    out = float(value)
    if not math.isfinite(out):
        raise ValueError(f"{name} must be a finite number.")
    return out


def as_optional_finite_float(value: Any, name: str) -> Optional[float]:
    """Return a finite float, or ``None`` when *value* is ``None``."""
    if value is None:
        return None
    return as_finite_float(value, name)


def as_positive(value: Any, name: str) -> float:
    """Return *value* as a strictly positive finite float or raise."""
    out = as_finite_float(value, name)
    if out <= 0:
        raise ValueError(f"{name} must be positive.")
    return out


def as_nonnegative(value: Any, name: str) -> float:
    """Return *value* as a non-negative finite float or raise."""
    out = as_finite_float(value, name)
    if out < 0:
        raise ValueError(f"{name} must be >= 0.")
    return out


def as_optional_positive(value: Any, name: str) -> Optional[float]:
    """Return a strictly positive finite float, or ``None``."""
    if value is None:
        return None
    return as_positive(value, name)


def as_probability(value: Any, name: str) -> float:
    """Return *value* as a float in the closed interval ``[0, 1]``."""
    out = as_finite_float(value, name)
    if not 0.0 <= out <= 1.0:
        raise ValueError(f"{name} must be between 0 and 1.")
    return out


def as_timestamp(value: Any, name: str = "timestamp") -> float:
    """Return a finite, non-negative epoch/relative timestamp.

    Rejects ``NaN``, infinities, and negative values, which are common
    signs of an un-synchronised or mis-parsed field clock. Zero is
    permitted so relative session clocks starting at ``t=0`` remain
    valid.
    """
    out = float(value)
    if not math.isfinite(out):
        raise ValueError(f"{name} must be finite (got {value!r}).")
    if out < 0:
        raise ValueError(f"{name} must be non-negative (got {out}).")
    return out


def as_qos(value: Any, name: str = "qos") -> int:
    """Return an MQTT-style quality-of-service level in ``{0, 1, 2}``."""
    out = int(value)
    if out not in {0, 1, 2}:
        raise ValueError(f"{name} must be 0, 1, or 2.")
    return out


def normalise_enum(value: Any, enum_cls: type[Enum], name: str) -> Enum:
    """Coerce *value* to a member of *enum_cls* by value or name."""
    if isinstance(value, enum_cls):
        return value
    text = str(value).strip().lower()
    for item in enum_cls:
        if text in {str(item.value).lower(), item.name.lower()}:
            return item
    allowed = ", ".join(str(item.value) for item in enum_cls)
    raise ValueError(f"{name} must be one of: {allowed}.")


def as_latitude(value: Any, name: str = "lat") -> Optional[float]:
    """Return a latitude in ``[-90, 90]`` degrees, or ``None``."""
    if value is None:
        return None
    out = as_finite_float(value, name)
    if not -90.0 <= out <= 90.0:
        raise ValueError(f"{name} must be within [-90, 90] degrees.")
    return out


def as_longitude(value: Any, name: str = "lon") -> Optional[float]:
    """Return a longitude in ``[-180, 360]`` degrees, or ``None``.

    Both ``[-180, 180]`` and ``[0, 360]`` conventions are accepted so
    that raw GPS strings and normalised longitudes both pass.
    """
    if value is None:
        return None
    out = as_finite_float(value, name)
    if not -180.0 <= out <= 360.0:
        raise ValueError(f"{name} must be within [-180, 360] degrees.")
    return out


def as_elevation(value: Any, name: str = "elevation") -> Optional[float]:
    """Return a finite elevation in metres, or ``None``."""
    if value is None:
        return None
    return as_finite_float(value, name)


def as_channel_list(values: Any) -> List[str]:
    """Return a list of lower-case, non-empty channel labels."""
    if values is None:
        return []
    if isinstance(values, str):
        values = [values]
    out: List[str] = []
    for value in list(values):
        text = str(value).strip().lower()
        if not text:
            raise ValueError("channel labels cannot be empty.")
        out.append(text)
    return out


def unique_preserving(values: Iterable[Any]) -> List[Any]:
    """Return *values* with duplicates removed, order preserved."""
    seen: set[Any] = set()
    out: List[Any] = []
    for value in values:
        if value not in seen:
            seen.add(value)
            out.append(value)
    return out


def utc_now() -> float:
    """Return the current UTC time as a POSIX timestamp (seconds)."""
    return float(time.time())
