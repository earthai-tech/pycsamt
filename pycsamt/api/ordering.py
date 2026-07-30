"""Package-wide site ordering policy.

Configure once and every loader/processor that normalizes through
``ensure_sites`` uses the same policy::

    from pycsamt.api import configure_ordering

    configure_ordering(mode="auto")

Per-call ``order_by=...`` arguments remain authoritative overrides.
"""

from __future__ import annotations

import copy
from collections.abc import Generator
from contextlib import contextmanager
from dataclasses import dataclass, fields
from typing import Any

__all__ = [
    "SiteOrderingConfig",
    "PYCSAMT_ORDERING",
    "configure_ordering",
    "reset_ordering",
]

_ALIASES = {
    "name": "station",
    "natural": "station",
    "lat": "latitude",
    "lon": "longitude",
    "profile": "chainage",
    "spatial": "chainage",
    "none": "input",
    "preserve": "input",
}
_MODES = {"auto", "chainage", "input", "station", "latitude", "longitude"}


@dataclass
class SiteOrderingConfig:
    """Global site-ordering strategy and automatic-line thresholds."""

    mode: str = "auto"
    min_linearity: float = 0.95
    max_cross_track_ratio: float = 0.15
    min_coordinate_fraction: float = 0.60

    def configure(self, **kw: Any) -> SiteOrderingConfig:
        valid = {field.name for field in fields(self)}
        unknown = set(kw) - valid
        if unknown:
            raise AttributeError(
                f"Unknown ordering config key(s) {sorted(unknown)}. "
                f"Valid keys: {sorted(valid)}"
            )
        values = {
            field.name: getattr(self, field.name) for field in fields(self)
        }
        values.update(kw)
        raw_mode = str(values["mode"]).strip().lower()
        mode = _ALIASES.get(raw_mode, raw_mode)
        if mode not in _MODES:
            raise ValueError(
                f"mode must be one of {sorted(_MODES)}, got {values['mode']!r}"
            )
        for key in (
            "min_linearity",
            "max_cross_track_ratio",
            "min_coordinate_fraction",
        ):
            value = float(values[key])
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"{key} must be between 0 and 1, got {value}")
            values[key] = value
        values["mode"] = mode
        for key, value in values.items():
            setattr(self, key, value)
        return self

    @contextmanager
    def context(self, **kw: Any) -> Generator[SiteOrderingConfig, None, None]:
        """Temporarily override ordering settings, then restore them."""

        snapshot = self.clone()
        try:
            self.configure(**kw)
            yield self
        finally:
            for field in fields(self):
                setattr(self, field.name, getattr(snapshot, field.name))

    def reset(self) -> None:
        """Restore package defaults."""

        defaults = SiteOrderingConfig()
        for field in fields(self):
            setattr(self, field.name, getattr(defaults, field.name))

    def clone(self) -> SiteOrderingConfig:
        """Return an independent copy."""

        return copy.deepcopy(self)


PYCSAMT_ORDERING = SiteOrderingConfig()


def configure_ordering(**kw: Any) -> SiteOrderingConfig:
    """Configure and return the global ordering singleton."""

    return PYCSAMT_ORDERING.configure(**kw)


def reset_ordering() -> None:
    """Reset the global ordering policy to ``mode='auto'``."""

    PYCSAMT_ORDERING.reset()
