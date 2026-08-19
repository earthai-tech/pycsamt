# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Reusable electromagnetic channel and site-layout metadata."""

from __future__ import annotations

from dataclasses import dataclass, field
from math import isfinite
from typing import Any

from ..api.property import PyCSAMTObject

__all__ = ["ChannelMeta", "SiteLayout"]

_FIELD_ALIASES = {
    "e": "electric",
    "electric": "electric",
    "h": "magnetic",
    "b": "magnetic",
    "magnetic": "magnetic",
    "other": "other",
}


@dataclass(repr=False)
class ChannelMeta(PyCSAMTObject):
    """Describe one electric, magnetic, or auxiliary EM channel.

    Coordinates describe the physical site layout and are independent of the
    orientation applied to transfer-function matrices.
    """

    name: str
    field_type: str
    orientation: float | None = None
    tilt: float | None = None
    x: float | None = None
    y: float | None = None
    z: float | None = None
    x2: float | None = None
    y2: float | None = None
    z2: float | None = None
    units: str | None = None
    reference: str | None = None
    sensor_id: str | None = None
    extra: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        self.name = str(self.name).strip()
        if not self.name:
            raise ValueError("channel name must be non-empty")

        raw_type = str(self.field_type).strip().lower()
        if raw_type not in _FIELD_ALIASES:
            raise ValueError(
                "field_type must describe electric, magnetic, or other data"
            )
        self.field_type = _FIELD_ALIASES[raw_type]

        for attr in (
            "orientation",
            "tilt",
            "x",
            "y",
            "z",
            "x2",
            "y2",
            "z2",
        ):
            value = getattr(self, attr)
            if value is None:
                continue
            numeric = float(value)
            if not isfinite(numeric):
                raise ValueError(f"{attr} must be finite")
            setattr(self, attr, numeric)

        for attr in ("units", "reference", "sensor_id"):
            value = getattr(self, attr)
            if value is not None:
                text = str(value).strip()
                setattr(self, attr, text or None)
        self.extra = dict(self.extra or {})

    @property
    def is_electric(self) -> bool:
        return self.field_type == "electric"

    @property
    def is_magnetic(self) -> bool:
        return self.field_type == "magnetic"

    @property
    def has_second_endpoint(self) -> bool:
        return any(value is not None for value in (self.x2, self.y2, self.z2))


@dataclass(repr=False)
class SiteLayout(PyCSAMTObject):
    """Original physical input/output channel geometry at a site."""

    input_channels: list[ChannelMeta] = field(default_factory=list)
    output_channels: list[ChannelMeta] = field(default_factory=list)
    input_units: str | None = None
    output_units: str | None = None
    input_reference: str | None = None
    output_reference: str | None = None
    extra: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        self.input_channels = self._validate_group(
            self.input_channels, "input"
        )
        self.output_channels = self._validate_group(
            self.output_channels, "output"
        )
        for attr in (
            "input_units",
            "output_units",
            "input_reference",
            "output_reference",
        ):
            value = getattr(self, attr)
            if value is not None:
                text = str(value).strip()
                setattr(self, attr, text or None)
        self.extra = dict(self.extra or {})

    @staticmethod
    def _validate_group(
        channels: list[ChannelMeta] | tuple[ChannelMeta, ...],
        role: str,
    ) -> list[ChannelMeta]:
        out = list(channels or [])
        for channel in out:
            if not isinstance(channel, ChannelMeta):
                raise TypeError(
                    f"{role}_channels must contain ChannelMeta instances"
                )
        names = [channel.name.lower() for channel in out]
        if len(names) != len(set(names)):
            raise ValueError(f"duplicate channel name in {role}_channels")
        return out

    @property
    def input_names(self) -> tuple[str, ...]:
        return tuple(channel.name for channel in self.input_channels)

    @property
    def output_names(self) -> tuple[str, ...]:
        return tuple(channel.name for channel in self.output_channels)

    @property
    def n_input(self) -> int:
        return len(self.input_channels)

    @property
    def n_output(self) -> int:
        return len(self.output_channels)

    def get_channel(
        self,
        name: str,
        *,
        role: str | None = None,
    ) -> ChannelMeta | None:
        """Return a channel by name, optionally restricted by role."""
        key = str(name).strip().lower()
        groups: list[list[ChannelMeta]]
        if role is None:
            groups = [self.input_channels, self.output_channels]
        else:
            normalized = str(role).strip().lower()
            if normalized not in {"input", "output"}:
                raise ValueError("role must be 'input', 'output', or None")
            groups = [
                self.input_channels
                if normalized == "input"
                else self.output_channels
            ]
        for group in groups:
            for channel in group:
                if channel.name.lower() == key:
                    return channel
        return None
