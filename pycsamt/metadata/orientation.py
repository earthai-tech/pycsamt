# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Transfer-function orientation metadata independent of file format."""

from __future__ import annotations

from dataclasses import dataclass, field
from math import isfinite
from typing import Any

from ..api.property import PyCSAMTObject

__all__ = ["OrientationMeta"]

_ALLOWED_MODES = frozenset({"orthogonal", "sitelayout"})
_MODE_ALIASES = {
    "site_layout": "sitelayout",
    "site-layout": "sitelayout",
    "layout": "sitelayout",
}


@dataclass(repr=False)
class OrientationMeta(PyCSAMTObject):
    """Describe the coordinate orientation of transfer-function data.

    ``OrientationMeta`` describes the orientation of the *data*, not the
    physical sensor geometry.  Physical channel geometry is retained by
    :class:`~pycsamt.metadata.channels.SiteLayout` so rotations remain
    reversible where sufficient information is available.

    Parameters
    ----------
    mode : {"orthogonal", "sitelayout"}, optional
        ``"orthogonal"`` means the TF matrix is expressed in an orthogonal
        coordinate frame. ``"sitelayout"`` means the TF follows the original
        site channel layout.
    angle_to_geographic_north : float, optional
        Clockwise angle in degrees of the orthogonal x-axis relative to
        geographic north. It is intentionally not fabricated when unknown.
    rotation_info : str, optional
        Human-readable record of known rotation history or ambiguity.
    extra : dict
        Additional orientation metadata.
    """

    mode: str | None = None
    angle_to_geographic_north: float | None = None
    rotation_info: str | None = None
    extra: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        if self.mode is not None:
            mode = str(self.mode).strip().lower()
            mode = _MODE_ALIASES.get(mode, mode)
            if mode not in _ALLOWED_MODES:
                raise ValueError(
                    "orientation mode must be one of "
                    f"{sorted(_ALLOWED_MODES)}"
                )
            self.mode = mode

        if self.angle_to_geographic_north is not None:
            angle = float(self.angle_to_geographic_north)
            if not isfinite(angle):
                raise ValueError("orientation angle must be finite")
            if self.mode == "sitelayout":
                raise ValueError(
                    "sitelayout orientation is defined by channel geometry; "
                    "angle_to_geographic_north must be None"
                )
            self.angle_to_geographic_north = angle

        if self.rotation_info is not None:
            value = str(self.rotation_info).strip()
            self.rotation_info = value or None
        self.extra = dict(self.extra or {})

    @property
    def is_orthogonal(self) -> bool:
        """Return whether the data are in an orthogonal coordinate frame."""
        return self.mode == "orthogonal"

    @property
    def follows_site_layout(self) -> bool:
        """Return whether data orientation follows the sensor layout."""
        return self.mode == "sitelayout"
