"""Global section-view controls for pyCSAMT figures.

This module centralizes layout choices used by 2-D section-like plots:
pseudo-sections, inversion sections, model slices, and compact dashboard
panels.  It complements :mod:`pycsamt.api.style`,
:mod:`pycsamt.api.station`, :mod:`pycsamt.api.control`, and
:mod:`pycsamt.api.plot`.

Examples
--------
Apply a publication section preset temporarily::

    from pycsamt.api.section import PYCSAMT_SECTION
    from pycsamt.tdem import plot_temavg_section

    with PYCSAMT_SECTION.context("publication"):
        ax = plot_temavg_section(avg, section="publication")

Use the data-aware dynamic layout for many stations::

    ax = plot_temavg_section(avg, section="dynamic")

Customize one preset globally::

    from pycsamt.api import configure_section
    configure_section(
        dynamic__figure__max_width=16.0,
        dynamic__colorbar__max_ticks=5,
    )
"""

from __future__ import annotations

import copy
from collections.abc import Generator, Sequence
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Any

from pycsamt.api.plot import add_colorbar
from pycsamt.api.station import (
    PYCSAMT_STATION_RENDERING,
    StationAxisStyle,
)

_PRESETS = {
    "compact",
    "dashboard",
    "dynamic",
    "inversion",
    "publication",
    "pseudosection",
}
_Y_DIRECTIONS = {"down", "up", "normal", "none"}
_ASPECTS = {"auto", "equal"}


@dataclass
class SectionFigureStyle:
    """Figure geometry for section-like plots."""

    figsize: tuple[float, float] | str = (9.5, 4.4)
    tight: bool = True
    constrained: bool = False
    title: bool = True
    min_width: float = 7.0
    max_width: float = 14.0
    min_height: float = 3.2
    max_height: float = 7.5
    station_width: float = 0.16
    y_sample_height: float = 0.030
    label_width: float = 0.055
    colorbar_width: float = 0.55

    def resolve(
        self,
        *,
        n_stations: int | None = None,
        n_y: int | None = None,
        max_label_len: int | None = None,
        colorbar: bool = True,
    ) -> tuple[float, float]:
        """Return a concrete figure size in inches."""
        if self.figsize != "dynamic":
            width, height = self.figsize
            return float(width), float(height)
        n_stations = max(1, int(n_stations or 1))
        n_y = max(1, int(n_y or 1))
        max_label_len = max(2, int(max_label_len or 4))
        label_room = max_label_len * float(self.label_width)
        width = 2.4 + n_stations * float(self.station_width) + label_room
        if colorbar:
            width += float(self.colorbar_width)
        height = 2.2 + n_y * float(self.y_sample_height)
        width = min(max(width, self.min_width), self.max_width)
        height = min(max(height, self.min_height), self.max_height)
        return float(width), float(height)


@dataclass
class SectionAxisStyle:
    """Axis behavior for section-like plots."""

    xlabel: str = "Station"
    ylabel: str | None = None
    y_direction: str = "down"
    aspect: str | float = "auto"
    grid: bool = False
    station_side: str = "top"
    title: bool = True
    xlabel_pad: float | None = None
    ylabel_pad: float | None = None

    def apply(
        self,
        ax: Any,
        *,
        xlabel: str | None = None,
        ylabel: str | None = None,
        title: str | None = None,
    ) -> None:
        """Apply axis labels, direction, aspect, and grid state."""
        y_direction = self.y_direction.lower()
        if y_direction not in _Y_DIRECTIONS:
            msg = (
                "section y_direction must be one of "
                f"{sorted(_Y_DIRECTIONS)}."
            )
            raise ValueError(msg)
        if isinstance(self.aspect, str) and self.aspect not in _ASPECTS:
            msg = "section aspect must be 'auto', 'equal', or a number."
            raise ValueError(msg)

        xlab = self.xlabel if xlabel is None else xlabel
        ylab = self.ylabel if ylabel is None else ylabel
        if xlab:
            if self.xlabel_pad is None:
                ax.set_xlabel(xlab)
            else:
                ax.set_xlabel(xlab, labelpad=self.xlabel_pad)
        if ylab:
            if self.ylabel_pad is None:
                ax.set_ylabel(ylab)
            else:
                ax.set_ylabel(ylab, labelpad=self.ylabel_pad)
        if title and self.title:
            ax.set_title(title)
        if y_direction == "down" and not ax.yaxis_inverted():
            ax.invert_yaxis()
        elif y_direction in {"up", "normal"} and ax.yaxis_inverted():
            ax.invert_yaxis()
        if self.aspect is not None:
            ax.set_aspect(self.aspect)
        ax.grid(bool(self.grid))


@dataclass
class SectionColorbarStyle:
    """Colorbar geometry and tick-density controls."""

    show: bool = True
    side: str = "right"
    size: str = "3.5%"
    pad: float = 0.04
    max_ticks: int | None = 6
    tick_format: str | None = None

    def add(
        self,
        mappable: Any,
        ax: Any,
        *,
        label: str | None = None,
        **kw: Any,
    ) -> Any | None:
        """Attach a section colorbar to ``ax``."""
        if not self.show:
            return None
        options = {
            "side": self.side,
            "size": self.size,
            "pad": self.pad,
            "max_ticks": self.max_ticks,
            "tick_format": self.tick_format,
        }
        options.update(kw)
        return add_colorbar(mappable, ax, label=label, **options)


@dataclass
class SectionStyle:
    """Complete visual-control bundle for one section preset."""

    figure: SectionFigureStyle = field(default_factory=SectionFigureStyle)
    axis: SectionAxisStyle = field(default_factory=SectionAxisStyle)
    colorbar: SectionColorbarStyle = field(
        default_factory=SectionColorbarStyle,
    )
    station_preset: str = "pseudosection"

    def copy(self, **kw: Any) -> "SectionStyle":
        """Return a modified deep copy of this section style."""
        new = copy.deepcopy(self)
        for key, value in kw.items():
            setattr(new, key, value)
        return new

    def figsize_for(
        self,
        *,
        n_stations: int | None = None,
        n_y: int | None = None,
        labels: Sequence[Any] | None = None,
        colorbar: bool = True,
    ) -> tuple[float, float]:
        """Return a concrete figure size for this section."""
        max_label_len = None
        if labels is not None:
            max_label_len = max((len(str(label)) for label in labels), default=4)
        return self.figure.resolve(
            n_stations=n_stations,
            n_y=n_y,
            max_label_len=max_label_len,
            colorbar=colorbar,
        )

    def apply_axis(self, ax: Any, **kw: Any) -> None:
        """Apply this style's axis policy to ``ax``."""
        self.axis.apply(ax, **kw)

    def add_colorbar(
        self,
        mappable: Any,
        ax: Any,
        *,
        label: str | None = None,
        **kw: Any,
    ) -> Any | None:
        """Attach a colorbar using this style's colorbar policy."""
        return self.colorbar.add(mappable, ax, label=label, **kw)

    def apply_stations(
        self,
        ax: Any,
        positions: Sequence[float],
        labels: Sequence[Any] | None = None,
        *,
        station_style: StationAxisStyle | None = None,
        xlim: tuple[float, float] | None = None,
    ) -> Any:
        """Apply station rendering for this section style."""
        style = copy.deepcopy(
            station_style
            or PYCSAMT_STATION_RENDERING.style_for(self.station_preset),
        )
        style.side = self.axis.station_side
        style.xlabel = self.axis.xlabel
        return style.apply(ax, positions, labels, xlim=xlim)


class PyCSAMTSection:
    """Package-wide section-view control container."""

    def __init__(self) -> None:
        self._init_defaults()

    def _init_defaults(self) -> None:
        self.pseudosection = SectionStyle(
            figure=SectionFigureStyle(figsize=(9.5, 4.4)),
            axis=SectionAxisStyle(
                xlabel="Station",
                y_direction="down",
                aspect="auto",
                grid=False,
                station_side="top",
            ),
            colorbar=SectionColorbarStyle(
                size="3.5%",
                pad=0.04,
                max_ticks=6,
            ),
            station_preset="pseudosection",
        )
        self.inversion = SectionStyle(
            figure=SectionFigureStyle(figsize=(10.5, 5.2)),
            axis=SectionAxisStyle(
                xlabel="Station",
                ylabel="Depth (m)",
                y_direction="down",
                aspect="auto",
                grid=False,
                station_side="top",
            ),
            colorbar=SectionColorbarStyle(
                size="3.2%",
                pad=0.035,
                max_ticks=6,
            ),
            station_preset="inversion",
        )
        self.publication = SectionStyle(
            figure=SectionFigureStyle(
                figsize=(8.2, 3.6),
                tight=True,
                title=False,
            ),
            axis=SectionAxisStyle(
                xlabel="Station",
                y_direction="down",
                aspect="auto",
                grid=False,
                station_side="top",
                title=False,
            ),
            colorbar=SectionColorbarStyle(
                size="2.8%",
                pad=0.035,
                max_ticks=5,
            ),
            station_preset="pseudosection",
        )
        self.compact = SectionStyle(
            figure=SectionFigureStyle(
                figsize=(7.0, 3.1),
                tight=True,
                title=False,
            ),
            axis=SectionAxisStyle(
                xlabel="Station",
                y_direction="down",
                aspect="auto",
                grid=False,
                station_side="top",
                title=False,
            ),
            colorbar=SectionColorbarStyle(
                size="2.8%",
                pad=0.03,
                max_ticks=4,
            ),
            station_preset="pseudosection",
        )
        self.dashboard = SectionStyle(
            figure=SectionFigureStyle(
                figsize=(6.0, 3.0),
                tight=False,
                title=True,
            ),
            axis=SectionAxisStyle(
                xlabel="Station",
                y_direction="down",
                aspect="auto",
                grid=False,
                station_side="top",
            ),
            colorbar=SectionColorbarStyle(
                size="3.0%",
                pad=0.03,
                max_ticks=4,
            ),
            station_preset="pseudosection",
        )
        self.dynamic = SectionStyle(
            figure=SectionFigureStyle(
                figsize="dynamic",
                min_width=7.5,
                max_width=15.0,
                min_height=3.4,
                max_height=7.2,
                station_width=0.15,
                y_sample_height=0.032,
                label_width=0.060,
                colorbar_width=0.55,
            ),
            axis=SectionAxisStyle(
                xlabel="Station",
                y_direction="down",
                aspect="auto",
                grid=False,
                station_side="top",
            ),
            colorbar=SectionColorbarStyle(
                size="3.2%",
                pad=0.035,
                max_ticks=6,
            ),
            station_preset="pseudosection",
        )

    def style_for(self, preset: str = "pseudosection") -> SectionStyle:
        """Return the section style associated with ``preset``."""
        preset = str(preset).lower()
        if preset not in _PRESETS:
            msg = f"section preset must be one of {sorted(_PRESETS)}."
            raise ValueError(msg)
        return getattr(self, preset)

    def configure(self, **kw: Any) -> None:
        """Configure section styles using ``section__attribute`` paths."""
        for path, value in kw.items():
            parts = path.split("__")
            obj = self
            for part in parts[:-1]:
                obj = getattr(obj, part)
            setattr(obj, parts[-1], value)

    @contextmanager
    def context(
        self,
        preset: str | None = None,
        **kw: Any,
    ) -> Generator["PyCSAMTSection", None, None]:
        """Temporarily override section styles, then restore them."""
        snapshot = self._snapshot()
        try:
            if preset:
                self.use_preset(preset)
            if kw:
                self.configure(**kw)
            yield self
        finally:
            self._restore(snapshot)

    def use_preset(self, preset: str) -> None:
        """Copy one preset into the pseudosection slot."""
        self.pseudosection = copy.deepcopy(self.style_for(preset))

    def reset(self) -> None:
        """Reset section controls to package defaults."""
        self._init_defaults()

    def summary(self) -> str:
        """Return a human-readable summary."""
        lines = ["PyCSAMTSection"]
        for name in sorted(_PRESETS):
            style = self.style_for(name)
            lines.append(
                "  "
                f"{name}: figsize={style.figure.figsize!r}, "
                f"y={style.axis.y_direction!r}, "
                f"stations={style.station_preset!r}, "
                f"cbar_ticks<={style.colorbar.max_ticks}",
            )
        return "\n".join(lines)

    def __repr__(self) -> str:  # noqa: D105
        return self.summary()

    def _snapshot(self) -> dict[str, SectionStyle]:
        return {
            name: copy.deepcopy(self.style_for(name))
            for name in _PRESETS
        }

    def _restore(self, snapshot: dict[str, SectionStyle]) -> None:
        for name, value in snapshot.items():
            setattr(self, name, value)


PYCSAMT_SECTION = PyCSAMTSection()


def configure_section(**kw: Any) -> None:
    """Configure :data:`PYCSAMT_SECTION`."""
    PYCSAMT_SECTION.configure(**kw)


def reset_section() -> None:
    """Reset :data:`PYCSAMT_SECTION` to package defaults."""
    PYCSAMT_SECTION.reset()


__all__ = [
    "PYCSAMT_SECTION",
    "PyCSAMTSection",
    "SectionAxisStyle",
    "SectionColorbarStyle",
    "SectionFigureStyle",
    "SectionStyle",
    "configure_section",
    "reset_section",
]
