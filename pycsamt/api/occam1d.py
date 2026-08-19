"""Visual style API for native and file-backed Occam1D plots.

The objects in this module contain presentation choices only.  They never
alter inversion data or numerical results.  A style may be selected globally,
passed to one plot, or copied with local overrides::

    from pycsamt.api import PYCSAMT_OCCAM1D

    style = PYCSAMT_OCCAM1D.default.copy(
        observed__marker="x",
        predicted__color="navy",
        model__linewidth=3.0,
        target__visible=False,
        response_legend=True,
    )

Double-underscore paths address fields in nested artist styles.
"""

from __future__ import annotations

import copy
from contextlib import contextmanager
from dataclasses import dataclass, field, fields
from typing import Any

__all__ = [
    "Occam1DArtistStyle",
    "Occam1DPlotStyle",
    "PyCSAMTOccam1D",
    "PYCSAMT_OCCAM1D",
    "configure_occam1d_style",
    "reset_occam1d_style",
    "resolve_occam1d_style",
    "use_occam1d_style",
]


@dataclass
class Occam1DArtistStyle:
    """Validated Matplotlib properties for one Occam1D visual layer."""

    visible: bool = True
    color: str = "#1f77b4"
    linestyle: str = "-"
    linewidth: float = 1.6
    marker: str | None = None
    markersize: float = 5.0
    markerfacecolor: str | None = None
    markeredgecolor: str | None = None
    alpha: float = 1.0
    zorder: float = 2.0
    label: str | None = None

    def kwargs(self, *, include_marker=True) -> dict[str, Any]:
        """Return non-null Matplotlib keyword arguments."""
        values = {
            "color": self.color,
            "linestyle": self.linestyle,
            "linewidth": self.linewidth,
            "alpha": self.alpha,
            "zorder": self.zorder,
            "label": self.label,
        }
        if include_marker:
            values.update(
                marker=self.marker,
                markersize=self.markersize,
                markerfacecolor=self.markerfacecolor,
                markeredgecolor=self.markeredgecolor,
            )
        return {
            key: value
            for key, value in values.items()
            if value is not None
        }


@dataclass
class Occam1DPlotStyle:
    """Complete visual policy for Occam1D model and diagnostic figures.

    Legends are independently controlled for each panel.  The model legend is
    disabled by default because a single step curve is self-explanatory.
    """

    observed: Occam1DArtistStyle = field(
        default_factory=lambda: Occam1DArtistStyle(
            color="black",
            linestyle="none",
            marker="o",
            markersize=5.0,
            markerfacecolor="black",
            markeredgecolor="black",
            label="Observed",
            zorder=4,
        )
    )
    predicted: Occam1DArtistStyle = field(
        default_factory=lambda: Occam1DArtistStyle(
            color="#b2182b", linewidth=1.8, label="Modeled", zorder=3
        )
    )
    phase_predicted: Occam1DArtistStyle = field(
        default_factory=lambda: Occam1DArtistStyle(
            color="#2166ac", linewidth=1.8, label="Modeled", zorder=3
        )
    )
    model: Occam1DArtistStyle = field(
        default_factory=lambda: Occam1DArtistStyle(
            color="#b2182b", linewidth=2.0, label="Layer model"
        )
    )
    iteration: Occam1DArtistStyle = field(
        default_factory=lambda: Occam1DArtistStyle(
            color="#2166ac",
            marker="o",
            markersize=4.0,
            linewidth=1.7,
            label="Normalized RMS",
        )
    )
    target: Occam1DArtistStyle = field(
        default_factory=lambda: Occam1DArtistStyle(
            color="#444444",
            linestyle="--",
            linewidth=1.1,
            label="Target",
            zorder=1,
        )
    )
    roughness: Occam1DArtistStyle = field(
        default_factory=lambda: Occam1DArtistStyle(
            color="#d97706",
            marker="s",
            markersize=3.0,
            linewidth=1.2,
            label="Roughness",
        )
    )
    current_iteration: Occam1DArtistStyle = field(
        default_factory=lambda: Occam1DArtistStyle(
            color="0.75", linestyle=":", linewidth=0.9, label=None
        )
    )
    errorbar_capsize: float = 2.0
    errorbar_linewidth: float = 0.8
    grid_alpha: float = 0.22
    grid_linewidth: float = 0.6
    model_legend: bool = False
    response_legend: bool = True
    convergence_legend: bool = True
    roughness_legend: bool = False
    legend_frame: bool = False
    legend_fontsize: float = 8.0
    legend_location: str = "best"

    def copy(self, **overrides: Any) -> Occam1DPlotStyle:
        """Return a deep copy modified by double-underscore field paths."""
        result = copy.deepcopy(self)
        for path, value in overrides.items():
            parts = path.split("__")
            target = result
            for name in parts[:-1]:
                if not hasattr(target, name):
                    raise ValueError(f"Unknown Occam1D style path {path!r}.")
                target = getattr(target, name)
            name = parts[-1]
            if not hasattr(target, name):
                raise ValueError(f"Unknown Occam1D style path {path!r}.")
            setattr(target, name, value)
        result.validate()
        return result

    def validate(self) -> Occam1DPlotStyle:
        """Validate booleans and finite non-negative numeric properties."""
        artists = (
            self.observed,
            self.predicted,
            self.phase_predicted,
            self.model,
            self.iteration,
            self.target,
            self.roughness,
            self.current_iteration,
        )
        for artist in artists:
            if not isinstance(artist, Occam1DArtistStyle):
                raise TypeError(
                    "Occam1D artist fields require Occam1DArtistStyle."
                )
            for name in ("linewidth", "markersize", "alpha", "zorder"):
                value = getattr(artist, name)
                if isinstance(value, bool) or not isinstance(
                    value, (int, float)
                ):
                    raise TypeError(f"{name} must be a real number.")
                if name != "zorder" and value < 0:
                    raise ValueError(f"{name} cannot be negative.")
        for item in fields(self):
            if item.name.endswith("_legend") and not isinstance(
                getattr(self, item.name), bool
            ):
                raise TypeError(f"{item.name} must be a bool.")
        return self


def _pycsamt_style() -> Occam1DPlotStyle:
    return Occam1DPlotStyle()


def _publication_style() -> Occam1DPlotStyle:
    return Occam1DPlotStyle().copy(
        observed__markerfacecolor="white",
        predicted__color="0.15",
        phase_predicted__color="0.15",
        model__color="0.1",
        iteration__color="0.1",
        roughness__color="0.5",
    )


def _minimal_style() -> Occam1DPlotStyle:
    return Occam1DPlotStyle().copy(
        response_legend=False,
        convergence_legend=False,
        roughness__visible=False,
        current_iteration__visible=False,
    )


class PyCSAMTOccam1D:
    """Registry and mutable default for Occam1D visual presets."""

    def __init__(self):
        self.pycsamt = _pycsamt_style()
        self.publication = _publication_style()
        self.minimal = _minimal_style()
        self.default = copy.deepcopy(self.pycsamt)

    def style_for(self, preset="pycsamt") -> Occam1DPlotStyle:
        """Return a deep copy of a named preset."""
        key = str(preset).strip().lower().replace("-", "_")
        if key not in {"pycsamt", "publication", "minimal", "default"}:
            raise ValueError(
                "Unknown Occam1D style preset. Available: default, minimal, "
                "publication, pycsamt."
            )
        return copy.deepcopy(getattr(self, key))

    def use(self, preset="pycsamt") -> Occam1DPlotStyle:
        """Replace the live default with one named preset."""
        self.default = self.style_for(preset)
        return copy.deepcopy(self.default)

    def configure(self, **overrides) -> Occam1DPlotStyle:
        """Modify the live default with double-underscore paths."""
        self.default = self.default.copy(**overrides)
        return copy.deepcopy(self.default)

    def reset(self) -> Occam1DPlotStyle:
        """Restore the pyCSAMT house style."""
        return self.use("pycsamt")

    @contextmanager
    def context(self, preset=None, **overrides):
        """Temporarily replace or customize the live default style."""
        previous = copy.deepcopy(self.default)
        try:
            base = self.style_for(preset) if preset else self.default
            self.default = base.copy(**overrides)
            yield self.default
        finally:
            self.default = previous


PYCSAMT_OCCAM1D = PyCSAMTOccam1D()


def resolve_occam1d_style(style=None, **overrides) -> Occam1DPlotStyle:
    """Resolve ``None``, a preset name, or an explicit style instance."""
    if style is None:
        base = PYCSAMT_OCCAM1D.default
    elif isinstance(style, str):
        base = PYCSAMT_OCCAM1D.style_for(style)
    elif isinstance(style, Occam1DPlotStyle):
        base = style
    else:
        raise TypeError(
            "style must be None, a preset name, or Occam1DPlotStyle."
        )
    return copy.deepcopy(base).copy(**overrides)


def configure_occam1d_style(**overrides) -> Occam1DPlotStyle:
    """Configure and return a copy of the live default style."""
    return PYCSAMT_OCCAM1D.configure(**overrides)


def use_occam1d_style(preset="pycsamt") -> Occam1DPlotStyle:
    """Select a named Occam1D style as the live default."""
    return PYCSAMT_OCCAM1D.use(preset)


def reset_occam1d_style() -> Occam1DPlotStyle:
    """Restore the default pyCSAMT Occam1D style."""
    return PYCSAMT_OCCAM1D.reset()
