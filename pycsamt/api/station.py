"""Global station-rendering controls for pyCSAMT figures."""

from __future__ import annotations

import copy
from collections.abc import Generator, Sequence
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Any

import numpy as np

_NICE_STEPS = (1, 2, 3, 4, 5, 8, 10, 12, 15, 20, 25, 50, 100)
_SIDES = {"top", "bottom"}
_PRESETS = {"pseudosection", "inversion", "survey"}


@dataclass
class StationMarkerStyle:
    """Marker style used to draw stations along a profile axis."""

    marker: str = "v"
    size: float = 34.0
    facecolor: str = "white"
    edgecolor: str = "black"
    linewidth: float = 0.9
    alpha: float = 1.0
    offset: float = 1.025
    zorder: int = 5

    def kwargs(self, **overrides: Any) -> dict[str, Any]:
        """Return Matplotlib scatter keyword arguments."""
        out = {
            "marker": self.marker,
            "s": self.size,
            "facecolors": self.facecolor,
            "edgecolors": self.edgecolor,
            "linewidths": self.linewidth,
            "alpha": self.alpha,
            "zorder": self.zorder,
            "clip_on": False,
        }
        out.update(overrides)
        return out


@dataclass
class StationAxisStyle:
    """Station-axis tick, label, and marker configuration."""

    side: str = "top"
    xlabel: str = "Station"
    every: int | str = "auto"
    max_labels: int = 14
    rotation: float = 90.0
    fontsize: int = 7
    tick_length: float = 3.0
    label_tick_length: float = 5.0
    label_pad: float = 8.0
    xlabel_pad: float = 6.0
    show_all_ticks: bool = True
    show_labels: bool = True
    show_markers: bool = True
    marker: StationMarkerStyle = field(default_factory=StationMarkerStyle)

    def compute_every(
        self,
        n: int,
        figwidth_in: float = 10.0,
        max_label_len: int = 4,
    ) -> int:
        """Return the interval between visible station labels."""
        if isinstance(self.every, int):
            return max(1, int(self.every))
        n = max(1, int(n))
        if n <= max(1, int(self.max_labels)):
            return 1
        width_factor = max(float(figwidth_in), 1.0) / 10.0
        label_budget = max(2.0, float(self.max_labels) * width_factor)
        needed = max(1.0, n / label_budget)
        if int(max_label_len) > 8:
            needed *= 1.25
        for step in _NICE_STEPS:
            if step >= needed:
                return int(step)
        return int(np.ceil(needed))

    def label_indices(
        self,
        labels: Sequence[Any],
        figwidth_in: float = 10.0,
    ) -> np.ndarray:
        """Return indices whose station labels should be visible."""
        n = len(labels)
        if n <= 0:
            return np.asarray([], dtype=int)
        max_len = max((len(str(label)) for label in labels), default=4)
        step = self.compute_every(n, figwidth_in, max_len)
        idx = np.arange(0, n, step, dtype=int)
        if n - 1 not in idx:
            idx = np.r_[idx, n - 1]
        return idx

    def apply(
        self,
        ax: Any,
        positions: Sequence[float],
        labels: Sequence[Any] | None = None,
        *,
        xlim: tuple[float, float] | None = None,
        topo_elev: np.ndarray | None = None,
    ) -> np.ndarray:
        """Apply station ticks, labels, and markers to ``ax``.

        Parameters
        ----------
        ax : Axes
        positions : (n,) station x positions
        labels : (n,) station name strings, optional
        xlim : (xmin, xmax) x-axis limits, optional
        topo_elev : (n,) terrain elevation in *data y-units*, optional
            When provided the markers are placed at the terrain surface
            (+ a small pad toward the viewer) instead of at the flat
            axis edge.  Labels are drawn inline above the markers;
            xticklabels at the axis edge are suppressed.

        Returns
        -------
        numpy.ndarray
            Indices of stations whose labels are visible.
        """
        side = self.side.lower()
        if side not in _SIDES:
            msg = f"station axis side must be one of {sorted(_SIDES)}."
            raise ValueError(msg)
        x = np.asarray(positions, dtype=float)
        if x.size == 0:
            return np.asarray([], dtype=int)
        raw = labels if labels is not None else x.tolist()
        raw_labels = [str(label) for label in raw]
        figwidth = ax.figure.get_figwidth() if ax.figure is not None else 10.0
        idx = self.label_indices(raw_labels, figwidth)

        if xlim is not None:
            ax.set_xlim(*xlim)

        # ── Topo-surface mode ─────────────────────────────────────────────
        # Markers are placed at the real terrain elevation (data coordinates)
        # instead of at the flat axis edge.  Labels are drawn inline above
        # the markers; the axis ticklabels at the edge are suppressed.
        if topo_elev is not None:
            elev = np.asarray(topo_elev, dtype=float)

            # Compute pad: move marker slightly above the topo surface
            # toward the viewer (upper part of the display).
            y0, y1 = ax.get_ylim()
            y_range = abs(y1 - y0)
            from pycsamt.topo.config import PYCSAMT_TOPO

            pad_frac = PYCSAMT_TOPO.marker_pad_fraction
            # Direction toward the top of the display (+1 if normal, -1 if inverted)
            toward_top = 1.0 if y1 >= y0 else -1.0
            pad = toward_top * pad_frac * y_range

            marker_y = elev + pad
            label_y = elev + pad * 3.5  # label further above the marker

            # Suppress axis ticklabels — labels are drawn inline instead
            ax.set_xticks(x, minor=True)
            ax.tick_params(axis="x", which="minor", length=self.tick_length, top=True)
            ax.set_xticks([])
            ax.set_xticklabels([])

            # All markers always visible
            if self.show_markers:
                ax.scatter(
                    x,
                    marker_y,
                    transform=ax.transData,
                    **self.marker.kwargs(),
                )

            # Labels only at thinned subset
            if self.show_labels:
                ha = "center"
                for i in idx:
                    ax.text(
                        x[i],
                        label_y[i],
                        raw_labels[i],
                        fontsize=self.fontsize,
                        rotation=self.rotation,
                        va="bottom" if toward_top > 0 else "top",
                        ha=ha,
                        clip_on=True,
                        zorder=self.marker.zorder + 1,
                    )
            return idx

        # ── Flat-datum mode (default) ─────────────────────────────────────
        if self.show_all_ticks:
            ax.set_xticks(x, minor=True)
            ax.tick_params(axis="x", which="minor", length=self.tick_length)
        ax.set_xticks(x[idx])
        if self.show_labels:
            # "right" anchoring makes a *diagonal* label hang naturally
            # from its tick (the classic ha="right" + rotation idiom),
            # but at rotation~90 (fully vertical) that same anchor,
            # combined with matplotlib's own top-axis default va
            # ("bottom"), drags the whole label a few pixels left of
            # its tick instead of centering it. Only use "right" in the
            # genuinely diagonal band; center it at 0 and at ~90.
            rot_mod = abs(float(self.rotation)) % 180.0
            ha = "right" if 20.0 < rot_mod < 70.0 else "center"
            ax.set_xticklabels(
                [raw_labels[i] for i in idx],
                rotation=self.rotation,
                fontsize=self.fontsize,
                ha=ha,
            )
        else:
            ax.set_xticklabels(["" for _ in idx])
        ax.tick_params(axis="x", which="major", length=self.label_tick_length)
        if side == "top":
            ax.xaxis.tick_top()
            ax.xaxis.set_label_position("top")
            ax.tick_params(
                axis="x",
                which="both",
                top=True,
                bottom=False,
                labeltop=self.show_labels,
                labelbottom=False,
                pad=self.label_pad,
            )
        else:
            ax.xaxis.tick_bottom()
            ax.xaxis.set_label_position("bottom")
            ax.tick_params(
                axis="x",
                which="both",
                top=False,
                bottom=True,
                labeltop=False,
                labelbottom=self.show_labels,
                pad=self.label_pad,
            )
        ax.set_xlabel(self.xlabel, labelpad=self.xlabel_pad)
        if self.show_markers:
            ax.scatter(
                x,
                np.full(x.size, self.marker.offset),
                transform=ax.get_xaxis_transform(),
                **self.marker.kwargs(),
            )
        return idx


class PyCSAMTStationRendering:
    """Package-wide station rendering control container."""

    def __init__(self) -> None:
        self._init_defaults()

    def _init_defaults(self) -> None:
        self.pseudosection = StationAxisStyle(
            side="top",
            label_pad=10.0,
            xlabel_pad=8.0,
            marker=StationMarkerStyle(
                marker="v",
                facecolor="white",
                edgecolor="black",
                size=34.0,
            ),
        )
        self.inversion = StationAxisStyle(
            side="top",
            label_pad=10.0,
            xlabel_pad=8.0,
            marker=StationMarkerStyle(
                marker="v",
                facecolor="black",
                edgecolor="black",
                size=38.0,
            ),
        )
        self.survey = StationAxisStyle(
            side="bottom",
            rotation=45.0,
            marker=StationMarkerStyle(
                marker="o",
                facecolor="white",
                edgecolor="black",
                size=28.0,
                offset=-0.055,
            ),
        )

    def style_for(self, preset: str = "pseudosection") -> StationAxisStyle:
        """Return the station-axis style associated with ``preset``."""
        preset = str(preset).lower()
        if preset not in _PRESETS:
            msg = f"station preset must be one of {sorted(_PRESETS)}."
            raise ValueError(msg)
        return getattr(self, preset)

    def apply(
        self,
        ax: Any,
        positions: Sequence[float],
        labels: Sequence[Any] | None = None,
        *,
        preset: str = "pseudosection",
        xlim: tuple[float, float] | None = None,
    ) -> np.ndarray:
        """Apply one named station-rendering preset to ``ax``."""
        return self.style_for(preset).apply(ax, positions, labels, xlim=xlim)

    def configure(self, **kw: Any) -> None:
        """Configure station rendering using ``section__attribute`` paths."""
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
    ) -> Generator[PyCSAMTStationRendering, None, None]:
        """Temporarily override station rendering, then restore it."""
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
        source = copy.deepcopy(self.style_for(preset))
        self.pseudosection = source

    def reset(self) -> None:
        """Reset station rendering to package defaults."""
        self._init_defaults()

    def summary(self) -> str:
        """Return a human-readable summary."""
        lines = ["PyCSAMTStationRendering"]
        for name in ("pseudosection", "inversion", "survey"):
            style = self.style_for(name)
            marker = style.marker
            lines.append(
                "  "
                f"{name}: side={style.side!r}, marker={marker.marker!r}, "
                f"face={marker.facecolor!r}, labels<={style.max_labels}",
            )
        return "\n".join(lines)

    def __repr__(self) -> str:  # noqa: D105
        return self.summary()

    def _snapshot(self) -> dict[str, Any]:
        return {
            "pseudosection": copy.deepcopy(self.pseudosection),
            "inversion": copy.deepcopy(self.inversion),
            "survey": copy.deepcopy(self.survey),
        }

    def _restore(self, snapshot: dict[str, Any]) -> None:
        self.pseudosection = snapshot["pseudosection"]
        self.inversion = snapshot["inversion"]
        self.survey = snapshot["survey"]


PYCSAMT_STATION_RENDERING = PyCSAMTStationRendering()


def configure_station_rendering(**kw: Any) -> None:
    """Configure :data:`PYCSAMT_STATION_RENDERING`."""
    PYCSAMT_STATION_RENDERING.configure(**kw)


def reset_station_rendering() -> None:
    """Reset :data:`PYCSAMT_STATION_RENDERING` to package defaults."""
    PYCSAMT_STATION_RENDERING.reset()


__all__ = [
    "PYCSAMT_STATION_RENDERING",
    "PyCSAMTStationRendering",
    "StationAxisStyle",
    "StationMarkerStyle",
    "configure_station_rendering",
    "reset_station_rendering",
]
