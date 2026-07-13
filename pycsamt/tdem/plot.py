"""Visualization helpers for time-domain EM data."""

from __future__ import annotations

import copy
from dataclasses import dataclass
from typing import Any

import numpy as np

from ..api.control import PYCSAMT_CONTROL
from ..api.plot import add_colorbar, save_fig
from ..api.property import PyCSAMTObject
from ..api.section import PYCSAMT_SECTION, SectionStyle
from ..api.station import PYCSAMT_STATION_RENDERING
from ..api.style import PYCSAMT_STYLE

__all__ = [
    "PlotDecayCurve",
    "PlotElevationProfile",
    "PlotGateProfile",
    "PlotSurveyMap",
    "PlotSurveyOverview",
    "PlotTEMDashboard",
    "PlotTEMAVGSection",
    "PlotTEMZSection",
    "PlotTransformedRho",
    "StationTickConfig",
    "TDEMPlotStyle",
    "plot_decay",
    "plot_elevation_profile",
    "plot_gate_profile",
    "plot_survey_map",
    "plot_survey_overview",
    "plot_tem_dashboard",
    "plot_tem_z_section",
    "plot_temavg_section",
    "plot_transformed_rho",
]


@dataclass(repr=False)
class TDEMPlotStyle(PyCSAMTObject):
    """Shared style values for TDEM figures."""

    primary: str = "#2166ac"
    secondary: str = "#d6604d"
    accent: str = "#1b7837"
    warning: str = "#b2182b"
    grid: str = "#ededed"
    text: str = "#1a1a1a"
    decay_cmap: str = "viridis"
    section_cmap: str = "RdYlBu_r"
    elevation_cmap: str = "terrain"
    figsize_single: tuple[float, float] = (7.0, 5.0)
    figsize_double: tuple[float, float] = (7.0, 8.0)
    figsize_wide: tuple[float, float] = (9.0, 4.5)
    multiline: object | None = None
    verbose: int = 0
    logger: object | None = None

    def colors(self, n: int) -> list[str]:
        """Return ``n`` line colors."""
        multiline = self.multiline or PYCSAMT_STYLE.multiline
        if multiline is not None and hasattr(multiline, "colors"):
            return list(multiline.colors(n))
        base = [
            self.primary,
            self.secondary,
            self.accent,
            self.warning,
            "#762a83",
            "#4d4d4d",
        ]
        return [base[i % len(base)] for i in range(n)]

    def line_kwargs(
        self,
        idx: int,
        n: int,
        **overrides: Any,
    ) -> dict[str, Any]:
        """Return line kwargs from the package multiline style."""
        multiline = self.multiline or PYCSAMT_STYLE.multiline
        if multiline is not None and hasattr(multiline, "line_kwargs"):
            return dict(multiline.line_kwargs(idx, n, **overrides))
        kwargs = {"color": self.colors(n)[idx], "lw": 1.5, "alpha": 0.9}
        kwargs.update(overrides)
        return kwargs


class TDEMPlotBase(PyCSAMTObject):
    """Base class shared by TDEM plot objects."""

    def __init__(
        self,
        *,
        style: TDEMPlotStyle | None = None,
        verbose: int = 0,
        logger: object | None = None,
    ) -> None:
        self.style = style or TDEMPlotStyle()
        self.verbose = int(verbose)
        self.logger = logger

    def _finish(self, fig, *, tight: bool = True):
        """Apply final figure layout and return the figure."""
        if tight:
            fig.tight_layout()
        return fig

    def save(self, fig_or_ax, path, **kwargs):
        """Save a TDEM plot with the global pyCSAMT plot config."""
        return save_fig(fig_or_ax, path, **kwargs)


@dataclass(repr=False)
class StationTickConfig(PyCSAMTObject):
    """Station-axis tick configuration for TDEM plots."""

    every: int | str = "auto"
    rotation: float = 45.0
    fontsize: int = 8
    fmt: str = "{:g}"
    max_ticks: int = 12
    preset: str = "pseudosection"
    side: str | None = None
    show_markers: bool = True
    use_shared_api: bool = True
    verbose: int = 0
    logger: object | None = None

    _nice_steps = (1, 2, 5, 10, 20, 25, 50, 100, 200, 250, 500, 1000)

    def compute_every(
        self,
        n: int,
        figwidth_in: float = 10.0,
        max_label_len: int = 4,
    ) -> int:
        """Return the station tick step."""
        if isinstance(self.every, int):
            return max(1, self.every)

        n = max(1, int(n))
        if n <= max(1, int(self.max_ticks)):
            return 1

        available = max(float(figwidth_in), 1.0) * 0.82
        space_per_station = available / n
        char_width = float(self.fontsize) * 0.006
        rotation_rad = np.deg2rad(abs(float(self.rotation)))
        label_width = char_width * max(int(max_label_len), 2)
        label_height = float(self.fontsize) * 0.014
        effective_width = (
            label_width * abs(np.sin(rotation_rad))
            + label_height * abs(np.cos(rotation_rad))
            + 0.025
        )
        needed = max(effective_width / max(space_per_station, 1e-9), 1.0)
        for step in self._nice_steps:
            if step >= needed:
                return step
        return int(np.ceil(needed))

    def apply(
        self,
        ax,
        stations: np.ndarray,
        *,
        labels: list[str] | None = None,
        xlabel: str = "Station",
        xlim: tuple[float, float] | None = None,
    ) -> None:
        """Apply station ticks to ``ax``."""
        stations = np.asarray(stations, dtype=float)
        if stations.size == 0:
            return
        raw_labels = labels or [self._format_label(v, v) for v in stations]
        if self.use_shared_api:
            style = copy.copy(
                PYCSAMT_STATION_RENDERING.style_for(self.preset),
            )
            style.every = self.every
            style.max_labels = int(self.max_ticks)
            style.rotation = float(self.rotation)
            style.fontsize = int(self.fontsize)
            style.show_markers = bool(self.show_markers)
            style.xlabel = xlabel
            if self.side is not None:
                style.side = self.side
            style.apply(
                ax,
                stations,
                raw_labels,
                xlim=xlim,
            )
            return
        figwidth = ax.figure.get_figwidth() if ax.figure is not None else 10.0
        max_len = max((len(str(lbl)) for lbl in raw_labels), default=4)
        step = self.compute_every(stations.size, figwidth, max_len)
        tick_labels = [
            self._format_label(raw_labels[i], stations[i])
            if i % step == 0 or i == stations.size - 1
            else ""
            for i in range(stations.size)
        ]
        ha = "right" if abs(float(self.rotation)) > 20.0 else "center"
        ax.set_xticks(stations)
        ax.set_xticklabels(
            tick_labels,
            rotation=self.rotation,
            fontsize=self.fontsize,
            ha=ha,
        )
        if xlim is not None:
            ax.set_xlim(*xlim)
        ax.set_xlabel(xlabel)

    def _format_label(self, label: Any, value: float) -> str:
        """Format a station label with numeric fallback."""
        try:
            return self.fmt.format(label)
        except (TypeError, ValueError):
            try:
                return self.fmt.format(value)
            except (TypeError, ValueError):
                return str(label)


class PlotDecayCurve(TDEMPlotBase):
    r"""Plot one or more TEM decay curves on log-log axes."""

    def __init__(
        self,
        soundings,
        *,
        title: str | None = None,
        figsize: tuple[float, float] | None = None,
        y_mode: str = "dBdt",
        show_error: bool = True,
        style: TDEMPlotStyle | None = None,
        verbose: int = 0,
        logger: object | None = None,
    ) -> None:
        super().__init__(style=style, verbose=verbose, logger=logger)
        from ._base import TEMSounding

        if isinstance(soundings, TEMSounding):
            soundings = [soundings]
        self.soundings = list(soundings)
        self.title = title
        self.figsize = figsize or self.style.figsize_single
        self.y_mode = y_mode
        self.show_error = bool(show_error)

    def plot(self, ax=None):
        """Draw the decay curves and return the axes."""
        import matplotlib.pyplot as plt

        if ax is None:
            _fig, ax = plt.subplots(figsize=self.figsize)
        for idx, sounding in enumerate(self.soundings):
            y = _sounding_decay_values(sounding, self.y_mode)
            t_ms = np.asarray(sounding.time_gates, dtype=float) * 1e3
            label = sounding.station_name or f"Sounding {idx + 1}"
            line_kw = self.style.line_kwargs(
                idx,
                len(self.soundings),
                marker="o",
                ms=4,
                lw=1.2,
                label=label,
            )
            if (
                self.show_error
                and sounding.error is not None
                and self.y_mode == "data"
            ):
                ax.errorbar(
                    t_ms,
                    np.abs(y),
                    yerr=np.abs(sounding.error),
                    fmt="o-",
                    capsize=2,
                    **line_kw,
                )
            else:
                ax.plot(
                    t_ms,
                    np.abs(y),
                    **line_kw,
                )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Time (ms)")
        ax.set_ylabel(_decay_ylabel(self.y_mode))
        ax.set_title(self.title or _default_decay_title(self.soundings))
        ax.grid(True, which="both", color=self.style.grid, linewidth=0.7)
        ax.legend(fontsize=8)
        return ax


class PlotTransformedRho(TDEMPlotBase):
    """Plot transformed apparent resistivity and optional phase."""

    def __init__(
        self,
        data,
        *,
        show_phase: bool = True,
        freq_convention: str = "skin_depth",
        phase_mode: str = "weidelt",
        use_control: bool = False,
        panel_height_ratios: tuple[float, float] = (3.0, 1.15),
        figsize: tuple[float, float] | None = None,
        style: TDEMPlotStyle | None = None,
        verbose: int = 0,
        logger: object | None = None,
    ) -> None:
        super().__init__(style=style, verbose=verbose, logger=logger)
        self.show_phase = bool(show_phase)
        self.freq_convention = freq_convention
        self.phase_mode = phase_mode
        self.use_control = bool(use_control)
        self.panel_height_ratios = panel_height_ratios
        self.figsize = figsize or (
            self.style.figsize_double
            if self.show_phase
            else self.style.figsize_single
        )
        self.results = _as_transform_results(
            data,
            freq_convention=freq_convention,
            phase_mode=phase_mode,
        )

    def plot(self, axes=None):
        """Draw apparent resistivity and phase panels."""
        import matplotlib.pyplot as plt

        n_panels = 2 if self.show_phase else 1
        created = axes is None
        if axes is None:
            if self.show_phase:
                fig, axes = plt.subplots(
                    2,
                    1,
                    figsize=self.figsize,
                    sharex=True,
                    gridspec_kw={
                        "height_ratios": self.panel_height_ratios,
                        "hspace": 0.08,
                    },
                )
            else:
                fig, ax = plt.subplots(figsize=self.figsize)
                axes = [ax]
        else:
            fig = axes[0].figure if hasattr(axes, "__len__") else axes.figure
            if not hasattr(axes, "__len__"):
                axes = [axes]

        ax_rho = axes[0]
        ax_phase = axes[1] if self.show_phase else None
        for idx, result in enumerate(self.results):
            freq = np.asarray(result["freq"], dtype=float)
            rho = np.asarray(result["rho_a"], dtype=float)
            x = freq
            y_rho = rho
            if self.use_control:
                x = PYCSAMT_CONTROL.x.transform(freq)
                y_rho = PYCSAMT_CONTROL.rho.transform(rho)
            label = result.get("station_name") or f"Site {idx + 1}"
            line_kw = self.style.line_kwargs(
                idx,
                len(self.results),
                marker="o",
                ms=4,
                lw=1.2,
                label=label,
            )
            if self.use_control:
                ax_rho.plot(x, y_rho, **line_kw)
            else:
                ax_rho.loglog(x, y_rho, **line_kw)
            phase = result.get("phase_xy")
            if ax_phase is not None and phase is not None:
                y_phase = phase
                if self.use_control:
                    y_phase = PYCSAMT_CONTROL.phase.transform(phase)
                phase_kw = self.style.line_kwargs(
                    idx,
                    len(self.results),
                    marker="s",
                    ms=4,
                    lw=1.0,
                    ls="--",
                )
                if self.use_control:
                    ax_phase.plot(x, y_phase, **phase_kw)
                else:
                    ax_phase.semilogx(x, y_phase, **phase_kw)
        if self.use_control:
            ax_rho.set_ylabel(PYCSAMT_CONTROL.rho.label())
            xlabel = PYCSAMT_CONTROL.x.label()
            if PYCSAMT_CONTROL.x.use_log_scale():
                ax_rho.set_xscale("log")
                if ax_phase is not None:
                    ax_phase.set_xscale("log")
        else:
            ax_rho.set_ylabel(r"$\rho_a$ ($\Omega\,\mathrm{m}$)")
            xlabel = "Pseudo-frequency (Hz)"
        ax_rho.set_title("TEM-derived apparent resistivity")
        ax_rho.grid(True, which="both", color=self.style.grid)
        ax_rho.legend(fontsize=8)
        if ax_phase is not None:
            ax_phase.set_xlabel(xlabel)
            ax_phase.set_ylabel(
                PYCSAMT_CONTROL.phase.label()
                if self.use_control
                else r"Phase $Z_{xy}$ ($^\circ$)"
            )
            if not self.use_control:
                ax_phase.set_ylim(0, 90)
            ax_phase.axhline(
                45,
                color="#777777",
                linestyle=":",
                linewidth=0.8,
            )
            ax_phase.grid(True, which="both", color=self.style.grid)
        else:
            ax_rho.set_xlabel(xlabel)
        if created:
            if self.show_phase:
                fig.subplots_adjust(
                    left=0.12,
                    right=0.97,
                    bottom=0.10,
                    top=0.95,
                    hspace=0.08,
                )
                return self._finish(fig, tight=False)
            return self._finish(fig)
        return tuple(axes[:n_panels])


class PlotTEMAVGSection(TDEMPlotBase):
    """Plot a TEMAVG pseudo-section from station-gate records."""

    def __init__(
        self,
        data,
        *,
        value: str = "ramp_app_res",
        y: str = "depth",
        log_value: bool = True,
        absolute: bool = False,
        cmap: str | None = None,
        title: str | None = None,
        section: str | SectionStyle = "dynamic",
        figsize: tuple[float, float] | None = None,
        station_ticks: StationTickConfig | None = None,
        style: TDEMPlotStyle | None = None,
        verbose: int = 0,
        logger: object | None = None,
    ) -> None:
        super().__init__(style=style, verbose=verbose, logger=logger)
        self.data = data
        self.value = value
        self.y = y
        self.log_value = bool(log_value)
        self.absolute = bool(absolute)
        self.cmap = cmap or self.style.section_cmap
        self.title = title
        self.section_style = _resolve_section_style(section)
        self.figsize = figsize
        self.station_ticks = station_ticks or StationTickConfig(
            preset=self.section_style.station_preset,
            side=self.section_style.axis.station_side,
        )

    def plot(self, ax=None, *, colorbar: bool = True):
        """Draw the pseudo-section and return the axes."""
        import matplotlib.pyplot as plt

        section = _records_to_section(
            _avg_records(self.data),
            x_key="station",
            y_key=self.y,
            value_key=self.value,
            absolute=self.absolute,
            log_value=self.log_value,
        )
        if ax is None:
            _fig, ax = plt.subplots(
                figsize=self.figsize
                or self.section_style.figsize_for(
                    n_stations=section["x"].size,
                    n_y=section["y"].size,
                    colorbar=colorbar,
                ),
            )
        mesh = ax.pcolormesh(
            section["x_edges"],
            section["y_edges"],
            section["values"],
            cmap=self.cmap,
            shading="auto",
        )
        self.section_style.apply_axis(
            ax,
            xlabel="Station",
            ylabel=_axis_label(self.y),
            title=self.title or _section_title(self.value),
        )
        self.station_ticks.apply(ax, section["x"], xlabel="Station")
        if colorbar:
            self.section_style.add_colorbar(
                mesh,
                ax,
                label=_value_label(self.value, self.log_value),
            )
        return ax


class PlotTEMZSection(TDEMPlotBase):
    """Plot a ZPLOT ``.Z`` pseudo-section."""

    def __init__(
        self,
        data,
        *,
        value: str = "magnitude",
        y: str = "time_ms",
        log_value: bool = True,
        absolute: bool = True,
        cmap: str | None = None,
        title: str | None = None,
        section: str | SectionStyle = "dynamic",
        figsize: tuple[float, float] | None = None,
        station_ticks: StationTickConfig | None = None,
        style: TDEMPlotStyle | None = None,
        verbose: int = 0,
        logger: object | None = None,
    ) -> None:
        super().__init__(style=style, verbose=verbose, logger=logger)
        self.data = data
        self.value = value
        self.y = y
        self.log_value = bool(log_value)
        self.absolute = bool(absolute)
        self.cmap = cmap or self.style.decay_cmap
        self.title = title
        self.section_style = _resolve_section_style(section)
        self.figsize = figsize
        self.station_ticks = station_ticks or StationTickConfig(
            preset=self.section_style.station_preset,
            side=self.section_style.axis.station_side,
        )

    def plot(self, ax=None, *, colorbar: bool = True):
        """Draw the ZPLOT pseudo-section and return the axes."""
        import matplotlib.pyplot as plt

        section = _records_to_section(
            _z_records(self.data),
            x_key="station",
            y_key=self.y,
            value_key=self.value,
            absolute=self.absolute,
            log_value=self.log_value,
        )
        if ax is None:
            _fig, ax = plt.subplots(
                figsize=self.figsize
                or self.section_style.figsize_for(
                    n_stations=section["x"].size,
                    n_y=section["y"].size,
                    colorbar=colorbar,
                ),
            )
        mesh = ax.pcolormesh(
            section["x_edges"],
            section["y_edges"],
            section["values"],
            cmap=self.cmap,
            shading="auto",
        )
        self.section_style.apply_axis(
            ax,
            xlabel="Station",
            ylabel=_axis_label(self.y),
            title=self.title or _section_title(self.value),
        )
        self.station_ticks.apply(ax, section["x"], xlabel="Station")
        if colorbar:
            self.section_style.add_colorbar(
                mesh,
                ax,
                label=_value_label(self.value, self.log_value),
            )
        return ax


class PlotSurveyMap(TDEMPlotBase):
    """Plot TEM station coordinates from a survey or coordinate table."""

    def __init__(
        self,
        data,
        *,
        color_by: str = "elevation",
        annotate: bool = False,
        contour: bool = False,
        contour_levels: int | list[float] = 8,
        contour_labels: bool = True,
        marker_preset: str = "survey",
        marker_size: float | None = None,
        marker_alpha: float | None = None,
        padding: float = 0.04,
        colorbar_size: str = "3.5%",
        colorbar_pad: float = 0.04,
        colorbar_max_ticks: int | None = 5,
        cmap: str | None = None,
        title: str | None = "TEM survey map",
        figsize: tuple[float, float] | None = None,
        style: TDEMPlotStyle | None = None,
        verbose: int = 0,
        logger: object | None = None,
    ) -> None:
        super().__init__(style=style, verbose=verbose, logger=logger)
        self.data = data
        self.color_by = color_by
        self.annotate = bool(annotate)
        self.contour = bool(contour)
        self.contour_levels = contour_levels
        self.contour_labels = bool(contour_labels)
        self.marker_preset = marker_preset
        self.marker_size = marker_size
        self.marker_alpha = marker_alpha
        self.padding = float(padding)
        self.colorbar_size = colorbar_size
        self.colorbar_pad = float(colorbar_pad)
        self.colorbar_max_ticks = colorbar_max_ticks
        self.cmap = cmap or self.style.elevation_cmap
        self.title = title
        self.figsize = figsize or self.style.figsize_single

    def plot(self, ax=None, *, colorbar: bool = True):
        """Draw the survey station map and return the axes."""
        import matplotlib.pyplot as plt

        rows = _coordinate_records(self.data)
        if not rows:
            msg = "No coordinate rows are available for survey map."
            raise ValueError(msg)
        x = np.asarray([row["x"] for row in rows], dtype=float)
        y = np.asarray([row["y"] for row in rows], dtype=float)
        c = _as_float_array([row.get(self.color_by, np.nan) for row in rows])
        if ax is None:
            _fig, ax = plt.subplots(figsize=self.figsize)
        marker = PYCSAMT_STATION_RENDERING.style_for(
            self.marker_preset,
        ).marker
        scatter_kw = {
            "marker": marker.marker,
            "s": self.marker_size or marker.size,
            "edgecolors": marker.edgecolor,
            "linewidths": marker.linewidth,
            "alpha": marker.alpha
            if self.marker_alpha is None
            else self.marker_alpha,
            "zorder": marker.zorder,
        }
        sc = ax.scatter(
            x,
            y,
            c=c if np.isfinite(c.astype(float)).any() else None,
            cmap=self.cmap,
            **scatter_kw,
        )
        if self.contour:
            self._add_contours(ax, x, y, c)
        if self.annotate:
            for row in rows:
                label = row.get("point", row.get("station", ""))
                ax.annotate(
                    str(label),
                    (row["x"], row["y"]),
                    xytext=(2.0, 2.0),
                    textcoords="offset points",
                    fontsize=7,
                )
        ax.set_aspect("equal", adjustable="box")
        _set_padded_limits(ax, x, y, self.padding)
        ax.set_xlabel("Relative X (m)")
        ax.set_ylabel("Relative Y (m)")
        if self.title:
            ax.set_title(self.title)
        ax.grid(True, color=self.style.grid)
        if colorbar and np.isfinite(c).any():
            add_colorbar(
                sc,
                ax,
                label=_axis_label(self.color_by),
                size=self.colorbar_size,
                pad=self.colorbar_pad,
                max_ticks=self.colorbar_max_ticks,
            )
        return ax

    def _add_contours(self, ax, x, y, c) -> None:
        """Overlay contours for the selected coordinate attribute."""
        finite = np.isfinite(x) & np.isfinite(y) & np.isfinite(c)
        if finite.sum() < 4:
            return
        values = c[finite]
        if np.nanmax(values) <= np.nanmin(values):
            return
        try:
            cs = ax.tricontour(
                x[finite],
                y[finite],
                values,
                levels=self.contour_levels,
                colors="#2f2f2f",
                linewidths=0.75,
                alpha=0.82,
                zorder=6,
            )
        except (RuntimeError, ValueError):
            return
        if self.contour_labels:
            ax.clabel(cs, inline=True, fontsize=7, fmt="%.0f")


class PlotElevationProfile(TDEMPlotBase):
    """Plot TEM station elevation along one or more survey profiles."""

    def __init__(
        self,
        data,
        *,
        profiles: float | list[float] | None = None,
        x: str = "point",
        station_ticks: StationTickConfig | None = None,
        title: str | None = "TEM elevation profile",
        figsize: tuple[float, float] | None = None,
        style: TDEMPlotStyle | None = None,
        verbose: int = 0,
        logger: object | None = None,
    ) -> None:
        super().__init__(style=style, verbose=verbose, logger=logger)
        self.data = data
        self.profiles = _as_profile_list(profiles)
        self.x = x
        self.station_ticks = station_ticks or StationTickConfig(
            every="auto",
            rotation=45.0,
            max_ticks=14,
            preset="survey",
            side="top",
            show_markers=False,
        )
        self.title = title
        self.figsize = figsize or self.style.figsize_wide

    def plot(self, ax=None):
        """Draw station elevation profiles and return the axes."""
        import matplotlib.pyplot as plt

        rows = _coordinate_records(self.data)
        grouped = _coordinate_profiles(rows, self.profiles)
        if not grouped:
            msg = "No coordinate rows are available for elevation profile."
            raise ValueError(msg)
        if ax is None:
            _fig, ax = plt.subplots(figsize=self.figsize)

        for idx, (profile, prows) in enumerate(grouped.items()):
            x_values, labels, xlabel = _profile_x_values(prows, self.x)
            elev = np.asarray(
                [float(row["elevation"]) for row in prows],
                dtype=float,
            )
            line_kw = self.style.line_kwargs(
                idx,
                len(grouped),
                marker="o",
                ms=3.0,
                lw=1.2,
                label=f"P{profile:g}",
            )
            ax.plot(x_values, elev, **line_kw)

        if len(grouped) == 1:
            self.station_ticks.apply(
                ax,
                x_values,
                labels=labels,
                xlabel=xlabel,
            )
        else:
            ax.set_xlabel(_profile_axis_label(self.x))
            ax.legend(title="Profile", fontsize=8)
        ax.set_ylabel("Elevation (m)")
        if self.title:
            ax.set_title(self.title)
        ax.grid(True, color=self.style.grid)
        return ax


class PlotSurveyOverview(TDEMPlotBase):
    """Plot a TEM survey map with a matched elevation profile panel."""

    def __init__(
        self,
        data,
        *,
        profile: float | None = None,
        profile_x: str = "point",
        map_kwargs: dict[str, Any] | None = None,
        profile_kwargs: dict[str, Any] | None = None,
        title: str | None = "TEM survey overview",
        figsize: tuple[float, float] = (10.0, 6.5),
        height_ratios: tuple[float, float] = (1.25, 1.0),
        style: TDEMPlotStyle | None = None,
        verbose: int = 0,
        logger: object | None = None,
    ) -> None:
        super().__init__(style=style, verbose=verbose, logger=logger)
        self.data = data
        self.profile = profile
        self.profile_x = profile_x
        self.map_kwargs = dict(map_kwargs or {})
        self.profile_kwargs = dict(profile_kwargs or {})
        self.title = title
        self.figsize = figsize
        self.height_ratios = height_ratios

    def plot(self):
        """Draw the survey overview and return the figure."""
        import matplotlib.pyplot as plt

        profile = self.profile
        if profile is None:
            profiles = _coordinate_profiles(
                _coordinate_records(self.data), None
            )
            if profiles:
                profile = next(iter(profiles))

        fig, axes = plt.subplots(
            2,
            1,
            figsize=self.figsize,
            gridspec_kw={
                "height_ratios": self.height_ratios,
                "hspace": 0.48,
            },
        )
        map_kwargs = {
            "color_by": "elevation",
            "contour": True,
            "marker_size": 8,
            "marker_alpha": 0.72,
            "padding": 0.025,
            "title": None,
            **self.map_kwargs,
        }
        PlotSurveyMap(
            self.data,
            style=self.style,
            **map_kwargs,
        ).plot(ax=axes[0])
        axes[0].set_xlabel("")

        profile_kwargs = {
            "profiles": profile,
            "x": self.profile_x,
            "title": None,
            **self.profile_kwargs,
        }
        PlotElevationProfile(
            self.data,
            style=self.style,
            **profile_kwargs,
        ).plot(ax=axes[1])
        if self.title:
            fig.suptitle(self.title, y=0.98)
        fig.subplots_adjust(top=0.90, hspace=0.48)
        return fig


class PlotGateProfile(TDEMPlotBase):
    """Plot selected TEMAVG windows as profiles along stations."""

    def __init__(
        self,
        data,
        *,
        windows: list[int] | None = None,
        value: str = "magnitude",
        absolute: bool = True,
        log_y: bool = True,
        title: str | None = None,
        figsize: tuple[float, float] | None = None,
        station_ticks: StationTickConfig | None = None,
        style: TDEMPlotStyle | None = None,
        verbose: int = 0,
        logger: object | None = None,
    ) -> None:
        super().__init__(style=style, verbose=verbose, logger=logger)
        self.data = data
        self.windows = windows
        self.value = value
        self.absolute = bool(absolute)
        self.log_y = bool(log_y)
        self.title = title
        self.figsize = figsize or self.style.figsize_wide
        self.station_ticks = station_ticks or StationTickConfig(
            rotation=0.0,
        )

    def plot(self, ax=None):
        """Draw selected gate profiles and return the axes."""
        import matplotlib.pyplot as plt

        rows = _avg_records(self.data)
        windows = self.windows or _representative_windows(
            row["window"] for row in rows
        )
        if ax is None:
            _fig, ax = plt.subplots(figsize=self.figsize)
        all_stations = np.asarray(
            sorted({float(row["station"]) for row in rows}),
            dtype=float,
        )
        for idx, window in enumerate(windows):
            sub = [row for row in rows if int(row["window"]) == int(window)]
            sub = sorted(sub, key=lambda row: float(row["station"]))
            stations = np.asarray(
                [float(row["station"]) for row in sub],
                dtype=float,
            )
            values = np.asarray(
                [float(row[self.value]) for row in sub],
                dtype=float,
            )
            if self.absolute:
                values = np.abs(values)
            line_kw = self.style.line_kwargs(
                idx,
                len(windows),
                marker="o",
                ms=3.5,
                lw=1.1,
                label=f"W{int(window)}",
            )
            ax.plot(
                stations,
                values,
                **line_kw,
            )
        if self.log_y:
            ax.set_yscale("log")
        self.station_ticks.apply(ax, all_stations, xlabel="Station")
        ax.set_ylabel(_value_label(self.value, False))
        ax.set_title(
            self.title or f"Gate profiles: {_axis_label(self.value)}"
        )
        ax.grid(True, color=self.style.grid)
        ax.legend(title="Window", fontsize=8, ncol=min(4, len(windows)))
        return ax


class PlotTEMDashboard(TDEMPlotBase):
    """Create a compact multi-panel TDEM real-data dashboard."""

    def __init__(
        self,
        avg,
        zplot,
        soundings,
        *,
        title: str = "TDEM profile dashboard",
        figsize: tuple[float, float] = (12.0, 9.0),
        station_ticks: StationTickConfig | None = None,
        style: TDEMPlotStyle | None = None,
        verbose: int = 0,
        logger: object | None = None,
    ) -> None:
        super().__init__(style=style, verbose=verbose, logger=logger)
        self.avg = avg
        self.zplot = zplot
        self.soundings = list(soundings)
        self.title = title
        self.figsize = figsize
        self.station_ticks = station_ticks or StationTickConfig()

    def plot(self):
        """Draw dashboard and return the figure."""
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=self.figsize)
        selected = _select_evenly(self.soundings, 4)
        PlotDecayCurve(
            selected,
            title="Selected station decays",
            style=self.style,
        ).plot(ax=axes[0, 0])
        PlotTransformedRho(
            selected,
            show_phase=False,
            style=self.style,
        ).plot(axes=axes[0, 1])
        PlotTEMAVGSection(
            self.avg,
            value="ramp_app_res",
            y="depth",
            title="AVG ramp apparent resistivity",
            station_ticks=self.station_ticks,
            style=self.style,
        ).plot(ax=axes[1, 0], colorbar=True)
        PlotTEMZSection(
            self.zplot,
            value="magnitude",
            y="time_ms",
            title="ZPLOT transient magnitude",
            station_ticks=self.station_ticks,
            style=self.style,
        ).plot(ax=axes[1, 1], colorbar=True)
        fig.suptitle(self.title)
        fig.tight_layout()
        return fig


def plot_decay(soundings, **kwargs):
    """Plot TEM decay curves."""
    return PlotDecayCurve(soundings, **kwargs).plot()


def plot_transformed_rho(data, **kwargs):
    """Plot transformed apparent resistivity and phase."""
    return PlotTransformedRho(data, **kwargs).plot()


def plot_gate_profile(data, **kwargs):
    """Plot selected TEMAVG gate profiles."""
    return PlotGateProfile(data, **kwargs).plot()


def plot_elevation_profile(data, **kwargs):
    """Plot station elevation along TEM survey profiles."""
    return PlotElevationProfile(data, **kwargs).plot()


def plot_temavg_section(data, **kwargs):
    """Plot a TEMAVG pseudo-section."""
    return PlotTEMAVGSection(data, **kwargs).plot()


def plot_tem_z_section(data, **kwargs):
    """Plot a TEMAVG ``.Z`` pseudo-section."""
    return PlotTEMZSection(data, **kwargs).plot()


def plot_survey_map(data, **kwargs):
    """Plot station coordinates from a TEM survey."""
    return PlotSurveyMap(data, **kwargs).plot()


def plot_survey_overview(data, **kwargs):
    """Plot a TEM survey map with an elevation-profile panel."""
    return PlotSurveyOverview(data, **kwargs).plot()


def plot_tem_dashboard(avg, zplot, soundings, **kwargs):
    """Plot a compact TDEM dashboard."""
    return PlotTEMDashboard(avg, zplot, soundings, **kwargs).plot()


def _as_transform_results(
    data,
    *,
    freq_convention: str,
    phase_mode: str,
) -> list[dict[str, Any]]:
    """Normalize transform plot inputs to result dictionaries."""
    from ._base import TEMSounding
    from .transform import LateTimeTransform

    if isinstance(data, (TEMSounding, dict)):
        data = [data]
    results: list[dict[str, Any]] = []
    transformer = LateTimeTransform(
        freq_convention=freq_convention,
        phase_mode=phase_mode,
    )
    for item in data:
        if isinstance(item, TEMSounding):
            results.append(transformer.transform(item))
        elif isinstance(item, dict) and "freq" in item:
            results.append(item)
        else:
            msg = (
                "Expected TEMSounding, transform result dict, or "
                f"sequence of those; got {type(item)!r}."
            )
            raise TypeError(msg)
    return results


def _resolve_section_style(section: str | SectionStyle) -> SectionStyle:
    """Return a copied section style for TDEM section plots."""
    if isinstance(section, SectionStyle):
        return section.copy()
    return PYCSAMT_SECTION.style_for(str(section)).copy()


def _sounding_decay_values(sounding, y_mode: str) -> np.ndarray:
    """Return decay values for a sounding."""
    mode = y_mode.lower()
    if mode == "dbdt":
        return sounding.dBdt()
    if mode == "data":
        return np.asarray(sounding.data, dtype=float)
    msg = "y_mode must be 'dBdt' or 'data'."
    raise ValueError(msg)


def _decay_ylabel(y_mode: str) -> str:
    """Return the y-axis label for a decay mode."""
    if y_mode.lower() == "dbdt":
        return r"$|\partial B_z / \partial t|$ (T/s)"
    return "Absolute transient data"


def _default_decay_title(soundings: list[Any]) -> str:
    """Return a default decay title."""
    if len(soundings) == 1:
        return soundings[0].station_name or "TEM decay"
    return "TEM decay curves"


def _select_evenly(values: list[Any], n: int) -> list[Any]:
    """Select up to ``n`` items evenly from a list."""
    if len(values) <= n:
        return list(values)
    idx = np.linspace(0, len(values) - 1, n).round().astype(int)
    return [values[int(i)] for i in idx]


def _representative_windows(windows) -> list[int]:
    """Return early, middle, and late windows for profile plots."""
    values = sorted({int(window) for window in windows})
    if len(values) <= 4:
        return values
    idx = np.linspace(0, len(values) - 1, 4).round().astype(int)
    return [values[int(i)] for i in idx]


def _avg_records(data) -> list[dict[str, Any]]:
    """Return TEMAVG-like records from supported inputs."""
    from .avg import TEMAVG
    from .survey import TEMSurvey

    if isinstance(data, TEMAVG):
        return data.to_records()
    if isinstance(data, TEMSurvey):
        return data.to_records()
    if isinstance(data, list):
        return data
    msg = f"Unsupported TEMAVG section input: {type(data)!r}."
    raise TypeError(msg)


def _z_records(data) -> list[dict[str, Any]]:
    """Return ZPLOT-like records from supported inputs."""
    from .survey import TEMSurvey
    from .zplot import TEMZPlot

    if isinstance(data, TEMZPlot):
        return data.to_records()
    if isinstance(data, TEMSurvey):
        rows: list[dict[str, Any]] = []
        for zplot in data.z_files.values():
            rows.extend(zplot.to_records())
        return rows
    if isinstance(data, list):
        return data
    msg = f"Unsupported TEM .Z section input: {type(data)!r}."
    raise TypeError(msg)


def _coordinate_records(data) -> list[dict[str, Any]]:
    """Return coordinate rows from supported map inputs."""
    from .coordinates import TEMCoordinateTable
    from .survey import TEMSurvey

    if isinstance(data, TEMCoordinateTable):
        return data.to_records()
    if isinstance(data, TEMSurvey):
        if data.coordinates is not None:
            return data.coordinates.to_records()
        rows = data.to_records()
        return [
            row
            for row in rows
            if {"x", "y"}.issubset(row) and row.get("x") is not None
        ]
    if isinstance(data, list):
        return data
    msg = f"Unsupported survey map input: {type(data)!r}."
    raise TypeError(msg)


def _as_profile_list(
    profiles: float | list[float] | None,
) -> list[float] | None:
    """Return profile selectors as floats."""
    if profiles is None:
        return None
    if isinstance(profiles, (str, int, float)):
        return [float(profiles)]
    return [float(profile) for profile in profiles]


def _coordinate_profiles(
    rows: list[dict[str, Any]],
    profiles: list[float] | None,
) -> dict[float, list[dict[str, Any]]]:
    """Group coordinate records by profile id."""
    wanted = set(profiles) if profiles is not None else None
    grouped: dict[float, list[dict[str, Any]]] = {}
    for row in rows:
        profile = float(row.get("profile", 0.0))
        if wanted is not None and profile not in wanted:
            continue
        grouped.setdefault(profile, []).append(row)
    return {
        profile: sorted(
            prows,
            key=lambda row: (
                float(row.get("point", row.get("station", 0.0))),
                float(row.get("x", 0.0)),
                float(row.get("y", 0.0)),
            ),
        )
        for profile, prows in sorted(grouped.items())
    }


def _profile_x_values(
    rows: list[dict[str, Any]],
    mode: str,
) -> tuple[np.ndarray, list[str], str]:
    """Return x values, station labels, and axis label for profile rows."""
    mode = mode.lower()
    labels = [
        f"{float(row.get('point', row.get('station', idx))):g}"
        for idx, row in enumerate(rows)
    ]
    if mode in {"point", "station"}:
        return (
            np.asarray(
                [
                    float(row.get("point", row.get("station", idx)))
                    for idx, row in enumerate(rows)
                ],
                dtype=float,
            ),
            labels,
            "Station point",
        )
    if mode in {"x", "relative_x"}:
        return (
            np.asarray([float(row["x"]) for row in rows], dtype=float),
            labels,
            "Relative X (m)",
        )
    if mode in {"y", "relative_y"}:
        return (
            np.asarray([float(row["y"]) for row in rows], dtype=float),
            labels,
            "Relative Y (m)",
        )
    if mode in {"distance", "dist", "chainage"}:
        return _profile_distance(rows), labels, "Distance along profile (m)"
    msg = "x must be 'point', 'station', 'x', 'y', or 'distance'."
    raise ValueError(msg)


def _profile_distance(rows: list[dict[str, Any]]) -> np.ndarray:
    """Return cumulative station-to-station distance from coordinates."""
    x = np.asarray([float(row["x"]) for row in rows], dtype=float)
    y = np.asarray([float(row["y"]) for row in rows], dtype=float)
    if x.size == 0:
        return x
    step = np.hypot(np.diff(x), np.diff(y))
    return np.r_[0.0, np.cumsum(step)]


def _profile_axis_label(mode: str) -> str:
    """Return the x-axis label for an elevation profile mode."""
    labels = {
        "point": "Station point",
        "station": "Station point",
        "x": "Relative X (m)",
        "relative_x": "Relative X (m)",
        "y": "Relative Y (m)",
        "relative_y": "Relative Y (m)",
        "distance": "Distance along profile (m)",
        "dist": "Distance along profile (m)",
        "chainage": "Distance along profile (m)",
    }
    return labels.get(mode.lower(), mode)


def _records_to_section(
    rows: list[dict[str, Any]],
    *,
    x_key: str,
    y_key: str,
    value_key: str,
    absolute: bool,
    log_value: bool,
) -> dict[str, np.ndarray]:
    """Convert sparse row records to a dense section matrix."""
    if not rows:
        msg = "No records are available for section plotting."
        raise ValueError(msg)
    x_vals = sorted({float(row[x_key]) for row in rows})
    y_vals = sorted({float(row[y_key]) for row in rows})
    x_index = {value: idx for idx, value in enumerate(x_vals)}
    y_index = {value: idx for idx, value in enumerate(y_vals)}
    matrix = np.full((len(y_vals), len(x_vals)), np.nan, dtype=float)
    for row in rows:
        x = float(row[x_key])
        y = float(row[y_key])
        value = float(row[value_key])
        if absolute:
            value = abs(value)
        if log_value:
            value = np.log10(value) if value > 0.0 else np.nan
        matrix[y_index[y], x_index[x]] = value
    return {
        "x": np.asarray(x_vals, dtype=float),
        "y": np.asarray(y_vals, dtype=float),
        "x_edges": _edges(x_vals),
        "y_edges": _edges(y_vals),
        "values": matrix,
    }


def _edges(values: list[float]) -> np.ndarray:
    """Return cell edges from sorted cell centres."""
    arr = np.asarray(values, dtype=float)
    if arr.size == 1:
        return np.asarray([arr[0] - 0.5, arr[0] + 0.5])
    delta = np.diff(arr)
    mids = arr[:-1] + delta / 2.0
    lo = arr[0] - delta[0] / 2.0
    hi = arr[-1] + delta[-1] / 2.0
    return np.concatenate([[lo], mids, [hi]])


def _as_float_array(values: list[Any]) -> np.ndarray:
    """Return numeric values, coercing failures to nan."""
    out = []
    for value in values:
        try:
            out.append(float(value))
        except (TypeError, ValueError):
            out.append(np.nan)
    return np.asarray(out, dtype=float)


def _set_padded_limits(ax, x, y, padding: float) -> None:
    """Set compact equal-aspect map limits around finite coordinates."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    finite = np.isfinite(x) & np.isfinite(y)
    if not finite.any():
        return
    xmin = float(np.nanmin(x[finite]))
    xmax = float(np.nanmax(x[finite]))
    ymin = float(np.nanmin(y[finite]))
    ymax = float(np.nanmax(y[finite]))
    xspan = max(xmax - xmin, 1.0)
    yspan = max(ymax - ymin, 1.0)
    pad = max(float(padding), 0.0)
    ax.set_xlim(xmin - xspan * pad, xmax + xspan * pad)
    ax.set_ylim(ymin - yspan * pad, ymax + yspan * pad)


def _axis_label(name: str) -> str:
    """Return a readable axis label from a field name."""
    labels = {
        "depth": "Depth (m)",
        "time_ms": "Time (ms)",
        "time_s": "Time (s)",
        "elevation": "Elevation (m)",
        "ramp_app_res": r"$\rho_a$ ($\Omega\,\mathrm{m}$)",
        "app_res": r"$\rho_a$ ($\Omega\,\mathrm{m}$)",
        "rho_a": r"$\rho_a$ ($\Omega\,\mathrm{m}$)",
        "phase": r"Phase ($^\circ$)",
        "x": "Relative X (m)",
        "y": "Relative Y (m)",
    }
    return labels.get(name, name.replace("_", " ").title())


def _value_label(name: str, log_value: bool) -> str:
    """Return a readable colorbar label."""
    label = _axis_label(name)
    return rf"$\log_{{10}}$ {label}" if log_value else label


def _section_title(value: str) -> str:
    """Return a default pseudo-section title."""
    return f"TDEM pseudo-section: {_axis_label(value)}"
