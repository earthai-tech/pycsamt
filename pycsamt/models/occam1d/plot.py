# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""Validated, style-aware visualization of Occam1D inversion results.

All functions plot physical units: resistivity in ohm metres, phase in
degrees, frequency in hertz, and depth in metres. Matplotlib is imported
lazily so non-plotting Occam1D workflows do not require it at import time.
"""

from __future__ import annotations

import copy
import math

import numpy as np

from ...api.occam1d import resolve_occam1d_style
from ...api.station import PYCSAMT_STATION_RENDERING
from .base import Occam1DBase

__all__ = [
    "PlotConvergence",
    "PlotModel",
    "PlotResponse",
    "PlotSummary",
]

_RHO_CODE = 103
_PHASE_CODE = 104


def _positive_real(name: str, value, *, allow_none=False) -> float | None:
    """Validate a finite positive plotting limit."""
    if value is None and allow_none:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{name} must be a positive real number.")
    result = float(value)
    if not math.isfinite(result) or result <= 0:
        raise ValueError(f"{name} must be finite and strictly positive.")
    return result


def _figure_size(value) -> tuple[float, float] | None:
    """Return a validated two-element figure size in inches."""
    if value is None:
        return None
    try:
        width, height = value
    except (TypeError, ValueError):
        raise ValueError("figsize must contain exactly two values.") from None
    return (
        _positive_real("figsize width", width),
        _positive_real("figsize height", height),
    )


class _PlotBase(Occam1DBase):
    """Common validated state for Occam1D figure builders.

    Parameters
    ----------
    result : Occam1DResult
        Fully loaded result aggregate.
    figsize : (float, float), optional
        Figure width and height in inches.
    dpi : int, default=120
        Positive raster resolution used when constructing the figure.
    grid : bool, default=True
        Draw light scientific reference grids.
    title : str or None, optional
        Per-figure title override. ``None`` uses an informative default.
    """

    def __init__(
        self,
        result,
        figsize=None,
        dpi=120,
        *,
        grid=True,
        title=None,
        style=None,
        style_overrides=None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        required = ("data", "model", "resistivity", "iteration")
        missing = [name for name in required if not hasattr(result, name)]
        if missing:
            names = ", ".join(missing)
            raise TypeError(f"result is missing required attributes: {names}.")
        if not isinstance(dpi, int) or isinstance(dpi, bool) or dpi < 1:
            raise ValueError("dpi must be a positive integer.")
        if not isinstance(grid, bool):
            raise TypeError("grid must be a bool.")
        if title is not None and not isinstance(title, str):
            raise TypeError("title must be a string or None.")
        self.result = result
        self.figsize = _figure_size(figsize)
        self.dpi = dpi
        self.grid = grid
        self.title = title
        if style_overrides is None:
            style_overrides = {}
        if not isinstance(style_overrides, dict):
            raise TypeError("style_overrides must be a mapping or None.")
        self.style = resolve_occam1d_style(style, **style_overrides)

    @staticmethod
    def _pyplot():
        """Import and return pyplot only when a plot is requested."""
        try:
            import matplotlib.pyplot as plt
        except ImportError as error:
            raise ImportError(
                "Occam1D plotting requires matplotlib. Install matplotlib "
                "to create model-image products."
            ) from error
        return plt

    def _new_subplots(self, *args, **kwargs):
        """Construct a figure using validated common configuration."""
        plt = self._pyplot()
        kwargs.setdefault("dpi", self.dpi)
        if self.figsize is not None:
            kwargs.setdefault("figsize", self.figsize)
        return plt.subplots(*args, **kwargs)


def _finite_xy(x, *values, positive_x=False):
    """Return aligned finite arrays without turning missing data into zero."""
    arrays = [np.asarray(x, dtype=float)] + [
        np.asarray(value, dtype=float) for value in values
    ]
    sizes = {array.size for array in arrays}
    if len(sizes) != 1:
        raise ValueError(
            "Plot arrays must have identical one-dimensional size."
        )
    mask = np.logical_and.reduce([np.isfinite(array) for array in arrays])
    if positive_x:
        mask &= arrays[0] > 0
    return tuple(array[mask] for array in arrays)


def _model_coordinates(result, depth_max=None):
    """Return finite step coordinates in ohm metres and metres."""
    depth = np.asarray(result.model.depth, dtype=float)
    thickness = np.asarray(result.model.thickness, dtype=float)
    rho = np.asarray(result.resistivity, dtype=float)
    if not (depth.size == thickness.size == rho.size) or depth.size == 0:
        raise ValueError("Model depth, thickness, and resistivity must align.")
    if np.any(~np.isfinite(depth)) or np.any(depth < 0):
        raise ValueError("Model depths must be finite and non-negative.")
    if np.any(~np.isfinite(rho)) or np.any(rho <= 0):
        raise ValueError("Model resistivities must be finite and positive.")
    finite_thickness = thickness[np.isfinite(thickness) & (thickness > 0)]
    terminal = finite_thickness[-1] if finite_thickness.size else 1.0
    terminal = max(terminal, depth[-1] * 0.35, 1.0)
    bottom = np.r_[depth[1:], depth[-1] + terminal]
    limit = _positive_real("depth_max", depth_max, allow_none=True)
    return np.repeat(rho, 2), np.column_stack((depth, bottom)).ravel(), (
        limit or bottom[-1]
    )


def _apply_station_marker(ax, station: str) -> None:
    """Draw the established filled inversion station marker."""
    style = copy.deepcopy(
        PYCSAMT_STATION_RENDERING.style_for("inversion").marker
    )
    left, right = ax.get_xlim()
    center = math.sqrt(left * right)
    ax.scatter(
        [center],
        [0.0],
        transform=ax.transData,
        **style.kwargs(),
    )
    ax.annotate(
        station,
        (center, 0.0),
        xytext=(0, 9),
        textcoords="offset points",
        ha="center",
        va="bottom",
        fontsize=8,
        clip_on=False,
        zorder=style.zorder + 1,
    )


def _draw_model(
    ax, result, style, *, depth_max=None, show_station=True, grid=True
):
    """Render the selected physical layer model on ``ax``."""
    rho, depth, limit = _model_coordinates(result, depth_max)
    artist = style.model
    if artist.visible:
        ax.plot(rho, depth, **artist.kwargs())
    ax.set_xscale("log")
    ax.set_xlabel(r"Resistivity ($\Omega\,m$)")
    ax.set_ylabel("Depth (m)")
    ax.set_ylim(limit, 0.0)
    ax.grid(
        grid,
        which="both",
        alpha=style.grid_alpha,
        linewidth=style.grid_linewidth,
    )
    if show_station:
        _apply_station_marker(ax, str(result.data.station))
    if style.model_legend and artist.visible and artist.label:
        ax.legend(
            frameon=style.legend_frame,
            fontsize=style.legend_fontsize,
            loc=style.legend_location,
        )


def _response_rows(result, code):
    """Return finite frequency, values, and physical errors for one type."""
    if result.response is None:
        raise RuntimeError(
            "No response is available for the selected iteration."
        )
    response = result.response
    mask = response.type_code == code
    if not np.any(mask):
        return tuple(np.empty(0, dtype=float) for _ in range(4))
    indices = response.frequency_index[mask] - 1
    valid = (indices >= 0) & (indices < result.data.frequency.size)
    observed, modeled = response.physical_values()
    errors = response.physical_errors()
    return _finite_xy(
        result.data.frequency[indices[valid]],
        observed[mask][valid],
        modeled[mask][valid],
        errors[mask][valid],
        positive_x=True,
    )


def _draw_response_axis(
    ax, result, code, style, *, grid=True, legend=False
):
    """Render one response type, including a clear no-data annotation."""
    frequency, observed, modeled, error = _response_rows(result, code)
    observed_style = style.observed
    predicted_style = (
        style.predicted if code == _RHO_CODE else style.phase_predicted
    )
    ylabel = (
        r"Apparent resistivity ($\Omega\,m$)"
        if code == _RHO_CODE
        else r"Phase ($^\circ$)"
    )
    if not frequency.size:
        ax.text(
            0.5,
            0.5,
            f"No {ylabel.split(' (')[0].lower()} rows",
            transform=ax.transAxes,
            ha="center",
            va="center",
            color="0.4",
        )
    else:
        error = np.abs(error)
        if observed_style.visible:
            observed_kw = observed_style.kwargs()
            observed_kw["linewidth"] = style.errorbar_linewidth
            ax.errorbar(
                frequency,
                observed,
                yerr=error,
                capsize=style.errorbar_capsize,
                **observed_kw,
            )
        order = np.argsort(frequency)
        if predicted_style.visible:
            ax.plot(
                frequency[order],
                modeled[order],
                **predicted_style.kwargs(),
            )
        if code == _RHO_CODE:
            positive = (observed > 0) & (modeled > 0)
            if not np.all(positive):
                raise ValueError(
                    "Apparent-resistivity plot values must be positive."
                )
            ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_ylabel(ylabel)
    ax.grid(
        grid,
        which="both",
        alpha=style.grid_alpha,
        linewidth=style.grid_linewidth,
    )
    if legend and frequency.size:
        ax.legend(
            frameon=style.legend_frame,
            fontsize=style.legend_fontsize,
            loc=style.legend_location,
        )


def _draw_response(rho_ax, phase_ax, result, style, *, grid=True):
    """Render paired apparent-resistivity and phase axes."""
    _draw_response_axis(
        rho_ax,
        result,
        _RHO_CODE,
        style,
        grid=grid,
        legend=style.response_legend,
    )
    _draw_response_axis(
        phase_ax,
        result,
        _PHASE_CODE,
        style,
        grid=grid,
    )
    rho_ax.tick_params(labelbottom=False)
    phase_ax.set_xlabel("Frequency (Hz)")


def _draw_convergence(
    ax,
    result,
    style,
    *,
    target=None,
    show_roughness=True,
    grid=True,
):
    """Render finite convergence metrics and return an optional twin axis."""
    if result.log is None or not result.log.n_iterations:
        raise RuntimeError("No convergence history is available.")
    target = (
        result.iteration_data.target_misfit
        if target is None
        else _positive_real("target", target)
    )
    iteration, rms = _finite_xy(result.log.iterations, result.log.rms)
    if not iteration.size:
        raise RuntimeError(
            "Convergence history contains no finite RMS values."
        )
    if style.iteration.visible:
        ax.plot(iteration, rms, **style.iteration.kwargs())
    if style.target.visible:
        target_kw = style.target.kwargs(include_marker=False)
        if style.target.label:
            target_kw["label"] = f"{style.target.label} ({target:g})"
        ax.axhline(target, **target_kw)
    if style.current_iteration.visible:
        ax.axvline(
            result.iteration,
            **style.current_iteration.kwargs(include_marker=False),
        )
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Normalized RMS", color=style.iteration.color)
    ax.tick_params(axis="y", labelcolor=style.iteration.color)
    ax.grid(
        grid,
        alpha=style.grid_alpha,
        linewidth=style.grid_linewidth,
    )
    if style.convergence_legend:
        ax.legend(
            frameon=style.legend_frame,
            fontsize=style.legend_fontsize,
            loc=style.legend_location,
        )
    twin = None
    if show_roughness and style.roughness.visible:
        it_rough, roughness = _finite_xy(
            result.log.iterations,
            result.log.roughness,
        )
        if roughness.size:
            twin = ax.twinx()
            twin.plot(it_rough, roughness, **style.roughness.kwargs())
            if np.all(roughness > 0):
                twin.set_yscale("log")
            twin.set_ylabel("Model roughness", color=style.roughness.color)
            twin.tick_params(axis="y", labelcolor=style.roughness.color)
            if style.roughness_legend:
                twin.legend(
                    frameon=style.legend_frame,
                    fontsize=style.legend_fontsize,
                    loc=style.legend_location,
                )
    return twin


class PlotModel(_PlotBase):
    """Plot one physical stepwise resistivity model against depth."""

    def plot(self, *, depth_max=None, show_station=True):
        """Return a model-profile :class:`matplotlib.figure.Figure`."""
        if not isinstance(show_station, bool):
            raise TypeError("show_station must be a bool.")
        figure, axis = self._new_subplots(figsize=self.figsize or (5.2, 7.0))
        _draw_model(
            axis,
            self.result,
            self.style,
            depth_max=depth_max,
            show_station=show_station,
            grid=self.grid,
        )
        axis.set_title(self.title or self._default_title())
        figure.tight_layout()
        return figure

    def _default_title(self) -> str:
        rms = self.result.final_rms
        metric = f"RMS {rms:.3g}" if math.isfinite(rms) else "RMS unavailable"
        return (
            f"{self.result.data.station} — Occam1D model\n"
            f"iteration {self.result.iteration}, {metric}"
        )


class PlotResponse(_PlotBase):
    """Plot physical observed and modeled resistivity and phase responses."""

    def plot(self):
        """Return a two-panel response-fit figure."""
        if self.result.response is None:
            raise RuntimeError(
                "No response is available for the selected iteration."
            )
        figure, axes = self._new_subplots(
            2,
            1,
            figsize=self.figsize or (7.0, 7.0),
            sharex=True,
        )
        _draw_response(
            axes[0], axes[1], self.result, self.style, grid=self.grid
        )
        figure.suptitle(
            self.title
            or f"{self.result.data.station} — observed and modeled response"
        )
        figure.tight_layout()
        return figure


class PlotConvergence(_PlotBase):
    """Plot normalized RMS, target misfit, and optional model roughness."""

    def plot(self, *, target=None, show_roughness=True):
        """Return a convergence-history figure."""
        if not isinstance(show_roughness, bool):
            raise TypeError("show_roughness must be a bool.")
        figure, axis = self._new_subplots(
            figsize=self.figsize or (7.0, 4.5)
        )
        _draw_convergence(
            axis,
            self.result,
            self.style,
            target=target,
            show_roughness=show_roughness,
            grid=self.grid,
        )
        axis.set_title(
            self.title
            or f"{self.result.data.station} — Occam1D convergence"
        )
        figure.tight_layout()
        return figure


class PlotSummary(_PlotBase):
    """Combine model, response fit, and convergence in one main image."""

    def plot(self, *, depth_max=None, target=None):
        """Return the main four-panel Occam1D summary figure."""
        if self.result.response is None:
            raise RuntimeError("Summary plotting requires a response file.")
        if self.result.log is None or not self.result.log.n_iterations:
            raise RuntimeError("Summary plotting requires a convergence log.")
        plt = self._pyplot()
        figure = plt.figure(
            figsize=self.figsize or (12.0, 8.5),
            dpi=self.dpi,
            constrained_layout=True,
        )
        grid = figure.add_gridspec(2, 2, width_ratios=(0.85, 1.5))
        model_axis = figure.add_subplot(grid[:, 0])
        response_grid = grid[0, 1].subgridspec(2, 1, hspace=0.06)
        rho_axis = figure.add_subplot(response_grid[0])
        phase_axis = figure.add_subplot(response_grid[1], sharex=rho_axis)
        convergence_axis = figure.add_subplot(grid[1, 1])
        _draw_model(
            model_axis,
            self.result,
            self.style,
            depth_max=depth_max,
            show_station=True,
            grid=self.grid,
        )
        _draw_response(
            rho_axis,
            phase_axis,
            self.result,
            self.style,
            grid=self.grid,
        )
        _draw_convergence(
            convergence_axis,
            self.result,
            self.style,
            target=target,
            show_roughness=True,
            grid=self.grid,
        )
        model_axis.set_title("Resistivity model")
        rho_axis.set_title("Response fit")
        convergence_axis.set_title("Convergence")
        rms = self.result.final_rms
        metric = f"RMS {rms:.3g}" if math.isfinite(rms) else "RMS unavailable"
        figure.suptitle(
            self.title
            or (
                f"{self.result.data.station} — Occam1D inversion summary "
                f"(iteration {self.result.iteration}, {metric})"
            ),
            fontsize=14,
        )
        return figure
