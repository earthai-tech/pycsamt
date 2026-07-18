# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Profile-view and pseudosection API for pyCSAMT surveys."""

from __future__ import annotations

from dataclasses import replace
from typing import Any

import numpy as np

from ._backends import (
    require_matplotlib_pyplot,
    require_plotly,
    require_plotly_subplots,
)
from ._core import (
    MapData,
    ensure_map_data,
    pseudosection_table,
)
from .config import ProfileMapOptions
from .styles import theme_colors, to_plotly_cmap


class ProfileMap:
    """Builder object for profile maps."""

    def __init__(
        self,
        data: Any,
        *,
        options: ProfileMapOptions | None = None,
        **ensure_kwargs: Any,
    ) -> None:
        self.data: MapData = ensure_map_data(
            data,
            **ensure_kwargs,
        )
        self.options = options or ProfileMapOptions()

    def with_options(self, **kwargs: Any) -> ProfileMap:
        """Return a copy with updated options."""
        new = object.__new__(ProfileMap)
        new.data = self.data
        new.options = replace(self.options, **kwargs)
        return new

    def with_quantity(self, quantity: str) -> ProfileMap:
        """Return a copy using another mapped quantity."""
        return self.with_options(quantity=quantity)

    def with_component(self, component: str) -> ProfileMap:
        """Return a copy using one impedance component."""
        return self.with_options(
            component=component,
            components=(component,),
        )

    def with_components(
        self,
        *components: str,
    ) -> ProfileMap:
        """Return a copy using multiple components."""
        if not components:
            components = self.options.components
        return self.with_options(
            component=components[0],
            components=tuple(components),
        )

    def figure(self) -> Any:
        """Build and return the profile figure."""
        return build_profile_map(self.data, self.options)

    def pseudosection(self) -> Any:
        """Build and return the pseudosection figure."""
        return build_pseudosection(self.data, self.options)


def plot_profile_map(
    data: Any,
    *,
    options: ProfileMapOptions | None = None,
    **kwargs: Any,
) -> Any:
    """Build a profile-view map."""
    return ProfileMap(
        data,
        options=options,
        **kwargs,
    ).figure()


def plot_pseudosection(
    data: Any,
    *,
    options: ProfileMapOptions | None = None,
    **kwargs: Any,
) -> Any:
    """Build a resistivity or phase pseudosection."""
    return ProfileMap(
        data,
        options=options,
        **kwargs,
    ).pseudosection()


def build_profile_map(
    data: MapData,
    options: ProfileMapOptions,
) -> Any:
    """Build the concrete profile-view figure."""
    return build_pseudosection(data, options)


def build_pseudosection(
    data: MapData,
    options: ProfileMapOptions,
) -> Any:
    """Build the concrete pseudosection figure."""
    if options.backend == "matplotlib":
        if options.by_line:
            return _matplotlib_line_grid_pseudosection(data, options)
        return _matplotlib_pseudosection(data, options)
    if options.backend != "plotly":
        msg = f"Unknown profile backend: {options.backend}"
        raise ValueError(msg)
    if options.by_line:
        return _line_grid_pseudosection(data, options)
    if len(options.components) > 1:
        return _multi_component_pseudosection(data, options)
    return _single_pseudosection(
        data,
        options,
        component=options.component,
    )


def _single_pseudosection(
    data: MapData,
    options: ProfileMapOptions,
    *,
    component: str,
):
    go = require_plotly()

    quantity = _quantity_name(options.quantity)
    table = pseudosection_table(
        data,
        quantity=quantity,
        component=component,
    )
    table = _filter_table(table, options, quantity)
    colors = theme_colors(options.theme)
    if table.empty:
        return _empty_profile_figure(colors)
    if options.x_axis == "distance":
        x_col = "distance"
    else:
        x_col = "station"
    piv = table.pivot_table(
        index="period",
        columns=x_col,
        values="value",
        aggfunc="median",
    ).sort_index()
    z_values = _profile_values_for_plot(
        piv.to_numpy(dtype=float),
        quantity,
        options,
    )
    if quantity == "rho":
        if options.log_rho:
            colorbar = "log₁₀ ρ (Ω·m)"
        else:
            colorbar = "ρ (Ω·m)"
        cmap = options.cmap or "Jet"
    else:
        colorbar = "φ (°)"
        cmap = options.cmap or "RdBu_r"
    y_values = np.log10(piv.index.to_numpy(dtype=float))
    fig = go.Figure(
        go.Heatmap(
            z=z_values,
            x=list(piv.columns),
            y=y_values,
            colorscale=to_plotly_cmap(cmap),
            zmin=_zmin(options),
            zmax=_zmax(options),
            colorbar=dict(title=dict(text=colorbar, side="right")),
            hovertemplate=(
                "<b>%{x}</b><br>"
                "log10(T): %{y:.2f}<br>"
                "value: %{z:.2f}<extra></extra>"
            ),
        )
    )
    title = f"{quantity.upper()} {component.upper()} pseudosection"
    _style_profile(fig, colors, title)
    if options.x_axis == "distance":
        fig.update_xaxes(title_text="Distance (km)")
    return fig


def _line_grid_pseudosection(data: MapData, options: ProfileMapOptions):
    """Build one pseudosection panel per line, laid out in a grid.

    A merged pseudosection concatenates every line's stations onto one
    x-axis, which reads as a single continuous section even though the
    lines are unrelated traverses. This renders each line as its own
    small-multiple panel instead, sharing one color scale so panels
    stay comparable.
    """
    go = require_plotly()
    make_subplots = require_plotly_subplots()

    quantity = _quantity_name(options.quantity)
    component = options.component
    table = pseudosection_table(
        data,
        quantity=quantity,
        component=component,
    )
    table = _filter_table(table, options, quantity)
    colors = theme_colors(options.theme)
    if table.empty:
        return _empty_profile_figure(colors)

    lines = sorted(str(v) for v in table["line"].dropna().unique())
    if len(lines) <= 1:
        return _single_pseudosection(data, options, component=component)

    x_col = "distance" if options.x_axis == "distance" else "station"
    if quantity == "rho":
        colorbar = "log₁₀ ρ (Ω·m)" if options.log_rho else "ρ (Ω·m)"
        cmap = options.cmap or "Jet"
    else:
        colorbar = "φ (°)"
        cmap = options.cmap or "RdBu_r"

    panels: list[tuple[list, np.ndarray, np.ndarray] | None] = []
    finite_values = []
    for line in lines:
        piv = (
            table[table["line"] == line]
            .pivot_table(
                index="period",
                columns=x_col,
                values="value",
                aggfunc="median",
            )
            .sort_index()
        )
        if piv.empty:
            panels.append(None)
            continue
        z = _profile_values_for_plot(
            piv.to_numpy(dtype=float),
            quantity,
            options,
        )
        y = np.log10(piv.index.to_numpy(dtype=float))
        panels.append((list(piv.columns), y, z))
        finite = z[np.isfinite(z)]
        if finite.size:
            finite_values.append(finite)

    zmin, zmax = _zmin(options), _zmax(options)
    if (zmin is None or zmax is None) and finite_values:
        merged = np.concatenate(finite_values)
        zmin = float(merged.min()) if zmin is None else zmin
        zmax = float(merged.max()) if zmax is None else zmax

    cols = max(
        1,
        int(options.line_cols) if options.line_cols else min(3, len(lines)),
    )
    rows = -(-len(lines) // cols)

    fig = make_subplots(
        rows=rows,
        cols=cols,
        subplot_titles=lines,
        shared_yaxes=True,
        horizontal_spacing=0.06,
        vertical_spacing=0.12,
    )
    scale_shown = False
    last_row_for_col: dict[int, int] = {}
    for idx, panel in enumerate(panels):
        if panel is None:
            continue
        row, col = divmod(idx, cols)
        x, y, z = panel
        fig.add_trace(
            go.Heatmap(
                z=z,
                x=x,
                y=y,
                colorscale=to_plotly_cmap(cmap),
                zmin=zmin,
                zmax=zmax,
                showscale=not scale_shown,
                colorbar=dict(title=dict(text=colorbar, side="right")),
                hovertemplate=(
                    "<b>%{x}</b><br>"
                    "log10(T): %{y:.2f}<br>"
                    "value: %{z:.2f}<extra></extra>"
                ),
            ),
            row=row + 1,
            col=col + 1,
        )
        scale_shown = True
        last_row_for_col[col] = max(last_row_for_col.get(col, 0), row + 1)

    fig.update_yaxes(autorange="reversed")
    fig.update_yaxes(title_text="log10 period (s)", col=1)
    x_title = "Distance (km)" if x_col == "distance" else "Station"
    for col, row in last_row_for_col.items():
        fig.update_xaxes(title_text=x_title, row=row, col=col + 1)

    fig.update_layout(
        title=f"{quantity.upper()} {component.upper()} pseudosection by line",
        paper_bgcolor=colors["paper"],
        plot_bgcolor=colors["plot"],
        font=dict(color=colors["text"]),
        height=max(320, int(options.height_per_panel) * rows),
    )
    return fig


def _matplotlib_pseudosection(
    data: MapData,
    options: ProfileMapOptions,
):
    plt = require_matplotlib_pyplot()

    quantity = _quantity_name(options.quantity)
    components = options.components or (options.component,)
    colors = theme_colors(options.theme)
    fig, axes = plt.subplots(
        len(components),
        1,
        figsize=(8, max(3, 2.8 * len(components))),
        squeeze=False,
    )
    fig.patch.set_facecolor(colors["paper"])
    for ax, comp in zip(axes.ravel(), components):
        _draw_mpl_panel(
            ax,
            data,
            options,
            quantity,
            comp,
            colors,
        )
    fig.tight_layout()
    return fig


def _draw_mpl_panel(
    ax,
    data,
    options,
    quantity,
    comp,
    colors,
):
    table = pseudosection_table(
        data,
        quantity=quantity,
        component=comp,
    )
    table = _filter_table(table, options, quantity)
    ax.set_facecolor(colors["plot"])
    ax.tick_params(colors=colors["text"])
    if table.empty:
        ax.text(
            0.5,
            0.5,
            "No profile data available",
            ha="center",
            va="center",
            color=colors["text"],
            transform=ax.transAxes,
        )
        return
    if options.x_axis == "distance":
        x_col = "distance"
    else:
        x_col = "station"
    piv = table.pivot_table(
        index="period",
        columns=x_col,
        values="value",
        aggfunc="median",
    ).sort_index()
    z_values = _profile_values_for_plot(
        piv.to_numpy(dtype=float),
        quantity,
        options,
    )
    im = ax.imshow(
        z_values,
        aspect="auto",
        origin="upper",
        cmap=options.cmap or _default_profile_cmap(quantity),
        vmin=_zmin(options),
        vmax=_zmax(options),
    )
    ax.set_title(
        f"{quantity.upper()} {comp.upper()}",
        color=colors["text"],
    )
    ax.set_ylabel("Period index", color=colors["text"])
    ax.set_xlabel(x_col.title(), color=colors["text"])
    ax.figure.colorbar(im, ax=ax)


def _matplotlib_line_grid_pseudosection(
    data: MapData,
    options: ProfileMapOptions,
):
    plt = require_matplotlib_pyplot()

    quantity = _quantity_name(options.quantity)
    component = options.component
    table = pseudosection_table(
        data,
        quantity=quantity,
        component=component,
    )
    table = _filter_table(table, options, quantity)
    colors = theme_colors(options.theme)
    if table.empty:
        fig, ax = plt.subplots(figsize=(6, 4))
        fig.patch.set_facecolor(colors["paper"])
        ax.set_facecolor(colors["plot"])
        ax.text(
            0.5,
            0.5,
            "No profile data available",
            ha="center",
            va="center",
            color=colors["text"],
            transform=ax.transAxes,
        )
        return fig

    lines = sorted(str(v) for v in table["line"].dropna().unique())
    if len(lines) <= 1:
        return _matplotlib_pseudosection(
            data,
            replace(options, components=(component,)),
        )

    x_col = "distance" if options.x_axis == "distance" else "station"
    pivots: list[Any] = []
    finite_values = []
    for line in lines:
        piv = (
            table[table["line"] == line]
            .pivot_table(
                index="period",
                columns=x_col,
                values="value",
                aggfunc="median",
            )
            .sort_index()
        )
        pivots.append(piv)
        if piv.empty:
            continue
        z = _profile_values_for_plot(
            piv.to_numpy(dtype=float),
            quantity,
            options,
        )
        finite = z[np.isfinite(z)]
        if finite.size:
            finite_values.append(finite)

    vmin, vmax = _zmin(options), _zmax(options)
    if (vmin is None or vmax is None) and finite_values:
        merged = np.concatenate(finite_values)
        vmin = float(merged.min()) if vmin is None else vmin
        vmax = float(merged.max()) if vmax is None else vmax

    cols = max(
        1,
        int(options.line_cols) if options.line_cols else min(3, len(lines)),
    )
    rows = -(-len(lines) // cols)
    fig, axes = plt.subplots(
        rows,
        cols,
        figsize=(4.2 * cols, 3.0 * rows),
        squeeze=False,
    )
    fig.patch.set_facecolor(colors["paper"])
    for idx, (line, piv) in enumerate(zip(lines, pivots)):
        row, col = divmod(idx, cols)
        ax = axes[row][col]
        ax.set_facecolor(colors["plot"])
        ax.tick_params(colors=colors["text"])
        if piv.empty:
            ax.text(
                0.5,
                0.5,
                "No data",
                ha="center",
                va="center",
                color=colors["text"],
                transform=ax.transAxes,
            )
            continue
        z_values = _profile_values_for_plot(
            piv.to_numpy(dtype=float),
            quantity,
            options,
        )
        im = ax.imshow(
            z_values,
            aspect="auto",
            origin="upper",
            cmap=options.cmap or _default_profile_cmap(quantity),
            vmin=vmin,
            vmax=vmax,
        )
        ax.set_title(line, color=colors["text"])
        ax.set_xlabel(x_col.title(), color=colors["text"])
        fig.colorbar(im, ax=ax)
    for idx in range(len(lines), rows * cols):
        row, col = divmod(idx, cols)
        axes[row][col].axis("off")
    axes[0][0].set_ylabel("Period index", color=colors["text"])
    fig.suptitle(
        f"{quantity.upper()} {component.upper()} pseudosection by line",
        color=colors["text"],
    )
    fig.tight_layout()
    return fig


def _multi_component_pseudosection(
    data: MapData,
    options: ProfileMapOptions,
):
    make_subplots = require_plotly_subplots()

    components = options.components or (options.component,)
    fig = make_subplots(
        rows=len(components),
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.08,
        subplot_titles=[
            (f"{_quantity_name(options.quantity).upper()} {c.upper()}")
            for c in components
        ],
    )
    colors = theme_colors(options.theme)
    any_data = False
    for row, comp in enumerate(components, start=1):
        sub = _single_pseudosection(
            data,
            options,
            component=comp,
        )
        if not sub.data:
            continue
        trace = sub.data[0]
        fig.add_trace(trace, row=row, col=1)
        any_data = True
    if not any_data:
        return _empty_profile_figure(colors)
    _style_profile(
        fig,
        colors,
        "Multi-component pseudosection",
    )
    # ``_style_profile`` sets the x-axis title on row 1's axis, which
    # ``shared_xaxes`` hides beneath the next subplot's title. Move it
    # onto the bottom row, where the (only) visible x-axis ticks are.
    fig.update_layout(xaxis_title=None)
    x_title = "Distance (km)" if options.x_axis == "distance" else "Station"
    fig.update_xaxes(title_text=x_title, row=len(components), col=1)
    fig.update_layout(
        height=max(
            320,
            int(options.height_per_panel) * len(components),
        )
    )
    return fig


def _filter_table(table, options, quantity: str):
    if table.empty:
        return table
    out = table.copy()
    if options.period_range:
        lo, hi = options.period_range
        mask = (out["period"] >= lo) & (out["period"] <= hi)
        out = out[mask]
    if quantity == "phase" and options.phase_range:
        lo, hi = options.phase_range
        out = out[(out["value"] >= lo) & (out["value"] <= hi)]
    if options.value_range:
        lo, hi = options.value_range
        out = out[(out["value"] >= lo) & (out["value"] <= hi)]
    return out


def _profile_values_for_plot(values, quantity, options):
    if quantity == "rho" and options.log_rho:
        return np.where(
            values > 0,
            np.log10(values),
            np.nan,
        )
    return values


def _default_profile_cmap(quantity: str) -> str:
    return "jet" if quantity == "rho" else "RdBu_r"


def _zmin(options) -> float | None:
    if not options.value_range:
        return None
    lo, _ = options.value_range
    if options.quantity.lower() in {"rho", "resistivity"}:
        if options.log_rho and lo > 0:
            return float(np.log10(lo))
    return float(lo)


def _zmax(options) -> float | None:
    if not options.value_range:
        return None
    _, hi = options.value_range
    if options.quantity.lower() in {"rho", "resistivity"}:
        if options.log_rho and hi > 0:
            return float(np.log10(hi))
    return float(hi)


def _quantity_name(quantity: str) -> str:
    value = quantity.lower()
    if value in {"phase", "phi"}:
        return "phase"
    return "rho"


def _style_profile(fig, colors, title: str) -> None:
    fig.update_layout(
        title=title,
        xaxis_title="Station",
        yaxis_title="log10 period (s)",
        yaxis_autorange="reversed",
        paper_bgcolor=colors["paper"],
        plot_bgcolor=colors["plot"],
        font=dict(color=colors["text"]),
    )


def _empty_profile_figure(colors):
    go = require_plotly()

    fig = go.Figure()
    fig.add_annotation(
        text="No profile data available",
        x=0.5,
        y=0.5,
        xref="paper",
        yref="paper",
        showarrow=False,
        font=dict(color=colors["text"]),
    )
    fig.update_layout(
        paper_bgcolor=colors["paper"],
        plot_bgcolor=colors["plot"],
    )
    return fig
