# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""3-D map API for fence, block, and depth-slice views."""

from __future__ import annotations

from dataclasses import replace
from typing import Any

import numpy as np

from ._backends import require_plotly
from ._core import (
    MapData,
    ensure_map_data,
    pseudosection_table,
)
from .config import VolumeMapOptions
from .geometry import (
    normalize_offsets,
    resolve_offset,
    survey_uv,
)
from .styles import theme_colors, to_plotly_cmap


class Map3D:
    """Builder object for 3-D survey maps."""

    def __init__(
        self,
        data: Any,
        *,
        options: VolumeMapOptions | None = None,
        **ensure_kwargs: Any,
    ) -> None:
        self.data: MapData = ensure_map_data(
            data,
            **ensure_kwargs,
        )
        self.options = options or VolumeMapOptions()

    def with_mode(self, mode: str) -> Map3D:
        """Return a copy with another 3-D mode."""
        opts = replace(self.options, mode=mode)  # type: ignore[arg-type]
        new = object.__new__(Map3D)
        new.data = self.data
        new.options = opts
        return new

    def with_options(self, **kwargs: Any) -> Map3D:
        """Return a copy with updated options."""
        new = object.__new__(Map3D)
        new.data = self.data
        new.options = replace(self.options, **kwargs)
        return new

    def with_quantity(self, quantity: str) -> Map3D:
        """Return a copy using another mapped quantity."""
        return self.with_options(quantity=quantity)

    def with_component(self, component: str) -> Map3D:
        """Return a copy using another component."""
        return self.with_options(component=component)

    def figure(self) -> Any:
        """Build and return the 3-D map figure."""
        return build_3d_map(self.data, self.options)


VolumeMap = Map3D


def plot_3d_map(
    data: Any,
    *,
    options: VolumeMapOptions | None = None,
    **kwargs: Any,
) -> Any:
    """Build a 3-D map."""
    return Map3D(data, options=options, **kwargs).figure()


def plot_volume_map(
    data: Any,
    *,
    options: VolumeMapOptions | None = None,
    **kwargs: Any,
) -> Any:
    """Build a 3-D volume map."""
    return plot_3d_map(data, options=options, **kwargs)


def build_3d_map(
    data: MapData,
    options: VolumeMapOptions,
) -> Any:
    """Build the concrete 3-D figure."""
    profiles = _profile_grids(data, options)
    colors = theme_colors(options.theme)
    if not profiles:
        return _empty_3d_figure(colors)
    builders = {
        "block": _block_figure,
        "depth": _depth_figure,
        "surface": _surface_figure,
        "fence": _fence_figure,
    }
    builder = builders.get(options.mode)
    if builder is None:
        msg = f"Unknown 3-D map mode: {options.mode}"
        raise ValueError(msg)
    fig = builder(profiles, options, colors)
    if options.show_stations:
        _add_station_markers(fig, profiles, options)
    _apply_axis_units(fig, options)
    return fig


def _unit_scale(unit) -> float:
    """Metres -> display factor ('m' -> 1, 'km' -> 1e-3)."""
    return 0.001 if str(unit).lower() == "km" else 1.0


def _apply_axis_units(fig, options) -> None:
    """Rescale trace coordinates from metres to the chosen display units.

    Profile (x) and line-offset (y) use ``x_unit``; depth (z) uses
    ``depth_unit``.  The axes are independent, so x can stay in metres
    while depth is shown in km without squashing the scene.
    """
    sx = _unit_scale(getattr(options, "x_unit", "m"))
    sz = _unit_scale(getattr(options, "depth_unit", "m"))
    if sx == 1.0 and sz == 1.0:
        return
    for trace in fig.data:
        for attr, scale in (("x", sx), ("y", sx), ("z", sz)):
            values = getattr(trace, attr, None)
            if values is None:
                continue
            arr = np.asarray(values, dtype=float)
            if arr.size:
                setattr(trace, attr, arr * scale)


def _profile_grids(
    data: MapData,
    options: VolumeMapOptions,
) -> dict[str, dict[str, np.ndarray]]:
    """Per-line (x, z, rho, ...) grids for the 3-D builders.

    Merges two sources: lines with a precomputed section (e.g. a
    ModEM 3-D volume already sliced per line — see
    :mod:`pycsamt.map.inversion`) are used directly; any other line
    falls back to the EDI pseudosection/skin-depth path. Both share
    one real-geometry projection (:func:`_dataset_uv`) so mixed
    EDI + inversion multi-line views stay consistently placed.
    """
    uv = _dataset_uv(data)
    out = _edi_profile_grids(data, options, uv)
    out.update(_inversion_profile_grids(data, options, uv))
    return out


def _edi_profile_grids(
    data: MapData,
    options: VolumeMapOptions,
    uv: dict[str, tuple[float, float]],
) -> dict[str, dict[str, np.ndarray]]:
    quantity = _volume_quantity(options.quantity)
    value_table = pseudosection_table(
        data,
        quantity=quantity,
        component=options.component,
    )
    rho_table = pseudosection_table(
        data,
        quantity="rho",
        component=options.component,
    )
    if value_table.empty or rho_table.empty:
        return {}
    station_line = {
        station.id: station.line or "line" for station in data.stations
    }
    station_elev = {
        station.id: (
            float(station.elevation)
            if station.elevation is not None
            else np.nan
        )
        for station in data.stations
    }
    value_table = _prepare_volume_table(
        value_table,
        station_line,
        options,
    )
    rho_table = _prepare_volume_table(
        rho_table,
        station_line,
        options,
    )
    out: dict[str, dict[str, np.ndarray]] = {}
    for line, group in value_table.groupby("line"):
        rho_group = rho_table[rho_table["line"] == line]
        if rho_group.empty:
            continue
        value_piv = group.pivot_table(
            index="period",
            columns="station",
            values="value",
            aggfunc="median",
        ).sort_index()
        rho_piv = rho_group.pivot_table(
            index="period",
            columns="station",
            values="value",
            aggfunc="median",
        ).sort_index()
        rho_piv = rho_piv.reindex_like(value_piv)
        if value_piv.empty or rho_piv.empty:
            continue
        periods = value_piv.index.to_numpy(dtype=float)
        values = value_piv.to_numpy(dtype=float)
        rho = rho_piv.to_numpy(dtype=float)
        rho = np.where(rho > 0, rho, np.nan)
        median_rho = np.nanmedian(rho, axis=1)
        depth = 503.0 * np.sqrt(median_rho * periods)
        depth = np.where(np.isfinite(depth), depth, periods)
        value_full = values  # before any depth-range clip -- see below
        if options.depth_range:
            lo, hi = options.depth_range
            keep = (depth >= lo) & (depth <= hi)
            depth = depth[keep]
            rho = rho[keep, :]
            values = values[keep, :]
            periods = periods[keep]
        if values.size == 0:
            continue
        stations = np.array(
            list(value_piv.columns),
            dtype=object,
        )
        x, y_offset = _station_uv(group, stations, uv)
        out[str(line)] = {
            "x": x,
            "y_offset": y_offset,
            "z": depth,
            "period": periods,
            "rho": rho,
            "value": values,
            # Full (pre depth-clip) values so the colour scale stays put
            # while the user scrubs the depth window / slice depth --
            # same principle as ignoring rho_range for cmin/cmax.
            "value_full": value_full,
            "elev": np.array(
                [station_elev.get(str(s), np.nan) for s in stations],
                dtype=float,
            ),
            "quantity": np.array([quantity], dtype=object),
            "stations": stations,
        }
    return out


def _inversion_profile_grids(
    data: MapData,
    options: VolumeMapOptions,
    uv: dict[str, tuple[float, float]],
) -> dict[str, dict[str, np.ndarray]]:
    """Per-line grids for lines with a precomputed section.

    Consumes ``data.metadata["sections"]`` (see
    :func:`pycsamt.map.inversion.load_modem_lines`) directly instead
    of the EDI pseudosection/skin-depth path — an inversion result is
    already a real ``(x, z, rho)`` section, not raw impedance spectra
    to be gridded.
    """
    sections = (data.metadata or {}).get("sections")
    if not sections:
        return {}
    # A line survives (via merge/filter/mask) in ``metadata["sections"]``
    # even after its stations are removed from ``data.stations`` — skip
    # it here so a filtered-out or masked-away line doesn't keep
    # rendering its precomputed section.
    active_lines = {s.line for s in data.stations if s.line}
    # Elevation is read live from ``data.stations`` (falling back to the
    # section's own snapshot for stations no longer present there) so
    # MapView.with_elevations()/the topo-apply flow — which only update
    # station records, not the precomputed sections — still take effect.
    station_elev = {
        s.id: s.elevation for s in data.stations if s.elevation is not None
    }
    quantity = _volume_quantity(options.quantity)
    out: dict[str, dict[str, np.ndarray]] = {}
    for line, section in sections.items():
        if str(line) not in active_lines:
            continue
        stations = np.asarray(section["stations"], dtype=object)
        rho = np.asarray(section["rho"], dtype=float)
        z = np.asarray(section["z"], dtype=float)
        if stations.size == 0 or rho.size == 0:
            continue
        rho_full = rho  # before any depth-range clip
        if options.depth_range:
            lo, hi = options.depth_range
            keep = (z >= lo) & (z <= hi)
            z = z[keep]
            rho = rho[keep, :]
        if rho.size == 0:
            continue
        x, y_offset = _section_uv(stations, uv)
        section_elev = section.get("elev", np.full(stations.size, np.nan))
        elev = np.array(
            [
                station_elev.get(str(sid), section_elev[i])
                for i, sid in enumerate(stations)
            ],
            dtype=float,
        )
        out[str(line)] = {
            "x": x,
            "y_offset": y_offset,
            "z": z,
            "period": np.full_like(z, np.nan),
            "rho": rho,
            "value": rho,
            "value_full": rho_full,
            "elev": elev,
            "quantity": np.array([quantity], dtype=object),
            "stations": stations,
        }
    return out


def _fence_figure(profiles, options, colors):
    go = require_plotly()

    fig = go.Figure()
    az = np.deg2rad(float(options.azimuth))
    unit = _line_offset_unit(profiles)
    real_offsets = _line_real_offsets(profiles)
    cmin, cmax = _full_value_crange(profiles, options)
    for idx, (name, grid) in enumerate(profiles.items()):
        x, z_pos, values, elev = _prepare_section(grid, options)
        z = -z_pos
        xx, zz = np.meshgrid(x, z)
        offset = _line_offset(name, idx, real_offsets, unit, options)
        yy = np.full_like(
            xx,
            offset * np.cos(az),
        )
        xx = xx + offset * np.sin(az)
        if options.topography:
            zz = zz + elev[np.newaxis, :]
        # Open real holes in the mesh wherever the section is masked
        # (resistivity-range selection or genuine data gap), so the
        # panel shows only the selected cells instead of a full surface
        # tinted by a NaN colour.
        hole = ~np.isfinite(np.asarray(values, dtype=float))
        if hole.any():
            zz = np.where(hole, np.nan, zz)
        fig.add_trace(
            go.Surface(
                x=xx,
                y=yy,
                z=zz,
                surfacecolor=_color_values(values, options),
                colorscale=to_plotly_cmap(options.cmap),
                cmin=cmin,
                cmax=cmax,
                opacity=float(options.opacity),
                name=name,
                showscale=idx == 0,
                colorbar=dict(
                    title=dict(text=_colorbar_title(options), side="right")
                ),
                contours=_surface_contours(options),
            )
        )
        if options.show_terrain and options.topography:
            fig.add_trace(_terrain_trace(xx[0], yy[0], elev, name))
        if options.show_labels:
            fig.add_trace(
                go.Scatter3d(
                    x=[float(np.nanmean(xx))],
                    y=[float(np.nanmean(yy))],
                    z=[0.0],
                    text=[name],
                    mode="text",
                    showlegend=False,
                )
            )
    _style_3d(fig, options, colors)
    return fig


def _prepare_section(grid, options):
    """Return ``(x, depth, values, elev)`` for one line's fence panel.

    Coordinates are sorted/de-duplicated so the section is a proper
    monotonic grid, then optionally resampled onto a denser regular
    grid (cubic/linear spline) so the panel reads as a smooth section
    instead of a handful of raw per-station stripes — matching the
    interpolated volume the web app's 3-D map view renders.
    """
    x = np.asarray(grid["x"], dtype=float)
    z = np.asarray(grid["z"], dtype=float)
    values = _filtered_values(grid, options)
    # Cells the resistivity-range selection hides (options.rho_range).
    # Tracked separately from genuine NaN data gaps so the smoothing
    # pass below can re-open exactly these — see the re-mask step.
    hidden = _rho_range_mask(grid, options)
    elev = _elev_for(grid, options)

    order = np.argsort(x)
    x, values, elev = x[order], values[:, order], elev[order]
    hidden = hidden[:, order]
    x, uniq = np.unique(x, return_index=True)
    values, elev = values[:, uniq], elev[uniq]
    hidden = hidden[:, uniq]

    zorder = np.argsort(z)
    z, values = z[zorder], values[zorder, :]
    hidden = hidden[zorder, :]
    z, zuniq = np.unique(z, return_index=True)
    values = values[zuniq, :]
    hidden = hidden[zuniq, :]

    if not getattr(options, "smooth_sections", True):
        return x, z, values, elev
    if x.size < 3 or z.size < 3:
        return x, z, values, elev

    filled = _fill_nan_2d(z, x, values)
    if filled is None:
        return x, z, values, elev
    try:
        from scipy.interpolate import RectBivariateSpline

        res = max(int(getattr(options, "section_res", 100)), 2)
        kx = min(3, z.size - 1)
        ky = min(3, x.size - 1)
        spline = RectBivariateSpline(z, x, filled, kx=kx, ky=ky)
        xi = np.linspace(x.min(), x.max(), max(res, x.size))
        zi = np.linspace(z.min(), z.max(), max(res, z.size))
        vi = spline(zi, xi)
        # cubic splines can overshoot past the source range near sharp
        # gradients; clip back so the colorscale isn't stretched by
        # spline artefacts rather than real data.
        vi = np.clip(vi, float(np.nanmin(filled)), float(np.nanmax(filled)))
        # _fill_nan_2d filled the resistivity-range holes so the spline
        # could fit; punch them back through on the dense grid via the
        # nearest source cell. Without this the fence stays whole and
        # merely recolours when a range is selected, instead of showing
        # only the cells inside the range and masking the rest.
        if hidden.any():
            zi_idx = np.abs(zi[:, None] - z[None, :]).argmin(axis=1)
            xi_idx = np.abs(xi[:, None] - x[None, :]).argmin(axis=1)
            vi[hidden[np.ix_(zi_idx, xi_idx)]] = np.nan
        elev_i = np.interp(xi, x, elev)
        return xi, zi, vi, elev_i
    except Exception:  # noqa: BLE001 - fall back to the raw grid
        return x, z, values, elev


def _fill_nan_2d(z, x, values):
    """Fill NaN gaps with nearest-neighbour values so a spline can fit."""
    values = np.asarray(values, dtype=float)
    good = np.isfinite(values)
    if good.sum() < 3:
        return None
    if good.all():
        return values
    zz, xx = np.meshgrid(z, x, indexing="ij")
    try:
        from scipy.interpolate import griddata

        filled = griddata(
            (zz[good], xx[good]),
            values[good],
            (zz, xx),
            method="nearest",
        )
        return np.where(good, values, filled)
    except Exception:  # noqa: BLE001
        return np.where(good, values, float(np.nanmedian(values[good])))


def _volume_point_cloud(profiles, options):
    """Flatten the dense ``(x, y, z)`` volume grid to *finite-only*
    1-D point arrays for ``go.Volume``/``go.Isosurface``, plus an
    iso-value window clamped to the data actually present.

    ``go.Volume``/``go.Isosurface`` render **nothing** when the value
    array still carries the ``NaN`` cells left by the topography drape
    / no-data gaps, or when ``[isomin, isomax]`` lands outside the
    values that remain (e.g. a resistivity-band selection with no cell
    that conductive) or collapses to a single value. Both used to
    happen here, so the block silently disappeared under topography or
    any resistivity filter. This mirrors the working
    ``pycsamt.app.web`` block builder.

    Returns ``(x, y, z, value, iso_lo, iso_hi, cmin, cmax, banded)``
    where the grid is kept whole (``go.Volume`` needs a regular lattice
    to reconstruct a shape). ``cmin/cmax`` is always the full-model
    colour range (never the band). When ``banded`` (a resistivity-range
    selection is active), every out-of-band cell is set to a
    far-below-``iso_lo`` sentinel so ``go.Isosurface`` renders *only*
    the in-band body, each surviving cell keeping the exact hue it had
    on the unfiltered colourbar -- the 3-D equivalent of what fence mode
    already does. Otherwise only genuine no-data cells get the sentinel.
    Returns ``None`` when nothing can be drawn.
    """
    grid = _dense_volume_grid(profiles, options)
    if grid is None:
        return None
    x_arr, y_arr, z_arr, rho_vol = grid
    X, Y, Z = np.meshgrid(x_arr, y_arr, z_arr, indexing="ij")
    finite = rho_vol[np.isfinite(rho_vol)]
    if finite.size == 0:
        return None
    data_lo, data_hi = float(finite.min()), float(finite.max())

    # Colour scale = the *full, unfiltered* model range, exactly the
    # scale fence/depth mode use (:func:`_full_value_crange`). A
    # resistivity-range selection must only decide *which* cells show,
    # never rescale the colourbar -- the surviving cells keep the same
    # hue they had before the filter, other cells are masked away.
    cmin, cmax = _full_value_crange(profiles, options)
    if cmin is None or cmax is None:
        cmin, cmax = data_lo, data_hi
    if cmax <= cmin:  # a near-uniform model -- keep a valid colour span
        pad = max(abs(cmin) * 0.05, 0.05)
        cmin, cmax = cmin - pad, cmax + pad

    iso_lo, iso_hi = _iso_range(finite, options)
    # A resistivity band that doesn't overlap the model at all -> honour
    # it (draw nothing), rather than silently widen back to full range.
    if options.rho_range and (iso_lo > data_hi or iso_hi < data_lo):
        return None
    # go.Volume/go.Isosurface draw nothing when [isomin, isomax] lands
    # outside the values present or collapses to a point -- clamp it
    # into the data range and widen a degenerate window back to full.
    iso_lo = max(float(iso_lo), data_lo)
    iso_hi = min(float(iso_hi), data_hi)
    if iso_hi <= iso_lo:
        if data_hi > data_lo:
            iso_lo, iso_hi = data_lo, data_hi
        else:  # a perfectly uniform model -- pad so the window is valid
            pad = max(abs(data_lo) * 0.05, 0.05)
            iso_lo, iso_hi = data_lo - pad, data_hi + pad

    banded = bool(options.rho_range)

    span = iso_hi - iso_lo
    sentinel = iso_lo - max(span, 1.0) * (10.0 if banded else 1.0)
    if banded:
        # ``isomin``/``isomax`` alone do NOT isolate a value band -- the
        # isomax iso-surface still wraps everything *above* it, so a
        # "conductive" selection would drape the whole resistive model.
        # Push every out-of-band cell (both sides) to the sentinel so
        # only the band itself has a renderable value.
        in_band = (
            np.isfinite(rho_vol)
            & (rho_vol >= iso_lo)
            & (rho_vol <= iso_hi)
        )
        value = np.where(in_band, rho_vol, sentinel)
        # Structure smoothing: feather the sentinel <-> in-band cliff so
        # go.Isosurface cuts a rounded envelope, not a voxel staircase.
        # (Non-banded block/iso already round off from the denser
        # lattice + grid blur in :func:`_dense_volume_grid`.)
        _, _, _, sigma_cells = _volume_smoothing_params(
            float(getattr(options, "volume_smoothing", 0.0) or 0.0)
        )
        if sigma_cells > 0:
            value = _feather_volume(
                value, sentinel,
                _axis_sigma(value.shape, sigma_cells + 0.6),
            )
    else:
        value = np.where(np.isfinite(rho_vol), rho_vol, sentinel)
    return (
        X.ravel(),
        Y.ravel(),
        Z.ravel(),
        value.ravel(),
        iso_lo,
        iso_hi,
        cmin,
        cmax,
        banded,
    )


def _annotate_3d(fig, text, colors):
    """Drop a centred note on an otherwise-empty 3-D figure."""
    fig.add_annotation(
        text=text,
        x=0.5,
        y=0.5,
        xref="paper",
        yref="paper",
        showarrow=False,
        font=dict(color=colors["text"], size=12),
    )
    return fig


_VOLUME_EMPTY_MSG = (
    "No 3-D block for the current filters — every cell is outside the "
    "selected depth or resistivity range."
)


def _block_figure(profiles, options, colors):
    go = require_plotly()

    fig = go.Figure()
    cloud = _volume_point_cloud(profiles, options)
    if cloud is None:
        _annotate_3d(fig, _VOLUME_EMPTY_MSG, colors)
        _style_3d(fig, options, colors)
        return fig
    vx, vy, vz, vv, iso_lo, iso_hi, cmin, cmax, banded = cloud
    colorbar = dict(
        title=dict(text=_colorbar_title(options), side="right")
    )
    if banded:
        # A resistivity-band selection is a "show me only this zone"
        # request -- render it as one closed iso-surface body with
        # filled caps, not a translucent whole-model volume.
        trace = go.Isosurface(
            x=vx,
            y=vy,
            z=vz,
            value=vv,
            isomin=iso_lo,
            isomax=iso_hi,
            surface_count=2,
            opacity=float(options.opacity),
            caps=dict(x_show=True, y_show=True, z_show=True),
            colorscale=to_plotly_cmap(options.cmap),
            cmin=cmin,
            cmax=cmax,
            showscale=True,
            colorbar=colorbar,
        )
    else:
        trace = go.Volume(
            x=vx,
            y=vy,
            z=vz,
            value=vv,
            isomin=iso_lo,
            isomax=iso_hi,
            # Fade rather than hard-cut so the block still reads as one
            # solid shape instead of a jagged threshold.
            opacityscale=[
                [0.0, 0.0],
                [0.2, 0.3],
                [0.5, 0.7],
                [1.0, 1.0],
            ],
            opacity=float(options.opacity),
            surface_count=max(2, int(options.surface_count)),
            caps=dict(x_show=False, y_show=False, z_show=False),
            colorscale=to_plotly_cmap(options.cmap),
            cmin=cmin,
            cmax=cmax,
            showscale=True,
            colorbar=colorbar,
        )
    fig.add_trace(trace)
    _style_3d(fig, options, colors)
    return fig


def _depth_figure(profiles, options, colors):
    go = require_plotly()

    depths = _slice_depths(profiles, options)
    az = np.deg2rad(float(options.azimuth))
    unit = _line_offset_unit(profiles)
    real_offsets = _line_real_offsets(profiles)
    fig = go.Figure()
    # Lines may hold different station counts; pad every row
    # to the widest so the per-depth surface grid is
    # rectangular. NaN gaps render as holes in Plotly.
    width = max(
        (np.asarray(g["x"]).size for g in profiles.values()),
        default=0,
    )
    if width == 0:
        return _empty_3d_figure(colors)
    cmin, cmax = _full_value_crange(profiles, options)

    # go.Surface interpolates between *adjacent* rows of the stacked
    # grid below -- lines must be vstacked in real cross-strike
    # (y-offset) order, not whatever order the profiles dict happens
    # to iterate in, or the reconstructed surface connects unrelated
    # lines and twists on itself. _line_offset's own idx argument (the
    # original insertion order) still has to stay unchanged, since
    # it's also the fallback spacing multiplier shared with
    # _add_station_markers/_dense_volume_grid -- only the row-stacking
    # order below is resorted.
    names = list(profiles.keys())
    offsets = [
        _line_offset(name, idx, real_offsets, unit, options)
        for idx, name in enumerate(names)
    ]
    stack_order = sorted(range(len(names)), key=lambda i: offsets[i])
    # Each line's own stations must be in real along-profile (x) order
    # too -- go.Surface (and the terrain polyline below) connects
    # *adjacent* columns, and grid["x"]'s column order only reflects
    # whatever order the source pivot table produced (e.g. lexical by
    # station name), not necessarily real spatial order. x itself is
    # already a correct per-station geometric distance regardless of
    # that order (see _station_uv/_station_x); this just makes the
    # column order agree with it, the same way _prepare_section
    # (fence) and _dense_volume_grid (block/surface) already do.
    x_orders = {
        i: np.argsort(np.asarray(profiles[names[i]]["x"], dtype=float))
        for i in range(len(names))
    }

    drawn_any = False
    # Deepest first so the shallower slices render on top (the usual
    # oblique / from-above viewing angle).
    for depth in sorted(depths, reverse=True):
        x_rows = []
        y_rows = []
        z_rows = []
        val_rows = []
        for i in stack_order:
            name = names[i]
            grid = profiles[name]
            offset = offsets[i]
            x_order = x_orders[i]
            values = _values_at_depth(
                grid,
                depth,
                options,
            )[x_order]
            x = np.asarray(grid["x"], dtype=float)[x_order]
            y = np.full_like(
                x,
                offset * np.cos(az),
            )
            x_rows.append(_pad_row(x + offset * np.sin(az), width))
            y_rows.append(_pad_row(y, width))
            z = _elev_for(grid, options)[x_order] - float(depth)
            z_rows.append(_pad_row(z, width))
            color = _color_values(values, options)
            val_rows.append(_pad_row(color, width))
        color_grid = np.vstack(val_rows)
        if not np.isfinite(color_grid).any():
            # This depth is outside every line's sampled range -- an
            # all-NaN ghost surface, skip it.
            continue
        fig.add_trace(
            go.Surface(
                x=np.vstack(x_rows),
                y=np.vstack(y_rows),
                z=np.vstack(z_rows),
                surfacecolor=color_grid,
                colorscale=to_plotly_cmap(options.cmap),
                cmin=cmin,
                cmax=cmax,
                opacity=float(options.opacity),
                showscale=not drawn_any,
                colorbar=dict(
                    title=dict(text=_colorbar_title(options), side="right")
                ),
                contours=_surface_contours(options),
                name=f"{depth:.0f} m",
                hovertemplate=(
                    f"depth {depth:.0f} m<br>"
                    "%{surfacecolor:.2f}<extra></extra>"
                ),
            )
        )
        drawn_any = True
    if not drawn_any:
        _annotate_3d(fig, _VOLUME_EMPTY_MSG, colors)
        _style_3d(fig, options, colors)
        return fig
    if options.show_terrain and options.topography:
        for i in stack_order:
            name = names[i]
            grid = profiles[name]
            offset = offsets[i]
            x_order = x_orders[i]
            x = np.asarray(grid["x"], dtype=float)[x_order]
            fig.add_trace(
                _terrain_trace(
                    x + offset * np.sin(az),
                    np.full_like(x, offset * np.cos(az)),
                    _elev_for(grid, options)[x_order],
                    name,
                )
            )
    _style_3d(fig, options, colors)
    return fig


def _surface_figure(profiles, options, colors):
    go = require_plotly()

    fig = go.Figure()
    cloud = _volume_point_cloud(profiles, options)
    if cloud is None:
        _annotate_3d(fig, _VOLUME_EMPTY_MSG, colors)
        _style_3d(fig, options, colors)
        return fig
    vx, vy, vz, vv, iso_lo, iso_hi, cmin, cmax, banded = cloud
    fig.add_trace(
        go.Isosurface(
            x=vx,
            y=vy,
            z=vz,
            value=vv,
            isomin=iso_lo,
            isomax=iso_hi,
            surface_count=2 if banded else max(2, int(options.surface_count)),
            opacity=float(options.opacity),
            colorscale=to_plotly_cmap(options.cmap),
            caps=dict(
                x_show=banded, y_show=banded, z_show=banded
            ),
            colorbar=dict(
                title=dict(text=_colorbar_title(options), side="right")
            ),
            cmin=cmin,
            cmax=cmax,
        )
    )
    _style_3d(fig, options, colors)
    return fig


def _prepare_volume_table(table, station_line, options):
    out = table.copy()
    out["line"] = out["station"].map(station_line).fillna("line")
    if options.period_range:
        lo, hi = options.period_range
        keep = (out["period"] >= lo) & (out["period"] <= hi)
        out = out[keep]
    return out


def _station_x(group, stations):
    x_map = group.groupby("station")["distance"].median().to_dict()
    x = np.array(
        [x_map.get(str(station), np.nan) for station in stations],
        dtype=float,
    )
    if not np.isfinite(x).all():
        # No real geometry: fall back to an index spaced ~100 m apart so
        # the profile axis stays comparable to depth (metres), not 0,1,2…
        return np.arange(len(stations), dtype=float) * 100.0
    # ``distance`` is in km; the 3-D depth axis is in metres — convert so
    # the profile (x) and depth (z) axes share a unit and render to scale.
    return x * 1000.0


def _line_offset_unit(profiles) -> float:
    """Base Y-offset unit for fence/depth lines, in the same (metre)
    unit as the profile axis, so lines separate proportionally.

    Only used as a fallback when stations carry no real coordinates —
    see :func:`_line_real_offsets` for the geometry-based placement.
    """
    spans = []
    for grid in profiles.values():
        x = np.asarray(grid["x"], dtype=float)
        x = x[np.isfinite(x)]
        if x.size > 1:
            spans.append(float(x.max() - x.min()))
    return max(spans) if spans else 1000.0


def _dataset_uv(data: MapData) -> dict[str, tuple[float, float]]:
    """Adapt :class:`MapData` stations to :func:`geometry.survey_uv`.

    See :mod:`pycsamt.map.geometry` for the shared projection math
    (also used by ``pycsamt.app.web``'s 3-D map) — this just extracts
    the plain id/lat/lon/line arrays it needs from ``data.stations``.
    """
    ids: list[str] = []
    lats: list[float] = []
    lons: list[float] = []
    lines: list[str] = []
    for station in data.stations:
        if station.latitude is None or station.longitude is None:
            continue
        ids.append(station.id)
        lats.append(float(station.latitude))
        lons.append(float(station.longitude))
        lines.append(station.line or "line")
    return survey_uv(ids, lats, lons, lines)


def _station_uv(group, stations, uv):
    """Return ``(x, y_offset)`` for one line.

    ``x`` is the per-station along-strike distance (metres) and
    ``y_offset`` is the line's cross-strike position (metres), both
    derived from :func:`_dataset_uv` when every station in the line
    has real coordinates. Falls back to the path-length/index scheme
    in :func:`_station_x` (with ``y_offset=nan``, signalling callers
    to use the synthetic per-panel spacing) otherwise.
    """
    ids = [str(s) for s in stations]
    u = np.array(
        [uv[i][0] if i in uv else np.nan for i in ids],
        dtype=float,
    )
    v = np.array(
        [uv[i][1] if i in uv else np.nan for i in ids],
        dtype=float,
    )
    if u.size and np.isfinite(u).all():
        return u, float(np.nanmedian(v))
    return _station_x(group, stations), float("nan")


def _section_uv(stations, uv):
    """Same as :func:`_station_uv`, for a precomputed section with no
    pseudosection ``group`` to fall back on — synthetic index spacing
    is used instead when coordinates are unavailable."""
    ids = [str(s) for s in stations]
    u = np.array(
        [uv[i][0] if i in uv else np.nan for i in ids],
        dtype=float,
    )
    v = np.array(
        [uv[i][1] if i in uv else np.nan for i in ids],
        dtype=float,
    )
    if u.size and np.isfinite(u).all():
        return u, float(np.nanmedian(v))
    return np.arange(len(ids), dtype=float) * 100.0, float("nan")


def _line_real_offsets(profiles) -> dict[str, float] | None:
    """Return ``{line: cross-strike offset (m)}`` when every line has
    real geometry, else ``None`` so callers fall back to synthetic,
    index-based spacing."""
    return normalize_offsets(
        {name: grid.get("y_offset") for name, grid in profiles.items()}
    )


def _line_offset(name, idx, real_offsets, unit, options) -> float:
    """Cross-strike placement for one line, in metres — see
    :func:`geometry.resolve_offset`."""
    return resolve_offset(
        name, idx, real_offsets, unit, float(options.line_spacing)
    )


def _rho_range_mask(grid, options):
    """Boolean grid, ``True`` where a cell falls outside
    ``options.rho_range`` and must be hidden. All-``False`` when no
    range is set."""
    value = np.asarray(grid["value"], dtype=float)
    if not options.rho_range:
        return np.zeros(value.shape, dtype=bool)
    lo, hi = options.rho_range
    rho = np.asarray(grid["rho"], dtype=float)
    return (rho < lo) | (rho > hi)


def _filtered_values(grid, options):
    out = np.asarray(grid["value"], dtype=float).copy()
    out[_rho_range_mask(grid, options)] = np.nan
    return out


def _color_values(values, options):
    arr = np.asarray(values, dtype=float)
    if _volume_quantity(options.quantity) == "rho" and options.log_color:
        return np.where(arr > 0, np.log10(arr), np.nan)
    return arr


def _full_value_crange(profiles, options) -> tuple[float | None, float | None]:
    """``(cmin, cmax)`` for the fence/depth-slice colour scale.

    Returns the user's own fixed ``options.value_range`` when set
    (via :func:`_crange`), else falls back to every profile's *full,
    unfiltered* value array, in the same colour space
    :func:`_color_values` renders (log10 rho when requested) --
    deliberately independent of ``options.rho_range``. The
    resistivity-range filter only decides which cells
    :func:`_filtered_values` masks to NaN (what is shown); it must
    never also renormalise cmin/cmax (the scale it is shown on), or a
    narrow selection (e.g. a "conductive" 1-100 ohm.m slice) gets
    stretched across the whole colourscale and reads as one near-solid
    colour instead of the same hues it had before filtering.
    """
    cmin, cmax = _crange(options)
    if cmin is not None and cmax is not None:
        return cmin, cmax
    finite_chunks = []
    for grid in profiles.values():
        # Prefer the pre-depth-clip values (``value_full``) so scrubbing
        # the depth window / slice depth never rescales the colourbar.
        raw = np.asarray(grid.get("value_full", grid["value"]), dtype=float)
        colored = _color_values(raw, options)
        finite = colored[np.isfinite(colored)]
        if finite.size:
            finite_chunks.append(finite)
    if not finite_chunks:
        return None, None
    combined = np.concatenate(finite_chunks)
    return float(combined.min()), float(combined.max())


def _crange(options) -> tuple[float | None, float | None]:
    """Return (cmin, cmax) in color space for the colorbar range."""
    if not options.value_range:
        return None, None
    lo, hi = options.value_range
    if _volume_quantity(options.quantity) == "rho" and options.log_color:
        lo = float(np.log10(lo)) if lo and lo > 0 else None
        hi = float(np.log10(hi)) if hi and hi > 0 else None
        return lo, hi
    return float(lo), float(hi)


def _elev_for(grid, options):
    """Return per-station surface elevation, or zeros when topo is off."""
    x = np.asarray(grid["x"], dtype=float)
    elev = np.asarray(grid.get("elev", []), dtype=float)
    if not options.topography or elev.size != x.size:
        return np.zeros_like(x)
    return np.where(np.isfinite(elev), elev, 0.0)


def _depth_interp_column(z, col_values, depth):
    """Interpolate one station column at ``depth``; NaN outside its
    sampled range or when it has < 1 finite sample."""
    good = np.isfinite(z) & np.isfinite(col_values)
    if not good.any():
        return np.nan
    return float(
        np.interp(
            float(depth),
            z[good],
            col_values[good],
            left=np.nan,
            right=np.nan,
        )
    )


def _values_at_depth(grid, depth, options):
    """Resistivity (or phase) on a horizontal slice at ``depth``.

    The resistivity-range filter is applied *after* the depth
    interpolation, against the resistivity actually interpolated at
    that depth -- not before it. Filtering first (NaN-ing whole cells)
    let ``np.interp`` bridge across an out-of-band layer and paint the
    slice a solid colour where it should be masked.
    """
    z = np.asarray(grid["z"], dtype=float)
    values = np.asarray(grid["value"], dtype=float)
    rho = np.asarray(grid["rho"], dtype=float)
    order = np.argsort(z)
    z = z[order]
    values = values[order, :]
    rho = rho[order, :]

    out = np.array(
        [
            _depth_interp_column(z, values[:, col], depth)
            for col in range(values.shape[1])
        ],
        dtype=float,
    )
    if options.rho_range:
        lo, hi = options.rho_range
        rho_here = np.array(
            [
                _depth_interp_column(z, rho[:, col], depth)
                for col in range(rho.shape[1])
            ],
            dtype=float,
        )
        drop = ~np.isfinite(rho_here) | (rho_here < lo) | (rho_here > hi)
        out[drop] = np.nan
    return out


def _marker_z_offset(profiles) -> float:
    """Small upward nudge (1% of the survey's own depth span) so
    station markers clear the panel surface instead of z-fighting it."""
    spans = []
    for grid in profiles.values():
        z = np.asarray(grid.get("z", []), dtype=float)
        z = z[np.isfinite(z)]
        if z.size:
            spans.append(float(z.max() - z.min()))
    span = max(spans) if spans else 0.0
    return max(span * 0.01, 1.0) if span > 0 else 1.0


def _ref_line_x_range(profiles) -> tuple[float, float]:
    """``(x_min, x_max)`` of the line :func:`_dense_volume_grid` uses as
    its shared x-reference (the one with the most stations) -- the exact
    x-extent the block / iso-surface is reconstructed over."""
    if not profiles:
        return (-np.inf, np.inf)
    ref = max(profiles.values(), key=lambda g: np.asarray(g["x"]).size)
    x = np.asarray(ref["x"], dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return (-np.inf, np.inf)
    return (float(x.min()), float(x.max()))


def _thin_indices(n: int, max_count: int | None) -> np.ndarray:
    """Evenly-spaced station indices, thinning ``n`` items down to at
    most ``max_count`` while always keeping the first and last one.

    Returns every index unchanged (``0..n-1``) when ``max_count`` is
    ``None`` or already ``>= n`` -- the default, so existing callers
    that never set :attr:`VolumeMapOptions.max_stations` keep showing
    every station exactly as before this option existed.
    """
    if max_count is None or n == 0 or max_count >= n:
        return np.arange(n)
    if max_count <= 1:
        return np.array([0], dtype=int)
    return np.unique(np.round(np.linspace(0, n - 1, max_count)).astype(int))


def _label_keep_set(drawn, stations, options) -> set:
    """Which of a line's *drawn* marker indices (``drawn``) also get a
    text label.

    ``station_label_names`` wins: label only those station ids.
    Otherwise keep an evenly-spaced ``station_label_fraction`` of the
    drawn markers (first and last always kept). Markers themselves are
    never affected -- this only thins the labels so a dense line stays
    legible.
    """
    if not drawn:
        return set()
    names = getattr(options, "station_label_names", None)
    if names:
        want = {str(n).strip() for n in names if str(n).strip()}
        return {j for j in drawn if str(stations[j]) in want}
    frac = float(getattr(options, "station_label_fraction", 1.0) or 1.0)
    if frac >= 1.0 or len(drawn) <= 2:
        return set(drawn)
    count = max(2, int(np.ceil(len(drawn) * max(frac, 0.0))))
    picks = _thin_indices(len(drawn), count)
    return {drawn[p] for p in picks}


# Plotly's Scatter3d marker symbol enum has no triangle variant at all
# (only circle/circle-open/cross/diamond/diamond-open/square/square-open/x
# — triangles exist only for 2-D go.Scatter). These two values are
# pycsamt's own extension: _add_station_markers detects them and draws
# real custom geometry (a Mesh3d for the filled case, a closed Scatter3d
# line loop for the open/outline case) instead of a sprite marker.
TRIANGLE_DOWN_SYMBOLS = ("triangle-down", "triangle-down-open")


def _triangle_marker_size(profiles, options) -> tuple[float, float]:
    """``(half_width, height)`` for a custom triangle-down glyph, in
    world (data) units — unlike a sprite marker's pixel ``size``, this
    scales with the survey's own geometry so the glyph reads the same
    relative size regardless of survey scale. ``options.station_size``
    still controls it (as a multiplier against its own default of 4),
    so the one existing size control keeps working for this symbol too.

    Width and height are deliberately scaled off *different* axes,
    each against a quantity that stays roughly constant as a survey
    grows: width against the median along-profile station spacing (not
    the total profile length, which can span anywhere from tens of
    metres to kilometres depending on station count with no change in
    how "big" one marker should look), height against the depth-axis
    span the same way :func:`_marker_z_offset` already does. Coupling
    height to the width instead (e.g. a fixed aspect ratio) breaks
    badly whenever the along-profile and depth extents differ by
    orders of magnitude, which is the common case, not an edge case.
    """
    spacings = []
    z_spans = []
    for grid in profiles.values():
        x = np.sort(np.asarray(grid.get("x", []), dtype=float))
        x = x[np.isfinite(x)]
        if x.size > 1:
            spacings.append(float(np.median(np.diff(x))))
        z = np.asarray(grid.get("z", []), dtype=float)
        z = z[np.isfinite(z)]
        if z.size:
            z_spans.append(float(z.max() - z.min()))
    spacing = abs(float(np.median(spacings))) if spacings else 100.0
    z_span = max(z_spans) if z_spans else spacing
    size_frac = max(float(options.station_size), 0.1) / 4.0
    # Capped independently of size_frac: a user cranking station_size
    # way up should get a visibly bigger marker, never one so big it
    # swallows the gap to the next station or dwarfs the whole panel.
    half_width = min(max(spacing * 0.15 * size_frac, 1.0), spacing * 0.45)
    # Height is a *small* glyph, not a structural element: keep it near
    # the along-profile spacing and hard-cap it at a few percent of the
    # depth span so a marker never floats visibly clear of the section /
    # block it is meant to sit on.
    height = min(
        max(z_span * 0.02 * size_frac, spacing * 0.35, 1.0),
        z_span * 0.06,
    )
    return half_width, height


def _triangle_down_mesh_trace(
    xs, ys, zs, half_width: float, height: float, color: str, hover: list[str]
):
    """One combined ``Mesh3d`` trace: a filled, downward-pointing
    triangle per station (apex at the station's own point, base
    extending ``height`` upward in +Z, ``2 * half_width`` wide in X) —
    lies in the vertical X-Z plane at each station's Y, the orientation
    that stays least foreshortened under Plotly's default oblique
    camera (an isometric-style eye vector with x == y == z)."""
    go = require_plotly()
    vx: list[float] = []
    vy: list[float] = []
    vz: list[float] = []
    vi: list[int] = []
    vj: list[int] = []
    vk: list[int] = []
    vhover: list[str] = []
    for x, y, z, h in zip(xs, ys, zs, hover):
        base = len(vx)
        vx.extend([x, x - half_width, x + half_width])
        vy.extend([y, y, y])
        vz.extend([z, z + height, z + height])
        vhover.extend([h, h, h])
        vi.append(base)
        vj.append(base + 1)
        vk.append(base + 2)
    return go.Mesh3d(
        x=vx,
        y=vy,
        z=vz,
        i=vi,
        j=vj,
        k=vk,
        color=color,
        opacity=1.0,
        flatshading=True,
        hovertext=vhover,
        hovertemplate="%{hovertext}<extra></extra>",
        name="stations",
        showlegend=False,
    )


def _triangle_down_outline_trace(
    xs, ys, zs, half_width: float, height: float, color: str, hover: list[str]
):
    """One combined ``Scatter3d`` line trace: the unfilled/outline
    variant of :func:`_triangle_down_mesh_trace` — each triangle drawn
    as a closed 4-point loop, with a ``None`` gap between triangles so
    they render as separate closed shapes in a single trace."""
    go = require_plotly()
    lx: list[float | None] = []
    ly: list[float | None] = []
    lz: list[float | None] = []
    lhover: list[str | None] = []
    for x, y, z, h in zip(xs, ys, zs, hover):
        loop_x = [x, x - half_width, x + half_width, x]
        loop_z = [z, z + height, z + height, z]
        lx.extend(loop_x + [None])
        ly.extend([y] * 4 + [None])
        lz.extend(loop_z + [None])
        lhover.extend([h] * 4 + [None])
    return go.Scatter3d(
        x=lx,
        y=ly,
        z=lz,
        mode="lines",
        line=dict(color=color, width=3),
        hovertext=lhover,
        hovertemplate="%{hovertext}<extra></extra>",
        name="stations",
        showlegend=False,
    )


#: Modes whose grid is built from :func:`_dense_volume_grid` -- one
#: literal, axis-aligned ``(x, y, z)`` lattice shared by both
#: ``go.Volume`` (block) and ``go.Isosurface`` (surface), neither of
#: which can be rotated by azimuth without shearing the reconstructed
#: shape (see ``test_block_mode_ignores_azimuth``). Markers must sit
#: in that same unrotated frame for both, not just block.
_AXIS_ALIGNED_MODES = {"block", "surface"}


def _add_station_markers(fig, profiles, options) -> None:
    """Overlay per-station markers at the survey surface across lines.

    Block and surface (isosurface) mode share one dense ``(x, y, z)``
    grid that is built axis-aligned and deliberately never rotated by
    azimuth -- ``go.Volume``/``go.Isosurface`` need one literal
    rectangular lattice to reconstruct a shape, and rotating it would
    shear that reconstruction. Markers must sit in that same unrotated
    frame for those two modes, or they drift outside the block/
    isosurface's own footprint whenever ``azimuth != 0`` -- fence/depth
    mode keep applying azimuth as before, since those panels rotate
    with it.
    """
    go = require_plotly()

    az = (
        0.0
        if options.mode in _AXIS_ALIGNED_MODES
        else np.deg2rad(float(options.azimuth))
    )
    unit = _line_offset_unit(profiles)
    real_offsets = _line_real_offsets(profiles)
    # A marker at exactly the same (x, y, z) as the panel's own top
    # edge z-fights with that surface in WebGL — same coordinates, so
    # the two triangles flicker/occlude depending on view angle. Lift
    # markers a hair above the surface, scaled to the survey's own
    # depth extent so it reads the same regardless of scale.
    z_offset = _marker_z_offset(profiles)
    # Block / iso-surface are reconstructed on ONE shared x-axis -- the
    # densest line's (see _dense_volume_grid). A station on a longer line
    # would otherwise plot its marker well outside the rendered block.
    x_lo, x_hi = (
        _ref_line_x_range(profiles)
        if options.mode in _AXIS_ALIGNED_MODES
        else (-np.inf, np.inf)
    )
    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []
    labels: list[str | None] = []
    hover: list[str] = []
    for line_idx, (name, grid) in enumerate(profiles.items()):
        x = np.asarray(grid["x"], dtype=float)
        elev = _elev_for(grid, options)
        offset = _line_offset(name, line_idx, real_offsets, unit, options)
        stations = grid.get("stations", np.arange(x.size))
        drawn = [
            j
            for j in _thin_indices(x.size, options.max_stations)
            if x_lo <= x[j] <= x_hi
        ]
        label_keep = _label_keep_set(drawn, stations, options)
        for j in drawn:
            xs.append(float(x[j] + offset * np.sin(az)))
            ys.append(float(offset * np.cos(az)))
            zs.append(float(elev[j]) + z_offset)
            labels.append(str(stations[j]) if j in label_keep else None)
            hover.append(f"{name} · {stations[j]}")
    if not xs:
        return
    show_labels = bool(options.station_labels)
    symbol = str(options.station_symbol)

    if symbol in TRIANGLE_DOWN_SYMBOLS:
        half_width, height = _triangle_marker_size(profiles, options)
        color = str(options.station_color)
        if symbol == "triangle-down":
            fig.add_trace(
                _triangle_down_mesh_trace(
                    xs, ys, zs, half_width, height, color, hover
                )
            )
        else:
            fig.add_trace(
                _triangle_down_outline_trace(
                    xs, ys, zs, half_width, height, color, hover
                )
            )
        if show_labels:
            _add_station_labels(
                fig, xs, ys, zs, labels, options,
                z_lift=height, clear_px=10,
            )
        return

    fig.add_trace(
        go.Scatter3d(
            x=xs,
            y=ys,
            z=zs,
            mode="markers",
            marker=dict(
                symbol=symbol,
                size=int(options.station_size),
                color=str(options.station_color),
                line=dict(width=0),
            ),
            hovertext=hover,
            hovertemplate="%{hovertext}<extra></extra>",
            name="stations",
            showlegend=False,
        )
    )
    if show_labels:
        _add_station_labels(
            fig, xs, ys, zs, labels, options,
            clear_px=int(options.station_size) + 9,
        )


def _add_station_labels(
    fig, xs, ys, zs, labels, options, *, z_lift=0.0, clear_px=10
) -> None:
    """Per-station labels as rotatable, screen-facing scene annotations.

    ``go.Scatter3d`` text has no rotation, so the labels ride on
    ``layout.scene.annotations`` instead -- which do support
    ``textangle`` (``options.station_label_angle``), letting a crowded
    line tilt its labels the way a 2-D pseudosection does.

    ``None`` entries in ``labels`` are markers the label-density /
    named-only filter dropped -- skipped here. Each label is anchored by
    its *bottom* edge and lifted ``clear_px`` screen pixels above the
    marker so a rotated label sits on top of the glyph instead of
    running through it.
    """
    angle = float(getattr(options, "station_label_angle", 0.0) or 0.0)
    color = str(options.station_color)
    notes = [
        dict(
            x=float(x),
            y=float(y),
            z=float(z) + z_lift,
            text=str(t),
            showarrow=False,
            textangle=angle,
            yanchor="bottom",
            yshift=clear_px,
            font=dict(size=9, color=color),
        )
        for x, y, z, t in zip(xs, ys, zs, labels)
        if t is not None
    ]
    if not notes:
        return
    existing = list(getattr(fig.layout.scene, "annotations", None) or [])
    fig.update_scenes(annotations=tuple(existing) + tuple(notes))


def _terrain_trace(x, y, elev, name):
    """Return a Scatter3d line tracing the terrain top for one line."""
    go = require_plotly()

    return go.Scatter3d(
        x=np.asarray(x, dtype=float),
        y=np.asarray(y, dtype=float),
        z=np.asarray(elev, dtype=float),
        mode="lines",
        line=dict(color="#8d6e63", width=5),
        name=f"{name} terrain",
        showlegend=False,
        hoverinfo="skip",
    )


def _pad_row(arr, width):
    """Pad a 1-D array with NaN for ragged stacking."""
    arr = np.asarray(arr, dtype=float)
    if arr.size >= width:
        return arr[:width]
    out = np.full(width, np.nan, dtype=float)
    out[: arr.size] = arr
    return out


# Structure smoothing ("Light / Medium / Strong / Very strong" in the
# mapview). Each level: (a) rebuilds the block/iso volume on a lattice
# ``factor`` times denser per axis (capped per axis, and by a global
# cell budget) so Plotly's marching-cubes has the resolution to draw a
# rounded surface -- this is what actually removes the facets; (b) runs
# a NaN-aware Gaussian of ``sigma_cells`` *refined cells* (small on
# purpose -- a wide blur on the already-dense lattice erodes thin bodies
# away entirely); and (c) -- see :func:`_feather_volume` -- ramps the
# in/out-of-band cliff so a resistivity-band body reads as one smooth
# shape. Strength 0 (the default) skips all of it. Earlier revisions
# blurred the raw *sparse* grid with the preset value straight as a
# sigma, which only blobbed bodies without ever rounding the surfaces.
_VOLUME_SMOOTHING_PRESETS = (
    # (strength_threshold, refine_factor, axis_cap, y_cap, sigma_cells)
    (3.5, 3.4, 170, 105, 3.0),  # very strong
    (2.5, 3.0, 140, 92, 2.2),   # strong
    (1.5, 2.5, 118, 82, 1.5),   # medium
    (0.0, 2.0, 92, 66, 0.9),    # light (any strength > 0)
)
# go.Volume / go.Isosurface run marching-cubes in the browser; past
# ~0.5M lattice cells the render janks badly or drops the trace
# entirely (which is what "Very strong makes the block vanish" was).
# Keep the refined lattice under this budget -- it is still ~10x the
# native station/depth grid, plenty for a facet-free surface once the
# Gaussian + envelope feather are applied. See the cube-root downscale.
_VOLUME_CELL_BUDGET = 500_000


def _volume_smoothing_params(strength):
    """Resolve the UI structure-smoothing strength (0 / 0.8 / 1.5 / 2.5
    / 4.0) to ``(refine_factor, axis_cap, y_cap, sigma_cells)``.

    ``strength <= 0`` -> ``(1.0, 0, 0, 0.0)``: no refinement, no blur --
    the default, untouched volume. Higher strength = denser lattice
    (the main lever) plus a slightly wider Gaussian.
    """
    if not strength or strength <= 0:
        return 1.0, 0, 0, 0.0
    s = float(strength)
    for thresh, factor, cap, y_cap, sigma_cells in _VOLUME_SMOOTHING_PRESETS:
        if s >= thresh:
            return factor, cap, y_cap, sigma_cells
    return _VOLUME_SMOOTHING_PRESETS[-1][1:]


def _axis_sigma(shape, sigma_cells):
    """Per-axis Gaussian sigma for a volume of ``shape``: ``sigma_cells``
    on every axis that has real extent (>= 5 samples), 0 on a thin axis
    (e.g. a 1-line survey's cross-line axis) so the blur never collapses
    it. Returns ``None`` when there is nothing to blur."""
    if not sigma_cells or sigma_cells <= 0:
        return None
    sig = tuple(float(sigma_cells) if n >= 5 else 0.0 for n in shape)
    return sig if any(s > 0 for s in sig) else None


def _smooth_volume(rho_vol, sigma):
    """NaN-aware Gaussian blur of the dense resistivity volume.

    ``sigma`` may be a scalar or a per-axis sequence; ``None`` / all
    non-positive returns the volume untouched (the default). The blur
    is normalised by a matching blur of the finite-mask so no-data
    cells don't drag values toward zero. Cells where the blurred mask
    is below ``0.4`` (i.e. mostly reaching into no-data) are dropped, so
    the smoothed body keeps essentially the original footprint but its
    boundary rounds off by ~one sigma instead of stair-stepping.
    """
    sig = np.atleast_1d(np.asarray(sigma, dtype=float)) if sigma is not None \
        else np.array([0.0])
    if sig.size == 0 or not np.any(sig > 0):
        return rho_vol
    try:
        from scipy.ndimage import gaussian_filter
    except Exception:  # noqa: BLE001 - scipy.ndimage is part of scipy
        return rho_vol
    finite = np.isfinite(rho_vol)
    if not finite.any():
        return rho_vol
    sig_use = float(sig[0]) if sig.size == 1 else tuple(sig.tolist())
    filled = np.where(finite, rho_vol, 0.0)
    num = gaussian_filter(filled, sigma=sig_use, mode="nearest")
    den = gaussian_filter(finite.astype(float), sigma=sig_use, mode="nearest")
    out = np.divide(num, den, out=np.full_like(num, np.nan), where=den > 1e-6)
    out[den < 0.3] = np.nan  # keep ~the original footprint, edge rounded
    if not np.isfinite(out).any():
        # A wide blur on a thin footprint can erase everything -- never
        # let smoothing turn the whole block into no-data.
        return rho_vol
    return out


def _feather_volume(value, sentinel, sigma):
    """Ramp the cliff between the ``sentinel`` fill and the live cells so
    ``go.Isosurface`` cuts a rounded envelope instead of a voxel
    staircase.

    Only the sentinel cells are moved: they are set to a Gaussian-
    blurred field that rises from a moderate floor up toward the band
    near the boundary. **Live cells keep their exact value**, so the
    body can never be blurred out of existence and hues are unchanged.

    ``sigma`` is a per-axis sequence (or scalar). Returns ``value``
    untouched when SciPy is missing or ``sigma`` is empty.
    """
    sig = np.atleast_1d(np.asarray(sigma, dtype=float)) if sigma is not None \
        else np.array([0.0])
    if sig.size == 0 or not np.any(sig > 0):
        return value
    try:
        from scipy.ndimage import gaussian_filter
    except Exception:  # noqa: BLE001
        return value
    sig_use = float(sig[0]) if sig.size == 1 else tuple(sig.tolist())
    live_mask = value > sentinel
    live = value[live_mask]
    if live.size == 0:
        return value
    # Clamp the sentinel to a moderate distance below the live range
    # first: blurring against an extreme outlier makes a needle-thin
    # ramp (still a hard edge). A gentle offset gives a real ramp.
    floor = float(live.min()) - max(float(np.ptp(live)), 1.0) * 0.5
    clamped = np.where(live_mask, value, floor)
    ramp = gaussian_filter(clamped, sigma=sig_use, mode="nearest")
    return np.where(live_mask, value, ramp)


def _dense_volume_grid(profiles, options):
    """Build a regular ``(x, y, z)`` grid + resistivity volume for
    ``go.Volume``/``go.Isosurface``.

    Those trace types need genuine 3-D grid density to reconstruct
    smooth iso-surfaces. Each line's own ``(x, z)`` axes are too
    sparse and mutually inconsistent for that directly — this
    regrids every line onto one shared ``(x, z)`` reference, then
    densifies the cross-line (``y``) axis by interpolation, the same
    technique used by the ``pycsamt.app.web`` 3-D map's block mode.

    Returns ``None`` when there isn't enough data (or SciPy isn't
    available) to build a grid.
    """
    if not profiles:
        return None
    try:
        from scipy.interpolate import RegularGridInterpolator
    except Exception:  # noqa: BLE001 - scipy is an optional extra
        return None

    names = list(profiles.keys())
    n_lines = len(names)
    unit = _line_offset_unit(profiles)
    real_offsets = _line_real_offsets(profiles)
    y_vals = np.array(
        [
            _line_offset(name, idx, real_offsets, unit, options)
            for idx, name in enumerate(names)
        ],
        dtype=float,
    )

    ref_name = max(names, key=lambda n: np.asarray(profiles[n]["x"]).size)
    x_ref = np.unique(np.asarray(profiles[ref_name]["x"], dtype=float))
    z_ref = np.unique(np.asarray(profiles[ref_name]["z"], dtype=float))
    if x_ref.size < 2 or z_ref.size < 2:
        return None

    # Structure smoothing: when it's on, reconstruct every line onto a
    # denser shared (x, z) lattice first so the block/iso surfaces round
    # off instead of faceting on the raw station/depth spacing.
    refine, axis_cap, y_cap, sigma_cells = _volume_smoothing_params(
        float(getattr(options, "volume_smoothing", 0.0) or 0.0)
    )
    n_y_target = None
    if refine > 1:
        n_x_t = min(axis_cap, max(x_ref.size, round(x_ref.size * refine)))
        n_z_t = min(axis_cap, max(z_ref.size, round(z_ref.size * refine)))
        n_y_t = (
            min(y_cap, max(20, n_lines * 6) * 2) if n_lines >= 2 else 1
        )
        # Keep the refined lattice inside a global cell budget so a big
        # survey at "Very strong" stays renderable -- shrink all three
        # target counts by the cube-root of the overshoot, but never
        # below the native resolution (coarsening would lose structure).
        _overshoot = max(n_x_t * n_y_t * n_z_t, 1)
        scale = (_VOLUME_CELL_BUDGET / _overshoot) ** (1 / 3)
        if scale < 1.0:
            n_x_t = max(int(n_x_t * scale), x_ref.size)
            n_z_t = max(int(n_z_t * scale), z_ref.size)
            n_y_t = max(int(n_y_t * scale), min(n_y_t, 20))
        if n_x_t > x_ref.size:
            x_ref = np.linspace(float(x_ref.min()), float(x_ref.max()), n_x_t)
        if n_z_t > z_ref.size:
            z_ref = np.linspace(float(z_ref.min()), float(z_ref.max()), n_z_t)
        n_y_target = n_y_t

    rho_vol = np.full((x_ref.size, n_lines, z_ref.size), np.nan, dtype=float)
    for idx, name in enumerate(names):
        grid = profiles[name]
        x = np.asarray(grid["x"], dtype=float)
        z = np.asarray(grid["z"], dtype=float)
        if x.size < 2 or z.size < 2:
            continue
        values = _color_values(np.asarray(grid["value"], dtype=float), options)
        xo, zo = np.argsort(x), np.argsort(z)
        x_sorted, z_sorted = x[xo], z[zo]
        if np.unique(x_sorted).size < 2 or np.unique(z_sorted).size < 2:
            continue
        vals_sorted = values[np.ix_(zo, xo)]
        try:
            interp = RegularGridInterpolator(
                (z_sorted, x_sorted),
                vals_sorted,
                bounds_error=False,
                fill_value=np.nan,
            )
        except ValueError:
            continue
        zz, xx = np.meshgrid(z_ref, x_ref, indexing="ij")
        sampled = interp(np.column_stack([zz.ravel(), xx.ravel()]))
        rho_vol[:, idx, :] = sampled.reshape(z_ref.size, x_ref.size).T

    elev_lines = None
    if options.topography:
        elev_lines = np.stack(
            [
                _line_elevation_on_x_ref(profiles[name], options, x_ref)
                for name in names
            ],
            axis=1,
        )  # (x_ref.size, n_lines)

    if n_lines >= 2:
        order = np.argsort(y_vals)
        y_sorted = y_vals[order]
        n_y_dense = max(20, n_lines * 6)
        if n_y_target is not None:
            # Cross-line data is genuinely sparse (a handful of lines);
            # densify it too, but modestly -- more than ~2x just invents
            # structure between lines that was never measured. Use the
            # budget-adjusted target from above (floored at one per line).
            n_y_dense = max(n_lines, n_y_target)
        y_dense = np.linspace(y_sorted[0], y_sorted[-1], n_y_dense)
        dense = np.empty((x_ref.size, n_y_dense, z_ref.size), dtype=float)
        for ix in range(x_ref.size):
            for iz in range(z_ref.size):
                dense[ix, :, iz] = _interp_finite_1d(
                    y_dense, y_sorted, rho_vol[ix, order, iz]
                )
        rho_vol, y_vals = dense, y_dense

        if elev_lines is not None:
            elev_dense = np.empty((x_ref.size, n_y_dense), dtype=float)
            for ix in range(x_ref.size):
                elev_dense[ix, :] = _interp_finite_1d(
                    y_dense, y_sorted, elev_lines[ix, order]
                )
            elev_lines = elev_dense

    rho_vol = _smooth_volume(rho_vol, _axis_sigma(rho_vol.shape, sigma_cells))

    z_depth = -np.abs(z_ref)
    if elev_lines is not None:
        rho_vol, z_out = _drape_dense_volume(rho_vol, z_depth, elev_lines)
        return x_ref, y_vals, z_out, rho_vol

    return x_ref, y_vals, z_depth, rho_vol


def _line_elevation_on_x_ref(grid, options, x_ref):
    """Real station elevation for one line, interpolated onto the
    shared ``x_ref`` axis -- same source (:func:`_elev_for`) the fence
    and depth-slice panels already drape by, densified across ``x``
    the same way :func:`_dense_volume_grid` already densifies
    resistivity."""
    x = np.asarray(grid["x"], dtype=float)
    elev = _elev_for(grid, options)
    if elev.size != x.size or x.size < 2:
        return np.zeros_like(x_ref)
    order = np.argsort(x)
    x_sorted, elev_sorted = x[order], elev[order]
    good = np.isfinite(x_sorted) & np.isfinite(elev_sorted)
    if good.sum() < 2:
        return np.zeros_like(x_ref)
    return np.interp(x_ref, x_sorted[good], elev_sorted[good])


def _drape_dense_volume(rho_vol, z_depth, elev_grid):
    """Resample a flat (depth-below-surface) volume onto one shared
    absolute-elevation lattice, terrain-following per ``(x, y)``
    column.

    ``rho_vol`` : (nx, ny, nz) array of values as a function of
        ``z_depth`` (depth below each column's own local surface,
        ``<= 0``, ``0`` at the surface) -- the flat datum
        :func:`_dense_volume_grid` builds before this call.
    ``elev_grid`` : (nx, ny) real elevation (m) at each ``(x, y)``
        column, from the same station topography the fence and
        depth-slice panels already drape by.

    Plotly's ``go.Volume``/``go.Isosurface`` need one shared, literal
    rectangular ``(x, y, z)`` lattice to reconstruct a shape -- unlike
    ``go.Surface``, they cannot be warped per index the way a fence
    panel's terrain-following ``z`` mesh can. So instead of warping
    the *coordinate* grid, ``z`` stays one shared absolute-elevation
    axis and the *values* are resampled into it: for each ``(x, y)``
    column, a query elevation is converted back to that column's own
    local depth below its real surface and sampled from the flat
    volume at that depth. A query above the local terrain (or deeper
    than any profile actually sampled) is left ``NaN`` -- real air or
    genuine no-data, never a fabricated value -- which is exactly what
    makes the rendered block's top surface follow real terrain instead
    of a flat datum.
    """
    nx, ny, nz = rho_vol.shape
    finite_elev = elev_grid[np.isfinite(elev_grid)]
    if finite_elev.size == 0:
        return rho_vol, z_depth
    top = float(finite_elev.max())
    bottom = float(finite_elev.min() + z_depth.min())
    z_abs = np.linspace(bottom, top, nz)

    order = np.argsort(z_depth)
    z_sorted = z_depth[order]

    draped = np.full((nx, ny, nz), np.nan, dtype=float)
    for ix in range(nx):
        for iy in range(ny):
            elev = elev_grid[ix, iy]
            if not np.isfinite(elev):
                continue
            col = rho_vol[ix, iy, order]
            good = np.isfinite(col)
            if good.sum() < 2:
                continue
            depth_query = z_abs - elev
            sampled = np.interp(
                depth_query,
                z_sorted[good],
                col[good],
                left=np.nan,
                right=np.nan,
            )
            # Above the local terrain is real air, not no-data below
            # the profile's own range -- mask it explicitly rather
            # than let np.interp's edge clamp show a fabricated value.
            sampled = np.where(depth_query > 0, np.nan, sampled)
            draped[ix, iy, :] = sampled
    return draped, z_abs


def _interp_finite_1d(x_new, x, values):
    """1-D interpolation that skips NaNs instead of propagating them."""
    x = np.asarray(x, dtype=float)
    values = np.asarray(values, dtype=float)
    good = np.isfinite(x) & np.isfinite(values)
    if good.sum() >= 2:
        return np.interp(x_new, x[good], values[good])
    if good.sum() == 1:
        return np.full_like(x_new, values[good][0], dtype=float)
    return np.full_like(x_new, np.nan, dtype=float)


def _iso_range(vals, options):
    lo, hi = _iso_color_range(options)
    if lo is not None and hi is not None:
        return lo, hi
    arr = np.asarray(vals, dtype=float)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return 0.0, 1.0
    return (
        float(np.nanpercentile(finite, 20)),
        float(np.nanpercentile(finite, 80)),
    )


def _iso_color_range(options) -> tuple[float | None, float | None]:
    """``(isomin, isomax)`` in colour space, from the user-facing
    resistivity-range filter (``options.rho_range``) so it masks the
    block/iso-surface the same way it already masks fence/depth-slice
    cell values — without pre-filtering the dense grid to NaN, which
    would starve ``go.Volume``/``go.Isosurface`` of the density they
    need to interpolate a smooth shape.
    """
    if options.rho_range:
        lo, hi = options.rho_range
        if _volume_quantity(options.quantity) == "rho" and options.log_color:
            lo = float(np.log10(lo)) if lo and lo > 0 else None
            hi = float(np.log10(hi)) if hi and hi > 0 else None
        return lo, hi
    if options.iso_range:
        return options.iso_range
    return None, None


def _volume_quantity(quantity: str) -> str:
    if quantity.lower() in {"phase", "phi"}:
        return "phase"
    return "rho"


def _colorbar_title(options) -> str:
    if _volume_quantity(options.quantity) == "phase":
        return "φ (°)"
    return "log₁₀ ρ (Ω·m)" if options.log_color else "ρ (Ω·m)"


def _slice_depths(profiles, options):
    """Depths (m below surface) to cut a horizontal resistivity slice at.

    Always kept inside the data's real z-extent -- a slice at a depth
    the model never sampled comes back all-NaN and renders as an empty
    ghost surface. ``n_slices == 1`` cuts through the *middle* of the
    requested window (so "0-200 m, 1 slice" gives the 100 m slice, not
    the 0 m one).
    """
    all_z = np.concatenate(
        [np.asarray(grid["z"], dtype=float) for grid in profiles.values()]
    )
    all_z = all_z[np.isfinite(all_z)]
    if all_z.size == 0:
        return np.array([0.0])
    z_lo, z_hi = float(np.nanmin(all_z)), float(np.nanmax(all_z))
    req_lo, req_hi = (
        (float(options.depth_range[0]), float(options.depth_range[1]))
        if options.depth_range
        else (z_lo, z_hi)
    )
    n = max(1, int(options.n_slices))
    if n == 1:
        # the middle of the *requested* window, clamped into real data
        mid = 0.5 * (req_lo + req_hi)
        return np.array([float(np.clip(mid, z_lo, z_hi))])
    lo = max(req_lo, z_lo)
    hi = min(req_hi, z_hi)
    if not (hi > lo):
        return np.array([float(np.clip(0.5 * (z_lo + z_hi), z_lo, z_hi))])
    return np.linspace(lo, hi, n)


def _surface_contours(options):
    go = require_plotly()

    return go.surface.Contours(
        z=go.surface.contours.Z(
            show=bool(options.show_contours),
            usecolormap=True,
            project_z=False,
        )
    )


def _style_3d(fig, options, colors) -> None:
    title = options.title or f"pyCSAMT 3-D {options.mode} map"
    xu = getattr(options, "x_unit", "m")
    zu = getattr(options, "depth_unit", "m")
    z_label = "Elevation − depth" if options.topography else "Depth"
    z_title = f"{z_label} ({zu})"
    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title=f"Profile distance ({xu})",
            yaxis_title=f"Line offset ({xu})",
            zaxis_title=z_title,
            bgcolor=colors["plot"],
            aspectmode=getattr(options, "aspectmode", "data") or "data",
        ),
        paper_bgcolor=colors["paper"],
        font=dict(color=colors["text"]),
        margin=dict(l=0, r=0, t=40, b=0),
    )


def _empty_3d_figure(colors):
    go = require_plotly()

    fig = go.Figure()
    fig.add_annotation(
        text="No 3-D map data available",
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
