# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
MapPanel — multi-type station map rendered in Matplotlib.

Map types
---------
station      Coloured scatter (default)
elevation    Elevation colour map with optional contours
depth        Skin-depth (Niblett-Bostick) map at user-chosen depth
resistivity  Apparent ρ map at user-chosen frequency / period

Basemap overlay
---------------
When a basemap provider is selected station coordinates are automatically
projected to Web Mercator (EPSG:3857) before plotting so that contextily can
add tiles WITHOUT using the `crs=` parameter.  This sidesteps the common
conda/pip PROJ database version conflict that causes a silent CRSError when
`crs="EPSG:4326"` is passed.

The pyproj WGS84 → Mercator transform is purely arithmetic and does NOT
require a valid proj.db, so it works even in mixed conda/pip environments.

Colourbar
---------
Always placed as an inset inside the map axes (ax.inset_axes), so it never
steals layout space and the map cannot shrink on repeated Refresh calls.

Note on fig.clf()
-----------------
Every _draw_map() call does fig.clf() + fig.add_subplot(111), giving a
completely clean figure on every refresh — no accumulated axes, no residual
layout adjustments.
"""

from __future__ import annotations

import importlib
import subprocess
import sys

import numpy as np
import pandas as pd
from PySide6.QtCore import QEvent, Qt, Signal
from PySide6.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas

# ── optional deps ─────────────────────────────────────────────────────────────


def _try_ctx():
    try:
        import contextily as _c

        return _c
    except ImportError:
        return None


def _try_griddata():
    try:
        from scipy.interpolate import griddata

        return griddata
    except ImportError:
        return None


# ── coordinate helpers ────────────────────────────────────────────────────────


def _project_to_merc(
    xs: np.ndarray, ys: np.ndarray, source_crs: str = "EPSG:4326"
):
    """
    Project coordinates from *source_crs* to Web Mercator (EPSG:3857).

    Works for any valid EPSG / proj string: WGS84, UTM zones, national grids
    (CGCS2000, RGF93, GDA2020 …) or any other EPSG the user supplies.

    The WGS84 → Mercator path is purely arithmetic and survives the common
    conda/pip PROJ database version conflict.  Other pairs use pyproj's
    built-in tables (EPSG codes are embedded in the pyproj wheel, so they
    work even when the system proj.db is stale).
    """
    from pyproj import Transformer

    t = Transformer.from_crs(source_crs, "EPSG:3857", always_xy=True)
    return t.transform(xs, ys)


def _is_geographic(crs_str: str) -> bool:
    """Return True when *crs_str* uses angular units (degrees, not metres)."""
    try:
        from pyproj import CRS

        return CRS.from_user_input(crs_str).is_geographic
    except Exception:
        # Safe fallback: treat only well-known geographic codes as geographic
        code = crs_str.upper().replace("EPSG:", "").strip()
        return code in ("4326", "4490", "4610", "4214", "4269", "4258")


def _merc_tick_formatter():
    """Return (lon_fmt, lat_fmt) callables for Mercator-axis tick labels."""
    from pyproj import Transformer

    inv = Transformer.from_crs("EPSG:3857", "EPSG:4326", always_xy=True)

    def fmt_lon(val, _pos):
        lon, _ = inv.transform(val, 0)
        return f"{lon:.4g}°"

    def fmt_lat(val, _pos):
        _, lat = inv.transform(0, val)
        return f"{lat:.4g}°"

    return fmt_lon, fmt_lat


# ── impedance helpers ─────────────────────────────────────────────────────────

_COMP_IDX = {"x": 0, "y": 1}


def _slice_comp(z: np.ndarray, comp: str) -> np.ndarray:
    a, b = _COMP_IDX[comp[0]], _COMP_IDX[comp[1]]
    return z[:, a, b]


def _rho_app(zc: np.ndarray, freqs: np.ndarray) -> np.ndarray:
    return 0.2 * np.abs(zc) ** 2 / (freqs + 1e-24)


def _reshape_z(z) -> np.ndarray | None:
    if z is None:
        return None
    z = np.asarray(z)
    if z.ndim == 3 and z.shape[1] == 2 and z.shape[2] == 2:
        return z
    if z.ndim == 2 and z.shape[1] >= 4:
        return z[:, :4].reshape(-1, 2, 2)
    return None


# ── contextily provider table ─────────────────────────────────────────────────

_BASEMAP_PROVIDERS: dict = {}


def _build_provider_table(ctx) -> dict:
    p = ctx.providers
    tbl = {}
    for name, getter in [
        ("OpenStreetMap", lambda: p.OpenStreetMap.Mapnik),
        ("Stamen Terrain", lambda: p.Stamen.Terrain),
        ("Satellite (ESRI)", lambda: p.Esri.WorldImagery),
        ("CartoDB Positron", lambda: p.CartoDB.Positron),
        ("CartoDB DarkMatter", lambda: p.CartoDB.DarkMatter),
    ]:
        try:
            tbl[name] = getter()
        except AttributeError:
            pass
    return tbl


# ── colourbar inset positions [x0, y0, w, h] in axes-fraction coords ─────────

_CBAR_INSET = {
    "vertical": [0.965, 0.08, 0.025, 0.45],
    "horizontal": [0.05, 0.03, 0.45, 0.030],
}


# ── hover-reveal pop-out button ───────────────────────────────────────────────


class _MapPopOutButton(QPushButton):
    """
    Transparent expand button that floats over the MapPanel canvas.

    Invisible by default; fades in when the mouse enters the MapPanel.
    Click opens a larger MapViewerWindow pre-loaded with the current data and
    all rendering settings so the user can explore the map in full detail.
    """

    _SS_HIDDEN = """QPushButton {
        background: transparent;
        border: none;
        color:       transparent;
    }"""
    _SS_SHOWN = """QPushButton {
        background:   rgba(18, 20, 32, 0.72);
        border:       1px solid rgba(255, 255, 255, 0.22);
        border-radius: 7px;
        color:        rgba(255, 255, 255, 0.90);
        font-size:    14px;
    }
    QPushButton:hover {
        background:   rgba(55, 60, 90, 0.90);
        border-color: rgba(255, 255, 255, 0.45);
        color:        white;
    }
    QPushButton:pressed {
        background:   rgba(80, 86, 130, 0.95);
    }"""

    def __init__(self, container: MapPanel) -> None:
        super().__init__("⤢", container)
        self._container = container
        self._shown = False

        self.setFixedSize(30, 30)
        self.setToolTip("Open in full map viewer with all controls")
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground)
        self.setStyleSheet(self._SS_HIDDEN)

        self.clicked.connect(self._on_click)
        container.installEventFilter(self)
        self._reposition()
        self.raise_()

    # ── event filter on the MapPanel container ────────────────────────────────

    def eventFilter(self, obj, event) -> bool:
        # ``self`` can be re-entered here while its own C++ object is being
        # torn down: destroying this button sends a ChildRemoved event to
        # ``container``, and since this button is registered as one of
        # ``container``'s filters, Qt calls back into this very instance --
        # possibly after CPython has already cleared its ``__dict__`` as
        # part of deallocation. ``getattr`` guards that reentrant case
        # instead of raising out of a virtual-method override.
        if obj is getattr(self, "_container", None):
            t = event.type()
            if t == QEvent.Type.Resize:
                self._reposition()
                self.raise_()
            elif t == QEvent.Type.Enter:
                self._set_shown(True)
            elif t == QEvent.Type.Leave:
                self._set_shown(False)
        return super().eventFilter(obj, event)

    def _set_shown(self, shown: bool) -> None:
        if shown == self._shown:
            return
        self._shown = shown
        self.setStyleSheet(self._SS_SHOWN if shown else self._SS_HIDDEN)
        if shown:
            self.raise_()

    def _reposition(self) -> None:
        m = 8
        self.move(self._container.width() - self.width() - m, m)

    # ── click → open full viewer ──────────────────────────────────────────────

    def _on_click(self) -> None:
        panel = self._container
        # Deferred import avoids circular dependency (map_window imports MapPanel)
        from pycsamt.app.desktop.windows.map_window import (
            MapDetailWindow,
        )

        win = MapDetailWindow(parent=panel.window())
        # Load data into the new window
        if not panel._df.empty:
            win.set_dataframe(panel._df)
        if panel._sites is not None:
            win.set_sites(panel._sites)
        win.set_dark_mode(panel._dark)
        # Sync all controls and render from the current embedded panel state
        win.sync_from_panel(panel)
        win.show()


# ── MapPanel ──────────────────────────────────────────────────────────────────


class MapPanel(QWidget):
    """
    Interactive station map.

    Signals
    -------
    station_selected : str
    freq_list_ready  : list[float]
    """

    station_selected = Signal(str)
    freq_list_ready = Signal(list)

    # defaults
    _map_type: str = "station"
    _color_by: str = "index"
    _cmap_name: str = "plasma"
    _marker_size: int = 8
    _marker_alpha: float = 0.9
    _show_labels: bool = True
    _show_profile: bool = True
    _show_grid: bool = True
    _provider: str = "None"
    _basemap_alpha: float = 0.7
    _contour_mode: str = "none"
    _contour_levels: int = 10
    _component: str = "xy"
    _target_freq_hz: float = 1.0
    _target_depth_m: float = 500.0
    _log_scale: bool = True
    _cbar_orient: str = "vertical"
    _show_cbar: bool = True
    # CRS of the input station coordinates (standard EDI = EPSG:4326)
    _source_crs: str = "EPSG:4326"

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._dark: bool = True
        self._df: pd.DataFrame = pd.DataFrame(
            columns=["ID", "Latitude", "Longitude"]
        )
        self._sites = None
        self._selected_id: str | None = None
        self._scatter = None
        self._annots: dict = {}
        self._build_ui()

    # ── construction ──────────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        bar = QHBoxLayout()
        bar.setContentsMargins(4, 2, 4, 2)
        self._info_label = QLabel("")
        self._info_label.setObjectName("InfoLabel")
        bar.addWidget(self._info_label)
        bar.addStretch()
        layout.addLayout(bar)

        self._canvas = MplCanvas(self, toolbar=True)
        layout.addWidget(self._canvas)

        # Hover-reveal pop-out button — floats over the top-right corner.
        # Keep a reference on self: PySide6 only owns the button on the C++
        # side via its `container` parent, so an unassigned Python wrapper
        # can be garbage-collected while the C++ object survives. Qt then
        # has to re-wrap it from the bare C++ pointer without going through
        # __init__, and the resulting instance is missing `_container`.
        self._pop_out_button = _MapPopOutButton(container=self)

        # Connect once — survive fig.clf() because they're on the canvas object
        self._canvas.figure.canvas.mpl_connect("pick_event", self._on_pick)
        self._canvas.figure.canvas.mpl_connect(
            "button_press_event", self._on_click
        )

    # ── public API ────────────────────────────────────────────────────────────

    def set_dataframe(self, df: pd.DataFrame) -> None:
        self._df = df.copy().dropna(subset=["Latitude", "Longitude"])
        self._selected_id = None
        self._draw_map()

    def set_sites(self, sites) -> None:
        self._sites = sites
        freqs = self._collect_frequencies()
        if freqs:
            self.freq_list_ready.emit(freqs)
        self._draw_map()

    def highlight_station(self, station_id: str) -> None:
        self._selected_id = station_id
        self._update_highlight()

    def set_dark_mode(self, dark: bool) -> None:
        self._dark = dark
        self._draw_map()

    def clear(self) -> None:
        self._df = pd.DataFrame(columns=["ID", "Latitude", "Longitude"])
        self._selected_id = None
        self._scatter = None
        self._annots = {}
        self._canvas.axes.cla()
        self._canvas.draw()

    def redraw(self, **kwargs) -> None:
        for k, v in kwargs.items():
            if hasattr(self, f"_{k}"):
                setattr(self, f"_{k}", v)
        self._draw_map()

    # ── frequency collection ──────────────────────────────────────────────────

    def _collect_frequencies(self) -> list:
        if self._sites is None:
            return []
        freqs_set: set = set()
        for site in self._sites:
            try:
                f = np.asarray(site.freq, float)
                freqs_set.update(f[np.isfinite(f) & (f > 0)].tolist())
            except Exception:
                pass
        return sorted(freqs_set, reverse=True)

    # ── impedance extraction ──────────────────────────────────────────────────

    def _z_arrays(self, site):
        try:
            from pycsamt.site.base import _extract_z_arrays

            arrs = _extract_z_arrays(site.edi)
            freqs = np.asarray(arrs["freq"], float)
            z = _reshape_z(arrs["z"])
            if freqs.size == 0 or z is None:
                return None, None
            mask = np.isfinite(freqs) & (freqs > 0)
            return freqs[mask], z[mask]
        except Exception:
            return None, None

    def _rho_comp(
        self, z3d: np.ndarray, freqs: np.ndarray, comp: str
    ) -> np.ndarray:
        if comp == "det":
            zc = np.sqrt(np.abs(z3d[:, 0, 1] * z3d[:, 1, 0]))
        else:
            zc = _slice_comp(z3d, comp)
        return _rho_app(zc, freqs)

    def _rho_at_freq(self, target_hz: float, comp: str = "xy") -> dict:
        out: dict = {}
        for site in self._sites or []:
            freqs, z3d = self._z_arrays(site)
            if freqs is None:
                continue
            try:
                rho = self._rho_comp(z3d, freqs, comp)
                idx = int(np.argmin(np.abs(freqs - target_hz)))
                val = float(rho[idx])
                if np.isfinite(val) and val > 0:
                    out[str(site.name)] = val
            except Exception:
                pass
        return out

    def _rho_at_depth(self, target_m: float, comp: str = "xy") -> dict:
        """Niblett-Bostick: δ (m) = 503 × √(ρ_app / f)."""
        out: dict = {}
        for site in self._sites or []:
            freqs, z3d = self._z_arrays(site)
            if freqs is None:
                continue
            try:
                rho = self._rho_comp(z3d, freqs, comp)
                skin = 503.0 * np.sqrt(
                    np.maximum(rho, 1e-10) / (freqs + 1e-24)
                )
                idx = int(np.argmin(np.abs(skin - target_m)))
                val = float(rho[idx])
                if np.isfinite(val) and val > 0:
                    out[str(site.name)] = val
            except Exception:
                pass
        return out

    # ── main drawing entry ────────────────────────────────────────────────────

    def _draw_map(self) -> None:
        # Completely reset figure → no accumulated axes, no layout drift
        fig = self._canvas.figure
        fig.clf()
        ax = fig.add_subplot(111)
        self._canvas.axes = ax
        self._scatter = None
        self._annots = {}

        if self._df.empty:
            ax.set_title("No stations loaded", fontsize=10, color="gray")
            self._canvas.draw()
            return

        # Decide coordinate system up front
        use_merc = self._provider not in (None, "None", "")
        raw_xs = self._df["Longitude"].values  # column name kept for compat
        raw_ys = self._df["Latitude"].values
        is_geo = _is_geographic(self._source_crs)

        if use_merc:
            try:
                xs, ys = _project_to_merc(raw_xs, raw_ys, self._source_crs)
            except Exception as exc:
                use_merc = False
                self._info_label.setText(f"CRS projection error: {exc}")
                xs, ys = raw_xs, raw_ys
        else:
            xs, ys = raw_xs, raw_ys

        ids = self._df["ID"].values
        mt = self._map_type

        if mt == "elevation":
            self._draw_typed_map(
                ax,
                xs,
                ys,
                ids,
                self._elevation_values(),
                cbar_label="Elevation (m)",
                title="Elevation Map",
                log_scale=False,
            )
        elif mt == "depth":
            vm = self._rho_at_depth(self._target_depth_m, self._component)
            self._draw_typed_map(
                ax,
                xs,
                ys,
                ids,
                vm,
                cbar_label=(
                    r"$\rho_a\ (\Omega\cdot m)$  depth ≈ "
                    f"{self._target_depth_m:.0f} m"
                ),
                title=(
                    f"Depth Map  ({self._component.upper()}, "
                    f"depth ≈ {self._target_depth_m:.0f} m)"
                ),
                log_scale=self._log_scale,
            )
        elif mt == "resistivity":
            T = 1.0 / self._target_freq_hz if self._target_freq_hz > 0 else 0
            vm = self._rho_at_freq(self._target_freq_hz, self._component)
            self._draw_typed_map(
                ax,
                xs,
                ys,
                ids,
                vm,
                cbar_label=(
                    r"$\rho_a\ (\Omega\cdot m)$  "
                    f"f={self._target_freq_hz:.4g} Hz  T={T:.4g} s"
                ),
                title=(
                    f"Resistivity Map  ({self._component.upper()}, "
                    f"{self._target_freq_hz:.4g} Hz)"
                ),
                log_scale=self._log_scale,
            )
        else:
            self._draw_station_map(ax, xs, ys, ids)

        # Basemap — always added last so it appears behind scatter
        if use_merc:
            self._add_basemap_mercator(ax)
            self._format_merc_ticks(ax)

        self._style_axes(ax, is_geo=is_geo, use_merc=use_merc)
        self._canvas.draw()

    # ── station map ───────────────────────────────────────────────────────────

    def _draw_station_map(
        self, ax, xs: np.ndarray, ys: np.ndarray, ids: np.ndarray
    ) -> None:
        cby = self._color_by.lower()
        if cby == "elevation" and "Elevation" in self._df.columns:
            colours = self._df["Elevation"].values.astype(float)
            cbar_label = "Elevation (m)"
        elif cby == "n_freq" and "N_freq" in self._df.columns:
            colours = self._df["N_freq"].values.astype(float)
            cbar_label = "N frequencies"
        elif cby == "tipper" and "Tipper" in self._df.columns:
            colours = self._df["Tipper"].astype(float).values
            cbar_label = "Has Tipper"
        else:
            colours = np.arange(len(xs), dtype=float)
            cbar_label = "Station index"

        self._scatter = ax.scatter(
            xs,
            ys,
            c=colours,
            cmap=self._cmap_name,
            s=self._marker_size**2,
            linewidths=0.6,
            edgecolors="white",
            picker=True,
            pickradius=6,
            alpha=self._marker_alpha,
            zorder=3,
        )

        self._add_colorbar(ax, self._scatter, cbar_label, log_scale=False)

        if self._show_profile:
            self._add_profile_lines(ax, xs, ys)
        if self._show_labels:
            self._add_labels(ax, xs, ys, ids)

        self._update_highlight()

    # ── typed map (elevation / depth / resistivity) ───────────────────────────

    def _draw_typed_map(
        self,
        ax,
        xs: np.ndarray,
        ys: np.ndarray,
        ids: np.ndarray,
        value_map: dict,
        cbar_label: str,
        title: str,
        log_scale: bool = False,
    ) -> None:

        raw = np.array([value_map.get(str(sid), np.nan) for sid in ids], float)
        finite = np.isfinite(raw)

        if not finite.any():
            ax.set_title(
                f"{title}\n(no valid data — check EDI elevation/impedance)",
                fontsize=9,
                color="gray",
            )
            ax.scatter(xs, ys, c="gray", s=self._marker_size**2, zorder=3)
            return

        plot_vals = raw.copy()
        if log_scale:
            valid = finite & (raw > 0)
            plot_vals = np.where(valid, np.log10(raw), np.nan)
            finite = np.isfinite(plot_vals)

        # Contour under scatter
        self._add_contour(
            ax, xs[finite], ys[finite], plot_vals[finite], log_scale
        )

        vmin = (
            float(np.nanpercentile(plot_vals, 5)) if finite.sum() > 2 else None
        )
        vmax = (
            float(np.nanpercentile(plot_vals, 95))
            if finite.sum() > 2
            else None
        )

        self._scatter = ax.scatter(
            xs,
            ys,
            c=plot_vals,
            cmap=self._cmap_name,
            s=self._marker_size**2,
            linewidths=0.6,
            edgecolors="white",
            picker=True,
            pickradius=6,
            alpha=self._marker_alpha,
            zorder=5,
            vmin=vmin,
            vmax=vmax,
        )

        self._add_colorbar(ax, self._scatter, cbar_label, log_scale=log_scale)

        if self._show_profile:
            self._add_profile_lines(ax, xs, ys)
        if self._show_labels:
            self._add_labels(ax, xs, ys, ids)

        ax.set_title(title, fontsize=9)
        self._update_highlight()

    # ── colourbar — always an inset (never resizes the map) ──────────────────

    def _add_colorbar(self, ax, mappable, label: str, log_scale: bool) -> None:
        if not self._show_cbar:
            return
        fig = self._canvas.figure
        tc = "white" if self._dark else "#222222"

        pos = _CBAR_INSET.get(self._cbar_orient, _CBAR_INSET["vertical"])
        cax = ax.inset_axes(pos)
        orient = self._cbar_orient

        cbar = fig.colorbar(mappable, cax=cax, orientation=orient)
        cbar.set_label(label, fontsize=7, color=tc)
        cbar.ax.tick_params(labelsize=6, colors=tc)
        for spine in cbar.ax.spines.values():
            spine.set_edgecolor(tc)
            spine.set_linewidth(0.5)

        if log_scale:
            try:
                import matplotlib.ticker as mticker

                fmt = mticker.FuncFormatter(lambda v, _: f"{10**v:.3g}")
                if orient == "vertical":
                    cbar.ax.yaxis.set_major_formatter(fmt)
                else:
                    cbar.ax.xaxis.set_major_formatter(fmt)
            except Exception:
                pass

    # ── contour ───────────────────────────────────────────────────────────────

    def _add_contour(
        self,
        ax,
        xs: np.ndarray,
        ys: np.ndarray,
        values: np.ndarray,
        log_scale: bool,
    ) -> None:
        mode = self._contour_mode
        if mode == "none" or len(xs) < 4:
            return
        griddata = _try_griddata()
        if griddata is None:
            self._info_label.setText("Install scipy for contour interpolation")
            return
        self._info_label.setText("")

        xi = np.linspace(xs.min(), xs.max(), 200)
        yi = np.linspace(ys.min(), ys.max(), 200)
        Xi, Yi = np.meshgrid(xi, yi)
        try:
            Zi = griddata((xs, ys), values, (Xi, Yi), method="cubic")
            Zi_lin = griddata((xs, ys), values, (Xi, Yi), method="linear")
            bad = ~np.isfinite(Zi)
            Zi[bad] = Zi_lin[bad]
        except Exception:
            return

        if not np.isfinite(Zi).any():
            return

        lvl = self._contour_levels
        cmap = self._cmap_name
        fmt = (lambda v, _: f"{10**v:.3g}") if log_scale else "%.2g"

        try:
            if mode == "lines":
                ax.contour(
                    Xi, Yi, Zi, levels=lvl, cmap=cmap, alpha=0.85, zorder=2
                )
            elif mode == "filled":
                ax.contourf(
                    Xi, Yi, Zi, levels=lvl, cmap=cmap, alpha=0.55, zorder=1
                )
                ax.contour(
                    Xi,
                    Yi,
                    Zi,
                    levels=lvl,
                    colors="k",
                    linewidths=0.3,
                    alpha=0.4,
                    zorder=2,
                )
            elif mode == "filled_labels":
                ax.contourf(
                    Xi, Yi, Zi, levels=lvl, cmap=cmap, alpha=0.55, zorder=1
                )
                cs = ax.contour(
                    Xi,
                    Yi,
                    Zi,
                    levels=lvl,
                    colors="k",
                    linewidths=0.4,
                    alpha=0.6,
                    zorder=2,
                )
                ax.clabel(cs, inline=True, fontsize=6, fmt=fmt)
        except Exception:
            pass

    # ── elevation ─────────────────────────────────────────────────────────────

    def _elevation_values(self) -> dict:
        out: dict = {}
        if "Elevation" not in self._df.columns:
            return out
        for _, row in self._df.iterrows():
            v = row["Elevation"]
            if pd.notna(v):
                out[str(row["ID"])] = float(v)
        return out

    # ── overlays ──────────────────────────────────────────────────────────────

    def _add_profile_lines(self, ax, xs: np.ndarray, ys: np.ndarray) -> None:
        if len(xs) < 2:
            return
        lc = "#bbbbbb" if self._dark else "#555555"
        ax.plot(xs, ys, "-", color=lc, linewidth=0.9, alpha=0.5, zorder=2)

    def _add_labels(
        self, ax, xs: np.ndarray, ys: np.ndarray, ids: np.ndarray
    ) -> None:
        tc = "white" if self._dark else "#222222"
        for x, y, sid in zip(xs, ys, ids):
            self._annots[str(sid)] = ax.annotate(
                sid,
                xy=(x, y),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=7,
                color=tc,
                clip_on=True,
                zorder=6,
            )

    # ── basemap (Web Mercator, no CRS conversion) ─────────────────────────────

    def _add_basemap_mercator(self, ax) -> None:
        """
        Add basemap tiles to axes that are already in EPSG:3857 (meters).

        contextily.add_basemap() is called WITHOUT the crs= parameter, which
        avoids the PROJ database lookup that fails in mixed conda/pip setups.
        """
        ctx = _try_ctx()
        if ctx is None:
            self._ask_install_contextily()
            return

        global _BASEMAP_PROVIDERS
        if not _BASEMAP_PROVIDERS:
            _BASEMAP_PROVIDERS = _build_provider_table(ctx)

        src = _BASEMAP_PROVIDERS.get(self._provider)
        if src is None:
            # Provider not available in this contextily version
            self._info_label.setText(
                f"'{self._provider}' not available in contextily "
                f"{getattr(ctx, '__version__', '?')} — try OpenStreetMap or CartoDB."
            )
            return

        try:
            ctx.add_basemap(
                ax,
                source=src,
                zoom="auto",
                alpha=self._basemap_alpha,
            )
            self._info_label.setText("")
        except Exception as exc:
            msg = str(exc)
            self._info_label.setText(f"Basemap error: {msg[:100]}")

    def _format_merc_ticks(self, ax) -> None:
        """Reformat Mercator-axes ticks to display lat/lon degree labels."""
        try:
            import matplotlib.ticker as mticker

            fmt_lon, fmt_lat = _merc_tick_formatter()
            ax.xaxis.set_major_formatter(mticker.FuncFormatter(fmt_lon))
            ax.yaxis.set_major_formatter(mticker.FuncFormatter(fmt_lat))
        except Exception:
            pass
        ax.set_xlabel("Longitude", fontsize=9)
        ax.set_ylabel("Latitude", fontsize=9)

    # ── install contextily ────────────────────────────────────────────────────

    def _ask_install_contextily(self) -> None:
        reply = QMessageBox.question(
            self,
            "Install basemap library?",
            "Basemap tiles require 'contextily', which is not installed.\n\n"
            "Install it now?  (pip install contextily — requires internet)",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
        )
        if reply == QMessageBox.StandardButton.Yes:
            self._install_contextily()

    def _install_contextily(self) -> None:
        try:
            res = subprocess.run(
                [sys.executable, "-m", "pip", "install", "contextily"],
                capture_output=True,
                text=True,
                timeout=180,
            )
            if res.returncode == 0:
                importlib.invalidate_caches()
                QMessageBox.information(
                    self,
                    "Installed",
                    "contextily installed.\nClick Refresh to load the basemap.",
                )
            else:
                QMessageBox.warning(
                    self,
                    "Install failed",
                    f"pip install contextily failed:\n{res.stderr[-600:]}",
                )
        except Exception as exc:
            QMessageBox.warning(
                self, "Install error", f"Could not run pip:\n{exc}"
            )

    # ── axes styling ──────────────────────────────────────────────────────────

    def _style_axes(
        self, ax, is_geo: bool = True, use_merc: bool = False
    ) -> None:
        bg = "#181825" if self._dark else "#eff1f5"
        fbg = "#1e1e2e" if self._dark else "#e6e9ef"
        tc = "white" if self._dark else "#222222"

        ax.set_facecolor(bg)
        self._canvas.figure.patch.set_facecolor(fbg)

        # Axis labels depend on coordinate system in use
        if use_merc:
            # _format_merc_ticks() already set "Longitude"/"Latitude"
            pass
        elif is_geo:
            ax.set_xlabel("Longitude", fontsize=9, color=tc)
            ax.set_ylabel("Latitude", fontsize=9, color=tc)
        else:
            ax.set_xlabel(
                f"Easting (m)  [{self._source_crs}]", fontsize=8, color=tc
            )
            ax.set_ylabel("Northing (m)", fontsize=9, color=tc)

        ax.tick_params(labelsize=8, colors=tc)
        for spine in ax.spines.values():
            spine.set_edgecolor(tc)
        ax.grid(self._show_grid, linewidth=0.4, alpha=0.25, color=tc)

    # ── highlight ─────────────────────────────────────────────────────────────

    def _update_highlight(self) -> None:
        if self._scatter is None:
            return
        ids = self._df["ID"].values
        sizes = np.full(len(ids), self._marker_size**2, dtype=float)
        widths = np.full(len(ids), 0.6)
        ec = ["white"] * len(ids)

        if self._selected_id is not None:
            for i, sid in enumerate(ids):
                if str(sid) == self._selected_id:
                    sizes[i] = (self._marker_size * 2.2) ** 2
                    widths[i] = 2.0
                    ec[i] = "#f5c542"

        self._scatter.set_sizes(sizes)
        self._scatter.set_linewidths(widths)
        self._scatter.set_edgecolors(ec)
        self._canvas.draw()

    # ── interaction ───────────────────────────────────────────────────────────

    def _on_pick(self, event) -> None:
        if event.artist is not self._scatter:
            return
        ind = event.ind
        if not len(ind):
            return
        clicked_id = str(self._df.iloc[ind[0]]["ID"])
        self._selected_id = clicked_id
        self._update_highlight()
        self.station_selected.emit(clicked_id)

    def _on_click(self, event) -> None:
        if event.inaxes is None:
            return
