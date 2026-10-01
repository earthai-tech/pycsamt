# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Topography and Map View hand-off for the PCSF 3-D window (Qt-free).

Topography
----------
A ``.pcsf``/``.pcsm`` model can carry elevations three ways: the station
table's ``z``, a per-station topography block, or a topography *raster*
(``TopographyRaster``: a grid of elevations in the model frame).
:meth:`pycsamt.map.MapView.from_pcsf` only reads the first two, so stations
the table leaves blank stay without ground.  :func:`file_elevations`
collects all three (the raster is sampled at each station's model x/y), and
:func:`fill_line_gaps` interpolates any station still missing along its
line, so the draped surface has no holes.

Map View hand-off
-----------------
:func:`mapview_state` turns the window's scene settings, elevations and
overlays into the JSON the Map View app applies on first load
(``pycsamt-mapview --state``): it opens straight on the 3-D view with the
same mode, colours, depth window, topography, overlays and spin.
"""

from __future__ import annotations

import math

import numpy as np

from pycsamt.api.colormaps import colormap_choices

__all__ = [
    "CMAPS",
    "DEPTH_PRESETS",
    "MODES",
    "file_elevations",
    "fill_line_gaps",
    "mapview_controls",
    "mapview_state",
]

# the shared catalogue, as Map View offers it (a hand-off always finds its
# colour map); RdYlBu_r stays this window's default
CMAPS = tuple(colormap_choices(first=("RdYlBu_r",)))
MODES = (("fence", "Fence"), ("block", "Block"), ("depth", "Depth slices"),
         ("surface", "Iso-surface"))
DEPTH_PRESETS = ((0.0, "Full model"), (500.0, "500 m"), (1000.0, "1 km"),
                 (2000.0, "2 km"), (5000.0, "5 km"), (10000.0, "10 km"),
                 (20000.0, "20 km"))


def _finite(v) -> bool:
    try:
        return v is not None and math.isfinite(float(v))
    except (TypeError, ValueError):
        return False


def file_elevations(path) -> tuple[dict[str, float], str]:
    """``({station: elevation}, how)`` from everything the file carries.

    *how* names the sources used ("station table", "topography block",
    "topography raster").  Zero-only tables (ModEM writes 0 when it has
    no topography) count as no topography.
    """
    from pycsamt.format.text import read_pcsf_or_pcsm

    pm = read_pcsf_or_pcsm(str(path))
    st = pm.stations
    if st is None:
        return {}, ""
    names = [str(n) for n in st.name]
    elev: dict[str, float] = {}
    used = []
    z = getattr(st, "z", None)
    if z is not None:
        for n, e in zip(names, np.asarray(z, dtype=float)):
            if _finite(e):
                elev[n] = float(e)
        if elev:
            used.append("station table")
    topo = getattr(pm, "topography", None)
    kind = type(topo).__name__ if topo is not None else ""
    if kind == "TopographyPerStation":
        before = len(elev)
        for n, e in zip(topo.station_id, np.asarray(topo.elevation, float)):
            if _finite(e) and str(n) not in elev:
                elev[str(n)] = float(e)
        if len(elev) > before:
            used.append("topography block")
    elif kind == "TopographyRaster":
        missing = [i for i, n in enumerate(names) if n not in elev]
        sx = np.asarray(getattr(st, "x", []), dtype=float)
        sy = np.asarray(getattr(st, "y", None) if getattr(st, "y", None)
                        is not None else np.zeros_like(sx), dtype=float)
        if missing and sx.size == len(names):
            got = _sample_raster(topo, sx[missing], sy[missing])
            for i, e in zip(missing, got):
                if _finite(e):
                    elev[names[i]] = float(e)
            if any(_finite(e) for e in got):
                used.append("topography raster")
    if elev and not any(abs(e) > 0 for e in elev.values()):
        return {}, ""  # an all-zero table is "no topography"
    return elev, " + ".join(used)


def _sample_raster(topo, x, y) -> np.ndarray:
    from scipy.interpolate import RegularGridInterpolator

    gx = np.asarray(topo.x, dtype=float)
    gy = np.asarray(topo.y, dtype=float)
    grid = np.asarray(topo.elevation, dtype=float)  # (ny, nx)
    if grid.shape != (gy.size, gx.size):
        grid = grid.T
    fn = RegularGridInterpolator((gy, gx), grid, bounds_error=False,
                                 fill_value=np.nan)
    return fn(np.column_stack([y, x]))


def fill_line_gaps(stations, elev: dict[str, float]) -> tuple[dict, int]:
    """Interpolate missing elevations along each line.

    *stations* are the view's station records (``id``, ``line``,
    ``latitude``/``longitude`` or ``x``/``y``).  Stations are ordered
    along each line by their principal direction; a line with fewer than
    two known elevations is left as it is.  Returns the completed mapping
    and how many stations were filled.
    """
    out = dict(elev)
    lines: dict[str, list] = {}
    for s in stations:
        lines.setdefault(str(getattr(s, "line", None) or "line"), []).append(s)
    filled = 0
    for recs in lines.values():
        pos = _along_line(recs)
        known = [(p, out[r.id]) for p, r in zip(pos, recs) if r.id in out]
        if len(known) < 2:
            continue
        kp, ke = zip(*sorted(known))
        for p, r in zip(pos, recs):
            if r.id not in out:
                out[r.id] = float(np.interp(p, kp, ke))
                filled += 1
    return out, filled


def _along_line(recs) -> np.ndarray:
    pts = []
    for r in recs:
        lat, lon = getattr(r, "latitude", None), getattr(r, "longitude", None)
        if _finite(lat) and _finite(lon):
            pts.append((float(lon) * math.cos(math.radians(float(lat))),
                        float(lat)))
        else:
            x, y = getattr(r, "x", None), getattr(r, "y", None)
            pts.append((float(x) if _finite(x) else np.nan,
                        float(y) if _finite(y) else 0.0))
    xy = np.asarray(pts, dtype=float)
    if len(xy) < 2 or not np.isfinite(xy).all():
        return np.arange(len(recs), dtype=float)
    c = xy - xy.mean(axis=0)
    _u, _s, vt = np.linalg.svd(c, full_matrices=False)
    return c @ vt[0]


_DENSITIES = (1.0, 0.75, 0.5, 0.25, 0.1, 0.05)
_VOL_SMOOTH = (0.0, 0.8, 1.5, 2.5, 4.0)
_SECTION_RES = (60, 100, 160, 240)


def _nearest(value, choices):
    return min(choices, key=lambda c: abs(float(c) - float(value)))


def _num_str(v) -> str:
    return f"{float(v):g}"


def mapview_controls(render: dict) -> dict:
    """Map View widget values for *render* (``map3d`` keyword names).

    The keys are the ones Map View's control gathering produces, and the
    values are in each widget's own format (select values are strings).
    """
    r = dict(render)
    c: dict = {}

    def put(key, value):
        if value is not None:
            c[key] = value

    put("mode3d", r.get("mode"))
    cmap = r.get("cmap")
    put("cmap", cmap if cmap in CMAPS else ("RdYlBu_r" if cmap else None))
    dr = r.get("depth_range")
    if dr:
        c["depth_lo"], c["depth_hi"] = float(dr[0]), float(dr[1])
    else:
        c["depth_lo"] = c["depth_hi"] = None
    for key in ("topography", "opacity", "show_stations", "station_labels",
                "station_symbol", "station_size", "station_color",
                "station_label_angle", "n_slices", "surface_count",
                "line_spacing", "azimuth", "x_unit", "depth_unit",
                "smooth_sections", "geology_legend", "geology_fill",
                "vertical_exaggeration"):
        put(key, r.get(key))
    put("terrain", r.get("show_terrain"))
    put("labels", r.get("show_labels"))
    put("contours", r.get("show_contours"))
    if "aspectmode" in r:
        c["aspect"] = r["aspectmode"]
    if "max_stations" in r:
        c["station_max"] = r["max_stations"]
    if r.get("station_label_fraction") is not None:
        c["station_label_density"] = _num_str(
            _nearest(r["station_label_fraction"], _DENSITIES))
    names = r.get("station_label_names")
    c["station_label_names"] = ", ".join(names) if names else ""
    if "log_color" in r:
        c["scale"] = "log" if r["log_color"] else "linear"
    vr = r.get("value_range")
    c["vmin"], c["vmax"] = (vr if vr else (None, None))
    cp = r.get("crange_percentile")
    if cp:
        c["crange_plo"], c["crange_phi"] = float(cp[0]), float(cp[1])
    rr = r.get("rho_range")
    c["rho_lo"], c["rho_hi"] = (rr if rr else (None, None))
    c["rho_cutoff"] = r.get("rho_display_max")
    if r.get("section_res") is not None:
        c["section_res"] = str(_nearest(r["section_res"], _SECTION_RES))
    if r.get("volume_smoothing") is not None:
        c["volume_smoothing"] = _num_str(
            _nearest(r["volume_smoothing"], _VOL_SMOOTH))
    return c


def mapview_state(render: dict, *, elevations: dict | None = None,
                  spin: bool = False, dark: bool = False,
                  camera: dict | None = None, geology=None,
                  boreholes=None, structure=None,
                  source: str = "") -> dict:
    """The Map View seed state for the desktop scene.

    *render* holds the scene's ``map3d`` options (mode, cmap, depth,
    stations, labels, resistivity, view, ...); see
    :func:`mapview_controls`.
    """
    state = {
        "version": 1,
        "source": source,
        "view": "map3d",
        "controls": mapview_controls(render),
        "elevations": {str(k): float(v) for k, v in (elevations or {})
                       .items() if _finite(v)},
        "spin": bool(spin),
        "theme": "dark" if dark else "light",
        "camera": camera or None,
    }
    if geology is not None:
        from pycsamt.app._geology import store_from_legend

        state["geo"] = store_from_legend(geology)
    if boreholes is not None:
        from pycsamt.format.borehole import pcbh_to_dict

        state["pcbh"] = {"filename": "boreholes.pcbh.json",
                         "document": pcbh_to_dict(boreholes),
                         "n_boreholes": len(boreholes.boreholes)}
    if structure is not None:
        from pycsamt.app._structure import store_from_structure

        state["struct"] = store_from_structure(structure)
    return state
