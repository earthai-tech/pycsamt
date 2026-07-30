"""Generate the static figures used by the map volume guide."""

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.patches import Patch
from scipy.interpolate import griddata
from skimage.measure import marching_cubes

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.map import VolumeMapOptions, load_lines
from pycsamt.map.volume import _profile_grids, _values_at_depth
from pycsamt.models.modem.results import InversionResult


IMAGES = ROOT / "docs/source/images/user_guide/map"
MODEM = ROOT / "data/modem/willy_27freq_watex_line02_sample"
EDI = ROOT / "data/AMT/WILLY_DATA"


def make_separated_depth_slices() -> None:
    """Render depth slices as independent panels with explicit depths."""
    data = load_lines(EDI, detect="folder", recursive=True)
    options = VolumeMapOptions(
        mode="depth",
        quantity="resistivity",
        component="xy",
        depth_range=(0.0, 3000.0),
        log_color=True,
    )
    profiles = _profile_grids(data, options)
    depths = (250.0, 750.0, 1500.0, 2500.0)
    finite = []
    prepared = []
    for depth in depths:
        rows = []
        for line, grid in profiles.items():
            values = _values_at_depth(grid, depth, options)
            rows.append((line, np.asarray(grid["x"]), values))
            finite.extend(values[np.isfinite(values)])
        prepared.append(rows)

    norm = colors.Normalize(
        vmin=np.log10(np.nanpercentile(finite, 2)),
        vmax=np.log10(np.nanpercentile(finite, 98)),
    )
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 7.2), constrained_layout=True)
    levels = np.linspace(norm.vmin, norm.vmax, 15)
    for ax, depth, rows in zip(axes.flat, depths, prepared):
        px, py, pv = [], [], []
        for line_index, (_line, distance, values) in enumerate(rows):
            valid = np.isfinite(values) & (values > 0)
            px.extend(distance[valid] / 1000.0)
            py.extend(np.full(valid.sum(), line_index))
            pv.extend(np.log10(values[valid]))
        px, py, pv = map(np.asarray, (px, py, pv))
        gx = np.linspace(px.min(), px.max(), 240)
        gy = np.linspace(py.min(), py.max(), 150)
        GX, GY = np.meshgrid(gx, gy)
        # Linear interpolation respects observed gradients; nearest values
        # only close small internal holes, while the convex-hull exterior
        # remains masked.
        GZ = griddata((px, py), pv, (GX, GY), method="linear")
        nearest = griddata((px, py), pv, (GX, GY), method="nearest")
        GZ = np.where(np.isfinite(GZ), GZ, nearest)
        ax.contourf(
            GX, GY, GZ, levels=levels, cmap="turbo", norm=norm,
            extend="both",
        )
        ax.contour(
            GX, GY, GZ, levels=levels[::2], colors="#172554",
            linewidths=0.45, alpha=0.55,
        )
        for line_index, (_line, distance, values) in enumerate(rows):
            valid = np.isfinite(values)
            ax.plot(
                distance[valid] / 1000.0,
                np.full(valid.sum(), line_index),
                color="white", linestyle=":", linewidth=1.15,
                alpha=0.95,
            )
            ax.scatter(
                distance[valid] / 1000.0,
                np.full(valid.sum(), line_index),
                s=5, color="#0f172a", alpha=0.65, linewidths=0,
            )
        ax.set_title(f"Depth = {depth:,.0f} m", loc="left", weight="bold")
        ax.set_xlabel("Profile distance (km)")
        ax.set_ylabel("Survey line")
        ax.set_yticks(range(len(rows)), [row[0] for row in rows])
        ax.grid(color="#94a3b8", linestyle=":", linewidth=0.7, alpha=0.45)
        ax.set_facecolor("#f8fafc")
    sm = plt.cm.ScalarMappable(norm=norm, cmap="turbo")
    cbar = fig.colorbar(sm, ax=axes, shrink=0.88, pad=0.02)
    cbar.set_label(r"$\log_{10}(\rho_a\;[\Omega\,m])$")
    fig.suptitle("Independent apparent-resistivity depth slices", weight="bold")
    fig.savefig(IMAGES / "map_volume_depth_slices_preview.png", dpi=190)
    plt.close(fig)


def _modem_model():
    result = InversionResult(str(MODEM), verbose=0)
    model = result.model_final
    rho = np.asarray(model.rho_linear, dtype=float)
    x = 0.5 * (np.asarray(model.x_nodes[:-1]) + np.asarray(model.x_nodes[1:]))
    y = 0.5 * (np.asarray(model.y_nodes[:-1]) + np.asarray(model.y_nodes[1:]))
    z = 0.5 * (np.asarray(model.z_nodes[:-1]) + np.asarray(model.z_nodes[1:]))
    return result, rho, x, y, z


def make_modem_threshold_blocks() -> None:
    """Show conductive and resistive cells from the final ModEM model."""
    result, rho, x, y, z = _modem_model()
    # Remove the laterally padded boundary and deepest coarse cells, then
    # decimate the long x direction to keep individual voxels legible.
    volume = rho[:28, 5:-5, 18:-18:4]
    xx = x[18:-18:4] / 1000.0
    yy = y[5:-5] / 1000.0
    zz = z[:28] / 1000.0
    X, Y, Z = np.meshgrid(xx, yy, zz, indexing="xy")
    station_lines = {}
    x_origin = 0.5 * (x[0] + x[-1]) / 1000.0
    y_origin = 0.5 * (y[0] + y[-1]) / 1000.0
    for name, (east, north, _height) in result.data_obs.site_coords.items():
        line = name.split("-")[1]
        station_lines.setdefault(line, []).append(
            (east / 1000.0 + x_origin, north / 1000.0 + y_origin)
        )

    fig = plt.figure(figsize=(13.0, 5.8), constrained_layout=True)
    cases = ((30.0, "Conductive cells: $\\rho \\leq 30$ $\\Omega$ m", "#06b6d4"),
             (1000.0, "Resistive cells: $\\rho \\geq 1000$ $\\Omega$ m", "#ef4444"))
    for index, (threshold, title, color) in enumerate(cases, start=1):
        ax = fig.add_subplot(1, 2, index, projection="3d")
        mask = volume <= threshold if index == 1 else volume >= threshold
        zi, yi, xi = np.nonzero(mask)
        values = volume[mask]
        keep = np.arange(values.size) % max(1, values.size // 3500 + 1) == 0
        ax.scatter(
            xx[xi[keep]], yy[yi[keep]], -zz[zi[keep]],
            c=color, s=7, marker="s", alpha=0.32, linewidths=0,
        )
        for line_index, (line, points) in enumerate(sorted(station_lines.items())):
            points = np.asarray(points)
            order = np.argsort(points[:, 0])
            points = points[order]
            ax.plot(
                points[:, 0], points[:, 1], np.full(points.shape[0], 0.025),
                color="#111827", linewidth=1.0, marker="o", markersize=1.8,
                alpha=0.9,
            )
            endpoint = -1 if line_index % 2 == 0 else 0
            ax.text(
                points[endpoint, 0], points[endpoint, 1],
                0.045 + 0.006 * line_index, f"L{line}",
                fontsize=7.5, weight="bold", color="#0f172a",
                bbox=dict(facecolor="white", alpha=0.72,
                          edgecolor="none", pad=1.0),
            )
        ax.set_title(title, weight="bold")
        ax.set_xlabel("East (km)")
        ax.set_ylabel("North (km)")
        ax.set_zlabel("Elevation relative to surface (km)")
        ax.grid(True, linestyle=":", linewidth=0.65, alpha=0.32)
        ax.xaxis.pane.set_alpha(0.0)
        ax.yaxis.pane.set_alpha(0.0)
        ax.zaxis.pane.set_alpha(0.0)
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis._axinfo["grid"].update(
                color=(0.39, 0.45, 0.55, 0.35),
                linestyle=":",
                linewidth=0.65,
            )
        ax.view_init(elev=25, azim=-58)
    fig.suptitle(
        f"Final ModEM inversion: thresholded model cells (RMS {result.final_rms:.3f})",
        weight="bold",
    )
    fig.savefig(IMAGES / "map_volume_modem_threshold_blocks.png", dpi=190)
    plt.close(fig)


def make_phase_fence() -> None:
    """Render phase curtains with station tracks and line identifiers."""
    data = load_lines(EDI, detect="folder", recursive=True)
    options = VolumeMapOptions(
        mode="fence", quantity="phase", component="xy",
        depth_range=(0.0, 3000.0), show_stations=True,
    )
    profiles = _profile_grids(data, options)
    norm = colors.Normalize(vmin=-180.0, vmax=180.0)
    cmap = plt.get_cmap("RdBu_r")
    fig = plt.figure(figsize=(10.8, 8.2), constrained_layout=True)
    ax = fig.add_subplot(111, projection="3d")
    for line_index, (line, grid) in enumerate(profiles.items()):
        distance = np.asarray(grid["x"], dtype=float)
        depth = np.asarray(grid["z"], dtype=float)
        phase = np.asarray(grid["value"], dtype=float)
        keep = depth <= 3000.0
        distance, depth, phase = distance, depth[keep], phase[keep]
        XX, ZZ = np.meshgrid(distance, -depth)
        offset = float(line_index * 200.0)
        YY = np.full_like(XX, offset)
        ax.plot_surface(
            XX, YY, ZZ, facecolors=cmap(norm(phase)),
            shade=False, alpha=0.82, linewidth=0,
        )
        ax.plot(
            distance, np.full(distance.size, offset), np.zeros(distance.size),
            color="#111827", linewidth=1.0, marker="o", markersize=2.4,
        )
        endpoint = -1 if line_index % 2 == 0 else 0
        ax.text(
            distance[endpoint], offset, 160.0 + 38.0 * line_index,
            line, fontsize=8, weight="bold", color="#0f172a",
            ha="left" if endpoint == -1 else "right",
            bbox=dict(facecolor="white", alpha=0.76,
                      edgecolor="none", pad=1.1),
        )
    ax.set_title("Phase fence: geometry held fixed", weight="bold")
    ax.set_xlabel("Profile distance (m)")
    ax.set_ylabel("Survey line")
    ax.set_yticks(
        [float(index * 200.0) for index in range(len(profiles))],
        list(profiles),
    )
    ax.set_zlabel("Pseudo-depth (m)")
    ax.set_zlim(-3000.0, 400.0)
    ax.view_init(elev=24, azim=-61)
    ax.grid(True, linestyle=":", linewidth=0.6, alpha=0.3)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_alpha(0.0)
        axis._axinfo["grid"].update(
            color=(0.39, 0.45, 0.55, 0.34),
            linestyle=":", linewidth=0.6,
        )
    scalar = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(scalar, ax=ax, shrink=0.65, pad=0.08)
    cbar.set_label("Phase (deg)")
    fig.savefig(IMAGES / "map_volume_phase_fence.png", dpi=190)
    plt.close(fig)


def make_modem_isosurfaces() -> None:
    """Extract conductive and resistive boundaries with marching cubes."""
    result, rho, x, y, z = _modem_model()
    log_rho = np.log10(rho[:32, 4:-4, 12:-12])
    spacing = (
        float(np.median(np.diff(z[:32]))) / 1000.0,
        float(np.median(np.diff(y[4:-4]))) / 1000.0,
        float(np.median(np.diff(x[12:-12]))) / 1000.0,
    )
    fig = plt.figure(figsize=(10.2, 7.2), constrained_layout=True)
    ax = fig.add_subplot(111, projection="3d")
    for level, color, label in (
        (np.log10(30.0), "#0891b2", r"30 $\Omega$ m conductor boundary"),
        (np.log10(1000.0), "#dc2626", r"1000 $\Omega$ m resistor boundary"),
    ):
        verts, faces, _, _ = marching_cubes(log_rho, level=level, spacing=spacing)
        # marching_cubes returns coordinates in (z, y, x) order.
        ax.plot_trisurf(
            verts[:, 2], verts[:, 1], -verts[:, 0],
            triangles=faces, color=color, alpha=0.28,
            linewidth=0.05, edgecolor=color,
        )
    x_origin = 0.5 * (x[0] + x[-1]) / 1000.0 - x[12] / 1000.0
    y_origin = 0.5 * (y[0] + y[-1]) / 1000.0 - y[4] / 1000.0
    station_lines = {}
    for name, (east, north, _height) in result.data_obs.site_coords.items():
        line = name.split("-")[1]
        station_lines.setdefault(line, []).append(
            (east / 1000.0 + x_origin, north / 1000.0 + y_origin)
        )
    for line_index, (line, points) in enumerate(sorted(station_lines.items())):
        points = np.asarray(points)
        order = np.argsort(points[:, 0])
        points = points[order]
        ax.plot(
            points[:, 0], points[:, 1], np.full(points.shape[0], 0.018),
            color="#111827", linewidth=1.05, marker="o", markersize=2.0,
            alpha=0.92,
        )
        endpoint = -1 if line_index % 2 == 0 else 0
        ax.text(
            points[endpoint, 0], points[endpoint, 1],
            0.035 + 0.009 * line_index, f"L{line}",
            fontsize=8, weight="bold", color="#0f172a",
            ha="left" if endpoint == -1 else "right",
            bbox=dict(facecolor="white", alpha=0.76,
                      edgecolor="none", pad=1.0),
        )
    ax.set_xlabel("East (km)")
    ax.set_ylabel("North (km)")
    ax.set_zlabel("Elevation relative to surface (km)")
    ax.set_title(
        f"Final ModEM inversion: selected resistivity isosurfaces\nRMS {result.final_rms:.3f}",
        weight="bold",
    )
    ax.grid(True, linestyle=":", linewidth=0.65, alpha=0.3)
    ax.xaxis.pane.set_alpha(0.0)
    ax.yaxis.pane.set_alpha(0.0)
    ax.zaxis.pane.set_alpha(0.0)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis._axinfo["grid"].update(
            color=(0.39, 0.45, 0.55, 0.35),
            linestyle=":",
            linewidth=0.65,
        )
    ax.view_init(elev=24, azim=-58)
    ax.legend(
        handles=[
            Patch(color="#0891b2", alpha=0.45,
                  label=r"30 $\Omega$ m conductor boundary"),
            Patch(color="#dc2626", alpha=0.45,
                  label=r"1000 $\Omega$ m resistor boundary"),
        ],
        loc="upper left",
    )
    fig.savefig(IMAGES / "map_volume_modem_isosurfaces.png", dpi=190)
    plt.close(fig)


def _light_3d_frame(ax) -> None:
    """Apply the transparent, dotted frame shared by static 3-D views."""
    ax.grid(True, linestyle=":", linewidth=0.55, alpha=0.3)
    ax.xaxis.pane.set_alpha(0.0)
    ax.yaxis.pane.set_alpha(0.0)
    ax.zaxis.pane.set_alpha(0.0)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis._axinfo["grid"].update(
            color=(0.39, 0.45, 0.55, 0.32),
            linestyle=":",
            linewidth=0.55,
        )
    ax.set_xlabel("East (km)", labelpad=2)
    ax.set_ylabel("North (km)", labelpad=2)
    ax.set_zlabel("Elevation (km)", labelpad=2)
    ax.view_init(elev=25, azim=-58)


def make_modem_modes_overview() -> None:
    """Compare four views of the same final ModEM inversion model."""
    result, rho, x, y, z = _modem_model()
    # Crop padded boundary cells and the deep, very coarse mesh.  Moderate
    # decimation keeps the four-panel figure readable without changing the
    # threshold definitions.
    xs = x[18:-18:4] / 1000.0
    ys = y[5:-5:2] / 1000.0
    zs = z[:28] / 1000.0
    vol = rho[:28, 5:-5:2, 18:-18:4]
    log_vol = np.log10(vol)
    norm = colors.Normalize(0.0, 4.0)
    cmap = plt.get_cmap("turbo")

    fig = plt.figure(figsize=(13.2, 10.0), constrained_layout=True)

    # Fence: three genuine north-indexed sections through the inversion mesh.
    ax = fig.add_subplot(2, 2, 1, projection="3d")
    curtain_ax = ax
    XX, ZZ = np.meshgrid(xs, -zs)
    line_positions = {}
    for name, (_east, north, _height) in result.data_obs.site_coords.items():
        line = name.split("-")[1]
        line_positions.setdefault(line, []).append(north / 1000.0)
    model_north_origin = 0.5 * (y[0] + y[-1]) / 1000.0
    for line_index, (line, northings) in enumerate(
        sorted(line_positions.items(), key=lambda item: item[0])
    ):
        line_y = float(np.mean(northings)) + model_north_origin
        yi = int(np.argmin(np.abs(ys - line_y)))
        YY = np.full_like(XX, ys[yi])
        ax.plot_surface(
            XX, YY, ZZ, facecolors=cmap(norm(log_vol[:, yi, :])),
            shade=False, alpha=0.88, linewidth=0,
        )
        label_x = xs.max() + 0.04 if line_index % 2 == 0 else xs.min() - 0.04
        ax.text(
            label_x, ys[yi], -zs[0] - 0.012 * line_index, f"L{line}",
            fontsize=8, weight="bold", color="#0f172a",
            ha="left" if line_index % 2 == 0 else "right",
            bbox=dict(facecolor="white", alpha=0.72, edgecolor="none", pad=1.2),
        )
    ax.set_title("Inversion curtains", weight="bold")
    _light_3d_frame(ax)

    # Block: retain only interpretable end members instead of an opaque cube.
    ax = fig.add_subplot(2, 2, 2, projection="3d")
    for mask, color, label in (
        (vol <= 30.0, "#06b6d4", r"$\rho\leq30$ $\Omega$ m"),
        (vol >= 1000.0, "#ef4444", r"$\rho\geq1000$ $\Omega$ m"),
    ):
        zi, yi, xi = np.nonzero(mask)
        stride = max(1, zi.size // 1800 + 1)
        ax.scatter(
            xs[xi[::stride]], ys[yi[::stride]], -zs[zi[::stride]],
            s=5, marker="s", color=color, alpha=0.28,
            linewidths=0, label=label,
        )
    ax.set_title("Thresholded model cells", weight="bold")
    ax.legend(loc="upper left", fontsize=8)
    _light_3d_frame(ax)

    # Depth: independent horizontal surfaces with a small vertical gap.
    ax = fig.add_subplot(2, 2, 3, projection="3d")
    depth_ax = ax
    XI, YI = np.meshgrid(xs, ys)
    for zi in (5, 12, 20):
        plane = log_vol[zi]
        ax.plot_surface(
            XI, YI, np.full_like(XI, -zs[zi]),
            facecolors=cmap(norm(plane)), shade=False,
            alpha=0.86, linewidth=0,
        )
        ax.text(
            xs.max() + 0.04, ys.max(), -zs[zi],
            f"depth {zs[zi]:.2f} km", fontsize=8, weight="bold",
            color="#0f172a",
            bbox=dict(facecolor="white", alpha=0.78, edgecolor="none", pad=1.2),
        )
    ax.set_title("Independent model-depth slices", weight="bold")
    _light_3d_frame(ax)

    # Surface: explicit conductive and resistive boundaries.
    ax = fig.add_subplot(2, 2, 4, projection="3d")
    spacing = (
        float(np.median(np.diff(zs))),
        float(np.median(np.diff(ys))),
        float(np.median(np.diff(xs))),
    )
    for level, color in ((np.log10(30.0), "#0891b2"),
                         (np.log10(1000.0), "#dc2626")):
        verts, faces, _, _ = marching_cubes(log_vol, level=level, spacing=spacing)
        ax.plot_trisurf(
            verts[:, 2] + xs.min(), verts[:, 1] + ys.min(), -verts[:, 0],
            triangles=faces, color=color, alpha=0.32,
            linewidth=0.04, edgecolor=color,
        )
    ax.set_title("30 and 1000 $\Omega$ m isosurfaces", weight="bold")
    ax.legend(
        handles=[
            Patch(color="#0891b2", alpha=0.45,
                  label=r"30 $\Omega$ m conductor boundary"),
            Patch(color="#dc2626", alpha=0.45,
                  label=r"1000 $\Omega$ m resistor boundary"),
        ],
        loc="upper left", fontsize=8,
    )
    _light_3d_frame(ax)

    scalar = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(
        scalar, ax=[curtain_ax, depth_ax], shrink=0.7,
        pad=0.02, location="left",
    )
    cbar.set_label(r"$\log_{10}(\rho\;[\Omega\,m])$")
    fig.suptitle(
        f"Final ModEM inversion rendered four ways (RMS {result.final_rms:.3f})",
        fontsize=17, weight="bold",
    )
    fig.savefig(IMAGES / "map_volume_modes_overview.png", dpi=190)
    plt.close(fig)


def make_modem_topography_overlay() -> None:
    """Drape observed elevations above bodies from the final ModEM model."""
    result, rho, x, y, z = _modem_model()
    edi = load_lines(EDI, detect="folder", recursive=True)
    elevation_by_id = {
        station.id: station.elevation
        for station in edi.stations
        if station.elevation is not None
    }
    modem_data = result.data_obs
    samples = []
    for name, (east, north, _height) in modem_data.site_coords.items():
        short_name = name.split("-", 1)[-1]
        elevation = elevation_by_id.get(short_name)
        if elevation is not None:
            samples.append((name, east / 1000.0, north / 1000.0,
                            elevation / 1000.0))
    names, sx, sy, sh = zip(*samples)
    sx, sy, sh = map(np.asarray, (sx, sy, sh))

    gx = np.linspace(sx.min(), sx.max(), 100)
    gy = np.linspace(sy.min(), sy.max(), 90)
    GX, GY = np.meshgrid(gx, gy)
    GH = griddata((sx, sy), sh, (GX, GY), method="linear")
    nearest = griddata((sx, sy), sh, (GX, GY), method="nearest")
    GH = np.where(np.isfinite(GH), GH, nearest)

    # ModEM nodes are stored from zero; station offsets are centred on the
    # model origin. Shift cell centres into that same coordinate frame.
    mx = (x - 0.5 * (x[0] + x[-1])) / 1000.0
    my = (y - 0.5 * (y[0] + y[-1])) / 1000.0
    xi = np.flatnonzero((mx >= sx.min() - 0.25) & (mx <= sx.max() + 0.25))
    yi = np.flatnonzero((my >= sy.min() - 0.25) & (my <= sy.max() + 0.25))
    zi = np.arange(min(31, len(z)))
    sub = np.log10(rho[np.ix_(zi, yi, xi)])
    spacing = (
        float(np.median(np.diff(z[zi]))) / 1000.0,
        float(np.median(np.diff(my[yi]))),
        float(np.median(np.diff(mx[xi]))),
    )

    fig = plt.figure(figsize=(12.2, 8.2), constrained_layout=True)
    ax = fig.add_subplot(111, projection="3d")
    terrain = ax.plot_surface(
        GX, GY, GH, cmap="terrain", alpha=0.58,
        linewidth=0, antialiased=True, shade=True,
    )
    surface_level = float(np.nanmedian(sh))
    for level, color, alpha in (
        (np.log10(30.0), "#06b6d4", 0.34),
        (np.log10(1000.0), "#dc2626", 0.18),
    ):
        verts, faces, _, _ = marching_cubes(sub, level=level, spacing=spacing)
        ax.plot_trisurf(
            verts[:, 2] + mx[xi[0]],
            verts[:, 1] + my[yi[0]],
            surface_level - verts[:, 0],
            triangles=faces, color=color, alpha=alpha,
            linewidth=0.04, edgecolor=color,
        )

    line_names = np.array([name.split("-")[1] for name in names])
    for line in sorted(set(line_names)):
        keep = line_names == line
        order = np.argsort(sx[keep])
        lx, ly, lh = sx[keep][order], sy[keep][order], sh[keep][order]
        ax.plot(lx, ly, lh + 0.006, color="#111827", linewidth=1.25,
                marker="o", markersize=2.3, zorder=20)
        ax.text(lx[-1], ly[-1], lh[-1] + 0.018, f"L{line}",
                fontsize=8, weight="bold")

    ax.set_title(
        "Observed topography over the final ModEM inversion\n"
        "cyan: $\\rho\leq30$ $\\Omega$ m; red: $\\rho\geq1000$ $\\Omega$ m",
        weight="bold",
    )
    ax.set_xlabel("Relative east (km)")
    ax.set_ylabel("Relative north (km)")
    ax.set_zlabel("Elevation relative to datum (km)")
    ax.view_init(elev=27, azim=-57)
    ax.grid(True, linestyle=":", linewidth=0.6, alpha=0.3)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_alpha(0.0)
        axis._axinfo["grid"].update(
            color=(0.39, 0.45, 0.55, 0.32), linestyle=":", linewidth=0.6,
        )
    cbar = fig.colorbar(terrain, ax=ax, shrink=0.62, pad=0.08)
    cbar.set_label("Station-interpolated elevation (km)")
    fig.savefig(IMAGES / "map_overlays_topography_surface.png", dpi=190)
    plt.close(fig)


if __name__ == "__main__":
    IMAGES.mkdir(parents=True, exist_ok=True)
    make_separated_depth_slices()
    make_modem_threshold_blocks()
    make_phase_fence()
    make_modem_isosurfaces()
    make_modem_modes_overview()
    make_modem_topography_overlay()
