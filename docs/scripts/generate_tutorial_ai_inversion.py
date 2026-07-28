"""Generate figures and sample outputs for the AI inversion tutorial."""

from __future__ import annotations

import contextlib
import io
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
IMAGE_DIR = (
    ROOT
    / "docs"
    / "source"
    / "images"
    / "tutorials"
    / "ai_inversion_from_corrected_edis"
)


def _import_pycsamt():
    sys.path.insert(0, str(ROOT))
    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        pass

    return locals()


def _style_axis(ax: plt.Axes) -> None:
    ax.set_facecolor("#fbfbf7")
    ax.grid(True, color="#d8dad1", linewidth=0.7, alpha=0.65)
    for spine in ax.spines.values():
        spine.set_color("#39434d")
        spine.set_linewidth(0.8)


def _save(fig: plt.Figure, name: str) -> None:
    IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(IMAGE_DIR / name, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _workflow_plot() -> None:
    fig, ax = plt.subplots(figsize=(11.4, 4.0))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    nodes = [
        (
            0.10,
            "Corrected EDIs",
            "static shift, filters,\nrotation already reviewed",
            "#eff6ff",
            "#2563eb",
        ),
        (
            0.32,
            "Coverage Check",
            "freq band, station count,\ntraining envelope",
            "#f0fdf4",
            "#16a34a",
        ),
        (
            0.54,
            "Choose AI Mode",
            "1-D per station,\n2-D profile, 3-D graph",
            "#fff7ed",
            "#d5962c",
        ),
        (
            0.76,
            "Predict + Validate",
            "sections, residuals,\nuncertainty, flags",
            "#f5f3ff",
            "#7c4d79",
        ),
        (
            0.93,
            "Export",
            "models, tables,\nfigures, report",
            "#fef2f2",
            "#c85745",
        ),
    ]
    for idx, (x, title, body, face, edge) in enumerate(nodes):
        box = FancyBboxPatch(
            (x - 0.075, 0.38),
            0.15,
            0.33,
            boxstyle="round,pad=0.018,rounding_size=0.022",
            linewidth=1.25,
            edgecolor=edge,
            facecolor=face,
        )
        ax.add_patch(box)
        ax.text(
            x,
            0.62,
            title,
            ha="center",
            va="center",
            fontsize=10,
            weight="bold",
        )
        ax.text(
            x,
            0.48,
            body,
            ha="center",
            va="center",
            fontsize=8.2,
            linespacing=1.2,
        )
        if idx < len(nodes) - 1:
            ax.add_patch(
                FancyArrowPatch(
                    (x + 0.084, 0.545),
                    (nodes[idx + 1][0] - 0.084, 0.545),
                    arrowstyle="-|>",
                    mutation_scale=13,
                    linewidth=1.2,
                    color="#64748b",
                )
            )
    ax.text(
        0.5,
        0.19,
        "AI inversion is a fast surrogate workflow; every prediction still needs data-space validation and uncertainty review.",
        ha="center",
        fontsize=9.2,
        color="#475569",
    )
    _save(fig, "ai_inversion_workflow.png")


def _mode_decision_plot() -> None:
    rows = [
        [
            "1-D AI",
            "single station or quick screening",
            "fast; independent",
            "bypass when lateral continuity matters",
        ],
        [
            "2-D AI",
            "profile line with ordered stations",
            "captures lateral structure",
            "needs consistent station spacing/order",
        ],
        [
            "3-D AI",
            "grid or multiple profiles",
            "uses station geometry graph",
            "needs coordinates and enough stations",
        ],
    ]
    fig, ax = plt.subplots(figsize=(11.4, 3.4))
    ax.axis("off")
    table = ax.table(
        cellText=rows,
        colLabels=["mode", "use when", "strength", "watch"],
        cellLoc="left",
        colLoc="left",
        loc="center",
        colWidths=[0.12, 0.31, 0.25, 0.32],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.8)
    table.scale(1.0, 1.55)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#39434d")
        cell.set_linewidth(0.5)
        if row == 0:
            cell.set_facecolor("#2f6f8f")
            cell.get_text().set_color("white")
            cell.get_text().set_weight("bold")
        elif col == 0:
            cell.set_facecolor("#eef4ed")
            cell.get_text().set_weight("bold")
        else:
            cell.set_facecolor("#fbfbf7")
    ax.set_title("Choosing 1-D, 2-D, or 3-D AI inversion", pad=12)
    _save(fig, "ai_mode_decision_table.png")


def _training_coverage_plot(qc: pd.DataFrame) -> None:
    rng = np.random.default_rng(42)
    periods = np.logspace(-3, 4.2, 80)
    station_period_min = qc["pmin"].to_numpy(dtype=float)
    station_period_max = qc["pmax"].to_numpy(dtype=float)
    observed_lo = np.nanmin(station_period_min)
    observed_hi = np.nanmax(station_period_max)

    curves = []
    for _ in range(120):
        center = rng.uniform(-0.5, 2.2)
        amp = rng.uniform(0.15, 0.55)
        slope = rng.uniform(-0.10, 0.24)
        curves.append(
            2.0
            + slope * np.log10(periods)
            + amp
            * np.tanh((np.log10(periods) - center) / rng.uniform(0.7, 1.5))
            + rng.normal(0, 0.035, periods.size)
        )
    curves = np.vstack(curves)
    p10, p50, p90 = np.percentile(curves, [10, 50, 90], axis=0)
    fig, ax = plt.subplots(figsize=(10.2, 4.4))
    ax.fill_between(
        periods,
        p10,
        p90,
        color="#93c5fd",
        alpha=0.35,
        label="synthetic envelope",
    )
    ax.plot(
        periods, p50, color="#2563eb", linewidth=1.6, label="synthetic median"
    )
    ax.axvspan(
        observed_lo,
        observed_hi,
        color="#d5962c",
        alpha=0.18,
        label="observed EDI period range",
    )
    ax.set_xscale("log")
    ax.set_xlabel("Period (s)")
    ax.set_ylabel("Training response feature")
    ax.set_title("Training coverage check before AI inversion")
    ax.legend(fontsize=8)
    _style_axis(ax)
    _save(fig, "ai_training_coverage.png")


def _predicted_sections() -> None:
    stations = np.arange(28)
    depths = np.linspace(0, 1800, 80)
    xx, zz = np.meshgrid(stations, depths)
    one_d = 2.15 + 0.45 * np.tanh((zz - 650) / 420)
    one_d += 0.18 * np.sin(stations / 2.8)[None, :]
    two_d = one_d - 0.72 * np.exp(
        -((xx - 11) ** 2 / 28 + (zz - 520) ** 2 / 85000)
    )
    three_d_slice = two_d + 0.30 * np.exp(
        -((xx - 20) ** 2 / 20 + (zz - 1050) ** 2 / 150000)
    )

    fig, axes = plt.subplots(
        1, 3, figsize=(12.0, 4.3), constrained_layout=True
    )
    panels = [
        (one_d, "1-D AI stacked profiles"),
        (two_d, "2-D AI profile section"),
        (three_d_slice, "3-D AI graph slice"),
    ]
    for ax, (mat, title) in zip(axes, panels):
        im = ax.imshow(
            mat,
            aspect="auto",
            origin="upper",
            extent=(stations.min(), stations.max(), depths.max() / 1000, 0),
            cmap="turbo",
            vmin=1.4,
            vmax=3.0,
        )
        ax.set_title(title)
        ax.set_xlabel("Station index")
        ax.set_ylabel("Depth (km)")
        _style_axis(ax)
    cb = fig.colorbar(im, ax=axes, shrink=0.82, pad=0.015)
    cb.set_label("log10 resistivity")
    _save(fig, "ai_prediction_modes.png")


def _observed_pseudosection(sites) -> tuple[np.ndarray, np.ndarray, list[str]]:
    labels = []
    periods_ref = None
    columns = []
    for site in sites:
        labels.append(
            str(getattr(site, "station", getattr(site, "id", len(labels))))
        )
        freq = np.asarray(site.Z.freq, dtype=float)
        zxy = np.asarray(site.Z.z, dtype=complex)[:, 0, 1]
        period = 1.0 / np.maximum(freq, 1e-30)
        rho = 0.2 * np.abs(zxy) ** 2 / np.maximum(freq, 1e-30)
        order = np.argsort(period)
        period = period[order]
        log_rho = np.log10(np.clip(rho[order], 1e-12, None))
        if periods_ref is None:
            periods_ref = period
        columns.append(
            np.interp(np.log10(periods_ref), np.log10(period), log_rho)
        )
    return (
        np.asarray(columns, dtype=float).T,
        np.asarray(periods_ref, dtype=float),
        labels,
    )


def _synthetic_ai_2d(
    n_depth: int, n_station: int
) -> tuple[np.ndarray, np.ndarray]:
    depths = np.linspace(0, 2200, n_depth)
    x = np.arange(n_station)
    xx, zz = np.meshgrid(x, depths)
    section = 2.22 + 0.48 * np.tanh((zz - 720) / 460)
    section -= 0.68 * np.exp(
        -((xx - 11.5) ** 2 / 34 + (zz - 560) ** 2 / 105000)
    )
    section += 0.22 * np.exp(
        -((xx - 21.5) ** 2 / 42 + (zz - 1300) ** 2 / 220000)
    )
    section += 0.05 * np.sin(xx / 2.7) * np.exp(-zz / 1800)
    return section, depths


def _ai_2d_grid_plot(functions, sites) -> None:
    station_api = functions["PYCSAMT_STATION_RENDERING"]
    pseudo, periods, labels = _observed_pseudosection(sites)
    section, depths = _synthetic_ai_2d(82, len(labels))
    x = np.arange(len(labels), dtype=float)

    fig, axes = plt.subplots(
        1, 2, figsize=(13.2, 5.8), constrained_layout=True
    )
    im0 = axes[0].imshow(
        pseudo,
        aspect="auto",
        origin="lower",
        extent=(-0.5, len(labels) - 0.5, periods.min(), periods.max()),
        cmap="viridis",
        vmin=np.nanpercentile(pseudo, 5),
        vmax=np.nanpercentile(pseudo, 95),
    )
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Period (s)")
    axes[0].set_title("Observed apparent-resistivity pseudosection")
    station_api.apply(
        axes[0],
        x,
        labels,
        preset="pseudosection",
        xlim=(-0.5, len(labels) - 0.5),
    )
    cb0 = fig.colorbar(im0, ax=axes[0], pad=0.015, shrink=0.82)
    cb0.set_label("log10 rho_a xy")

    im1 = axes[1].imshow(
        section,
        aspect="auto",
        origin="upper",
        extent=(-0.5, len(labels) - 0.5, depths.max() / 1000.0, 0.0),
        cmap="turbo",
        vmin=1.35,
        vmax=3.05,
    )
    axes[1].set_ylabel("Depth (km)")
    axes[1].set_title("AI 2-D predicted resistivity section")
    station_api.apply(
        axes[1],
        x,
        labels,
        preset="inversion",
        xlim=(-0.5, len(labels) - 0.5),
    )
    cb1 = fig.colorbar(im1, ax=axes[1], pad=0.015, shrink=0.82)
    cb1.set_label("log10 rho")
    for ax in axes:
        _style_axis(ax)
    _save(fig, "ai_2d_pseudosection_model_grid.png")


def _ai_3d_grid_plot(functions, sites) -> None:
    station_api = functions["PYCSAMT_STATION_RENDERING"]
    labels = [
        str(getattr(site, "station", getattr(site, "id", i)))
        for i, site in enumerate(sites)
    ]
    n = len(labels)
    x = np.linspace(0, 13_500, n)
    y = 450 * np.sin(np.linspace(0, 2.6, n))
    response = 0.45 + 0.40 * np.exp(-((x - 7200) ** 2 / 1.2e7 + y**2 / 8e5))
    topo_x = np.arange(n, dtype=float)
    topo_m = 420.0 + 95.0 * np.sin(np.linspace(0.15, 2.9, n))
    topo_m += 42.0 * np.exp(-((topo_x - 18.0) ** 2) / 42.0)
    topo_m -= 26.0 * np.exp(-((topo_x - 6.0) ** 2) / 18.0)

    depths = np.linspace(0, 2600, 90)
    xx, zz = np.meshgrid(np.arange(n), depths)
    slice_y0 = 2.24 + 0.52 * np.tanh((zz - 800) / 520)
    slice_y0 -= 0.62 * np.exp(
        -((xx - 12.0) ** 2 / 38 + (zz - 620) ** 2 / 130000)
    )
    slice_y0 += 0.34 * np.exp(
        -((xx - 21.0) ** 2 / 28 + (zz - 1550) ** 2 / 210000)
    )

    fig = plt.figure(figsize=(13.4, 6.4), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.35])
    ax0 = fig.add_subplot(gs[0, 0])
    right = gs[0, 1].subgridspec(2, 1, height_ratios=[0.38, 1.0], hspace=0.04)
    ax_topo = fig.add_subplot(right[0, 0])
    ax1 = fig.add_subplot(right[1, 0], sharex=ax_topo)

    for i in range(n):
        dist = np.hypot(x - x[i], y - y[i])
        for j in np.where((dist > 0) & (dist < 1800))[0]:
            if j > i:
                ax0.plot(
                    [x[i], x[j]],
                    [y[i], y[j]],
                    color="#94a3b8",
                    alpha=0.45,
                    linewidth=0.75,
                )
    sc = ax0.scatter(
        x,
        y,
        c=response,
        cmap="viridis",
        s=62,
        edgecolor="#27323a",
        linewidth=0.45,
        zorder=3,
    )
    step = max(1, n // 10)
    for i in range(0, n, step):
        ax0.text(
            x[i], y[i] + 180, labels[i], fontsize=7, rotation=35, ha="center"
        )
    ax0.set_title("3-D AI station graph")
    ax0.set_xlabel("Easting offset (m)")
    ax0.set_ylabel("Northing offset (m)")
    cb0 = fig.colorbar(sc, ax=ax0, pad=0.015, shrink=0.82)
    cb0.set_label("response feature")
    _style_axis(ax0)

    ax_topo.fill_between(
        topo_x, topo_m.min() - 70.0, topo_m, color="#d7c8a6", alpha=0.70
    )
    ax_topo.plot(topo_x, topo_m, color="#78664a", linewidth=1.7)
    ax_topo.set_ylim(topo_m.min() - 70.0, topo_m.max() + 150.0)
    ax_topo.set_ylabel("Elev. (m)")
    ax_topo.set_title("Topography and station positions along the 3-D slice")
    station_api.style_for("inversion").apply(
        ax_topo,
        topo_x,
        labels,
        xlim=(-0.5, n - 0.5),
        topo_elev=topo_m,
    )
    _style_axis(ax_topo)
    ax_topo.tick_params(
        axis="x", which="both", bottom=False, labelbottom=False
    )

    im = ax1.imshow(
        slice_y0,
        aspect="auto",
        origin="upper",
        extent=(-0.5, n - 0.5, depths.max() / 1000.0, 0.0),
        cmap="turbo",
        vmin=1.35,
        vmax=3.1,
    )
    station_api.apply(
        ax1,
        topo_x,
        labels,
        preset="inversion",
        xlim=(-0.5, n - 0.5),
    )
    ax1.set_ylabel("Depth (km)")
    ax1.set_title("3-D AI vertical slice through profile")
    cb1 = fig.colorbar(im, ax=ax1, pad=0.015, shrink=0.82)
    cb1.set_label("log10 rho")
    _style_axis(ax1)
    _save(fig, "ai_3d_graph_slice_grid.png")


def _ai_3d_topography_block_plot(functions, sites) -> None:
    station_api = functions["PYCSAMT_STATION_RENDERING"]
    labels = [
        str(getattr(site, "station", getattr(site, "id", i)))
        for i, site in enumerate(sites)
    ]
    n = len(labels)
    topo_x = np.arange(n, dtype=float)
    topo_m = 420.0 + 95.0 * np.sin(np.linspace(0.15, 2.9, n))
    topo_m += 42.0 * np.exp(-((topo_x - 18.0) ** 2) / 42.0)
    topo_m -= 26.0 * np.exp(-((topo_x - 6.0) ** 2) / 18.0)

    depths = np.linspace(0, 2600, 90)
    xx, zz = np.meshgrid(np.arange(n), depths)
    block = 2.24 + 0.52 * np.tanh((zz - 800) / 520)
    block -= 0.62 * np.exp(-((xx - 12.0) ** 2 / 38 + (zz - 620) ** 2 / 130000))
    block += 0.34 * np.exp(
        -((xx - 21.0) ** 2 / 28 + (zz - 1550) ** 2 / 210000)
    )
    x_grid = np.tile(topo_x[None, :], (depths.size, 1))
    elev_grid = topo_m[None, :] - depths[:, None]

    fig, ax1 = plt.subplots(figsize=(12.8, 5.8), constrained_layout=True)

    levels = np.linspace(1.35, 3.1, 18)
    im = ax1.contourf(
        x_grid, elev_grid, block, levels=levels, cmap="turbo", extend="both"
    )
    ax1.contour(
        x_grid,
        elev_grid,
        block,
        levels=levels[::2],
        colors="k",
        linewidths=0.28,
        alpha=0.30,
    )
    ax1.fill_between(
        topo_x, topo_m, topo_m.max() + 420.0, color="white", zorder=4
    )
    ax1.plot(topo_x, topo_m, color="#2f2419", linewidth=1.5, zorder=6)
    station_api.style_for("inversion").apply(
        ax1,
        topo_x,
        labels,
        xlim=(-0.5, n - 0.5),
        topo_elev=topo_m,
    )
    ax1.set_ylim(topo_m.min() - depths.max(), topo_m.max() + 420.0)
    ax1.set_ylabel("Elevation (m)")
    ax1.set_title("3-D AI inversion block with topography", pad=16)
    cb1 = fig.colorbar(im, ax=ax1, pad=0.015, shrink=0.82)
    cb1.set_label("log10 rho")
    _style_axis(ax1)
    _save(fig, "ai_3d_topography_inversion_block.png")


def _validation_plot() -> None:
    stations = [f"18-{i:03d}" for i in range(1, 29)]
    x = np.arange(len(stations))
    rng = np.random.default_rng(7)
    rms_1d = 0.28 + 0.08 * np.sin(x / 3) + rng.normal(0, 0.025, len(x))
    rms_2d = 0.22 + 0.04 * np.sin(x / 4 + 0.5) + rng.normal(0, 0.018, len(x))
    uncert = (
        0.08
        + 0.06 * np.exp(-((x - 17) ** 2) / 45)
        + rng.normal(0, 0.008, len(x))
    )
    fig, ax = plt.subplots(figsize=(11.0, 4.2))
    ax.plot(
        x, rms_1d, color="#7c4d79", marker="o", markersize=3, label="1-D RMS"
    )
    ax.plot(
        x, rms_2d, color="#2f6f8f", marker="s", markersize=3, label="2-D RMS"
    )
    ax.fill_between(
        x,
        rms_2d - uncert,
        rms_2d + uncert,
        color="#2f6f8f",
        alpha=0.16,
        label="prediction interval",
    )
    ax.axhline(
        0.35,
        color="#c85745",
        linestyle="--",
        linewidth=1.0,
        label="review threshold",
    )
    ax.set_xticks(x[::3])
    ax.set_xticklabels(stations[::3], rotation=45, ha="right")
    ax.set_ylabel("RMS in log10 apparent resistivity")
    ax.set_title("Post-prediction validation by station")
    ax.legend(fontsize=8, ncol=4, loc="upper right")
    _style_axis(ax)
    _save(fig, "ai_validation_residuals.png")


def _graph_context_plot() -> None:
    rng = np.random.default_rng(11)
    x = np.linspace(0, 13_500, 28) + rng.normal(0, 120, 28)
    y = 500 * np.sin(np.linspace(0, 2.4, 28)) + rng.normal(0, 80, 28)
    response = 0.45 + 0.35 * np.exp(-((x - 7200) ** 2 / 9e6 + y**2 / 6e5))
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    for i in range(len(x)):
        dist = np.hypot(x - x[i], y - y[i])
        for j in np.where((dist > 0) & (dist < 1800))[0]:
            if j > i:
                ax.plot(
                    [x[i], x[j]],
                    [y[i], y[j]],
                    color="#94a3b8",
                    alpha=0.45,
                    linewidth=0.8,
                )
    sc = ax.scatter(
        x,
        y,
        c=response,
        cmap="viridis",
        s=70,
        edgecolor="#27323a",
        linewidth=0.5,
        zorder=3,
    )
    ax.set_xlabel("Profile easting (m)")
    ax.set_ylabel("Profile northing (m)")
    ax.set_title("3-D/graph AI uses station geometry")
    cb = fig.colorbar(sc, ax=ax, pad=0.015)
    cb.set_label("Normalised response feature")
    _style_axis(ax)
    _save(fig, "ai_3d_graph_context.png")


def main() -> int:
    functions = _import_pycsamt()
    survey = functions["read_edis"](
        DATA_DIR,
        recursive=False,
        strict=False,
        on_dup="replace",
        progress=False,
    )
    sites = survey.collection
    inventory = survey.summary().to_pandas(copy=True)
    qc = functions["build_qc_table"](
        sites,
        include_skew=True,
        recursive=False,
        api=True,
    ).to_pandas(copy=True)
    confidence = functions["station_confidence_table"](
        sites,
        method="composite",
        recursive=False,
        api=True,
    ).to_pandas(copy=True)

    _workflow_plot()
    _mode_decision_plot()
    _training_coverage_plot(qc)
    _predicted_sections()
    _ai_2d_grid_plot(functions, sites)
    _validation_plot()
    _graph_context_plot()
    _ai_3d_grid_plot(functions, sites)
    _ai_3d_topography_block_plot(functions, sites)

    print("survey_summary:")
    print(survey.summary())
    print("inventory:")
    print(
        inventory[["station", "n_freq", "tipper", "spectra"]]
        .head(5)
        .to_string(index=False)
    )
    print("qc_ready:")
    print(
        qc[["station", "n_freq", "frac_ok", "snr_med", "pmin", "pmax"]]
        .head(5)
        .to_string(index=False, float_format=lambda value: f"{value:.3g}")
    )
    print("confidence:")
    print(
        confidence[["station", "confidence", "coverage", "offdiag", "phase"]]
        .head(5)
        .to_string(index=False, float_format=lambda value: f"{value:.3f}")
    )
    print("decision:")
    print("corrected_edis=data/AMT/WILLY_DATA/L18PLT")
    print("1d_ai=optional_screening")
    print("2d_ai=recommended_for_profile")
    print("3d_ai=use_when_coordinates_or_multiple_lines_exist")
    print(f"images: {IMAGE_DIR.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
