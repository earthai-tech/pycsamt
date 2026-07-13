"""Generate figures and sample outputs for the Occam2D preparation tutorial."""

from __future__ import annotations

import contextlib
import io
import shutil
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
IMAGE_DIR = (
    ROOT
    / "docs"
    / "source"
    / "images"
    / "tutorials"
    / "prepare_occam2d_inversion"
)
BUILD_DIR = ROOT / "docs" / "_generated" / "tutorial_occam2d"


def _import_pycsamt():
    sys.path.insert(0, str(ROOT))
    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        from pycsamt.api import read_edis
        from pycsamt.emtools.qc import station_confidence_table
        from pycsamt.models.occam2d import InputBuilder, OccamConfig
        from pycsamt.models.occam2d.validation import detect_file_type

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


def _data_coverage_plot(builder) -> None:
    data = builder.data
    blocks = data.data_blocks
    sites = list(data.sites)
    codes = blocks[:, 2].astype(int)
    mode_masks = {
        "TE rho": codes == 1,
        "TE phase": codes == 2,
        "TM rho": codes == 5,
        "TM phase": codes == 6,
    }
    x = np.arange(len(sites))
    width = 0.2
    fig, ax = plt.subplots(figsize=(11.3, 4.2))
    colors = ["#2f6f8f", "#82a6b1", "#c85745", "#e2a36c"]
    for i, (label, mask) in enumerate(mode_masks.items()):
        counts = np.zeros(len(sites), dtype=int)
        for site_idx in blocks[mask, 0].astype(int):
            counts[site_idx - 1] += 1
        ax.bar(
            x + (i - 1.5) * width,
            counts,
            width=width,
            label=label,
            color=colors[i],
            edgecolor="#27323a",
            linewidth=0.35,
            alpha=0.88,
        )
    ax.set_xticks(x[::2])
    ax.set_xticklabels([s.replace("23-", "") for s in sites[::2]], rotation=45, ha="right")
    ax.set_ylabel("Rows written per station")
    ax.set_xlabel("Station")
    ax.set_title("OccamDataFile rows by station and data type")
    ax.legend(fontsize=8, ncol=4, loc="upper right")
    _style_axis(ax)
    _save(fig, "occam_data_rows_by_station.png")


def _frequency_mesh_plot(builder) -> None:
    data = builder.data
    mesh = builder.mesh
    fig, axes = plt.subplots(1, 2, figsize=(11.3, 4.4), constrained_layout=True)

    freqs = np.asarray(data.frequencies, dtype=float)
    periods = 1.0 / freqs
    axes[0].plot(
        np.arange(1, len(freqs) + 1),
        periods,
        color="#2f6f8f",
        marker="o",
        markersize=4,
        linewidth=1.4,
    )
    axes[0].set_yscale("log")
    axes[0].set_xlabel("Occam frequency index")
    axes[0].set_ylabel("Period (s)")
    axes[0].set_title("Selected frequency band")
    _style_axis(axes[0])

    x_nodes = np.asarray(mesh.x_nodes, dtype=float)
    z_nodes = np.asarray(mesh.z_nodes, dtype=float)
    for x in x_nodes:
        axes[1].plot([x, x], [z_nodes[0], z_nodes[-1]], color="#7e8b93", linewidth=0.35)
    for z in z_nodes:
        axes[1].plot([x_nodes[0], x_nodes[-1]], [z, z], color="#7e8b93", linewidth=0.35)
    axes[1].invert_yaxis()
    axes[1].set_xlabel("Profile mesh coordinate (m)")
    axes[1].set_ylabel("Depth mesh coordinate (m)")
    axes[1].set_title(f"Mesh skeleton: {mesh.n_xcells} x {mesh.n_zcells} cells")
    axes[1].grid(False)
    axes[1].set_facecolor("#fbfbf7")
    for spine in axes[1].spines.values():
        spine.set_color("#39434d")
    _save(fig, "occam_frequency_mesh_summary.png")


def _file_size_plot(workdir: Path, cfg) -> None:
    files = [cfg.data_file, cfg.mesh_file, cfg.model_file, cfg.startup_file]
    sizes = [((workdir / name).stat().st_size / 1024.0) for name in files]
    fig, ax = plt.subplots(figsize=(8.6, 3.8))
    ax.bar(
        files,
        sizes,
        color=["#2f6f8f", "#7c4d79", "#d5962c", "#3f8f61"],
        edgecolor="#27323a",
        linewidth=0.5,
    )
    ax.set_ylabel("File size (KiB)")
    ax.set_title("Native Occam2D files written by InputBuilder")
    ax.tick_params(axis="x", rotation=20)
    _style_axis(ax)
    _save(fig, "occam_native_file_sizes.png")


def main() -> int:
    functions = _import_pycsamt()
    if BUILD_DIR.exists():
        shutil.rmtree(BUILD_DIR)
    workdir = BUILD_DIR / "L18PLT_occam2d_native"
    tm_workdir = BUILD_DIR / "L18PLT_occam2d_tm"
    workdir.mkdir(parents=True, exist_ok=True)

    survey = functions["read_edis"](
        DATA_DIR,
        recursive=False,
        strict=False,
        progress=False,
    )
    sites = survey.collection
    confidence = functions["station_confidence_table"](
        sites,
        method="composite",
        api=True,
    ).to_pandas(copy=True)

    cfg = functions["OccamConfig"](
        modes=["TE", "TM"],
        freq_min=0.1,
        freq_max=1000.0,
        error_floor_rho=0.05,
        error_floor_phase=0.5,
        n_layers=32,
        n_airlayers=4,
        cell_size_horizontal=100.0,
        cell_size_vertical_top=25.0,
        depth_scale=1.15,
        target_misfit=1.0,
        max_iterations=80,
        initial_rho=100.0,
    )
    cfg.to_template(workdir / "occam2d.yml")
    builder = functions["InputBuilder"](
        sites,
        workdir=workdir,
        config=cfg,
        verbose=0,
    ).build(title="L18PLT pyCSAMT v2 Occam2D preparation")

    tm_builder = functions["InputBuilder"](
        sites,
        workdir=tm_workdir,
        config=functions["OccamConfig"](),
        verbose=0,
    ).build(
        modes=["TM"],
        freq_min=0.2,
        freq_max=500.0,
        error_floor_rho=0.07,
        error_floor_phase=1.0,
        n_layers=28,
        cell_size=150.0,
        title="L18PLT TM-only Occam2D first trial",
    )

    _data_coverage_plot(builder)
    _frequency_mesh_plot(builder)
    _file_size_plot(workdir, cfg)

    print("confidence_head:")
    print(confidence[["station", "confidence", "coverage"]].head(5).to_string(index=False))
    print("builder_summary:")
    print(builder.summary())
    print("object_counts:")
    print(
        f"sites={builder.data.n_sites} frequencies={builder.data.n_frequencies} "
        f"data_rows={builder.data.n_data} mesh={builder.mesh.n_xcells}x{builder.mesh.n_zcells} "
        f"parameters={builder.model.n_params} startup={cfg.startup_file}"
    )
    print("file_types:")
    for filename in (cfg.data_file, cfg.mesh_file, cfg.model_file, cfg.startup_file):
        path = workdir / filename
        kind = functions["detect_file_type"](path)
        print(f"{path.name:18s} -> {kind:7s} {path.stat().st_size:7d} bytes")
    print("expected_files:")
    for path in sorted(workdir.iterdir()):
        print(path.name)
    print("tm_summary:")
    print(tm_builder.summary())
    print(f"images: {IMAGE_DIR.relative_to(ROOT)}")
    shutil.rmtree(BUILD_DIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
