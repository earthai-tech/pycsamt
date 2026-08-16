"""Build the two synthetic Occam2D lines used by the two-line interpretation
tutorial pair (build_two_line_occam2d_survey.rst /
interpret_two_line_occam2d_survey.rst), and render their figures.

This script is documentation-only glue: every real step -- the true earth
model, the forward physics, the site packaging, the display topography --
is real, tested pycsamt library code
(:meth:`pycsamt.forward.grid2d.Grid2D.layered_with_fault`,
:func:`pycsamt.models.occam2d.sites_from_response`,
:func:`pycsamt.topo.synthetic_elevation_profile`,
:meth:`pycsamt.interp.ResistivityModel.clip_to_stations`), not duplicated
here. Only figure rendering lives in this file.
"""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.forward.em2d import MT2DForward  # noqa: E402
from pycsamt.forward.grid2d import Grid2D  # noqa: E402
from pycsamt.forward.synthetic import LayeredModel  # noqa: E402
from pycsamt.interp import ResistivityModel  # noqa: E402
from pycsamt.models.occam2d import (  # noqa: E402
    InputBuilder,
    InversionResult,
    OccamConfig,
    OccamRunner,
    sites_from_response,
)
from pycsamt.topo import plot_topo_section, synthetic_elevation_profile  # noqa: E402

IMAGES = ROOT / "docs/source/images/tutorials/build_two_line_occam2d_survey"
RUN_ROOT = ROOT / "runs"
OCCAM_BINARY = (ROOT / "pycsamt/models/occam2d/_source/Occam2D.exe").resolve()

FREQS_HZ = np.array(
    [2000.0, 1000.0, 500.0, 200.0, 100.0, 50.0, 20.0, 10.0, 5.0, 2.0, 1.0, 0.5, 0.25]
)

# One shared true layer stack: overburden / weathered zone / basement.
TRUE_LAYERS = LayeredModel(resistivity=[80.0, 300.0, 3000.0], thickness=[60.0, 260.0])

LINE_A = dict(
    label="Line A",
    fault_x_m=1200.0,
    apparent_dip_deg=65.0,
    downthrown_side="right",
    throw_m=700.0,
)
LINE_B = dict(
    label="Line B",
    fault_x_m=1400.0,
    apparent_dip_deg=47.0,
    downthrown_side="left",
    throw_m=700.0,
)

GRID_KWARGS = dict(nx=60, nz=50, x_max=2400.0, z_max=1600.0, n_stations=25)

OCCAM_CFG_KWARGS = dict(
    modes=["TE"],
    n_layers=32,
    n_airlayers=0,
    cell_size_horizontal=60.0,
    cell_size_vertical_top=15.0,
    depth_scale=1.16,
    initial_rho=200.0,
    error_floor_rho=0.10,
    error_floor_phase=2.0,
    target_misfit=1.5,
    max_iterations=80,
)

LINE_TOPO_KWARGS = {
    "Line A": dict(base_m=120.0, amplitude_m=35.0, phase_m=0.0),
    "Line B": dict(base_m=95.0, amplitude_m=28.0, phase_m=650.0),
}


def build_true_grid(line: dict) -> Grid2D:
    """Build one line's true earth model on its own solver grid."""
    kwargs = {k: v for k, v in line.items() if k != "label"}
    return Grid2D.layered_with_fault(TRUE_LAYERS, name=line["label"], **kwargs, **GRID_KWARGS)


def build_line(line: dict, workdir_name: str):
    """Run the complete real pipeline for one line and return
    ``(grid, resp, workdir, result)`` with *result* a real, loaded
    :class:`~pycsamt.models.occam2d.InversionResult`.
    """
    grid = build_true_grid(line)
    resp = MT2DForward(FREQS_HZ, grid, verbose=False).run()

    names = [f"L{i:02d}" for i in range(grid.n_stations)]
    sites = sites_from_response(resp, grid.x_stations, names)

    workdir = RUN_ROOT / workdir_name
    cfg = OccamConfig(**OCCAM_CFG_KWARGS)
    builder = InputBuilder(sites, workdir=workdir, config=cfg, verbose=0)
    builder.build(title=f"Synthetic {line['label']}")

    runner = OccamRunner(workdir=workdir, binary_path=OCCAM_BINARY, verbose=0)
    runner.run(max_iter=cfg.max_iterations, target_misfit=cfg.target_misfit, auto_compile=False)

    result = InversionResult(workdir)
    return grid, resp, workdir, result


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def make_true_models_figure():
    """Write the two-panel true-model figure (both lines' Grid2D), draped
    over the same real terrain the inversion-result and structural
    figures use, via :func:`pycsamt.topo.plot_topo_section`'s raw
    ``(x_centers, z_centers, rho_2d)`` array-triple input -- Grid2D is
    not a ResistivityModel, so this is the generic entry point rather
    than the ResistivityModel-specific one."""
    fig, axes = plt.subplots(2, 1, figsize=(9, 8))
    for ax, line in zip(axes, (LINE_A, LINE_B)):
        grid = build_true_grid(line)
        names = [f"L{i:02d}" for i in range(grid.n_stations)]
        station_x = grid.x_stations - grid.core_x_offset
        elevation = synthetic_elevation_profile(station_x, **LINE_TOPO_KWARGS[line["label"]])
        plot_topo_section(
            (grid.x_centers - grid.core_x_offset, grid.z_centers, np.log10(grid.resistivity)),
            ax=ax,
            elevation=elevation,
            station_x=station_x,
            station_names=names,
            depth_max=1600.0,
            cmap="turbo_r",
            vmin=1.5,
            vmax=3.6,
            title=f"{line['label']} -- true model",
        )
        # Grid2D.layered_with_fault pads its mesh far past the station-
        # carrying core (matching Occam2D's own padding convention), and
        # plot_topo_section's auto xlim spans that full padded mesh.
        # ResistivityModel.clip_to_stations() crops this for the real
        # inversion figures; Grid2D has no such method, so it is cropped
        # here for display instead, purely visually -- no data is dropped.
        ax.set_xlim(-0.1, 2.5)
    fig.tight_layout()
    fig.savefig(IMAGES / "true_models.png", dpi=170, bbox_inches="tight")
    plt.close(fig)


def make_topography_figure():
    """Write the synthetic display-topography figure for both lines."""
    fig, ax = plt.subplots(figsize=(9, 3.2))
    chainage = np.linspace(0, 2400, 200)
    for line, color in zip((LINE_A, LINE_B), ("#B9770E", "#1F618D")):
        elev = synthetic_elevation_profile(chainage, **LINE_TOPO_KWARGS[line["label"]])
        ax.plot(chainage, elev, color=color, label=line["label"])
    ax.set_xlabel("Profile position x (m)")
    ax.set_ylabel("Elevation (m)")
    ax.set_title("Synthetic display topography (not used by the forward solve or Occam2D)")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(IMAGES / "topography.png", dpi=170, bbox_inches="tight")
    plt.close(fig)


def make_inversion_results_figure(result_a, result_b):
    """Write the two-panel real recovered-resistivity figure, draped over
    the synthetic display topography via :func:`pycsamt.topo.plot_topo_section`."""
    fig, axes = plt.subplots(2, 1, figsize=(9, 8))
    for ax, result, line in zip(axes, (result_a, result_b), (LINE_A, LINE_B)):
        model = ResistivityModel.from_occam2d(result).clip_to_stations()
        elevation = synthetic_elevation_profile(model.station_x, **LINE_TOPO_KWARGS[line["label"]])
        plot_topo_section(
            model,
            ax=ax,
            elevation=elevation,
            chainage=model.station_x / 1000.0,
            topo_source="array",
            depth_max=1600.0,
            cmap="turbo_r",
            vmin=1.5,
            vmax=4.0,
            title=f"{line['label']} -- real Occam2D (TE), RMS={model.rms:.2f}",
        )
    fig.tight_layout()
    fig.savefig(IMAGES / "inversion_results.png", dpi=170, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    IMAGES.mkdir(parents=True, exist_ok=True)
    RUN_ROOT.mkdir(parents=True, exist_ok=True)

    make_true_models_figure()
    make_topography_figure()

    _, _, _, result_a = build_line(LINE_A, "line_a_occam2d")
    _, _, _, result_b = build_line(LINE_B, "line_b_occam2d")
    print("Line A final RMS", result_a.final_rms)
    print("Line B final RMS", result_b.final_rms)

    make_inversion_results_figure(result_a, result_b)
