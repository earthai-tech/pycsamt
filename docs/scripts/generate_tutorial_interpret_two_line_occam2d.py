"""Interpret the two real Occam2D lines built by
build_two_line_occam2d_survey.rst / generate_tutorial_two_line_occam2d.py:
calibrate against boreholes, classify into lithology, record structural
evidence, and validate against held-out boreholes.

As in the companion script, every reusable step is real, tested pycsamt
library code (:class:`pycsamt.geology.Borehole`,
:class:`pycsamt.interp.ModelCalibrator`,
:class:`pycsamt.geology.FaultTrace`); only figure rendering and this
tutorial's specific borehole positions live here.
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
sys.path.insert(0, str(ROOT / "docs" / "scripts"))

from generate_tutorial_two_line_occam2d import (  # noqa: E402
    LINE_A,
    LINE_B,
    LINE_TOPO_KWARGS,
    RUN_ROOT,
    build_true_grid,
)
from pycsamt.geology import (  # noqa: E402
    Borehole,
    FaultTrace,
    Interval,
    RockDatabase,
    RockEntry,
    StructuralModel,
)
from pycsamt.interp import ModelCalibrator, ResistivityModel  # noqa: E402
from pycsamt.interp import plot as iplot  # noqa: E402
from pycsamt.models.occam2d import InversionResult  # noqa: E402
from pycsamt.topo import (  # noqa: E402
    interp_elev,
    plot_topo_array,
    plot_topo_section,
    synthetic_elevation_profile,
)

IMAGES = ROOT / "docs/source/images/tutorials/interpret_two_line_occam2d_survey"

# A project-specific table, not the generic built-in one -- see
# geology/rock_database.rst for why classifying against the wrong table
# changes the answer even though the resistivity does not.
REGIONAL_DB = RockDatabase([
    RockEntry(name="Overburden", rho_min=30, rho_max=180, color="#D4AC0D",
              source="Site-specific, this survey"),
    RockEntry(name="Weathered basement", rho_min=180, rho_max=900, color="#A9780C",
              source="Site-specific, this survey"),
    RockEntry(name="Fresh basement", rho_min=900, rho_max=8000, color="#4A4A4A",
              source="Site-specific, this survey"),
])

# (name, x_m, max_depth_m) -- positions chosen clear of the fault's own
# lateral migration with depth (see the tutorial text for why that zone
# is excluded rather than merely avoided by luck).
LINE_A_CALIB_BOREHOLES = [("A1", 300.0, 450.0), ("A2", 2100.0, 1100.0)]
LINE_A_VALID_BOREHOLES = [("A3", 900.0, 450.0), ("A4", 1800.0, 1100.0)]
LINE_B_CALIB_BOREHOLES = [("B1", 300.0, 1100.0), ("B2", 1900.0, 450.0)]
LINE_B_VALID_BOREHOLES = [("B3", 100.0, 1100.0), ("B4", 2100.0, 450.0)]


def true_intervals(grid, x_m: float, max_depth_m: float, dz: float = 1.0):
    """Real interval boundaries at along-profile *x_m*, read directly off
    the true grid with :meth:`~pycsamt.forward.grid2d.Grid2D.value_at` --
    the same technique borehole ground truth uses throughout this page,
    guaranteed consistent with the grid actually solved and inverted."""
    zs = np.arange(0.0, max_depth_m + dz, dz)
    vals = np.array([grid.value_at(x_m, z, chainage=True) for z in zs])
    intervals = []
    start = 0
    for i in range(1, len(zs)):
        if vals[i] != vals[i - 1]:
            intervals.append((float(zs[start]), float(zs[i]), float(vals[i - 1])))
            start = i
    intervals.append((float(zs[start]), float(zs[-1]), float(vals[-1])))
    return intervals


def borehole_from_true_grid(grid, name: str, x_m: float, max_depth_m: float) -> Borehole:
    """Build a real :class:`~pycsamt.geology.Borehole` from the true grid's
    own resistivity -- what a real drilling program would have measured
    at this position, classified with the same regional table used
    throughout this page."""
    intervals = [
        Interval(top=top, bottom=bottom, lithology=REGIONAL_DB.classify(rho).name,
                  resistivity=rho)
        for top, bottom, rho in true_intervals(grid, x_m, max_depth_m)
    ]
    return Borehole(name, x=x_m, intervals=intervals)


def build_boreholes(line: dict, calib_spec, valid_spec):
    grid = build_true_grid(line)
    calib = [borehole_from_true_grid(grid, *spec) for spec in calib_spec]
    valid = [borehole_from_true_grid(grid, *spec) for spec in valid_spec]
    return calib, valid


def load_calibrated_model(workdir_name: str, calib_boreholes):
    result = InversionResult(RUN_ROOT / workdir_name)
    model = ResistivityModel.from_occam2d(result).clip_to_stations()
    cal = ModelCalibrator(ptol=0.10, db=REGIONAL_DB, verbose=False).fit(model, calib_boreholes)
    return model, cal


def fault_trace_for_line(line: dict) -> FaultTrace:
    return FaultTrace(
        x=line["fault_x_m"], dip_deg=line["apparent_dip_deg"],
        downthrown_side=line["downthrown_side"], throw_m=line["throw_m"],
        sense="normal", evidence="true model (synthetic control)",
    )


def elevation_for_line(line: dict, station_x: np.ndarray) -> np.ndarray:
    """Real terrain elevation at each station, real Line A/B config --
    the same :func:`~pycsamt.topo.synthetic_elevation_profile` /
    ``LINE_TOPO_KWARGS`` combination the topography figure in Part 1
    uses, kept in one place so every figure on this page drapes exactly
    the same terrain."""
    return synthetic_elevation_profile(station_x, **LINE_TOPO_KWARGS[line["label"]])


def validation_table(model: ResistivityModel, cal: ModelCalibrator, valid_boreholes):
    """Compare each validation borehole's true top-layer boundary against
    the classified station log nearest to it -- station spacing here (100 m)
    means every borehole above lands exactly on a station, so "nearest" is
    exact, not approximate."""
    logs = {log.station_name: log for log in cal.stratigraphic_logs()}
    rows = []
    for bh in valid_boreholes:
        ix = int(np.argmin(np.abs(model.station_x - bh.x)))
        name = model.station_names[ix]
        log = logs[name]
        rows.append({
            "borehole": bh.name,
            "x_m": bh.x,
            "station": name,
            "true_overburden_base_m": bh.intervals[0].bottom,
            "classified_overburden_base_m": log.layers[0].bottom if log.layers else float("nan"),
            "true_basement_top_m": bh.intervals[-1].top,
            "classified_basement_top_m": (
                next((l.top for l in log.layers if l.lithology == "Fresh basement"), float("nan"))
            ),
        })
    return rows


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def make_fence_diagrams():
    """Classified fence diagram for every station, with a real terrain
    strip and the shared ``inversion`` station marker drawn above the
    panels -- the same station-index convention
    pycsamt.topo.draw_topo_strip uses for pseudosections, applied here
    via PlotFenceDiagram's own elevation_m parameter."""
    for line, calib_spec, valid_spec, workdir in (
        (LINE_A, LINE_A_CALIB_BOREHOLES, LINE_A_VALID_BOREHOLES, "line_a_occam2d"),
        (LINE_B, LINE_B_CALIB_BOREHOLES, LINE_B_VALID_BOREHOLES, "line_b_occam2d"),
    ):
        calib, _ = build_boreholes(line, calib_spec, valid_spec)
        model, cal = load_calibrated_model(workdir, calib)
        logs = cal.stratigraphic_logs()
        elev = elevation_for_line(line, model.station_x)
        fig = iplot.PlotFenceDiagram(
            logs, max_depth=1600.0, elevation_m=elev,
            title=f"{line['label']} -- classified fence diagram",
        ).plot()
        fig.savefig(IMAGES / f"fence_{line['label'].split()[-1].lower()}.png", dpi=170, bbox_inches="tight")
        plt.close(fig)


def make_borehole_fence_figures():
    """Compare every field borehole on a line side by side, ordered by
    profile position -- the raw ground truth the calibrator and the
    validation step both draw on, before any inversion model enters
    the picture."""
    for line, calib_spec, valid_spec in (
        (LINE_A, LINE_A_CALIB_BOREHOLES, LINE_A_VALID_BOREHOLES),
        (LINE_B, LINE_B_CALIB_BOREHOLES, LINE_B_VALID_BOREHOLES),
    ):
        calib, valid = build_boreholes(line, calib_spec, valid_spec)
        boreholes = sorted(calib + valid, key=lambda bh: bh.x)
        fig = iplot.PlotBoreholeFence(
            boreholes, db=REGIONAL_DB,
            title=f"{line['label']} -- field boreholes (ground truth)",
        ).plot()
        letter = line["label"].split()[-1].lower()
        fig.savefig(IMAGES / f"borehole_fence_{letter}.png", dpi=170, bbox_inches="tight")
        plt.close(fig)


def make_calibration_effect_figures():
    """CRM vs calibrated NM vs the G (%) misfit map, draped over real
    topography with the shared inversion station marker, for both
    lines -- the direct visual record of where ModelCalibrator's
    soft-replace step actually changed the raw inversion, and by how
    much.

    All three panels are clipped to the top 1600 m -- Occam2D's own
    mesh reaches almost 10 km depth for numerical-boundary reasons,
    which would otherwise squeeze the entire zone of geological
    interest into a sliver of the figure. CRM/NM use
    pycsamt.topo.plot_topo_section directly (they are real
    ResistivityModel objects); the misfit map is not a resistivity
    grid, so it uses the generic pycsamt.topo.plot_topo_array instead.
    misfit_map()'s per-column G (%) runs 6-13% almost everywhere on
    this survey (autolayer reclassification touches every column
    outside max_borehole_distance of a borehole, not just the
    calibrated ones), so vmax=15 is used instead of the 0-10% range
    that would suit a typical, mostly-untouched line.

    CRM and NM share one colour scale, so they share one colorbar
    rather than two identical ones. The three panels are laid out on a
    GridSpec with a dedicated, fixed-width colorbar column instead of
    letting each plotting call attach its own -- plot_topo_section's
    colorbar (mpl_toolkits' make_axes_locatable, a fixed percentage of
    its own axes) and a bare fig.colorbar (matplotlib's default
    fraction/pad) do not reserve the same width, which otherwise leaves
    the misfit panel a different width from CRM/NM.
    """
    for line, calib_spec, valid_spec, workdir in (
        (LINE_A, LINE_A_CALIB_BOREHOLES, LINE_A_VALID_BOREHOLES, "line_a_occam2d"),
        (LINE_B, LINE_B_CALIB_BOREHOLES, LINE_B_VALID_BOREHOLES, "line_b_occam2d"),
    ):
        calib, _ = build_boreholes(line, calib_spec, valid_spec)
        model, cal = load_calibrated_model(workdir, calib)
        nm = cal.calibrated_model()
        mm = cal.misfit_map()
        elev = elevation_for_line(line, model.station_x)

        fig = plt.figure(figsize=(12, 12))
        gs = fig.add_gridspec(
            3, 2, width_ratios=[30, 1], height_ratios=[1, 1, 1.08],
            hspace=0.65, wspace=0.06,
        )
        ax_crm = fig.add_subplot(gs[0, 0])
        ax_nm = fig.add_subplot(gs[1, 0], sharex=ax_crm)
        ax_g = fig.add_subplot(gs[2, 0], sharex=ax_crm)
        cax_rho = fig.add_subplot(gs[0:2, 1])
        cax_g = fig.add_subplot(gs[2, 1])

        common = dict(
            elevation=elev, station_x=model.station_x,
            station_names=model.station_names, depth_max=1600.0,
        )
        plot_topo_section(
            model, ax=ax_crm, cmap="jet", vmin=1.0, vmax=5.0,
            colorbar=False, title="CRM -- inversion result", **common,
        )
        plot_topo_section(
            nm, ax=ax_nm, cmap="jet", vmin=1.0, vmax=5.0,
            colorbar=False, title="NM -- calibrated model", **common,
        )
        plot_topo_array(
            model.x_centers, model.z_centers, mm, ax=ax_g,
            cmap="RdYlBu_r", vmin=0.0, vmax=15.0, colorbar=False,
            title="Misfit G (%)", **common,
        )
        fig.colorbar(
            ax_crm.collections[0], cax=cax_rho,
            label=r"$\log_{10}\rho$ ($\Omega\cdot$m)",
        )
        fig.colorbar(ax_g.collections[0], cax=cax_g, label="G (%)")

        # Only the bottom panel needs the shared x-axis label/ticks --
        # repeating "Profile distance (km)" on all three, this close
        # together, just crowds the titles above each lower panel.
        ax_crm.set_xlabel("")
        ax_nm.set_xlabel("")
        ax_crm.tick_params(labelbottom=False)
        ax_nm.tick_params(labelbottom=False)

        fig.suptitle(
            f"{line['label']} -- CRM vs calibrated NM, real topography (top 1600 m)",
            fontweight="bold", y=0.995,
        )
        letter = line["label"].split()[-1].lower()
        fig.savefig(IMAGES / f"calibration_effect_{letter}.png", dpi=170, bbox_inches="tight")
        plt.close(fig)


def make_structural_section_figure():
    """Draw both lines' calibrated sections draped over real topography,
    with the true fault trace and every borehole overlaid in the same
    terrain-following (elevation, not flat depth) coordinates the
    background section itself uses.

    pycsamt.topo.plot_topo_section does the draping and draws the
    station markers/labels; the fault trace and borehole traces are
    plotted afterwards on the same axes, converted from
    (x, depth)-space into (x, elevation)-space with the same
    pycsamt.topo.interp_elev this page's boreholes are built against.
    """
    fig, axes = plt.subplots(2, 1, figsize=(11, 10))
    specs = (
        (LINE_A, LINE_A_CALIB_BOREHOLES, LINE_A_VALID_BOREHOLES, "line_a_occam2d"),
        (LINE_B, LINE_B_CALIB_BOREHOLES, LINE_B_VALID_BOREHOLES, "line_b_occam2d"),
    )
    for ax, (line, calib_spec, valid_spec, workdir) in zip(axes, specs):
        calib, valid = build_boreholes(line, calib_spec, valid_spec)
        model, cal = load_calibrated_model(workdir, calib)
        nm = cal.calibrated_model()
        elev = elevation_for_line(line, model.station_x)
        chain_km = model.station_x / 1000.0
        elev_km = elev / 1000.0

        plot_topo_section(
            nm, ax=ax, elevation=elev, station_x=model.station_x,
            station_names=model.station_names, depth_max=1600.0,
            cmap="turbo_r", vmin=1.5, vmax=4.0,
            title=f"{line['label']} -- calibrated model, true fault, and boreholes",
        )

        def elev_at_km(x_m):
            return float(interp_elev(chain_km, elev_km, np.array([x_m / 1000.0]))[0])

        fault = fault_trace_for_line(line)
        z = np.linspace(0.0, 1600.0, 60)
        direction = 1.0 if fault.downthrown_side == "right" else -1.0
        x_line_m = fault.x + direction * z / np.tan(np.deg2rad(fault.dip_deg))
        y_line_km = np.array(
            [elev_at_km(xi) - zi / 1000.0 for xi, zi in zip(x_line_m, z)]
        )
        ax.plot(x_line_m / 1000.0, y_line_km, color="white", linewidth=2.5,
                linestyle="--", zorder=6,
                label=f"true fault ({fault.dip_deg:.0f} deg apparent dip)")

        for bh, bold, ls, lw in (
            *((bh, True, "-", 2.0) for bh in calib),
            *((bh, False, ":", 1.2) for bh in valid),
        ):
            e0_km = elev_at_km(bh.x)
            y_top_km = e0_km
            y_bot_km = e0_km - bh.max_depth / 1000.0
            ax.plot([bh.x / 1000.0, bh.x / 1000.0], [y_top_km, y_bot_km],
                    color="black", linewidth=lw, linestyle=ls, zorder=6)
            ax.text(bh.x / 1000.0, y_bot_km - 0.03, bh.name, ha="center",
                    va="top", fontsize=8, fontweight="bold" if bold else None,
                    style=None if bold else "italic", zorder=7)

        ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(IMAGES / "structural_overlay.png", dpi=170, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    IMAGES.mkdir(parents=True, exist_ok=True)
    make_borehole_fence_figures()
    make_calibration_effect_figures()
    make_fence_diagrams()
    make_structural_section_figure()
    print("figures written")
