"""Build and invert the Gabbs Valley Occam1D example, and render its figures.

This script is documentation-only glue for
``docs/source/user_guide/models/occam1d.rst``: the actual inversion runs
through the same public :mod:`pycsamt.models.occam1d` API shown on the page
(:class:`~pycsamt.models.occam1d.Occam1DBatch`,
:class:`~pycsamt.models.occam1d.Occam1DInversion`,
:meth:`~pycsamt.models.occam1d.Occam1DInversion.save_main_images`); only the
multi-station comparison figure is custom, since no dedicated
"compare independent Occam1D stations" plot function exists in the package.
"""

from __future__ import annotations

from pathlib import Path
import shutil
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.models.occam1d import (  # noqa: E402
    Occam1DBatch,
    Occam1DConfig,
    Occam1DInversion,
)
from pycsamt.seg.edi import EDIFile  # noqa: E402

DATA_DIR = ROOT / "data" / "gv_data" / "gv_final_edi"
IMAGES = ROOT / "docs/source/images/user_guide/models"
BUILD_DIR = ROOT / "docs/_generated/tutorial_occam1d"

# gv100/gv130/gv163: the low, middle, and high end of the tracked station
# subset (see data/gv_data/README.md) -- spread across the survey area
# rather than three adjacent soundings. The full gv_final_edi/ directory may
# hold all 59 stations on a machine that downloaded the complete USGS
# release (see data/gv_data/README.md); this script loads the three
# tracked EDIs explicitly rather than the whole directory, both to stay
# reproducible from a fresh checkout and to keep this example's run time
# independent of how many extra local-only EDIs happen to be present.
STATIONS = ("gv100", "gv130", "gv163")

_STATION_COLORS = {"gv100": "#2f6f8f", "gv130": "#c85745", "gv163": "#3f8f61"}


def build_and_invert():
    """Build native Occam1D inputs and run every station's inversion once.

    Returns the ready :class:`Occam1DBatch` and a
    ``{station: (Occam1DInversion, Occam1DInversionResult)}`` mapping so the
    rest of this script (and, separately,
    :meth:`~pycsamt.models.occam1d.Occam1DBatch.invert_all`) never repeats a
    real inversion just to get a different view of the same result.
    """
    config = Occam1DConfig(
        mode="determinant",
        n_layers=30,
        first_thickness=5.0,
        depth_max=5000.0,
        starting_resistivity=100.0,
        target_misfit=1.5,
        max_iterations=25,
    )
    sources = [EDIFile(DATA_DIR / f"{station}.edi") for station in STATIONS]
    batch = Occam1DBatch(
        sources, BUILD_DIR / "gv_occam1d", config=config, verbose=0
    ).build_all()

    runs = {}
    for builder in batch.builders:
        station = builder.data.station
        if station not in STATIONS:
            continue
        inversion = Occam1DInversion(
            builder.data, builder.model, config=builder.config,
            startup=builder.startup, verbose=0,
        )
        result = inversion.run()
        runs[station] = (inversion, result)
    return batch, runs


def make_station_summary_figure(inversion, result) -> None:
    """Write one station's combined model/response/convergence figure.

    Thin wrapper around the public
    :meth:`~pycsamt.models.occam1d.Occam1DInversion.save_main_images`, kept
    here (rather than shown inline) only so the page can present it as one
    call alongside the multi-station figure below.
    """
    paths = inversion.save_main_images(IMAGES, result, dpi=200)
    target = IMAGES / "occam1d_gv100_summary.png"
    shutil.move(str(paths["summary"]), str(target))
    for kind in ("model", "response", "convergence"):
        paths[kind].unlink(missing_ok=True)


def make_valley_resistivity_comparison_figure(runs) -> None:
    """Write the three-station recovered resistivity-vs-depth comparison.

    Each :class:`~pycsamt.models.occam1d.Occam1DInversion` result is
    independent -- Occam1D has no lateral coupling between stations -- so
    this lays the three final layered-earth models out as one row, one
    column per station, sharing a single depth axis. A shared axis (rather
    than three independently scaled panels) is what makes it possible to
    compare the actual depth of a conductive zone across stations at a
    glance, which a per-station summary image cannot show on its own.
    """
    fig, axes = plt.subplots(
        1, len(STATIONS), figsize=(9.6, 6.0), sharey=True,
    )
    for ax, station in zip(axes, STATIONS):
        inversion, result = runs[station]
        depth = inversion.model.depth
        resistivity = 10.0 ** result.final.parameters
        ax.step(
            resistivity, depth, where="post",
            color=_STATION_COLORS[station], linewidth=1.8,
        )
        ax.set_xscale("log")
        ax.set_xlabel(r"Resistivity ($\Omega\,\mathrm{m}$)")
        ax.set_title(station)
        ax.grid(True, which="both", axis="x", alpha=0.3)
    axes[0].set_ylim(5000.0, 0.0)
    axes[0].set_ylabel("Depth (m)")
    fig.suptitle("Gabbs Valley -- independent Occam1D models, 3 stations")
    fig.tight_layout()
    fig.savefig(IMAGES / "occam1d_gv_valley_models.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    IMAGES.mkdir(parents=True, exist_ok=True)
    if BUILD_DIR.exists():
        shutil.rmtree(BUILD_DIR)

    batch, runs = build_and_invert()

    print("build_summary:")
    print(batch.summary())
    print("per_station_results:")
    for station in STATIONS:
        _, result = runs[station]
        print(
            f"{station}: status={result.convergence.value} "
            f"iterations={result.n_iterations} "
            f"initial_rms={result.initial.rms:.2f} "
            f"final_rms={result.final.rms:.3f}"
        )

    # Same three stations, run again through the batch API to demonstrate
    # invert_all's own return shape (a second, deliberately separate,
    # inversion pass -- see the "Batch inversion" section below).
    outcome = batch.invert_all(n_jobs=1, export_text=False, export_images=False)
    print("invert_all_results:")
    for station, summary in sorted(outcome["results"].items()):
        print(
            f"{station}: status={summary['status']} "
            f"iterations={summary['iterations']} "
            f"final_rms={summary['final_rms']:.3f}"
        )

    make_station_summary_figure(*runs["gv100"])
    make_valley_resistivity_comparison_figure(runs)

    shutil.rmtree(BUILD_DIR)
    print(f"images: {IMAGES.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
