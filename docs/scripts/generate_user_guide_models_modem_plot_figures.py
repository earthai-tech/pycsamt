"""Figures for the ModEM plotting sections of the models/modem guide.

The bundled ``willy_27freq_watex_line02_sample`` run stalled at RMS ~3
and its station geometry is a single dense grid, which makes the
geo-referenced map and section plotters hard to read.  These figures
therefore use the bundled 21-station **Broken Hill** MT survey and its
converged 3-D ModEM model (``data/MT/broken-hill``; see that folder's
``README.md`` for the data source and citation).
"""

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.models.modem import InversionResult  # noqa: E402
from pycsamt.models.modem.plot import PlotDepthMap, PlotSection  # noqa: E402

IMAGES = ROOT / "docs/source/images/user_guide/models"
WORKDIR = ROOT / "data/MT/broken-hill/final-models"

# True grid-centre lat/lon.  The ModEM data file reduces longitude by
# 100 degrees (its header reads 41.53481).
ORIGIN_LAT = -31.95556
ORIGIN_LON = 141.53481


def make_depth_map_figures() -> None:
    """Depth-slice grid, hull-masked crop, and a smooth conductance map."""
    result = InversionResult(WORKDIR)

    fig = PlotDepthMap(
        result,
        depths={"(a) 1 km": 1000, "(b) 2 km": 2000,
                "(c) 3 km": 3000, "(d) 5 km": 5000},
        origin_lat=ORIGIN_LAT,
        origin_lon=ORIGIN_LON,
        rho_range=(1, 10000),
        mask_outside_hull=True,
        contours=[1000],
        scalebar=True,
        north_arrow=True,
        title="Broken Hill ModEM -- resistivity depth slices",
    ).plot()
    fig.savefig(
        IMAGES / "modem_depthmap_slices.png", dpi=200, bbox_inches="tight"
    )
    plt.close(fig)

    fig = PlotDepthMap(
        result,
        quantity="conductance",
        conductance_window=(1000, 6000),
        origin_lat=ORIGIN_LAT,
        origin_lon=ORIGIN_LON,
        render="gouraud",
        cmap="magma",
        norm="linear",
        rho_range=(10, 50),
        smooth_sigma=0.8,
        mask_outside_hull=True,
        station_color="white",
        cbar_orientation="horizontal",
        title="Broken Hill -- conductance, 1-6 km",
    ).plot()
    fig.savefig(
        IMAGES / "modem_depthmap_conductance.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(fig)


def make_section_figure() -> None:
    """One arbitrary-azimuth vertical resistivity section."""
    result = InversionResult(WORKDIR)
    fig = PlotSection(
        result=result,
        start_point=(-31.860, 141.490),   # NW
        end_point=(-32.020, 141.610),     # SE, across strike
        use_latlon=True,
        origin_lat=ORIGIN_LAT,
        origin_lon=ORIGIN_LON,
        n_samples=320,
        depth_max=8000.0,
        rho_min=1.0,
        rho_max=5000.0,
        show_station_names=False,
        station_tol=2500.0,
        figsize=(11.5, 4.6),
        title="Broken Hill ModEM -- section B-B' (NW -> SE)",
    ).plot()
    fig.savefig(
        IMAGES / "modem_section_broken_hill.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(fig)


if __name__ == "__main__":
    make_depth_map_figures()
    make_section_figure()
