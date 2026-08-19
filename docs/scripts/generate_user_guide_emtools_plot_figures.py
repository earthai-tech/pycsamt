"""Generate the multi-figure report-bundle image for the emtools plot guide."""

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.emtools import (  # noqa: E402
    ensure_sites,
    plot_raw_sites_1d,
    plot_response_overview,
    plot_sites_compare,
    plot_sites_fit_grid,
    plot_sites_panels,
    smooth_mavg,
)

IMAGES = ROOT / "docs/source/images/user_guide/emtools"
L18PLT = ROOT / "data/AMT/WILLY_DATA/L18PLT"
GV_DATA = ROOT / "data/gv_data/gv_final_edi"


def make_plot_report_bundle() -> None:
    """Write the four report figures used by the "Save A Multi-Figure
    Plot Bundle" section: a station-panel grid, raw full-tensor panels,
    a raw-vs-smoothed comparison, and a measured-vs-predicted fit grid,
    all built from the same L18PLT survey and station subset.
    """
    survey = ensure_sites(L18PLT, strict=True)
    stations = ["18-001A", "18-007U", "18-016A", "18-018A"]
    processed = smooth_mavg(survey, k=5)

    fig = plot_sites_panels(
        survey,
        stations=stations,
        components=("xy", "yx"),
        ncols=4,
        show_legend=True,
    )
    fig.savefig(IMAGES / "user-guide-emtools-plot-15-01.png", dpi=200)
    plt.close(fig)

    fig = plot_raw_sites_1d(
        survey,
        stations=stations[:3],
        components=("xx", "xy", "yx", "yy"),
        ncols_groups=3,
    )
    fig.savefig(IMAGES / "user-guide-emtools-plot-15-02.png", dpi=200)
    plt.close(fig)

    fig = plot_sites_compare(
        survey,
        processed,
        stations=stations[:3],
        components=("xy", "yx"),
        labels=("raw", "smoothed"),
        ncols_groups=3,
    )
    fig.savefig(IMAGES / "user-guide-emtools-plot-15-03.png", dpi=200)
    plt.close(fig)

    fig = plot_sites_fit_grid(
        survey,
        processed,
        stations=stations[:2],
        components=("xy", "yx"),
        ncols_groups=2,
    )
    fig.savefig(IMAGES / "user-guide-emtools-plot-15-04.png", dpi=200)
    plt.close(fig)


def make_response_overview_figures() -> None:
    """Write the two figures used by the "Full-Response Station Overview"
    section: the default single-station overview (Gabbs Valley gv100)
    and the MTpy-style phimin ellipse-colouring variant (gv130).
    """
    survey = ensure_sites(GV_DATA, strict=True)

    fig = plot_response_overview(survey, station="gv100")
    fig.savefig(IMAGES / "user-guide-emtools-plot-16.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    fig = plot_response_overview(
        survey,
        station="gv130",
        c_by="phimin_deg",
        clim=(0.0, 90.0),
        cmap="turbo",
    )
    fig.savefig(IMAGES / "user-guide-emtools-plot-17.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    IMAGES.mkdir(parents=True, exist_ok=True)
    make_plot_report_bundle()
    make_response_overview_figures()
