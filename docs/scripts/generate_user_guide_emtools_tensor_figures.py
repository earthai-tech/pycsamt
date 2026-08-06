"""Generate the multi-panel figure bundles for the emtools tensor guide."""

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.emtools import (  # noqa: E402
    ensure_sites,
    plot_dimensionality_grid,
    plot_dimensionality_psection,
    plot_ellipticity_psection,
    plot_phase_tensor_psection,
    plot_phase_tensor_rose,
    plot_phase_tensor_skewmap,
    plot_phase_tensor_summary,
    plot_skew_ellipt_density,
    plot_theta_rose_grid,
    plot_theta_stability_stripe,
    plot_theta_vs_period,
)

IMAGES = ROOT / "docs/source/images/user_guide/emtools"
L18PLT = ROOT / "data/AMT/WILLY_DATA/L18PLT"


def make_simple_views_bundle() -> None:
    """Write the five quick-look diagnostic figures used by the "Simple
    Phase-Tensor Views" section: strike scatter, skew and ellipticity
    heatmaps, and the two dimensionality classifications.
    """
    sites = ensure_sites(L18PLT, recursive=True)

    plot_theta_vs_period(sites, recursive=False)
    plt.gcf().savefig(IMAGES / "user-guide-emtools-tensor-06-01.png", dpi=200, bbox_inches="tight")
    plt.close()

    plot_phase_tensor_skewmap(sites, recursive=False, axis_y="logperiod")
    plt.gcf().savefig(IMAGES / "user-guide-emtools-tensor-06-02.png", dpi=200, bbox_inches="tight")
    plt.close()

    plot_ellipticity_psection(sites, recursive=False)
    plt.gcf().savefig(IMAGES / "user-guide-emtools-tensor-06-03.png", dpi=200, bbox_inches="tight")
    plt.close()

    plot_dimensionality_psection(
        sites,
        skew_th=3.0,
        ellipt_th=0.2,
        recursive=False,
    )
    plt.gcf().savefig(IMAGES / "user-guide-emtools-tensor-06-04.png", dpi=200, bbox_inches="tight")
    plt.close()

    plot_dimensionality_grid(
        sites,
        skew_th=3.0,
        ellipt_th=0.2,
        recursive=False,
    )
    plt.gcf().savefig(IMAGES / "user-guide-emtools-tensor-06-05.png", dpi=200, bbox_inches="tight")
    plt.close()


def make_rose_stability_bundle() -> None:
    """Write the three rose/stability figures used by the "Rose And
    Stability Plots" section: the overall axial rose, the six-band rose
    grid, and the hue/saturation stability stripe.
    """
    sites = ensure_sites(L18PLT, recursive=True)

    plot_phase_tensor_rose(
        sites,
        band=(0.001, 10.0),
        bins=36,
        recursive=False,
    )
    plt.gcf().savefig(IMAGES / "user-guide-emtools-tensor-09-01.png", dpi=200, bbox_inches="tight")
    plt.close()

    plot_theta_rose_grid(
        sites,
        n_bands=6,
        bins=24,
        recursive=False,
    )
    plt.gcf().savefig(IMAGES / "user-guide-emtools-tensor-09-02.png", dpi=200, bbox_inches="tight")
    plt.close()

    plot_theta_stability_stripe(
        sites,
        win=5,
        recursive=False,
    )
    plt.gcf().savefig(IMAGES / "user-guide-emtools-tensor-09-03.png", dpi=200, bbox_inches="tight")
    plt.close()


def run_recommended_workflow() -> None:
    """Build the core pseudosection, rose grid, density plot, and summary
    figure used by the "Recommended Interpretation Workflow" section,
    all with the same period band and thresholds.
    """
    sites = ensure_sites(L18PLT, recursive=True)

    period_range = (0.001, 10.0)
    skew_threshold = 3.0
    ellipt_threshold = 0.2

    plot_phase_tensor_psection(
        sites,
        period_range=period_range,
        c_by="skew",
        skew_threshold=skew_threshold,
        mark_3d=True,
        recursive=False,
    )
    plt.gcf().savefig(IMAGES / "user-guide-emtools-tensor-24-01.png", dpi=200, bbox_inches="tight")
    plt.close()

    plot_theta_rose_grid(
        sites,
        n_bands=6,
        recursive=False,
    )
    plt.gcf().savefig(IMAGES / "user-guide-emtools-tensor-24-02.png", dpi=200, bbox_inches="tight")
    plt.close()

    plot_skew_ellipt_density(
        sites,
        band=period_range,
        recursive=False,
    )
    plt.gcf().savefig(IMAGES / "user-guide-emtools-tensor-24-03.png", dpi=200, bbox_inches="tight")
    plt.close()

    fig = plot_phase_tensor_summary(
        sites,
        period_range=period_range,
        skew_threshold=skew_threshold,
        ellipt_threshold=ellipt_threshold,
        recursive=False,
    )
    fig.savefig(IMAGES / "user-guide-emtools-tensor-24-04.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    make_simple_views_bundle()
    make_rose_stability_bundle()
    run_recommended_workflow()
