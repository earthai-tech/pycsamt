"""Generate the recommended-workflow figures for the emtools strike guide."""

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.emtools import (  # noqa: E402
    ensure_sites,
    estimate_strike_consensus,
    plot_strike_profile,
    plot_strike_rose,
    rotate_to_strike,
)

IMAGES = ROOT / "docs/source/images/user_guide/emtools"
L18PLT = ROOT / "data/AMT/WILLY_DATA/L18PLT"


def run_strike_workflow() -> None:
    """Estimate, visualize, and conditionally rotate L18PLT onto its
    consensus strike, saving the profile and rose figures used by the
    "Recommended Workflow" section.
    """
    sites = ensure_sites(L18PLT, recursive=True)
    band = (0.001, 10.0)

    consensus = estimate_strike_consensus(
        sites,
        band=band,
        w_sweep=0.4,
        w_pt=0.6,
    )

    ax = plot_strike_profile(
        sites,
        method="consensus",
        band=band,
        sort_by="lat",
    )
    ax.figure.savefig(
        IMAGES / "user-guide-emtools-strike-18-01.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(ax.figure)

    fig = plot_strike_rose(
        sites,
        method="consensus",
        band=band,
        weight="inv_iqr",
    )
    fig.savefig(
        IMAGES / "user-guide-emtools-strike-18-02.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(fig)

    stable = consensus["iqr"].median() < 45.0
    if stable:
        rotate_to_strike(
            sites,
            method="consensus",
            band=band,
            inplace=False,
        )


if __name__ == "__main__":
    run_strike_workflow()
