"""Generate the figures for the emtools ztem guide.

Reads the two synthetic ZTEM sample surveys committed under
``data/ZTEM/`` (see that directory's generator,
``scripts/airborne/generate_synthetic_survey_xml.py``, for the full
synthetic-data notice): ``gold_springs_nv`` (a real 7-line, 105-
station block survey, loosely inspired by Legault et al. 2012's Gold
Springs epithermal-gold ZTEM survey) and ``forrestania_wa`` (a single
15-station line, loosely inspired by Sattel and Witherly 2012's
Forrestania test-range survey). Neither is a vendor delivery or a
reproduction of the cited papers' actual field data.

Single-line figures (profile/pseudosection) all use
``gold_springs_nv``'s line 4 -- the survey's own central line, sitting
directly over the synthetic target's along-strike peak (target
amplitude follows a Gaussian centred on line 4, tapering toward lines
1 and 7 at the block's edges) -- rather than line 1, so the crossover
these figures demonstrate is the clean, strong one the target
actually produces, not the much weaker one an off-target edge line
shows.
"""

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.airborne.site import ensure_asites  # noqa: E402
from pycsamt.emtools.ztem import (  # noqa: E402
    plot_ztem_band_mask_psection,
    plot_ztem_divergence_profile,
    plot_ztem_divergence_psection_grid,
    plot_ztem_flight_lines,
    plot_ztem_map,
    plot_ztem_phase_rotation_profile,
    plot_ztem_tipper_profile,
)

IMAGES = ROOT / "docs/source/images/user_guide/emtools"
GOLD_SPRINGS = ROOT / "data/ZTEM/gold_springs_nv"
FORRESTANIA = ROOT / "data/ZTEM/forrestania_wa"

_TARGET_LINE_ID = "gold_springs_nv_L4"  # the survey's own central,
# strongest-response line -- see the module docstring for why this
# (not L1) is used for every single-line figure below.


def _load_sites():
    gs_all = ensure_asites(GOLD_SPRINGS)
    gs_l4 = gs_all.select(
        predicate=lambda s: (
            s.emtf.metadata.get("notes", {}).get("ZTEM", {}).get("LineId")
            == _TARGET_LINE_ID
        )
    )
    forrestania = ensure_asites(FORRESTANIA)
    return gs_all, gs_l4, forrestania


def make_tipper_profile(gs_l4) -> None:
    ax = plot_ztem_tipper_profile(gs_l4, frequency_hz=90.0)
    ax.figure.savefig(
        IMAGES / "user-guide-emtools-ztem-01.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(ax.figure)


def make_divergence_and_phase_rotation(gs_l4) -> None:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.0, 4.2))
    _ = plot_ztem_divergence_profile(gs_l4, frequency_hz=90.0, ax=ax1)
    _ = plot_ztem_phase_rotation_profile(
        gs_l4, frequency_hz=90.0, ax=ax2,
    )
    fig.tight_layout()
    fig.savefig(
        IMAGES / "user-guide-emtools-ztem-02.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(fig)


def make_divergence_psection_grid(gs_all) -> None:
    fig = plot_ztem_divergence_psection_grid(
        gs_all, max_lines=6, n_cols=3, cmap="seismic",
    )
    fig.savefig(
        IMAGES / "user-guide-emtools-ztem-03.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(fig)


def make_flight_lines(gs_all) -> None:
    ax = plot_ztem_flight_lines(gs_all)
    ax.figure.savefig(
        IMAGES / "user-guide-emtools-ztem-04.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(ax.figure)


def make_maps(gs_all) -> None:
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(13.0, 6.0), sharey=True,
    )
    _ = plot_ztem_map(
        gs_all, quantity="tipper", frequency_hz=90.0, ax=ax1,
    )
    _ = plot_ztem_map(
        gs_all, quantity="divergence", frequency_hz=90.0, ax=ax2,
    )
    # sharey=True already ties the two y-axes together; only the
    # left panel needs to keep its latitude tick labels/label.
    ax2.set_ylabel("")
    ax2.tick_params(axis="y", labelleft=False)
    fig.tight_layout()
    fig.savefig(
        IMAGES / "user-guide-emtools-ztem-05.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(fig)


def make_band_mask(forrestania) -> None:
    fig = plot_ztem_band_mask_psection(forrestania, band_hz=(40.0, 300.0))
    fig.savefig(
        IMAGES / "user-guide-emtools-ztem-06.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(fig)


def run_emtools_ztem_workflow() -> None:
    gs_all, gs_l4, forrestania = _load_sites()
    make_tipper_profile(gs_l4)
    make_divergence_and_phase_rotation(gs_l4)
    make_divergence_psection_grid(gs_all)
    make_flight_lines(gs_all)
    make_maps(gs_all)
    make_band_mask(forrestania)


if __name__ == "__main__":
    IMAGES.mkdir(parents=True, exist_ok=True)
    run_emtools_ztem_workflow()
    print("done ->", IMAGES)
