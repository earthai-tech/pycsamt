"""Generate the figures for the emtools inspect-and-QC user guide.

Currently covers the MTPy-style apparent-resistivity / phase
pseudo-section (:class:`~pycsamt.emtools.PlotResPhasePseudoSection`).
Run from the repository root::

    python docs/scripts/generate_user_guide_emtools_inspect_figures.py
"""

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.emtools import (  # noqa: E402
    PlotResPhasePseudoSection,
    ensure_sites,
)

IMAGES = ROOT / "docs/source/images/user_guide/emtools"
L18PLT = ROOT / "data/AMT/WILLY_DATA/L18PLT"
L22PLT = ROOT / "data/AMT/WILLY_DATA/L22PLT"
KAP03 = ROOT / "data/MT/kap03lmt_edis"


def _save(fig, name: str, dpi: int = 150) -> None:
    fig.savefig(IMAGES / name, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def run_resphase_pseudosection() -> None:
    """Save the pseudo-section figures used by the
    "MTPy-Style Resistivity and Phase Pseudo-Sections" section.
    """
    line18 = ensure_sites(L18PLT, recursive=True, strict=True)
    line22 = ensure_sites(L22PLT, recursive=True, strict=True)
    kap03 = ensure_sites(KAP03, recursive=True, strict=True)

    # 15 -- one survey, every component auto-detected (4 columns),
    #       station markers on top (the default)
    _save(
        PlotResPhasePseudoSection(
            line18,
            title="L18PLT -- apparent resistivity and phase",
        ).plot(),
        "user-guide-emtools-inspect-15.png",
        dpi=160,
    )

    # 16 -- two AMT lines stacked, one shared period axis and colour
    #       scale, lettered panels
    _save(
        PlotResPhasePseudoSection(
            {"L18PLT": line18, "L22PLT": line22},
            components=["xy", "yx"],
            panel_labels=True,
        ).plot(),
        "user-guide-emtools-inspect-16.png",
        dpi=160,
    )

    # 17 -- KAP03 long-period line: resistivity twice the height of
    #       phase, station markers above each panel
    _save(
        PlotResPhasePseudoSection(
            kap03,
            components=["xy", "yx"],
            res_phase_ratio=2.0,
            station_side="top",
            title="KAP03 (LMT) -- apparent resistivity and phase",
        ).plot(),
        "user-guide-emtools-inspect-17.png",
        dpi=160,
    )

    # 18 -- an LMT line above an AMT line: different bands, so each
    #       group keeps its own period window
    _save(
        PlotResPhasePseudoSection(
            {"KAP03 (LMT)": kap03, "L18PLT (AMT)": line18},
            components=["xy", "yx"],
            res_phase_ratio=2.0,
            share_period=False,
            panel_labels=True,
        ).plot(),
        "user-guide-emtools-inspect-18.png",
        dpi=160,
    )


if __name__ == "__main__":
    run_resphase_pseudosection()
