"""Generate the figures for the emtools mobilemt guide.

Reads the two synthetic MobileMT sample surveys committed under
``data/mobileMT/`` (see ``scripts/airborne/generate_synthetic_survey_xml.py``
for the full synthetic-data notice): ``flammefjeld_greenland`` (a single
12-station line, loosely inspired by Zhdanov et al. 2024's Climax-style
porphyry Mo-Cu breccia-pipe survey in East Greenland) and
``timiskaming_kimberlite_on`` (a single 12-station line, loosely
inspired by Prikhodko et al. 2022's KL-22 kimberlite-pipe survey at
Lake Timiskaming, Ontario). Neither is a vendor delivery or a
reproduction of the cited papers' actual field data.
"""

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.emtools.mobilemt import (  # noqa: E402
    admittance_determinant_table,
    admittance_table,
    ensure_mobilemt_dataset,
    mask_outside_mobilemt_band,
    plot_mobilemt_admittance_profile,
    plot_mobilemt_conductivity_psection,
    plot_mobilemt_skew_profile,
)

IMAGES = ROOT / "docs/source/images/user_guide/emtools"
FLAMMEFJELD = ROOT / "data/mobileMT/flammefjeld_greenland"
TIMISKAMING = ROOT / "data/mobileMT/timiskaming_kimberlite_on"


def _load_datasets():
    flammefjeld = ensure_mobilemt_dataset(FLAMMEFJELD)
    timiskaming = ensure_mobilemt_dataset(TIMISKAMING)
    return flammefjeld, timiskaming


def make_admittance_profile(flammefjeld) -> None:
    ax = plot_mobilemt_admittance_profile(
        flammefjeld, component="det", part="abs", frequency_hz=390.0,
    )
    ax.figure.savefig(
        IMAGES / "user-guide-emtools-mobilemt-01.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(ax.figure)


def make_theoretical_vs_native(flammefjeld) -> None:
    # Shared color limits so the two panels are directly comparable --
    # they plot the same physical quantity from two independent
    # sources (pyCSAMT's own determinant formula vs. the
    # vendor-delivered native field recovered from the XML notes).
    df_th = admittance_determinant_table(flammefjeld)
    vals = np.concatenate(
        [
            df_th["theoretical_sigma_a_Sm"].to_numpy(),
            df_th["apparent_conductivity_native_Sm"].dropna().to_numpy(),
        ]
    )
    clim = tuple(np.percentile(vals, [2.0, 98.0]))

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(15.0, 5.2), sharey=True,
    )
    _ = plot_mobilemt_conductivity_psection(
        flammefjeld, source="theoretical", clim=clim, ax=ax1,
    )
    _ = plot_mobilemt_conductivity_psection(
        flammefjeld, source="native", clim=clim, ax=ax2,
    )
    ax2.set_ylabel("")
    ax2.tick_params(axis="y", labelleft=False)
    fig.tight_layout()
    fig.savefig(
        IMAGES / "user-guide-emtools-mobilemt-02.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(fig)


def make_skew_profile(timiskaming) -> None:
    ax = plot_mobilemt_skew_profile(timiskaming, frequency_hz=770.0)
    ax.figure.savefig(
        IMAGES / "user-guide-emtools-mobilemt-03.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(ax.figure)


def make_band_mask_comparison(timiskaming) -> None:
    masked = mask_outside_mobilemt_band(timiskaming, band_hz=(100.0, 3000.0))

    df_full = admittance_table(timiskaming)
    df_masked = admittance_table(masked)
    vals = np.concatenate(
        [
            df_full["apparent_conductivity_native_Sm"].dropna().to_numpy(),
            df_masked["apparent_conductivity_native_Sm"].dropna().to_numpy(),
        ]
    )
    clim = tuple(np.percentile(vals, [2.0, 98.0]))

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(15.0, 5.2), sharey=True,
    )
    _ = plot_mobilemt_conductivity_psection(
        timiskaming, source="native", clim=clim, ax=ax1,
    )
    ax1.set_title("Full band (25-12,000 Hz)", fontsize=10)
    _ = plot_mobilemt_conductivity_psection(
        masked, source="native", clim=clim, ax=ax2,
    )
    ax2.set_ylabel("")
    ax2.tick_params(axis="y", labelleft=False)
    ax2.set_title("Masked to 100-3,000 Hz", fontsize=10)
    fig.tight_layout()
    fig.savefig(
        IMAGES / "user-guide-emtools-mobilemt-04.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(fig)


def run_emtools_mobilemt_workflow() -> None:
    flammefjeld, timiskaming = _load_datasets()
    make_admittance_profile(flammefjeld)
    make_theoretical_vs_native(flammefjeld)
    make_skew_profile(timiskaming)
    make_band_mask_comparison(timiskaming)


if __name__ == "__main__":
    IMAGES.mkdir(parents=True, exist_ok=True)
    run_emtools_mobilemt_workflow()
    print("done ->", IMAGES)
