"""Generate the figures for the emtools afmag guide.

Reuses the bundled KAP03 long-period MT profile (real tipper, real
station coordinates) purely as a source of a genuine complex
``(Tx, Ty)`` tipper -- the same object type AFMAG's classical tilt-angle
readout maps onto. KAP03 is a ground MT survey, not an airborne AFMAG
survey; the motion-coupling physics section below uses a separate,
explicitly synthetic attitude time series, since no bundled dataset
carries INS/attitude data.
"""

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.emtools import (  # noqa: E402
    afmag_tilt_angles,
    ensure_sites,
    flag_motion_susceptible_band,
    motion_coupling_angle,
    motion_coupling_cosine,
    motion_susceptibility_table,
    plot_afmag_correction_comparison,
    plot_afmag_tilt_polar,
    plot_afmag_tilt_profile,
    plot_afmag_tilt_psection,
    plot_motion_coupling,
    plot_motion_susceptibility_map,
)

IMAGES = ROOT / "docs/source/images/user_guide/emtools"
KAP03 = ROOT / "data/MT/kap03lmt_edis"

# Illustrative mid-southern-latitude geomagnetic geometry for the KAP03
# footprint (roughly southern South Africa/Namibia/Botswana, ~-22 to
# -32 deg latitude). These are representative round values, not an
# IGRF lookup for a specific station and epoch -- for real survey work,
# obtain the true local inclination/declination from a proper
# geomagnetic reference model instead.
KAP03_INCLINATION_DEG = -62.0
KAP03_DECLINATION_DEG = -20.0

# Nominal platform attitude amplitude for the QC demonstration, chosen
# to match the small-angle regime Liu et al. (2018) themselves used for
# their literature comparison ("angular range of 10 degrees").
ROLL_AMPLITUDE_DEG = 8.0
PITCH_AMPLITUDE_DEG = 5.0


def _load_sites():
    return ensure_sites(KAP03, recursive=True, strict=False)


def make_tilt_profile(sites) -> None:
    ax = plot_afmag_tilt_profile(sites, component="real", period_s=650.0)
    ax.figure.savefig(
        IMAGES / "user-guide-emtools-afmag-01.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(ax.figure)


def make_tilt_psection(sites) -> None:
    ax = plot_afmag_tilt_psection(sites, component="resultant")
    ax.figure.savefig(
        IMAGES / "user-guide-emtools-afmag-02.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(ax.figure)


def make_tilt_polar(sites, station: str) -> None:
    ax = plot_afmag_tilt_polar(sites, station=station, component="real")
    ax.figure.savefig(
        IMAGES / "user-guide-emtools-afmag-03.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(ax.figure)


#: Inclination values swept in the sensitivity-grid figure, matching
#: the style of Liu et al. (2018) Fig. 2 (their own sweep is likewise
#: entirely synthetic/simulated, not field data).
SENSITIVITY_INCLINATIONS_DEG = (-30.0, -50.0, -70.0, -90.0)


def _build_motion_coupling_sensitivity_figure() -> plt.Figure:
    """Reproduce the qualitative structure of Liu et al. (2018) Fig. 2:
    theta(t) and cos(theta(t)) swept over roll, pitch, and yaw at
    several inclination values, with declination held at zero.
    """
    angles = np.linspace(-90.0, 90.0, 181)
    cmap = plt.get_cmap("viridis")
    colors = [
        cmap(x)
        for x in np.linspace(0.15, 0.9, len(SENSITIVITY_INCLINATIONS_DEG))
    ]

    fig, axes = plt.subplots(
        2, 3, figsize=(12.0, 6.5), sharex="col", sharey="row"
    )
    labels = ["roll", "pitch", "yaw"]
    for col, label in enumerate(labels):
        for inc, color in zip(SENSITIVITY_INCLINATIONS_DEG, colors):
            zeros = np.zeros_like(angles)
            if label == "roll":
                yaw_arr, pitch_arr, roll_arr = zeros, zeros, angles
            elif label == "pitch":
                yaw_arr, pitch_arr, roll_arr = zeros, angles, zeros
            else:
                yaw_arr, pitch_arr, roll_arr = angles, zeros, zeros
            theta = motion_coupling_angle(
                yaw_arr, pitch_arr, roll_arr, inc, 0.0
            )
            cos_theta = motion_coupling_cosine(
                yaw_arr, pitch_arr, roll_arr, inc, 0.0
            )
            axes[0, col].plot(
                angles, theta, color=color, lw=1.6, label=f"I={inc:g}°"
            )
            axes[1, col].plot(angles, cos_theta, color=color, lw=1.6)
        axes[0, col].set_title(f"{label} sweep", fontsize=10)
        axes[1, col].set_xlabel(f"{label} angle (deg)")
        for row in (0, 1):
            axes[row, col].grid(True, ls=":", alpha=0.35)

    axes[0, 0].set_ylabel(r"$\theta(t)$ (deg)")
    axes[1, 0].set_ylabel(r"$\cos\theta(t)$")
    axes[0, 0].legend(fontsize=7, loc="upper right", ncol=2)
    fig.suptitle(
        "Motion-coupling sensitivity to attitude and inclination\n"
        "(declination = 0°, cf. Liu et al. 2018 Fig. 2)"
    )
    fig.tight_layout()
    return fig


def make_motion_coupling_sensitivity(
    filename: str = "user-guide-emtools-afmag-04.png",
) -> None:
    fig = _build_motion_coupling_sensitivity_figure()
    fig.savefig(IMAGES / filename, dpi=200, bbox_inches="tight")
    plt.close(fig)


def make_motion_coupling() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(42)
    t = np.linspace(0.0, 20.0, 400)
    # Slow platform oscillation plus a little independent jitter per
    # channel, at the paper's own reported small-angle test amplitude.
    yaw = 15.0 * np.sin(0.15 * t) + rng.normal(0.0, 0.5, t.size)
    pitch = PITCH_AMPLITUDE_DEG * np.sin(0.31 * t + 0.6) + rng.normal(
        0.0, 0.3, t.size
    )
    roll = ROLL_AMPLITUDE_DEG * np.sin(0.47 * t + 1.1) + rng.normal(
        0.0, 0.3, t.size
    )
    ax = plot_motion_coupling(
        yaw,
        pitch,
        roll,
        KAP03_INCLINATION_DEG,
        KAP03_DECLINATION_DEG,
        x=t,
        xlabel="time (s)",
    )
    ax.figure.savefig(
        IMAGES / "user-guide-emtools-afmag-05.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(ax.figure)
    return yaw, pitch, roll


def make_susceptibility_map(sites) -> None:
    ax = plot_motion_susceptibility_map(
        sites,
        inclination=KAP03_INCLINATION_DEG,
        declination=KAP03_DECLINATION_DEG,
        roll_amplitude_deg=ROLL_AMPLITUDE_DEG,
        pitch_amplitude_deg=PITCH_AMPLITUDE_DEG,
    )
    ax.figure.savefig(
        IMAGES / "user-guide-emtools-afmag-06.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(ax.figure)


def make_correction_comparison(sites) -> None:
    fig = plot_afmag_correction_comparison(
        sites,
        inclination=KAP03_INCLINATION_DEG,
        declination=KAP03_DECLINATION_DEG,
        roll_amplitude_deg=ROLL_AMPLITUDE_DEG,
        pitch_amplitude_deg=PITCH_AMPLITUDE_DEG,
        threshold=0.1,
        band_hz=(0.005, 0.04),
    )
    fig.savefig(
        IMAGES / "user-guide-emtools-afmag-07.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(fig)


def run_afmag_workflow() -> None:
    """Run the recommended AFMAG diagnostic sequence on KAP03 end to
    end: tilt profile, tilt pseudosection, tilt polar view, motion-
    coupling sensitivity grid, one illustrative motion-coupling time
    series, motion-susceptibility map, and before/after/delta
    correction comparison. Each step reproduces one of the figures
    already shown individually above; this function is documented here
    as a single, self-contained script rather than as a separate set
    of figures.
    """
    sites = _load_sites()
    strongest_station = "kap151"

    plot_afmag_tilt_profile(sites, component="real", period_s=650.0)
    plt.close()

    plot_afmag_tilt_psection(sites, component="resultant")
    plt.close()

    plot_afmag_tilt_polar(sites, station=strongest_station, component="real")
    plt.close()

    plt.close(_build_motion_coupling_sensitivity_figure())

    rng = np.random.default_rng(42)
    t = np.linspace(0.0, 20.0, 400)
    yaw = 15.0 * np.sin(0.15 * t) + rng.normal(0.0, 0.5, t.size)
    pitch = PITCH_AMPLITUDE_DEG * np.sin(0.31 * t + 0.6) + rng.normal(
        0.0, 0.3, t.size
    )
    roll = ROLL_AMPLITUDE_DEG * np.sin(0.47 * t + 1.1) + rng.normal(
        0.0, 0.3, t.size
    )
    plot_motion_coupling(
        yaw,
        pitch,
        roll,
        KAP03_INCLINATION_DEG,
        KAP03_DECLINATION_DEG,
        x=t,
        xlabel="time (s)",
    )
    plt.close()

    plot_motion_susceptibility_map(
        sites,
        inclination=KAP03_INCLINATION_DEG,
        declination=KAP03_DECLINATION_DEG,
        roll_amplitude_deg=ROLL_AMPLITUDE_DEG,
        pitch_amplitude_deg=PITCH_AMPLITUDE_DEG,
    )
    plt.close()

    plot_afmag_correction_comparison(
        sites,
        inclination=KAP03_INCLINATION_DEG,
        declination=KAP03_DECLINATION_DEG,
        roll_amplitude_deg=ROLL_AMPLITUDE_DEG,
        pitch_amplitude_deg=PITCH_AMPLITUDE_DEG,
        threshold=0.1,
        band_hz=(0.005, 0.04),
    )
    plt.close()


if __name__ == "__main__":
    sites = _load_sites()
    make_tilt_profile(sites)
    make_tilt_psection(sites)
    make_tilt_polar(sites, "kap151")
    make_motion_coupling_sensitivity()
    make_motion_coupling()
    make_susceptibility_map(sites)
    make_correction_comparison(sites)
    run_afmag_workflow()
