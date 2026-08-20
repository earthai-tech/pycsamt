"""Generate the figures for the airborne site user_guide page.

Reads four of the committed synthetic airborne sample surveys directly
through :func:`pycsamt.airborne.site.ensure_asites` -- one representative
site per response family (``afmag_original``, ``afmag_airmt``, ``ztem``,
``mobilemt``) for the response-family comparison, and three surveys
(one per technology) for the diagnostic panel, each rendered with its
own technology's real ``emtools`` plotting function. None of these are
vendor deliveries; see ``scripts/airborne/generate_synthetic_survey_xml.py``
for the synthetic-data notice.
"""

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.airborne.site import ensure_asites  # noqa: E402

IMAGES = ROOT / "docs/source/images/user_guide/airborne"
IMAGES.mkdir(parents=True, exist_ok=True)

AFMAG_ORIGINAL = ROOT / "data/AFMAG/abitibi_on"
AFMAG_AIRMT = ROOT / "data/AFMAG/yulong_belt_cn"
ZTEM_SURVEY = ROOT / "data/ZTEM/gold_springs_nv"
MOBILEMT_SURVEY = ROOT / "data/mobileMT/flammefjeld_greenland"


def _tensor_frobenius(arr: np.ndarray) -> np.ndarray:
    """Per-frequency Frobenius norm of a (nf, ..., ...) complex tensor."""
    axes = tuple(range(1, arr.ndim))
    return np.sqrt(np.sum(np.abs(arr) ** 2, axis=axes))


def make_response_family_grid() -> None:
    afmag_o = ensure_asites(AFMAG_ORIGINAL)[0]
    afmag_a = ensure_asites(AFMAG_AIRMT)[4]
    ztem = ensure_asites(ZTEM_SURVEY).get("GO_L4_008")
    mmt = ensure_asites(MOBILEMT_SURVEY)[0]

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.2))

    ax = axes[0, 0]
    ax.plot(afmag_o.freq, afmag_o.afmag_tilt_deg, "o", ms=8, color="C0")
    ax.axhline(0.0, color="0.6", lw=0.8)
    ax.set_title(f"afmag_original -- {afmag_o.name}")
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("Tilt angle [deg]")

    ax = axes[0, 1]
    mag = _tensor_frobenius(afmag_a.interstation_tensor)
    ax.semilogx(afmag_a.freq, mag, "o-", color="C1")
    ax.set_title(f"afmag_airmt -- {afmag_a.name}")
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel(r"$\|\mathbf{T}\|_F$ (interstation tensor)")

    ax = axes[1, 0]
    mag = _tensor_frobenius(ztem.tipper)
    ax.plot(ztem.freq, mag, "o-", color="C2")
    ax.set_title(f"ztem -- {ztem.name}")
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel(r"$\|\mathbf{T}\|_F$ (tipper)")

    ax = axes[1, 1]
    mag = _tensor_frobenius(mmt.admittance)
    ax.loglog(mmt.freq, mag, "o-", color="C3")
    ax.set_title(f"mobilemt -- {mmt.name}")
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel(r"$\|\mathbf{Y}\|_F$ (admittance)")

    fig.tight_layout()
    fig.savefig(
        IMAGES / "user-guide-airborne-site-01.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(fig)


def make_technology_diagnostic_panel() -> None:
    """One publication-style diagnostic per technology, side by side.

    Each panel calls the exact same ``plot_*`` function its own
    ``emtools`` guide page uses (:func:`~pycsamt.emtools.ztem.plot_ztem_map`,
    :func:`~pycsamt.emtools.afmag.plot_original_afmag_dual_frequency_profile`,
    :func:`~pycsamt.emtools.mobilemt.plot_mobilemt_conductivity_psection`)
    directly on the ``AirborneSites`` objects this page builds -- the
    point being that the container needs no conversion to reach a
    real, literature-style figure.
    """
    from pycsamt.emtools.afmag import plot_original_afmag_dual_frequency_profile
    from pycsamt.emtools.mobilemt import plot_mobilemt_conductivity_psection
    from pycsamt.emtools.ztem import plot_ztem_map

    ztem = ensure_asites(ZTEM_SURVEY)
    afmag_o = ensure_asites(AFMAG_ORIGINAL)
    mmt = ensure_asites(MOBILEMT_SURVEY)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16.5, 5.6))
    plot_ztem_map(ztem, ax=ax1)
    plot_original_afmag_dual_frequency_profile(afmag_o, ax=ax2)
    plot_mobilemt_conductivity_psection(mmt, ax=ax3)
    for ax in (ax1, ax2, ax3):
        ax.title.set_fontsize(10)

    fig.tight_layout(w_pad=2.5)
    fig.savefig(
        IMAGES / "user-guide-airborne-site-02.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(fig)


if __name__ == "__main__":
    make_response_family_grid()
    make_technology_diagnostic_panel()
