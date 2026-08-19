"""Generate the recommended-workflow figures for the emtools tf guide."""

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.emtools import (  # noqa: E402
    ensure_sites,
    plot_induction_convention,
    plot_induction_map,
    plot_induction_multiperiod_map,
    plot_induction_rose,
    plot_induction_section,
    plot_tipper_hodograms,
    plot_tipper_polar,
)

IMAGES = ROOT / "docs/source/images/user_guide/emtools"
KAP03 = ROOT / "data/MT/kap03lmt_edis"


def run_tf_workflow() -> None:
    """Run station, period, and convention diagnostics for KAP03 in the
    order recommended on the tf.rst page: hodogram, polar view, map,
    convention comparison, short/long-period roses, period section, and
    a multi-period map. Each step reproduces one of the figures already
    shown individually above; this function is documented here as a
    single, self-contained script rather than as a separate set of
    figures.
    """
    sites = ensure_sites(KAP03, recursive=True)

    strongest_station = "kap151"
    periods = [25.0, 650.0, 2000.0, 17000.0]

    plot_tipper_hodograms(
        sites,
        station=strongest_station,
        bands=[(25.0, 200.0), (200.0, 2000.0), (2000.0, 17000.0)],
    )
    plt.close()

    plot_tipper_polar(
        sites,
        station=strongest_station,
        component="real",
    )
    plt.close()

    plot_induction_map(
        sites,
        period=2000.0,
        convention="park",
        show_real=True,
        show_imag=True,
        scale=4.0,
    )
    plt.close()

    plot_induction_convention(
        sites,
        period=650.0,
        station_labels=False,
    )
    plt.close()

    plot_induction_rose(
        sites,
        component="real",
        pband=(25.0, 200.0),
        title="Short-period induction azimuths",
    )
    plt.close()

    plot_induction_rose(
        sites,
        component="real",
        pband=(2000.0, 17000.0),
        title="Long-period induction azimuths",
    )
    plt.close()

    plot_induction_section(
        sites,
        component="abs",
        n_periods=30,
    )
    plt.close()

    plot_induction_multiperiod_map(
        sites,
        periods=periods,
        convention="park",
        arrow_scale=6.0,
        show_background_cbar=False,
        station_labels=False,
        xlabel="Longitude",
        ylabel="Latitude",
    )
    plt.close()


if __name__ == "__main__":
    run_tf_workflow()
