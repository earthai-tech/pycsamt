"""Generate the figure for the airborne data_model user_guide page.

Builds one small, deterministic synthetic flight line (a hill under a
radar-altimeter clearance hold) purely from
:class:`pycsamt.airborne.NavigationTrack` -- no committed sample survey
is read here, since this page's point is the container mechanics
themselves, not reading real data (that is :doc:`site`'s job).
"""

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.airborne import NavigationTrack  # noqa: E402

IMAGES = ROOT / "docs/source/images/user_guide/airborne"
IMAGES.mkdir(parents=True, exist_ok=True)


def make_navigation_profile() -> None:
    n = 20
    x = np.arange(n, dtype=float) * 50.0
    terrain = 1500.0 + 80.0 * np.exp(-((x - 500.0) / 200.0) ** 2)
    rng = np.random.default_rng(7)
    target_clearance = 90.0
    platform = terrain + target_clearance + rng.normal(0.0, 3.0, size=n)
    nav = NavigationTrack(
        sample_ids=tuple(f"S{i:02d}" for i in range(n)),
        easting=x,
        northing=np.zeros(n),
        terrain_elevation=terrain,
        platform_elevation=platform,
        heading=np.full(n, 90.0),
    )
    cv = nav.clearance_values

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(9.0, 6.2), sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1]},
    )
    ax1.fill_between(x, 0, terrain, color="#c9b58a", alpha=0.6, label="Terrain")
    ax1.plot(x, terrain, color="#8a6d3b", lw=1.5)
    ax1.plot(x, platform, color="#1f6feb", lw=1.8, label="Platform")
    ax1.set_ylim(terrain.min() - 20, platform.max() + 20)
    ax1.set_ylabel("Elevation [m]")
    ax1.legend(loc="upper right", fontsize=9)
    ax1.set_title("Terrain-following flight profile, target clearance 90 m")

    ax2.plot(x, cv, color="0.2", lw=1.2, marker="o", ms=3)
    ax2.axhline(target_clearance, color="0.5", ls="--", lw=1.0, label="Target clearance")
    ax2.fill_between(x, cv, target_clearance, color="#1f6feb", alpha=0.15)
    ax2.set_ylabel("Clearance [m]")
    ax2.set_xlabel("Along-line distance [m]")
    ax2.legend(loc="upper right", fontsize=9)

    fig.tight_layout()
    fig.savefig(
        IMAGES / "user-guide-airborne-data-model-01.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(fig)


if __name__ == "__main__":
    make_navigation_profile()
