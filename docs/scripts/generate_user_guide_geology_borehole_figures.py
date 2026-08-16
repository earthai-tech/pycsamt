"""Generate the borehole figures for the geology borehole guide."""

from pathlib import Path
import sys
import tempfile

import matplotlib

matplotlib.use("Agg")

import matplotlib.cm as cm  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.geology import Borehole, Interval  # noqa: E402

IMAGES = ROOT / "docs/source/images/user_guide/geology"

# Same minimal LAS 2.0 fixture shown on the geology/borehole.rst page.
_BH02_LAS = """\
~VERSION INFORMATION
VERS.   2.0 : CWLS LOG ASCII STANDARD - VERSION 2.0
WRAP.    NO : ONE LINE PER DEPTH STEP
~WELL INFORMATION
WELL.        BH02 : WELL NAME
NULL.    -9999.25 : NULL VALUE
~CURVE INFORMATION
DEPT.M       : DEPTH
RESD.OHMM    : DEEP RESISTIVITY
LITH.        : LITHOLOGY CODE
~ASCII
0.0     450.0   1
2.0     430.0   1
4.0     410.0   1
6.0      60.0   2
8.0      45.0   2
10.0     40.0   2
12.0    190.0   3
14.0    180.0   3
"""

_BH02_RAW_DEPTHS = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0]
_BH02_RAW_RESD = [450.0, 430.0, 410.0, 60.0, 45.0, 40.0, 190.0, 180.0]

_BH01_COLORS = ["#B5651D", "#C9A66B", "#8E8E8E", "#4A4A4A"]
_BH02_COLORS = ["#2E86AB", "#A23B72", "#F18F01"]


def _make_bh01() -> Borehole:
    """BH01, built directly from ``Interval`` objects (no CSV/LAS)."""
    intervals = [
        Interval(top=0.0, bottom=8.0, lithology="lateritic soil", resistivity=450.0),
        Interval(top=8.0, bottom=31.0, lithology="clayey sand", resistivity=42.0),
        Interval(top=31.0, bottom=67.0, lithology="weathered granite", resistivity=185.0),
        Interval(top=67.0, bottom=120.0, lithology="fresh granite", resistivity=3200.0),
    ]
    return Borehole("BH01", x=500.0, intervals=intervals, collar_elevation=238.4)


def _make_bh02() -> Borehole:
    """BH02, read from the minimal LAS 2.0 fixture shown on the page."""
    with tempfile.TemporaryDirectory() as tmp:
        las_path = Path(tmp) / "BH02.las"
        las_path.write_text(_BH02_LAS, encoding="utf-8")
        return Borehole.from_las(las_path, x=750.0)


def make_borehole_log_comparison() -> None:
    """Write the two-panel BH01/BH02 well-log figure.

    Left panel: BH01, drawn as a classic resistivity-vs-depth log track
    colored by lithology. Right panel: BH02, showing the raw 2 m samples
    alongside the three intervals ``from_las`` grouped them into.
    """
    bh01 = _make_bh01()
    bh02 = _make_bh02()

    fig, axes = plt.subplots(1, 2, figsize=(9, 6))

    ax = axes[0]
    for iv, color in zip(bh01.intervals, _BH01_COLORS):
        ax.fill_betweenx(
            [iv.top, iv.bottom], 0, iv.resistivity,
            color=color, edgecolor="0.2", linewidth=0.7,
        )
        ax.text(
            iv.resistivity * 1.15, (iv.top + iv.bottom) / 2,
            f"{iv.lithology}\n{iv.resistivity:.0f} Ohm.m",
            va="center", fontsize=8,
        )
    ax.set_xscale("log")
    ax.set_xlim(1, 2e4)
    ax.set_ylim(bh01.max_depth, 0.0)
    ax.set_xlabel(r"Resistivity ($\Omega\,\mathrm{m}$, log scale)")
    ax.set_ylabel("Depth (m)")
    ax.set_title(f"{bh01.name} -- from_csv / direct construction")
    ax.grid(axis="x", which="both", alpha=0.3)

    ax = axes[1]
    ax.scatter(
        _BH02_RAW_RESD, _BH02_RAW_DEPTHS,
        s=30, color="0.15", zorder=3, label="raw 2 m samples",
        edgecolor="white", linewidth=0.6,
    )
    for iv, color in zip(bh02.intervals, _BH02_COLORS):
        ax.fill_betweenx(
            [iv.top, iv.bottom], 0, iv.resistivity,
            color=color, alpha=0.55, edgecolor="0.2", linewidth=0.7,
        )
        ax.text(
            iv.resistivity * 1.15, (iv.top + iv.bottom) / 2,
            f"code {iv.lithology}\nmean {iv.resistivity:.0f}",
            va="center", fontsize=8,
        )
    ax.set_xscale("log")
    ax.set_xlim(1, 2e3)
    ax.set_ylim(bh02.max_depth, 0.0)
    ax.set_xlabel(r"Resistivity ($\Omega\,\mathrm{m}$, log scale)")
    ax.set_title(f"{bh02.name} -- from_las grouping")
    ax.grid(axis="x", which="both", alpha=0.3)
    ax.legend(loc="lower right", fontsize=8)

    fig.tight_layout()
    fig.savefig(IMAGES / "borehole_log_comparison.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def make_borehole_profile_positions() -> None:
    """Write the three-borehole profile-position figure.

    Three independently constructed boreholes (BH01, plus a shallow
    BH03 and a deeper BH04), drawn as vertical log strips at their own
    ``x`` -- the picture a set of ground-truth control points along a
    survey line actually looks like before any of them reach a
    calibrator.
    """
    bh01 = _make_bh01()
    bh03 = Borehole("BH03", x=150.0, intervals=[
        Interval(top=0.0, bottom=6.0, lithology="topsoil", resistivity=300.0),
        Interval(top=6.0, bottom=35.0, lithology="saprolite", resistivity=120.0),
        Interval(top=35.0, bottom=60.0, lithology="fresh basement", resistivity=8000.0),
    ])
    bh04 = Borehole("BH04", x=900.0, intervals=[
        Interval(top=0.0, bottom=15.0, lithology="alluvium", resistivity=30.0),
        Interval(top=15.0, bottom=45.0, lithology="aquifer sand", resistivity=90.0),
        Interval(top=45.0, bottom=60.0, lithology="clay aquitard", resistivity=8.0),
        Interval(top=60.0, bottom=90.0, lithology="bedrock", resistivity=4000.0),
    ])
    boreholes = [bh03, bh01, bh04]

    all_liths: list[str] = []
    for bh in boreholes:
        for iv in bh.intervals:
            if iv.lithology not in all_liths:
                all_liths.append(iv.lithology)
    palette = cm.get_cmap("tab20").colors
    lith_color = {name: palette[i % len(palette)] for i, name in enumerate(all_liths)}

    fig, ax = plt.subplots(figsize=(8, 6))
    width = 60.0
    for bh in boreholes:
        for iv in bh.intervals:
            ax.fill_betweenx(
                [iv.top, iv.bottom], bh.x - width / 2, bh.x + width / 2,
                color=lith_color[iv.lithology], edgecolor="0.25", linewidth=0.6,
            )
        ax.text(
            bh.x, -3.0, f"{bh.name}\nx={bh.x:.0f} m",
            ha="center", va="bottom", fontsize=9, fontweight="bold",
        )

    max_depth = max(bh.max_depth for bh in boreholes)
    ax.set_xlim(0, 1050)
    ax.set_ylim(max_depth + 5, -12)
    ax.set_xlabel("Profile position x (m)")
    ax.set_ylabel("Depth (m)")
    ax.set_title("Three boreholes positioned along one survey profile")
    ax.grid(axis="both", alpha=0.25)

    handles = [plt.Rectangle((0, 0), 1, 1, color=lith_color[name]) for name in all_liths]
    ax.legend(
        handles, all_liths, loc="upper center", bbox_to_anchor=(0.5, -0.12),
        ncol=4, fontsize=8, frameon=False,
    )

    fig.tight_layout()
    fig.savefig(
        IMAGES / "borehole_profile_positions.png", dpi=200, bbox_inches="tight"
    )
    plt.close(fig)


def make_borehole_tres_sampling() -> None:
    """Write the ``tres_column()`` sampling figure for BH01.

    Draws BH01's log as a resistivity-vs-depth step curve and marks
    where ``tres_column(z)`` samples it for a hypothetical model's
    depth cells, distinguishing a found value from ``nan`` below
    ``max_depth`` -- what a calibrator sees when it reads this log
    cell by cell.
    """
    bh01 = _make_bh01()
    z = np.array([4.0, 20.0, 50.0, 100.0, 150.0])
    tres = bh01.tres_column(z)

    fig, ax = plt.subplots(figsize=(7, 5))

    for iv in bh01.intervals:
        ax.hlines(iv.resistivity, iv.top, iv.bottom, color="#2E4053", linewidth=3, zorder=2)
    for i in range(len(bh01.intervals) - 1):
        top = bh01.intervals[i + 1].top
        r0 = bh01.intervals[i].resistivity
        r1 = bh01.intervals[i + 1].resistivity
        ax.vlines(top, min(r0, r1), max(r0, r1), color="#2E4053", linewidth=1.2,
                  linestyle="--", zorder=2)
        ax.vlines(top, 1, min(r0, r1), color="0.75", linewidth=0.8, zorder=1, linestyle=":")

    found = np.isfinite(tres)
    ax.scatter(z[found], tres[found], color="#C0392B", zorder=4, s=50,
               label="tres_column() -- found")
    ax.scatter(z[~found], np.full((~found).sum(), 1.5), color="#C0392B",
               zorder=4, s=60, marker="x", linewidth=2,
               label="tres_column() -- nan (below max_depth)")

    for zi, ti in zip(z, tres):
        label = f"z={zi:.0f}\n{'nan' if not np.isfinite(ti) else f'{ti:.0f}'}"
        y = 1.5 if not np.isfinite(ti) else ti
        ax.annotate(label, (zi, y), textcoords="offset points", xytext=(0, 10),
                    ha="center", fontsize=8)

    ax.axvline(bh01.max_depth, color="0.4", linestyle="-.", linewidth=1)
    ax.text(bh01.max_depth + 2, 2000, f"max_depth={bh01.max_depth:.0f} m",
            fontsize=8, color="0.3")

    ax.set_yscale("log")
    ax.set_ylim(1, 9000)
    ax.set_xlim(-5, 160)
    ax.set_xlabel("Depth z (m)")
    ax.set_ylabel(r"Resistivity ($\Omega\,\mathrm{m}$, log scale)")
    ax.set_title("BH01.tres_column(z) sampled against the logged intervals", pad=14)
    ax.grid(alpha=0.25)
    ax.legend(loc="upper left", fontsize=8)

    fig.tight_layout()
    fig.savefig(IMAGES / "borehole_tres_sampling.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    IMAGES.mkdir(parents=True, exist_ok=True)
    make_borehole_log_comparison()
    make_borehole_profile_positions()
    make_borehole_tres_sampling()
