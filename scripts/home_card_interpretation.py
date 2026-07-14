# -*- coding: utf-8 -*-
"""Home-page card image: "Interpretation & mapping".

Left panel: an apparent-resistivity pseudosection of the real L18PLT demo
line (the same sites the corrections gallery processes) — the map-view
product MapView renders interactively.

Right panel: interpretation — :class:`pycsamt.interp.HydroInterpreter`
classifies the shared synthetic section into stratigraphic logs, drawn as a
fence across the profile with each layer's own lithology colour.

Output: docs/source/_static/images/home/card-interpretation.png

Usage (any cwd):
    python scripts/home_card_interpretation.py
"""
import os
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

os.environ.setdefault("PYCSAMT_DOCS_BUILD", "1")
os.environ.setdefault("PYCSAMT_DOCS_REPO_ROOT", str(ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "docs" / "examples" / "corrections"))
sys.path.insert(0, str(ROOT / "docs" / "examples" / "interp"))

from _corr_data import curves, demo_line  # noqa: E402
from _interp_data import demo_model  # noqa: E402

from pycsamt.interp import HydroInterpreter  # noqa: E402

INK = "#1c2941"
FMIN = 2.5e-3

STYLE = {
    "font.size": 9,
    "axes.linewidth": 0.8,
    "axes.edgecolor": "#c8d0dc",
    "xtick.color": "#5c677d",
    "ytick.color": "#5c677d",
    "axes.labelcolor": INK,
    "axes.titlesize": 10.5,
    "axes.titleweight": "bold",
    "axes.titlecolor": INK,
}


def edges(centers):
    c = np.asarray(centers, dtype=float)
    mid = np.sqrt(c[:-1] * c[1:]) if (c > 0).all() else (c[:-1] + c[1:]) / 2
    first = c[0] * c[0] / mid[0] if (c > 0).all() else 2 * c[0] - mid[0]
    last = c[-1] * c[-1] / mid[-1] if (c > 0).all() else 2 * c[-1] - mid[-1]
    return np.concatenate([[first], mid, [last]])


def pseudosection():
    """Station x frequency grid of apparent resistivity for L18PLT."""
    raw = curves(demo_line("L18PLT"), "rho")
    stations = sorted(raw)
    f0 = np.asarray(raw[stations[0]][0], dtype=float)
    keep = f0 >= FMIN
    f = f0[keep]
    grid = np.full((keep.sum(), len(stations)), np.nan)
    for j, s in enumerate(stations):
        fs, rs = np.asarray(raw[s][0], float), np.asarray(raw[s][1], float)
        # interpolate in log-log onto the shared frequency axis
        m = np.isfinite(rs) & (rs > 0) & np.isfinite(fs)
        order = np.argsort(fs[m])
        grid[:, j] = 10 ** np.interp(
            np.log10(f),
            np.log10(fs[m][order]),
            np.log10(rs[m][order]),
        )
    return stations, f, grid


def main():
    plt.rcParams.update(STYLE)
    fig, (ax_l, ax_r) = plt.subplots(
        1,
        2,
        figsize=(6.4, 3.0),
        dpi=200,
        gridspec_kw={"width_ratios": [1.3, 1.0]},
    )
    fig.patch.set_facecolor("white")

    # -- mapping: apparent-resistivity pseudosection ------------------------
    stations, f, grid = pseudosection()
    x_edges = np.arange(len(stations) + 1) - 0.5
    mesh = ax_l.pcolormesh(
        x_edges,
        edges(f),
        grid,
        norm=mcolors.LogNorm(vmin=10.0, vmax=1e4),
        cmap="viridis",
        rasterized=True,
    )
    ax_l.set_yscale("log")
    ax_l.set_ylim(f.min(), f.max())  # high frequency (shallow) at the top
    ax_l.set_title("mapping — ρ$_a$ pseudosection", pad=6)
    ax_l.set_xlabel("station (L18)", fontsize=8.5)
    ax_l.set_ylabel("frequency (Hz)", fontsize=8.5)
    ax_l.set_xticks(np.arange(0, len(stations), 6))
    ax_l.tick_params(labelsize=7.5, length=2.5)
    cbar = fig.colorbar(mesh, ax=ax_l, fraction=0.055, pad=0.03)
    cbar.set_label(r"$\rho_a$ ($\Omega\cdot$m)", fontsize=7.5)
    cbar.ax.tick_params(labelsize=6.5, length=2)

    # -- interpretation: classified stratigraphic fence ----------------------
    hydro = HydroInterpreter(
        water_table_depth=20.0,
        aquifer_range=(30.0, 300.0),
        clay_max=20.0,
        min_zone_thickness=14.0,
    ).fit(demo_model())
    logs = hydro.logs

    DEPTH_MAX = 118.0  # below this the synthetic basement classifies noisily
    width = 0.62 * (logs[1].station_x - logs[0].station_x)
    thickness_by_lith = defaultdict(float)
    color_by_lith = {}
    for log in logs:
        for layer in log.layers:
            top = min(layer.top, DEPTH_MAX)
            bot = min(layer.bottom, DEPTH_MAX)
            if bot <= top:
                continue
            ax_r.add_patch(
                mpatches.Rectangle(
                    (log.station_x - width / 2, top),
                    width,
                    bot - top,
                    facecolor=layer.color,
                    edgecolor="white",
                    linewidth=0.4,
                )
            )
            thickness_by_lith[layer.lithology] += bot - top
            color_by_lith[layer.lithology] = layer.color

    ax_r.set_xlim(logs[0].station_x - width, logs[-1].station_x + width)
    ax_r.set_ylim(DEPTH_MAX, 0)
    ax_r.set_title("interpretation — strat. logs", pad=6)
    ax_r.set_xlabel("distance (m)", fontsize=8.5)
    ax_r.set_ylabel("depth (m)", fontsize=8.5)
    ax_r.tick_params(labelsize=7.5, length=2.5)

    top_liths = sorted(
        thickness_by_lith, key=thickness_by_lith.get, reverse=True
    )[:4]
    ax_r.legend(
        handles=[
            mpatches.Patch(color=color_by_lith[t], label=t)
            for t in top_liths
        ],
        fontsize=6.2,
        loc="lower right",
        frameon=True,
        framealpha=0.85,
        edgecolor="none",
    )

    fig.tight_layout(pad=0.6, w_pad=1.2)
    out = ROOT / "docs/source/_static/images/home/card-interpretation.png"
    fig.savefig(out, dpi=200, facecolor="white", bbox_inches="tight")
    print(f"saved: {out} ({out.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
