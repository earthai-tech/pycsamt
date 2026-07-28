"""Home-page card image: "Forward modelling".

Builds three :class:`pycsamt.forward.LayeredModel` scenarios, runs the real
1-D MT forward (:class:`pycsamt.forward.MT1DForward`) on each, and renders a
model-to-response figure sized for the landing-page feature card.

Output: docs/source/_static/images/home/card-forward.png

Usage (any cwd):
    python scripts/home_card_forward.py
"""

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

os.environ.setdefault("PYCSAMT_DOCS_BUILD", "1")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(ROOT))

from pycsamt.forward import LayeredModel, MT1DForward  # noqa: E402

BLUE, ORANGE, GOLD = "#3e65b0", "#f15a29", "#d99114"
INK = "#1c2941"

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

MODELS = [
    (
        LayeredModel(
            [1_000.0, 20.0, 5.0, 300.0],
            [200.0, 600.0, 1_500.0],
            name="sedimentary",
        ),
        BLUE,
    ),
    (
        LayeredModel(
            [500.0, 8.0, 250.0, 3_000.0],
            [100.0, 400.0, 2_500.0],
            name="geothermal",
        ),
        ORANGE,
    ),
    (
        LayeredModel(
            [800.0, 8_000.0, 600.0], [2_000.0, 15_000.0], name="crystalline"
        ),
        GOLD,
    ),
]

FREQS = np.logspace(-3, 4, 40)
DEPTH_MAX_KM = 5.0


def step_profile(model):
    """Resistivity-vs-depth step arrays (depth in km) for plotting."""
    rho = list(model.resistivity)
    tops = [0.0] + list(np.cumsum(model.thickness))
    xs, ys = [], []
    for r, top, bot in zip(rho, tops, tops[1:] + [DEPTH_MAX_KM * 1_000.0 * 2]):
        xs += [r, r]
        ys += [top / 1_000.0, bot / 1_000.0]
    return np.array(xs), np.array(ys)


def main():
    plt.rcParams.update(STYLE)
    fig, (ax_m, ax_r) = plt.subplots(
        1,
        2,
        figsize=(6.4, 3.0),
        dpi=200,
        gridspec_kw={"width_ratios": [1.0, 1.55]},
    )
    fig.patch.set_facecolor("white")

    for model, color in MODELS:
        xs, ys = step_profile(model)
        ax_m.plot(xs, ys, "-", lw=2.4, color=color, zorder=3)

        resp = MT1DForward(FREQS).run(model)
        ax_r.loglog(
            resp.freqs,
            resp.rho_a,
            "-",
            lw=2.4,
            color=color,
            label=model.name,
            zorder=3,
        )

    ax_m.set_xscale("log")
    ax_m.set_xlim(1.0, 3e4)
    ax_m.set_ylim(DEPTH_MAX_KM, 0)
    ax_m.set_title("layered-earth models", pad=6)
    ax_m.set_xlabel(r"$\rho$  ($\Omega\cdot$m)", fontsize=8.5)
    ax_m.set_ylabel("depth (km)", fontsize=8.5)

    ax_r.set_title("responses — MT1DForward", pad=6)
    ax_r.set_xlabel("frequency (Hz)", fontsize=8.5)
    ax_r.set_ylabel(r"$\rho_a$  ($\Omega\cdot$m)", fontsize=8.5)
    ax_r.invert_xaxis()
    ax_r.legend(fontsize=7.5, loc="lower left", frameon=False)

    for ax in (ax_m, ax_r):
        ax.grid(True, which="both", ls=":", lw=0.4, color="#c8d0dc", alpha=0.7)
        ax.tick_params(labelsize=7.5, length=2.5)

    fig.tight_layout(pad=0.6, w_pad=1.4)
    out = ROOT / "docs/source/_static/images/home/card-forward.png"
    fig.savefig(out, dpi=200, facecolor="white", bbox_inches="tight")
    print(f"saved: {out} ({out.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
