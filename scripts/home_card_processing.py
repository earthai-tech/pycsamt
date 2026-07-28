"""Home-page card image: "Processing & corrections".

Runs a real pyCSAMT noise-removal pass (:func:`pycsamt.emtools.smooth_logfreq`)
on the L18PLT demo line from the corrections gallery and renders a compact
raw-vs-denoised figure sized for the landing-page feature card.

Output: docs/source/_static/images/home/card-processing.png

Usage (any cwd):
    python scripts/home_card_processing.py
"""

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

os.environ.setdefault("PYCSAMT_DOCS_BUILD", "1")
os.environ.setdefault("PYCSAMT_DOCS_REPO_ROOT", str(ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "docs" / "examples" / "corrections"))

from _corr_data import curves, demo_line  # noqa: E402

from pycsamt.emtools import smooth_logfreq  # noqa: E402

BLUE, ORANGE, GOLD = "#3e65b0", "#f15a29", "#d99114"
GRAY = "#a8b0bd"
INK = "#1c2941"
FMIN = 2.5e-3  # keep the coherent band; drop the noisy dead-band tail

STYLE = {
    "font.size": 9,
    "axes.linewidth": 0.8,
    "axes.edgecolor": "#c8d0dc",
    "xtick.color": "#5c677d",
    "ytick.color": "#5c677d",
    "axes.labelcolor": INK,
    "axes.titlesize": 11,
    "axes.titleweight": "bold",
    "axes.titlecolor": INK,
}


def clip(f, r):
    f, r = np.asarray(f), np.asarray(r)
    m = f >= FMIN
    return f[m], r[m]


def main():
    S = demo_line("L18PLT")
    raw = curves(S, "rho")
    clean = curves(
        smooth_logfreq(S, win=7, kind="tri", recursive=False), "rho"
    )

    # three moderately ragged stations at distinct resistivity levels
    def roughness(s):
        _, r = clip(*raw[s])
        return np.nanstd(np.diff(np.log10(r)))

    levels = {s: np.log10(np.nanmedian(raw[s][1])) for s in raw}
    ranked = sorted(raw, key=roughness, reverse=True)[4:14]
    ranked.sort(key=lambda s: levels[s])
    pick = [ranked[0], ranked[len(ranked) // 2], ranked[-1]]

    plt.rcParams.update(STYLE)
    fig, ax = plt.subplots(figsize=(6.4, 3.0), dpi=200)
    fig.patch.set_facecolor("white")

    for s, color in zip(pick, (BLUE, ORANGE, GOLD)):
        f, r = clip(*raw[s])
        ax.loglog(f, r, "-", lw=1.0, color=GRAY, alpha=0.95, zorder=2)
        ax.loglog(f, r, ".", ms=3.2, color=GRAY, alpha=0.9, zorder=2)
        fc, rc = clip(*clean[s])
        ax.loglog(fc, rc, "-", lw=2.6, color=color, zorder=3)

    ax.loglog([], [], ".-", color=GRAY, lw=1.0, ms=4, label="raw station data")
    ax.loglog([], [], "-", color=BLUE, lw=2.6, label="denoised (emtools)")
    ax.legend(fontsize=8, loc="upper right", frameon=False)

    ax.set_title("Noise removal — emtools.smooth_logfreq", pad=7)
    ax.grid(True, which="both", ls=":", lw=0.4, color="#c8d0dc", alpha=0.7)
    ax.tick_params(labelsize=7.5, length=2.5)
    ax.set_xlabel("frequency (Hz)", fontsize=8.5)
    ax.set_ylabel(r"$\rho_a$  ($\Omega\cdot$m)", fontsize=8.5)
    ax.invert_xaxis()

    fig.tight_layout(pad=0.6)
    out = ROOT / "docs/source/_static/images/home/card-processing.png"
    fig.savefig(out, dpi=200, facecolor="white", bbox_inches="tight")
    print(
        f"saved: {out} ({out.stat().st_size / 1024:.0f} KB)  stations: {pick}"
    )


if __name__ == "__main__":
    main()
