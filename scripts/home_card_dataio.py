"""Home-page card image: "Data I/O & quality control".

Loads the bundled kap03 MT line (26 real EDI files) with
:func:`pycsamt.api.read_edis`, audits it with
:func:`pycsamt.emtools.build_qc_table`, and renders the audit as the card
figure: per-station period coverage coloured by SNR, and the SNR ranking
with noisy stations flagged.

Output: docs/source/_static/images/home/card-dataio.png

Usage (any cwd):
    python scripts/home_card_dataio.py
"""

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

os.environ.setdefault("PYCSAMT_DOCS_BUILD", "1")

import matplotlib

matplotlib.use("Agg")
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(ROOT))

from pycsamt.api import read_edis  # noqa: E402
from pycsamt.emtools import build_qc_table  # noqa: E402

BLUE, ORANGE = "#3e65b0", "#f15a29"
INK = "#1c2941"
SNR_FLAG = 20.0  # stations below this median SNR get flagged

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


def main():
    sites = read_edis(str(ROOT / "data" / "MT" / "kap03lmt_edis"))
    qc = build_qc_table(sites)
    qc = qc.sort_values("station").reset_index(drop=True)

    plt.rcParams.update(STYLE)
    fig, (ax_l, ax_r) = plt.subplots(
        1,
        2,
        figsize=(6.4, 3.0),
        dpi=200,
        gridspec_kw={"width_ratios": [1.35, 1.0]},
    )
    fig.patch.set_facecolor("white")

    # -- audit: period coverage per station, coloured by SNR ----------------
    norm = mcolors.LogNorm(
        vmin=max(qc["snr_med"].min(), 1.0), vmax=qc["snr_med"].max()
    )
    cmap = matplotlib.colormaps["viridis"]
    for i, row in qc.iterrows():
        ax_l.hlines(
            i,
            row["pmin"],
            row["pmax"],
            colors=cmap(norm(row["snr_med"])),
            lw=3.4,
            zorder=3,
        )
    ax_l.set_xscale("log")
    ax_l.set_ylim(-1, len(qc))
    ax_l.set_yticks(np.arange(0, len(qc), 4))
    ax_l.set_yticklabels(qc["station"][::4], fontsize=6.8)
    ax_l.set_xlabel("period (s)", fontsize=8.5)
    ax_l.set_title(f"audit — coverage · {len(qc)} EDIs", pad=6)
    ax_l.grid(
        True,
        axis="x",
        which="both",
        ls=":",
        lw=0.4,
        color="#c8d0dc",
        alpha=0.7,
    )
    ax_l.tick_params(labelsize=7, length=2.5)
    cbar = fig.colorbar(
        cm.ScalarMappable(norm=norm, cmap=cmap),
        ax=ax_l,
        fraction=0.05,
        pad=0.03,
    )
    cbar.set_label("median SNR", fontsize=7.5)
    cbar.ax.tick_params(labelsize=6.5, length=2)

    # -- QC: SNR ranking with flagged stations -------------------------------
    ranked = qc.sort_values("snr_med", ascending=False).reset_index(drop=True)
    colors = [ORANGE if s < SNR_FLAG else BLUE for s in ranked["snr_med"]]
    ax_r.bar(
        np.arange(len(ranked)),
        ranked["snr_med"],
        color=colors,
        width=0.8,
        zorder=3,
    )
    ax_r.axhline(SNR_FLAG, ls="--", lw=1.2, color=INK, alpha=0.6, zorder=4)
    ax_r.text(
        len(ranked) - 0.5,
        SNR_FLAG * 1.15,
        f"flag < {SNR_FLAG:.0f}",
        ha="right",
        fontsize=7,
        color=INK,
        alpha=0.8,
    )
    n_flag = int((ranked["snr_med"] < SNR_FLAG).sum())
    ax_r.set_yscale("log")
    ax_r.set_xticks([])
    ax_r.set_xlabel("stations (ranked)", fontsize=8.5)
    ax_r.set_ylabel("median SNR", fontsize=8.5)
    ax_r.set_title(f"QC — {n_flag} stations flagged", pad=6)
    ax_r.grid(
        True,
        axis="y",
        which="both",
        ls=":",
        lw=0.4,
        color="#c8d0dc",
        alpha=0.7,
    )
    ax_r.tick_params(labelsize=7.5, length=2.5)

    fig.tight_layout(pad=0.6, w_pad=1.2)
    out = ROOT / "docs/source/_static/images/home/card-dataio.png"
    fig.savefig(out, dpi=200, facecolor="white", bbox_inches="tight")
    print(f"saved: {out} ({out.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
