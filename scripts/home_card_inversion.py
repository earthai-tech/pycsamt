"""Home-page card image: "Inversion — classical & AI".

Left panel: a real ModEM 3-D result (bundled Willy line-02 sample) sliced
into the central N-S resistivity section, exactly like the inversion
gallery does with :class:`pycsamt.models.modem.InversionResult`.

Right panel: the AI side — :class:`pycsamt.ai.inversion.inv1d.EMInverter1D`
trained on a :func:`pycsamt.forward.batch.generate_dataset` synthetic set,
shown as predicted-vs-true layer resistivities on the held-out test split.

Output: docs/source/_static/images/home/card-inversion.png

Usage (any cwd):
    python scripts/home_card_inversion.py
"""

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

os.environ.setdefault("PYCSAMT_DOCS_BUILD", "1")

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(ROOT))

from pycsamt.ai.inversion.inv1d import EMInverter1D  # noqa: E402
from pycsamt.forward.batch import generate_dataset  # noqa: E402
from pycsamt.models.modem import InversionResult  # noqa: E402

BLUE, ORANGE, GOLD, SLATE = "#3e65b0", "#f15a29", "#d99114", "#5c677d"
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

SAMPLE = ROOT / "data" / "modem" / "willy_27freq_watex_line02_sample"
DEPTH_MAX_KM = 1.6
DIST_MAX_KM = 2.0


def modem_section():
    """Most-structured N-S section (distance, depth, rho) from the sample."""
    result = InversionResult(SAMPLE, load_control=False, load_covariance=False)
    key = "final" if "final" in result.models else sorted(result.models)[-1]
    model = result.models[key]

    # pick the slice the inversion actually updated: highest log-rho spread
    log_rho = np.log10(model.rho_linear)
    spread = np.array(
        [np.nanstd(log_rho[:, j, :]) for j in range(log_rho.shape[1])]
    )
    y_index = int(np.argmax(spread))
    rho = model.rho_linear[:, y_index, :]

    x_nodes = np.asarray(model.x_nodes, dtype=float)
    dist_km = (x_nodes - x_nodes[-1] / 2.0) / 1_000.0
    depth_km = np.asarray(model.z_nodes, dtype=float) / 1_000.0

    nz = int(np.searchsorted(depth_km, DEPTH_MAX_KM))
    keep = np.abs(dist_km) <= DIST_MAX_KM
    cols = np.where(keep[:-1] & keep[1:])[0]
    return (
        dist_km[cols[0] : cols[-1] + 2],
        depth_km[: nz + 1],
        rho[:nz, cols[0] : cols[-1] + 1],
        key,
    )


def main():
    plt.rcParams.update(STYLE)
    fig, (ax_l, ax_r) = plt.subplots(
        1,
        2,
        figsize=(6.4, 3.0),
        dpi=200,
        gridspec_kw={"width_ratios": [1.5, 1.0]},
    )
    fig.patch.set_facecolor("white")

    # -- classical: ModEM section ------------------------------------------
    dist, depth, rho, key = modem_section()
    mesh = ax_l.pcolormesh(
        dist,
        depth,
        rho,
        norm=mcolors.LogNorm(vmin=1.0, vmax=3_000.0),
        cmap="turbo_r",
        shading="flat",
        rasterized=True,
    )
    ax_l.set_ylim(DEPTH_MAX_KM, 0)
    ax_l.set_title("classical — ModEM section", pad=6)
    ax_l.set_xlabel("distance (km)", fontsize=8.5)
    ax_l.set_ylabel("depth (km)", fontsize=8.5)
    cbar = fig.colorbar(mesh, ax=ax_l, fraction=0.055, pad=0.03)
    cbar.set_label(r"$\rho$ ($\Omega\cdot$m)", fontsize=7.5)
    cbar.ax.tick_params(labelsize=6.5, length=2)

    # -- AI: learned 1-D inversion accuracy --------------------------------
    freqs = np.logspace(-1, 3, 24)
    ds = generate_dataset(
        solver="mt1d",
        n_samples=1_200,
        freqs=freqs,
        n_layers=4,
        noise_level=0.05,
        seed=0,
        verbose=False,
    )
    train, _val, test = ds.split()
    inv = EMInverter1D(arch="cnn1d", n_layers=4, solver="mt1d")
    inv.fit(train, epochs=45, batch_size=64, verbose=False)
    pred = inv.predict(test.X)

    lims = (0.0, 4.0)
    ax_r.plot(lims, lims, "--", lw=1.2, color=INK, alpha=0.6, zorder=2)
    for k, color in enumerate((BLUE, ORANGE, GOLD, SLATE)):
        ax_r.plot(
            test.y[:, k],
            pred[:, k],
            "o",
            ms=3.2,
            mec="none",
            color=color,
            alpha=0.75,
            label=f"layer {k + 1}",
            zorder=3,
        )
    ax_r.set_xlim(lims)
    ax_r.set_ylim(lims)
    ax_r.set_title("AI — EMInverter1D", pad=6)
    ax_r.set_xlabel(r"true $\log_{10}\rho$", fontsize=8.5)
    ax_r.set_ylabel(r"predicted $\log_{10}\rho$", fontsize=8.5)
    ax_r.legend(
        fontsize=6.8, loc="upper left", frameon=False, handletextpad=0.2
    )
    ax_r.grid(True, ls=":", lw=0.4, color="#c8d0dc", alpha=0.7)
    ax_r.tick_params(labelsize=7.5, length=2.5)
    ax_r.set_aspect("equal")

    fig.tight_layout(pad=0.6, w_pad=1.2)
    out = ROOT / "docs/source/_static/images/home/card-inversion.png"
    fig.savefig(out, dpi=200, facecolor="white", bbox_inches="tight")
    print(f"saved: {out} ({out.stat().st_size / 1024:.0f} KB)  model: {key}")


if __name__ == "__main__":
    main()
