"""Shared helpers for the ``Data corrections`` example gallery.

Not an example itself (sphinx-gallery skips ``_``-prefixed files). Loads the
survey line every correction wave works on, extracts apparent-resistivity /
phase curves from a ``Sites`` object, and draws the before/after comparison
panels the waves reuse. Unique name to avoid the shared-``sys.modules``
collision that same-named gallery helpers would hit.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np

from pycsamt.emtools._core import (
    _iter_items,
    _name,
    ensure_sites,
)

_IJ = {"xx": (0, 0), "xy": (0, 1), "yx": (1, 0), "yy": (1, 1)}


def demo_line(line: str = "L18PLT"):
    """Load one WILLY_DATA line as ``Sites`` (L18PLT carries real error tensors,
    which the static-shift and QC steps need)."""
    root = os.environ.get("PYCSAMT_DOCS_REPO_ROOT", ".")
    return ensure_sites(
        str(Path(root) / "data" / "AMT" / "WILLY_DATA" / line)
    )


def curves(sites, quantity: str = "rho", component: str = "xy"):
    """Return ``{station: (period_s, values)}`` for ``rho`` or ``phase``."""
    i, j = _IJ[component]
    out = {}
    for k, ed in enumerate(_iter_items(sites)):
        period = 1.0 / np.asarray(ed.freq, dtype=float)
        arr = np.asarray(
            ed.rho if quantity == "rho" else ed.phase, dtype=float
        )
        out[_name(ed, k)] = (period, arr[:, i, j])
    return out


def plot_before_after(
    before,
    after,
    stations,
    *,
    quantity="rho",
    labels=("raw", "corrected"),
    colors=("#b0b7c3", "#3e65b0"),
    title="",
    figsize=None,
):
    """Overlay raw vs corrected curves for several stations, one panel each."""
    import matplotlib.pyplot as plt

    n = len(stations)
    fig, axes = plt.subplots(
        1,
        n,
        figsize=figsize or (4.0 * n, 4.2),
        sharey=True,
        constrained_layout=True,
    )
    if n == 1:
        axes = [axes]
    for ax, st in zip(axes, stations):
        pb, vb = before[st]
        pa, va = after[st]
        if quantity == "rho":
            ax.loglog(pb, vb, ".", ms=3, color=colors[0], label=labels[0])
            ax.loglog(pa, va, "-", lw=1.7, color=colors[1], label=labels[1])
        else:
            ax.semilogx(pb, vb, ".", ms=3, color=colors[0], label=labels[0])
            ax.semilogx(pa, va, "-", lw=1.7, color=colors[1], label=labels[1])
        ax.set_title(st, fontsize=9)
        ax.set_xlabel("period (s)")
        ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.6)
    ylab = (
        r"$\rho_a$  ($\Omega\cdot$m)" if quantity == "rho" else "phase (deg)"
    )
    axes[0].set_ylabel(ylab)
    axes[0].legend(fontsize=8, framealpha=0.85)
    if title:
        fig.suptitle(title, fontsize=12)
    return fig
