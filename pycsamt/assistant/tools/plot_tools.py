# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.tools.plot_tools
==================================

Line-aware plotting tools — produce a named section figure for a survey
line (or path) without running a full workflow agent.

Each entry in :data:`PLOTTERS` wraps a real ``pycsamt.emtools`` function
that takes a ``Sites`` object, so the assistant can ask for a specific
figure ("plot the SNR section for L22PLT") and get a saved image back.
"""

from __future__ import annotations

import importlib
from pathlib import Path
from typing import Any

__all__ = ["PLOTTERS", "list_plots", "make_plot"]

# plot kind → (module, function). Every function accepts ``sites`` first.
PLOTTERS: dict[str, tuple[str, str]] = {
    "confidence": (
        "pycsamt.emtools.qc",
        "plot_frequency_confidence_psection",
    ),
    "apparent_depth": (
        "pycsamt.emtools.frequency",
        "plot_apparent_depth_psection",
    ),
    "snr": ("pycsamt.emtools.advanced", "plot_snr_section"),
    "composite": (
        "pycsamt.emtools.advanced",
        "plot_mt_composite_section",
    ),
    "strike": (
        "pycsamt.emtools.advanced",
        "plot_strike_stability_bands",
    ),
}


def list_plots() -> list[str]:
    """Available plot kinds."""
    return sorted(PLOTTERS)


def _to_figure(obj: Any) -> Any | None:
    """Coerce a plotter return (Figure / Axes / tuple) to a Figure."""
    from matplotlib.figure import Figure

    if isinstance(obj, Figure):
        return obj
    if hasattr(obj, "get_figure"):
        try:
            return obj.get_figure()
        except Exception:  # noqa: BLE001
            return None
    if isinstance(obj, (list, tuple)):
        for x in obj:
            fig = _to_figure(x)
            if fig is not None:
                return fig
    return None


def make_plot(
    kind: str,
    line_or_path: str,
    *,
    output_dir: str | None = None,
    registry: Any | None = None,
    fmt: str = "png",
    dpi: int = 150,
) -> dict[str, Any]:
    """Render *kind* for *line_or_path* and save it.

    Returns ``{"kind", "line", "path", "figure_path", "status"}``.
    Raises ``ValueError`` for an unknown kind and ``FileNotFoundError``
    when the resolved EDI directory is missing.
    """
    if kind not in PLOTTERS:
        raise ValueError(f"Unknown plot kind {kind!r}. Available: {list_plots()}")

    # headless: never open a window
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    from .workflow_tools import resolve_target

    target = resolve_target(line_or_path, registry=registry)
    edi_dir = target.get("edi_dir") or target.get("path")
    if not target.get("exists"):
        raise FileNotFoundError(f"EDI directory not found: {edi_dir!r}")

    out = Path(output_dir or target.get("output_root") or "pycsamt_plots")
    out.mkdir(parents=True, exist_ok=True)

    from pycsamt.emtools._core import ensure_sites

    sites = ensure_sites(edi_dir, recursive=True, verbose=0)

    mod_name, fn_name = PLOTTERS[kind]
    fn = getattr(importlib.import_module(mod_name), fn_name)
    result = fn(sites)
    fig = _to_figure(result)

    fig_path: str | None = None
    status = "success"
    if fig is not None:
        line = target.get("line") or "section"
        fig_path = str(out / f"{kind}_{line}.{fmt}")
        try:
            fig.savefig(fig_path, dpi=dpi, bbox_inches="tight")
        except Exception:  # noqa: BLE001
            status = "saved_failed"
        plt.close(fig)
    else:
        status = "no_figure"

    return {
        "kind": kind,
        "line": target.get("line"),
        "path": edi_dir,
        "figure_path": fig_path,
        "status": status,
    }
