# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Style helpers for pyCSAMT map renderers."""

from __future__ import annotations

from .config import MapTheme

PLOTLY_CMAP_REMAP: dict[str, str] = {
    "coolwarm": "balance",
    "seismic": "rdbu",
    "terrain": "earth",
    "YlOrRd": "ylorrd",
    "RdBu_r": "rdbu_r",
    "copper": "thermal",
    "gnuplot2": "turbid",
    "tab10": "plotly3",
    "tab20": "plasma",
}

THEME_COLORS: dict[MapTheme, dict[str, str]] = {
    "dark": {
        "paper": "#1e1e2e",
        "plot": "#181825",
        "text": "#cdd6f4",
        "grid": "#313244",
        "accent": "#89b4fa",
    },
    "light": {
        "paper": "#eff1f5",
        "plot": "#ffffff",
        "text": "#4c4f69",
        "grid": "#ccd0da",
        "accent": "#1e66f5",
    },
    "publication": {
        "paper": "#ffffff",
        "plot": "#ffffff",
        "text": "#111111",
        "grid": "#cccccc",
        "accent": "#1565c0",
    },
}


def to_plotly_cmap(
    cmap: str | None,
    fallback: str = "plasma",
) -> str:
    """Return a Plotly-compatible colorscale name."""
    if not cmap:
        return fallback
    return PLOTLY_CMAP_REMAP.get(cmap, cmap)


def theme_colors(theme: MapTheme = "light") -> dict[str, str]:
    """Return a copied color palette for *theme*."""
    colors = THEME_COLORS.get(theme, THEME_COLORS["light"])
    return dict(colors)


def is_dark_theme(theme: MapTheme = "light") -> bool:
    """Return True for dark map styling."""
    return theme == "dark"
