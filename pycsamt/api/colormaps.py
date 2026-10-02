# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Colormap catalogue shared by every colormap picker.

The desktop studios, the API-configuration pages, the web app and Map View
all offer colormaps; each used to keep its own short list, and ``jet_r`` --
the classic resistivity colormap (red = conductive, blue = resistive) --
was missing from most of them.  :data:`COLORMAPS` is the one list; the
pickers call :func:`colormap_choices`.

Examples
--------
>>> from pycsamt.api.colormaps import COLORMAPS, colormap_choices
>>> COLORMAPS[:3]
('jet_r', 'jet', 'RdYlBu_r')
>>> colormap_choices("Blues")[-1]    # a valid current value is kept
'Blues'
"""

from __future__ import annotations

__all__ = [
    "COLORMAPS",
    "PLOTLY_COLORMAPS",
    "colormap_choices",
    "colormap_label",
    "colormap_options",
    "is_colormap",
]

#: Matplotlib colormaps offered by pickers, most useful for EM first.
COLORMAPS: tuple[str, ...] = (
    # resistivity / conductivity sections
    "jet_r", "jet", "RdYlBu_r", "RdYlBu", "Spectral_r", "Spectral",
    "turbo", "turbo_r", "rainbow_r", "rainbow",
    # perceptually uniform
    "viridis", "viridis_r", "plasma", "plasma_r", "inferno", "magma",
    "cividis",
    # diverging (anomalies, residuals, skew)
    "RdBu_r", "RdBu", "seismic", "coolwarm", "bwr", "PiYG", "BrBG",
    # terrain / single hue
    "terrain", "gist_earth", "hot", "hot_r", "YlOrRd", "copper",
    "Greys", "Greys_r", "Blues",
)

#: Plotly colour-scale names (lower-case; ``_r`` reverses any of them).
PLOTLY_COLORMAPS: tuple[str, ...] = (
    "jet_r", "jet", "rdylbu_r", "rdylbu", "spectral_r", "spectral",
    "turbo", "turbo_r", "rainbow_r", "rainbow",
    "viridis", "viridis_r", "plasma", "plasma_r", "inferno", "magma",
    "cividis", "rdbu_r", "rdbu", "balance", "earth", "hot", "ylorrd",
    "thermal", "greys",
)

_LABELS = {
    "jet_r": "jet_r (classic resistivity)",
    "RdYlBu_r": "RdYlBu_r (pyCSAMT default)",
}


def is_colormap(name) -> bool:
    """``True`` when matplotlib knows the colormap *name*."""
    if not name or not isinstance(name, str):
        return False
    try:
        import matplotlib

        matplotlib.colormaps[name]
    except Exception:
        return False
    return True


def colormap_choices(current: str | None = None, *, first=(),
                     plotly: bool = False) -> list[str]:
    """The catalogue as a list, keeping *current* when it is not in it.

    *first* names go to the top (a view's own default); *plotly* gives the
    Plotly colour-scale names instead of matplotlib's.
    """
    base = PLOTLY_COLORMAPS if plotly else COLORMAPS
    out = list(dict.fromkeys([*first, *base]))
    if current and current not in out and (plotly or is_colormap(current)):
        out.append(current)
    return out


def colormap_label(name: str) -> str:
    """Display text for *name* (a short hint for the key colormaps)."""
    return _LABELS.get(name, name)


def colormap_options(current: str | None = None, *, first=(),
                     plotly: bool = False) -> list[dict]:
    """:func:`colormap_choices` as ``[{"label", "value"}]`` (dropdowns)."""
    return [{"label": colormap_label(c), "value": c}
            for c in colormap_choices(current, first=first, plotly=plotly)]
