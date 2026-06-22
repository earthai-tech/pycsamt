# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Package-wide topography configuration for 2-D section displays.

Follows the same singleton + configure/reset/context pattern as the
rest of :mod:`pycsamt.api`.  Import the singleton directly::

    from pycsamt.topo import PYCSAMT_TOPO, configure_topo

Or via the public API interface::

    from pycsamt.api.topo import configure_topo, reset_topo

Examples
--------
Enable topography globally (applies to every 2-D plot until reset)::

    from pycsamt.topo import configure_topo
    configure_topo(enabled=True)

Temporarily override inside a ``with`` block::

    from pycsamt.topo import PYCSAMT_TOPO
    with PYCSAMT_TOPO.context(enabled=True, exaggeration=3.0):
        fig = section.plot()
"""

from __future__ import annotations

import copy
from contextlib import contextmanager
from dataclasses import dataclass, field, fields
from typing import Any, Generator

__all__ = [
    "TopoConfig",
    "PYCSAMT_TOPO",
    "Y_DEPTH_TYPES",
    "Y_FREQ_TYPES",
    "configure_topo",
    "reset_topo",
]

# Y-axis types where real elevation exists and topo rendering is meaningful.
Y_DEPTH_TYPES: frozenset[str] = frozenset(
    {"depth", "elevation", "skin_depth", "elev", "z"}
)

# Y-axis types that represent EM frequency or period — topo is not applicable.
Y_FREQ_TYPES: frozenset[str] = frozenset(
    {"period", "frequency", "freq", "per", "t", "f"}
)


@dataclass
class TopoConfig:
    """Package-wide topography rendering policy for 2-D section plots.

    Once configured, the settings apply to every 2-D section plot
    (pseudosections, inversion sections, interpretation panels) until
    :meth:`reset` is called.

    Attributes
    ----------
    enabled : bool
        Master switch.  ``False`` → stations appear at a flat z = 0
        datum (default).  ``True`` → terrain-following geometry active.
    source : {"sites", "file", "array"}
        Where to read elevation data from.
        ``"sites"`` — read ``.elev`` from each Site's EDI HEAD (default).
        ``"file"``  — load from ``elev_file`` (CSV / parquet).
        ``"array"`` — use the ``elev_array`` ndarray directly.
    elev_array : array-like or None
        External elevation values (m a.s.l.) when ``source="array"``.
        Must have one value per station in profile order.
    elev_file : str or None
        Path to a tabular file when ``source="file"``.  Expected columns:
        ``station, elevation``  (or lat / lon / elevation).
    interp_method : {"linear", "cubic", "nearest"}
        Method for interpolating station elevations to intermediate x
        positions on the section grid.
    exaggeration : float
        Vertical exaggeration factor applied to the elevation profile
        during rendering.  1.0 = true scale.
    fill_color : str
        Matplotlib color for the above-surface terrain fill polygon.
    fill_alpha : float
        Opacity of the terrain fill (0–1).
    line_color : str
        Color of the terrain surface polyline.
    line_width : float
        Line width of the terrain surface polyline.
    show_surface_line : bool
        Draw the surface polyline on top of the terrain fill.
    clip_below_surface : bool
        Mask model cells that lie above the terrain surface (set to NaN).
    station_pins_at_surface : bool
        Place station marker triangles at the real terrain elevation
        instead of at the flat z = 0 datum.
    show_topo_strip : bool
        For period-vs-station pseudosections: add a thin elevation
        profile strip above the main colour image.
    strip_height_ratio : float
        Height of the topo strip relative to the main axes (0–1).
    """

    # ── Master switch ─────────────────────────────────────────────────
    enabled: bool = False

    # ── Elevation source ──────────────────────────────────────────────
    source: str = "sites"
    elev_array: Any = field(default=None, repr=False)
    elev_file: str | None = None

    # ── Interpolation ─────────────────────────────────────────────────
    interp_method: str = "linear"
    exaggeration: float = 1.0

    # ── Terrain rendering style ───────────────────────────────────────
    fill_color: str = "#a89070"
    fill_alpha: float = 0.40
    line_color: str = "#6b4e2a"
    line_width: float = 1.2
    show_surface_line: bool = True

    # ── Data masking + marker behaviour ──────────────────────────────
    clip_below_surface: bool = True
    station_pins_at_surface: bool = True

    # ── Pseudosection topo strip ──────────────────────────────────────
    show_topo_strip: bool = True
    strip_height_ratio: float = 0.18

    # ── Station marker placement ───────────────────────────────────────
    # Fraction of the visible y-range added above the topo surface to
    # separate the station marker from the terrain line.
    marker_pad_fraction: float = 0.015

    # ── Public API ────────────────────────────────────────────────────

    def configure(self, **kw: Any) -> "TopoConfig":
        """Set one or more configuration attributes by keyword."""
        for k, v in kw.items():
            if not hasattr(self, k):
                valid = [f.name for f in fields(self)]
                raise AttributeError(
                    f"Unknown topo config key {k!r}.  Valid keys: {valid}"
                )
            setattr(self, k, v)
        return self

    @contextmanager
    def context(self, **kw: Any) -> Generator["TopoConfig", None, None]:
        """Temporarily override config values, then restore them."""
        snapshot = {k: getattr(self, k) for k in kw}
        try:
            self.configure(**kw)
            yield self
        finally:
            for k, v in snapshot.items():
                setattr(self, k, v)

    def reset(self) -> None:
        """Restore all attributes to package defaults."""
        _defaults = TopoConfig()
        for f in fields(self):
            setattr(self, f.name, getattr(_defaults, f.name))

    def clone(self) -> "TopoConfig":
        """Return a deep copy of this config."""
        return copy.deepcopy(self)

    def is_active_for(self, y_type: str) -> bool:
        """Return True when topo rendering should fire for a given y-axis type.

        Topo only makes physical sense when the y-axis represents a real
        spatial quantity (depth, elevation, skin depth).  Period and
        frequency pseudosections carry no elevation information — topo
        rendering is silently skipped for those.

        Parameters
        ----------
        y_type : str
            String tag describing what the y-axis represents.  Recognised
            depth-like values: ``"depth"``, ``"elevation"``, ``"elev"``,
            ``"skin_depth"``, ``"z"``.  Recognised freq-like values:
            ``"period"``, ``"frequency"``, ``"freq"``, ``"per"``,
            ``"t"``, ``"f"``.

        Returns
        -------
        bool
        """
        if not self.enabled:
            return False
        return y_type.lower() in Y_DEPTH_TYPES

    def summary(self) -> str:
        """One-line status string."""
        src = self.source if self.enabled else "disabled"
        return f"TopoConfig(enabled={self.enabled}, source={src!r}, exag={self.exaggeration})"

    def __repr__(self) -> str:
        lines = ["TopoConfig("]
        for f in fields(self):
            lines.append(f"  {f.name}={getattr(self, f.name)!r},")
        lines.append(")")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Singleton + convenience helpers
# ---------------------------------------------------------------------------

PYCSAMT_TOPO: TopoConfig = TopoConfig()


def configure_topo(**kw: Any) -> None:
    """Configure the global :data:`PYCSAMT_TOPO` singleton.

    Parameters
    ----------
    **kw
        Keyword arguments forwarded to :meth:`TopoConfig.configure`.

    Examples
    --------
    >>> from pycsamt.topo import configure_topo
    >>> configure_topo(enabled=True, exaggeration=2.0)
    """
    PYCSAMT_TOPO.configure(**kw)


def reset_topo() -> None:
    """Reset :data:`PYCSAMT_TOPO` to package defaults."""
    PYCSAMT_TOPO.reset()
