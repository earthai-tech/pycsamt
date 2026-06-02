"""
pycsamt.api.style
=================

Central visual-style system for the entire pyCSAMT package.

All plot functions that respect the style system read from the
module-level singleton :data:`PYCSAMT_STYLE`.  Users configure it once
and every subsequent figure inherits those choices.

Quick start
-----------
**Apply a named preset**::

    from pycsamt.api.style import use_style
    use_style("publication")          # grayscale publication look
    use_style("dark")                 # dark-background figures

**Override specific attributes**::

    from pycsamt.api.style import configure_style
    configure_style(
        rose__compass_labels = "degrees",   # dotted path: section__attr
        multiline__mode      = "gradient",
        multiline__base_color = "red",
        mt__xy__color        = "#003f88",
    )

**Temporary changes with a context manager**::

    from pycsamt.api.style import PYCSAMT_STYLE
    with PYCSAMT_STYLE.context("publication"):
        fig = plot_phase_tensor_rose(sites)
    # reverts to previous state automatically

**Build a fully custom style and share it**::

    from pycsamt.api.style import PyCSAMTStyle, MultilineStyle
    s = PyCSAMTStyle()
    s.multiline = MultilineStyle(mode="gradient", base_color="teal")
    s.mt.xy.color = "#004e89"
    # pass to functions that accept style=
    plot_phase_tensor_rose(sites, style=s.rose)

Style sections
--------------
* :class:`RoseStyle` — rose diagram visuals (bars, rings, compass, mean line)
* :class:`MultilineStyle` — gradient or cycle coloring for multi-line plots
* :class:`MTComponentStyle` — consistent per-component colors (XY/YX/XX/YY/TE/TM)
"""
from __future__ import annotations

import copy
import colorsys
from contextlib import contextmanager
from dataclasses import dataclass, field, fields as dc_fields
from typing import Any, Dict, Generator, List, Optional, Sequence, Tuple

import numpy as np

# ── Import RoseStyle from emtools (single source of truth) ────────────────────
from pycsamt.emtools._rose_style import (
    RoseStyle,
    resolve_rose_style,
    _PRESETS as _ROSE_PRESETS,
)


# ═══════════════════════════════════════════════════════════════════════════════
# MultilineStyle
# ═══════════════════════════════════════════════════════════════════════════════

# Maps common color names → sequential matplotlib cmaps (dark-first variant)
# Map base_color → sequential cmap where HIGH fraction = DARK shade.
# We use non-reversed cmaps so linspace(dark, light) = dark-first.
_COLOR_TO_CMAP: Dict[str, str] = {
    "blue":    "Blues",    "b":      "Blues",
    "red":     "Reds",     "r":      "Reds",
    "green":   "Greens",   "g":      "Greens",
    "orange":  "Oranges",
    "purple":  "Purples",
    "grey":    "Greys",    "gray":   "Greys",
    "teal":    "GnBu",
    "brown":   "YlOrBr",
    "navy":    "Blues",
    "crimson": "Reds",
}


@dataclass
class MultilineStyle:
    """Gradient-first coloring for multi-line (multi-station / multi-profile) plots.

    When *mode* is ``"gradient"``, the *n* lines returned by :meth:`colors`
    are shades of the same hue arranged from dark (line 0) to light (line
    n-1) — or reversed.  This preserves visual coherence while still making
    individual lines distinguishable.

    Parameters
    ----------
    mode : {"gradient", "cycle", "matplotlib"}
        ``"gradient"`` — sequential shades from *base_color* / *cmap*.
        ``"cycle"`` — rotate through a fixed palette (default: ``tab10``).
        ``"matplotlib"`` — fall back to matplotlib's default ``"C0"``,
        ``"C1"``, … colour cycle.
    base_color : str
        Base hue for ``mode="gradient"``.  If *cmap* is ``None``, a
        sequential colormap is auto-selected from this colour name (e.g.
        ``"blue"`` → ``"Blues_r"``).  Any matplotlib colour spec works;
        unrecognised names fall back to *cmap* if set, else ``"Blues_r"``.
    cmap : str or None
        Explicit colormap override.  Takes precedence over *base_color* for
        cmap selection.
    dark : float
        Fraction (0–1) of the cmap for the darkest (first) line.
        Default ``0.85`` avoids pure black.
    light : float
        Fraction for the lightest (last) line.  Default ``0.25`` avoids
        near-white.
    reverse : bool
        If ``True``, light lines come first and dark lines last.
    lw : float
        Default line width for multi-line plots.
    alpha : float
        Default line alpha.
    cycle_palette : list of str
        Colour list used when *mode* is ``"cycle"``.
    """

    mode:           str            = "gradient"
    base_color:     str            = "blue"
    cmap:           Optional[str]  = None
    dark:           float          = 0.85
    light:          float          = 0.25
    reverse:        bool           = False
    lw:             float          = 1.5
    alpha:          float          = 0.85
    cycle_palette:  List[str]      = field(default_factory=lambda: [
        "#1f77b4", "#d62728", "#2ca02c", "#ff7f0e",
        "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
        "#bcbd22", "#17becf",
    ])

    # ── public API ────────────────────────────────────────────────────────────

    def colors(self, n: int) -> List[Any]:
        """Return *n* colours for a multi-line plot.

        Parameters
        ----------
        n : int
            Number of lines / profiles.

        Returns
        -------
        list
            RGBA tuples (for ``"gradient"``) or colour strings
            (for ``"cycle"`` / ``"matplotlib"``).

        Examples
        --------
        >>> ms = MultilineStyle(mode="gradient", base_color="red")
        >>> colors = ms.colors(5)   # 5 shades of red, dark → light
        """
        if n <= 0:
            return []
        if self.mode == "matplotlib":
            return [f"C{i}" for i in range(n)]
        if self.mode == "cycle":
            pal = self.cycle_palette
            return [pal[i % len(pal)] for i in range(n)]
        # gradient
        import matplotlib.pyplot as plt
        cm_name = self.cmap or _COLOR_TO_CMAP.get(
            self.base_color.lower(), "Blues_r"
        )
        try:
            cm = plt.get_cmap(cm_name)
        except ValueError:
            cm = plt.get_cmap("Blues_r")
        fracs = np.linspace(self.dark, self.light, max(n, 2))[:n]
        if self.reverse:
            fracs = fracs[::-1]
        return [cm(float(f)) for f in fracs]

    def line_kwargs(self, idx: int, n: int, **overrides: Any) -> Dict[str, Any]:
        """Return ``ax.plot`` keyword arguments for line *idx* of *n*.

        Parameters
        ----------
        idx : int
            Zero-based line index.
        n : int
            Total number of lines.
        **overrides
            Keyword arguments that win over the style defaults.

        Examples
        --------
        >>> kw = ms.line_kwargs(0, 5)
        >>> ax.plot(x, y, **kw)
        """
        c = self.colors(n)[idx]
        d: Dict[str, Any] = dict(color=c, lw=self.lw, alpha=self.alpha)
        d.update(overrides)
        return d

    def copy(self, **kw: Any) -> "MultilineStyle":
        """Return a modified copy."""
        new = copy.copy(self)
        for k, v in kw.items():
            setattr(new, k, v)
        return new


# ═══════════════════════════════════════════════════════════════════════════════
# MTComponentStyle
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class _MTComp:
    """Visual properties for a single MT impedance component or mode.

    Attributes
    ----------
    color : str
        Line / marker colour.
    ls : str
        Line style (``"-"``, ``"--"``, ``":"``).
    lw : float
        Line width.
    marker : str
        Marker style.
    ms : float
        Marker size.
    mfc : str
        Marker face colour.  ``"white"`` gives open (hollow) markers.
    mew : float
        Marker edge width.
    alpha : float
        Overall opacity.
    label : str
        Default legend label.  Empty string → no label.
    elinewidth : float
        Error-bar cap line width.
    capsize : float
        Error-bar cap size in points.
    """

    color:       str   = "#1f77b4"
    ls:          str   = "-"
    lw:          float = 1.5
    marker:      str   = "o"
    ms:          float = 4.0
    mfc:         str   = "white"
    mew:         float = 1.2
    alpha:       float = 0.90
    label:       str   = ""
    elinewidth:  float = 0.8
    capsize:     float = 2.5

    def plot_kwargs(self, **overrides: Any) -> Dict[str, Any]:
        """Keyword arguments for :func:`matplotlib.axes.Axes.plot`.

        Examples
        --------
        >>> ax.plot(period, rho_xy, **style.mt.xy.plot_kwargs())
        """
        d: Dict[str, Any] = dict(
            color=self.color, ls=self.ls, lw=self.lw,
            marker=self.marker, ms=self.ms,
            mfc=self.mfc, mew=self.mew, alpha=self.alpha,
        )
        if self.label:
            d["label"] = self.label
        d.update(overrides)
        return d

    def errorbar_kwargs(self, **overrides: Any) -> Dict[str, Any]:
        """Keyword arguments for :func:`matplotlib.axes.Axes.errorbar`.

        The error bars are drawn in a 50 %-transparent version of the
        component colour.

        Examples
        --------
        >>> ax.errorbar(period, rho_xy, rho_err, **style.mt.xy.errorbar_kwargs())
        """
        import matplotlib.colors as mc
        rgba = list(mc.to_rgba(self.color))
        rgba[3] = 0.50          # semi-transparent error bars
        d = self.plot_kwargs()
        d.update(
            ecolor=tuple(rgba),
            elinewidth=self.elinewidth,
            capsize=self.capsize,
        )
        d.update(overrides)
        return d

    def copy(self, **kw: Any) -> "_MTComp":
        """Return a modified copy."""
        new = copy.copy(self)
        for k, v in kw.items():
            setattr(new, k, v)
        return new

    def __repr__(self) -> str:  # noqa: D105
        return (f"_MTComp(color={self.color!r}, ls={self.ls!r}, "
                f"lw={self.lw}, marker={self.marker!r}, label={self.label!r})")


@dataclass
class MTComponentStyle:
    """Consistent colours and markers for MT impedance components and modes.

    Attributes
    ----------
    xy : _MTComp
        Off-diagonal XY component (TE-like).  Default: solid blue circles.
    yx : _MTComp
        Off-diagonal YX component (TM-like).  Default: solid red squares.
    xx : _MTComp
        Diagonal XX component.  Default: dashed green upward triangles.
    yy : _MTComp
        Diagonal YY component.  Default: dashed purple downward triangles.
    te : _MTComp
        TE mode (when TE/TM decomposition is used).  Default: blue circles.
    tm : _MTComp
        TM mode.  Default: red squares.
    det : _MTComp
        Determinant / rotationally invariant quantity.  Default: grey circles.

    Examples
    --------
    Spread kwargs into a matplotlib call::

        ax.errorbar(per, rho, drho, **style.mt.xy.errorbar_kwargs())

    Override for one call only::

        kw = style.mt.yx.plot_kwargs(lw=2.5, alpha=1.0)

    Change the XY colour package-wide::

        PYCSAMT_STYLE.mt.xy.color = "#003f88"
    """

    xy:  _MTComp = field(default_factory=lambda: _MTComp(
        color="#1f77b4", ls="-",  lw=1.5, marker="o", ms=4.0,
        mfc="white", label="XY",
    ))
    yx:  _MTComp = field(default_factory=lambda: _MTComp(
        color="#d62728", ls="-",  lw=1.5, marker="s", ms=4.0,
        mfc="white", label="YX",
    ))
    xx:  _MTComp = field(default_factory=lambda: _MTComp(
        color="#2ca02c", ls="--", lw=1.2, marker="^", ms=3.5,
        mfc="white", label="XX",
    ))
    yy:  _MTComp = field(default_factory=lambda: _MTComp(
        color="#9467bd", ls="--", lw=1.2, marker="v", ms=3.5,
        mfc="white", label="YY",
    ))
    te:  _MTComp = field(default_factory=lambda: _MTComp(
        color="#1f77b4", ls="-",  lw=1.5, marker="o", ms=4.0,
        mfc="white", label="TE",
    ))
    tm:  _MTComp = field(default_factory=lambda: _MTComp(
        color="#d62728", ls="-",  lw=1.5, marker="s", ms=4.0,
        mfc="white", label="TM",
    ))
    det: _MTComp = field(default_factory=lambda: _MTComp(
        color="#7f7f7f", ls="-",  lw=1.3, marker="D", ms=3.5,
        mfc="white", label="det",
    ))

    # ── convenience ───────────────────────────────────────────────────────────

    def component(self, key: str) -> _MTComp:
        """Return the :class:`_MTComp` for *key* (case-insensitive).

        Parameters
        ----------
        key : str
            One of ``"xy"``, ``"yx"``, ``"xx"``, ``"yy"``,
            ``"te"``, ``"tm"``, ``"det"``.

        Raises
        ------
        KeyError
            If *key* is not a recognised component name.
        """
        try:
            return getattr(self, key.lower())
        except AttributeError:
            valid = [f.name for f in dc_fields(self)]
            raise KeyError(
                f"{key!r} is not a recognised MT component. "
                f"Valid keys: {valid}"
            ) from None

    def copy(self) -> "MTComponentStyle":
        """Return a deep copy."""
        return copy.deepcopy(self)


# ═══════════════════════════════════════════════════════════════════════════════
# PyCSAMTStyle — main container
# ═══════════════════════════════════════════════════════════════════════════════

class PyCSAMTStyle:
    """Global visual-style container for pyCSAMT.

    Holds three style sub-objects:

    * :attr:`rose` — :class:`~pycsamt.emtools._rose_style.RoseStyle`
    * :attr:`multiline` — :class:`MultilineStyle`
    * :attr:`mt` — :class:`MTComponentStyle`

    The module-level singleton :data:`PYCSAMT_STYLE` is the recommended
    entry point; import it and mutate it in your script or notebook before
    creating figures.

    Parameters
    ----------
    preset : str or None
        Named preset to initialise from.  ``None`` → ``"pycsamt"``.

    Examples
    --------
    Inspect the current style::

        from pycsamt.api.style import PYCSAMT_STYLE
        print(PYCSAMT_STYLE)

    Apply a named preset::

        PYCSAMT_STYLE.use("publication")

    Configure individual attributes with dotted paths::

        PYCSAMT_STYLE.configure(
            rose__compass_labels="degrees",
            multiline__base_color="red",
            mt__xy__color="#003f88",
        )

    Temporary style change (context manager)::

        with PYCSAMT_STYLE.context("dark"):
            fig = plot_phase_tensor_rose(sites)
    """

    #: Named full-package presets registered at import time.
    _PRESETS: Dict[str, Dict[str, Any]] = {}

    def __init__(self, preset: Optional[str] = None) -> None:
        self._init_defaults()
        if preset is not None:
            self.use(preset)

    def _init_defaults(self) -> None:
        self.rose      = RoseStyle()
        self.multiline = MultilineStyle()
        self.mt        = MTComponentStyle()

    # ── named presets ─────────────────────────────────────────────────────────

    def use(self, preset: str) -> None:
        """Apply a named full-package style preset.

        Parameters
        ----------
        preset : str
            ``"pycsamt"``  — default: gradient bars, NESW compass, blue/red MT.
            ``"publication"`` — grayscale, degree compass, open markers, tight.
            ``"dark"``     — dark-background optimised colours.

        Raises
        ------
        ValueError
            If *preset* is not registered.
        """
        key = preset.lower().strip()
        if key not in self._PRESETS:
            avail = ", ".join(f"{k!r}" for k in self._PRESETS)
            raise ValueError(
                f"Unknown style preset {preset!r}. Available: {avail}"
            )
        spec = self._PRESETS[key]
        # rose
        if "rose" in spec:
            self.rose = resolve_rose_style(spec["rose"])
        # multiline
        if "multiline" in spec:
            self.multiline = MultilineStyle(**spec["multiline"])
        # mt
        if "mt" in spec:
            self._apply_mt(spec["mt"])

    def _apply_mt(self, mt_spec: Dict[str, Any]) -> None:
        for comp_name, attrs in mt_spec.items():
            comp = getattr(self.mt, comp_name, None)
            if comp is None:
                continue
            for k, v in attrs.items():
                setattr(comp, k, v)

    # ── configure with dotted paths ───────────────────────────────────────────

    def configure(self, **kw: Any) -> None:
        """Configure style attributes using double-underscore dotted paths.

        Parameters
        ----------
        **kw
            Each key uses ``__`` as a section separator, e.g.
            ``rose__compass_labels="degrees"``,
            ``mt__xy__color="#003f88"``,
            ``multiline__base_color="red"``.

        Examples
        --------
        >>> PYCSAMT_STYLE.configure(
        ...     rose__compass_labels="degrees",
        ...     rose__show_secondary=False,
        ...     multiline__mode="gradient",
        ...     multiline__base_color="red",
        ...     mt__xy__color="#003f88",
        ...     mt__yx__lw=2.0,
        ... )
        """
        for path, value in kw.items():
            parts = path.split("__")
            obj = self
            for part in parts[:-1]:
                obj = getattr(obj, part)
            setattr(obj, parts[-1], value)

    # ── context manager ───────────────────────────────────────────────────────

    @contextmanager
    def context(
        self,
        preset: Optional[str] = None,
        **kw: Any,
    ) -> Generator["PyCSAMTStyle", None, None]:
        """Temporarily apply a preset and/or overrides, then restore.

        Parameters
        ----------
        preset : str or None
            Named preset to activate inside the context.
        **kw
            Additional dotted-path overrides (same syntax as
            :meth:`configure`).

        Yields
        ------
        PyCSAMTStyle
            The (temporarily modified) style object.

        Examples
        --------
        >>> with PYCSAMT_STYLE.context("publication", mt__xy__lw=2.0):
        ...     fig = plot_phase_tensor_rose(sites)
        """
        saved_rose      = copy.deepcopy(self.rose)
        saved_multiline = copy.deepcopy(self.multiline)
        saved_mt        = copy.deepcopy(self.mt)
        try:
            if preset is not None:
                self.use(preset)
            if kw:
                self.configure(**kw)
            yield self
        finally:
            self.rose      = saved_rose
            self.multiline = saved_multiline
            self.mt        = saved_mt

    # ── reset ─────────────────────────────────────────────────────────────────

    def reset(self) -> None:
        """Reset all style sections to package defaults (``"pycsamt"`` preset)."""
        self._init_defaults()

    # ── display ───────────────────────────────────────────────────────────────

    def summary(self) -> str:
        """Return a human-readable summary of the current style."""
        lines = ["PyCSAMTStyle"]
        lines.append(f"  rose.bar_style        = {self.rose.bar_style!r}")
        lines.append(f"  rose.cmap             = {self.rose.cmap!r}")
        lines.append(f"  rose.compass_labels   = {self.rose.compass_labels!r}")
        lines.append(f"  rose.show_mean        = {self.rose.show_mean}")
        lines.append(f"  rose.show_secondary   = {self.rose.show_secondary}")
        lines.append(f"  multiline.mode        = {self.multiline.mode!r}")
        lines.append(f"  multiline.base_color  = {self.multiline.base_color!r}")
        lines.append(f"  multiline.dark/light  = {self.multiline.dark}/{self.multiline.light}")
        lines.append(f"  mt.xy  color={self.mt.xy.color!r}  marker={self.mt.xy.marker!r}")
        lines.append(f"  mt.yx  color={self.mt.yx.color!r}  marker={self.mt.yx.marker!r}")
        lines.append(f"  mt.te  color={self.mt.te.color!r}  ls={self.mt.te.ls!r}")
        lines.append(f"  mt.tm  color={self.mt.tm.color!r}  ls={self.mt.tm.ls!r}")
        return "\n".join(lines)

    def __repr__(self) -> str:  # noqa: D105
        return self.summary()


# ═══════════════════════════════════════════════════════════════════════════════
# Full-package presets
# ═══════════════════════════════════════════════════════════════════════════════

def _register_preset(name: str, spec: Dict[str, Any]) -> None:
    PyCSAMTStyle._PRESETS[name.lower()] = spec


# "pycsamt" — default package look
_register_preset("pycsamt", {
    "rose": "pycsamt",              # string → resolved via _ROSE_PRESETS
    "multiline": dict(
        mode="gradient", base_color="blue", cmap=None,
        dark=0.85, light=0.25, reverse=False, lw=1.5, alpha=0.85,
    ),
    "mt": {
        "xy":  dict(color="#1f77b4", ls="-",  lw=1.5, marker="o", ms=4.0,
                    mfc="white", label="XY"),
        "yx":  dict(color="#d62728", ls="-",  lw=1.5, marker="s", ms=4.0,
                    mfc="white", label="YX"),
        "xx":  dict(color="#2ca02c", ls="--", lw=1.2, marker="^", ms=3.5,
                    mfc="white", label="XX"),
        "yy":  dict(color="#9467bd", ls="--", lw=1.2, marker="v", ms=3.5,
                    mfc="white", label="YY"),
        "te":  dict(color="#1f77b4", ls="-",  lw=1.5, marker="o", ms=4.0,
                    mfc="white", label="TE"),
        "tm":  dict(color="#d62728", ls="-",  lw=1.5, marker="s", ms=4.0,
                    mfc="white", label="TM"),
        "det": dict(color="#7f7f7f", ls="-",  lw=1.3, marker="D", ms=3.5,
                    mfc="white", label="det"),
    },
})

# "publication" — grayscale / BW-friendly for journals
_register_preset("publication", {
    "rose": "publication",
    "multiline": dict(
        mode="gradient", base_color="grey", cmap="Greys",
        dark=0.80, light=0.20, reverse=False, lw=1.5, alpha=0.90,
    ),
    "mt": {
        "xy":  dict(color="black",  ls="-",  lw=1.5, marker="o", ms=4.0,
                    mfc="white", label="XY"),
        "yx":  dict(color="black",  ls="-",  lw=1.5, marker="s", ms=4.0,
                    mfc="black", label="YX"),
        "xx":  dict(color="#555555", ls="--", lw=1.2, marker="^", ms=3.5,
                    mfc="white", label="XX"),
        "yy":  dict(color="#555555", ls="--", lw=1.2, marker="v", ms=3.5,
                    mfc="#555555", label="YY"),
        "te":  dict(color="black",  ls="-",  lw=1.5, marker="o", ms=4.0,
                    mfc="white", label="TE"),
        "tm":  dict(color="black",  ls="--", lw=1.5, marker="s", ms=4.0,
                    mfc="black", label="TM"),
        "det": dict(color="#777777", ls=":",  lw=1.3, marker="D", ms=3.5,
                    mfc="white", label="det"),
    },
})

# "dark" — optimised for dark-background Jupyter / presentation slides
_register_preset("dark", {
    "rose": "pycsamt",              # keep the default rose preset
    "multiline": dict(
        mode="gradient", base_color="blue", cmap="cool",
        dark=0.95, light=0.45, reverse=False, lw=2.0, alpha=0.90,
    ),
    "mt": {
        "xy":  dict(color="#74b9ff", ls="-",  lw=2.0, marker="o", ms=4.5,
                    mfc="#74b9ff", label="XY"),
        "yx":  dict(color="#ff7675", ls="-",  lw=2.0, marker="s", ms=4.5,
                    mfc="#ff7675", label="YX"),
        "xx":  dict(color="#55efc4", ls="--", lw=1.6, marker="^", ms=4.0,
                    mfc="#55efc4", label="XX"),
        "yy":  dict(color="#a29bfe", ls="--", lw=1.6, marker="v", ms=4.0,
                    mfc="#a29bfe", label="YY"),
        "te":  dict(color="#74b9ff", ls="-",  lw=2.0, marker="o", ms=4.5,
                    mfc="#74b9ff", label="TE"),
        "tm":  dict(color="#ff7675", ls="-",  lw=2.0, marker="s", ms=4.5,
                    mfc="#ff7675", label="TM"),
        "det": dict(color="#dfe6e9", ls="-",  lw=1.5, marker="D", ms=4.0,
                    mfc="#dfe6e9", label="det"),
    },
})


# ═══════════════════════════════════════════════════════════════════════════════
# Module-level singleton + convenience functions
# ═══════════════════════════════════════════════════════════════════════════════

#: Package-level singleton.  Import and mutate this object to configure
#: all pycsamt figures globally.
PYCSAMT_STYLE: PyCSAMTStyle = PyCSAMTStyle()


def use_style(preset: str) -> None:
    """Apply a named full-package style preset to :data:`PYCSAMT_STYLE`.

    Parameters
    ----------
    preset : str
        One of ``"pycsamt"``, ``"publication"``, ``"dark"``.

    Examples
    --------
    >>> from pycsamt.api.style import use_style
    >>> use_style("publication")
    """
    PYCSAMT_STYLE.use(preset)


def reset_style() -> None:
    """Reset :data:`PYCSAMT_STYLE` to package defaults.

    Examples
    --------
    >>> from pycsamt.api.style import reset_style
    >>> reset_style()
    """
    PYCSAMT_STYLE.reset()


def configure_style(**kw: Any) -> None:
    """Configure :data:`PYCSAMT_STYLE` with dotted-path keyword arguments.

    Parameters
    ----------
    **kw
        Dotted-path attributes using ``__`` as separator, e.g.
        ``rose__compass_labels="degrees"``,
        ``multiline__base_color="red"``,
        ``mt__xy__color="#003f88"``.

    Examples
    --------
    >>> from pycsamt.api.style import configure_style
    >>> configure_style(
    ...     rose__compass_labels="degrees",
    ...     multiline__mode="gradient",
    ...     multiline__base_color="teal",
    ...     mt__xy__color="#003f88",
    ...     mt__yx__color="#9b0000",
    ... )
    """
    PYCSAMT_STYLE.configure(**kw)


__all__ = [
    # main objects
    "PyCSAMTStyle",
    "MultilineStyle",
    "MTComponentStyle",
    "_MTComp",
    "RoseStyle",
    # singleton + helpers
    "PYCSAMT_STYLE",
    "use_style",
    "reset_style",
    "configure_style",
    # re-export for convenience
    "resolve_rose_style",
]
