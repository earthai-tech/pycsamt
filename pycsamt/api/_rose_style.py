"""
pycsamt.emtools._rose_style
===========================

Shared visual-style system for all rose diagram functions
(:func:`~pycsamt.emtools.strike.plot_strike_rose`,
:func:`~pycsamt.emtools.tensor.plot_phase_tensor_rose`, and any future
rose plot added to emtools).

Usage
-----
**Named preset** (recommended default)::

    ax = plot_phase_tensor_rose(sites, style="pycsamt")

**Override individual attributes** while keeping the rest of the preset::

    ax = plot_phase_tensor_rose(
        sites,
        style="pycsamt",
        compass_labels="degrees",  # switch NESW → degree values
        mean_ls="--",  # dashed mean line
        show_secondary=False,  # hide conjugate line
    )

**Build a custom style** and reuse it across multiple figures::

    from pycsamt.emtools import RoseStyle

    my_style = RoseStyle(
        bar_style="solid",
        bar_color="#2c7fb8",
        cmap="Blues",
        compass_labels="degrees",
        show_mean=True,
        mean_color="navy",
        show_secondary=False,
    )
    ax1 = plot_phase_tensor_rose(sites1, style=my_style)
    ax2 = plot_strike_rose(sites2, style=my_style)

**Named presets available**

- ``"pycsamt"`` / ``"pycsamt-rose"`` — default: gradient bars, bold outer
  ring, NESW compass, crimson mean line, dashed conjugate, annotation box.
- ``"minimal"`` — solid bars, sparse rings, degree compass, no mean lines,
  no annotation.
- ``"publication"`` — grayscale gradient, degree compass, solid/dotted mean
  lines, tight annotation.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass
from typing import Any

# ── Sentinel ────────────────────────────────────────────────────────────────────
# Used as the default for visual kwargs in rose-plot functions.
# Allows the caller to distinguish "I didn't pass this" (= use style value)
# from "I explicitly set this to X" (= override the style).
_UNSET = object()


# ── RoseStyle dataclass ─────────────────────────────────────────────────────────


@dataclass
class RoseStyle:
    """Visual style bundle shared by all pycsamt rose diagram functions.

    Every attribute maps 1-to-1 to a keyword argument accepted by
    :func:`~pycsamt.emtools.strike.plot_strike_rose` and
    :func:`~pycsamt.emtools.tensor.plot_phase_tensor_rose`.  Pass an
    instance via the ``style=`` parameter of either function; individual
    kwargs still override the style.

    Attributes
    ----------
    bar_style : {"gradient", "solid", "bands"}
        ``"gradient"`` — bar height mapped to *cmap*;
        ``"solid"`` — uniform *bar_color*;
        ``"bands"`` — one colour per period sub-band (stacked).
    bar_color : str
        Bar fill colour for ``bar_style="solid"``.
    bar_alpha : float
        Bar opacity (0–1).
    bar_edgecolor : str
        Bar edge colour.  ``"none"`` disables edges.
    bar_edgelw : float
        Bar edge line-width.
    cmap : str
        Matplotlib colormap name for ``bar_style="gradient"``.
    outer_ring_lw : float
        Line-width of the bold outer bounding circle.
    outer_ring_color : str
        Colour of the outer bounding circle.
    n_rings : int
        Number of concentric reference rings inside the plot.
    ring_color : str
        Colour of concentric rings.
    ring_ls : str
        Line-style of concentric rings (``":"`` dotted, ``"--"`` dashed, …).
    ring_lw : float
        Line-width of concentric rings.
    ring_labels : list[float] or None
        Explicit count values to annotate on the rings
        (e.g. ``[25, 50, 75, 100]``).  ``None`` → evenly auto-spaced.
    ring_label_angle : float
        Clockwise degrees from North where count labels appear.
    ring_label_fontsize : float
        Font size for ring count annotations.
    ring_label_color : str
        Colour for ring count annotations.
    ring_label_fmt : str
        Python format string for ring labels (e.g. ``"{:.0f}"``).
    spoke_every : float
        Angular spacing (degrees) between radial spokes / tick marks.
    spoke_color : str
        Colour of radial spokes.
    spoke_ls : str
        Line-style of radial spokes.
    spoke_lw : float
        Line-width of radial spokes.
    compass_labels : {"NESW", "degrees", "none"}
        Perimeter labels.  ``"NESW"`` → N/E/S/W;
        ``"degrees"`` → 0°/45°/…; ``"none"`` → hidden.
    compass_fontsize : float
        Font size for compass / degree perimeter labels.
    compass_color : str
        Colour for perimeter labels.
    compass_fontweight : str
        Font weight for perimeter labels (``"bold"``, ``"normal"``, …).
    show_mean : bool
        Draw the mean direction as a diameter line through the centre.
    mean_color : str
        Colour of the mean direction line.
    mean_lw : float
        Line-width of the mean direction line.
    mean_ls : str
        Line-style of the mean direction line.
    show_secondary : bool
        Draw the 180°-conjugate (axial-symmetry) mean line.
    secondary_color : str or None
        Colour of the conjugate line.  ``None`` → same as *mean_color*.
    secondary_ls : str
        Line-style of the conjugate line.
    secondary_lw : float or None
        Line-width of the conjugate line.  ``None`` → same as *mean_lw*.
    show_annotation : bool
        Show the text box with mean angle and count.
    annotation_pos : (float, float)
        Axes-fraction ``(x, y)`` position of the annotation box.
    annotation_fontsize : float
        Font size of the annotation text.
    annotation_bg : str
        Background colour of the annotation box.
    annotation_ec : str
        Edge colour of the annotation box.
    show_n : bool
        Append ``n = N`` (data-point count) to the annotation text.
    """

    # ── bars ──────────────────────────────────────────────────────────────────
    bar_style: str = "gradient"
    bar_color: str = "#e53935"
    bar_alpha: float = 0.88
    bar_edgecolor: str = "none"
    bar_edgelw: float = 0.4
    cmap: str = "YlOrRd"

    # ── outer ring ────────────────────────────────────────────────────────────
    outer_ring_lw: float = 2.5
    outer_ring_color: str = "0.12"

    # ── concentric rings ──────────────────────────────────────────────────────
    n_rings: int = 3
    ring_color: str = "0.75"
    ring_ls: str = ":"
    ring_lw: float = 0.7
    ring_labels: list[float] | None = None  # None → auto-spaced
    ring_label_angle: float = 22.5
    ring_label_fontsize: float = 7.0
    ring_label_color: str = "0.30"
    ring_label_fmt: str = "{:.0f}"

    # ── spokes ────────────────────────────────────────────────────────────────
    spoke_every: float = 45.0
    spoke_color: str = "0.72"
    spoke_ls: str = ":"
    spoke_lw: float = 0.7

    # ── compass ───────────────────────────────────────────────────────────────
    compass_labels: str = "NESW"  # "NESW" | "degrees" | "none"
    compass_fontsize: float = 8.5
    compass_color: str = "0.15"
    compass_fontweight: str = "bold"

    # ── mean direction ────────────────────────────────────────────────────────
    show_mean: bool = True
    mean_color: str = "crimson"
    mean_lw: float = 2.2
    mean_ls: str = "-"
    show_secondary: bool = True
    secondary_color: str | None = None  # None → same as mean_color
    secondary_ls: str = "--"
    secondary_lw: float | None = None  # None → same as mean_lw

    # ── annotation box ────────────────────────────────────────────────────────
    show_annotation: bool = True
    annotation_pos: tuple = (0.05, 0.93)
    annotation_fontsize: float = 8.0
    annotation_bg: str = "white"
    annotation_ec: str = "0.25"
    show_n: bool = True

    # ── helpers ───────────────────────────────────────────────────────────────

    def copy(self, **overrides: Any) -> RoseStyle:
        """Return a shallow copy with *overrides* applied.

        Parameters
        ----------
        **overrides
            Any :class:`RoseStyle` attribute name and its new value.

        Returns
        -------
        RoseStyle
            A new :class:`RoseStyle` instance.

        Raises
        ------
        ValueError
            If an unknown attribute name is passed.

        Examples
        --------
        >>> rs = RoseStyle()
        >>> rs2 = rs.copy(compass_labels="degrees", show_secondary=False)
        """
        new = copy.copy(self)
        for k, v in overrides.items():
            if not hasattr(new, k):
                raise ValueError(
                    f"RoseStyle has no attribute {k!r}. "
                    f"Valid attributes: {sorted(vars(new).keys())}"
                )
            setattr(new, k, v)
        return new

    def __repr__(self) -> str:  # noqa: D105
        pairs = ", ".join(
            f"{k}={v!r}"
            for k, v in vars(self).items()
            if not k.startswith("_")
        )
        return f"RoseStyle({pairs})"


# ── Named presets ───────────────────────────────────────────────────────────────

_PRESETS: dict[str, RoseStyle] = {}


def _reg(name: str, **kw: Any) -> None:
    _PRESETS[name] = RoseStyle(**kw)


# "pycsamt" / "pycsamt-rose" — default look matching the paper style
_reg(
    "pycsamt",
    bar_style="gradient",
    cmap="YlOrRd",
    bar_alpha=0.88,
    bar_edgecolor="none",
    outer_ring_lw=2.5,
    outer_ring_color="0.12",
    n_rings=3,
    ring_color="0.75",
    ring_ls=":",
    ring_lw=0.7,
    ring_labels=None,
    ring_label_angle=22.5,
    ring_label_fontsize=7.0,
    ring_label_color="0.30",
    ring_label_fmt="{:.0f}",
    spoke_every=45.0,
    spoke_color="0.72",
    spoke_ls=":",
    spoke_lw=0.7,
    compass_labels="NESW",
    compass_fontsize=8.5,
    compass_color="0.15",
    compass_fontweight="bold",
    show_mean=True,
    mean_color="crimson",
    mean_lw=2.2,
    mean_ls="-",
    show_secondary=True,
    secondary_ls="--",
    show_annotation=True,
    annotation_fontsize=8.0,
    annotation_bg="white",
    annotation_ec="0.25",
    show_n=True,
)
_PRESETS["pycsamt-rose"] = _PRESETS["pycsamt"]  # alias

# "minimal" — solid bars, sparse decorations, degree compass, no mean
_reg(
    "minimal",
    bar_style="solid",
    bar_color="#4e79a7",
    bar_alpha=0.80,
    outer_ring_lw=1.2,
    outer_ring_color="0.30",
    n_rings=2,
    ring_color="0.72",
    ring_ls=":",
    ring_lw=0.5,
    spoke_every=90.0,
    spoke_color="0.65",
    spoke_ls=":",
    spoke_lw=0.5,
    compass_labels="degrees",
    compass_fontsize=7.5,
    compass_fontweight="normal",
    show_mean=False,
    show_secondary=False,
    show_annotation=False,
    show_n=False,
)

# "publication" — grayscale gradient, degree compass, solid+dotted mean
_reg(
    "publication",
    bar_style="gradient",
    cmap="Greys",
    bar_alpha=0.85,
    bar_edgecolor="0.4",
    bar_edgelw=0.3,
    outer_ring_lw=2.0,
    outer_ring_color="0.0",
    n_rings=4,
    ring_color="0.60",
    ring_ls=":",
    ring_lw=0.6,
    ring_label_angle=15.0,
    ring_label_fontsize=6.5,
    spoke_every=45.0,
    spoke_color="0.60",
    spoke_ls=":",
    spoke_lw=0.6,
    compass_labels="degrees",
    compass_fontsize=7.5,
    compass_color="0.0",
    compass_fontweight="bold",
    show_mean=True,
    mean_color="black",
    mean_lw=1.8,
    mean_ls="-",
    show_secondary=True,
    secondary_ls=":",
    secondary_lw=1.2,
    show_annotation=True,
    annotation_fontsize=7.5,
    annotation_bg="white",
    annotation_ec="0.0",
    show_n=True,
)


# ── Public resolver ─────────────────────────────────────────────────────────────


def resolve_rose_style(
    style: str | RoseStyle | None,
    **overrides: Any,
) -> RoseStyle:
    """Resolve *style* to a :class:`RoseStyle`, then apply *overrides*.

    Parameters
    ----------
    style : str, RoseStyle, or None
        Named preset string, a :class:`RoseStyle` instance, or ``None``
        (falls back to ``"pycsamt"``).
    **overrides
        Any :class:`RoseStyle` attribute to override after resolving the
        base preset, e.g. ``compass_labels="degrees"``.

    Returns
    -------
    RoseStyle

    Raises
    ------
    ValueError
        If *style* is a string that does not match a known preset.
    TypeError
        If *style* is not a str, RoseStyle, or None.

    Examples
    --------
    >>> rs = resolve_rose_style("pycsamt", compass_labels="degrees")
    >>> rs.compass_labels
    'degrees'
    """
    if style is None:
        base = _PRESETS["pycsamt"]
    elif isinstance(style, str):
        key = style.lower().replace("_", "-")
        if key not in _PRESETS:
            avail = ", ".join(f"{k!r}" for k in _PRESETS)
            raise ValueError(
                f"Unknown rose style {style!r}. Available presets: {avail}"
            )
        base = _PRESETS[key]
    elif isinstance(style, RoseStyle):
        base = style
    else:
        raise TypeError(
            f"style must be a str, RoseStyle, or None; got {type(style).__name__!r}"
        )
    return base.copy(**overrides) if overrides else base


__all__ = ["RoseStyle", "resolve_rose_style", "_UNSET"]
