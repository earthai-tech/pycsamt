"""
pycsamt.api.plot
================

Global plot-export configuration for pyCSAMT.

Every pyCSAMT plot function returns a ``Figure`` or ``Axes``.  Saving that
object through :func:`save_fig` (or :meth:`PlotConfig.save`) automatically
honours whatever format, DPI, and output directory you have configured
once — instead of having to repeat the same ``savefig`` kwargs in every
script.

Quick start
-----------
**Single format** (default is PNG)::

    from pycsamt.api.plot import save_fig, set_fmt, set_dpi

    ax = plot_phase_tensor_psection(sites)
    save_fig(ax, "output/fig_pt")          # → fig_pt.png

**Switch the global format**::

    set_fmt("svg")
    save_fig(ax, "output/fig_pt")          # → fig_pt.svg

**Additive** — keep the default PNG *plus* add SVG::

    set_fmt("+svg")
    save_fig(ax, "output/fig_pt")          # → fig_pt.png  +  fig_pt.svg

**Multiple explicit formats** (no PNG unless listed)::

    set_fmt("svg", "pdf")
    save_fig(ax, "output/fig_pt")          # → fig_pt.svg  +  fig_pt.pdf

**Additive stack** — PNG + SVG + PDF::

    set_fmt("+svg", "+pdf")
    # equivalent:
    PLOT_CONFIG.fmt = ["+svg", "+pdf"]

**Per-call override** (does not change the global config)::

    save_fig(ax, "output/fig_pt", fmt="eps")

**Context manager** (reverts on exit)::

    from pycsamt.api.plot import PLOT_CONFIG
    with PLOT_CONFIG.context(fmt=["+svg"], dpi=300):
        save_fig(fig, "publication/fig_pt")   # → .png + .svg @ 300 DPI
    # back to previous config

**Set a global output directory**::

    set_savedir("~/figures")
    save_fig(ax, "fig_pt")     # → ~/figures/fig_pt.png

**Global DPI**::

    set_dpi(300)

External control
----------------
Any attribute can be set **without modifying Python code** via:

1. **Environment variables** (highest precedence)::

    PYCSAMT_FMT="+svg,+pdf"  PYCSAMT_DPI=300  python make_figures.py

   ``PYCSAMT_FMT`` is a comma-separated list of format tokens (``+`` prefix
   supported).  Other variables: ``PYCSAMT_DPI``, ``PYCSAMT_SAVEDIR``,
   ``PYCSAMT_TRANSPARENT``, ``PYCSAMT_BBOX``.

2. **Config file** (``pycsamt_plot.cfg`` in the working directory, or
   ``~/.pycsamt/plot.cfg`` as a user-global fallback)::

    [plot]
    fmt        = +svg, +pdf
    dpi        = 300
    savedir    = ~/my_figures
    transparent = false
    bbox_inches = tight

   Local file takes precedence over the home-dir file.  Both are read at
   import time and overridden by environment variables.

3. **Python API** (programmatic, in scripts or notebooks)::

    from pycsamt.api.plot import PLOT_CONFIG, set_fmt, set_dpi
    PLOT_CONFIG.dpi = 300
    set_fmt("+svg")

Format token rules
------------------
Each token is a format extension (without the leading dot):

* ``"png"``    → save as PNG only
* ``"svg"``    → save as SVG only
* ``"+svg"``   → save as *base* (default ``"png"``) **and** SVG
* ``"+pdf"``   → save as base **and** PDF
* A comma-separated string is split: ``"png,svg"`` → ``["png", "svg"]``
* The base format can be changed via ``PLOT_CONFIG.base_fmt = "tiff"``

Supported formats: ``png``, ``svg``, ``pdf``, ``eps``, ``tiff``,
``jpg``/``jpeg``.

PlotConfig attributes
---------------------
.. autosummary::

   PlotConfig.fmt
   PlotConfig.base_fmt
   PlotConfig.dpi
   PlotConfig.bbox_inches
   PlotConfig.transparent
   PlotConfig.facecolor
   PlotConfig.savedir
   PlotConfig.close_after_save
   PlotConfig.verbose
"""
from __future__ import annotations

import configparser
import copy
import os
from collections.abc import Generator, Sequence
from contextlib import contextmanager
from pathlib import Path
from typing import (
    Any,
)

# ── supported format extensions ───────────────────────────────────────────────
_VALID_FMTS = frozenset({"png", "svg", "pdf", "eps", "tiff", "jpg", "jpeg"})

# ── where to look for the config file ─────────────────────────────────────────
_CFG_LOCAL  = Path("pycsamt_plot.cfg")
_CFG_GLOBAL = Path.home() / ".pycsamt" / "plot.cfg"


# ═══════════════════════════════════════════════════════════════════════════════
# Internal helpers
# ═══════════════════════════════════════════════════════════════════════════════

def _normalise_token(tok: str) -> str:
    """Strip whitespace and leading dot from a single format token."""
    return tok.strip().lstrip(".")


def _parse_fmt_tokens(raw: str | Sequence) -> list[str]:
    """Return a list of format tokens (preserving ``+`` prefix) from *raw*.

    Accepts a string (comma-separated), a single token, or an iterable.
    """
    if isinstance(raw, str):
        parts = [p.strip() for p in raw.split(",") if p.strip()]
    else:
        # flatten: each element may itself be a comma-separated string
        parts = []
        for item in raw:
            parts.extend(p.strip() for p in str(item).split(",") if p.strip())
    return [_normalise_token(p) for p in parts]


def _resolve_formats(tokens: list[str], base_fmt: str) -> list[str]:
    """Expand additive tokens and return a deduplicated ordered format list.

    Parameters
    ----------
    tokens : list[str]
        Raw tokens, possibly prefixed with ``+``.
    base_fmt : str
        The base format included when any additive token is present.

    Returns
    -------
    list[str]
        Clean, deduplicated list of format strings (no ``+``, no dot).

    Examples
    --------
    >>> _resolve_formats(["png"], "png")
    ['png']
    >>> _resolve_formats(["+svg"], "png")
    ['png', 'svg']
    >>> _resolve_formats(["+svg", "+pdf"], "png")
    ['png', 'svg', 'pdf']
    >>> _resolve_formats(["svg", "pdf"], "png")
    ['svg', 'pdf']
    >>> _resolve_formats(["png", "+svg"], "png")
    ['png', 'svg']
    """
    base_fmt = _normalise_token(base_fmt)
    result: list[str] = []
    has_additive = any(t.startswith("+") for t in tokens)
    if has_additive:
        result.append(base_fmt)
    for t in tokens:
        if t.startswith("+"):
            fmt = t.lstrip("+")
        else:
            fmt = t
        if fmt and fmt not in result:
            result.append(fmt)
    return result if result else [base_fmt]


def _validate_fmt(fmts: list[str]) -> None:
    """Warn about unrecognised format names (does not raise)."""
    for f in fmts:
        if f.lower() not in _VALID_FMTS:
            import warnings
            warnings.warn(
                f"pycsamt PlotConfig: format {f!r} is not in the known list "
                f"{sorted(_VALID_FMTS)}.  matplotlib may still accept it.",
                stacklevel=4,
            )


def _load_config_file() -> dict:
    """Read ``pycsamt_plot.cfg`` (CWD) or ``~/.pycsamt/plot.cfg``."""
    cp = configparser.ConfigParser()
    for candidate in (_CFG_LOCAL, _CFG_GLOBAL):
        if candidate.exists():
            cp.read(candidate)
            break
    if "plot" not in cp:
        return {}
    s = cp["plot"]
    out: dict = {}
    if "fmt" in s:
        out["fmt"] = [t.strip() for t in s["fmt"].split(",") if t.strip()]
    if "base_fmt" in s:
        out["base_fmt"] = s["base_fmt"].strip()
    if "dpi" in s:
        out["dpi"] = int(s["dpi"])
    if "bbox_inches" in s:
        out["bbox_inches"] = s["bbox_inches"].strip()
    if "transparent" in s:
        out["transparent"] = s.getboolean("transparent")
    if "facecolor" in s:
        out["facecolor"] = s["facecolor"].strip()
    if "savedir" in s:
        out["savedir"] = os.path.expanduser(s["savedir"].strip())
    if "close_after_save" in s:
        out["close_after_save"] = s.getboolean("close_after_save")
    if "verbose" in s:
        out["verbose"] = s.getboolean("verbose")
    return out


def _load_env() -> dict:
    """Read ``PYCSAMT_*`` environment variables."""
    out: dict = {}
    if v := os.environ.get("PYCSAMT_FMT"):
        out["fmt"] = [t.strip() for t in v.split(",") if t.strip()]
    if v := os.environ.get("PYCSAMT_BASE_FMT"):
        out["base_fmt"] = v.strip()
    if v := os.environ.get("PYCSAMT_DPI"):
        out["dpi"] = int(v)
    if v := os.environ.get("PYCSAMT_SAVEDIR"):
        out["savedir"] = os.path.expanduser(v)
    if v := os.environ.get("PYCSAMT_TRANSPARENT"):
        out["transparent"] = v.lower() in ("1", "true", "yes")
    if v := os.environ.get("PYCSAMT_BBOX"):
        out["bbox_inches"] = v.strip()
    return out


# ═══════════════════════════════════════════════════════════════════════════════
# PlotConfig
# ═══════════════════════════════════════════════════════════════════════════════

class PlotConfig:
    """Global plot-export configuration singleton for pyCSAMT.

    The module-level :data:`PLOT_CONFIG` instance is the recommended entry
    point.  Mutate it once (at the top of your script or notebook) and every
    subsequent :func:`save_fig` call inherits those settings.

    Parameters
    ----------
    fmt : str or list[str], default ``"png"``
        Export format token(s).  See :ref:`format-token-rules` above.
    base_fmt : str, default ``"png"``
        The reference format for additive ``+`` tokens.
    dpi : int, default ``150``
        Raster resolution.  Use ``300``–``600`` for publication.
    bbox_inches : str, default ``"tight"``
        Passed to :func:`~matplotlib.figure.Figure.savefig`.
    transparent : bool, default ``False``
        Transparent background.
    facecolor : str, default ``"white"``
        Figure background colour.
    savedir : str or Path or None
        If set, relative *path* arguments passed to :func:`save_fig` are
        resolved under this directory.
    close_after_save : bool, default ``False``
        Call ``plt.close(fig)`` after saving.
    verbose : bool, default ``True``
        Print a line for each file written.

    Examples
    --------
    Configure once::

        from pycsamt.api.plot import PLOT_CONFIG, save_fig
        PLOT_CONFIG.fmt = ["+svg"]
        PLOT_CONFIG.dpi = 300
        PLOT_CONFIG.savedir = "~/paper/figures"

    Use the :meth:`configure` helper (dotted-style not needed — direct
    attribute access is cleaner for this class)::

        PLOT_CONFIG.configure(fmt="+svg", dpi=300)
    """

    def __init__(
        self,
        *,
        fmt: str | list[str] = "png",
        base_fmt: str = "png",
        dpi: int = 150,
        bbox_inches: str = "tight",
        transparent: bool = False,
        facecolor: str = "white",
        savedir: str | Path | None = None,
        close_after_save: bool = False,
        verbose: bool = True,
    ) -> None:
        self.fmt             = fmt
        self.base_fmt        = base_fmt
        self.dpi             = dpi
        self.bbox_inches     = bbox_inches
        self.transparent     = transparent
        self.facecolor       = facecolor
        self.savedir         = savedir
        self.close_after_save = close_after_save
        self.verbose         = verbose

    # ── format resolution ─────────────────────────────────────────────────────

    def resolve_formats(
        self,
        fmt: str | list[str] | None = None,
    ) -> list[str]:
        """Return the resolved list of format strings for this call.

        Parameters
        ----------
        fmt : str, list[str], or None
            Per-call override.  ``None`` → use :attr:`fmt`.

        Returns
        -------
        list[str]
            Clean format strings (no ``+``, no leading dot).

        Examples
        --------
        >>> cfg = PlotConfig(fmt="+svg")
        >>> cfg.resolve_formats()
        ['png', 'svg']
        >>> cfg.resolve_formats(fmt="pdf")
        ['pdf']
        >>> cfg.resolve_formats(fmt=["+svg", "+pdf"])
        ['png', 'svg', 'pdf']
        """
        raw = fmt if fmt is not None else self.fmt
        tokens = _parse_fmt_tokens(raw)
        result = _resolve_formats(tokens, self.base_fmt)
        _validate_fmt(result)
        return result

    # ── save ─────────────────────────────────────────────────────────────────

    def save(
        self,
        fig_or_ax: Any,
        path: str | Path,
        *,
        fmt: str | list[str] | None = None,
        dpi: int | None = None,
        bbox_inches: str | None = None,
        transparent: bool | None = None,
        facecolor: str | None = None,
        **savefig_kw: Any,
    ) -> list[Path]:
        """Save *fig_or_ax* to one or more format files.

        Parameters
        ----------
        fig_or_ax : Figure or Axes
            Accepts either object type — the parent figure is extracted
            automatically from an :class:`~matplotlib.axes.Axes`.
        path : str or Path
            Base output path.  The file extension is *replaced* (or added)
            per format.  If relative and :attr:`savedir` is set, the path is
            resolved under :attr:`savedir`.
        fmt : str, list[str], or None
            Per-call format override.  ``None`` → :attr:`fmt`.
        dpi : int or None
            Per-call DPI override.  ``None`` → :attr:`dpi`.
        bbox_inches : str or None
            Per-call override.
        transparent : bool or None
            Per-call override.
        facecolor : str or None
            Per-call override.
        **savefig_kw
            Additional keyword arguments forwarded to
            :meth:`~matplotlib.figure.Figure.savefig`.

        Returns
        -------
        list[Path]
            Paths of every file written (one per resolved format).

        Examples
        --------
        >>> paths = PLOT_CONFIG.save(ax, "fig_pt")
        >>> paths = PLOT_CONFIG.save(fig, "fig_pt", fmt="+svg", dpi=300)
        """
        import matplotlib.pyplot as plt

        # resolve figure
        try:
            import matplotlib.figure as _mfig
            fig = (fig_or_ax if isinstance(fig_or_ax, _mfig.Figure)
                   else fig_or_ax.get_figure())
        except AttributeError:
            raise TypeError(
                f"Expected a matplotlib Figure or Axes, got {type(fig_or_ax)!r}."
            ) from None

        # resolve output path
        path = Path(path)
        if not path.is_absolute() and self.savedir is not None:
            path = Path(os.path.expanduser(str(self.savedir))) / path
        path.parent.mkdir(parents=True, exist_ok=True)

        # resolve per-call overrides
        _dpi         = dpi         if dpi         is not None else self.dpi
        _bbox        = bbox_inches if bbox_inches is not None else self.bbox_inches
        _transparent = transparent if transparent is not None else self.transparent
        _facecolor   = facecolor   if facecolor  is not None else self.facecolor

        fmts = self.resolve_formats(fmt)
        saved: list[Path] = []
        for f in fmts:
            out = path.with_suffix(f".{f}")
            fig.savefig(
                out,
                dpi=_dpi,
                bbox_inches=_bbox,
                transparent=_transparent,
                facecolor=_facecolor,
                format=f,
                **savefig_kw,
            )
            saved.append(out)
            if self.verbose:
                print(f"  ✔  {out}")

        if self.close_after_save:
            plt.close(fig)

        return saved

    # ── configure ─────────────────────────────────────────────────────────────

    def configure(self, **kw: Any) -> None:
        """Set multiple attributes at once.

        Parameters
        ----------
        **kw
            Any writable :class:`PlotConfig` attribute by name.

        Examples
        --------
        >>> PLOT_CONFIG.configure(fmt="+svg", dpi=300, savedir="~/figures")
        """
        for k, v in kw.items():
            if not hasattr(self, k):
                raise AttributeError(
                    f"PlotConfig has no attribute {k!r}.  "
                    f"Valid: {list(self._fields())}"
                )
            setattr(self, k, v)

    # ── reset ─────────────────────────────────────────────────────────────────

    def reset(self) -> None:
        """Reset to package defaults (ignores config file and env vars)."""
        self.fmt              = "png"
        self.base_fmt         = "png"
        self.dpi              = 150
        self.bbox_inches      = "tight"
        self.transparent      = False
        self.facecolor        = "white"
        self.savedir          = None
        self.close_after_save = False
        self.verbose          = True

    # ── context manager ───────────────────────────────────────────────────────

    @contextmanager
    def context(self, **kw: Any) -> Generator[PlotConfig, None, None]:
        """Temporarily override config attributes, then restore.

        Parameters
        ----------
        **kw
            Any writable :class:`PlotConfig` attribute by name.

        Yields
        ------
        PlotConfig
            The (temporarily modified) config object.

        Examples
        --------
        >>> with PLOT_CONFIG.context(fmt="pdf", dpi=600):
        ...     save_fig(fig, "publication/fig")
        >>> # config reverts here
        """
        saved = self._snapshot()
        try:
            self.configure(**kw)
            yield self
        finally:
            self._restore(saved)

    # ── serialise / summarise ─────────────────────────────────────────────────

    def to_dict(self) -> dict:
        """Return current config as a plain dict."""
        return {f: getattr(self, f) for f in self._fields()}

    def summary(self) -> str:
        """Return a human-readable summary of the current config."""
        resolved = self.resolve_formats()
        lines = [
            "PlotConfig",
            f"  fmt              = {self.fmt!r}",
            f"  resolved formats = {resolved}",
            f"  base_fmt         = {self.base_fmt!r}",
            f"  dpi              = {self.dpi}",
            f"  bbox_inches      = {self.bbox_inches!r}",
            f"  transparent      = {self.transparent}",
            f"  facecolor        = {self.facecolor!r}",
            f"  savedir          = {self.savedir!r}",
            f"  close_after_save = {self.close_after_save}",
            f"  verbose          = {self.verbose}",
        ]
        return "\n".join(lines)

    def __repr__(self) -> str:  # noqa: D105
        return self.summary()

    # ── helpers ───────────────────────────────────────────────────────────────

    def _fields(self) -> list[str]:
        return [
            "fmt", "base_fmt", "dpi", "bbox_inches",
            "transparent", "facecolor", "savedir",
            "close_after_save", "verbose",
        ]

    def _snapshot(self) -> dict:
        return copy.deepcopy(self.to_dict())

    def _restore(self, snap: dict) -> None:
        for k, v in snap.items():
            setattr(self, k, v)


# ═══════════════════════════════════════════════════════════════════════════════
# Build the module-level singleton
# ═══════════════════════════════════════════════════════════════════════════════

def _build_singleton() -> PlotConfig:
    """Construct PLOT_CONFIG from defaults → config file → env vars."""
    cfg = PlotConfig()

    # layer 1: config file (CWD overrides home-dir)
    file_kw = _load_config_file()
    for k, v in file_kw.items():
        setattr(cfg, k, v)

    # layer 2: environment variables (highest precedence)
    env_kw = _load_env()
    for k, v in env_kw.items():
        setattr(cfg, k, v)

    return cfg


#: Package-level singleton.  Import and mutate to configure all
#: :func:`save_fig` calls globally.
PLOT_CONFIG: PlotConfig = _build_singleton()


# ═══════════════════════════════════════════════════════════════════════════════
# Convenience functions
# ═══════════════════════════════════════════════════════════════════════════════

def save_fig(
    fig_or_ax: Any,
    path: str | Path,
    *,
    fmt: str | list[str] | None = None,
    dpi: int | None = None,
    **kw: Any,
) -> list[Path]:
    """Save a pyCSAMT figure using the global :data:`PLOT_CONFIG` settings.

    Parameters
    ----------
    fig_or_ax : Figure or Axes
        The object to save.
    path : str or Path
        Base output path (extension auto-managed per format).
    fmt : str, list[str], or None
        Per-call format override.  ``None`` → :attr:`PLOT_CONFIG.fmt`.
    dpi : int or None
        Per-call DPI override.  ``None`` → :attr:`PLOT_CONFIG.dpi`.
    **kw
        Forwarded to :meth:`PlotConfig.save`.

    Returns
    -------
    list[Path]
        File paths written.

    Examples
    --------
    Default (PNG)::

        save_fig(ax, "output/fig_pt")

    Override format for this call::

        save_fig(ax, "output/fig_pt", fmt="svg")

    Two formats for this call::

        save_fig(ax, "output/fig_pt", fmt=["+svg"])
    """
    return PLOT_CONFIG.save(fig_or_ax, path, fmt=fmt, dpi=dpi, **kw)


def add_colorbar(
    mappable: Any,
    ax: Any,
    *,
    label: str | None = None,
    side: str = "right",
    size: str = "3.5%",
    pad: float | str = 0.06,
    max_ticks: int | None = 6,
    tick_format: str | None = None,
    **colorbar_kw: Any,
) -> Any:
    """Attach an axes-aligned colorbar with smart tick density.

    The helper creates a dedicated colorbar axes whose height matches the
    plotted axes.  This avoids overly tall colorbars on equal-aspect maps
    and gives package plots a consistent colorbar geometry.
    """
    from matplotlib.ticker import (
        FormatStrFormatter,
        MaxNLocator,
    )
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    side = side.lower()
    if side not in {"right", "left", "top", "bottom"}:
        msg = "colorbar side must be 'right', 'left', 'top', or 'bottom'."
        raise ValueError(msg)

    divider = make_axes_locatable(ax)
    orientation = (
        "vertical" if side in {"right", "left"} else "horizontal"
    )
    cax = divider.append_axes(side, size=size, pad=pad)
    cbar = ax.figure.colorbar(
        mappable,
        cax=cax,
        orientation=orientation,
        **colorbar_kw,
    )
    if label:
        cbar.set_label(label)
    if max_ticks is not None:
        locator = MaxNLocator(nbins=max(2, int(max_ticks)))
        cbar.locator = locator
    if tick_format:
        formatter = FormatStrFormatter(tick_format)
        cbar.formatter = formatter
    cbar.update_ticks()
    return cbar


def add_polar_colorbar(
    mappable: Any,
    ax: Any,
    *,
    label: str | None = None,
    pad: float = 0.10,
    shrink: float = 0.72,
    aspect: float = 20,
    max_ticks: int | None = 5,
    tick_format: str | None = None,
    **colorbar_kw: Any,
) -> Any:
    """Attach a compact colorbar beside a polar axes.

    Polar plots do not work well with the axes-divider geometry used by
    :func:`add_colorbar`, especially when figures are saved with tight
    bounding boxes.  This helper keeps the same tick-density policy while
    using Matplotlib's polar-friendly colorbar placement.
    """
    from matplotlib.ticker import (
        FormatStrFormatter,
        MaxNLocator,
    )

    cbar = ax.figure.colorbar(
        mappable,
        ax=ax,
        pad=pad,
        shrink=shrink,
        aspect=aspect,
        **colorbar_kw,
    )
    if label:
        cbar.set_label(label)
    if max_ticks is not None:
        cbar.locator = MaxNLocator(nbins=max(2, int(max_ticks)))
    if tick_format:
        cbar.formatter = FormatStrFormatter(tick_format)
    cbar.update_ticks()
    return cbar


def set_fmt(*formats: str) -> None:
    """Set the global export format(s) on :data:`PLOT_CONFIG`.

    Parameters
    ----------
    *formats : str
        One or more format tokens.  A single ``"png"`` or ``"svg"`` sets
        that format exclusively.  A ``"+svg"`` token adds SVG to the base
        format.

    Examples
    --------
    >>> set_fmt("svg")          # only SVG
    >>> set_fmt("+svg")         # PNG + SVG
    >>> set_fmt("+svg", "+pdf") # PNG + SVG + PDF
    >>> set_fmt("svg", "pdf")   # SVG + PDF (no PNG)
    """
    PLOT_CONFIG.fmt = list(formats) if len(formats) > 1 else formats[0]


def set_dpi(dpi: int) -> None:
    """Set the global raster DPI on :data:`PLOT_CONFIG`.

    Parameters
    ----------
    dpi : int
        Dots per inch.  Typical values: ``150`` (screen), ``300`` (print),
        ``600`` (high-res print).

    Examples
    --------
    >>> set_dpi(300)
    """
    PLOT_CONFIG.dpi = int(dpi)


def set_savedir(path: str | Path) -> None:
    """Set the global output directory on :data:`PLOT_CONFIG`.

    Parameters
    ----------
    path : str or Path
        Directory where :func:`save_fig` writes files when a relative
        *path* is given.  Tilde (``~``) is expanded.

    Examples
    --------
    >>> set_savedir("~/paper/figures")
    >>> save_fig(ax, "fig_pt")   # → ~/paper/figures/fig_pt.png
    """
    PLOT_CONFIG.savedir = os.path.expanduser(str(path))


def reset_plot_config() -> None:
    """Reset :data:`PLOT_CONFIG` to package defaults.

    Does **not** re-read config files or environment variables.

    Examples
    --------
    >>> reset_plot_config()
    """
    PLOT_CONFIG.reset()


def load_plot_config(path: str | Path | None = None) -> None:
    """(Re-)load a config file into :data:`PLOT_CONFIG`.

    Parameters
    ----------
    path : str, Path, or None
        Path to an INI config file with a ``[plot]`` section.  ``None``
        searches the default locations (``pycsamt_plot.cfg`` in CWD then
        ``~/.pycsamt/plot.cfg``).

    Raises
    ------
    FileNotFoundError
        If *path* is given explicitly and does not exist.

    Examples
    --------
    Load a specific file::

        load_plot_config("project/plot_settings.cfg")

    Reload defaults::

        load_plot_config()
    """
    if path is not None:
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"Plot config file not found: {p}")
        cp = configparser.ConfigParser()
        cp.read(p)
        kw: dict = {}
        if "plot" in cp:
            s = cp["plot"]
            if "fmt" in s:
                kw["fmt"] = [t.strip() for t in s["fmt"].split(",")]
            if "base_fmt" in s:
                kw["base_fmt"] = s["base_fmt"].strip()
            if "dpi" in s:
                kw["dpi"] = int(s["dpi"])
            if "bbox_inches" in s:
                kw["bbox_inches"] = s["bbox_inches"].strip()
            if "transparent" in s:
                kw["transparent"] = s.getboolean("transparent")
            if "facecolor" in s:
                kw["facecolor"] = s["facecolor"].strip()
            if "savedir" in s:
                kw["savedir"] = os.path.expanduser(s["savedir"].strip())
            if "close_after_save" in s:
                kw["close_after_save"] = s.getboolean("close_after_save")
        PLOT_CONFIG.configure(**kw)
    else:
        kw = _load_config_file()
        PLOT_CONFIG.configure(**kw)
    # always let env vars win
    for k, v in _load_env().items():
        setattr(PLOT_CONFIG, k, v)


def write_default_config(path: str | Path = "pycsamt_plot.cfg") -> Path:
    """Write a template config file with the current :data:`PLOT_CONFIG` values.

    Useful for users who want to tweak settings without touching Python code.

    Parameters
    ----------
    path : str or Path
        Destination file path.  Default ``"pycsamt_plot.cfg"`` writes to the
        current working directory.

    Returns
    -------
    Path
        The file that was written.

    Examples
    --------
    >>> write_default_config()          # creates pycsamt_plot.cfg in CWD
    >>> write_default_config("~/.pycsamt/plot.cfg")   # user-global config
    """
    p = Path(os.path.expanduser(str(path)))
    p.parent.mkdir(parents=True, exist_ok=True)
    cfg = PLOT_CONFIG
    fmt_str = (
        ", ".join(cfg.fmt) if isinstance(cfg.fmt, list) else str(cfg.fmt)
    )
    content = f"""\
# pyCSAMT plot-export configuration
# Generated by pycsamt.api.plot.write_default_config()
#
# Format tokens
#   png / svg / pdf / eps / tiff / jpg   → save in that format only
#   +svg / +pdf                          → add to the base format (see base_fmt)
#   Multiple formats: fmt = +svg, +pdf
#
# Environment variable overrides (always take precedence):
#   PYCSAMT_FMT=+svg,+pdf  PYCSAMT_DPI=300  PYCSAMT_SAVEDIR=~/figures

[plot]
fmt             = {fmt_str}
base_fmt        = {cfg.base_fmt}
dpi             = {cfg.dpi}
bbox_inches     = {cfg.bbox_inches}
transparent     = {str(cfg.transparent).lower()}
facecolor       = {cfg.facecolor}
savedir         ={"" if cfg.savedir is None else " " + str(cfg.savedir)}
close_after_save = {str(cfg.close_after_save).lower()}
"""
    p.write_text(content)
    print(f"  ✔  config written → {p}")
    return p


__all__ = [
    # main class
    "PlotConfig",
    # singleton
    "PLOT_CONFIG",
    # convenience setters
    "add_colorbar",
    "add_polar_colorbar",
    "save_fig",
    "set_fmt",
    "set_dpi",
    "set_savedir",
    "reset_plot_config",
    # config-file helpers
    "load_plot_config",
    "write_default_config",
]
