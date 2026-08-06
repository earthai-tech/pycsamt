from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PathLike = Union[str, Path]


def _as_1d(a: Any) -> np.ndarray:
    v = np.asarray(a, dtype=float).ravel()
    m = np.isfinite(v) & (v > 0.0)
    return v[m]


def _prep_curve(
    misfit: Any,
    rough: Any,
    lam: Any | None,
    *,
    sort: str = "auto",  # auto|x|lambda
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = _as_1d(rough)
    y = _as_1d(misfit)
    n = int(min(x.size, y.size))
    x, y = x[:n], y[:n]
    if lam is None:
        l = np.arange(n, dtype=float)
    else:
        l = np.asarray(lam, dtype=float).ravel()[:n]
    if sort == "x":
        idx = np.argsort(x)
    elif sort == "lambda":
        idx = np.argsort(l)
    else:
        mono = np.all(np.diff(l) >= 0) or np.all(np.diff(l) <= 0)
        idx = np.argsort(l) if mono else np.argsort(x)
    return x[idx], y[idx], l[idx]


def _movavg(v: np.ndarray, w: int) -> np.ndarray:
    w = int(max(1, w))
    if v.size == 0 or w == 1:
        return v
    k = np.ones(w, dtype=float) / float(w)
    return np.convolve(v, k, mode="same")


def _curvature_kappa(
    x: np.ndarray,
    y: np.ndarray,
    *,
    smooth: int = 3,
) -> np.ndarray:
    lx = np.log10(np.maximum(x, 1e-300))
    ly = np.log10(np.maximum(y, 1e-300))
    if smooth > 1:
        lx = _movavg(lx, smooth)
        ly = _movavg(ly, smooth)
    t = np.arange(lx.size, dtype=float)
    x1 = np.gradient(lx, t, edge_order=2)
    y1 = np.gradient(ly, t, edge_order=2)
    x2 = np.gradient(x1, t, edge_order=2)
    y2 = np.gradient(y1, t, edge_order=2)
    den = (x1 * x1 + y1 * y1 + 1e-24) ** 1.5
    k = np.abs(x1 * y2 - y1 * x2) / den
    return k


def _maxdist_knee(
    x: np.ndarray,
    y: np.ndarray,
) -> np.ndarray:
    lx = np.log10(np.maximum(x, 1e-300))
    ly = np.log10(np.maximum(y, 1e-300))
    p0 = np.array([lx[0], ly[0]], dtype=float)
    p1 = np.array([lx[-1], ly[-1]], dtype=float)
    v = p1 - p0
    den = np.hypot(v[0], v[1]) + 1e-24
    # signed distance to end-point line
    d = np.abs(v[1] * (lx - p0[0]) - v[0] * (ly - p0[1])) / den
    return d


def _pick_corner(
    x: np.ndarray,
    y: np.ndarray,
    *,
    method: str = "curvature",  # curvature|maxdist
    smooth: int = 3,
    skip: int = 1,
) -> tuple[int, np.ndarray]:
    n = int(min(x.size, y.size))
    if n <= 2:
        return 0, np.zeros(n)
    if method == "maxdist":
        s = _maxdist_knee(x, y)
    else:
        s = _curvature_kappa(x, y, smooth=smooth)
    lo = int(max(0, skip))
    hi = int(max(lo + 1, n - skip))
    j = int(lo + np.nanargmax(s[lo:hi]))
    return j, s


def lcurve_table(
    misfit: Any,
    rough: Any,
    lam: Any | None = None,
    *,
    sort: str = "auto",
    method: str = "curvature",
    smooth: int = 3,
    skip: int = 1,
    return_dict: bool = False,
):
    x, y, l = _prep_curve(misfit, rough, lam, sort=sort)
    if x.size == 0:
        if return_dict:
            return dict(rough=x, misfit=y, lam=l, curv=x, slope=x, corner=None)
        if pd is None:
            return None
        df = pd.DataFrame(columns=["rough", "misfit", "lam", "curv", "slope"])
        df.attrs["corner_idx"] = None
        return df
    j, s = _pick_corner(x, y, method=method, smooth=smooth, skip=skip)
    lx = np.log10(x)
    ly = np.log10(y)
    # local slope d log(misfit) / d log(rough)
    sl = np.gradient(ly) / (np.gradient(lx) + 1e-24)
    if return_dict:
        return dict(rough=x, misfit=y, lam=l, curv=s, slope=sl, corner=j)
    df = pd.DataFrame(dict(rough=x, misfit=y, lam=l, curv=s, slope=sl))
    df.attrs["corner_idx"] = j
    return df


# ------------------------- real-inversion adapters ----------------------- #


@dataclass
class LCurveData:
    """Misfit/roughness/lambda sweep extracted from a real inversion log.

    ``misfit[i]``, ``rough[i]``, ``lam[i]``, and ``iterations[i]`` all
    describe the same iteration, so the arrays can be passed straight
    into :func:`lcurve_table` or :func:`plot_lcurve` -- which is exactly
    what :meth:`table` and :meth:`plot` do.
    """

    misfit: np.ndarray
    rough: np.ndarray
    lam: np.ndarray
    iterations: np.ndarray
    backend: str
    source: Path

    def table(self, **kwargs: Any):
        """Call :func:`lcurve_table` on this sweep."""
        return lcurve_table(self.misfit, self.rough, self.lam, **kwargs)

    def plot(self, **kwargs: Any):
        """Call :func:`plot_lcurve` on this sweep."""
        kwargs.setdefault("labels", [self.backend])
        return plot_lcurve(self.misfit, self.rough, self.lam, **kwargs)


def _finite_positive_mask(*arrays: np.ndarray) -> np.ndarray:
    m = np.ones(arrays[0].shape, dtype=bool)
    for a in arrays:
        m &= np.isfinite(a) & (a > 0.0)
    return m


def lcurve_from_occam2d(path: PathLike, **kwargs: Any) -> LCurveData:
    """Build an :class:`LCurveData` sweep from an Occam2D convergence log.

    Reads the per-iteration ``LogFile.logfile`` written by the Occam2D
    Fortran binary with :class:`pycsamt.models.occam2d.log.OccamLog`,
    which already parses accepted RMS misfit, roughness, and Lagrange
    multiplier for every completed iteration. Rows with a non-finite or
    non-positive misfit or roughness (for example a final iteration that
    stopped on "Convergence problems" before writing ``ROUGHNESS IS``)
    are dropped, since :func:`lcurve_table` requires strictly positive
    values for its log-log scoring.

    Parameters
    ----------
    path : path-like
        Path to an Occam2D log file, typically ``LogFile.logfile``.
    **kwargs
        Forwarded to ``OccamLog.read`` (for example ``verbose``).

    Returns
    -------
    LCurveData
        ``rough`` holds Occam's reported ``ROUGHNESS IS`` values,
        ``misfit`` the accepted RMS, and ``lam`` the accepted Lagrange
        multiplier (linear scale, not log10).

    Examples
    --------
    >>> from pycsamt.emtools.lcurve import lcurve_from_occam2d
    >>> sweep = lcurve_from_occam2d("data/occam2D/LogFile.logfile")
    >>> sweep.backend
    'occam2d'
    >>> ax = sweep.plot()
    """
    from pycsamt.models.occam2d.log import OccamLog

    log = OccamLog.read(path, **kwargs)
    mask = _finite_positive_mask(log.rms, log.roughness)
    return LCurveData(
        misfit=log.rms[mask],
        rough=log.roughness[mask],
        lam=log.lagrange[mask],
        iterations=log.iterations[mask],
        backend="occam2d",
        source=Path(path),
    )


def lcurve_from_modem(path: PathLike, **kwargs: Any) -> LCurveData:
    """Build an :class:`LCurveData` sweep from a ModEM NLCG log.

    Reads the per-iteration ``Completed NLCG iteration`` blocks written
    by ModEM with :class:`pycsamt.models.modem.log.ModEmLog`. ModEM
    reports the model-regularization term directly as ``m2`` rather than
    a roughness norm, so ``rough`` here is ``m2`` -- the same quantity
    :math:`\\lambda\\,\\Phi_m(m)` trades off against data misfit in
    ModEM's objective function, and it plays the identical x-axis role
    on an L-curve.

    Parameters
    ----------
    path : path-like
        Path to a ModEM log file, typically ``Modular_NLCG.log``.
    **kwargs
        Forwarded to ``ModEmLog.read`` (for example ``verbose``).

    Returns
    -------
    LCurveData
        ``rough`` holds ModEM's ``m2`` model-norm term, ``misfit`` the
        reported ``rms``, and ``lam`` the damping parameter ``lambda``.

    Examples
    --------
    >>> from pycsamt.emtools.lcurve import lcurve_from_modem
    >>> sweep = lcurve_from_modem(
    ...     "data/modem/willy_27freq_watex_line02_sample/Modular_NLCG.log"
    ... )
    >>> sweep.backend
    'modem'
    >>> ax = sweep.plot()
    """
    from pycsamt.models.modem.log import ModEmLog

    log = ModEmLog.read(path, **kwargs)
    mask = _finite_positive_mask(log.rms, log.model_norm)
    return LCurveData(
        misfit=log.rms[mask],
        rough=log.model_norm[mask],
        lam=log.lagrange[mask],
        iterations=log.iterations[mask],
        backend="modem",
        source=Path(path),
    )


def lcurve_from_mare2dem(path: PathLike, **kwargs: Any) -> LCurveData:
    """Build an :class:`LCurveData` sweep from a MARE2DEM convergence log.

    Reads the per-iteration ``** Iteration N **`` blocks written by
    MARE2DEM with :class:`pycsamt.models.mare2dem.log.Mare2DEMLog`,
    which already parses ``Model Misfit``, ``Roughness``, and
    ``Optimal Mu`` for every completed iteration. MARE2DEM reports
    ``Optimal Mu`` as :math:`\\log_{10}\\mu`, so it is converted back to
    linear scale here to match the convention used by the other two
    adapters.

    Parameters
    ----------
    path : path-like
        Path to a MARE2DEM log file, typically ``*.logfile``.
    **kwargs
        Accepted for interface symmetry with the other adapters;
        ``Mare2DEMLog`` takes no extra keyword arguments.

    Returns
    -------
    LCurveData
        ``rough`` holds MARE2DEM's reported roughness, ``misfit`` the
        model misfit, and ``lam`` the optimal mu converted to linear
        scale (``10 ** log10_mu``).

    Examples
    --------
    >>> from pycsamt.emtools.lcurve import lcurve_from_mare2dem
    >>> sweep = lcurve_from_mare2dem(
    ...     "data/mare2dem/demo_mt_inversion/demo.logfile"
    ... )
    >>> sweep.backend
    'mare2dem'
    >>> ax = sweep.plot()
    """
    from pycsamt.models.mare2dem.log import Mare2DEMLog

    del kwargs  # accepted for interface symmetry only
    log = Mare2DEMLog(path)
    misfit = np.array([r.rms for r in log.iterations], dtype=float)
    rough = np.array([r.roughness for r in log.iterations], dtype=float)
    lam = 10.0 ** np.array([r.lambda_ for r in log.iterations], dtype=float)
    iterations = np.array([r.iteration for r in log.iterations], dtype=int)
    mask = _finite_positive_mask(misfit, rough)
    return LCurveData(
        misfit=misfit[mask],
        rough=rough[mask],
        lam=lam[mask],
        iterations=iterations[mask],
        backend="mare2dem",
        source=Path(path),
    )


# ------------------------------- plotting ------------------------------- #


def plot_lcurve(
    misfit: Any | Sequence[Any],
    rough: Any | Sequence[Any],
    lam: Any | Sequence[Any] | None = None,
    *,
    labels: Sequence[str] | None = None,
    colors: Sequence[str] | None = None,
    cmap: str = "viridis",
    marker: str = "o",
    ms: float = 3.0,
    lw: float = 1.4,
    alpha: float = 0.9,
    show_points: bool = True,
    show_path: bool = True,
    arrow_every: int = 0,  # 0 disables
    method: str = "curvature",
    smooth: int = 3,
    skip: int = 1,
    show_corner: bool = True,
    corner_style: dict[str, Any] | None = None,
    show_inset: bool = True,
    inset_loc: tuple[float, float, float, float] = (0.62, 0.12, 0.32, 0.32),
    label_every: int = 0,  # 0 disables per-point lambda/tau labels
    label_prefix: str = "",
    label_fontsize: float = 7.0,
    target_misfit: float | None = None,
    target_label: str = "target misfit",
    figsize: tuple[float, float] = (6.0, 4.6),
    ax: plt.Axes | None = None,
):
    # normalize to list of curves
    if isinstance(misfit, (list, tuple)):
        Ms = list(misfit)
        Rs = list(rough)  # type: ignore
        Ls = list(lam) if isinstance(lam, (list, tuple)) else [None] * len(Ms)
        labs = (
            list(labels)
            if labels is not None
            else [f"C{i}" for i in range(len(Ms))]
        )
    else:
        Ms = [misfit]
        Rs = [rough]
        Ls = [lam]
        labs = [labels[0]] if labels else ["curve"]
    if colors is None:
        cols = [None] * len(Ms)
    else:
        cols = list(colors)

    # main axes
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # store corners for legend text
    leg = []
    # Shared inset axes: created once and reused across curves, rather
    # than once per curve at the identical `inset_loc` (which used to
    # silently stack N insets on top of each other, leaving only the
    # last curve's score visible).
    axins = None

    for i, (m, r, la) in enumerate(zip(Ms, Rs, Ls)):
        x, y, l = _prep_curve(m, r, la, sort="auto")
        if x.size == 0:
            continue
        j, score = _pick_corner(x, y, method=method, smooth=smooth, skip=skip)
        c = cols[i]
        if c is None:
            c = plt.get_cmap(cmap)(0.15 + 0.65 * i / max(1, len(Ms) - 1))
        # path and points
        ax.set_xscale("log")
        ax.set_yscale("log")
        if show_path:
            ax.plot(x, y, "-", color=c, lw=lw, alpha=alpha)
        if show_points:
            ax.scatter(
                x,
                y,
                c=np.linspace(0, 1, x.size),
                s=12 * ms,
                cmap=cmap,
                edgecolors="none",
                alpha=0.85,
            )
        # optional arrows showing λ direction
        if arrow_every and x.size > 2:
            step = int(max(1, arrow_every))
            for k in range(0, x.size - step, step):
                ax.annotate(
                    "",
                    xy=(x[k + step], y[k + step]),
                    xytext=(x[k], y[k]),
                    arrowprops=dict(
                        arrowstyle="->", lw=0.8, color=c, shrinkA=0, shrinkB=0
                    ),
                )
        # per-point lambda/tau labels, e.g. "τ=100" at every Nth point
        if label_every and x.size > 0:
            step = int(max(1, label_every))
            for k in range(0, x.size, step):
                ax.annotate(
                    f"{label_prefix}{l[k]:.3g}",
                    xy=(x[k], y[k]),
                    xytext=(4, 4),
                    textcoords="offset points",
                    fontsize=label_fontsize,
                    color=c,
                )

        # corner highlight
        if show_corner:
            sty = dict(marker="*", ms=9, mec="k", mfc=c, mew=0.8)
            if corner_style:
                sty.update(corner_style)
            ax.plot([x[j]], [y[j]], **sty)
            leg.append(f"{labs[i]}  λ*≈{l[j]:.3g}")
        else:
            leg.append(labs[i])

        # inset curvature / distance -- one shared inset for all curves
        if show_inset:
            if axins is None:
                axins = fig.add_axes(inset_loc)
                axins.set_xticks([])
                axins.set_yticks([])
                axins.set_title(
                    "curv" if method == "curvature" else "knee",
                    fontsize=8,
                )
            idx = np.arange(score.size)
            axins.plot(idx, score, "-", color=c, lw=1.0)
            axins.axvline(j, color=c, ls="--", lw=1.0)

    if target_misfit is not None:
        ax.axhline(
            target_misfit, color="0.2", ls=":", lw=1.2, label=target_label
        )
        leg.append(target_label)

    ax.grid(True, alpha=0.25, which="both")
    ax.set_xlabel("||Lm|| (model roughness)")
    ax.set_ylabel("||Gm−d|| (data misfit)")
    if leg:
        ax.legend(leg, loc="best", frameon=False, fontsize=9)
    return ax
