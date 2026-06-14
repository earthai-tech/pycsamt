
from __future__ import annotations

from typing import Dict, Optional, Tuple, List, Any
import numpy as np

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.patches import Ellipse
from mpl_toolkits.axes_grid1 import make_axes_locatable


from ..site.base import Sites
from ..site import edit as _edit
from ..site import compute as _compute
from ..z import utils as zutils

from ._core import (
    ensure_sites,
    _iter_items,
    _name,
    _get_z_block,
    _get_t_block,
    _apply_each
)
from ..api._rose_style import RoseStyle, resolve_rose_style, _UNSET
from ..api.style import PYCSAMT_STYLE

_BACKWARD_SINCE = "2.0.0"
_BACKWARD_REMOVE = "2.17.0"

#
def rotate_to_strike(
    sites,
    *,
    method: str = "swift",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Sites:
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    def _one(Si: Sites) -> Sites:
        ang = 0.0
        try:
            r = _compute.strike_estimate(Si, method=method)
            if isinstance(r, (int, float)):
                ang = float(r)
            elif isinstance(r, dict):
                ang = float(r.get("angle", 0.0))
            else:
                ang = float(getattr(r, "angle", 0.0))
        except Exception:
            ang = 0.0
        return _edit.rotate(
            Si,
            ang,
            inplace=inplace,
        )

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# --- 2) rotate by angle or map ------------------------------------------- #

def rotate(
    sites,
    angle: float,
    *,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Sites:

    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    return _edit.rotate(S, angle=angle, inplace=inplace)


def rotate_by_map(
    sites,
    angle_by_station: dict,
    *,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Sites:

    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    def _one(Si: Sites) -> Sites:
        ed = next(_iter_items(Si))
        name = getattr(ed, "station", None) or getattr(ed, "name", None)
        ang = float(angle_by_station.get(name, 0.0))
        return _edit.rotate(Si, angle=ang, inplace=inplace)

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# --- 3) enforce off-diagonal antisymmetry -------------------------------- #

def antisymmetrize(
    sites,
    *,
    how: str = "rms",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Sites:

    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    def _one(Si: Sites) -> Sites:
        ed = next(_iter_items(Si))
        Z = getattr(ed, "Z", None) or getattr(ed, "z", None)
        if Z is None:
            return Si
        z = getattr(Z, "z", None)
        ze = getattr(Z, "z_err", None)
        if z is None:
            return Si
        z2, ze2 = zutils.enforce_offdiag_antisymmetry(
            z, ze, how=how
        )
        try:
            Z.z = z2
            if ze is not None and ze2 is not None:
                Z.z_err = ze2
        except Exception:
            pass
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# --- 4) invert impedance -------------------------------------------------- #

def invert(
    sites,
    *,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Sites:

    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    def _one(Si: Sites) -> Sites:
        ed = next(_iter_items(Si))
        Z = getattr(ed, "Z", None) or getattr(ed, "z", None)
        if Z is None:
            return Si
        z = getattr(Z, "z", None)
        ze = getattr(Z, "z_err", None)
        if z is None:
            return Si
        z2, ze2 = zutils.invert_z(z, ze)
        try:
            Z.z = z2
            if ze is not None and ze2 is not None:
                Z.z_err = ze2
        except Exception:
            pass
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# --- 5) sensor-orientation correction ------------------------------------ #

def orient_from_sensors(
    sites,
    ex: float,
    ey: float,
    bx: float,
    by: float,
    *,
    degrees: bool = True,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Sites:

    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    def _one(Si: Sites) -> Sites:
        ed = next(_iter_items(Si))
        Z = getattr(ed, "Z", None) or getattr(ed, "z", None)
        if Z is None:
            return Si
        z = getattr(Z, "z", None)
        ze = getattr(Z, "z_err", None)
        if z is None:
            return Si
        z2, ze2 = zutils.correct_for_sensor_orientation(
            z,
            ex=ex, ey=ey, bx=bx, by=by,
            degrees=degrees,
            z_err=ze,
        )
        try:
            Z.z = z2
            if ze is not None and ze2 is not None:
                Z.z_err = ze2
        except Exception:
            pass
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# --- 6) sigma-clip outliers in Z ----------------------------------------- #

def sigma_clip_z(
    sites,
    *,
    sigma: float = 3.0,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Sites:

    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    def _one(Si: Sites) -> Sites:
        ed = next(_iter_items(Si))
        Z = getattr(ed, "Z", None) or getattr(ed, "z", None)
        if Z is None:
            return Si
        z = getattr(Z, "z", None)
        if z is None:
            return Si
        mask = zutils.sigma_clip_mask(z, sigma=sigma)
        z2 = z.copy()
        z2[~mask] = np.nan
        try:
            Z.z = z2
        except Exception:
            pass
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# --- 7) balance off-diagonal magnitudes ---------------------------------- #

def balance_offdiag(
    sites,
    *,
    mode: str = "avgabs",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Sites:

    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    def _one(Si: Sites) -> Sites:
        ed = next(_iter_items(Si))
        Z = getattr(ed, "Z", None) or getattr(ed, "z", None)
        if Z is None:
            return Si
        z = getattr(Z, "z", None)
        if z is None:
            return Si
        zx = z[:, 0, 1]
        zy = z[:, 1, 0]
        if mode == "avgabs":
            m = 0.5 * (np.abs(zx) + np.abs(zy)) + 1e-12
            sx = m / (np.abs(zx) + 1e-12)
            sy = m / (np.abs(zy) + 1e-12)
            z2 = z.copy()
            z2[:, 0, 1] = zx * sx
            z2[:, 1, 0] = zy * sy
        else:
            z2 = z
        try:
            Z.z = z2
        except Exception:
            pass
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)




def _get_z(ed: Any) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    Zobj = getattr(ed, "Z", None) or getattr(ed, "z", None)
    if Zobj is None:
        return None, None
    z = getattr(Zobj, "z", None)
    fr = getattr(Zobj, "freq", None)
    if isinstance(z, np.ndarray) and isinstance(fr, np.ndarray):
        return z, fr
    if isinstance(Zobj, np.ndarray) and Zobj.ndim == 3:
        fr2 = getattr(ed, "freq", None)
        if isinstance(fr2, np.ndarray):
            return Zobj, fr2
    return None, None


def _period(fr: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        p = 1.0 / fr
    return p


def _angles_deg(a: np.ndarray, b: np.ndarray,
                c: np.ndarray, d: np.ndarray) -> Tuple[
                    np.ndarray, np.ndarray]:
    num_b = b + c
    den_b = a - d
    beta = 0.5 * np.degrees(np.arctan2(num_b, den_b))
    num_a = -(b - c)
    den_a = a + d
    alpha = 0.5 * np.degrees(np.arctan2(num_a, den_a))
    return alpha, beta


def _phi_from_z(z: np.ndarray) -> np.ndarray:
    # z: (n, 2, 2) complex
    X = z.real
    Y = z.imag
    n = z.shape[0]
    phi = np.empty((n, 2, 2), dtype=float)
    for i in range(n):
        Xi = X[i]
        Yi = Y[i]
        try:
            invX = np.linalg.inv(Xi)
        except np.linalg.LinAlgError:
            invX = np.linalg.pinv(Xi)
        phi[i] = invX @ Yi
    return phi


def _svd_axes(phi: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    # returns s (n,2) and angle of major axis (deg), shape (n,)
    n = phi.shape[0]
    s = np.empty((n, 2), dtype=float)
    ang = np.empty(n, dtype=float)
    for i in range(n):
        U, S, VT = np.linalg.svd(phi[i])
        s[i] = S
        vx = U[0, 0]
        vy = U[1, 0]
        ang[i] = np.degrees(np.arctan2(vy, vx))
    return s, ang


def _ellipticity(s: np.ndarray) -> np.ndarray:
    num = s[:, 0] - s[:, 1]
    den = s[:, 0] + s[:, 1] + 1e-12
    return num / den


def _df_phase_tensor_single(
    st: str,
    fr: np.ndarray,
    phi: np.ndarray,
    s: np.ndarray,
    theta: np.ndarray,
    alpha: np.ndarray,
    beta: np.ndarray,
) -> pd.DataFrame:
    p = _period(fr)
    e = _ellipticity(s)
    rec = dict(
        station=[st] * len(fr),
        freq=fr,
        period=p,
        s1=s[:, 0],
        s2=s[:, 1],
        theta=theta,
        alpha=alpha,
        beta=beta,
        skew=beta,
        ellipt=e,
    )
    return pd.DataFrame(rec)


# ------------------------ public computation ----------------------------- #

def build_phase_tensor_table(
    sites: Any,
    *,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    dfs: List[pd.DataFrame] = []
    for i, ed in enumerate(_iter_items(S)):
        nm = _name(ed, i)
        z, fr = _get_z(ed)
        if z is None or fr is None:
            continue
        phi = _phi_from_z(z)
        a = phi[:, 0, 0]
        b = phi[:, 0, 1]
        c = phi[:, 1, 0]
        d = phi[:, 1, 1]
        alpha, beta = _angles_deg(a, b, c, d)
        s, theta = _svd_axes(phi)
        df = _df_phase_tensor_single(
            nm, fr, phi, s, theta, alpha, beta
        )
        dfs.append(df)
    if not dfs:
        return pd.DataFrame(
            columns=[
                "station", "freq", "period", "s1", "s2",
                "theta", "alpha", "beta", "skew", "ellipt",
            ]
        )
    return pd.concat(dfs, ignore_index=True)


# ─── Plotting helpers ─────────────────────────────────────────────────────────

# Human-readable labels for each c_by specifier
_C_LABEL: Dict[str, str] = {
    "skew":     "Skew β (°)",
    "beta":     "Skew β (°)",
    "alpha":    "α (°)",
    "theta":    "Strike θ (°)",
    "ellipt":   "Ellipticity",
    "s1":       "φ_max (s.v.)",
    "s2":       "φ_min (s.v.)",
    "|skew|":   "|β| (°)",
    "|beta|":   "|β| (°)",
    "|theta|":  "|θ| (°)",
    "|alpha|":  "|α| (°)",
    "phi_mean": "Φ̄ = (φ_max+φ_min)/2",
    "phi_max":  "φ_max (tan units)",
    "phi_min":  "φ_min (tan units)",
}

# Symmetric c_by quantities (vmin = -vmax enforced by default)
_SYMMETRIC_C: frozenset = frozenset(("skew", "beta", "alpha", "|skew|", "|beta|"))


def _resolve_cvals(df: pd.DataFrame, c_by: str) -> Tuple[np.ndarray, str]:
    """Return ``(colour_array, colorbar_label)`` for a *c_by* specifier."""
    if c_by.startswith("|") and c_by.endswith("|"):
        inner = c_by[1:-1]
        col = (df[inner] if inner in df.columns else df["skew"]).to_numpy(float)
        return np.abs(col), _C_LABEL.get(c_by, f"|{inner}|")
    if c_by == "phi_mean":
        return (0.5 * (df["s1"].to_numpy(float) + df["s2"].to_numpy(float)),
                _C_LABEL["phi_mean"])
    if c_by == "phi_max":
        return df["s1"].to_numpy(float), _C_LABEL["phi_max"]
    if c_by == "phi_min":
        return df["s2"].to_numpy(float), _C_LABEL["phi_min"]
    col_name = c_by if c_by in df.columns else "skew"
    return df[col_name].to_numpy(float), _C_LABEL.get(c_by, c_by)


def _attach_cbar(
    ax: plt.Axes,
    mappable,
    label: str,
    *,
    size: str = "3%",
    pad: float = 0.06,
    labelsize: int = 9,
    ticksize: int = 8,
) -> None:
    """Attach a right-side colorbar to *ax* via make_axes_locatable."""
    div = make_axes_locatable(ax)
    cax = div.append_axes("right", size=size, pad=pad)
    cb = ax.get_figure().colorbar(mappable, cax=cax)
    cb.set_label(label, fontsize=labelsize)
    cb.ax.tick_params(labelsize=ticksize)


# ------------------------------ plotting --------------------------------- #

def plot_phase_tensor_psection(
    sites: Any,
    *,
    # ── data selection ────────────────────────────────────────────────────
    stations: Optional[List[str]] = None,
    period_range: Optional[Tuple[float, float]] = None,
    # ── y-axis ────────────────────────────────────────────────────────────
    axis_y: str = "logperiod",
    period_up: bool = True,
    # ── ellipse sizing ────────────────────────────────────────────────────
    scale          = _UNSET,   # default: PYCSAMT_STYLE.pt_ellipse.scale
    normalise_by   = _UNSET,   # default: PYCSAMT_STYLE.pt_ellipse.normalise_by
    s1_ref: Optional[float] = None,
    # ── colour ────────────────────────────────────────────────────────────
    c_by           = _UNSET,   # default: PYCSAMT_STYLE.pt_ellipse.c_by
    cmap           = _UNSET,   # default: auto from c_by via resolve_cmap()
    clim: Optional[Tuple[float, float]] = None,
    clim_pct       = _UNSET,   # default: PYCSAMT_STYLE.pt_ellipse.clim_pct
    symmetric_clim = _UNSET,   # default: PYCSAMT_STYLE.pt_ellipse.resolve_symmetric_clim()
    # ── aesthetics ────────────────────────────────────────────────────────
    edgecolor      = _UNSET,   # default: PYCSAMT_STYLE.pt_ellipse.edgecolor
    linewidth      = _UNSET,   # default: PYCSAMT_STYLE.pt_ellipse.linewidth
    alpha          = _UNSET,   # default: PYCSAMT_STYLE.pt_ellipse.alpha
    # ── overlays ─────────────────────────────────────────────────────────
    skew_threshold = _UNSET,   # default: PYCSAMT_STYLE.pt_ellipse.skew_threshold
    mark_3d        = _UNSET,   # default: PYCSAMT_STYLE.pt_ellipse.mark_3d
    ref_ellipse    = _UNSET,   # default: PYCSAMT_STYLE.pt_ellipse.show_ref
    # ── labels & layout ───────────────────────────────────────────────────
    title: str = "",
    xlabel: str = "",
    ylabel: str = "",
    tick_label_rotation: float = 45.0,
    figsize: Tuple[float, float] = (10.0, 5.5),
    # ── standard emtools ──────────────────────────────────────────────────
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """
    Phase-tensor ellipse pseudo-section (Caldwell et al. 2004 style).

    Each cell ``(station × period)`` is drawn as an ellipse whose shape,
    orientation, and fill colour encode phase-tensor invariants:

    * **major axis** ∝ φ_max  (s1 eigenvalue — mean phase level)
    * **minor axis** ∝ φ_min  (s2 eigenvalue — minimum phase)
    * **aspect ratio** = φ_min / φ_max (1 = isotropic / 1-D)
    * **rotation** = θ (geoelectric strike, CCW from E)
    * **fill colour** controlled by *c_by* (default: skew angle β)

    Parameters
    ----------
    sites : any
        EDI paths, glob pattern, :class:`~pycsamt.site.base.Sites`, or any
        input accepted by :func:`~pycsamt.emtools.ensure_sites`.
    stations : list of str or None
        Restrict to a subset of station names.
    period_range : (T_min, T_max) in seconds or None
        Restrict the period range plotted.
    axis_y : ``"logperiod"`` | ``"logfreq"``
        Y-axis quantity.  ``"logperiod"`` is the MT convention.
    period_up : bool, default ``True``
        When ``True`` (MT convention) long period (low frequency) is at the
        top of the figure; ``False`` flips the y-axis.
    scale : float, default ``0.85``
        Fraction of each cell occupied by the *reference* ellipse
        (the one with s1 = s1_ref).  Values above 1 cause overlap.
    normalise_by : ``"cell"`` | ``"unity"`` | ``"abs"``
        Sizing strategy:

        ``"cell"``
            Sizes are normalised so the 90th-percentile s1 fills *scale*
            of its cell.  Preserves relative size information.
        ``"unity"``
            s1_ref = 1.0 (tan 45° = 1-D half-space reference).  A 45°
            phase tensor is drawn as a circle filling *scale* of its cell.
        ``"abs"``
            Raw: width = scale × s1, height = scale × s2 in data units
            (legacy behaviour).
    s1_ref : float or None
        Override the reference s1 for ``"cell"`` and ``"unity"`` modes.
    c_by : str, default ``"skew"``
        Column name or derived quantity to map to fill colour.  Supported
        values: ``"skew"``, ``"beta"``, ``"alpha"``, ``"theta"``,
        ``"ellipt"``, ``"s1"``, ``"s2"``, ``"|skew|"``, ``"|beta|"``,
        ``"|theta|"``, ``"phi_mean"``, ``"phi_max"``, ``"phi_min"``.
    cmap : str, default ``"RdBu_r"``
        Matplotlib colormap name.
    clim : (vmin, vmax) or None
        Explicit colour limits.  When ``None``, derived from *clim_pct*.
    clim_pct : (lo, hi), default ``(5.0, 95.0)``
        Percentile limits used when *clim* is ``None``.
    symmetric_clim : bool, default ``True``
        Enforce vmin = −vmax for skew-like quantities.
    edgecolor : str, default ``"k"``
        Ellipse border colour.  Set to ``"none"`` to suppress borders.
    linewidth : float, default ``0.2``
        Ellipse border width (pts).  3-D cells receive 3 × this width when
        *mark_3d* is ``True``.
    alpha : float, default ``0.92``
        Ellipse fill opacity.
    skew_threshold : float or None, default ``3.0``
        |β| threshold (degrees) separating 1-D/2-D from 3-D structure.
        Used to cap *clim* symmetrically when *clim* is ``None``, and to
        highlight 3-D cells with a thicker border when *mark_3d* is ``True``.
    mark_3d : bool, default ``True``
        Draw a thicker border on ellipses where |β| > *skew_threshold*.
    ref_ellipse : bool, default ``True``
        Draw a labelled reference circle (φ_max = φ_min = s1_ref, β = 0°)
        in the upper-right corner as a scale indicator.
    title, xlabel, ylabel : str
        Axes title and axis labels.  Sensible defaults are used when empty.
    tick_label_rotation : float, default ``45.0``
        Station-name tick rotation in degrees.
    figsize : (width, height), default ``(10.0, 5.5)``
        Figure size (ignored when *ax* is provided).
    recursive, on_dup, strict, verbose : see :func:`ensure_sites`.
    ax : :class:`~matplotlib.axes.Axes` or None
        If provided, draw into this axes and return it; otherwise create a
        new figure.

    Returns
    -------
    ax : :class:`~matplotlib.axes.Axes`

    See Also
    --------
    build_phase_tensor_table : Underlying data computation.
    plot_phase_tensor_summary : 3-panel figure combining ellipses,
        dimensionality, and skew distribution.
    plot_phase_tensor_map : Geographic map view.
    """
    if isinstance(sites, pd.DataFrame):
        df = sites
    else:
        df = build_phase_tensor_table(
            sites,
            recursive=recursive,
            on_dup=on_dup,
            strict=strict,
            verbose=verbose,
        )

    # ── resolve visual style from PYCSAMT_STYLE.pt_ellipse ──────────────
    _es = PYCSAMT_STYLE.pt_ellipse
    if scale          is _UNSET: scale          = _es.scale
    if normalise_by   is _UNSET: normalise_by   = _es.normalise_by
    if c_by           is _UNSET: c_by           = _es.c_by
    if cmap           is _UNSET: cmap           = _es.copy(c_by=c_by).resolve_cmap()
    if clim_pct       is _UNSET: clim_pct       = _es.clim_pct
    if symmetric_clim is _UNSET: symmetric_clim = _es.copy(c_by=c_by).resolve_symmetric_clim()
    if edgecolor      is _UNSET: edgecolor      = _es.edgecolor
    if linewidth      is _UNSET: linewidth      = _es.linewidth
    if alpha          is _UNSET: alpha          = _es.alpha
    if skew_threshold is _UNSET: skew_threshold = _es.skew_threshold
    if mark_3d        is _UNSET: mark_3d        = _es.mark_3d
    if ref_ellipse    is _UNSET: ref_ellipse    = _es.show_ref

    # ── data selection ────────────────────────────────────────────────────
    if not df.empty and stations is not None:
        df = df[df["station"].isin(stations)]
    if not df.empty and period_range is not None:
        lo, hi = float(period_range[0]), float(period_range[1])
        df = df[(df["period"] >= lo) & (df["period"] <= hi)]

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    if df.empty:
        ax.text(0.5, 0.5, "no phase tensor data", ha="center", va="center",
                transform=ax.transAxes)
        return ax

    # ── station ordering & y-axis ─────────────────────────────────────────
    st_list = list(dict.fromkeys(df["station"].tolist()))   # preserve order
    st_map: Dict[str, int] = {s: i for i, s in enumerate(st_list)}
    n_st = len(st_list)

    if axis_y == "logperiod":
        df = df.copy()
        df["_y"] = np.log10(df["period"].to_numpy(float))
        ylab = ylabel or r"$\log_{10}(T)$ (s)"
    else:
        df = df.copy()
        df["_y"] = np.log10(np.maximum(df["freq"].to_numpy(float), 1e-30))
        ylab = ylabel or r"$\log_{10}(f)$ (Hz)"

    # ── cell dimensions in data coordinates ───────────────────────────────
    y_unique = np.unique(df["_y"].to_numpy(float))
    cell_y = float(np.nanmedian(np.abs(np.diff(y_unique)))) if len(y_unique) > 1 else 1.0
    cell_x = 1.0  # one station spacing

    # ── ellipse-size normalisation reference ──────────────────────────────
    s1_vals = df["s1"].to_numpy(float)
    if normalise_by == "cell":
        if s1_ref is not None:
            ref = float(s1_ref)
        else:
            finite_s1 = s1_vals[np.isfinite(s1_vals)]
            ref = float(np.nanpercentile(finite_s1, 90)) if len(finite_s1) else 1.0
        ref = max(ref, 1e-9)
    elif normalise_by == "unity":
        ref = float(s1_ref) if s1_ref is not None else 1.0
        ref = max(ref, 1e-9)
    else:  # "abs" — legacy: w = scale × s1 in raw data units
        ref = None

    # ── colour values & limits ────────────────────────────────────────────
    cvals, cbar_label = _resolve_cvals(df, c_by)
    finite_c = cvals[np.isfinite(cvals)]

    if clim is not None:
        vmin, vmax = float(clim[0]), float(clim[1])
    elif skew_threshold is not None and c_by in ("skew", "beta"):
        vmin, vmax = -float(skew_threshold), float(skew_threshold)
    else:
        lo_pct, hi_pct = clim_pct
        vmin = float(np.nanpercentile(finite_c, lo_pct)) if len(finite_c) else -1.0
        vmax = float(np.nanpercentile(finite_c, hi_pct)) if len(finite_c) else  1.0
        if symmetric_clim and c_by in _SYMMETRIC_C:
            vlim = max(abs(vmin), abs(vmax))
            vmin, vmax = -vlim, vlim

    norm = Normalize(vmin=vmin, vmax=vmax)
    cm   = plt.get_cmap(cmap)

    # ── extract arrays for the drawing loop ───────────────────────────────
    x_arr  = df["station"].map(st_map).to_numpy(float)
    y_arr  = df["_y"].to_numpy(float)
    s1_arr = df["s1"].to_numpy(float)
    s2_arr = df["s2"].to_numpy(float)
    th_arr = df["theta"].to_numpy(float)
    # skew values for mark_3d (always needed even if c_by != "skew")
    skew_arr = df["skew"].to_numpy(float) if "skew" in df.columns else np.zeros(len(df))

    # ── draw ellipses ─────────────────────────────────────────────────────
    for xi, yi, s1i, s2i, thi, ci, beti in zip(
            x_arr, y_arr, s1_arr, s2_arr, th_arr, cvals, skew_arr):
        if not np.all(np.isfinite([xi, yi, s1i, s2i, thi, ci])):
            continue

        if ref is not None:
            w = min(scale * cell_x * s1i / ref, 0.99 * cell_x)
            h = min(scale * cell_y * s2i / ref, 0.99 * cell_y)
        else:
            w = scale * s1i
            h = scale * s2i

        is_3d = (skew_threshold is not None and mark_3d
                 and np.isfinite(beti) and abs(beti) > skew_threshold)
        lw_i = linewidth * 3.0 if is_3d else linewidth

        ax.add_patch(Ellipse(
            (xi, yi),
            width=w, height=h, angle=thi,
            facecolor=cm(norm(ci)),
            edgecolor=edgecolor,
            linewidth=lw_i,
            alpha=alpha,
            zorder=3,
        ))

    # ── axes limits and ticks ─────────────────────────────────────────────
    ax.set_xlim(-0.6, n_st - 0.4)

    y_all  = y_arr[np.isfinite(y_arr)]
    y_lo   = float(np.nanmin(y_all)) if len(y_all) else 0.0
    y_hi   = float(np.nanmax(y_all)) if len(y_all) else 1.0
    margin = max(0.5 * cell_y, 0.02 * (y_hi - y_lo + 1e-9))

    if period_up:
        ax.set_ylim(y_lo - margin, y_hi + margin)  # long period at top ✓
    else:
        ax.set_ylim(y_hi + margin, y_lo - margin)  # short period at top

    # x-ticks: station names
    ax.set_xticks(np.arange(n_st))
    ax.set_xticklabels(st_list, rotation=tick_label_rotation,
                       ha="right", fontsize=8)

    # y-ticks: integer log10 values
    y_int = np.arange(int(np.floor(y_lo)), int(np.ceil(y_hi)) + 1)
    ax.set_yticks(y_int)
    ax.set_yticklabels([str(int(v)) for v in y_int], fontsize=8)

    ax.set_xlabel(xlabel or "Station", fontsize=9)
    ax.set_ylabel(ylab, fontsize=9)
    if title:
        ax.set_title(title, fontsize=10)
    ax.grid(True, ls=":", lw=0.4, color="#cccccc", alpha=0.7, zorder=0)

    # ── colorbar ─────────────────────────────────────────────────────────
    sm = ScalarMappable(cmap=cm, norm=norm)
    sm.set_array([])
    _attach_cbar(ax, sm, cbar_label)

    # ── optional overlays ─────────────────────────────────────────────────
    # Annotation boxes: 1-D/2-D vs 3-D text
    if skew_threshold is not None and c_by in ("skew", "beta"):
        ax.text(0.02, 0.98, f"|β| < {skew_threshold:.0f}° → 1-D/2-D",
                transform=ax.transAxes, fontsize=7.5, va="top",
                color="#2166ac",
                bbox=dict(fc="white", ec="none", alpha=0.75))
        ax.text(0.02, 0.91, f"|β| ≥ {skew_threshold:.0f}° → 3-D",
                transform=ax.transAxes, fontsize=7.5, va="top",
                color="#b2182b",
                bbox=dict(fc="white", ec="none", alpha=0.75))

    # Reference ellipse (circle with s1 = s1_ref = "45° 1-D response")
    if ref_ellipse and ref is not None:
        rx = n_st - 1 - 0.05 * (n_st - 1)   # near right edge
        ry = y_hi - 1.5 * cell_y             # near top
        rw = scale * cell_x                  # reference circle width
        rh = scale * cell_y                  # reference circle height
        ax.add_patch(Ellipse(
            (rx, ry), width=rw, height=rh, angle=0.0,
            fill=False, edgecolor="k", linewidth=0.8, zorder=5,
        ))
        ax.text(rx, ry - rh * 0.6,
                f"φ = s₁_ref", ha="center", va="top", fontsize=6.5)

    return ax


def plot_phase_tensor_skewmap(
    sites: Any,
    *,
    axis_y: str = "logperiod",
    agg: str = "median",
    figsize: Tuple[float, float] = (9.0, 4.0),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    df = build_phase_tensor_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if df.empty:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no phase tensor",
                ha="center", va="center")
        return ax
    if axis_y == "logperiod":
        df = df.copy()
        df["y"] = np.log10(df["period"].to_numpy())
        ylab = "LogPeriod (s)"
    else:
        df = df.copy()
        df["y"] = df["period"].to_numpy()
        ylab = "Period (s)"
    piv = df.pivot_table(
        index="y",
        columns="station",
        values="skew",
        aggfunc=agg,
    )
    piv = piv.sort_index()
    Z = piv.to_numpy(dtype=float)        # (n_logp, n_st) — no transpose
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(
        Z,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
    )
    ax.set_xlabel("Station")
    ax.set_ylabel(ylab)
    ax.set_xticks(np.arange(len(piv.columns)))
    ax.set_xticklabels(piv.columns, rotation=90)
    yt = np.linspace(0, Z.shape[0] - 1, num=min(8, Z.shape[0]))
    yvals = np.linspace(piv.index.min(), piv.index.max(),
                        num=min(8, len(piv.index)))
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{v:.3g}" for v in yvals])
    cb = plt.colorbar(im, ax=ax)
    cb.set_label("skew")
    return ax

# --- more plots ----------------------------------------------------------- #

def _coords_of_sites(sites) -> dict:
    m = {}
    for i, ed in enumerate(_iter_items(sites)):
        nm = _name(ed, i)
        lat = getattr(ed, "lat", None) or getattr(ed, "latitude", None)
        lon = getattr(ed, "lon", None) or getattr(ed, "longitude", None)
        try:
            if lat is None or lon is None:
                continue
            m[nm] = (float(lat), float(lon))
        except Exception:
            continue
    return m


def plot_theta_vs_period(
    sites: Any,
    *,
    figsize: Tuple[float, float] = (8.0, 4.0),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    df = build_phase_tensor_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    if df.empty:
        ax.text(0.5, 0.5, "no phase tensor",
                ha="center", va="center")
        return ax
    for st, sdf in df.groupby("station"):
        x = sdf["period"].to_numpy()
        y = sdf["theta"].to_numpy()
        ax.plot(x, y, ".", ms=3, label=st)
    ax.set_xscale("log")
    ax.set_xlabel("Period (s)")
    ax.set_ylabel("theta (deg)")
    ax.legend(ncols=2, fontsize=8)
    return ax


def plot_ellipticity_psection(
    sites: Any,
    *,
    figsize: Tuple[float, float] = (8.5, 4.0),
    agg: str = "median",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    df = build_phase_tensor_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    if df.empty:
        ax.text(0.5, 0.5, "no phase tensor",
                ha="center", va="center")
        return ax
    df = df.copy()
    df["logp"] = np.log10(df["period"].to_numpy())
    piv = df.pivot_table(
        index="logp",
        columns="station",
        values="ellipt",
        aggfunc=agg,
    )
    piv = piv.sort_index()
    Z = piv.to_numpy(dtype=float)        # (n_logp, n_st) — no transpose
    im = ax.imshow(
        Z,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
    )
    ax.set_xlabel("Station")
    ax.set_ylabel("LogPeriod (s)")
    ax.set_xticks(np.arange(Z.shape[1]))          # shape[1] = n_stations
    ax.set_xticklabels(piv.columns, rotation=90)
    yt = np.linspace(0, Z.shape[0] - 1, num=min(8, Z.shape[0]))  # shape[0] = n_logp
    yvals = np.linspace(piv.index.min(), piv.index.max(),
                        num=min(8, len(piv.index)))
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{v:.2g}" for v in yvals])
    cb = plt.colorbar(im, ax=ax)
    cb.set_label("ellipticity")
    return ax


def plot_dimensionality_psection(
    sites: Any,
    *,
    skew_th: float = 3.0,
    ellipt_th: float = 0.2,
    figsize: Tuple[float, float] = (8.5, 4.0),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    df = build_phase_tensor_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    if df.empty:
        ax.text(0.5, 0.5, "no phase tensor",
                ha="center", va="center")
        return ax
    df = df.copy()
    df["logp"] = np.log10(df["period"].to_numpy())
    # 0=1D, 1=2D, 2=3D (simple rule)
    lab = np.zeros(len(df), dtype=int)
    lab[(np.abs(df["skew"]) <= skew_th) &
        (df["ellipt"].abs() <= ellipt_th)] = 0
    lab[(np.abs(df["skew"]) <= skew_th) &
        (df["ellipt"].abs() > ellipt_th)] = 1
    lab[(np.abs(df["skew"]) > skew_th)] = 2
    df["dim"] = lab
    piv = df.pivot_table(
        index="logp",
        columns="station",
        values="dim",
        aggfunc="median",
    )
    piv = piv.sort_index()
    Z = piv.to_numpy(dtype=float)        # (n_logp, n_st) — no transpose
    im = ax.imshow(
        Z,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
    )
    ax.set_xlabel("Station")
    ax.set_ylabel("LogPeriod (s)")
    ax.set_xticks(np.arange(Z.shape[1]))          # shape[1] = n_stations
    ax.set_xticklabels(piv.columns, rotation=90)
    yt = np.linspace(0, Z.shape[0] - 1, num=min(8, Z.shape[0]))  # shape[0] = n_logp
    yvals = np.linspace(piv.index.min(), piv.index.max(),
                        num=min(8, len(piv.index)))
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{v:.2g}" for v in yvals])
    cb = plt.colorbar(im, ax=ax)
    cb.set_label("dim (0=1D,1=2D,2=3D)")
    return ax


def plot_phase_tensor_rose(
    sites: Any,
    *,
    # ── visual style ─────────────────────────────────────────────────────
    style: "str | RoseStyle | None" = "pycsamt",
    # ── data selection ───────────────────────────────────────────────────
    band: Optional[Tuple[float, float]] = None,
    freq_bands: Optional[List[Tuple[float, float]]] = None,
    band_labels: Optional[List[str]] = None,
    band_colors: Optional[List] = None,
    # ── histogram ────────────────────────────────────────────────────────
    bins: int = 36,
    # ── visual overrides (sentinel → taken from *style*) ─────────────────
    bar_style=_UNSET,
    bar_color=_UNSET,
    bar_alpha=_UNSET,
    bar_edgecolor=_UNSET,
    bar_edgelw=_UNSET,
    cmap=_UNSET,
    outer_ring_lw=_UNSET,
    outer_ring_color=_UNSET,
    n_rings=_UNSET,
    ring_color=_UNSET,
    ring_ls=_UNSET,
    ring_lw=_UNSET,
    ring_labels=_UNSET,
    ring_label_angle=_UNSET,
    ring_label_fontsize=_UNSET,
    ring_label_color=_UNSET,
    ring_label_fmt=_UNSET,
    spoke_every=_UNSET,
    spoke_color=_UNSET,
    spoke_ls=_UNSET,
    spoke_lw=_UNSET,
    compass_labels=_UNSET,
    compass_fontsize=_UNSET,
    compass_color=_UNSET,
    compass_fontweight=_UNSET,
    show_mean=_UNSET,
    mean_color=_UNSET,
    mean_lw=_UNSET,
    mean_ls=_UNSET,
    show_secondary=_UNSET,
    secondary_ls=_UNSET,
    secondary_lw=_UNSET,
    secondary_color=_UNSET,
    show_annotation=_UNSET,
    annotation_pos=_UNSET,
    annotation_fontsize=_UNSET,
    annotation_bg=_UNSET,
    annotation_ec=_UNSET,
    show_n=_UNSET,
    # ── layout ───────────────────────────────────────────────────────────
    figsize: Tuple[float, float] = (5.5, 5.5),
    title: str = "",
    title_fontsize: float = 10.0,
    # ── core ─────────────────────────────────────────────────────────────
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """Publication-quality phase-tensor θ rose diagram.

    Bars are drawn with axial symmetry — 0–180° mirrored to 180–360° —
    reflecting the inherent 180° ambiguity of the phase-tensor principal
    axis direction.

    Parameters
    ----------
    sites : any
        EDI paths, objects, or collection accepted by
        :func:`~pycsamt.emtools._core.ensure_sites`.
    band : (lo_s, hi_s) or None
        Period window in seconds.  ``None`` uses all available periods.
    freq_bands : list of (lo_s, hi_s), optional
        Sub-bands for ``bar_style="bands"`` — one stacked bar colour per
        band.
    band_labels : list[str], optional
        Legend labels for each entry in *freq_bands*.
    band_colors : list, optional
        Colours for each *freq_bands* entry.  Defaults to ``tab10``.
    bins : int
        Number of bins over 0–180°, mirrored to 360°.  Default ``36``
        gives 5° bins.
    bar_style : {"gradient", "solid", "bands"}
        ``"gradient"`` colours bars by height via *cmap*;
        ``"solid"`` uses *bar_color* uniformly;
        ``"bands"`` stacks one colour per period sub-band.
    bar_color : str
        Bar fill colour for ``bar_style="solid"``.
    bar_edgecolor : str
        Bar edge colour (``"none"`` → no edge).
    bar_edgelw : float
        Bar edge line-width.
    bar_alpha : float
        Bar opacity (0–1).
    cmap : str
        Colormap for ``bar_style="gradient"``.
    outer_ring_lw : float
        Line-width of the bold outer bounding circle.
    outer_ring_color : str
        Colour of the outer circle.
    n_rings : int
        Number of concentric reference rings.
    ring_color : str
        Colour of concentric rings and spokes.
    ring_ls : str
        Line-style of concentric rings.
    ring_lw : float
        Line-width of concentric rings.
    ring_labels : list[float], optional
        Explicit count values to annotate on the rings (e.g.
        ``[25, 50, 75, 100]``).  ``None`` → evenly spaced from
        ``rmax / n_rings`` to ``rmax``.
    ring_label_angle : float
        Clockwise angle from North (degrees) at which ring count labels
        are placed.  Default ``22.5``.
    ring_label_fontsize : float
        Font size for ring count labels.
    ring_label_color : str
        Colour for ring count labels.
    ring_label_fmt : str
        Format string for ring labels, e.g. ``"{:.0f}"``.
    spoke_every : float
        Angular spacing (degrees) between radial spokes.
    spoke_color : str
        Colour of radial spokes.
    spoke_ls : str
        Line-style of radial spokes.
    spoke_lw : float
        Line-width of radial spokes.
    compass_labels : {"NESW", "degrees", "none"}
        Perimeter labels.  ``"NESW"`` shows cardinal directions;
        ``"degrees"`` shows degree values; ``"none"`` hides all.
    compass_fontsize : float
        Font size for perimeter labels.
    compass_color : str
        Colour for perimeter labels.
    compass_fontweight : str
        Font weight for perimeter labels.
    show_mean : bool
        Draw the axial mean direction as a line through the centre.
    mean_color : str
        Colour of the mean-direction line.
    mean_lw : float
        Line-width of the mean-direction line.
    mean_ls : str
        Line-style of the mean-direction line.
    show_secondary : bool
        Draw the 180°-conjugate mean line.
    secondary_ls : str
        Line-style of the conjugate line.
    secondary_lw : float, optional
        Line-width of the conjugate line; defaults to *mean_lw*.
    secondary_color : str, optional
        Colour of the conjugate line; defaults to *mean_color*.
    show_annotation : bool
        Show a text box with the mean θ and station/pair count.
    annotation_pos : (float, float)
        Axes-fraction ``(x, y)`` for the annotation box.
    annotation_fontsize : float
        Font size of the annotation text.
    annotation_bg : str
        Background colour of the annotation box.
    annotation_ec : str
        Edge colour of the annotation box.
    show_n : bool
        Append ``n = N`` to the annotation text.
    figsize : (float, float)
        Figure size in inches.
    title : str
        Axes title (set via ``ax.set_title``).
    title_fontsize : float
        Font size for *title*.
    recursive, on_dup, strict, verbose
        Passed to :func:`~pycsamt.emtools._core.ensure_sites`.
    ax : matplotlib.axes.Axes, optional
        Pre-existing polar axes to draw into.  Created when ``None``.

    Returns
    -------
    matplotlib.axes.Axes
        The polar axes containing the rose diagram.

    Examples
    --------
    Default gradient style, all periods:

    >>> from pycsamt.emtools import plot_phase_tensor_rose
    >>> ax = plot_phase_tensor_rose("path/to/edis/", figsize=(6, 6))

    Frequency-band decomposition (stacked):

    >>> ax = plot_phase_tensor_rose(
    ...     sites,
    ...     bar_style="bands",
    ...     freq_bands=[(1e-4, 1e-2), (1e-2, 1e0)],
    ...     band_labels=["Short period", "Long period"],
    ... )

    Custom ring count labels:

    >>> ax = plot_phase_tensor_rose(
    ...     sites,
    ...     ring_labels=[25, 50, 75, 100],
    ...     ring_label_angle=15.0,
    ... )
    """
    # ── resolve style → fill _UNSET visual params ─────────────────────────
    rs = resolve_rose_style(style)
    def _v(val, attr): return getattr(rs, attr) if val is _UNSET else val
    bar_style        = _v(bar_style,        "bar_style")
    bar_color        = _v(bar_color,        "bar_color")
    bar_alpha        = _v(bar_alpha,        "bar_alpha")
    bar_edgecolor    = _v(bar_edgecolor,    "bar_edgecolor")
    bar_edgelw       = _v(bar_edgelw,       "bar_edgelw")
    cmap             = _v(cmap,             "cmap")
    outer_ring_lw    = _v(outer_ring_lw,    "outer_ring_lw")
    outer_ring_color = _v(outer_ring_color, "outer_ring_color")
    n_rings          = _v(n_rings,          "n_rings")
    ring_color       = _v(ring_color,       "ring_color")
    ring_ls          = _v(ring_ls,          "ring_ls")
    ring_lw          = _v(ring_lw,          "ring_lw")
    ring_labels      = _v(ring_labels,      "ring_labels")
    ring_label_angle = _v(ring_label_angle, "ring_label_angle")
    ring_label_fontsize = _v(ring_label_fontsize, "ring_label_fontsize")
    ring_label_color = _v(ring_label_color, "ring_label_color")
    ring_label_fmt   = _v(ring_label_fmt,   "ring_label_fmt")
    spoke_every      = _v(spoke_every,      "spoke_every")
    spoke_color      = _v(spoke_color,      "spoke_color")
    spoke_ls         = _v(spoke_ls,         "spoke_ls")
    spoke_lw         = _v(spoke_lw,         "spoke_lw")
    compass_labels   = _v(compass_labels,   "compass_labels")
    compass_fontsize = _v(compass_fontsize, "compass_fontsize")
    compass_color    = _v(compass_color,    "compass_color")
    compass_fontweight = _v(compass_fontweight, "compass_fontweight")
    show_mean        = _v(show_mean,        "show_mean")
    mean_color       = _v(mean_color,       "mean_color")
    mean_lw          = _v(mean_lw,          "mean_lw")
    mean_ls          = _v(mean_ls,          "mean_ls")
    show_secondary   = _v(show_secondary,   "show_secondary")
    secondary_color  = _v(secondary_color,  "secondary_color")
    secondary_ls     = _v(secondary_ls,     "secondary_ls")
    secondary_lw     = _v(secondary_lw,     "secondary_lw")
    show_annotation  = _v(show_annotation,  "show_annotation")
    annotation_pos   = _v(annotation_pos,   "annotation_pos")
    annotation_fontsize = _v(annotation_fontsize, "annotation_fontsize")
    annotation_bg    = _v(annotation_bg,    "annotation_bg")
    annotation_ec    = _v(annotation_ec,    "annotation_ec")
    show_n           = _v(show_n,           "show_n")

    # ── axial mean helper (local, avoids import cycle with strike.py) ──────
    def _axial_mean(a_deg: np.ndarray) -> float:
        """Circular mean for axial (0–180°) data, returned in [0,180°)."""
        rad = np.radians(2.0 * a_deg)
        mu = np.degrees(np.arctan2(np.nanmean(np.sin(rad)),
                                   np.nanmean(np.cos(rad)))) / 2.0
        return float(mu % 180.0)

    # ── build phase-tensor table ───────────────────────────────────────────
    df = build_phase_tensor_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    if ax is None:
        _, ax = plt.subplots(subplot_kw=dict(polar=True), figsize=figsize)

    if df.empty:
        ax.text(0.5, 0.5, "no phase tensor data",
                ha="center", va="center", transform=ax.transAxes)
        return ax

    # ── period filtering ───────────────────────────────────────────────────
    if band is not None:
        lo, hi = float(band[0]), float(band[1])
        sel = (df["period"] >= lo) & (df["period"] <= hi)
        band_label = f"{lo:.4g}–{hi:.4g} s"
    else:
        sel = np.ones(len(df), dtype=bool)
        p = df["period"]
        band_label = f"{p.min():.4g}–{p.max():.4g} s"

    # ── histogram over 0–180° ─────────────────────────────────────────────
    bins_ = int(max(12, bins))
    edges_deg = np.linspace(0.0, 180.0, bins_ + 1)
    dw = np.radians(180.0 / bins_)

    use_bands = bar_style == "bands" and bool(freq_bands)

    if use_bands:
        n_fb = len(freq_bands)  # type: ignore[arg-type]
        _bc: List = (list(band_colors) if band_colors is not None
                     else list(plt.get_cmap("tab10")(np.linspace(0, 0.8, n_fb))))
        _bl: List[str] = (list(band_labels) if band_labels is not None
                          else [f"{lo_:.4g}–{hi_:.4g} s"
                                for lo_, hi_ in freq_bands])  # type: ignore
        hists: List[np.ndarray] = []
        for fb in freq_bands:  # type: ignore[union-attr]
            lo_, hi_ = float(fb[0]), float(fb[1])
            sel_fb = (df["period"] >= lo_) & (df["period"] <= hi_)
            th_fb = df.loc[sel_fb, "theta"].to_numpy(float) % 180.0
            h, _ = np.histogram(th_fb, bins=edges_deg)
            hists.append(h)
        hist_total = np.sum(hists, axis=0)
    else:
        th_raw = df.loc[sel, "theta"].to_numpy(float) % 180.0
        hist_total, _ = np.histogram(th_raw, bins=edges_deg)

    rmax = float(hist_total.max()) if hist_total.max() > 0 else 1.0
    ang_mid = np.radians(0.5 * (edges_deg[1:] + edges_deg[:-1]))

    # ── polar setup ───────────────────────────────────────────────────────
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    # hide default radial ticks/labels; we draw our own
    ax.set_yticklabels([])
    ax.set_yticks([])
    ax.yaxis.grid(False)
    ax.xaxis.grid(False)
    ax.set_frame_on(False)

    # ── draw bars (0–180° and mirror 180–360°) ────────────────────────────
    def _draw_bars(ang_c: np.ndarray, heights: np.ndarray, color,
                   bottom: float = 0.0) -> None:
        for a, h in zip(ang_c, heights):
            for a_plot in (a, a + np.pi):      # axial symmetry
                ax.bar(
                    a_plot, h, width=dw,
                    bottom=bottom,
                    color=color,
                    edgecolor=bar_edgecolor,
                    linewidth=bar_edgelw,
                    alpha=bar_alpha,
                )

    if use_bands:
        bottoms = np.zeros(bins_)
        for h, col in zip(hists, _bc):
            _draw_bars(ang_mid, h, col, bottom=0.0)
            # stacked: accumulate actual bottom for next band
            bottoms += h
        # legend
        from matplotlib.patches import Patch as _Patch
        ax.legend(
            handles=[_Patch(fc=c, label=l, alpha=bar_alpha)
                     for c, l in zip(_bc, _bl)],
            loc="lower left",
            bbox_to_anchor=(1.05, 0.0),
            fontsize=max(6, annotation_fontsize - 1),
            framealpha=0.9,
        )
    elif bar_style == "gradient":
        cm_ = plt.get_cmap(cmap)
        for a, h in zip(ang_mid, hist_total):
            col = cm_(h / rmax) if rmax > 0 else cm_(0.5)
            for a_plot in (a, a + np.pi):
                ax.bar(
                    a_plot, h, width=dw,
                    color=col,
                    edgecolor=bar_edgecolor,
                    linewidth=bar_edgelw,
                    alpha=bar_alpha,
                )
    else:  # solid
        _draw_bars(ang_mid, hist_total, bar_color)

    # ── concentric rings ──────────────────────────────────────────────────
    if ring_labels is not None:
        r_levels = [float(v) for v in ring_labels]
    else:
        step = rmax / n_rings
        r_levels = [step * k for k in range(1, n_rings + 1)]

    theta_full = np.linspace(0, 2 * np.pi, 360)
    for rv in r_levels:
        ax.plot(theta_full, np.full_like(theta_full, rv),
                color=ring_color, ls=ring_ls, lw=ring_lw, zorder=0)

    # ring count annotations
    lbl_theta = np.radians(ring_label_angle)
    for rv in r_levels:
        ax.text(
            lbl_theta, rv,
            ring_label_fmt.format(rv),
            ha="center", va="center",
            fontsize=ring_label_fontsize,
            color=ring_label_color,
            zorder=5,
        )

    # ── radial spokes ─────────────────────────────────────────────────────
    n_spokes = int(round(360.0 / spoke_every))
    spoke_angles = np.radians(np.arange(n_spokes) * spoke_every)
    for sa in spoke_angles:
        ax.plot([sa, sa], [0, rmax],
                color=spoke_color, ls=spoke_ls, lw=spoke_lw, zorder=0)

    # ── outer ring ────────────────────────────────────────────────────────
    ax.spines["polar"].set_linewidth(outer_ring_lw)
    ax.spines["polar"].set_edgecolor(outer_ring_color)
    ax.set_ylim(0, rmax * 1.08)

    # ── compass labels ────────────────────────────────────────────────────
    n_spk = int(round(360.0 / spoke_every))
    spoke_degs = np.arange(n_spk) * spoke_every
    if compass_labels == "NESW":
        _card = {0: "N", 90: "E", 180: "S", 270: "W"}
        lbl_list = [_card.get(int(d) % 360, "") for d in spoke_degs]
    elif compass_labels == "degrees":
        lbl_list = [f"{int(d)}°" for d in spoke_degs]
    else:
        lbl_list = [""] * n_spk

    ax.set_thetagrids(
        spoke_degs, labels=lbl_list,
        fontsize=compass_fontsize,
    )
    for lbl in ax.get_xticklabels():
        lbl.set_color(compass_color)
        lbl.set_fontweight(compass_fontweight)

    # ── mean direction ────────────────────────────────────────────────────
    th_all = df.loc[sel, "theta"].to_numpy(float) % 180.0
    mu = _axial_mean(th_all) if len(th_all) > 0 else 0.0
    sec_lw  = secondary_lw  if secondary_lw  is not None else mean_lw
    sec_col = secondary_color if secondary_color is not None else mean_color

    if show_mean and len(th_all) > 0:
        mu_r = np.radians(mu)
        # draw full diameter: two radii (center→rim) in opposite directions
        for ang_pt in (mu_r, mu_r + np.pi):
            ax.plot([ang_pt, ang_pt], [0, rmax],
                    color=mean_color, lw=mean_lw, ls=mean_ls,
                    solid_capstyle="round", zorder=10)

        if show_secondary:
            mu90 = np.radians(mu + 90.0)
            for ang_pt in (mu90, mu90 + np.pi):
                ax.plot([ang_pt, ang_pt], [0, rmax],
                        color=sec_col, lw=sec_lw, ls=secondary_ls,
                        solid_capstyle="round", zorder=9)

    # ── annotation box ────────────────────────────────────────────────────
    if show_annotation and len(th_all) > 0:
        txt = f"θ̄ = {mu:.1f}°"
        if show_n:
            txt += f"\nn = {len(th_all)}"
        ax.text(
            annotation_pos[0], annotation_pos[1], txt,
            transform=ax.transAxes,
            fontsize=annotation_fontsize,
            va="top", ha="left",
            bbox=dict(
                boxstyle="round,pad=0.3",
                fc=annotation_bg,
                ec=annotation_ec,
                lw=0.8,
            ),
            zorder=20,
        )

    # ── title ─────────────────────────────────────────────────────────────
    if title:
        ax.set_title(title, fontsize=title_fontsize, pad=14)
    else:
        ax.set_title(f"Phase-tensor θ rose  ({band_label})",
                     fontsize=title_fontsize, pad=14)

    return ax

def plot_phase_tensor_map(
    sites: Any,
    *,
    period: float = 10.0,
    # ── ellipse style (honours PhaseTensorEllipseStyle) ───────────────
    scale          = _UNSET,
    normalise_by   = _UNSET,
    s1_ref         = _UNSET,
    c_by           = _UNSET,
    cmap           = _UNSET,
    clim: Optional[Tuple[float, float]]            = None,
    clim_pct       = _UNSET,
    symmetric_clim = _UNSET,
    alpha          = _UNSET,
    edgecolor      = _UNSET,
    linewidth      = _UNSET,
    skew_threshold = _UNSET,
    mark_3d        = _UNSET,
    lw_3d_factor   = _UNSET,
    ref_ellipse    = _UNSET,
    # ── tipper arrows ─────────────────────────────────────────────────
    show_tipper: bool = True,
    tipper_convention: str = "parkinson",
    tipper_component: str = "real",
    tipper_scale: Optional[float] = None,
    tipper_color: str = "k",
    tipper_lw: float = 1.4,
    # ── optional background grid ──────────────────────────────────────
    bg_grid: Optional[Dict[str, Any]] = None,
    # ── map decoration ────────────────────────────────────────────────
    station_labels: bool = True,
    station_marker: str = "v",
    station_ms: float = 5.0,
    station_color: str = "k",
    label_fontsize: float = 7.0,
    title: str = "",
    colorbar_label: Optional[str] = None,
    # ── coordinate override (station → (lat, lon)) ────────────────────
    coords: Optional[Dict[str, Tuple[float, float]]] = None,
    # ── layout ────────────────────────────────────────────────────────
    figsize: Tuple[float, float] = (9.0, 7.0),
    # ── core ──────────────────────────────────────────────────────────
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """Publication-quality phase-tensor map at a single period.

    Renders each station as a phase-tensor ellipse positioned at its
    geographic coordinates.  The ellipse shape encodes the phase-tensor
    principal values (φ_max, φ_min) and orientation (θ); fill colour
    encodes an additional scalar (default: skewness β).  Induction arrows
    can be overlaid when tipper data are available, and an optional
    background field (gravity, resistivity, …) may be drawn behind the
    ellipses.

    Parameters
    ----------
    sites : any
        EDI paths, objects, or collection accepted by
        :func:`~pycsamt.emtools._core.ensure_sites`.
    period : float
        Target period (s).  The nearest available period in the data is
        used for each station independently.
    scale : float or _UNSET
        Maximum ellipse semi-axis in geographic units (degrees).
        ``_UNSET`` → auto-derived from the median inter-station spacing.
    normalise_by : ``"cell"`` | ``"unity"`` | ``"abs"``
        Size normalisation strategy, see
        :class:`~pycsamt.api.style.PhaseTensorEllipseStyle`.
    s1_ref : float or None
        Manual reference s1 used with ``normalise_by="cell"``/``"unity"``.
    c_by : str
        Scalar quantity to map to fill colour.  Recognised values:
        ``"skew"``, ``"beta"``, ``"|skew|"``, ``"theta"``,
        ``"ellipt"``, ``"phi_mean"``, ``"phi_max"``, ``"phi_min"``,
        ``"s1"``, ``"s2"``, ``"alpha"``.  Defaults to
        ``PYCSAMT_STYLE.pt_ellipse.c_by``.
    cmap : str or None
        Matplotlib colormap.  Auto-selected from *c_by* when ``None``
        (same logic as :class:`~pycsamt.api.style.PhaseTensorEllipseStyle`).
    clim : (vmin, vmax) or None
        Explicit colour limits.  ``None`` → derived from *clim_pct*.
    clim_pct : (lo, hi)
        Percentile window for automatic colour limits.
    symmetric_clim : bool
        Force vmin = −vmax (useful for diverging quantities).
    alpha : float
        Ellipse fill opacity.
    edgecolor : str
        Ellipse border colour.
    linewidth : float
        Normal ellipse border width (pts).
    skew_threshold : float or None
        |β| above which a cell is flagged as 3-D.
    mark_3d : bool
        Draw thicker borders on 3-D flagged ellipses.
    lw_3d_factor : float
        Line-width multiplier for 3-D flagged borders.
    ref_ellipse : bool
        Draw a reference circle in the lower-left corner as a scale bar.
    show_tipper : bool
        Overlay induction arrows when tipper data are found.
    tipper_convention : ``"parkinson"`` | ``"wiese"``
        Parkinson: arrow toward anomaly (negated real T);
        Wiese: arrow along real T.
    tipper_component : ``"real"`` | ``"imag"`` | ``"both"``
        Which tipper component to draw.
    tipper_scale : float or None
        Arrow length = tipper_scale × *scale*.  Auto-derived when ``None``.
    tipper_color : str
        Arrow colour (real component).  Imaginary component uses
        ``mcolors.to_rgba(tipper_color, 0.55)``.
    tipper_lw : float
        Arrow line-width.
    bg_grid : dict or None
        Optional background field drawn behind ellipses::

            bg_grid = dict(
                lons   = 1D or 2D longitude array,
                lats   = 1D or 2D latitude array,
                values = 2D array (shape matches meshgrid of lons × lats),
                cmap   = "RdYlGn",     # colormap
                clim   = (vmin, vmax), # or None → auto
                alpha  = 0.55,         # opacity
                label  = "Gravity (gu)",
            )

    station_labels : bool
        Annotate each station position with its name.
    station_marker : str
        Marker style for station positions (default ``"v"``).
    station_ms : float
        Marker size.
    station_color : str
        Marker and label colour.
    label_fontsize : float
        Font size for station labels.
    title : str
        Axes title.
    colorbar_label : str or None
        Override the automatic colorbar label derived from *c_by*.
    coords : dict[str, (lat, lon)] or None
        Explicit station coordinates.  When ``None`` the function reads
        ``ed.coords`` (a ``(lat, lon, elev)`` tuple) from each Site object.
    figsize : (float, float)
        Figure size in inches.
    recursive, on_dup, strict, verbose
        Passed to :func:`~pycsamt.emtools._core.ensure_sites`.
    ax : matplotlib.axes.Axes or None
        Pre-existing axes to draw into.

    Returns
    -------
    matplotlib.axes.Axes

    Examples
    --------
    Default skew-coloured map:

    >>> from pycsamt.emtools import plot_phase_tensor_map
    >>> ax = plot_phase_tensor_map(sites, period=10.0)

    Ellipticity coloured, with tipper arrows:

    >>> ax = plot_phase_tensor_map(
    ...     sites, period=10.0,
    ...     c_by="ellipt", cmap="viridis",
    ...     show_tipper=True, tipper_convention="parkinson",
    ... )

    With a gravity background:

    >>> ax = plot_phase_tensor_map(
    ...     sites, period=10.0,
    ...     bg_grid=dict(lons=g_lon, lats=g_lat, values=g_bouguer,
    ...                  cmap="RdYlGn", alpha=0.45, label="Gravity (gu)"),
    ... )
    """
    from ..api.plot import add_colorbar
    import matplotlib.colors as _mc

    # ── 1. resolve style ──────────────────────────────────────────────────
    _es = PYCSAMT_STYLE.pt_ellipse
    def _sv(v, attr): return getattr(_es, attr) if v is _UNSET else v
    normalise_by   = _sv(normalise_by,   "normalise_by")
    s1_ref_        = _sv(s1_ref,         "s1_ref")
    c_by           = _sv(c_by,           "c_by")
    cmap_          = _sv(cmap,           "cmap")
    clim_pct       = _sv(clim_pct,       "clim_pct")
    symmetric_clim = _sv(symmetric_clim, "symmetric_clim")
    alpha_         = _sv(alpha,          "alpha")
    edgecolor_     = _sv(edgecolor,      "edgecolor")
    linewidth_     = _sv(linewidth,      "linewidth")
    skew_threshold_= _sv(skew_threshold, "skew_threshold")
    mark_3d_       = _sv(mark_3d,        "mark_3d")
    lw_3d_factor_  = _sv(lw_3d_factor,   "lw_3d_factor")
    ref_ellipse_   = _sv(ref_ellipse,     "show_ref")

    # auto-resolve cmap from c_by when None
    if cmap_ is None:
        cmap_ = _es.copy(c_by=c_by).resolve_cmap()
    # honour natural symmetric_clim for diverging quantities
    sym = symmetric_clim and _es.copy(c_by=c_by).resolve_symmetric_clim()

    # ── 2. load data ──────────────────────────────────────────────────────
    S = ensure_sites(sites, recursive=recursive, on_dup=on_dup,
                     strict=strict, verbose=verbose)
    df = build_phase_tensor_table(
        S, recursive=False, on_dup=on_dup, strict=False, verbose=verbose,
    )
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    if df.empty:
        ax.text(0.5, 0.5, "no phase tensor data",
                ha="center", va="center", transform=ax.transAxes)
        return ax

    # ── 3. station coordinates ────────────────────────────────────────────
    if coords is None:
        coords = {}
        for i, ed in enumerate(_iter_items(S)):
            st = _name(ed, i)
            # preferred: ed.coords → (lat, lon, elev)
            c = getattr(ed, "coords", None)
            if c is not None and len(c) >= 2 and c[0] is not None and c[1] is not None:
                coords[st] = (float(c[0]), float(c[1]))
                continue
            # fallback: ed.meta
            meta = getattr(ed, "meta", {}) or {}
            lat = meta.get("lat") or meta.get("latitude")
            lon = meta.get("lon") or meta.get("longitude") or meta.get("long")
            if lat is not None and lon is not None:
                coords[st] = (float(lat), float(lon))

    if not coords:
        ax.text(0.5, 0.5, "no geographic coordinates in EDI data\n"
                "pass coords={station: (lat, lon)} to override",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=9, style="italic", color="0.45")
        return ax

    # ── 4. select nearest period per station ─────────────────────────────
    rows: List[Any] = []
    for st, sdf in df.groupby("station"):
        if st not in coords:
            continue
        p = sdf["period"].to_numpy()
        idx = int(np.nanargmin(np.abs(p - period)))
        rows.append(sdf.iloc[idx])

    if not rows:
        ax.text(0.5, 0.5, f"no data near T = {period:g} s",
                ha="center", va="center", transform=ax.transAxes)
        return ax

    st_names = [r["station"] for r in rows]
    lats = np.array([coords[s][0] for s in st_names], float)
    lons = np.array([coords[s][1] for s in st_names], float)
    s1s  = np.array([float(r["s1"]) for r in rows], float)
    s2s  = np.array([float(r["s2"]) for r in rows], float)
    ths  = np.array([float(r["theta"]) for r in rows], float)
    betas = np.array([float(r.get("skew", r.get("beta", np.nan))) for r in rows], float)

    # ── 5. auto scale from grid spacing ──────────────────────────────────
    if scale is _UNSET:
        if len(lons) > 1:
            # Use unique grid spacings (cross-profile > within-profile)
            u_lons = np.unique(np.round(lons, 3))
            u_lats = np.unique(np.round(lats, 3))
            d_lon = (float(np.median(np.diff(u_lons)))
                     if len(u_lons) > 1 else float(np.ptp(lons) or 0.01))
            d_lat = (float(np.median(np.diff(u_lats)))
                     if len(u_lats) > 1 else float(np.ptp(lats) or 0.01))
            scale_ = 0.55 * max(abs(d_lon) + 1e-9, abs(d_lat) + 1e-9)
        else:
            scale_ = 0.01
    else:
        scale_ = float(scale)

    # ── 6. ellipse size normalisation ────────────────────────────────────
    # Clip extreme s1 outliers (robust 95th pct) before normalising
    s1s_clip = np.clip(s1s, 0.0, float(np.nanpercentile(s1s, 95)))
    if normalise_by == "unity":
        _s1ref = s1_ref_ if s1_ref_ is not None else float(np.tan(np.radians(45.0)))
    elif normalise_by == "cell":
        _s1ref = (s1_ref_ if s1_ref_ is not None
                  else float(np.nanpercentile(s1s_clip, 90)))
    else:  # "abs"
        _s1ref = s1_ref_ if s1_ref_ is not None else 1.0
    _s1ref = max(_s1ref, 1e-9)

    # ── 7. colour mapping ─────────────────────────────────────────────────
    c_vals = np.array([float(r.get(c_by, np.nan)) for r in rows], float)
    if clim is not None:
        v0, v1 = float(clim[0]), float(clim[1])
    else:
        lo_pct, hi_pct = float(clim_pct[0]), float(clim_pct[1])
        v0 = float(np.nanpercentile(c_vals, lo_pct))
        v1 = float(np.nanpercentile(c_vals, hi_pct))
    if sym and v0 != v1:
        vmax = max(abs(v0), abs(v1))
        v0, v1 = -vmax, vmax
    v0 = v0 if v0 != v1 else v0 - 1.0
    cm_obj = plt.get_cmap(cmap_)
    norm_c = plt.Normalize(vmin=v0, vmax=v1)

    # ── 8. optional background grid ───────────────────────────────────────
    bg_mappable = None
    if bg_grid is not None:
        g_lon = np.asarray(bg_grid.get("lons", []), float)
        g_lat = np.asarray(bg_grid.get("lats", []), float)
        g_val = np.asarray(bg_grid.get("values", []), float)
        g_cm  = bg_grid.get("cmap", "RdYlGn")
        g_cl  = bg_grid.get("clim", None)
        g_al  = float(bg_grid.get("alpha", 0.5))
        g_lbl = bg_grid.get("label", "")
        if g_val.size > 0:
            if g_val.ndim == 1:
                g_val = g_val.reshape(len(g_lat), len(g_lon))
            gv0 = float(np.nanpercentile(g_val, 5)) if g_cl is None else float(g_cl[0])
            gv1 = float(np.nanpercentile(g_val, 95)) if g_cl is None else float(g_cl[1])
            bg_mappable = ax.pcolormesh(
                g_lon, g_lat, g_val,
                cmap=g_cm, vmin=gv0, vmax=gv1,
                alpha=g_al, shading="auto", zorder=0,
            )
            if g_lbl:
                add_colorbar(bg_mappable, ax, label=g_lbl,
                             side="left", size="3.5%", pad=0.08)

    # ── 9. draw ellipses ──────────────────────────────────────────────────
    for row, lat, lon, s1, s2, th in zip(rows, lats, lons, s1s_clip, s2s, ths):
        # clamp to 2×scale_ so outlier sites don't swamp the map
        w = min(2.0 * scale_ * (s1 / _s1ref), 2.0 * scale_)
        h = min(2.0 * scale_ * (s2 / _s1ref), 2.0 * scale_)
        cv = float(row.get(c_by, np.nan))
        fc = cm_obj(norm_c(cv)) if np.isfinite(cv) else (0.8, 0.8, 0.8, 0.8)
        # 3-D flagging
        beta = float(row.get("skew", row.get("beta", 0.0)))
        is_3d = (skew_threshold_ is not None and
                 mark_3d_ and
                 abs(beta) > skew_threshold_)
        lw = linewidth_ * lw_3d_factor_ if is_3d else linewidth_
        ec = edgecolor_
        e = Ellipse(
            xy=(lon, lat),
            width=w, height=h,
            angle=th,
            facecolor=fc, edgecolor=ec,
            linewidth=lw, alpha=alpha_,
            zorder=2,
        )
        ax.add_patch(e)

    # ── 10. tipper arrows ─────────────────────────────────────────────────
    if show_tipper:
        _tip_data: Dict[str, Tuple[float, float]] = {}
        for i, ed in enumerate(_iter_items(S)):
            st = _name(ed, i)
            if st not in coords:
                continue
            _T, t, fr = _get_t_block(ed)
            if t is None or fr is None:
                continue
            per_arr = 1.0 / np.where(fr == 0, np.nan, fr)
            j = int(np.nanargmin(np.abs(per_arr - period)))
            tx, ty = t[j, 0], t[j, 1]
            if tipper_component == "real":
                u, v = float(np.real(ty)), float(np.real(tx))  # E, N
            elif tipper_component == "imag":
                u, v = float(np.imag(ty)), float(np.imag(tx))
            else:
                u, v = float(np.real(ty)), float(np.real(tx))
            if tipper_convention == "parkinson":
                u, v = -u, -v
            _tip_data[st] = (u, v)

        if _tip_data:
            # normalise arrow lengths
            magnitudes = np.array([np.hypot(u, v) for u, v in _tip_data.values()], float)
            mag_max = float(np.nanpercentile(magnitudes, 95)) + 1e-12
            t_scale = (tipper_scale if tipper_scale is not None
                       else scale_ * 1.6)

            t_lons, t_lats, t_u, t_v = [], [], [], []
            for st_t, (u, v) in _tip_data.items():
                lat_t, lon_t = coords[st_t]
                t_lons.append(lon_t)
                t_lats.append(lat_t)
                t_u.append(u / mag_max * t_scale)
                t_v.append(v / mag_max * t_scale)

            ax.quiver(
                np.array(t_lons), np.array(t_lats),
                np.array(t_u), np.array(t_v),
                color=tipper_color,
                scale=1.0, scale_units="xy", angles="xy",
                width=0.003, headwidth=4, headlength=5,
                linewidth=tipper_lw,
                zorder=4,
            )

            # imaginary component if "both"
            if tipper_component == "both":
                _tip_im: Dict[str, Tuple[float, float]] = {}
                for i, ed in enumerate(_iter_items(S)):
                    st = _name(ed, i)
                    if st not in coords:
                        continue
                    _T, t, fr = _get_t_block(ed)
                    if t is None or fr is None:
                        continue
                    per_arr = 1.0 / np.where(fr == 0, np.nan, fr)
                    j = int(np.nanargmin(np.abs(per_arr - period)))
                    tx, ty = t[j, 0], t[j, 1]
                    u_i, v_i = float(np.imag(ty)), float(np.imag(tx))
                    if tipper_convention == "parkinson":
                        u_i, v_i = -u_i, -v_i
                    _tip_im[st] = (u_i, v_i)
                if _tip_im:
                    im_lons, im_lats, im_u, im_v = [], [], [], []
                    for st_t, (u, v) in _tip_im.items():
                        lat_t, lon_t = coords[st_t]
                        im_lons.append(lon_t)
                        im_lats.append(lat_t)
                        im_u.append(u / mag_max * t_scale)
                        im_v.append(v / mag_max * t_scale)
                    ax.quiver(
                        np.array(im_lons), np.array(im_lats),
                        np.array(im_u), np.array(im_v),
                        color=_mc.to_rgba(tipper_color, 0.55),
                        scale=1.0, scale_units="xy", angles="xy",
                        width=0.002, headwidth=3, headlength=4,
                        linestyle="dashed",
                        zorder=3,
                    )

    # ── 11. station markers + labels ──────────────────────────────────────
    ax.scatter(lons, lats, marker=station_marker, c=station_color,
               s=station_ms**2, zorder=5, linewidths=0)
    if station_labels:
        offset_y = scale_ * 0.8
        for st, lat, lon in zip(st_names, lats, lons):
            ax.annotate(
                st, (lon, lat),
                xytext=(0, offset_y * 111000 * 0.00001),
                textcoords="offset points",
                fontsize=label_fontsize, color=station_color,
                ha="center", va="bottom", zorder=6,
            )

    # ── 12. reference ellipse ─────────────────────────────────────────────
    if ref_ellipse_:
        pad = scale_ * 0.7
        ref_lon = float(lons.min()) - pad
        ref_lat = float(lats.min()) - pad
        ref_w = 2.0 * scale_
        ref_h = 2.0 * scale_
        ref_e = Ellipse(
            (ref_lon, ref_lat), width=ref_w, height=ref_h,
            facecolor="none", edgecolor=_es.ref_edgecolor,
            linewidth=_es.ref_lw, zorder=5,
        )
        ax.add_patch(ref_e)
        ax.text(ref_lon, ref_lat - scale_ * 1.3,
                "ref", ha="center", va="top",
                fontsize=_es.ref_fontsize, color=_es.ref_edgecolor)

    # ── 13. colorbar ──────────────────────────────────────────────────────
    sm = plt.cm.ScalarMappable(cmap=cm_obj, norm=norm_c)
    sm.set_array([])
    cb_label = colorbar_label if colorbar_label is not None else c_by
    add_colorbar(sm, ax, label=cb_label, side="right", size="3.5%", pad=0.06)

    # ── 14. axes decoration ───────────────────────────────────────────────
    from matplotlib.ticker import ScalarFormatter
    pad_xy = scale_ * 2.2
    ax.set_xlim(lons.min() - pad_xy, lons.max() + pad_xy)
    ax.set_ylim(lats.min() - pad_xy, lats.max() + pad_xy)
    # geographic aspect ratio: 1° lon = cos(lat) × 1° lat in physical distance
    lat_mid = float(np.mean(lats))
    ax.set_aspect(1.0 / max(np.cos(np.radians(lat_mid)), 1e-6))
    ax.set_xlabel("Longitude", fontsize=9)
    ax.set_ylabel("Latitude",  fontsize=9)
    ax.tick_params(labelsize=8)
    # suppress matplotlib offset notation (e.g. "+1.1912e2")
    for _axis in (ax.xaxis, ax.yaxis):
        _fmt = ScalarFormatter(useOffset=False)
        _fmt.set_scientific(False)
        _axis.set_major_formatter(_fmt)
    ax.grid(True, alpha=0.18, linewidth=0.5)
    period_txt = (f"{period*1e3:.4g} ms" if period < 1.0
                  else f"{period:.4g} s")
    ax.set_title(
        title or f"Phase Tensor Map for {period_txt}",
        fontsize=10, pad=8, fontweight="bold",
    )
    return ax



def plot_phase_tensor_summary(
    sites: Any,
    *,
    stations: Optional[List[str]] = None,
    period_range: Optional[Tuple[float, float]] = None,
    scale          = _UNSET,   # default: PYCSAMT_STYLE.pt_ellipse.scale
    c_by           = _UNSET,   # default: PYCSAMT_STYLE.pt_ellipse.c_by
    cmap           = _UNSET,   # default: auto from c_by via resolve_cmap()
    clim: Optional[Tuple[float, float]] = None,
    skew_threshold = _UNSET,   # default: PYCSAMT_STYLE.pt_ellipse.skew_threshold
    ellipt_threshold: float = 0.2,
    figsize: Tuple[float, float] = (12.0, 9.0),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> plt.Figure:
    """
    Three-panel phase-tensor summary figure.

    Panels
    ------
    **(a)** Ellipse pseudo-section — the full :func:`plot_phase_tensor_psection`
        with all cell-size normalisations.
    **(b)** Dimensionality classification grid — 1-D (white), 2-D (grey),
        3-D (red), derived from *skew_threshold* and *ellipt_threshold*.
    **(c)** Skew–ellipticity density hexbin — joint distribution over all
        ``(station, period)`` cells, with contours at 20 / 50 / 80 % of the
        peak density.

    Parameters
    ----------
    sites : any
        Passed directly to :func:`build_phase_tensor_table`.
    stations : list of str or None
        Restrict to a station subset (applied to all panels).
    period_range : (T_min, T_max) or None
        Period window in seconds (applied to all panels).
    scale : float, default ``0.85``
        Passed to panel (a) — fraction of each cell to fill.
    c_by : str, default ``"skew"``
        Colour quantity for panel (a); see :func:`plot_phase_tensor_psection`.
    cmap : str, default ``"RdBu_r"``
        Colourmap for panel (a).
    clim : (vmin, vmax) or None
        Explicit colour limits for panel (a).
    skew_threshold : float, default ``3.0``
        |β| threshold (°) for 1-D/2-D vs 3-D in panels (a) and (b).
    ellipt_threshold : float, default ``0.2``
        Ellipticity threshold for 1-D vs 2-D in panel (b).
    figsize : (width, height), default ``(12.0, 9.0)``
    recursive, on_dup, strict, verbose : see :func:`ensure_sites`.

    Returns
    -------
    fig : :class:`~matplotlib.figure.Figure`
    """
    # ── resolve visual style from PYCSAMT_STYLE.pt_ellipse ──────────────
    _es = PYCSAMT_STYLE.pt_ellipse
    if scale          is _UNSET: scale          = _es.scale
    if c_by           is _UNSET: c_by           = _es.c_by
    if cmap           is _UNSET: cmap           = _es.copy(c_by=c_by).resolve_cmap()
    if skew_threshold is _UNSET: skew_threshold = _es.skew_threshold

    df = build_phase_tensor_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    # Apply filters once; pass the already-filtered data to sub-functions
    if not df.empty and stations is not None:
        df = df[df["station"].isin(stations)]
    if not df.empty and period_range is not None:
        lo, hi = float(period_range[0]), float(period_range[1])
        df = df[(df["period"] >= lo) & (df["period"] <= hi)]

    # ── layout: 2 rows × 2 cols; row-0 spans both columns ────────────────
    import matplotlib.gridspec as gridspec
    fig = plt.figure(figsize=figsize)
    gs  = gridspec.GridSpec(
        2, 2,
        figure=fig,
        height_ratios=[1.5, 1.0],
        hspace=0.38, wspace=0.35,
    )
    ax_ell  = fig.add_subplot(gs[0, :])   # (a) full-width ellipse section
    ax_dim  = fig.add_subplot(gs[1, 0])   # (b) dimensionality grid
    ax_den  = fig.add_subplot(gs[1, 1])   # (c) skew-ellipticity density

    # ── panel (a): ellipse pseudo-section ────────────────────────────────
    if df.empty:
        ax_ell.text(0.5, 0.5, "no phase tensor data",
                    ha="center", va="center", transform=ax_ell.transAxes)
    else:
        plot_phase_tensor_psection(
            df,                     # pass pre-filtered DataFrame directly
            scale=scale,
            c_by=c_by,
            cmap=cmap,
            clim=clim,
            skew_threshold=skew_threshold,
            mark_3d=True,
            ref_ellipse=True,
            title="(a) Phase-tensor ellipse pseudo-section",
            recursive=False,
            ax=ax_ell,
        )

    # ── panel (b): dimensionality grid ────────────────────────────────────
    if df.empty:
        ax_dim.text(0.5, 0.5, "no data", ha="center", va="center",
                    transform=ax_dim.transAxes)
    else:
        _draw_dim_grid(df, ax_dim, skew_threshold=skew_threshold,
                       ellipt_threshold=ellipt_threshold)
        ax_dim.set_title("(b) Dimensionality (1-D/2-D/3-D)", fontsize=9)

    # ── panel (c): skew–ellipticity density ───────────────────────────────
    if df.empty:
        ax_den.text(0.5, 0.5, "no data", ha="center", va="center",
                    transform=ax_den.transAxes)
    else:
        _draw_skew_ellipt_density(df, ax_den, skew_threshold=skew_threshold,
                                   ellipt_threshold=ellipt_threshold)
        ax_den.set_title("(c) Skew–ellipticity distribution", fontsize=9)

    return fig


def _draw_dim_grid(
    df: pd.DataFrame,
    ax: plt.Axes,
    *,
    skew_threshold: float = 3.0,
    ellipt_threshold: float = 0.2,
) -> None:
    """Internal helper: draw the dimensionality heatmap onto *ax*."""
    df = df.copy()
    df["_logp"] = np.log10(df["period"].to_numpy(float))
    a = np.abs(df["skew"].to_numpy(float))
    e = np.abs(df["ellipt"].to_numpy(float))
    lab = np.full(len(df), 2, dtype=int)          # default: 3-D
    lab[(a <= skew_threshold) & (e <= ellipt_threshold)] = 0   # 1-D
    lab[(a <= skew_threshold) & (e  > ellipt_threshold)] = 1   # 2-D
    df["_dim"] = lab

    st_list = list(dict.fromkeys(df["station"].tolist()))
    st_map  = {s: i for i, s in enumerate(st_list)}

    piv = df.pivot_table(
        index="_logp", columns="station", values="_dim", aggfunc="median",
    )
    piv = piv.reindex(columns=st_list).sort_index()
    # piv has shape (n_logp, n_stations); do NOT transpose so that
    # imshow maps x → stations and y → log-periods.
    Z = piv.to_numpy(dtype=float)            # (n_logp, n_st)

    cmap_d = mcolors.ListedColormap(["#f7f7f7", "#969696", "#d6604d"])
    bounds = [-0.5, 0.5, 1.5, 2.5]
    norm_d = mcolors.BoundaryNorm(bounds, cmap_d.N)

    ax.imshow(Z, aspect="auto", origin="lower", interpolation="nearest",
              cmap=cmap_d, norm=norm_d)
    ax.set_xlabel("Station", fontsize=8)
    ax.set_ylabel(r"$\log_{10}(T)$ (s)", fontsize=8)
    n_st = len(st_list)
    ax.set_xticks(np.arange(n_st))
    ax.set_xticklabels(st_list, rotation=45, ha="right", fontsize=7)
    n_logp = Z.shape[0]   # shape[0] = n_logp now that .T is removed
    y_ticks = np.linspace(0, n_logp - 1, num=min(6, n_logp))
    y_vals  = np.linspace(piv.index.min(), piv.index.max(), num=min(6, n_logp))
    ax.set_yticks(y_ticks)
    ax.set_yticklabels([f"{v:.1f}" for v in y_vals], fontsize=7)

    from matplotlib.patches import Patch
    ax.legend(
        handles=[
            Patch(fc="#f7f7f7", ec="k", lw=0.5, label="1-D"),
            Patch(fc="#969696", ec="k", lw=0.5, label="2-D"),
            Patch(fc="#d6604d", ec="k", lw=0.5, label="3-D"),
        ],
        fontsize=7, framealpha=0.9, loc="lower right",
        title=f"|β|≤{skew_threshold:.0f}°,ε≤{ellipt_threshold:.2f}", title_fontsize=6,
    )


def _draw_skew_ellipt_density(
    df: pd.DataFrame,
    ax: plt.Axes,
    *,
    skew_threshold: float = 3.0,
    ellipt_threshold: float = 0.2,
    gridsize: int = 30,
) -> None:
    """Internal helper: draw |β|–ellipticity hexbin density onto *ax*."""
    a = np.abs(df["skew"].to_numpy(float))
    e = np.abs(df["ellipt"].to_numpy(float))
    fin = np.isfinite(a) & np.isfinite(e)
    if not fin.any():
        ax.text(0.5, 0.5, "no data", ha="center", va="center",
                transform=ax.transAxes)
        return

    ax.hexbin(a[fin], e[fin], gridsize=gridsize, mincnt=1,
              cmap="YlOrRd", linewidths=0.0)
    # contour at 20 / 50 / 80 % of peak density
    nx = ny = gridsize
    H, xe, ye = np.histogram2d(a[fin], e[fin], bins=[nx, ny])
    if H.max() > 0:
        Xc = 0.5 * (xe[1:] + xe[:-1])
        Yc = 0.5 * (ye[1:] + ye[:-1])
        levs = [H.max() * f for f in (0.20, 0.50, 0.80)]
        ax.contour(Xc, Yc, H.T, levels=levs, colors="k", linewidths=0.6)

    # threshold lines
    ax.axvline(skew_threshold, color="#b2182b", lw=1.0, ls="--",
               label=f"|β|={skew_threshold:.0f}°")
    ax.axhline(ellipt_threshold, color="#2166ac", lw=1.0, ls="--",
               label=f"ε={ellipt_threshold:.2f}")
    ax.set_xlabel("|β| (°)", fontsize=8)
    ax.set_ylabel("Ellipticity", fontsize=8)
    ax.legend(fontsize=7, framealpha=0.9)


def phase_tensor_legend(
    *,
    size: float = 1.0,
    figsize: Tuple[float, float] = (2.5, 2.5),
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    e = Ellipse((0.0, 0.0), width=size, height=size,
                angle=0.0, fill=False, lw=1.0)
    ax.add_patch(e)
    ax.plot([0, 0], [0, size * 0.6], "-", lw=1.0)
    ax.set_xlim(-size, size)
    ax.set_ylim(-size, size)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    return ax



# -------------------- small utilities (local) --------------------------- #

def _rho_det_from_z(z: np.ndarray, fr: np.ndarray) -> np.ndarray:
    zx = z[:, 0, 1]
    zy = z[:, 1, 0]
    rx = 0.2 * (np.abs(zx) ** 2) / (fr + 1e-24)
    ry = 0.2 * (np.abs(zy) ** 2) / (fr + 1e-24)
    rdet = np.sqrt(rx * ry)
    return rdet


def _phi_det_from_z(z: np.ndarray) -> np.ndarray:
    d = z[:, 0, 0] * z[:, 1, 1] - z[:, 0, 1] * z[:, 1, 0]
    return np.degrees(np.angle(d))


def _attach_det_metrics(
    df: pd.DataFrame,
    S,
) -> pd.DataFrame:
    if df.empty:
        return df
    df = df.copy()
    df["rho_det"] = np.nan
    df["phi_det"] = np.nan
    for i, ed in enumerate(_iter_items(S)):
        st = _name(ed, i)
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            continue
        rho = _rho_det_from_z(z, fr)
        phi = _phi_det_from_z(z)
        per = 1.0 / fr
        m = df["station"] == st
        if not m.any():
            continue
        p = df.loc[m, "period"].to_numpy()
        idx = np.searchsorted(per, p)
        idx = np.clip(idx, 0, len(per) - 1)
        df.loc[m, "rho_det"] = rho[idx]
        df.loc[m, "phi_det"] = phi[idx]
    return df


def _color_from_hsv(h: np.ndarray, s: np.ndarray,
                    v: np.ndarray) -> np.ndarray:
    hsv = np.stack([h, s, v], axis=-1)
    rgb = mcolors.hsv_to_rgb(hsv)
    return rgb

# ---------------- 2) dimensionality grid (1D/2D/3D) --------------------- #
def plot_dimensionality_grid(
    sites: Any,
    *,
    skew_th: float = 3.0,
    ellipt_th: float = 0.2,
    figsize: Tuple[float, float] = (8.5, 4.0),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    df = build_phase_tensor_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if df.empty:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no phase tensor",
                ha="center", va="center")
        return ax
    df = df.copy()
    df["logp"] = np.log10(df["period"].to_numpy())
    lab = np.zeros(len(df), dtype=int)
    a = np.abs(df["skew"].to_numpy())
    e = np.abs(df["ellipt"].to_numpy())
    lab[:] = 2
    lab[(a <= skew_th) & (e <= ellipt_th)] = 0
    lab[(a <= skew_th) & (e > ellipt_th)] = 1
    df["dim"] = lab
    piv = df.pivot_table(
        index="logp",
        columns="station",
        values="dim",
        aggfunc="median",
    )
    piv = piv.sort_index()
    Z = piv.to_numpy(dtype=float)        # (n_logp, n_st) — no transpose
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(
        Z,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
    )
    ax.set_xlabel("Station")
    ax.set_ylabel("LogPeriod (s)")
    ax.set_xticks(np.arange(len(piv.columns)))
    ax.set_xticklabels(piv.columns, rotation=90)
    yt = np.linspace(0, Z.shape[0] - 1, num=min(8, Z.shape[0]))  # shape[0] = n_logp
    yv = np.linspace(piv.index.min(), piv.index.max(),
                     num=min(8, len(piv.index)))
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{v:.2g}" for v in yv])
    cb = plt.colorbar(im, ax=ax)
    cb.set_label("dim (0=1D,1=2D,2=3D)")
    return ax


# ---------------- 3) theta stability stripe (HSV image) ----------------- #

def plot_theta_stability_stripe(
    sites: Any,
    *,
    win: int = 5,
    figsize: Tuple[float, float] = (9.0, 4.0),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    df = build_phase_tensor_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if df.empty:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no phase tensor",
                ha="center", va="center")
        return ax
    df = df.copy()
    df["logp"] = np.log10(df["period"].to_numpy())
    sts = list(df["station"].unique())
    X = []
    H = []
    for st in sts:
        s = df[df["station"] == st].sort_values("logp")
        th = s["theta"].to_numpy()
        lp = s["logp"].to_numpy()
        # hue from theta in [0,180) → [0,1)
        h = ((th % 180.0) / 180.0)
        # variance in sliding window
        k = max(3, int(win))
        if len(h) < k:
            v = np.full_like(h, np.nan)
        else:
            v = np.convolve(
                (h - np.nanmean(h)) ** 2,
                np.ones(k) / k,
                mode="same",
            )
        # saturation from 1 - norm(var)
        v0 = np.nanpercentile(v[np.isfinite(v)], 5) if np.isfinite(
            v
        ).any() else 0.0
        v1 = np.nanpercentile(v[np.isfinite(v)], 95) if np.isfinite(
            v
        ).any() else 1.0
        s_sat = 1.0 - np.clip((v - v0) / (v1 - v0 + 1e-12), 0, 1)
        H.append(np.vstack([h, s_sat, np.ones_like(h)]))
        X.append(lp)
    # align to common y-grid for imshow
    yall = np.unique(np.concatenate(X))
    Himg = np.zeros((len(yall), len(sts), 3))
    for j, (lp, hs) in enumerate(zip(X, H)):
        i = np.searchsorted(yall, lp)
        i = np.clip(i, 0, len(yall) - 1)
        h, s, v = hs
        rgb = _color_from_hsv(h, s, v)
        for r, k in enumerate(i):
            Himg[k, j, :] = rgb[r]
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    ax.imshow(
        Himg,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
    )
    ax.set_xlabel("Station")
    ax.set_ylabel("LogPeriod (s)")
    ax.set_xticks(np.arange(len(sts)))
    ax.set_xticklabels(sts, rotation=90)
    yt = np.linspace(0, len(yall) - 1, num=min(8, len(yall)))
    yv = np.linspace(yall.min(), yall.max(),
                     num=min(8, len(yall)))
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{v:.2g}" for v in yv])
    return ax


# ---------------- 4) skew–ellipticity density (hexbin + contours) ------- #

def plot_skew_ellipt_density(
    sites: Any,
    *,
    band: Optional[Tuple[float, float]] = None,
    gridsize: int = 40,
    figsize: Tuple[float, float] = (6.5, 5.5),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    df = build_phase_tensor_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if df.empty:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no phase tensor",
                ha="center", va="center")
        return ax
    a = np.abs(df["beta"].to_numpy())
    e = np.abs(df["ellipt"].to_numpy())
    if band is not None:
        lo, hi = band
        m = (df["period"] >= lo) & (df["period"] <= hi)
        a = a[m]
        e = e[m]
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    ax.hexbin(
        a, e, gridsize=gridsize, mincnt=1, linewidths=0.0
    )
    ax.set_xlabel("|beta| (deg)")
    ax.set_ylabel("|ellipticity|")
    # contours from 2D hist
    nx = ny = gridsize
    H, xedges, yedges = np.histogram2d(a, e, bins=[nx, ny])
    Xc = 0.5 * (xedges[1:] + xedges[:-1])
    Yc = 0.5 * (yedges[1:] + yedges[:-1])
    lev = np.linspace(H.max() * 0.2, H.max() * 0.9, 4)
    ax.contour(Xc, Yc, H.T, levels=lev, colors="k", linewidths=0.6)
    return ax


# ---------------- 5) theta rose per band (small multiples) -------------- #
def plot_theta_rose_grid(
    sites: Any,
    *,
    n_bands: int = 6,
    figsize: Tuple[float, float] = (13.0, 3.8),
    bins: int = 24,
    style: "str | RoseStyle | None" = "pycsamt",
    panel_title_fontsize: float = 7.5,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> plt.Figure:
    """Phase-tensor θ rose grid — one pycsamt-styled rose per frequency decade.

    Draws *n_bands* polar rose diagrams side by side, each covering one
    equal-log-width period band.  Each rose applies axial symmetry (0–180°
    mirrored to 180–360°) and is rendered using the active
    :class:`~pycsamt.api._rose_style.RoseStyle`.

    Parameters
    ----------
    sites : any
        Input accepted by :func:`~pycsamt.emtools._core.ensure_sites`.
    n_bands : int, default ``6``
        Number of equal-log-width period bands.
    figsize : (width, height), default ``(13.0, 3.8)``
        Figure size in inches.  Constrained layout is used internally so
        an external ``fig.suptitle`` fits without producing blank space.
    bins : int, default ``24``
        Number of bins over 0–180° (mirrored to 360°).
    style : str, :class:`~pycsamt.api._rose_style.RoseStyle`, or None
        Rose visual style.  Strings resolved via
        :func:`~pycsamt.api._rose_style.resolve_rose_style`.
    panel_title_fontsize : float, default ``7.5``
        Font size for the period-band label above each panel.
    recursive, on_dup, strict, verbose
        Passed to :func:`~pycsamt.emtools._core.ensure_sites`.

    Returns
    -------
    :class:`matplotlib.figure.Figure`
    """
    # ── resolve rose style ────────────────────────────────────────────────
    rs = resolve_rose_style(style)

    # ── axial mean helper ─────────────────────────────────────────────────
    def _axial_mean(a_deg: np.ndarray) -> float:
        rad = np.radians(2.0 * a_deg)
        mu = np.degrees(np.arctan2(np.nanmean(np.sin(rad)),
                                   np.nanmean(np.cos(rad)))) / 2.0
        return float(mu % 180.0)

    df = build_phase_tensor_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if df.empty:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, polar=True)
        ax.text(0.5, 0.5, "no phase tensor", ha="center", va="center")
        return fig

    p = df["period"].to_numpy()
    lo_g = float(np.nanmin(p[p > 0])) if np.any(p > 0) else 1e-6
    hi_g = float(np.nanmax(p))
    edges = np.logspace(np.log10(lo_g), np.log10(hi_g), n_bands + 1)

    bins_ = int(max(12, bins))
    edges_deg = np.linspace(0.0, 180.0, bins_ + 1)
    dw = np.radians(180.0 / bins_)
    ang_mid = np.radians(0.5 * (edges_deg[1:] + edges_deg[:-1]))

    # ── constrained_layout avoids blank space when suptitle is added ─────
    fig = plt.figure(figsize=figsize, constrained_layout=True)

    for i in range(n_bands):
        ax = fig.add_subplot(1, n_bands, i + 1, polar=True)
        m = (p >= edges[i]) & (p < edges[i + 1])
        th_deg = df.loc[m, "theta"].to_numpy(float) % 180.0
        n_obs = int(np.sum(m))

        # histogram — 0–180°, mirrored for axial symmetry
        hist, _ = np.histogram(th_deg, bins=edges_deg)
        rmax = float(hist.max()) if hist.max() > 0 else 1.0

        # ── polar setup ───────────────────────────────────────────────────
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_yticklabels([])
        ax.set_yticks([])
        ax.yaxis.grid(False)
        ax.xaxis.grid(False)
        ax.set_frame_on(False)

        # ── bars with axial symmetry (0–180° + mirror 180–360°) ──────────
        if rs.bar_style == "gradient":
            cm_ = plt.get_cmap(rs.cmap)
            for a, h in zip(ang_mid, hist):
                col = cm_(h / rmax) if rmax > 0 else cm_(0.5)
                for a_plot in (a, a + np.pi):
                    ax.bar(a_plot, h, width=dw,
                           color=col,
                           edgecolor=rs.bar_edgecolor,
                           linewidth=rs.bar_edgelw,
                           alpha=rs.bar_alpha)
        else:  # solid
            for a, h in zip(ang_mid, hist):
                for a_plot in (a, a + np.pi):
                    ax.bar(a_plot, h, width=dw,
                           color=rs.bar_color,
                           edgecolor=rs.bar_edgecolor,
                           linewidth=rs.bar_edgelw,
                           alpha=rs.bar_alpha)

        # ── concentric rings ──────────────────────────────────────────────
        step = rmax / rs.n_rings
        r_levels = [step * k for k in range(1, rs.n_rings + 1)]
        theta_full = np.linspace(0, 2 * np.pi, 360)
        for rv in r_levels:
            ax.plot(theta_full, np.full_like(theta_full, rv),
                    color=rs.ring_color, ls=rs.ring_ls,
                    lw=rs.ring_lw, zorder=0)

        # ring count labels (outermost ring only to avoid clutter)
        lbl_theta = np.radians(rs.ring_label_angle)
        ax.text(lbl_theta, r_levels[-1],
                rs.ring_label_fmt.format(r_levels[-1]),
                ha="center", va="center",
                fontsize=rs.ring_label_fontsize,
                color=rs.ring_label_color, zorder=5)

        # ── radial spokes ─────────────────────────────────────────────────
        n_spokes = int(round(360.0 / rs.spoke_every))
        spoke_angles = np.radians(np.arange(n_spokes) * rs.spoke_every)
        for sa in spoke_angles:
            ax.plot([sa, sa], [0, rmax],
                    color=rs.spoke_color, ls=rs.spoke_ls,
                    lw=rs.spoke_lw, zorder=0)

        # ── outer bold ring ────────────────────────────────────────────────
        ax.spines["polar"].set_linewidth(rs.outer_ring_lw)
        ax.spines["polar"].set_edgecolor(rs.outer_ring_color)
        ax.set_ylim(0, rmax * 1.08)

        # ── compass / degree labels ───────────────────────────────────────
        spoke_degs = np.arange(n_spokes) * rs.spoke_every
        if rs.compass_labels == "NESW":
            _card = {0: "N", 90: "E", 180: "S", 270: "W"}
            lbl_list = [_card.get(int(d) % 360, "") for d in spoke_degs]
        elif rs.compass_labels == "degrees":
            lbl_list = [f"{int(d)}°" for d in spoke_degs]
        else:
            lbl_list = [""] * n_spokes
        ax.set_thetagrids(spoke_degs, labels=lbl_list,
                          fontsize=rs.compass_fontsize)
        for lbl in ax.get_xticklabels():
            lbl.set_color(rs.compass_color)
            lbl.set_fontweight(rs.compass_fontweight)

        # ── mean direction line ───────────────────────────────────────────
        if rs.show_mean and n_obs > 0:
            mu = _axial_mean(th_deg)
            mu_r = np.radians(mu)
            sec_lw = rs.mean_lw if rs.secondary_lw is None else rs.secondary_lw
            sec_col = rs.mean_color if rs.secondary_color is None else rs.secondary_color
            for ang_pt in (mu_r, mu_r + np.pi):
                ax.plot([ang_pt, ang_pt], [0, rmax],
                        color=rs.mean_color, lw=rs.mean_lw,
                        ls=rs.mean_ls, zorder=6)
            if rs.show_secondary:
                for ang_pt in (mu_r + np.pi / 2, mu_r - np.pi / 2):
                    ax.plot([ang_pt, ang_pt], [0, rmax],
                            color=sec_col, lw=sec_lw,
                            ls=rs.secondary_ls, zorder=6)

        # ── panel annotation: period range + n ────────────────────────────
        if rs.show_annotation and n_obs > 0:
            mu_str = ""
            if rs.show_mean and n_obs > 0:
                mu = _axial_mean(th_deg)
                mu_str = f"θ̄={mu:.0f}°"
            n_str = f"  n={n_obs}" if rs.show_n else ""
            ax.text(
                *rs.annotation_pos,
                mu_str + n_str,
                transform=ax.transAxes,
                fontsize=rs.annotation_fontsize,
                ha="left", va="top",
                bbox=dict(fc=rs.annotation_bg,
                          ec=rs.annotation_ec,
                          alpha=0.80, pad=1.5,
                          boxstyle="round,pad=0.3"),
                zorder=10,
            )

        # ── panel title: period band ───────────────────────────────────────
        ax.set_title(
            f"[{edges[i]:.2g}, {edges[i+1]:.2g}]s",
            fontsize=panel_title_fontsize,
            pad=8,
        )

    return fig

