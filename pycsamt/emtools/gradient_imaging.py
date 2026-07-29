"""
Gradient-based apparent resistivity pseudo-sections for CSAMT.

Implements the frequency-gradient, spatial-gradient, and joint
frequency-spatial gradient apparent resistivity formulations proposed by:

  zhang2021 : Zhang, Farquharson & Liu (2021), "Improved CSAMT apparent
              resistivity pseudo sections based on the frequency and
              frequency-spatial gradients of electromagnetic fields",
              *Geophysical Prospecting*, doi:10.1111/1365-2478.13059.

The three key gradient quantities are:

Transverse (along-line) gradient :math:`\\Delta\\rho_a^x`
    Finite difference of :math:`\\rho_a` between adjacent stations at each
    frequency.

Vertical (frequency) gradient :math:`\\Delta\\rho_a^z`
    Log-frequency finite difference of :math:`\\rho_a` at each station.

Joint vertical-transverse gradient :math:`\\Delta\\rho_a^{zx}`
    Frequency difference of the spatial gradient.  This is the principal
    product of zhang2021: it suppresses background interference caused by
    the non-zero spatial gradient over homogeneous half-spaces and
    improves boundary delineation in pseudo-section images.

All three quantities are derived from the standard impedance-based
:math:`\\rho_a`, so no source moment or source–receiver geometry is
required.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from ._core import (
    _get_z_block,
    _iter_items,
    _name,
    _station_positions,
    ensure_sites,
)

__all__ = [
    "rho_spatial_gradient",
    "rho_frequency_gradient",
    "rho_joint_gradient",
    "plot_gradient_section",
]

# ─────────────────────────────────────────────────────────────────────────────
# Column name constants
# ─────────────────────────────────────────────────────────────────────────────

_SPATIAL_COLS = [
    "station_a",
    "station_b",
    "x_m",
    "dx_m",
    "freq_hz",
    "period_s",
    "depth_m",
    "rho_a_ohmm",
    "delta_rho_x",
]

_FREQ_COLS = [
    "station",
    "x_m",
    "freq_hz",
    "period_s",
    "depth_m",
    "rho_a_ohmm",
    "delta_rho_z",
]

_JOINT_COLS = [
    "station_a",
    "station_b",
    "x_m",
    "dx_m",
    "freq_hz",
    "period_s",
    "depth_m",
    "delta_rho_zx",
]

# ─────────────────────────────────────────────────────────────────────────────
# Private helpers
# ─────────────────────────────────────────────────────────────────────────────


def _unwrap(ed: Any) -> Any:
    """Unwrap a Sites-level Site to its underlying EDI-like object."""
    edi = getattr(ed, "edi", None)
    if edi is not None and hasattr(edi, "Z"):
        return edi
    return ed


def _rho_a_comp(z: np.ndarray, fr: np.ndarray, comp: str) -> np.ndarray:
    """Apparent resistivity [Ω·m] for the requested impedance component."""
    if comp == "xy":
        return 0.2 * np.abs(z[:, 0, 1]) ** 2 / np.maximum(fr, 1e-24)
    if comp == "yx":
        return 0.2 * np.abs(z[:, 1, 0]) ** 2 / np.maximum(fr, 1e-24)
    rxy = 0.2 * np.abs(z[:, 0, 1]) ** 2 / np.maximum(fr, 1e-24)
    ryx = 0.2 * np.abs(z[:, 1, 0]) ** 2 / np.maximum(fr, 1e-24)
    return np.sqrt(np.maximum(rxy * ryx, 1e-12))


def _skin_depth_m(rho: float, freq: float) -> float:
    """Apparent depth δ = 503 √(ρ/f)  [m]."""
    return 503.0 * float(np.sqrt(max(rho, 1e-6) / max(freq, 1e-12)))


def _build_rho_grid(
    sites: Any,
    comp: str,
    spacing_m: float,
    recursive: bool,
    on_dup: str,
    strict: bool,
    verbose: int,
) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray]:
    """
    Build a 2-D apparent resistivity grid from a ``Sites`` object.

    Parameters
    ----------
    sites : Sites | list
    comp : {"det", "xy", "yx"}
    spacing_m : float
        Fall-back station spacing when no coordinate metadata is found.

    Returns
    -------
    stations : list of str
        Station names sorted by position along the survey line.
    x_pos : ndarray, shape (N,)
        Station positions [m].
    freqs : ndarray, shape (F,)
        Sorted unique frequencies [Hz].
    rho_grid : ndarray, shape (N, F)
        Apparent resistivity [Ω·m].  ``NaN`` where data are missing.
    """
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    items = list(_iter_items(S))
    if not items:
        return [], np.array([]), np.array([]), np.zeros((0, 0))

    raw_list: list[Any] = []
    names: list[str] = []
    freq_rho: dict[str, tuple[np.ndarray, np.ndarray]] = {}

    for i, ed in enumerate(items):
        raw = _unwrap(ed)
        nm = _name(raw, i)
        _, z, fr = _get_z_block(raw)
        if z is None or fr is None or fr.size == 0:
            continue
        raw_list.append(raw)
        names.append(nm)
        freq_rho[nm] = (fr, _rho_a_comp(z, fr, comp))

    if not names:
        return [], np.array([]), np.array([]), np.zeros((0, 0))

    x_pos = _station_positions(raw_list, spacing_m=spacing_m)

    # Sort stations by position along the line
    order = np.argsort(x_pos)
    names = [names[k] for k in order]
    x_pos = x_pos[order]

    # Common frequency axis (union of all frequencies)
    all_f: set = set()
    for fr, _ in freq_rho.values():
        all_f.update(fr.tolist())
    freqs = np.array(sorted(all_f), dtype=float)
    if freqs.size == 0:
        return names, x_pos, freqs, np.full((len(names), 0), np.nan)

    f_idx = {float(f): k for k, f in enumerate(freqs)}
    rho_grid = np.full((len(names), freqs.size), np.nan)

    for i, nm in enumerate(names):
        fr, rho = freq_rho[nm]
        for j in range(fr.size):
            k = f_idx.get(float(fr[j]))
            if k is not None:
                rho_grid[i, k] = rho[j]

    return names, x_pos, freqs, rho_grid


# ─────────────────────────────────────────────────────────────────────────────
# Public: rho_spatial_gradient
# ─────────────────────────────────────────────────────────────────────────────


def rho_spatial_gradient(
    sites: Any,
    spacing_m: float = 200.0,
    *,
    comp: str = "det",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    r"""
    Transverse (along-line) apparent resistivity gradient.

    Computes the first-order finite difference of :math:`\rho_a` between
    adjacent stations at each frequency (eq. 11 of zhang2021):

    .. math::

        \Delta\rho_a^x(j,\,f)
            \approx \rho_a(j,\,f) - \rho_a(j-1,\,f)

    where station indices are ordered by their position along the survey
    line.  The result is assigned to the spatial midpoint between the two
    stations.

    Parameters
    ----------
    sites : Sites | list
        EDI-like objects or a ``Sites`` container.
    spacing_m : float, default 200
        Fall-back inter-station spacing [m] used when no coordinate
        metadata is available.
    comp : {"det", "xy", "yx"}, default "det"
        Impedance component for :math:`\rho_a`.
        ``"det"`` uses the geometric-mean determinant.
    recursive, on_dup, strict, verbose
        Forwarded to :func:`~pycsamt.emtools._core.ensure_sites`.

    Returns
    -------
    pandas.DataFrame
        One row per (station-pair, frequency).  Columns:

        ``station_a``, ``station_b``
            Names of the left and right stations of each pair.
        ``x_m``
            Midpoint position along the survey line [m].
        ``dx_m``
            Spacing between the two stations [m].
        ``freq_hz``, ``period_s``
            Frequency [Hz] and period [s].
        ``depth_m``
            Skin depth :math:`\delta = 503\,\sqrt{\rho_a/f}` [m] at the
            midpoint :math:`\rho_a` and the given frequency.
        ``rho_a_ohmm``
            Mean :math:`\rho_a` of the pair at that frequency [Ω·m].
        ``delta_rho_x``
            :math:`\Delta\rho_a^x` [Ω·m].

    References
    ----------
    Zhang et al. (2021), eq. (11).
    """
    if comp not in ("det", "xy", "yx"):
        raise ValueError(f"comp must be 'det', 'xy', or 'yx'; got {comp!r}")

    names, x_pos, freqs, rho_grid = _build_rho_grid(
        sites, comp, spacing_m, recursive, on_dup, strict, verbose
    )

    N, F = rho_grid.shape
    if N < 2 or F == 0:
        return pd.DataFrame(columns=_SPATIAL_COLS)

    # Δρ_a^x[pair_idx, freq_idx] = rho_grid[j, f] - rho_grid[j-1, f]
    delta_x = rho_grid[1:, :] - rho_grid[:-1, :]  # (N-1, F)
    x_mid = 0.5 * (x_pos[1:] + x_pos[:-1])  # (N-1,)
    dx = x_pos[1:] - x_pos[:-1]  # (N-1,)
    rho_mid = 0.5 * (rho_grid[1:, :] + rho_grid[:-1, :])  # (N-1, F)

    rows: list[dict] = []
    for j in range(N - 1):
        for k in range(F):
            rho_m = float(rho_mid[j, k])
            rows.append(
                {
                    "station_a": names[j],
                    "station_b": names[j + 1],
                    "x_m": float(x_mid[j]),
                    "dx_m": float(dx[j]),
                    "freq_hz": float(freqs[k]),
                    "period_s": 1.0 / max(float(freqs[k]), 1e-12),
                    "depth_m": _skin_depth_m(rho_m, float(freqs[k])),
                    "rho_a_ohmm": rho_m,
                    "delta_rho_x": float(delta_x[j, k]),
                }
            )

    if not rows:
        return pd.DataFrame(columns=_SPATIAL_COLS)
    return pd.DataFrame(rows, columns=_SPATIAL_COLS)


# ─────────────────────────────────────────────────────────────────────────────
# Public: rho_frequency_gradient
# ─────────────────────────────────────────────────────────────────────────────


def rho_frequency_gradient(
    sites: Any,
    *,
    comp: str = "det",
    spacing_m: float = 200.0,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    r"""
    Vertical (log-frequency) apparent resistivity gradient.

    Computes the first-order finite difference of :math:`\rho_a` between
    adjacent frequencies at each station (eq. 12 of zhang2021):

    .. math::

        \Delta\rho_a^z(j,\,f_k)
            \approx \rho_a(j,\,f_k) - \rho_a(j,\,f_{k-1})

    where :math:`f_k > f_{k-1}` (ascending frequency, ascending
    skin-depth index).  Because different frequencies probe different
    depths, :math:`\Delta\rho_a^z` is associated with vertical changes in
    the subsurface.

    Parameters
    ----------
    sites : Sites | list
        EDI-like objects or a ``Sites`` container.
    comp : {"det", "xy", "yx"}, default "det"
        Impedance component for :math:`\rho_a`.
    spacing_m : float, default 200
        Fall-back inter-station spacing [m].
    recursive, on_dup, strict, verbose
        Forwarded to :func:`~pycsamt.emtools._core.ensure_sites`.

    Returns
    -------
    pandas.DataFrame
        One row per (station, adjacent-frequency pair).  Columns:

        ``station``
            Station name.
        ``x_m``
            Station position along the survey line [m].
        ``freq_hz``, ``period_s``
            Higher frequency :math:`f_k` of the pair [Hz] and its
            corresponding period [s].
        ``depth_m``
            Skin depth at the mean :math:`\rho_a` of the pair [m].
        ``rho_a_ohmm``
            Mean :math:`\rho_a` of the two adjacent frequencies [Ω·m].
        ``delta_rho_z``
            :math:`\Delta\rho_a^z` [Ω·m].

    References
    ----------
    Zhang et al. (2021), eq. (12).
    """
    if comp not in ("det", "xy", "yx"):
        raise ValueError(f"comp must be 'det', 'xy', or 'yx'; got {comp!r}")

    names, x_pos, freqs, rho_grid = _build_rho_grid(
        sites, comp, spacing_m, recursive, on_dup, strict, verbose
    )

    N, F = rho_grid.shape
    if N == 0 or F < 2:
        return pd.DataFrame(columns=_FREQ_COLS)

    # Δρ_a^z[station_idx, freq_pair_idx] = rho_grid[i, k] - rho_grid[i, k-1]
    # freq_pair_idx k corresponds to (freqs[k-1], freqs[k]) pair → assign to freqs[k]
    delta_z = rho_grid[:, 1:] - rho_grid[:, :-1]  # (N, F-1)
    rho_mean = 0.5 * (rho_grid[:, 1:] + rho_grid[:, :-1])  # (N, F-1)
    freq_k = freqs[1:]  # (F-1,) upper frequencies

    rows: list[dict] = []
    for i in range(N):
        for k in range(F - 1):
            rho_m = float(rho_mean[i, k])
            rows.append(
                {
                    "station": names[i],
                    "x_m": float(x_pos[i]),
                    "freq_hz": float(freq_k[k]),
                    "period_s": 1.0 / max(float(freq_k[k]), 1e-12),
                    "depth_m": _skin_depth_m(rho_m, float(freq_k[k])),
                    "rho_a_ohmm": rho_m,
                    "delta_rho_z": float(delta_z[i, k]),
                }
            )

    if not rows:
        return pd.DataFrame(columns=_FREQ_COLS)
    return pd.DataFrame(rows, columns=_FREQ_COLS)


# ─────────────────────────────────────────────────────────────────────────────
# Public: rho_joint_gradient
# ─────────────────────────────────────────────────────────────────────────────


def rho_joint_gradient(
    sites: Any,
    spacing_m: float = 200.0,
    *,
    comp: str = "det",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    r"""
    Joint vertical-transverse apparent resistivity gradient.

    Computes the frequency difference of the spatial gradient
    (eq. 13 of zhang2021), which simultaneously resolves lateral and
    vertical boundaries while suppressing the spurious background
    interference in the spatial gradient over homogeneous regions:

    .. math::

        \Delta\rho_a^{zx}(j,\,f_k)
            = \Delta\rho_a^x(j,\,f_k) - \Delta\rho_a^x(j,\,f_{k-1})

    Expanding in terms of :math:`\rho_a`:

    .. math::

        \Delta\rho_a^{zx}(j,\,f_k)
            = \bigl[\rho_a(j,\,f_k) - \rho_a(j-1,\,f_k)\bigr]
              - \bigl[\rho_a(j,\,f_{k-1}) - \rho_a(j-1,\,f_{k-1})\bigr]

    where station *j* is to the right of station *j-1* along the survey
    line and :math:`f_k > f_{k-1}`.

    The result is non-zero only where :math:`\rho_a` changes in **both**
    the lateral and vertical directions simultaneously, making it a
    sensitive indicator of target boundaries.

    Parameters
    ----------
    sites : Sites | list
        EDI-like objects or a ``Sites`` container.
    spacing_m : float, default 200
        Fall-back inter-station spacing [m].
    comp : {"det", "xy", "yx"}, default "det"
        Impedance component for :math:`\rho_a`.
    recursive, on_dup, strict, verbose
        Forwarded to :func:`~pycsamt.emtools._core.ensure_sites`.

    Returns
    -------
    pandas.DataFrame
        One row per (station-pair, adjacent-frequency pair).  Columns:

        ``station_a``, ``station_b``
            Names of the left and right stations of each pair.
        ``x_m``
            Midpoint position along the survey line [m].
        ``dx_m``
            Spacing between the two stations [m].
        ``freq_hz``, ``period_s``
            Upper frequency :math:`f_k` of the pair [Hz] and period [s].
        ``depth_m``
            Skin depth [m] estimated from the median :math:`\rho_a` of
            the four surrounding cells and the frequency :math:`f_k`.
        ``delta_rho_zx``
            :math:`\Delta\rho_a^{zx}` [Ω·m].

    References
    ----------
    Zhang et al. (2021), eqs. (1), (11), (13).
    """
    if comp not in ("det", "xy", "yx"):
        raise ValueError(f"comp must be 'det', 'xy', or 'yx'; got {comp!r}")

    names, x_pos, freqs, rho_grid = _build_rho_grid(
        sites, comp, spacing_m, recursive, on_dup, strict, verbose
    )

    N, F = rho_grid.shape
    if N < 2 or F < 2:
        return pd.DataFrame(columns=_JOINT_COLS)

    # spatial gradient: delta_x[j, k] = rho_grid[j+1,k] - rho_grid[j,k]  (N-1, F)
    delta_x = rho_grid[1:, :] - rho_grid[:-1, :]

    # joint gradient: delta_zx[j, k] = delta_x[j, k] - delta_x[j, k-1]  (N-1, F-1)
    delta_zx = delta_x[:, 1:] - delta_x[:, :-1]

    x_mid = 0.5 * (x_pos[1:] + x_pos[:-1])  # (N-1,)
    dx = x_pos[1:] - x_pos[:-1]  # (N-1,)
    freq_k = freqs[1:]  # (F-1,) upper frequencies

    # For depth: use median rho_a of the 4 surrounding corners
    # corners: rho_grid[j, k-1], rho_grid[j, k], rho_grid[j+1, k-1], rho_grid[j+1, k]
    rho_corners = np.nanmedian(
        np.stack(
            [
                rho_grid[:-1, :-1],  # (j,   k-1)
                rho_grid[:-1, 1:],  # (j,   k  )
                rho_grid[1:, :-1],  # (j+1, k-1)
                rho_grid[1:, 1:],  # (j+1, k  )
            ],
            axis=0,
        ),
        axis=0,
    )  # (N-1, F-1)

    rows: list[dict] = []
    for j in range(N - 1):
        for k in range(F - 1):
            rho_m = float(rho_corners[j, k])
            rows.append(
                {
                    "station_a": names[j],
                    "station_b": names[j + 1],
                    "x_m": float(x_mid[j]),
                    "dx_m": float(dx[j]),
                    "freq_hz": float(freq_k[k]),
                    "period_s": 1.0 / max(float(freq_k[k]), 1e-12),
                    "depth_m": _skin_depth_m(rho_m, float(freq_k[k])),
                    "delta_rho_zx": float(delta_zx[j, k]),
                }
            )

    if not rows:
        return pd.DataFrame(columns=_JOINT_COLS)
    return pd.DataFrame(rows, columns=_JOINT_COLS)


# ─────────────────────────────────────────────────────────────────────────────
# Public: plot_gradient_section
# ─────────────────────────────────────────────────────────────────────────────


def plot_gradient_section(
    sites: Any,
    quantity: str = "joint",
    spacing_m: float = 200.0,
    *,
    comp: str = "det",
    period_axis: bool = True,
    log_y: bool = True,
    figsize: tuple[float, float] = (10.0, 5.0),
    cmap: str = "RdBu_r",
    vlim: tuple[float, float] | None = None,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Any | None = None,
) -> Any:
    r"""
    Pseudo-section of a gradient apparent resistivity quantity.

    Produces a colour-coded pseudo-section (position × period or
    frequency) for one of three gradient quantities from zhang2021:

    * ``"spatial"`` (:math:`\Delta\rho_a^x`) — lateral boundaries.
    * ``"frequency"`` (:math:`\Delta\rho_a^z`) — vertical boundaries.
    * ``"joint"`` (:math:`\Delta\rho_a^{zx}`) — combined; best boundary
      delineation and suppressed background interference (default).

    A diverging colour-map centred at zero is used so that positive and
    negative gradient values (resistivity increase vs. decrease) can be
    distinguished at a glance.

    Parameters
    ----------
    sites : Sites | list
        EDI-like objects or a ``Sites`` container.
    quantity : {"joint", "spatial", "frequency"}, default "joint"
        Which gradient pseudo-section to plot.
    spacing_m : float, default 200
        Fall-back inter-station spacing [m].
    comp : {"det", "xy", "yx"}, default "det"
        Impedance component for :math:`\rho_a`.
    period_axis : bool, default True
        Show period [s] on the y-axis when *True*; frequency [Hz] when
        *False*.
    log_y : bool, default True
        Use a logarithmic y-axis.
    figsize : (float, float), default (10, 5)
    cmap : str, default ``"RdBu_r"``
        Diverging Matplotlib colour-map.
    vlim : (vmin, vmax) or None
        Colour-scale limits [Ω·m].  If *None*, a symmetric range centred
        at zero is derived from the data.
    ax : matplotlib.axes.Axes or None
        Draw on existing axes; a new figure is created if *None*.
    recursive, on_dup, strict, verbose
        Forwarded to :func:`~pycsamt.emtools._core.ensure_sites`.

    Returns
    -------
    matplotlib.axes.Axes

    References
    ----------
    Zhang et al. (2021), *Geophysical Prospecting*,
    doi:10.1111/1365-2478.13059.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize, TwoSlopeNorm

    _Q = {
        "joint": "joint",
        "zx": "joint",
        "spatial": "spatial",
        "x": "spatial",
        "frequency": "frequency",
        "z": "frequency",
    }
    q = _Q.get(quantity.lower())
    if q is None:
        raise ValueError(
            f"quantity must be 'joint', 'spatial', or 'frequency'; " f"got {quantity!r}"
        )

    if q == "joint":
        df = rho_joint_gradient(
            sites,
            spacing_m,
            comp=comp,
            recursive=recursive,
            on_dup=on_dup,
            strict=strict,
            verbose=verbose,
        )
        val_col = "delta_rho_zx"
        x_col = "x_m"
        title = r"$\Delta\rho_a^{zx}$ — joint gradient (zhang2021)"
    elif q == "spatial":
        df = rho_spatial_gradient(
            sites,
            spacing_m,
            comp=comp,
            recursive=recursive,
            on_dup=on_dup,
            strict=strict,
            verbose=verbose,
        )
        val_col = "delta_rho_x"
        x_col = "x_m"
        title = r"$\Delta\rho_a^x$ — spatial gradient (zhang2021)"
    else:
        df = rho_frequency_gradient(
            sites,
            comp=comp,
            spacing_m=spacing_m,
            recursive=recursive,
            on_dup=on_dup,
            strict=strict,
            verbose=verbose,
        )
        val_col = "delta_rho_z"
        x_col = "x_m"
        title = r"$\Delta\rho_a^z$ — frequency gradient (zhang2021)"

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    if df.empty:
        ax.text(
            0.5,
            0.5,
            "no data",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        ax.set_title(title)
        return ax

    y_col = "period_s" if period_axis else "freq_hz"

    all_x = np.sort(df[x_col].unique())
    all_y = np.sort(df[y_col].unique())

    x_idx = {float(v): k for k, v in enumerate(all_x)}
    y_idx = {float(v): k for k, v in enumerate(all_y)}

    grid = np.full((len(all_y), len(all_x)), np.nan)
    for _, row in df.iterrows():
        xi = x_idx.get(float(row[x_col]))
        yi = y_idx.get(float(row[y_col]))
        if xi is not None and yi is not None:
            grid[yi, xi] = row[val_col]

    valid = grid[np.isfinite(grid)]
    if vlim is None:
        vabs = float(np.nanmax(np.abs(valid))) if len(valid) else 1.0
        vabs = max(vabs, 1e-6)
        vlim = (-vabs, vabs)
    vmin, vmax = vlim

    try:
        norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=max(vmax, 1e-12))
    except Exception:
        norm = Normalize(vmin=vmin, vmax=vmax)

    X_idx, Y_idx = np.meshgrid(np.arange(len(all_x)), np.arange(len(all_y)))
    pcm = ax.pcolormesh(X_idx, Y_idx, grid, cmap=cmap, norm=norm, shading="nearest")
    plt.colorbar(pcm, ax=ax, label="Ω·m")

    # x-axis: station positions
    n_xtick = min(10, len(all_x))
    x_step = max(1, len(all_x) // n_xtick)
    xt_idx = np.arange(0, len(all_x), x_step)
    ax.set_xticks(xt_idx)
    ax.set_xticklabels(
        [f"{all_x[k]:.0f}" for k in xt_idx],
        rotation=45,
        ha="right",
        fontsize=8,
    )
    ax.set_xlabel("Position (m)")

    # y-axis: period or frequency
    n_ytick = min(8, len(all_y))
    y_step = max(1, len(all_y) // n_ytick)
    yt_idx = np.arange(0, len(all_y), y_step)
    ax.set_yticks(yt_idx)
    ax.set_yticklabels([f"{all_y[k]:.3g}" for k in yt_idx], fontsize=8)
    ax.set_ylabel("Period (s)" if period_axis else "Frequency (Hz)")

    ax.set_title(title, fontsize=10)
    return ax
