# pycsamt/emtools/spectra.py
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.emtools.spectra
=======================

Analysis and visualisation of MT cross-spectra stored in
:class:`~pycsamt.seg.spectra.Spectra` objects.

All plotting functions read visual defaults from the package-wide
API singletons at call time, so a single
:func:`~pycsamt.api.style.configure_style` or
:func:`~pycsamt.api.style.use_style` call propagates here too.

Analysis
--------
.. autosummary::

   coherence_matrix
   psd_table
   coherence_table
   snr_table
   band_select
   mask_low_coherence
   spectra_summary

Visualisation
-------------
.. autosummary::

   plot_psd
   plot_coherence
   plot_spectra_matrix
   plot_z_from_spectra
   plot_tipper_from_spectra
   plot_psd_section
   plot_coherence_section

Quick start
-----------
::

    from pycsamt.seg.spectra import Spectra
    from pycsamt.emtools.spectra import (
        plot_psd, plot_coherence, plot_z_from_spectra,
        plot_tipper_from_spectra, coherence_table,
    )

    sp = Spectra.from_file("site.edi")

    # power spectral density per channel
    plot_psd(sp, title="HBH05-03 PSD")

    # inter-channel coherence
    plot_coherence(sp)

    # apparent resistivity + phase recovered from spectra
    plot_z_from_spectra(sp)
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import (
    Any,
)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from ..api._rose_style import _UNSET
from ..api.control import PYCSAMT_CONTROL
from ..api.plot import add_colorbar
from ..api.section import PYCSAMT_SECTION, SectionStyle
from ..api.style import PYCSAMT_STYLE
from ..api.view import maybe_wrap_frame
from ._core import _axes_list

__all__ = [
    # analysis
    "coherence_matrix",
    "psd_table",
    "coherence_table",
    "snr_table",
    "band_select",
    "mask_low_coherence",
    "spectra_summary",
    # visualisation
    "plot_psd",
    "plot_coherence",
    "plot_spectra_matrix",
    "plot_z_from_spectra",
    "plot_tipper_from_spectra",
    "plot_psd_section",
    "plot_coherence_section",
]

# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────


def _x_vals(fr: np.ndarray) -> np.ndarray:
    return PYCSAMT_CONTROL.x.transform(fr)


def _x_label() -> str:
    return PYCSAMT_CONTROL.x.label()


def _x_log() -> bool:
    return PYCSAMT_CONTROL.x.use_log_scale()


def _spine_style(ax: Axes) -> None:
    ax.grid(True, which="both", ls=":", lw=0.4, color="0.75", zorder=0)
    ax.set_axisbelow(True)


def _check_spectra(sp: Any) -> None:
    if sp.n_freq == 0:
        raise ValueError("Spectra object has no frequency blocks.")
    if sp.n_chan == 0:
        raise ValueError("Spectra object has no channels.")


def _resolve_pairs(
    sp: Any,
    pairs: list[tuple[int, int]] | None,
) -> list[tuple[int, int]]:
    """Return (i, j) index pairs; default = all upper-triangle pairs."""
    nc = sp.n_chan
    if pairs is not None:
        return list(pairs)
    return [(i, j) for i in range(nc) for j in range(i + 1, nc)]


def _chan_label(sp: Any, idx: int) -> str:
    """Return a human-readable channel label for index *idx*."""
    ids = getattr(sp, "chan_ids", None)
    id_to_cht = getattr(sp, "id_to_chtype", {}) or {}
    if ids and idx < len(ids):
        raw = ids[idx]
        cht = id_to_cht.get(str(raw), "")
        return f"{cht}({raw})" if cht else str(raw)
    return f"ch{idx}"


def _sp_to_dict(
    sp_input: Any,
) -> dict[str, Any]:
    """Normalise single Spectra / list / dict to ``{name: Spectra}``."""
    # avoid importing Spectra at module level to prevent circular deps
    from ..seg.spectra import (
        Spectra as _Spectra,  # noqa: PLC0415
    )

    if isinstance(sp_input, _Spectra):
        name = getattr(sp_input, "name", None) or "site"
        return {str(name): sp_input}
    if isinstance(sp_input, dict):
        return {str(k): v for k, v in sp_input.items()}
    if isinstance(sp_input, (list, tuple)):
        out = {}
        for i, s in enumerate(sp_input):
            nm = getattr(s, "name", None) or f"site{i + 1}"
            out[str(nm)] = s
        return out
    raise TypeError(
        f"Expected Spectra, list[Spectra], or dict[str, Spectra]; "
        f"got {type(sp_input).__name__!r}"
    )


def _resolve_section_style(section: str | SectionStyle) -> SectionStyle:
    if isinstance(section, SectionStyle):
        return section.copy()
    return PYCSAMT_SECTION.style_for(str(section)).copy()


# ─────────────────────────────────────────────────────────────────────────────
# Analysis — returns arrays / DataFrames
# ─────────────────────────────────────────────────────────────────────────────


def coherence_matrix(sp: Any) -> np.ndarray:
    r"""Compute the inter-channel coherence matrix.

    The squared coherence between channels *i* and *j* is:

    .. math::

        \gamma^2_{ij}(f) = \frac{|S_{ij}(f)|^2}{S_{ii}(f)\,S_{jj}(f)}

    Parameters
    ----------
    sp : Spectra
        Cross-spectra container.

    Returns
    -------
    coh : ndarray, shape ``(n_freq, n_chan, n_chan)``
        Real-valued coherence matrix per frequency, values in [0, 1].
    """
    _check_spectra(sp)
    S = sp.S  # (nf, nc, nc)
    S_diag = np.real(np.diagonal(S, axis1=1, axis2=2))  # (nf, nc)
    denom = (
        S_diag[:, :, np.newaxis] * S_diag[:, np.newaxis, :]
    )  # (nf, nc, nc)
    with np.errstate(divide="ignore", invalid="ignore"):
        coh = np.abs(S) ** 2 / np.maximum(denom, 1e-40)
    coh = np.clip(np.real(coh), 0.0, 1.0)
    return coh


def psd_table(
    sp_input: Any,
    *,
    normalize: bool = False,
    api: bool | None = None,
) -> Any:
    """Power spectral density per channel as a tidy DataFrame.

    Parameters
    ----------
    sp_input : Spectra or list/dict of Spectra
        One or more cross-spectra containers.
    normalize : bool
        If ``True``, normalise each channel's PSD by its maximum value.

    Returns
    -------
    pd.DataFrame
        Columns: ``station``, ``freq``, ``period``, ``channel``, ``psd``.
    """
    rows = []
    for name, sp in _sp_to_dict(sp_input).items():
        _check_spectra(sp)
        S_diag = np.real(np.diagonal(sp.S, axis1=1, axis2=2))  # (nf, nc)
        for ch in range(sp.n_chan):
            psd = S_diag[:, ch].copy()
            if normalize:
                mx = np.nanmax(psd)
                if mx > 0:
                    psd = psd / mx
            lab = _chan_label(sp, ch)
            for fi, f in enumerate(sp.freq):
                rows.append(
                    {
                        "station": name,
                        "freq": float(f),
                        "period": float(1.0 / f) if f > 0 else np.nan,
                        "channel": lab,
                        "psd": float(psd[fi]),
                    }
                )
    df = pd.DataFrame(
        rows, columns=["station", "freq", "period", "channel", "psd"]
    )

    return maybe_wrap_frame(
        df,
        api=api,
        name="psd_table",
        kind="emtools.spectra.psd",
        source=sp_input,
        description="Power spectral density by station, frequency, and channel.",
    )


def coherence_table(
    sp_input: Any,
    *,
    pairs: list[tuple[int, int]] | None = None,
    api: bool | None = None,
) -> Any:
    """Inter-channel squared coherence as a tidy DataFrame.

    Parameters
    ----------
    sp_input : Spectra or list/dict of Spectra
    pairs : list of (i, j), optional
        Channel index pairs.  Default = all upper-triangle pairs.

    Returns
    -------
    pd.DataFrame
        Columns: ``station``, ``freq``, ``period``,
        ``ch_i``, ``ch_j``, ``pair``, ``coherence``.
    """
    rows = []
    for name, sp in _sp_to_dict(sp_input).items():
        _check_spectra(sp)
        coh = coherence_matrix(sp)  # (nf, nc, nc)
        prs = _resolve_pairs(sp, pairs)
        for i, j in prs:
            li = _chan_label(sp, i)
            lj = _chan_label(sp, j)
            label = f"{li}-{lj}"
            for fi, f in enumerate(sp.freq):
                rows.append(
                    {
                        "station": name,
                        "freq": float(f),
                        "period": float(1.0 / f) if f > 0 else np.nan,
                        "ch_i": li,
                        "ch_j": lj,
                        "pair": label,
                        "coherence": float(coh[fi, i, j]),
                    }
                )
    df = pd.DataFrame(
        rows,
        columns=[
            "station",
            "freq",
            "period",
            "ch_i",
            "ch_j",
            "pair",
            "coherence",
        ],
    )

    return maybe_wrap_frame(
        df,
        api=api,
        name="coherence_table",
        kind="emtools.spectra.coherence",
        source=sp_input,
        description="Squared coherence by station, frequency, and channel pair.",
    )


def snr_table(
    sp_input: Any,
    *,
    pairs: list[tuple[int, int]] | None = None,
    api: bool | None = None,
) -> Any:
    r"""Signal-to-noise ratio estimated from squared coherence.

    Uses the coherence-based estimator:

    .. math::

        \text{SNR} = \frac{\gamma^2}{1 - \gamma^2}

    with the dB version ``SNR_dB = 10 log₁₀(SNR)``.

    Parameters
    ----------
    sp_input : Spectra or list/dict of Spectra
    pairs : list of (i, j), optional

    Returns
    -------
    pd.DataFrame
        Columns: ``station``, ``freq``, ``period``, ``pair``,
        ``coherence``, ``snr``, ``snr_db``.
    """
    df = coherence_table(sp_input, pairs=pairs)
    g2 = df["coherence"].to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        snr = g2 / np.maximum(1.0 - g2, 1e-12)
        snr_db = 10.0 * np.log10(np.maximum(snr, 1e-12))
    df = df.copy()
    df["snr"] = snr
    df["snr_db"] = snr_db

    return maybe_wrap_frame(
        df,
        api=api,
        name="snr_table",
        kind="emtools.spectra.snr",
        source=sp_input,
        description="Coherence-derived signal-to-noise ratio table.",
    )


def band_select(
    sp: Any,
    f_min: float,
    f_max: float,
) -> Any:
    """Return a new :class:`Spectra` restricted to the frequency band
    ``[f_min, f_max]`` Hz.

    Parameters
    ----------
    sp : Spectra
    f_min, f_max : float
        Frequency limits in Hz (inclusive).

    Returns
    -------
    Spectra
        A shallow copy with arrays sliced to the band.
    """
    from ..seg.spectra import (
        Spectra as _Spectra,  # noqa: PLC0415
    )

    _check_spectra(sp)
    mask = (sp.freq >= float(f_min)) & (sp.freq <= float(f_max))
    if not mask.any():
        raise ValueError(
            f"No frequencies in [{f_min}, {f_max}] Hz.  "
            f"Available range: {sp.freq[-1]:.4g}–{sp.freq[0]:.4g} Hz."
        )
    idx = np.where(mask)[0]
    out = _Spectra(name=sp.name)
    out._freq = sp.freq[idx].copy()
    out._S = sp.S[idx].copy()
    out.bw = sp.bw[idx].copy()
    out.avgt = sp.avgt[idx].copy()
    out.avgf = sp.avgf[idx].copy()
    out.rotspec = sp.rotspec[idx].copy()
    out.segnum = sp.segnum[idx].copy()
    out.band = [sp.band[k] for k in idx.tolist()]
    out.chan_ids = list(sp.chan_ids)
    out.id_to_chtype = dict(sp.id_to_chtype)
    return out


def mask_low_coherence(
    sp: Any,
    *,
    pairs: list[tuple[int, int]] | None = None,
    threshold: float = 0.5,
    require_all: bool = False,
) -> np.ndarray:
    """Boolean mask of frequencies with *sufficient* coherence.

    Parameters
    ----------
    sp : Spectra
    pairs : list of (i, j), optional
        Pairs to evaluate.  Default = all upper-triangle pairs.
    threshold : float
        Minimum coherence required.  Default 0.5.
    require_all : bool
        If ``True``, all requested pairs must exceed *threshold*.
        If ``False`` (default), at least one pair suffices.

    Returns
    -------
    mask : ndarray of bool, shape ``(n_freq,)``
        ``True`` where coherence is acceptable.
    """
    _check_spectra(sp)
    coh = coherence_matrix(sp)  # (nf, nc, nc)
    prs = _resolve_pairs(sp, pairs)
    flags = np.zeros((sp.n_freq, len(prs)), dtype=bool)
    for k, (i, j) in enumerate(prs):
        flags[:, k] = coh[:, i, j] >= threshold
    if require_all:
        return flags.all(axis=1)
    return flags.any(axis=1)


def spectra_summary(sp: Any, *, api: bool | None = None) -> Any:
    """Compact per-frequency summary table.

    Columns: ``freq``, ``period``, ``bw``, ``avgt``, ``rotspec``,
    plus the diagonal PSD and mean off-diagonal coherence for each
    channel.

    Returns
    -------
    pd.DataFrame
    """
    _check_spectra(sp)
    coh = coherence_matrix(sp)  # (nf, nc, nc)
    psd = np.real(np.diagonal(sp.S, axis1=1, axis2=2))  # (nf, nc)

    rows = []
    for fi, f in enumerate(sp.freq):
        row: dict = {
            "freq": float(f),
            "period": float(1.0 / f) if f > 0 else np.nan,
            "bw": float(sp.bw[fi]),
            "avgt": float(sp.avgt[fi]),
            "rotspec": float(sp.rotspec[fi]),
        }
        for ch in range(sp.n_chan):
            lab = _chan_label(sp, ch)
            row[f"psd_{lab}"] = float(psd[fi, ch])
        # mean coherence over all pairs at this frequency
        nc = sp.n_chan
        if nc > 1:
            pairs_val = [
                coh[fi, i, j] for i in range(nc) for j in range(i + 1, nc)
            ]
            row["mean_coherence"] = float(np.nanmean(pairs_val))
        rows.append(row)

    df = pd.DataFrame(rows)

    return maybe_wrap_frame(
        df,
        api=api,
        name="spectra_summary",
        kind="emtools.spectra.summary",
        source=sp,
        description="Compact per-frequency spectra summary.",
    )


# ─────────────────────────────────────────────────────────────────────────────
# Visualisation — plot functions
# ─────────────────────────────────────────────────────────────────────────────


def plot_psd(
    sp: Any,
    *,
    channels: Sequence[int] | None = None,
    log_psd: bool = True,
    lw: float = _UNSET,
    alpha: float = _UNSET,
    title: str = "",
    figsize: tuple[float, float] = (9, 5),
    ax: Axes | None = None,
) -> Axes:
    """Plot the power spectral density per channel.

    The x-axis follows :data:`~pycsamt.api.control.PYCSAMT_CONTROL`
    (log₁₀ period by default).  Colours come from
    :attr:`~pycsamt.api.style.PYCSAMT_STYLE.multiline`.

    Parameters
    ----------
    sp : Spectra
    channels : sequence of int, optional
        Channel indices to plot.  Default = all channels.
    log_psd : bool
        Display log₁₀(PSD) when ``True``.
    lw : float
        Line width.  Default: ``PYCSAMT_STYLE.multiline.lw``.
    alpha : float
        Line alpha.  Default: ``PYCSAMT_STYLE.multiline.alpha``.
    title : str
    figsize : (float, float)
    ax : Axes or None

    Returns
    -------
    ax : Axes
    """
    _check_spectra(sp)
    _ml = PYCSAMT_STYLE.multiline
    if lw is _UNSET:
        lw = _ml.lw
    if alpha is _UNSET:
        alpha = _ml.alpha

    nc = sp.n_chan
    chs = list(channels) if channels is not None else list(range(nc))
    cols = _ml.colors(len(chs))
    psd = np.real(np.diagonal(sp.S, axis1=1, axis2=2))  # (nf, nc)
    x = _x_vals(sp.freq)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    for ki, chi in enumerate(chs):
        y = psd[:, chi]
        if log_psd:
            y = np.log10(np.maximum(y, 1e-40))
        lab = _chan_label(sp, chi)
        ax.plot(x, y, color=cols[ki], lw=lw, alpha=alpha, label=lab)

    if _x_log():
        ax.set_xscale("log")
    ax.set_xlabel(_x_label(), fontsize=9)
    ylabel = r"$\log_{10}$PSD" if log_psd else "PSD"
    ax.set_ylabel(ylabel, fontsize=9)
    ax.legend(fontsize=8, framealpha=0.8)
    ax.set_title(
        title or f"Power spectral density  —  {sp.name or 'site'}",
        fontsize=10,
        pad=6,
    )
    _spine_style(ax)
    return ax


def plot_coherence(
    sp: Any,
    *,
    pairs: list[tuple[int, int]] | None = None,
    threshold: float = 0.5,
    show_threshold: bool = True,
    lw: float = _UNSET,
    alpha: float = _UNSET,
    title: str = "",
    axes=None,
    figsize: tuple[float, float] | None = None,
) -> np.ndarray:
    """Plot squared coherence for selected channel pairs.

    Each pair gets its own sub-axis arranged in a grid.  A dashed
    horizontal line marks *threshold*.

    Parameters
    ----------
    sp : Spectra
    pairs : list of (i, j), optional
        Default = all upper-triangle pairs.
    threshold : float
        Quality threshold drawn as a dashed line.  Default 0.5.
    show_threshold : bool
    lw, alpha : float
        Style defaults from ``PYCSAMT_STYLE.multiline``.
    title : str
    figsize : (float, float) or None
        Auto-computed from the number of pairs when ``None``.

    Returns
    -------
    axes : ndarray of Axes, shape ``(n_pairs,)``
    """
    _check_spectra(sp)
    _ml = PYCSAMT_STYLE.multiline
    if lw is _UNSET:
        lw = _ml.lw
    if alpha is _UNSET:
        alpha = _ml.alpha

    coh = coherence_matrix(sp)
    prs = _resolve_pairs(sp, pairs)
    n = len(prs)
    if n == 0:
        raise ValueError("No channel pairs available.")

    ncols = min(n, 3)
    nrows = (n + ncols - 1) // ncols
    if figsize is None:
        figsize = (4.0 * ncols, 3.5 * nrows)

    axes_given = _axes_list(axes, n) if axes is not None else None
    if axes_given is None:
        fig, axs_raw = plt.subplots(
            nrows, ncols, figsize=figsize, constrained_layout=True
        )
        axs = np.asarray(axs_raw).ravel()
    else:
        axs = np.asarray(axes_given, dtype=object)
        fig = axs[0].figure
    cols = _ml.colors(n)
    x = _x_vals(sp.freq)

    for ki, (i, j) in enumerate(prs):
        ax = axs[ki]
        li = _chan_label(sp, i)
        lj = _chan_label(sp, j)
        y = coh[:, i, j]
        ax.plot(x, y, color=cols[ki], lw=lw, alpha=alpha, label=f"{li}-{lj}")
        if show_threshold:
            ax.axhline(
                threshold,
                ls="--",
                lw=0.8,
                color="0.5",
                label=f"thr={threshold}",
            )
        ax.set_ylim(-0.05, 1.05)
        if _x_log():
            ax.set_xscale("log")
        ax.set_xlabel(_x_label(), fontsize=8)
        ax.set_ylabel(r"$\gamma^2$", fontsize=8)
        ax.set_title(f"{li} – {lj}", fontsize=9, pad=4)
        ax.legend(fontsize=7, framealpha=0.7)
        _spine_style(ax)

    # hide surplus axes
    for ak in axs[n:]:
        ak.set_visible(False)

    fig.suptitle(
        title or f"Coherence  —  {sp.name or 'site'}",
        fontsize=10,
        y=1.01,
    )
    return axs[:n]


def plot_spectra_matrix(
    sp: Any,
    *,
    freq_idx: int = 0,
    quantity: str = "abs",
    cmap: str = _UNSET,
    log_scale: bool = True,
    title: str = "",
    ax=None,
    figsize: tuple[float, float] = (7, 6),
) -> Figure:
    """Visualise the full cross-spectral density matrix at one frequency.

    The matrix is drawn as a colour image.  The diagonal shows auto-
    spectra (PSD); the upper and lower triangles show the magnitude (or
    real/imaginary part) of the cross-spectra.

    Parameters
    ----------
    sp : Spectra
    freq_idx : int
        Index into ``sp.freq`` for the frequency slice.
    quantity : {'abs', 'real', 'imag', 'phase'}
        Quantity to colour.
    cmap : str or _UNSET
        Colour map.  Defaults to ``"viridis"`` (abs) or ``"RdBu_r"``
        (real/imag/phase).
    log_scale : bool
        Apply log₁₀ to the absolute value when ``quantity="abs"``.
    title : str
    figsize : (float, float)

    Returns
    -------
    fig : Figure
    """
    _check_spectra(sp)
    nc = sp.n_chan
    M = sp.S[freq_idx]  # (nc, nc) complex
    freq = sp.freq[freq_idx]

    if quantity == "abs":
        data = np.abs(M)
        if log_scale:
            data = np.log10(np.maximum(data, 1e-40))
            cb_label = r"$\log_{10}|S_{ij}|$"
        else:
            cb_label = r"$|S_{ij}|$"
        if cmap is _UNSET:
            cmap = "viridis"
    elif quantity == "real":
        data = np.real(M)
        cb_label = r"Re($S_{ij}$)"
        if cmap is _UNSET:
            cmap = "RdBu_r"
    elif quantity == "imag":
        data = np.imag(M)
        cb_label = r"Im($S_{ij}$)"
        if cmap is _UNSET:
            cmap = "RdBu_r"
    elif quantity == "phase":
        data = np.angle(M, deg=True)
        cb_label = r"Phase($S_{ij}$)  (°)"
        if cmap is _UNSET:
            cmap = "hsv"
    else:
        raise ValueError(
            f"quantity must be 'abs','real','imag','phase'; got {quantity!r}"
        )

    labels = [_chan_label(sp, k) for k in range(nc)]

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    else:
        fig = ax.figure
    im = ax.imshow(data, cmap=cmap, aspect="equal", origin="upper")
    add_colorbar(im, ax, label=cb_label, size="4%", pad=0.04)

    ax.set_xticks(range(nc))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_yticks(range(nc))
    ax.set_yticklabels(labels, fontsize=8)

    # annotate cells
    for r in range(nc):
        for c in range(nc):
            ax.text(
                c,
                r,
                f"{data[r, c]:.2g}",
                ha="center",
                va="center",
                fontsize=6.5,
                color="0.2",
            )

    ax.set_title(
        title
        or (
            f"Spectral matrix  —  {sp.name or 'site'}  "
            f"f = {freq:.4g} Hz  [{quantity}]"
        ),
        fontsize=10,
        pad=6,
    )
    return fig


def plot_z_from_spectra(
    sp: Any,
    *,
    e_labels: tuple[str, str] = ("EX", "EY"),
    h_labels: tuple[str, str] = ("HX", "HY"),
    ridge: float | None = None,
    estimate_error: bool = False,
    show_error: bool = True,
    title: str = "",
    axes=None,
    figsize: tuple[float, float] = (10, 5),
) -> Figure:
    """Plot apparent resistivity and phase recovered from spectra.

    Calls :meth:`~pycsamt.seg.spectra.Spectra.to_Z` internally and
    renders the result with the standard MT component styling from
    :attr:`~pycsamt.api.style.PYCSAMT_STYLE`.

    Parameters
    ----------
    sp : Spectra
    e_labels, h_labels : tuple of str
        Channel type labels used for the E and H blocks in
        :meth:`~pycsamt.seg.spectra.Spectra.to_Z`.
    ridge : float or None
        Tikhonov regularization for S_HH inversion.
    estimate_error : bool
        Estimate 1-σ errors in :meth:`to_Z`.  Default ``False``
        (avoids DoF warnings when metadata is incomplete).
    show_error : bool
        Shade error envelope when errors are available.
    title : str
    figsize : (float, float)

    Returns
    -------
    fig : Figure
    """
    _check_spectra(sp)
    z_obj, _ = sp.to_Z(
        e_labels=e_labels,
        h_labels=h_labels,
        ridge=ridge,
        estimate_error=estimate_error,
    )

    _st = PYCSAMT_STYLE.mt
    freqs = z_obj.freq
    x = _x_vals(freqs)
    rho = z_obj.resistivity  # (nf, 2, 2)
    phi = z_obj.phase  # (nf, 2, 2)
    z_err = z_obj.z_err  # None or (nf, 2, 2)

    axes_given = _axes_list(axes, 2)
    if axes_given is None:
        fig, (ax_r, ax_p) = plt.subplots(
            1, 2, figsize=figsize, constrained_layout=True
        )
    else:
        ax_r, ax_p = axes_given
        fig = ax_r.figure

    pairs = [
        ("xy", (0, 1), _st.xy),
        ("yx", (1, 0), _st.yx),
    ]
    for _comp, (r, c), sty in pairs:
        rho_c = np.log10(np.maximum(rho[:, r, c], 1e-12))
        phi_c = phi[:, r, c]
        kw = sty.plot_kwargs()
        ax_r.plot(x, rho_c, **kw)
        ax_p.plot(x, phi_c, **kw)

        if show_error and z_err is not None:
            # propagate |Z| error to log-rho error ≈ 2Δ|Z|/|Z| / ln10
            z_c = z_obj.z[:, r, c]
            z_err_c = z_err[:, r, c]
            rel = np.abs(z_err_c) / (np.abs(z_c) + 1e-24)
            drho = 2.0 * rel / np.log(10.0)
            ax_r.fill_between(
                x,
                rho_c - drho,
                rho_c + drho,
                color=sty.color,
                alpha=0.15,
                linewidth=0,
            )

    for ax in (ax_r, ax_p):
        _spine_style(ax)
        if _x_log():
            ax.set_xscale("log")
        ax.set_xlabel(_x_label(), fontsize=9)

    ax_r.set_ylabel(r"$\log_{10}\rho_a$  ($\Omega\cdot$m)", fontsize=9)
    ax_r.set_title(r"Apparent resistivity  $\rho_a$", fontsize=9, pad=6)
    ax_p.set_ylabel(r"Phase  (°)", fontsize=9)
    ax_p.set_title("Impedance phase", fontsize=9, pad=6)

    # shared legend on rho panel
    ax_r.legend(fontsize=8, framealpha=0.8)

    fig.suptitle(
        title or f"Z from spectra  —  {sp.name or 'site'}",
        fontsize=11,
        y=1.02,
    )
    return fig


def plot_tipper_from_spectra(
    sp: Any,
    *,
    h_labels: tuple[str, str] = ("HX", "HY"),
    ridge: float | None = None,
    estimate_error: bool = False,
    show_error: bool = True,
    title: str = "",
    axes=None,
    figsize: tuple[float, float] = (10, 5),
) -> np.ndarray:
    """Plot the induction tipper magnitude and phase from spectra.

    Displays the real and imaginary parts of T_x and T_y as well as
    their magnitudes on a two-panel figure (amplitude | phase).

    Parameters
    ----------
    sp : Spectra
    h_labels : tuple of str
    ridge : float or None
    estimate_error : bool
    show_error : bool
    title : str
    figsize : (float, float)

    Returns
    -------
    axes : ndarray of Axes, shape (2,)
        ``[ax_amp, ax_phase]``
    """
    _check_spectra(sp)
    _, tip = sp.to_Z(
        h_labels=h_labels,
        ridge=ridge,
        estimate_error=estimate_error,
    )

    if tip is None:
        axes_given = _axes_list(axes, 2) if axes is not None else None
        if axes_given is None:
            _, axs = plt.subplots(
                1, 2, figsize=figsize, constrained_layout=True
            )
        else:
            axs = axes_given
        for ax in axs:
            ax.text(
                0.5,
                0.5,
                "No tipper (HZ not found)",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=10,
                color="0.5",
            )
            _spine_style(ax)
        return np.array(axs)

    x = _x_vals(tip.freq)
    T = tip.tipper[:, 0, :]  # (nf, 2)  — Tx, Ty
    T_err = tip.tipper_err[:, 0, :] if tip.tipper_err is not None else None

    axes_given = _axes_list(axes, 2)
    if axes_given is None:
        fig, (ax_a, ax_p) = plt.subplots(
            1, 2, figsize=figsize, constrained_layout=True
        )
    else:
        ax_a, ax_p = axes_given
        fig = ax_a.figure

    _st = PYCSAMT_STYLE.mt
    specs = [
        ("T_x", T[:, 0], _st.xy),
        ("T_y", T[:, 1], _st.yx),
    ]
    for lab, Tc, sty in specs:
        amp = np.abs(Tc)
        phase = np.angle(Tc, deg=True)
        kw = sty.plot_kwargs(label=lab)
        ax_a.plot(x, amp, **kw)
        ax_p.plot(x, phase, **kw)

        if show_error and T_err is not None:
            idx = 0 if lab == "T_x" else 1
            err = np.abs(T_err[:, idx])
            ax_a.fill_between(
                x,
                amp - err,
                amp + err,
                color=sty.color,
                alpha=0.15,
                linewidth=0,
            )

    for ax in (ax_a, ax_p):
        _spine_style(ax)
        if _x_log():
            ax.set_xscale("log")
        ax.set_xlabel(_x_label(), fontsize=9)
        ax.legend(fontsize=8, framealpha=0.8)

    ax_a.set_ylabel("|T|", fontsize=9)
    ax_a.set_title("Tipper magnitude", fontsize=9, pad=6)
    ax_p.set_ylabel("Phase  (°)", fontsize=9)
    ax_p.set_title("Tipper phase", fontsize=9, pad=6)

    fig.suptitle(
        title or f"Tipper from spectra  —  {sp.name or 'site'}",
        fontsize=11,
        y=1.02,
    )
    return np.array([ax_a, ax_p])


def plot_psd_section(
    sp_input: Any,
    *,
    channel: int = 0,
    log_psd: bool = True,
    cmap: str = "viridis",
    vmin: float | None = None,
    vmax: float | None = None,
    section: str | SectionStyle = "pseudosection",
    title: str = "",
    figsize: tuple[float, float] | None = None,
    ax: Axes | None = None,
) -> Axes:
    """Pseudo-section of PSD across stations (station × period).

    Interpolates all Spectra objects to a common log-spaced frequency
    grid before assembling the 2-D colour map.

    Parameters
    ----------
    sp_input : Spectra or list/dict of Spectra
    channel : int
        Channel index to display.  Default 0.
    log_psd : bool
        Colour log₁₀(PSD) when ``True``.
    cmap : str
    vmin, vmax : float or None
    section : str or SectionStyle
        Layout preset from :data:`~pycsamt.api.section.PYCSAMT_SECTION`.
    title : str
    figsize : (float, float) or None
    ax : Axes or None

    Returns
    -------
    ax : Axes
    """
    sp_dict = _sp_to_dict(sp_input)
    names = list(sp_dict.keys())
    sps = list(sp_dict.values())

    # common log-frequency grid (intersection)
    all_freqs = [s.freq for s in sps]
    f_min = max(float(f.min()) for f in all_freqs)
    f_max = min(float(f.max()) for f in all_freqs)
    n_f = max(len(f) for f in all_freqs)
    f_grid = np.logspace(np.log10(f_min), np.log10(f_max), n_f)

    # build PSD matrix (n_stations, n_freqs)
    matrix = np.full((len(sps), n_f), np.nan)
    for si, sp in enumerate(sps):
        psd_raw = np.real(np.diagonal(sp.S, axis1=1, axis2=2))[:, channel]
        matrix[si] = np.interp(
            np.log10(f_grid),
            np.log10(sp.freq[::-1]),
            psd_raw[::-1],
        )

    data = np.log10(np.maximum(matrix, 1e-40)) if log_psd else matrix
    y = np.log10(1.0 / f_grid)  # log10(period)

    sty = _resolve_section_style(section)
    if figsize is None:
        figsize = sty.figsize_for(n_stations=len(sps), n_y=n_f)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    else:
        ax.get_figure()

    # pixel edges
    st_x = np.arange(len(sps), dtype=float)
    dx_half = 0.5
    x_edges = np.r_[st_x[0] - dx_half, st_x + dx_half]
    if len(y) > 1:
        dy = np.abs(np.diff(y)) / 2.0
        sgn = np.sign(np.diff(y))
        y_edges = np.r_[
            y[0] - dy[0], y[:-1] + sgn * dy, y[-1] + sgn[-1] * dy[-1]
        ]
    else:
        y_edges = np.r_[y[0] - 0.2, y[0] + 0.2]

    pc = ax.pcolormesh(
        x_edges,
        y_edges,
        data.T,
        cmap=cmap,
        shading="flat",
        vmin=vmin,
        vmax=vmax,
    )
    cb_label = r"$\log_{10}$PSD" if log_psd else "PSD"
    sty.add_colorbar(pc, ax, label=cb_label)

    sty.apply_axis(ax, xlabel="Station", ylabel=r"$\log_{10}T$ (s)")
    if y[0] > y[-1]:
        ax.invert_yaxis()

    sty.apply_stations(ax, st_x, names)

    chan_lab = _chan_label(sps[0], channel)
    ax.set_title(
        title or f"PSD pseudo-section  —  channel {chan_lab}",
        fontsize=10,
        pad=6,
    )
    _spine_style(ax)
    return ax


def plot_coherence_section(
    sp_input: Any,
    *,
    pair: tuple[int, int] | None = None,
    threshold: float = 0.5,
    show_threshold: bool = True,
    cmap: str = "RdYlGn",
    section: str | SectionStyle = "pseudosection",
    title: str = "",
    figsize: tuple[float, float] | None = None,
    ax: Axes | None = None,
) -> Axes:
    """Pseudo-section of coherence across stations (station × period).

    Parameters
    ----------
    sp_input : Spectra or list/dict of Spectra
    pair : (int, int) or None
        Single channel pair ``(i, j)`` to display.  When ``None`` the
        mean over all upper-triangle pairs is shown.
    threshold : float
        Value shown by the shared colorbar.  Default 0.5.
    show_threshold : bool
        Add a contour at *threshold* when ``True``.
    cmap : str
        Default ``"RdYlGn"`` — red=low, green=high coherence.
    section : str or SectionStyle
    title : str
    figsize : (float, float) or None
    ax : Axes or None

    Returns
    -------
    ax : Axes
    """
    sp_dict = _sp_to_dict(sp_input)
    names = list(sp_dict.keys())
    sps = list(sp_dict.values())

    # common grid
    all_freqs = [s.freq for s in sps]
    f_min = max(float(f.min()) for f in all_freqs)
    f_max = min(float(f.max()) for f in all_freqs)
    n_f = max(len(f) for f in all_freqs)
    f_grid = np.logspace(np.log10(f_min), np.log10(f_max), n_f)

    matrix = np.full((len(sps), n_f), np.nan)
    for si, sp in enumerate(sps):
        coh_m = coherence_matrix(sp)  # (nf, nc, nc)
        nc = sp.n_chan
        if pair is not None:
            i, j = pair
            coh_1d = coh_m[:, i, j]
        else:
            prs = [(i, j) for i in range(nc) for j in range(i + 1, nc)]
            coh_1d = np.nanmean([coh_m[:, i, j] for i, j in prs], axis=0)

        matrix[si] = np.interp(
            np.log10(f_grid),
            np.log10(sp.freq[::-1]),
            coh_1d[::-1],
        )

    y = np.log10(1.0 / f_grid)
    sty = _resolve_section_style(section)
    if figsize is None:
        figsize = sty.figsize_for(n_stations=len(sps), n_y=n_f)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize, constrained_layout=True)
        ax.get_figure()
    else:
        ax.get_figure()

    st_x = np.arange(len(sps), dtype=float)
    dx_half = 0.5
    x_edges = np.r_[st_x[0] - dx_half, st_x + dx_half]
    if len(y) > 1:
        dy = np.abs(np.diff(y)) / 2.0
        sgn = np.sign(np.diff(y))
        y_edges = np.r_[
            y[0] - dy[0], y[:-1] + sgn * dy, y[-1] + sgn[-1] * dy[-1]
        ]
    else:
        y_edges = np.r_[y[0] - 0.2, y[0] + 0.2]

    pc = ax.pcolormesh(
        x_edges,
        y_edges,
        matrix.T,
        cmap=cmap,
        shading="flat",
        vmin=0.0,
        vmax=1.0,
    )
    sty.add_colorbar(pc, ax, label=r"$\gamma^2$")

    if show_threshold and matrix.shape[1] > 1:
        ax.contour(
            st_x,
            y,
            matrix.T,
            levels=[threshold],
            colors=["k"],
            linewidths=0.8,
            linestyles=["--"],
        )

    sty.apply_axis(ax, xlabel="Station", ylabel=r"$\log_{10}T$ (s)")
    if y[0] > y[-1]:
        ax.invert_yaxis()

    sty.apply_stations(ax, st_x, names)

    pair_str = (
        f"ch{pair[0]}-ch{pair[1]}" if pair is not None else "mean all pairs"
    )
    ax.set_title(
        title or f"Coherence pseudo-section  —  {pair_str}",
        fontsize=10,
        pad=6,
    )
    _spine_style(ax)
    return ax
