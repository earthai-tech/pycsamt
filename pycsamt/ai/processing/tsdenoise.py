# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
TimeSeriesDenoiser -- MMF-SVM-K-SVD raw field-time-series denoising.

Every other estimator in :mod:`pycsamt.ai.processing` works on the
*frequency-domain* impedance tensor (:class:`~pycsamt.ai.processing.\
denoise.EMDenoiser`) or on station-level tables derived from it. This
module works one stage further upstream: on the raw
:class:`~pycsamt.ts.TSData` sample stream (Ex, Ey, Hx, Hy, Hz), the
same object :func:`pycsamt.ts.readers.read_ts` produces and
:func:`pycsamt.ts.process.ts_to_spectra` FFTs into cross-power spectra.
It is meant to run *before* :func:`pycsamt.ts.process.preprocess`, as
a strong-interference-suppression pass on the field record itself,
complementing rather than duplicating that function's light gap-fill
and the Huber-robust segment weighting already built into
:func:`~pycsamt.ts.process.cross_spectra`.

The method is a direct implementation of Gui et al. (2024) [1]_ --
mathematical-morphological filtering (MMF), a support-vector-machine
(SVM) signal/noise classifier, and K-SVD dictionary-learning
denoising -- extending the authors' earlier MMF-K-SVD method
(Gui et al., 2021) [2]_ with the SVM triage stage that prevents
high-quality low-amplitude MT signal from being treated as noise.

Pipeline
--------
1. :func:`mmf_split` separates each channel into a low-frequency part
   (kept untouched -- MT signal below ~1 Hz is exactly what a
   mathematical morphological filter with a sub-second structuring
   element preserves) and a high-frequency residual that may carry
   strong interference.
2. The high-frequency residual is cut into ``win_seconds`` windows;
   :class:`SignalQualityClassifier` (a linear-kernel SVM trained on
   four complexity features -- sample entropy, fuzzy entropy,
   approximate entropy, box-counting fractal dimension, see
   :func:`compute_entropy_features`) labels each window clean or
   noisy, mirroring the paper's Fig. 4-6.
3. Contiguous runs of noisy windows are handed to
   :class:`KSVDDenoiser`, which learns a small dictionary from
   overlapping patches of the run itself (self-supervised -- no
   pre-training needed) and subtracts the atoms' sparse
   reconstruction, which -- because transient interference is far
   easier to sparsely represent than the genuinely random background
   signal -- approximates the noise profile, not the signal.
4. Low-frequency part, untouched clean windows, and denoised noisy
   runs are merged back into one channel.

:class:`TimeSeriesDenoiser` orchestrates all three stages and accepts
/ returns a :class:`~pycsamt.ts.TSData` record directly via
:meth:`TimeSeriesDenoiser.apply`, the same TSData-in / TSData-out
convention :meth:`~pycsamt.ai.processing.denoise.EMDenoiser.apply`
uses for site collections.

References
----------
.. [1] Gui, T., Deng, J., Li, G., Chen, H., Yu, H., Feng, M., 2024.
       De-noising magnetotelluric data based on machine learning.
       *J. Appl. Geophys.* 230, 105538.
.. [2] Gui, T., Deng, J., Li, G., Liu, X., Chen, H., He, Z., 2021.
       De-noising magnetotelluric data based on mathematical
       morphology and K-SVD dictionary learning. *Chin. J. Nonferrous
       Met.* 31(12), 3713-3729 (in Chinese).
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from .._base import BaseEMProcessor

__all__ = [
    "mmf_split",
    "sample_entropy",
    "fuzzy_entropy",
    "approx_entropy",
    "box_dimension",
    "compute_entropy_features",
    "generate_synthetic_library",
    "SignalQualityClassifier",
    "omp",
    "omp_batch",
    "ksvd_dictionary",
    "KSVDDenoiser",
    "TimeSeriesDenoiser",
    "snr_db",
    "time_domain_ncc",
]

_FEATURE_NAMES = ("sample_entropy", "fuzzy_entropy", "approx_entropy",
                   "box_dimension")


# ─────────────────────────────────────────────────────────────────────────────
# NaN handling (self-contained -- no pycsamt.ts import at array level)
# ─────────────────────────────────────────────────────────────────────────────


def _fill_nan_linear(x: np.ndarray) -> np.ndarray:
    """Linear-interpolate interior NaNs; edge-fill leading/trailing runs."""
    x = np.asarray(x, dtype=float).copy()
    isnan = np.isnan(x)
    if not isnan.any():
        return x
    finite = np.flatnonzero(~isnan)
    if finite.size == 0:
        return np.zeros_like(x)
    x[: finite[0]] = x[finite[0]]
    x[finite[-1] + 1 :] = x[finite[-1]]
    idx = np.arange(x.size)
    still = np.isnan(x)
    x[still] = np.interp(idx[still], idx[~still], x[~still])
    return x


# ─────────────────────────────────────────────────────────────────────────────
# 1. Mathematical morphological filtering (MMF)
# ─────────────────────────────────────────────────────────────────────────────


def mmf_split(
    x: np.ndarray,
    *,
    size: int | None = None,
    dt: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    r"""
    Split a time series into a low-frequency part and a high-frequency
    residual using mathematical morphological filtering (MMF).

    Follows Tang et al. (2012) / Li et al. (2020): the low-frequency
    envelope is the average of the "open-then-close" and
    "close-then-open" morphological filters with a flat 1-D
    structuring element of length *size*,

    .. math::

        x_{low} = \tfrac{1}{2}\left[
            \mathrm{OC}(x) + \mathrm{CO}(x)
        \right], \qquad
        \mathrm{OC}(x) = \delta_g(\epsilon_g(x)),\;\;
        \mathrm{CO}(x) = \epsilon_g(\delta_g(x))

    where :math:`\epsilon_g`/:math:`\delta_g` are grayscale
    erosion/dilation. Averaging the two orderings removes both
    positive and negative transients symmetrically -- a spike that
    survives opening-then-closing is suppressed by closing-then-opening
    and vice versa -- which is what makes MMF, unlike CEEMD, VMD, or
    wavelet decomposition, remain reliable at extracting the
    low-frequency signal *even when* strong interference is present
    (see the module docstring and Fig. 7 of Gui et al., 2024).

    Parameters
    ----------
    x : ndarray, shape (n,)
        Input channel. ``NaN`` samples are linearly interpolated
        before filtering and restored as ``NaN`` in both outputs.
    size : int or None
        Structuring-element length in samples (odd, >= 3). When
        ``None``, defaults to the number of samples spanning
        approximately 2 s (rounded to the nearest odd integer, using
        *dt*) -- a spike/transient shorter than this is treated as
        "high frequency"; MT signal below roughly 1 Hz survives as
        "low frequency".

        This is the parameter that matters most, and the paper itself
        treats it as data-dependent (Li et al., 2020's own "detailed
        principles ... on selecting the corresponding structural
        element sizes"). *size* must exceed the widest interference
        pulse the record contains, or that pulse's plateau survives
        erosion/dilation and leaks straight into *low* instead of
        being rejected into *high* -- e.g. on a real 15 Hz field
        record with ~3-4 s wide square-wave bursts, the 2 s default
        left visible spikes in *low* and only cut noisy-window SNR by
        a few dB; raising *size* to ~8 s (121 samples) rejected the
        bursts cleanly and raised the same window's SNR by ~17 dB.
        Inspect ``low`` against a plot of *x* and widen *size* until
        the plateaus/steps disappear from it.
    dt : float or None
        Sampling interval (s), used only to derive the default *size*.
        Required when *size* is ``None``.

    Returns
    -------
    low, high : ndarray, shape (n,)
        ``low + high == x`` (up to floating point). *low* is the
        MMF envelope; *high* is the residual carrying any strong
        interference plus genuine high-frequency signal.

    See Also
    --------
    TimeSeriesDenoiser
        Uses this as its first stage.
    """
    from scipy.ndimage import grey_closing, grey_dilation, grey_erosion

    x_arr = np.asarray(x, dtype=float)
    nanmask = np.isnan(x_arr)
    xf = _fill_nan_linear(x_arr)

    if size is None:
        if dt is None:
            raise ValueError(
                "mmf_split: provide either `size` (samples) or `dt` "
                "(seconds) so a ~1s default structuring element can be "
                "derived."
            )
        size = int(round(2.0 / float(dt)))
    size = max(3, int(size))
    size = min(size, xf.size - 1 if xf.size > 3 else 3)
    size = max(3, size)
    if size % 2 == 0:
        size -= 1
    size = max(3, size)

    footprint = np.ones(size, dtype=bool)
    opening = grey_dilation(grey_erosion(xf, footprint=footprint),
                             footprint=footprint)
    closing = grey_erosion(grey_dilation(xf, footprint=footprint),
                            footprint=footprint)
    oc = grey_closing(opening, footprint=footprint)
    co = grey_dilation(grey_erosion(closing, footprint=footprint),
                        footprint=footprint)
    low = 0.5 * (oc + co)
    high = xf - low

    low = np.where(nanmask, np.nan, low)
    high = np.where(nanmask, np.nan, high)
    return low, high


# ─────────────────────────────────────────────────────────────────────────────
# 2. Complexity / entropy features (SE, FE, AE, BD)
# ─────────────────────────────────────────────────────────────────────────────


def _embed(x: np.ndarray, m: int) -> np.ndarray:
    """Time-delay embedding, shape (n - m + 1, m)."""
    n = x.size
    return np.stack([x[i : n - m + 1 + i] for i in range(m)], axis=1)


def _chebyshev_dist(vecs: np.ndarray) -> np.ndarray:
    """Pairwise Chebyshev (max-norm) distance matrix."""
    diff = np.abs(vecs[:, None, :] - vecs[None, :, :])
    return diff.max(axis=-1)


def sample_entropy(x: np.ndarray, *, m: int = 2, r: float | None = None
                    ) -> float:
    """
    Sample entropy (Richman & Moorman, 2000).

    Larger for more irregular, less self-similar signals -- the
    high-quality MT background -- and small (near 0) for the
    stereotyped, self-similar shapes of strong-interference
    transients (Fig. 4a of Gui et al., 2024).

    Parameters
    ----------
    x : ndarray, shape (n,)
    m : int, default 2
        Embedding dimension.
    r : float or None
        Tolerance. Defaults to ``0.2 * std(x)`` (Pincus's convention),
        the value used throughout the paper.

    Returns
    -------
    float
        ``np.nan`` when the window is too short or degenerate
        (``std(x) == 0``) to estimate.
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    if n < m + 2:
        return np.nan
    sd = np.std(x)
    if not np.isfinite(sd) or sd == 0:
        return np.nan
    tol = float(r) if r is not None else 0.2 * sd

    def _phi_count(mm: int) -> float:
        vecs = _embed(x, mm)
        d = _chebyshev_dist(vecs)
        np.fill_diagonal(d, np.inf)
        return float((d <= tol).sum())

    b = _phi_count(m)
    a = _phi_count(m + 1)
    if b == 0 or a == 0:
        return np.nan
    return float(-np.log(a / b))


def approx_entropy(x: np.ndarray, *, m: int = 2, r: float | None = None
                    ) -> float:
    """
    Approximate entropy (Pincus, 1995).

    Same interpretation as :func:`sample_entropy` (larger for
    high-quality signal, Fig. 4c), but counts self-matches, which
    biases it slightly and makes it less consistent for short windows
    -- kept for parity with the paper's four-feature vector.

    Parameters
    ----------
    x : ndarray, shape (n,)
    m : int, default 2
    r : float or None
        Defaults to ``0.2 * std(x)``.

    Returns
    -------
    float
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    if n < m + 2:
        return np.nan
    sd = np.std(x)
    if not np.isfinite(sd) or sd == 0:
        return np.nan
    tol = float(r) if r is not None else 0.2 * sd

    def _phi(mm: int) -> float:
        vecs = _embed(x, mm)
        nv = vecs.shape[0]
        d = _chebyshev_dist(vecs)
        counts = (d <= tol).sum(axis=1)  # includes self-match (diag==0)
        counts = np.maximum(counts, 1)
        return float(np.mean(np.log(counts / nv)))

    return float(_phi(m) - _phi(m + 1))


def fuzzy_entropy(
    x: np.ndarray, *, m: int = 2, r: float | None = None, n_exp: float = 2.0
) -> float:
    """
    Fuzzy entropy (Chen et al., 2007).

    Replaces sample entropy's hard ``d <= r`` match with a smooth
    exponential membership :math:`\\exp(-(d/r)^{n})`, making it less
    sensitive to the tolerance choice; same clean-vs-noise separation
    as :func:`sample_entropy` (Fig. 4b).

    Parameters
    ----------
    x : ndarray, shape (n,)
    m : int, default 2
    r : float or None
        Defaults to ``0.2 * std(x)``.
    n_exp : float, default 2.0
        Membership-function steepness.

    Returns
    -------
    float
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    if n < m + 2:
        return np.nan
    sd = np.std(x)
    if not np.isfinite(sd) or sd == 0:
        return np.nan
    tol = float(r) if r is not None else 0.2 * sd

    def _embed_centered(mm: int) -> np.ndarray:
        vecs = _embed(x, mm)
        return vecs - vecs.mean(axis=1, keepdims=True)

    def _phi(mm: int) -> float:
        vecs = _embed_centered(mm)
        d = _chebyshev_dist(vecs)
        np.fill_diagonal(d, np.inf)
        sim = np.exp(-((d / tol) ** n_exp))
        sim[~np.isfinite(sim)] = 0.0
        denom = max(vecs.shape[0] - 1, 1)
        return float(sim.sum() / (vecs.shape[0] * denom))

    b = _phi(m)
    a = _phi(m + 1)
    if b <= 0 or a <= 0:
        return np.nan
    return float(-np.log(a / b))


def box_dimension(x: np.ndarray, *, n_scales: int = 8) -> float:
    """
    Box-counting fractal dimension of the ``(t, x(t))`` graph.

    The time axis and amplitude axis are both normalised to ``[0, 1]``
    before counting, then a log(box count) vs. log(1/box size)
    regression slope is fit over dyadic box sizes
    (:math:`1/2, 1/4, \\dots`). Larger for the ragged, high-quality MT
    background than for a smooth or blocky interference shape
    (Fig. 4d).

    Parameters
    ----------
    x : ndarray, shape (n,)
    n_scales : int, default 8
        Number of dyadic box sizes to regress over (capped so the
        finest box size still spans >= 2 samples).

    Returns
    -------
    float
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    if n < 4:
        return np.nan
    t = np.linspace(0.0, 1.0, n)
    xmin, xmax = np.nanmin(x), np.nanmax(x)
    if not np.isfinite(xmin) or not np.isfinite(xmax) or xmax <= xmin:
        return np.nan
    xn = (x - xmin) / (xmax - xmin)

    max_k = int(np.floor(np.log2(max(n // 2, 2))))
    ks = np.arange(1, min(n_scales, max_k) + 1)
    if ks.size < 2:
        return np.nan

    log_inv_eps, log_counts = [], []
    for k in ks:
        n_boxes = 2**k
        eps = 1.0 / n_boxes
        cx = np.minimum((t / eps).astype(int), n_boxes - 1)
        cy = np.minimum((xn / eps).astype(int), n_boxes - 1)
        occupied = len(set(zip(cx.tolist(), cy.tolist())))
        if occupied > 0:
            log_inv_eps.append(np.log(1.0 / eps))
            log_counts.append(np.log(occupied))

    if len(log_inv_eps) < 2:
        return np.nan
    slope, _ = np.polyfit(log_inv_eps, log_counts, 1)
    return float(slope)


def compute_entropy_features(
    segments: np.ndarray,
    *,
    m: int = 2,
    r_frac: float = 0.2,
    fuzzy_n: float = 2.0,
) -> np.ndarray:
    """
    Batch feature extraction: (SE, FE, AE, BD) for each row of
    *segments*.

    Parameters
    ----------
    segments : ndarray, shape (n_segments, win_len)
    m : int, default 2
        Embedding dimension shared by SE/FE/AE.
    r_frac : float, default 0.2
        Tolerance as a fraction of each segment's own std.
    fuzzy_n : float, default 2.0
        Fuzzy-membership steepness (see :func:`fuzzy_entropy`).

    Returns
    -------
    F : ndarray, shape (n_segments, 4)
        Columns ``[sample_entropy, fuzzy_entropy, approx_entropy,
        box_dimension]``. Degenerate/too-short rows are ``NaN``.
    """
    segments = np.atleast_2d(np.asarray(segments, dtype=float))
    n = segments.shape[0]
    F = np.full((n, 4), np.nan, dtype=float)
    for i in range(n):
        seg = segments[i]
        sd = np.nanstd(seg)
        r = r_frac * sd if np.isfinite(sd) and sd > 0 else None
        F[i, 0] = sample_entropy(seg, m=m, r=r)
        F[i, 1] = fuzzy_entropy(seg, m=m, r=r, n_exp=fuzzy_n)
        F[i, 2] = approx_entropy(seg, m=m, r=r)
        F[i, 3] = box_dimension(seg)
    return F


# ─────────────────────────────────────────────────────────────────────────────
# Synthetic training-library generators (Fig. 3 noise morphologies)
# ─────────────────────────────────────────────────────────────────────────────


def _synth_clean(win_len: int, rng: np.random.Generator,
                  amp: float = 1.0) -> np.ndarray:
    """Weak-amplitude, strongly-random background -- colored noise."""
    white = rng.normal(size=win_len)
    spec = np.fft.rfft(white)
    freqs = np.fft.rfftfreq(win_len)
    freqs[0] = freqs[1] if win_len > 1 else 1.0
    spec = spec / (freqs**0.6)
    x = np.fft.irfft(spec, n=win_len)
    x = x / (np.std(x) + 1e-12) * amp
    return x


def _synth_charge_discharge(win_len: int, rng: np.random.Generator,
                             amp: float = 8.0) -> np.ndarray:
    t = np.arange(win_len, dtype=float)
    x = np.zeros(win_len)
    n_pulse = rng.integers(1, 4)
    for _ in range(n_pulse):
        c = rng.integers(0, win_len)
        tau = rng.uniform(2, max(3, win_len / 6))
        sign = rng.choice([-1.0, 1.0])
        env = sign * amp * np.exp(-np.abs(t - c) / tau)
        env[t < c] = 0.0
        x += env
    return x


def _synth_square(win_len: int, rng: np.random.Generator,
                   amp: float = 8.0) -> np.ndarray:
    x = np.zeros(win_len)
    n_step = rng.integers(1, 3)
    edges = sorted(rng.integers(0, win_len, size=2 * n_step))
    level = 0.0
    for i in range(len(edges) - 1):
        level = amp * rng.choice([-1.0, 1.0]) if i % 2 == 0 else 0.0
        x[edges[i] : edges[i + 1]] = level
    return x


def _synth_pulse(win_len: int, rng: np.random.Generator,
                  amp: float = 10.0) -> np.ndarray:
    x = np.zeros(win_len)
    n_spike = rng.integers(1, 6)
    width = max(1, win_len // 40)
    for _ in range(n_spike):
        c = rng.integers(0, win_len)
        lo, hi = max(0, c - width), min(win_len, c + width)
        x[lo:hi] += amp * rng.choice([-1.0, 1.0]) * rng.uniform(0.5, 1.0)
    return x


def _synth_triangle(win_len: int, rng: np.random.Generator,
                     amp: float = 8.0) -> np.ndarray:
    x = np.zeros(win_len)
    c = rng.integers(win_len // 4, 3 * win_len // 4)
    half = rng.integers(2, max(3, win_len // 4))
    sign = rng.choice([-1.0, 1.0])
    lo, hi = max(0, c - half), min(win_len, c + half)
    ramp = amp * sign * (1.0 - np.abs(np.arange(lo, hi) - c) / half)
    x[lo:hi] = ramp
    return x


_NOISE_GENERATORS = (
    _synth_charge_discharge,
    _synth_square,
    _synth_pulse,
    _synth_triangle,
)


def generate_synthetic_library(
    win_len: int,
    *,
    n_per_class: int = 500,
    clean_amp: float = 1.0,
    noise_amp: float = 8.0,
    seed: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Build a synthetic clean/noisy training library, mirroring Fig. 3
    of Gui et al. (2024).

    Each noisy sample is a weak-amplitude colored-noise background
    (the same generator as the clean class) with one strong,
    dominant-amplitude transient superimposed -- charge/discharge,
    square-wave, pulse, or triangle-wave, chosen uniformly at random
    -- matching the "strong interference signal samples have very
    obvious shapes and amplitudes compared to high-quality signal
    samples" description in the paper. Used by
    :meth:`SignalQualityClassifier.from_synthetic` when no measured
    labelled library is available.

    Parameters
    ----------
    win_len : int
        Segment length in samples (e.g. ``round(win_seconds / dt)``).
    n_per_class : int, default 500
        Number of clean and noisy samples each (paper's own default).
    clean_amp : float, default 1.0
        Background amplitude scale (arbitrary units -- only the
        clean/noisy amplitude *ratio* matters to the entropy
        features).
    noise_amp : float, default 8.0
        Transient amplitude scale, relative to *clean_amp*.
    seed : int or None

    Returns
    -------
    segments : ndarray, shape (2 * n_per_class, win_len)
    labels : int ndarray, shape (2 * n_per_class,)
        ``1`` = high-quality/clean, ``0`` = noisy.
    """
    rng = np.random.default_rng(seed)
    segs, labels = [], []
    for _ in range(n_per_class):
        segs.append(_synth_clean(win_len, rng, amp=clean_amp))
        labels.append(1)
    for _ in range(n_per_class):
        bg = _synth_clean(win_len, rng, amp=clean_amp)
        gen = _NOISE_GENERATORS[rng.integers(0, len(_NOISE_GENERATORS))]
        bg = bg + gen(win_len, rng, amp=noise_amp)
        segs.append(bg)
        labels.append(0)
    return np.asarray(segs, dtype=float), np.asarray(labels, dtype=int)


# ─────────────────────────────────────────────────────────────────────────────
# 2b. SignalQualityClassifier (linear-kernel SVM)
# ─────────────────────────────────────────────────────────────────────────────


class SignalQualityClassifier(BaseEMProcessor):
    """
    Linear-kernel SVM classifying 10-s time-series windows as
    high-quality or noisy, from four complexity features.

    Faithful to section 2.2 / 3.1 of Gui et al. (2024): a small,
    high-margin classifier on ``[sample_entropy, fuzzy_entropy,
    approx_entropy, box_dimension]`` (see :func:`compute_entropy_features`)
    is enough to separate strong interference (low complexity, very
    consistent shape) from genuine MT signal (high complexity, weak
    amplitude, strong randomness) -- no deep network is needed for a
    4-feature, high-margin problem like this one.

    Parameters
    ----------
    kernel : str, default "linear"
        Passed to :class:`sklearn.svm.SVC` -- ``"linear"`` matches the
        paper's own choice (Fig. 5).
    C : float, default 1.0
    gamma : str or float, default "scale"
    m : int, default 2
        Embedding dimension shared by SE/FE/AE (see
        :func:`compute_entropy_features`).
    r_frac : float, default 0.2
    fuzzy_n : float, default 2.0
    random_state : int or None

    Examples
    --------
    >>> from pycsamt.ai.processing.tsdenoise import SignalQualityClassifier
    >>> clf = SignalQualityClassifier.from_synthetic(
    ...     win_len=150, random_state=0,
    ... )
    >>> clf.predict(segments)  # doctest: +SKIP
    """

    def __init__(
        self,
        kernel: str = "linear",
        C: float = 1.0,
        gamma: str | float = "scale",
        m: int = 2,
        r_frac: float = 0.2,
        fuzzy_n: float = 2.0,
        random_state: int | None = None,
    ) -> None:
        self.kernel = kernel
        self.C = float(C)
        self.gamma = gamma
        self.m = int(m)
        self.r_frac = float(r_frac)
        self.fuzzy_n = float(fuzzy_n)
        self.random_state = random_state

        self._svc: Any = None
        self._x_mean: np.ndarray | None = None
        self._x_std: np.ndarray | None = None
        self._is_fitted: bool = False

    # ─── factory ───────────────────────────────────────────────────────────

    @classmethod
    def from_synthetic(
        cls,
        win_len: int,
        *,
        n_per_class: int = 500,
        clean_amp: float = 1.0,
        noise_amp: float = 8.0,
        random_state: int | None = None,
        **kwargs,
    ) -> SignalQualityClassifier:
        """
        Build and train from :func:`generate_synthetic_library`.

        Parameters
        ----------
        win_len : int
        n_per_class, clean_amp, noise_amp
            Passed to :func:`generate_synthetic_library`.
        random_state : int or None
        **kwargs
            Passed to the constructor.

        Returns
        -------
        SignalQualityClassifier
        """
        segs, labels = generate_synthetic_library(
            win_len,
            n_per_class=n_per_class,
            clean_amp=clean_amp,
            noise_amp=noise_amp,
            seed=random_state,
        )
        obj = cls(random_state=random_state, **kwargs)
        obj.fit(segs, labels)
        return obj

    # ─── BaseEMProcessor interface ────────────────────────────────────────

    def fit(self, X: np.ndarray, y: np.ndarray, **kwargs
            ) -> SignalQualityClassifier:
        """
        Train on raw segments.

        Parameters
        ----------
        X : ndarray, shape (n_segments, win_len)
            Raw time-series windows (not pre-computed features).
        y : int ndarray, shape (n_segments,)
            ``1`` = high-quality, ``0`` = noisy.

        Returns
        -------
        self
        """
        try:
            from sklearn.svm import SVC
        except ImportError as exc:
            raise ImportError(
                "scikit-learn is required for SignalQualityClassifier"
            ) from exc

        feat = compute_entropy_features(
            X, m=self.m, r_frac=self.r_frac, fuzzy_n=self.fuzzy_n
        )
        y_arr = np.asarray(y, dtype=int)
        valid = np.all(np.isfinite(feat), axis=1)
        if not valid.any():
            raise ValueError(
                "No valid (finite) feature rows -- segments may be too "
                "short or constant."
            )
        feat, y_arr = feat[valid], y_arr[valid]

        self._x_mean = feat.mean(axis=0, keepdims=True)
        self._x_std = feat.std(axis=0, keepdims=True) + 1e-8
        Xn = (feat - self._x_mean) / self._x_std

        self._svc = SVC(
            kernel=self.kernel,
            C=self.C,
            gamma=self.gamma,
            probability=True,
            random_state=self.random_state,
        )
        self._svc.fit(Xn, y_arr)
        self._is_fitted = True
        return self

    def transform(self, X: np.ndarray) -> np.ndarray:
        """
        Quality score for each segment.

        Parameters
        ----------
        X : ndarray, shape (n_segments, win_len)

        Returns
        -------
        scores : ndarray, shape (n_segments,)
            Probability of being high-quality, in ``[0, 1]``.
            ``NaN`` rows (degenerate segments) score ``0.0``
            (treated as noisy/untrustworthy).
        """
        if not self._is_fitted:
            raise RuntimeError("Call fit() before transform().")
        feat = compute_entropy_features(
            X, m=self.m, r_frac=self.r_frac, fuzzy_n=self.fuzzy_n
        )
        scores = np.zeros(feat.shape[0], dtype=float)
        valid = np.all(np.isfinite(feat), axis=1)
        if valid.any():
            Xn = (feat[valid] - self._x_mean) / self._x_std
            pos_col = list(self._svc.classes_).index(1)
            scores[valid] = self._svc.predict_proba(Xn)[:, pos_col]
        return scores

    def predict(self, X: np.ndarray, *, threshold: float = 0.5
                ) -> np.ndarray:
        """Predict labels (1=high-quality, 0=noisy)."""
        return (self.transform(X) >= threshold).astype(int)

    # ─── serialisation ────────────────────────────────────────────────────

    def _get_params(self) -> dict[str, Any]:
        return {
            "kernel": self.kernel,
            "C": self.C,
            "gamma": self.gamma,
            "m": self.m,
            "r_frac": self.r_frac,
            "fuzzy_n": self.fuzzy_n,
            "random_state": self.random_state,
        }

    def _get_weights(self) -> dict[str, np.ndarray]:
        if self._svc is None:
            return {}
        import io
        import pickle

        buf = io.BytesIO()
        pickle.dump(self._svc, buf)
        return {
            "_svc_pickle": np.frombuffer(buf.getvalue(), dtype=np.uint8),
            "_x_mean": self._x_mean,
            "_x_std": self._x_std,
        }

    def _load_weights(self, weights: dict[str, np.ndarray]) -> None:
        self._x_mean = weights.get("_x_mean")
        self._x_std = weights.get("_x_std")
        blob = weights.get("_svc_pickle")
        if blob is not None:
            import io
            import pickle

            self._svc = pickle.load(io.BytesIO(bytes(blob)))
            self._is_fitted = True

    def __repr__(self) -> str:
        status = "fitted" if self._is_fitted else "unfitted"
        return f"SignalQualityClassifier(kernel={self.kernel!r}, {status})"


# ─────────────────────────────────────────────────────────────────────────────
# 3. K-SVD dictionary learning
# ─────────────────────────────────────────────────────────────────────────────


def _normalize_columns(D: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(D, axis=0)
    norms = np.where(norms < 1e-12, 1.0, norms)
    return D / norms


def omp(D: np.ndarray, y: np.ndarray, sparsity: int) -> np.ndarray:
    """
    Orthogonal matching pursuit for one signal.

    Parameters
    ----------
    D : ndarray, shape (L, K)
        Column-normalized dictionary.
    y : ndarray, shape (L,)
    sparsity : int
        Maximum number of active atoms.

    Returns
    -------
    code : ndarray, shape (K,)
    """
    L, K = D.shape
    residual = y.astype(float).copy()
    idx_set: list[int] = []
    code = np.zeros(K, dtype=float)
    coef = np.array([])
    for _ in range(max(1, int(sparsity))):
        proj = D.T @ residual
        proj[idx_set] = 0.0
        j = int(np.argmax(np.abs(proj)))
        if j in idx_set:
            break
        idx_set.append(j)
        Dsub = D[:, idx_set]
        coef, *_ = np.linalg.lstsq(Dsub, y, rcond=None)
        residual = y - Dsub @ coef
        if np.linalg.norm(residual) < 1e-10:
            break
    if idx_set:
        code[idx_set] = coef
    return code


def omp_batch(D: np.ndarray, Y: np.ndarray, sparsity: int) -> np.ndarray:
    """
    Orthogonal matching pursuit for a batch of signals (patches).

    Parameters
    ----------
    D : ndarray, shape (L, K)
    Y : ndarray, shape (L, N)
    sparsity : int

    Returns
    -------
    codes : ndarray, shape (K, N)
    """
    K = D.shape[1]
    N = Y.shape[1]
    X = np.zeros((K, N), dtype=float)
    for i in range(N):
        X[:, i] = omp(D, Y[:, i], sparsity)
    return X


def ksvd_dictionary(
    Y: np.ndarray,
    n_atoms: int,
    sparsity: int,
    *,
    n_iter: int = 10,
    seed: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    K-SVD dictionary learning (Aharon, Elad & Bruckstein, 2006).

    Alternates OMP sparse coding with a per-atom SVD update -- for
    each atom, the rank-1 SVD of the residual restricted to the
    patches currently using that atom simultaneously updates the atom
    and its coefficients, which is what distinguishes K-SVD from a
    plain alternating-least-squares dictionary update (Fig. 2 of
    Gui et al., 2024). Atoms left unused after a coding pass are
    replaced by the training patch with the largest current
    reconstruction error, the standard K-SVD safeguard against
    dictionary collapse when there are few training patches relative
    to *n_atoms*.

    Parameters
    ----------
    Y : ndarray, shape (L, N)
        Training patches, one per column.
    n_atoms : int
        Dictionary size K (capped to ``N`` if larger).
    sparsity : int
        Maximum active atoms per patch (OMP budget).
    n_iter : int, default 10
    seed : int or None

    Returns
    -------
    D : ndarray, shape (L, K)
        Learned, column-normalized dictionary.
    X : ndarray, shape (K, N)
        Final sparse codes.
    """
    rng = np.random.default_rng(seed)
    L, N = Y.shape
    K = min(int(n_atoms), N)

    idx0 = rng.choice(N, size=K, replace=False)
    D = _normalize_columns(Y[:, idx0].copy())

    X = np.zeros((K, N), dtype=float)
    for _ in range(max(1, int(n_iter))):
        X = omp_batch(D, Y, sparsity)
        for k in range(K):
            using = np.nonzero(X[k, :])[0]
            if using.size == 0:
                errs = np.sum((Y - D @ X) ** 2, axis=0)
                worst = int(np.argmax(errs))
                atom = Y[:, worst].copy()
                nrm = np.linalg.norm(atom)
                D[:, k] = atom / nrm if nrm > 1e-12 else rng.normal(size=L)
                continue
            D[:, k] = 0.0
            Ek = Y[:, using] - D @ X[:, using]
            U, S, Vt = np.linalg.svd(Ek, full_matrices=False)
            D[:, k] = U[:, 0]
            X[k, using] = S[0] * Vt[0, :]
        D = _normalize_columns(D)
    return D, X


def _extract_patches(
    x: np.ndarray, patch_len: int, step: int
) -> tuple[np.ndarray, np.ndarray]:
    n = x.size
    patch_len = min(patch_len, n)
    last_start = max(n - patch_len, 0)
    starts = np.arange(0, last_start + 1, max(1, step))
    if starts.size == 0 or starts[-1] != last_start:
        starts = np.append(starts, last_start)
    P = np.stack([x[s : s + patch_len] for s in starts], axis=1)
    return P, starts


def _reconstruct_from_patches(
    patches: np.ndarray, starts: np.ndarray, patch_len: int, n: int
) -> np.ndarray:
    out = np.zeros(n, dtype=float)
    weight = np.zeros(n, dtype=float)
    for i, s in enumerate(starts):
        out[s : s + patch_len] += patches[:, i]
        weight[s : s + patch_len] += 1.0
    weight[weight == 0] = 1.0
    return out / weight


class KSVDDenoiser(BaseEMProcessor):
    """
    K-SVD dictionary-learning denoiser for a single time-series
    excursion (typically the noisy windows :class:`SignalQualityClassifier`
    flagged, concatenated over one contiguous run).

    Self-supervised: :meth:`fit` learns the dictionary directly from
    overlapping patches of the (noisy) input itself -- no pre-training
    or labelled data needed, matching section 2.3 of Gui et al. (2024).
    Because a short, low-sparsity dictionary can represent a stereotyped
    transient far better than the genuinely random background signal,
    the sparse reconstruction of the noisy patches approximates the
    *noise*, not the signal, and :meth:`transform` returns the noisy
    input minus that reconstruction.

    Parameters
    ----------
    patch_len : int, default 32
        Patch length in samples.
    n_atoms : int, default 64
        Dictionary size (capped to the number of extracted patches).
    sparsity : int, default 4
        OMP budget per patch.
    n_iter : int, default 10
        K-SVD alternations.
    overlap : float, default 0.5
        Fractional patch overlap in ``[0, 1)``.
    random_state : int or None

    Examples
    --------
    >>> from pycsamt.ai.processing.tsdenoise import KSVDDenoiser
    >>> den = KSVDDenoiser(random_state=0).fit(noisy_run)
    >>> clean_run = den.transform(noisy_run)
    """

    def __init__(
        self,
        patch_len: int = 32,
        n_atoms: int = 64,
        sparsity: int = 4,
        n_iter: int = 10,
        overlap: float = 0.5,
        random_state: int | None = None,
    ) -> None:
        self.patch_len = int(patch_len)
        self.n_atoms = int(n_atoms)
        self.sparsity = int(sparsity)
        self.n_iter = int(n_iter)
        self.overlap = float(overlap)
        self.random_state = random_state

        self._D: np.ndarray | None = None
        self._patch_len_used: int | None = None
        self._is_fitted: bool = False

    def fit(self, X: np.ndarray, **kwargs) -> KSVDDenoiser:
        """
        Learn the dictionary from patches of *X* itself.

        Parameters
        ----------
        X : ndarray, shape (n,)

        Returns
        -------
        self
        """
        x = np.asarray(X, dtype=float)
        patch_len = max(2, min(self.patch_len, x.size))
        step = max(1, int(round(patch_len * (1.0 - self.overlap))))
        P, _ = _extract_patches(x, patch_len, step)
        self._D, _ = ksvd_dictionary(
            P, self.n_atoms, self.sparsity,
            n_iter=self.n_iter, seed=self.random_state,
        )
        self._patch_len_used = patch_len
        self._is_fitted = True
        return self

    def noise_profile(self, X: np.ndarray) -> np.ndarray:
        """
        Sparse reconstruction of *X* against the learned dictionary --
        the extracted noise contour (Fig. 13c of Gui et al., 2024).
        """
        if not self._is_fitted:
            raise RuntimeError("Call fit() before noise_profile().")
        x = np.asarray(X, dtype=float)
        patch_len = self._patch_len_used
        step = max(1, int(round(patch_len * (1.0 - self.overlap))))
        P, starts = _extract_patches(x, patch_len, step)
        codes = omp_batch(self._D, P, self.sparsity)
        rec = self._D @ codes
        return _reconstruct_from_patches(rec, starts, patch_len, x.size)

    def transform(self, X: np.ndarray) -> np.ndarray:
        """
        Denoise *X*: ``X - noise_profile(X)``.

        Parameters
        ----------
        X : ndarray, shape (n,)

        Returns
        -------
        ndarray, shape (n,)
        """
        return np.asarray(X, dtype=float) - self.noise_profile(X)

    def _get_params(self) -> dict[str, Any]:
        return {
            "patch_len": self.patch_len,
            "n_atoms": self.n_atoms,
            "sparsity": self.sparsity,
            "n_iter": self.n_iter,
            "overlap": self.overlap,
            "random_state": self.random_state,
        }

    def _get_weights(self) -> dict[str, np.ndarray]:
        if self._D is None:
            return {}
        return {
            "_D": self._D,
            "_patch_len_used": np.array([self._patch_len_used]),
        }

    def _load_weights(self, weights: dict[str, np.ndarray]) -> None:
        if "_D" in weights:
            self._D = weights["_D"]
            self._patch_len_used = int(weights["_patch_len_used"][0])
            self._is_fitted = True

    def __repr__(self) -> str:
        status = "fitted" if self._is_fitted else "unfitted"
        return f"KSVDDenoiser(n_atoms={self.n_atoms}, {status})"


# ─────────────────────────────────────────────────────────────────────────────
# 4. Orchestration: TimeSeriesDenoiser
# ─────────────────────────────────────────────────────────────────────────────


def _contiguous_runs(mask: np.ndarray) -> list[tuple[int, int]]:
    """Start/stop (exclusive) index pairs of contiguous ``True`` runs."""
    runs = []
    n = mask.size
    i = 0
    while i < n:
        if mask[i]:
            j = i
            while j < n and mask[j]:
                j += 1
            runs.append((i, j))
            i = j
        else:
            i += 1
    return runs


class TimeSeriesDenoiser(BaseEMProcessor):
    """
    MMF-SVM-K-SVD field time-series denoiser (Gui et al., 2024).

    Orchestrates :func:`mmf_split`, :class:`SignalQualityClassifier`,
    and :class:`KSVDDenoiser` into the three-stage pipeline described
    in the module docstring, and works directly on a
    :class:`~pycsamt.ts.TSData` record via :meth:`apply` -- the same
    place :func:`pycsamt.ts.process.preprocess` sits, just upstream of
    it.

    Parameters
    ----------
    mmf_size : int or None
        Structuring-element length (samples) for :func:`mmf_split`.
        ``None`` derives it from each channel's own ``dt`` (~2 s) --
        the single most important knob to tune: it must exceed the
        widest interference pulse in the record, see
        :func:`mmf_split`'s docstring for a worked before/after.
    win_seconds : float, default 10.0
        SVM classification window length -- the paper's own choice.
    quality_clf : SignalQualityClassifier or None
        Pre-trained classifier to reuse (e.g. across stations sharing
        the same sampling rate *and* amplitude scale). When ``None``,
        one is trained automatically the first time :meth:`transform`
        needs it, with :func:`generate_synthetic_library`'s
        ``clean_amp`` set to a robust (MAD-based) estimate of that
        channel's own high-frequency background level -- required
        because raw channel amplitude varies enormously across
        acquisition systems (raw ADC counts, mV/km, nT, ...); training
        the classifier at a fixed, unrelated amplitude scale leaves
        every real feature vector far outside the training envelope, a
        regime where a linear SVM's decision boundary no longer
        reflects "more extreme = more noisy" and can misclassify
        confidently.
    noise_amp_ratio : float, default 8.0
        Synthetic interference amplitude as a multiple of the
        auto-detected (or supplied) ``clean_amp`` -- matches the
        roughly 8x contrast Fig. 3 of Gui et al. (2024) shows between
        their strong-interference and high-quality signal samples.
    patch_len, n_atoms, sparsity, ksvd_iter, ksvd_overlap
        Passed to each per-run :class:`KSVDDenoiser`.
    quality_threshold : float, default 0.5
        Windows scoring below this are treated as noisy.
    random_state : int or None

    Examples
    --------
    >>> from pycsamt.ts import read_ts
    >>> from pycsamt.ai.processing.tsdenoise import TimeSeriesDenoiser
    >>> ts = read_ts("data/MT/TS/kap103as.ts/kap103as.ts")  # doctest: +SKIP
    >>> clean_ts = TimeSeriesDenoiser(random_state=0).apply(ts)  # doctest: +SKIP
    """

    def __init__(
        self,
        mmf_size: int | None = None,
        win_seconds: float = 10.0,
        quality_clf: SignalQualityClassifier | None = None,
        noise_amp_ratio: float = 8.0,
        patch_len: int = 32,
        n_atoms: int = 64,
        sparsity: int = 4,
        ksvd_iter: int = 10,
        ksvd_overlap: float = 0.5,
        quality_threshold: float = 0.5,
        random_state: int | None = None,
    ) -> None:
        self.mmf_size = mmf_size
        self.win_seconds = float(win_seconds)
        self.noise_amp_ratio = float(noise_amp_ratio)
        self.patch_len = int(patch_len)
        self.n_atoms = int(n_atoms)
        self.sparsity = int(sparsity)
        self.ksvd_iter = int(ksvd_iter)
        self.ksvd_overlap = float(ksvd_overlap)
        self.quality_threshold = float(quality_threshold)
        self.random_state = random_state

        self._clf = quality_clf
        self._diagnostics: list[pd.DataFrame] = []
        self._is_fitted: bool = quality_clf is not None

    # ─── BaseEMProcessor interface ────────────────────────────────────────

    def fit(
        self,
        X: Any = None,
        y: Any = None,
        *,
        dt: float | None = None,
        clean_amp: float = 1.0,
        n_per_class: int = 500,
        **kwargs,
    ) -> TimeSeriesDenoiser:
        """
        Train (or attach) the internal :class:`SignalQualityClassifier`.

        Optional -- :meth:`transform`/:meth:`apply` train one
        automatically on first use if none is attached, with
        *clean_amp* set from the channel's own data (see
        :attr:`quality_clf`'s docstring entry above). Call explicitly
        to control *dt* / *clean_amp* / *n_per_class* directly, or to
        reuse the same classifier for several :meth:`apply` calls
        without retraining each time.

        Parameters
        ----------
        X, y : ignored
            Present for :class:`~pycsamt.ai._base.BaseEMProcessor`
            interface compatibility; training data always comes from
            :func:`generate_synthetic_library`.
        dt : float
            Sampling interval (s), needed to size the classification
            window (``round(win_seconds / dt)`` samples).
        clean_amp : float, default 1.0
            Background amplitude scale for the synthetic training
            library -- must be in the same units, and roughly the same
            order of magnitude, as the channel's actual high-frequency
            background level (see :attr:`quality_clf`'s docstring
            entry). Left at the default ``1.0`` only makes sense for
            data already normalised to that scale.
        n_per_class : int, default 500
            Passed to :func:`generate_synthetic_library`.

        Returns
        -------
        self
        """
        if dt is None:
            raise ValueError("TimeSeriesDenoiser.fit() requires `dt`.")
        win_len = max(8, int(round(self.win_seconds / float(dt))))
        self._clf = SignalQualityClassifier.from_synthetic(
            win_len,
            n_per_class=n_per_class,
            clean_amp=clean_amp,
            noise_amp=clean_amp * self.noise_amp_ratio,
            random_state=self.random_state,
        )
        self._is_fitted = True
        return self

    def transform(
        self, X: np.ndarray, *, dt: float, channel: str = ""
    ) -> np.ndarray:
        """
        Denoise a single 1-D channel array.

        Parameters
        ----------
        X : ndarray, shape (n,)
        dt : float
            Sampling interval (s).
        channel : str
            Label recorded in :attr:`diagnostics_` (purely
            informational).

        Returns
        -------
        ndarray, shape (n,)
        """
        x = np.asarray(X, dtype=float)
        nanmask = np.isnan(x)
        xf = _fill_nan_linear(x)

        low, high = mmf_split(xf, size=self.mmf_size, dt=dt)

        if not self._is_fitted or self._clf is None:
            # Robust (MAD-based) background level -- deliberately not a
            # plain std, which the very interference being screened for
            # would inflate and defeat the point of the estimate.
            mad = np.median(np.abs(high - np.median(high)))
            clean_amp = max(1.4826 * mad, 1e-12)
            self.fit(dt=dt, clean_amp=clean_amp)

        win_len = max(8, int(round(self.win_seconds / float(dt))))
        n = high.size
        n_win = max(1, n // win_len)
        edges = np.linspace(0, n_win * win_len, n_win + 1, dtype=int)
        if edges[-1] < n:
            edges = np.append(edges, n)
        edges = np.minimum(edges, n)

        windows = [high[edges[i] : edges[i + 1]] for i in range(len(edges) - 1)]
        win_lens = np.array([w.size for w in windows])
        long_enough = win_lens >= max(8, self.patch_len // 2)

        scores = np.ones(len(windows), dtype=float)
        if long_enough.any():
            padded = np.stack(
                [
                    np.pad(w, (0, win_len - w.size), mode="edge")
                    if w.size < win_len else w[:win_len]
                    for w, ok in zip(windows, long_enough) if ok
                ],
                axis=0,
            )
            scores[long_enough] = self._clf.transform(padded)

        noisy = long_enough & (scores < self.quality_threshold)

        out_high = high.copy()
        for start, stop in _contiguous_runs(noisy):
            s0, s1 = int(edges[start]), int(edges[stop])
            run = high[s0:s1]
            if run.size < 4:
                continue
            den = KSVDDenoiser(
                patch_len=self.patch_len,
                n_atoms=self.n_atoms,
                sparsity=self.sparsity,
                n_iter=self.ksvd_iter,
                overlap=self.ksvd_overlap,
                random_state=self.random_state,
            ).fit(run)
            out_high[s0:s1] = den.transform(run)

        denoised = low + out_high
        denoised = np.where(nanmask, np.nan, denoised)

        self._diagnostics.append(
            pd.DataFrame(
                {
                    "channel": channel,
                    "win_start": edges[:-1] * float(dt),
                    "win_stop": edges[1:] * float(dt),
                    "score": scores,
                    "label": np.where(noisy, "noisy", "clean"),
                }
            )
        )
        return denoised

    # ─── TSData-in / TSData-out ────────────────────────────────────────────

    def apply(
        self,
        ts: Any,
        *,
        channels: list[str] | None = None,
        inplace: bool = False,
    ) -> Any:
        """
        Denoise a :class:`~pycsamt.ts.TSData` record.

        Parameters
        ----------
        ts : TSData
        channels : list of str or None
            Channels to denoise. ``None`` (default) processes every
            channel present.
        inplace : bool, default ``False``
            When ``False`` (the default), *ts* is left untouched and
            a denoised copy is returned. When ``True``, *ts* is
            modified in place and returned.

        Returns
        -------
        TSData
            Same station/metadata, with each requested channel
            replaced by its denoised reconstruction.

        Examples
        --------
        >>> clean_ts = TimeSeriesDenoiser(
        ...     random_state=0,
        ... ).apply(ts)  # doctest: +SKIP
        """
        if ts.dt is None:
            raise ValueError(
                "TimeSeriesDenoiser.apply() requires ts.dt to be set."
            )
        chans = list(channels) if channels is not None else ts.channels()

        out = ts if inplace else ts.copy_meta()
        if not inplace:
            for cid in ts.channels():
                out.add_channel(cid, ts.get(cid))

        self._diagnostics = []
        for cid in chans:
            den = self.transform(ts.get(cid), dt=ts.dt, channel=cid)
            out.add_channel(cid, den)
        return out

    @property
    def diagnostics_(self) -> pd.DataFrame:
        """
        Per-window classification diagnostics from the last
        :meth:`apply`/:meth:`transform` call.

        Returns
        -------
        DataFrame
            Columns: ``channel``, ``win_start``, ``win_stop`` (s),
            ``score`` (probability of high-quality), ``label``
            (``"clean"``/``"noisy"``).
        """
        if not self._diagnostics:
            return pd.DataFrame(
                columns=["channel", "win_start", "win_stop", "score",
                         "label"]
            )
        return pd.concat(self._diagnostics, ignore_index=True)

    def _get_params(self) -> dict[str, Any]:
        return {
            "mmf_size": self.mmf_size,
            "win_seconds": self.win_seconds,
            "noise_amp_ratio": self.noise_amp_ratio,
            "patch_len": self.patch_len,
            "n_atoms": self.n_atoms,
            "sparsity": self.sparsity,
            "ksvd_iter": self.ksvd_iter,
            "ksvd_overlap": self.ksvd_overlap,
            "quality_threshold": self.quality_threshold,
            "random_state": self.random_state,
        }

    def _get_weights(self) -> dict[str, np.ndarray]:
        if self._clf is None:
            return {}
        return {f"clf__{k}": v for k, v in self._clf._get_weights().items()}

    def _load_weights(self, weights: dict[str, np.ndarray]) -> None:
        clf_weights = {
            k[len("clf__") :]: v
            for k, v in weights.items()
            if k.startswith("clf__")
        }
        if clf_weights:
            self._clf = SignalQualityClassifier(random_state=self.random_state)
            self._clf._load_weights(clf_weights)
            self._is_fitted = True

    def __repr__(self) -> str:
        status = "fitted" if self._is_fitted else "unfitted"
        return f"TimeSeriesDenoiser(win_seconds={self.win_seconds}, {status})"


# ─────────────────────────────────────────────────────────────────────────────
# 5. Validation metrics (paper's Eqs. 7-8 / Table 1)
# ─────────────────────────────────────────────────────────────────────────────


def snr_db(clean: np.ndarray, denoised: np.ndarray) -> float:
    r"""
    Signal-to-noise ratio, Eq. (7) of Gui et al. (2024).

    .. math::

        SNR = 20 \log_{10} \frac{\|m(s)\|_2}{\|m(s) - g(s)\|_2}

    Parameters
    ----------
    clean : ndarray
        Reference (original high-quality) signal :math:`m(s)`.
    denoised : ndarray
        Reconstructed signal :math:`g(s)`.

    Returns
    -------
    float
        Decibels; higher is better.
    """
    clean = np.asarray(clean, dtype=float)
    denoised = np.asarray(denoised, dtype=float)
    num = np.linalg.norm(clean)
    den = np.linalg.norm(clean - denoised)
    if den < 1e-24:
        return float("inf")
    return float(20.0 * np.log10(num / den))


def time_domain_ncc(clean: np.ndarray, denoised: np.ndarray) -> float:
    r"""
    Normalized cross-correlation, Eq. (8) of Gui et al. (2024).

    .. math::

        NCC = \frac{\sum_s m(s) g(s)}
                    {\sqrt{\sum_s m(s)^2 \sum_s g(s)^2}}

    Parameters
    ----------
    clean, denoised : ndarray

    Returns
    -------
    float
        In ``[-1, 1]``; closer to 1 means higher similarity.
    """
    clean = np.asarray(clean, dtype=float)
    denoised = np.asarray(denoised, dtype=float)
    denom = np.sqrt(np.sum(clean**2) * np.sum(denoised**2))
    if denom < 1e-24:
        return float("nan")
    return float(np.sum(clean * denoised) / denom)
