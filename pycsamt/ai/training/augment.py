# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Data augmentation for EM training datasets.

Augmenters are callable objects that accept ``(X, y)`` numpy pairs
and return modified copies.  They are designed to be composed with
:class:`Compose` and used inside a training loop or as a preprocessing
step on :class:`~pycsamt.ai.training.dataset.EMDataset`.

Available augmenters
--------------------
:class:`AugmentNoise`
    Additive log-Gaussian noise on feature amplitudes.
:class:`AugmentStaticShift`
    Random static shift — multiply all log-amplitude features by a
    site-specific constant drawn from a log-uniform distribution.
:class:`AugmentFreqDrop`
    Randomly zero-out a fraction of frequency channels.
:class:`AugmentMixup`
    Convex combination of two training samples (Zhang et al. 2018).
:class:`Compose`
    Sequential chain of augmenters.
:class:`RandomApply`
    Apply an augmenter with probability *p*.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Sequence

import numpy as np

__all__ = [
    "AugmentNoise",
    "AugmentStaticShift",
    "AugmentFreqDrop",
    "AugmentMixup",
    "Compose",
    "RandomApply",
]


# ─────────────────────────────────────────────────────────────────────────────
# Base class
# ─────────────────────────────────────────────────────────────────────────────


class _BaseAugmenter(ABC):
    """Abstract base for all EM data augmenters."""

    def __call__(
        self,
        X: np.ndarray,
        y: np.ndarray | None = None,
        *,
        rng: np.random.Generator | None = None,
    ) -> tuple[np.ndarray, np.ndarray | None]:
        """
        Apply the augmentation.

        Parameters
        ----------
        X : ndarray (n_samples, n_features)
        y : ndarray (n_samples, n_params) or None
        rng : numpy Generator or None
            Reproducibility.  When ``None`` a default Generator is used.

        Returns
        -------
        X_aug : ndarray
        y_aug : ndarray or None
        """
        if rng is None:
            rng = np.random.default_rng()
        return self._apply(X, y, rng)

    @abstractmethod
    def _apply(
        self,
        X: np.ndarray,
        y: np.ndarray | None,
        rng: np.random.Generator,
    ) -> tuple[np.ndarray, np.ndarray | None]:
        """Subclass implementation."""


# ─────────────────────────────────────────────────────────────────────────────
# AugmentNoise
# ─────────────────────────────────────────────────────────────────────────────


class AugmentNoise(_BaseAugmenter):
    """
    Add multiplicative log-Gaussian noise to feature amplitudes.

    For each sample a noise vector

    .. math::

        \\boldsymbol{\\varepsilon} \\sim \\mathcal{N}(\\mathbf{0},\\, \\sigma^2 \\mathbf{I})

    is added to the feature vector.  If ``phase_sigma`` is given a
    separate (typically smaller) noise level is used on phase features
    (assumed to occupy the second half of the feature vector when
    ``include_phase=True``).

    Parameters
    ----------
    sigma : float, default 0.05
        Noise standard deviation applied to all features.
    phase_sigma : float or None
        Separate noise level for phase features.  ``None`` → use *sigma*.
    clip : float or None
        Clip noise magnitude to ``[-clip, +clip]`` standard deviations.
    """

    def __init__(
        self,
        sigma: float = 0.05,
        phase_sigma: float | None = None,
        clip: float | None = 3.0,
    ) -> None:
        self.sigma = float(sigma)
        self.phase_sigma = (
            float(phase_sigma) if phase_sigma is not None else sigma
        )
        self.clip = clip

    def _apply(self, X, y, rng):
        X_out = X.copy()
        noise = rng.normal(0.0, self.sigma, X_out.shape).astype(X_out.dtype)
        if self.clip is not None:
            noise = np.clip(
                noise, -self.clip * self.sigma, self.clip * self.sigma
            )
        X_out += noise
        return X_out, y

    def __repr__(self) -> str:
        return f"AugmentNoise(sigma={self.sigma})"


# ─────────────────────────────────────────────────────────────────────────────
# AugmentStaticShift
# ─────────────────────────────────────────────────────────────────────────────


class AugmentStaticShift(_BaseAugmenter):
    """
    Apply a random static shift to amplitude features.

    Static shift is a common MT artefact: all apparent-resistivity
    values at one site are scaled by a constant factor due to local
    near-surface heterogeneity.  This augmenter simulates the effect
    by multiplying the first ``n_amp_features`` feature columns by a
    site-specific log-uniform random factor.

    Parameters
    ----------
    shift_range : (lo, hi), default (0.3, 3.0)
        Range of the multiplicative factor in linear scale.
    n_amp_features : int or None
        Number of amplitude columns at the start of the feature vector.
        ``None`` → shift all features (use when features are log₁₀ amplitude
        only, without phase).
    per_sample : bool
        If ``True`` draw an independent shift for each sample; otherwise
        draw one shift per batch.
    """

    def __init__(
        self,
        shift_range: float | tuple[float, float] = (0.3, 3.0),
        n_amp_features: int | None = None,
        per_sample: bool = True,
    ) -> None:
        if np.isscalar(shift_range):
            # Interpret scalar s as ±s log-decades: linear range (10^-s, 10^+s).
            s = abs(float(shift_range))
            shift_range = (10.0**-s, 10.0**s)
        self.shift_range = tuple(float(v) for v in shift_range)
        self.n_amp_features = n_amp_features
        self.per_sample = bool(per_sample)

    def _apply(self, X, y, rng):
        lo, hi = self.shift_range
        X_out = X.copy()
        n_feat = X_out.shape[1] if X_out.ndim > 1 else len(X_out)
        n_amp = (
            self.n_amp_features if self.n_amp_features is not None else n_feat
        )

        if X_out.ndim == 1:
            shift = rng.uniform(np.log10(lo), np.log10(hi))
            X_out[:n_amp] += shift
        else:
            n = len(X_out)
            if self.per_sample:
                shifts = rng.uniform(np.log10(lo), np.log10(hi), size=(n, 1))
            else:
                shifts = np.full(
                    (n, 1), rng.uniform(np.log10(lo), np.log10(hi))
                )
            X_out[:, :n_amp] += shifts.astype(X_out.dtype)

        return X_out, y

    def __repr__(self) -> str:
        return f"AugmentStaticShift(shift_range={self.shift_range})"


# ─────────────────────────────────────────────────────────────────────────────
# AugmentFreqDrop
# ─────────────────────────────────────────────────────────────────────────────


class AugmentFreqDrop(_BaseAugmenter):
    """
    Randomly zero-out a fraction of frequency channels.

    Simulates missing data at individual frequencies (bad data periods,
    cultural noise bursts, dead-band contamination).

    Parameters
    ----------
    drop_rate : float, default 0.1
        Fraction of feature columns to zero-out per sample.
    contiguous : bool, default False
        When ``True`` drop a contiguous block of channels (simulates
        a dead band in frequency).
    fill_value : float, default 0.0
        Replacement value for dropped channels.
    """

    def __init__(
        self,
        drop_rate: float = 0.1,
        contiguous: bool = False,
        fill_value: float = 0.0,
    ) -> None:
        self.drop_rate = float(drop_rate)
        self.contiguous = bool(contiguous)
        self.fill_value = float(fill_value)

    def _apply(self, X, y, rng):
        X_out = X.copy()
        n_feat = X_out.shape[-1]
        n_drop = max(1, int(n_feat * self.drop_rate))

        if X_out.ndim == 1:
            idx = self._drop_idx(n_feat, n_drop, rng)
            X_out[idx] = self.fill_value
        else:
            for i in range(len(X_out)):
                idx = self._drop_idx(n_feat, n_drop, rng)
                X_out[i, idx] = self.fill_value

        return X_out, y

    def _drop_idx(self, n_feat, n_drop, rng):
        if self.contiguous:
            start = rng.integers(0, max(1, n_feat - n_drop))
            return np.arange(start, start + n_drop)
        return rng.choice(n_feat, size=n_drop, replace=False)

    def __repr__(self) -> str:
        return f"AugmentFreqDrop(drop_rate={self.drop_rate})"


# ─────────────────────────────────────────────────────────────────────────────
# AugmentMixup
# ─────────────────────────────────────────────────────────────────────────────


class AugmentMixup(_BaseAugmenter):
    """
    Mixup augmentation (Zhang et al. 2018).

    Randomly pairs samples within the batch and creates convex
    combinations:

    .. math::

        \\tilde{\\mathbf{x}} = \\lambda \\mathbf{x}_i + (1-\\lambda) \\mathbf{x}_j,
        \\quad
        \\tilde{\\mathbf{y}} = \\lambda \\mathbf{y}_i + (1-\\lambda) \\mathbf{y}_j

    where :math:`\\lambda \\sim \\mathrm{Beta}(\\alpha, \\alpha)`.

    Parameters
    ----------
    alpha : float, default 0.2
        Beta distribution shape parameter.  Small values (≈0.1) produce
        near-original samples; ``alpha=1.0`` is uniform mixing.
    """

    def __init__(self, alpha: float = 0.2) -> None:
        self.alpha = float(alpha)

    def _apply(self, X, y, rng):
        if X.ndim == 1 or len(X) < 2:
            return X, y

        n = len(X)
        lam = rng.beta(self.alpha, self.alpha, size=n).astype(X.dtype)
        perm = rng.permutation(n)

        lam_x = lam[:, np.newaxis]
        X_out = lam_x * X + (1 - lam_x) * X[perm]

        if y is not None:
            lam_y = lam[:, np.newaxis] if y.ndim > 1 else lam
            y_out = lam_y * y + (1 - lam_y) * y[perm]
        else:
            y_out = None

        return X_out, y_out

    def __repr__(self) -> str:
        return f"AugmentMixup(alpha={self.alpha})"


# ─────────────────────────────────────────────────────────────────────────────
# Compose
# ─────────────────────────────────────────────────────────────────────────────


class Compose:
    """
    Sequentially apply a list of augmenters.

    Parameters
    ----------
    augmenters : list of _BaseAugmenter
        Applied in order.
    seed : int or None
        Seed for the shared random generator.

    Examples
    --------
    >>> aug = Compose([
    ...     AugmentStaticShift(shift_range=(0.5, 2.0)),
    ...     AugmentNoise(sigma=0.03),
    ...     AugmentFreqDrop(drop_rate=0.05),
    ... ])
    >>> X_aug, y_aug = aug(X_train, y_train)
    """

    def __init__(
        self,
        augmenters: Sequence[_BaseAugmenter],
        seed: int | None = None,
    ) -> None:
        self.augmenters = list(augmenters)
        self._rng = np.random.default_rng(seed)

    def __call__(
        self,
        X: np.ndarray,
        y: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray | None]:
        for aug in self.augmenters:
            # Only _BaseAugmenter subclasses accept the rng kwarg;
            # wrappers (RandomApply, Compose) manage their own rngs.
            if isinstance(aug, _BaseAugmenter):
                X, y = aug(X, y, rng=self._rng)
            else:
                X, y = aug(X, y)
        return X, y

    def __repr__(self) -> str:
        inner = ", ".join(repr(a) for a in self.augmenters)
        return f"Compose([{inner}])"


# ─────────────────────────────────────────────────────────────────────────────
# RandomApply
# ─────────────────────────────────────────────────────────────────────────────


class RandomApply:
    """
    Apply an augmenter with probability *p*.

    Parameters
    ----------
    augmenter : _BaseAugmenter
    p : float, default 0.5
    seed : int or None

    Examples
    --------
    >>> aug = RandomApply(AugmentMixup(alpha=0.3), p=0.4)
    >>> X_aug, y_aug = aug(X_train, y_train)
    """

    def __init__(
        self,
        augmenter: _BaseAugmenter,
        p: float = 0.5,
        seed: int | None = None,
    ) -> None:
        self.augmenter = augmenter
        self.p = float(p)
        self._rng = np.random.default_rng(seed)

    def __call__(
        self,
        X: np.ndarray,
        y: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray | None]:
        if self._rng.random() < self.p:
            return self.augmenter(X, y, rng=self._rng)
        return X, y

    def __repr__(self) -> str:
        return f"RandomApply({self.augmenter!r}, p={self.p})"
