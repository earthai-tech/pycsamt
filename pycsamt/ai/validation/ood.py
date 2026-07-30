# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Out-of-distribution checks against the training realization set.

These checks support the "OOD sensitivity" column of the
Uncertainty row of the validation matrix, and the AI-inversion
plan's requirement that predictions outside training support are
rejected or flagged rather than returned as confident maps.  Inputs
are feature vectors, e.g. survey/realization descriptors, shaped
``(n_samples, n_features)``; a common source is per-realization
summary statistics of the geological priors in
:mod:`pycsamt.ai.geology`.

All functions operate on plain NumPy arrays so the module stays
importable without an optional deep-learning backend.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.spatial.distance import cdist

__all__ = [
    "OODReport",
    "ood_score",
    "flag_out_of_distribution",
]

_METHODS = ("mahalanobis", "knn")


def _readonly(value: Any, dtype: Any | None = None) -> np.ndarray:
    """Return a read-only copy of *value* as an ndarray."""
    array = np.array(value, dtype=dtype, copy=True)
    array.setflags(write=False)
    return array


def _validate_2d(name: str, value: Any) -> np.ndarray:
    """Return *value* as a finite float ``(n_samples, n_features)``
    array.
    """
    array = np.asarray(value, dtype=float)
    if array.ndim != 2:
        raise ValueError(
            f"{name} must be a 2-D (n_samples, n_features) array."
        )
    if array.shape[0] < 1 or array.shape[1] < 1:
        raise ValueError(f"{name} must not be empty.")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite.")
    return array


def _mahalanobis_scores(x: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """Return the Mahalanobis distance of each row of x from the
    reference set's mean and covariance.
    """
    n_samples, n_features = reference.shape
    if n_samples <= n_features:
        raise ValueError(
            "reference must contain more samples than features "
            f"({n_samples} <= {n_features}) to estimate a full-rank "
            "covariance for the mahalanobis method."
        )
    mean = reference.mean(axis=0)
    covariance = np.atleast_2d(np.cov(reference, rowvar=False))
    try:
        inverse_covariance = np.linalg.inv(covariance)
    except np.linalg.LinAlgError as error:
        raise ValueError(
            "reference covariance is singular; the mahalanobis "
            "method requires linearly independent features."
        ) from error
    difference = x - mean
    squared = np.einsum(
        "ij,jk,ik->i", difference, inverse_covariance, difference
    )
    return np.sqrt(np.clip(squared, 0.0, None))


def _knn_scores(x: np.ndarray, reference: np.ndarray, k: int) -> np.ndarray:
    """Return each row of x's distance to its k-th nearest neighbour
    in reference.
    """
    if (
        not isinstance(k, int)
        or isinstance(k, bool)
        or k < 1
        or k > reference.shape[0]
    ):
        raise ValueError(f"k must be an integer in [1, {reference.shape[0]}].")
    distances = cdist(x, reference)
    return np.sort(distances, axis=1)[:, k - 1]


def _self_knn_scores(reference: np.ndarray, k: int) -> np.ndarray:
    """Return each reference row's distance to its k-th nearest
    neighbour among the other reference rows, excluding itself.
    """
    upper = reference.shape[0] - 1
    if not isinstance(k, int) or isinstance(k, bool) or k < 1 or k > upper:
        raise ValueError(
            f"k must be an integer in [1, {upper}] for "
            "self-referential k-NN scoring."
        )
    distances = cdist(reference, reference)
    np.fill_diagonal(distances, np.inf)
    return np.sort(distances, axis=1)[:, k - 1]


def ood_score(
    x: Any,
    reference: Any,
    *,
    method: str = "mahalanobis",
    k: int = 5,
) -> np.ndarray:
    """Score how far new inputs fall from the training distribution.

    Parameters
    ----------
    x : array-like, shape (n_samples, n_features)
        Inputs to score, e.g. new survey feature vectors.
    reference : array-like, shape (n_reference, n_features)
        Training-set feature vectors defining the support region.
    method : {"mahalanobis", "knn"}, default="mahalanobis"
        ``"mahalanobis"`` measures deviation from the reference
        mean/covariance and requires ``n_reference > n_features``.
        ``"knn"`` measures Euclidean distance to the ``k``-th
        nearest reference point and makes no distributional
        assumption. If ``x`` shares exact points with ``reference``,
        their k-NN distance to those points is zero; use
        :func:`flag_out_of_distribution` for a leave-one-out
        self-score instead of passing ``reference`` as ``x`` here.
    k : int, default=5
        Neighbour rank used by ``method="knn"``. Ignored otherwise.

    Returns
    -------
    ndarray, shape (n_samples,)
        Higher values indicate inputs farther from the training
        support.

    Examples
    --------
    >>> import numpy as np
    >>> reference = np.array(
    ...     [
    ...         [0.0, 0.0],
    ...         [1.0, 0.0],
    ...         [0.0, 1.0],
    ...         [-1.0, 0.0],
    ...         [0.0, -1.0],
    ...         [0.5, 0.5],
    ...     ]
    ... )
    >>> x = np.array([[0.0, 0.0], [50.0, 50.0]])
    >>> scores = ood_score(x, reference, method="knn", k=2)
    >>> scores[0] < scores[1]
    True
    """
    if method not in _METHODS:
        raise ValueError(f"method must be one of {_METHODS}.")
    x_arr = _validate_2d("x", x)
    ref_arr = _validate_2d("reference", reference)
    if x_arr.shape[1] != ref_arr.shape[1]:
        raise ValueError(
            "x and reference must share the same number of features."
        )
    if method == "mahalanobis":
        return _mahalanobis_scores(x_arr, ref_arr)
    return _knn_scores(x_arr, ref_arr, k)


@dataclass(frozen=True)
class OODReport:
    """Immutable out-of-distribution screening result.

    Parameters
    ----------
    scores : ndarray
        Per-sample distance from :func:`ood_score`.
    threshold : float
        Score above which a sample is flagged.
    flagged : ndarray of bool
        Whether each sample exceeds ``threshold``, same shape as
        ``scores``.
    method : {"mahalanobis", "knn"}
        Distance measure used to compute ``scores``.
    quantile : float or None
        Quantile of the reference self-scores used to derive
        ``threshold``, or ``None`` when an explicit ``threshold``
        was supplied instead.
    n_reference, n_features : int
        Size of the reference set used to define support.

    Examples
    --------
    >>> import numpy as np
    >>> reference = np.array(
    ...     [
    ...         [0.0, 0.0],
    ...         [1.0, 0.0],
    ...         [0.0, 1.0],
    ...         [-1.0, 0.0],
    ...         [0.0, -1.0],
    ...         [0.5, 0.5],
    ...     ]
    ... )
    >>> x = np.array([[0.0, 0.0], [50.0, 50.0]])
    >>> report = flag_out_of_distribution(
    ...     x, reference, method="knn", k=2, quantile=0.5
    ... )
    >>> report.flagged.tolist()
    [False, True]
    """

    scores: np.ndarray
    threshold: float
    flagged: np.ndarray
    method: str
    quantile: float | None
    n_reference: int
    n_features: int

    def __post_init__(self) -> None:
        if self.method not in _METHODS:
            raise ValueError(f"method must be one of {_METHODS}.")
        scores = np.asarray(self.scores, dtype=float)
        if scores.ndim != 1 or scores.size == 0:
            raise ValueError("scores must be a non-empty 1-D array.")
        flagged = np.asarray(self.flagged, dtype=bool)
        if flagged.shape != scores.shape:
            raise ValueError("flagged must have the same shape as scores.")
        threshold = float(self.threshold)
        if not np.isfinite(threshold):
            raise ValueError("threshold must be finite.")
        quantile = None if self.quantile is None else float(self.quantile)
        if quantile is not None and not (0.0 < quantile < 1.0):
            raise ValueError("quantile must be within (0, 1).")
        if self.n_reference < 1 or self.n_features < 1:
            raise ValueError("n_reference and n_features must be positive.")
        object.__setattr__(self, "scores", _readonly(scores))
        object.__setattr__(self, "flagged", _readonly(flagged, bool))
        object.__setattr__(self, "threshold", threshold)
        object.__setattr__(self, "quantile", quantile)
        object.__setattr__(self, "n_reference", int(self.n_reference))
        object.__setattr__(self, "n_features", int(self.n_features))

    @property
    def fraction_flagged(self) -> float:
        """Return the fraction of scored samples flagged as OOD.

        Examples
        --------
        >>> import numpy as np
        >>> report = OODReport(
        ...     scores=np.array([0.1, 5.0]),
        ...     threshold=1.0,
        ...     flagged=np.array([False, True]),
        ...     method="knn",
        ...     quantile=None,
        ...     n_reference=10,
        ...     n_features=2,
        ... )
        >>> report.fraction_flagged
        0.5
        """
        return float(np.mean(self.flagged))


def flag_out_of_distribution(
    x: Any,
    reference: Any,
    *,
    method: str = "mahalanobis",
    k: int = 5,
    quantile: float = 0.99,
    threshold: float | None = None,
) -> OODReport:
    """Score inputs and flag those outside the training support.

    When ``threshold`` is not supplied, it is derived as the
    requested ``quantile`` of the reference set's own leave-one-out
    (``"knn"``) or full-sample (``"mahalanobis"``) self-scores, i.e.
    "how unusual is a typical reference point".

    Parameters
    ----------
    x : array-like, shape (n_samples, n_features)
        Inputs to score.
    reference : array-like, shape (n_reference, n_features)
        Training-set feature vectors defining the support region.
    method : {"mahalanobis", "knn"}, default="mahalanobis"
        Distance measure, as in :func:`ood_score`.
    k : int, default=5
        Neighbour rank used by ``method="knn"``. Ignored otherwise.
    quantile : float, default=0.99
        Quantile in ``(0, 1)`` of the reference self-scores used to
        derive ``threshold``. Ignored when ``threshold`` is given.
    threshold : float or None, optional
        Explicit score threshold. Overrides ``quantile`` when given.

    Returns
    -------
    OODReport
        Scores, threshold, and per-sample OOD flags.

    Examples
    --------
    >>> import numpy as np
    >>> reference = np.array(
    ...     [
    ...         [0.0, 0.0],
    ...         [1.0, 0.0],
    ...         [0.0, 1.0],
    ...         [-1.0, 0.0],
    ...         [0.0, -1.0],
    ...         [0.5, 0.5],
    ...     ]
    ... )
    >>> x = np.array([[0.0, 0.0], [50.0, 50.0]])
    >>> report = flag_out_of_distribution(x, reference, k=1)
    >>> report.method, report.n_reference
    ('mahalanobis', 6)
    """
    if method not in _METHODS:
        raise ValueError(f"method must be one of {_METHODS}.")
    x_arr = _validate_2d("x", x)
    ref_arr = _validate_2d("reference", reference)
    if x_arr.shape[1] != ref_arr.shape[1]:
        raise ValueError(
            "x and reference must share the same number of features."
        )
    scores = (
        _mahalanobis_scores(x_arr, ref_arr)
        if method == "mahalanobis"
        else _knn_scores(x_arr, ref_arr, k)
    )
    if threshold is not None:
        resolved_threshold = float(threshold)
        if not np.isfinite(resolved_threshold):
            raise ValueError("threshold must be finite.")
        resolved_quantile = None
    else:
        if not (0.0 < quantile < 1.0):
            raise ValueError("quantile must be within (0, 1).")
        self_scores = (
            _mahalanobis_scores(ref_arr, ref_arr)
            if method == "mahalanobis"
            else _self_knn_scores(ref_arr, k)
        )
        resolved_threshold = float(np.quantile(self_scores, quantile))
        resolved_quantile = float(quantile)
    return OODReport(
        scores=scores,
        threshold=resolved_threshold,
        flagged=scores > resolved_threshold,
        method=method,
        quantile=resolved_quantile,
        n_reference=ref_arr.shape[0],
        n_features=ref_arr.shape[1],
    )
