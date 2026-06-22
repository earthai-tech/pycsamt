# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Abstract base classes for all pycsamt AI/ML estimators.

Every neural-network–based estimator in :mod:`pycsamt.ai` inherits
from :class:`BaseEMNet` and every ML-based data processor inherits
from :class:`BaseEMProcessor`.  Both follow the sklearn interface
(``fit`` / ``predict`` / ``score``) augmented with EM-specific helpers
(``save`` / ``load`` / ``from_pretrained``).

No deep-learning framework is imported here — torch / tensorflow are
lazily imported inside concrete subclasses so that ``import pycsamt``
never fails for users who have not installed them.

Design principles
-----------------
* **Framework-agnostic** — the base declares the interface; subclasses
  choose PyTorch or TensorFlow.
* **Sklearn-compatible** — ``fit(X, y)``, ``predict(X)``, ``score(X, y)``
  work with numpy arrays so pipelines and CV are easy.
* **EM-aware** — helpers to accept :class:`~pycsamt.z.z.Z` objects and
  :class:`~pycsamt.forward.em1d.ForwardResponse` directly.
* **Serialisable** — ``save`` / ``load`` handle both model weights and
  hyperparameters in a single checkpoint file.
"""
from __future__ import annotations

import json
import warnings
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Tuple, Union

import numpy as np

from ..api.property import PyCSAMTObject

__all__ = [
    "BaseEMNet",
    "BaseEMProcessor",
    "BasePINNInverter",
    "BaseHybridInverter",
    "EMCheckpoint",
]

# ─────────────────────────────────────────────────────────────────────────────
# Checkpoint container
# ─────────────────────────────────────────────────────────────────────────────

class EMCheckpoint:
    """
    Minimal checkpoint for saving/loading EM network state.

    Holds hyperparameters (serialised as JSON) alongside framework-
    specific weight data (numpy arrays or a torch state dict path).

    Parameters
    ----------
    params : dict
        Hyperparameters (must be JSON-serialisable).
    weights : dict or None
        Framework weight dict (numpy arrays).
    history : dict or None
        Training history: ``{'train_loss': [...], 'val_loss': [...]}``.
    meta : dict or None
        Extra metadata (solver, freq range, training date, …).
    """

    def __init__(
        self,
        params: Dict[str, Any],
        weights: Optional[Dict[str, np.ndarray]] = None,
        history: Optional[Dict[str, list]] = None,
        meta: Optional[Dict[str, Any]] = None,
    ):
        self.params = params
        self.weights = weights or {}
        self.history = history or {}
        self.meta = meta or {}

    def save(self, path: Union[str, Path]) -> None:
        """Save to a ``.npz`` checkpoint file."""
        path = Path(path)
        arrays = {f"w_{k}": v for k, v in self.weights.items()}
        arrays["_params_json"] = np.array(json.dumps(self.params))
        arrays["_history_json"] = np.array(json.dumps(self.history))
        arrays["_meta_json"] = np.array(json.dumps(self.meta))
        np.savez_compressed(path, **arrays)

    @classmethod
    def load(cls, path: Union[str, Path]) -> "EMCheckpoint":
        """Load from a ``.npz`` checkpoint file."""
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Checkpoint not found: {path}")
        d = np.load(str(path), allow_pickle=True)
        params = json.loads(str(d["_params_json"]))
        history = json.loads(str(d["_history_json"]))
        meta = json.loads(str(d["_meta_json"]))
        weights = {
            k[2:]: d[k] for k in d.files
            if k.startswith("w_")
        }
        return cls(params=params, weights=weights, history=history, meta=meta)


# ─────────────────────────────────────────────────────────────────────────────
# BaseEMNet — neural-network–based inversion / processing
# ─────────────────────────────────────────────────────────────────────────────

class BaseEMNet(ABC):
    """
    Abstract base class for EM neural network estimators.

    Subclasses must implement :meth:`fit`, :meth:`predict`, and
    :meth:`_build_network`.

    All ``fit`` / ``predict`` calls accept either numpy arrays
    (the standard path) **or** EM domain objects (Z tensors,
    ForwardResponse, SiteCollection) which are coerced to arrays via
    the :meth:`_coerce_input` hook.

    Parameters
    ----------
    arch : str
        Architecture identifier (e.g. ``'cnn1d'``, ``'resnet'``).
    n_layers : int
        Number of earth layers to invert for.
    solver : str
        The EM method this network targets (``'mt1d'``, ``'tem1d'``,
        ``'csamt1d'``).
    device : str or None
        Compute device (``'cpu'``, ``'cuda'``, ``'mps'``).  If ``None``
        uses CUDA when available.
    """

    def __init__(
        self,
        arch: str = "cnn1d",
        n_layers: int = 5,
        solver: str = "mt1d",
        device: Optional[str] = None,
    ):
        self.arch = arch
        self.n_layers = n_layers
        self.solver = solver
        self.device = device
        self._network = None
        self._is_fitted: bool = False
        self._history: Dict[str, list] = {}
        self._params: Dict[str, Any] = {}
        self._meta: Dict[str, Any] = {}

    # ─── abstract interface ───────────────────────────────────────────────

    @abstractmethod
    def _build_network(self) -> Any:
        """
        Instantiate and return the neural network.

        Called once at the start of :meth:`fit`.
        Must *not* import torch/tensorflow at module level.
        """

    @abstractmethod
    def fit(
        self,
        X: Union[np.ndarray, str, "ForwardDataset"],
        y: Optional[np.ndarray] = None,
        **kwargs,
    ) -> "BaseEMNet":
        """
        Train the network.

        Parameters
        ----------
        X : ndarray (n_samples, n_features) or path-like or ForwardDataset
            Training data.  A string is treated as a path to a ``.npz``
            file produced by :func:`~pycsamt.forward.batch.generate_dataset`.
        y : ndarray (n_samples, n_params) or None
            Targets.  If *X* is a dataset object or path, *y* is ignored.
        **kwargs
            Solver-specific keyword arguments (e.g. ``epochs``,
            ``batch_size``, ``lr``).

        Returns
        -------
        self
        """

    @abstractmethod
    def predict(
        self,
        X: Union[np.ndarray, "Z", Sequence],
    ) -> np.ndarray:
        """
        Predict subsurface model parameters from data.

        Parameters
        ----------
        X : ndarray (n_samples, n_features) or Z object or list of Z
            Input data in array form or as EM domain objects.

        Returns
        -------
        y_pred : ndarray (n_samples, n_params)
            Predicted log₁₀(ρ) and thickness values.
        """

    # ─── sklearn-compatible helpers ──────────────────────────────────────

    def score(
        self,
        X: np.ndarray,
        y: np.ndarray,
        *,
        metric: str = "rmse",
    ) -> float:
        """
        Evaluate prediction quality.

        Parameters
        ----------
        X : ndarray
            Input features.
        y : ndarray
            True model parameters.
        metric : {'rmse', 'r2', 'mae', 'relative_rmse'}
            Evaluation metric.

        Returns
        -------
        score : float
            Lower is better for RMSE/MAE; higher is better for R².
        """
        y_pred = self.predict(X)
        return _compute_metric(y, y_pred, metric)

    def fit_predict(
        self,
        X: np.ndarray,
        y: np.ndarray,
        **kwargs,
    ) -> np.ndarray:
        """Fit then predict on the same data (useful for diagnostics)."""
        self.fit(X, y, **kwargs)
        return self.predict(X)

    # ─── serialisation ────────────────────────────────────────────────────

    def save(self, path: Union[str, Path]) -> None:
        """
        Save the model weights and hyperparameters to *path*.

        Parameters
        ----------
        path : str or Path
            Destination file path (a ``.pt`` extension is used for
            PyTorch backends; ``.npz`` for numpy-only models).
        """
        ckpt = EMCheckpoint(
            params=self._get_params(),
            weights=self._get_weights(),
            history=self._history,
            meta=self._meta,
        )
        ckpt.save(path)

    @classmethod
    def load(cls, path: Union[str, Path]) -> "BaseEMNet":
        """
        Load a model saved with :meth:`save`.

        Parameters
        ----------
        path : str or Path

        Returns
        -------
        estimator : instance of the calling class
        """
        ckpt = EMCheckpoint.load(path)
        obj = cls.__new__(cls)
        obj.__init__(**ckpt.params)
        obj._history = ckpt.history
        obj._meta = ckpt.meta
        if ckpt.weights:
            obj._load_weights(ckpt.weights)
            obj._is_fitted = True
        return obj

    @classmethod
    def from_pretrained(cls, name: str) -> "BaseEMNet":
        """
        Load a pre-trained model from the pycsamt model zoo.

        .. note::
            Pre-trained weights are hosted separately and will be
            downloaded on first call.  This feature is scheduled for
            Phase 5 of the pycsamt AI roadmap.

        Parameters
        ----------
        name : str
            Model identifier, e.g. ``'mt1d-resnet-v1'``.

        Raises
        ------
        NotImplementedError
            Until Phase 5 weights are published.
        """
        raise NotImplementedError(
            "Pre-trained model zoo is scheduled for Phase 5 of the "
            "pycsamt AI roadmap.  To use a custom model, call "
            "cls.load(path) with a checkpoint file."
        )

    # ─── EM domain input coercion ─────────────────────────────────────────

    def _coerce_input(
        self,
        X: Any,
        *,
        include_phase: bool = True,
        log_rho: bool = True,
    ) -> np.ndarray:
        """
        Convert EM domain objects to a numpy feature matrix.

        Accepts:
        * ndarray — returned as-is (cast to float32)
        * :class:`~pycsamt.z.z.Z` — single Z object
        * list of Z — multiple sites
        * :class:`~pycsamt.forward.em1d.ForwardResponse` — single response
        """
        if isinstance(X, np.ndarray):
            return X.astype(np.float32)

        # Single or list of Z objects
        try:
            from pycsamt.z.z import Z as ZClass
            if isinstance(X, ZClass):
                X = [X]
            if isinstance(X, (list, tuple)) and len(X) > 0 and isinstance(X[0], ZClass):
                return _z_list_to_array(X, include_phase=include_phase, log_rho=log_rho)
        except ImportError:
            pass

        # ForwardResponse
        try:
            from pycsamt.forward.em1d import ForwardResponse
            if isinstance(X, ForwardResponse):
                return X.to_array(
                    log_rho=log_rho, include_phase=include_phase
                ).reshape(1, -1).astype(np.float32)
        except ImportError:
            pass

        raise TypeError(
            f"Cannot coerce input of type {type(X)!r} to a numpy array.  "
            "Provide an ndarray, Z object, list of Z, or ForwardResponse."
        )

    # ─── hooks for subclasses ─────────────────────────────────────────────

    def _get_params(self) -> Dict[str, Any]:
        """Return JSON-serialisable hyperparameter dict."""
        return {
            "arch": self.arch,
            "n_layers": self.n_layers,
            "solver": self.solver,
            "device": self.device,
        }

    def _get_weights(self) -> Dict[str, np.ndarray]:
        """Return a dict of numpy weight arrays.  Override in subclass."""
        return {}

    def _load_weights(self, weights: Dict[str, np.ndarray]) -> None:
        """Restore weights from a dict of numpy arrays.  Override in subclass."""
        pass

    def _resolve_device(self) -> str:
        """Pick compute device via the active backend."""
        from ._backend_utils import resolve_device
        return resolve_device(self.device)

    # ─── dunder ───────────────────────────────────────────────────────────

    def __repr__(self) -> str:
        fitted = "fitted" if self._is_fitted else "unfitted"
        return (
            f"{self.__class__.__name__}("
            f"arch={self.arch!r}, n_layers={self.n_layers}, "
            f"solver={self.solver!r}, {fitted})"
        )


# ─────────────────────────────────────────────────────────────────────────────
# BaseEMProcessor — ML-based data processing (denoising, QC, …)
# ─────────────────────────────────────────────────────────────────────────────

class BaseEMProcessor(ABC):
    """
    Abstract base class for ML-based EM data processing tasks.

    Unlike :class:`BaseEMNet`, processors transform data to data rather
    than data to model parameters.

    Subclasses in :mod:`pycsamt.ai.processing` include
    :class:`~pycsamt.ai.processing.denoise.EMDenoiser`,
    :class:`~pycsamt.ai.processing.qc.EMQCScorer`, and
    :class:`~pycsamt.ai.processing.anomaly.AnomalyDetector`.
    """

    @abstractmethod
    def fit(self, X: np.ndarray, **kwargs) -> "BaseEMProcessor":
        """Fit the processor to training data."""

    @abstractmethod
    def transform(self, X: np.ndarray) -> np.ndarray:
        """Apply the learned transformation."""

    def fit_transform(self, X: np.ndarray, **kwargs) -> np.ndarray:
        """Fit then transform in one call."""
        return self.fit(X, **kwargs).transform(X)

    def save(self, path: Union[str, Path]) -> None:
        """Save to a checkpoint file."""
        ckpt = EMCheckpoint(
            params=self._get_params(),
            weights=self._get_weights(),
        )
        ckpt.save(path)

    @classmethod
    def load(cls, path: Union[str, Path]) -> "BaseEMProcessor":
        ckpt = EMCheckpoint.load(path)
        obj = cls.__new__(cls)
        obj.__init__(**ckpt.params)
        if ckpt.weights:
            obj._load_weights(ckpt.weights)
        return obj

    @classmethod
    def load_pretrained(cls, name: str) -> "BaseEMProcessor":
        raise NotImplementedError(
            "Pre-trained processor weights are scheduled for Phase 5."
        )

    def _get_params(self) -> Dict[str, Any]:
        return {}

    def _get_weights(self) -> Dict[str, np.ndarray]:
        return {}

    def _load_weights(self, weights: Dict[str, np.ndarray]) -> None:
        pass

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"


# ─────────────────────────────────────────────────────────────────────────────
# BasePINNInverter — physics-informed parameter optimisation
# ─────────────────────────────────────────────────────────────────────────────


class BasePINNInverter(PyCSAMTObject, ABC):
    r"""
    Abstract base for physics-informed EM inverters.

    Unlike :class:`BaseEMNet`, which wraps a supervised
    neural network, PINN inverters minimise a
    physics-informed loss directly in model-parameter
    space — no network weights are trained.

    Both PyTorch and TensorFlow are supported via the
    :mod:`pycsamt.backends` dispatch layer.  The active
    backend is selected automatically; override with
    ``PYCSAMT_AI_BACKEND`` or
    :func:`pycsamt.backends.set_backend`.

    Subclasses must implement :meth:`fit` and
    :meth:`residuals`.

    Parameters
    ----------
    n_layers : int
        Number of earth layers.
    depth_max : float
        Maximum investigation depth in metres.
    device : str or None
        Compute device (``'cpu'``, ``'cuda'``,
        ``'/GPU:0'``, …).  Auto-detected if ``None``.
    """

    def __init__(
        self,
        n_layers: int = 5,
        depth_max: float = 1000.0,
        device: Optional[str] = None,
    ) -> None:
        self.n_layers = int(n_layers)
        self.depth_max = float(depth_max)
        self.device = device
        self._is_fitted: bool = False
        self._history: list = []

    # ── abstract interface ────────────────────────────

    @abstractmethod
    def fit(self, **kwargs) -> "BasePINNInverter":
        """
        Minimise the physics-informed loss.

        Returns
        -------
        self
        """

    @abstractmethod
    def residuals(self) -> Any:
        """
        Return observed vs predicted data.

        Returns
        -------
        pandas.DataFrame
        """

    # ── common helpers ────────────────────────────────

    def convergence_curve(self) -> list:
        """
        Return Adam loss history.

        Returns
        -------
        list of float
        """
        self._check_fitted()
        return list(self._history)

    def _resolve_device(self) -> str:
        """Pick compute device via active backend."""
        from ._backend_utils import resolve_device
        return resolve_device(self.device)

    def _require_backend(self) -> str:
        """
        Return active backend name; raise if none.

        Returns
        -------
        name : str
            ``'torch'`` or ``'tensorflow'``.

        Raises
        ------
        ImportError
            When no DL framework is installed.
        """
        from ._backend_utils import active_backend
        name = active_backend()
        if name == "none":
            raise ImportError(
                "PINNInverter requires PyTorch or "
                "TensorFlow.\n"
                "Install one with:\n"
                "    pip install torch\n"
                "    pip install tensorflow"
            )
        return name

    def _check_fitted(self) -> None:
        """Raise RuntimeError if not fitted."""
        if not self._is_fitted:
            raise RuntimeError(
                "Call fit() before accessing results."
            )


# ─────────────────────────────────────────────────────────────────────────────
# BaseHybridInverter — AI warm-start + physics refinement
# ─────────────────────────────────────────────────────────────────────────────


class BaseHybridInverter(BasePINNInverter):
    r"""
    Abstract base for two-stage hybrid inverters.

    Stage 1 applies a pre-trained AI inverter to
    produce an initial Earth model.  Stage 2 refines
    that model with the same physics-informed Adam
    optimisation used by :class:`BasePINNInverter`.

    Subclasses must implement :meth:`_load_ai_inverter`,
    :meth:`_run_stage1`, :meth:`fit`, and
    :meth:`residuals`.

    Parameters
    ----------
    n_layers : int
        Number of earth layers.
    depth_max : float
        Maximum investigation depth in metres.
    device : str or None
        Compute device.  Auto-detected if ``None``.
    """

    def __init__(
        self,
        n_layers: int = 5,
        depth_max: float = 1000.0,
        device: Optional[str] = None,
    ) -> None:
        super().__init__(
            n_layers=n_layers,
            depth_max=depth_max,
            device=device,
        )
        self._ai_inv: Any = None
        self._stage1_fitted: bool = False

    @abstractmethod
    def _load_ai_inverter(
        self, ai_inverter: Any
    ) -> Any:
        """
        Validate and load the AI inverter.

        Returns
        -------
        ready-to-use inverter object
        """

    @abstractmethod
    def _run_stage1(self, *, verbose: bool) -> Any:
        """Apply the AI inverter to all stations."""


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

def _compute_metric(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    metric: str,
) -> float:
    """Compute a scalar evaluation metric, ignoring NaN entries."""
    mask = np.isfinite(y_true) & np.isfinite(y_pred)
    yt = y_true[mask]
    yp = y_pred[mask]
    if len(yt) == 0:
        return np.nan
    m = metric.lower()
    if m == "rmse":
        return float(np.sqrt(np.mean((yt - yp) ** 2)))
    if m == "mae":
        return float(np.mean(np.abs(yt - yp)))
    if m == "r2":
        ss_res = np.sum((yt - yp) ** 2)
        ss_tot = np.sum((yt - np.mean(yt)) ** 2)
        return float(1.0 - ss_res / (ss_tot + 1e-12))
    if m == "relative_rmse":
        rel = (yt - yp) / (np.abs(yt) + 1e-12)
        return float(np.sqrt(np.mean(rel ** 2)))
    raise ValueError(
        f"Unknown metric {metric!r}. "
        "Use 'rmse', 'mae', 'r2', or 'relative_rmse'."
    )


def _z_list_to_array(
    z_list: list,
    *,
    include_phase: bool = True,
    log_rho: bool = True,
) -> np.ndarray:
    """Convert a list of Z objects to a (n_sites, n_features) numpy array."""
    rows = []
    for z in z_list:
        rho_xy = z.resistivity_xy
        phs_xy = z.phase_xy
        if rho_xy is None or phs_xy is None:
            raise ValueError(
                "Z object has no computed resistivity/phase.  "
                "Call z.compute_resistivity_phase() first."
            )
        rho = np.log10(np.maximum(rho_xy, 1e-12)) if log_rho else rho_xy
        row = np.concatenate([rho, phs_xy]) if include_phase else rho
        rows.append(row)
    return np.vstack(rows).astype(np.float32)
