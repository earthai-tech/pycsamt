# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Abstract NeuralBackend protocol.

All concrete backends (:class:`~pycsamt.backends._torch.TorchBackend`
and :class:`~pycsamt.backends._tensorflow.TensorFlowBackend`) must
implement this interface.  The interface is intentionally minimal —
pycsamt only needs to:

1. Build a neural network from an architecture spec.
2. Run a training loop.
3. Run inference.
4. Serialise / deserialise weights to numpy dicts.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

__all__ = ["NeuralBackend", "BackendSpec"]

# Architecture specification — a JSON-serialisable dict describing a network.
BackendSpec = Dict[str, Any]


class NeuralBackend(ABC):
    """
    Abstract protocol for AI/ML backend implementations.

    Every method that touches a framework object is abstract; all
    utility logic (input validation, history formatting) lives here.
    """

    name: str  # 'torch' | 'tensorflow'

    # ─── availability ─────────────────────────────────────────────────────

    @classmethod
    @abstractmethod
    def is_available(cls) -> bool:
        """Return ``True`` if the framework is importable."""

    @classmethod
    def require(cls) -> None:
        """Raise :exc:`ImportError` with install instructions if unavailable."""
        if not cls.is_available():
            pkg = "torch" if cls.name == "torch" else "tensorflow"
            raise ImportError(
                f"Backend {cls.name!r} requires {pkg}.\n"
                f"Install it with:  pip install {pkg}"
            )

    # ─── device ───────────────────────────────────────────────────────────

    @abstractmethod
    def resolve_device(self, device: Optional[str] = None) -> str:
        """
        Return a concrete device string.

        For PyTorch: ``'cuda'``, ``'mps'``, or ``'cpu'``.
        For TensorFlow: ``'/GPU:0'``, ``'/CPU:0'``, etc.
        """

    # ─── model building ───────────────────────────────────────────────────

    @abstractmethod
    def build(self, spec: BackendSpec) -> Any:
        """
        Instantiate a neural network from *spec*.

        Parameters
        ----------
        spec : dict
            Must contain ``'arch'`` key (``'cnn1d'``, ``'resnet1d'``,
            ``'fcn1d'``, ``'unet2d'``, ``'drcnn'``).
            Architecture-specific hyper-parameters are passed as
            additional keys (``n_features``, ``n_out``, ``channels``,
            ``dropout``, …).

        Returns
        -------
        model : framework-native model object
        """

    # ─── training ─────────────────────────────────────────────────────────

    @abstractmethod
    def train(
        self,
        model: Any,
        X_train: np.ndarray,
        y_train: np.ndarray,
        X_val: np.ndarray,
        y_val: np.ndarray,
        *,
        epochs: int = 100,
        batch_size: int = 256,
        lr: float = 1e-3,
        patience: int = 20,
        grad_clip: Optional[float] = None,
        device: Optional[str] = None,
        verbose: bool = True,
        loss: str = "mse",
    ) -> Dict[str, List[float]]:
        """
        Run a training loop.

        Parameters
        ----------
        model : framework-native model
            Must have been created by :meth:`build`.
        X_train, y_train : ndarray
            Training data.
        X_val, y_val : ndarray
            Validation data.
        epochs, batch_size, lr, patience, grad_clip, device, verbose
            Standard hyper-parameters.
        loss : {'mse', 'masked_mse', 'cross_entropy'}

        Returns
        -------
        history : dict
            Keys ``'train_loss'`` and ``'val_loss'`` (lists of floats).
        """

    # ─── inference ────────────────────────────────────────────────────────

    @abstractmethod
    def predict(
        self,
        model: Any,
        X: np.ndarray,
        *,
        device: Optional[str] = None,
        batch_size: int = 256,
    ) -> np.ndarray:
        """
        Run inference.

        Parameters
        ----------
        model : framework-native model
        X : ndarray (n_samples, n_features)
        device : str or None
        batch_size : int

        Returns
        -------
        y_pred : ndarray (n_samples, n_out)
        """

    # ─── serialisation ────────────────────────────────────────────────────

    @abstractmethod
    def get_weights(self, model: Any) -> Dict[str, np.ndarray]:
        """
        Return model weights as a ``name → numpy`` array dict.

        Must be invertible via :meth:`set_weights`.
        """

    @abstractmethod
    def set_weights(self, model: Any, weights: Dict[str, np.ndarray]) -> None:
        """Restore weights from a ``name → numpy`` array dict."""

    # ─── dunder ───────────────────────────────────────────────────────────

    def __repr__(self) -> str:
        ok = "available" if self.is_available() else "NOT installed"
        return f"{self.__class__.__name__}({ok})"
