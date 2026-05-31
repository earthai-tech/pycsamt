# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
End-to-end 1-D EM inversion workflow.

:class:`EMInverter1D` is the primary user-facing estimator for 1-D
neural-network inversion of MT, CSAMT, and TEM data.  It:

1. Loads a :class:`~pycsamt.forward.batch.ForwardDataset` (or a
   ``.npz`` path to one).
2. Normalises features and targets with z-score normalisation
   (log₁₀ thickness transform applied automatically).
3. Instantiates the requested architecture
   (``'cnn1d'``, ``'resnet'``, ``'fcn'``).
4. Trains with :class:`~pycsamt.ai.training.trainer.EMTrainer`
   (early stopping, LR scheduling, masked MSE loss).
5. Predicts on new data as:
   * numpy arrays
   * lists of :class:`~pycsamt.z.z.Z` objects
   * :class:`~pycsamt.forward.em1d.ForwardResponse` objects

The sklearn-compatible interface (``fit`` / ``predict`` / ``score``)
from :class:`~pycsamt.ai._base.BaseEMNet` is preserved.

Quick start
-----------
>>> import numpy as np
>>> from pycsamt.forward.batch import generate_dataset
>>> from pycsamt.ai.inversion.inv1d import EMInverter1D

# Generate 2 000 training samples
>>> ds = generate_dataset(n_samples=2_000, seed=0, n_layers=5)

# Train
>>> inv = EMInverter1D(arch="resnet", n_layers=5, solver="mt1d")
>>> inv.fit(ds, epochs=30, batch_size=128, verbose=True)

# Predict on a new forward response
>>> from pycsamt.forward import MT1DForward, LayeredModel
>>> m_new = LayeredModel.random(n_layers=5, seed=99)
>>> resp = MT1DForward(np.logspace(-3, 4, 30)).run(m_new)
>>> y_pred = inv.predict_response(resp)  # → LayeredModel
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Union

import numpy as np

from .._base import BaseEMNet, EMCheckpoint

__all__ = ["EMInverter1D"]

# Architecture registry: name → factory callable (n_features, n_out, **kw)
_ARCH_REGISTRY: Dict[str, Any] = {}


# ─────────────────────────────────────────────────────────────────────────────
# from_pretrained helper (module-level to avoid class-level redefinition)
# ─────────────────────────────────────────────────────────────────────────────

def _load_pretrained(name: str, cache_dir: Optional[str] = None) -> "EMInverter1D":
    """Download (if needed) and load a pre-trained checkpoint."""
    from .._zoo import download_checkpoint, get_pretrained_info
    info = get_pretrained_info(name)
    fpath = download_checkpoint(name, cache_dir=cache_dir)
    obj = EMInverter1D.load(fpath)
    return obj


def _register(name, fn):
    _ARCH_REGISTRY[name.lower()] = fn


def _lazy_register():
    """Populate registry on first use (avoids eager torch import)."""
    if _ARCH_REGISTRY:
        return

    def _cnn1d(n_features, n_out, **kw):
        from ..nets.cnn1d import CNN1DNet
        return CNN1DNet(n_features, n_out, **kw).build()

    def _resnet(n_features, n_out, **kw):
        from ..nets.resnet import ResNet1DNet
        return ResNet1DNet(n_features, n_out, **kw).build()

    def _fcn(n_features, n_out, **kw):
        from ..nets.fcn import FCN1DNet
        return FCN1DNet(n_features, n_out, **kw).build()

    _register("cnn1d", _cnn1d)
    _register("cnn", _cnn1d)
    _register("resnet", _resnet)
    _register("resnet1d", _resnet)
    _register("fcn", _fcn)
    _register("fcn1d", _fcn)


# ─────────────────────────────────────────────────────────────────────────────
# EMInverter1D
# ─────────────────────────────────────────────────────────────────────────────

class EMInverter1D(BaseEMNet):
    """
    1-D EM neural-network inverter.

    Supports MT, CSAMT (far-field), and TEM step-off data.

    Parameters
    ----------
    arch : {'resnet', 'cnn1d', 'fcn'}
        Network architecture.  ``'resnet'`` (Liu 2021 style) gives the
        best accuracy on typical MT datasets.  ``'cnn1d'`` (Puzyrev
        style) is faster to train.  ``'fcn'`` (Moghadas style) handles
        variable input length.
    n_layers : int
        Number of earth layers to invert for.
    solver : {'mt1d', 'csamt1d', 'tem1d'}
        EM method this inverter targets (determines default frequency grid
        for Z-object input coercion).
    device : str or None
        Compute device.  Auto-detects CUDA/MPS/CPU if ``None``.
    log_thickness : bool
        Apply log₁₀ to thickness in training targets.  Strongly
        recommended when thicknesses span > 2 orders of magnitude.
    include_phase : bool
        Include impedance phase in the input feature vector for MT/CSAMT.
        Setting to ``False`` halves the input size.
    augment_noise : float
        On-the-fly noise level added to training inputs each epoch.
    net_kwargs : dict
        Extra keyword arguments forwarded to the network constructor
        (e.g. ``channels=(64, 128, 256)`` for ResNet).
    """

    def __init__(
        self,
        arch: str = "resnet",
        n_layers: int = 5,
        solver: str = "mt1d",
        *,
        device: Optional[str] = None,
        log_thickness: bool = True,
        include_phase: bool = True,
        augment_noise: float = 0.02,
        **net_kwargs,
    ):
        super().__init__(arch=arch, n_layers=n_layers, solver=solver, device=device)
        self.log_thickness = log_thickness
        self.include_phase = include_phase
        self.augment_noise = augment_noise
        self._net_kwargs = net_kwargs

        # Set after fit
        self._em_dataset = None    # EMDataset (holds normalizers)
        self._n_features: Optional[int] = None
        self._n_out: int = 2 * n_layers - 1

    # ─── fit ──────────────────────────────────────────────────────────────

    def fit(
        self,
        X,
        y=None,
        *,
        epochs: int = 100,
        batch_size: int = 256,
        lr: float = 1e-3,
        patience: int = 20,
        val_frac: float = 0.1,
        grad_clip: Optional[float] = 1.0,
        seed: Optional[int] = None,
        verbose: bool = True,
    ) -> "EMInverter1D":
        """
        Train the inverter on a :class:`~pycsamt.forward.batch.ForwardDataset`
        or a ``.npz`` file path.

        Parameters
        ----------
        X : ForwardDataset or str or Path
            Training data.  If *y* is also given, *X* and *y* are treated
            as raw feature / target numpy arrays.
        y : ndarray or None
            Raw targets (only when X is an ndarray).
        epochs : int
            Maximum training epochs.
        batch_size : int
        lr : float
            Initial learning rate.
        patience : int
            Early-stopping patience.
        val_frac : float
            Fraction of training data used for validation.
        grad_clip : float or None
            Gradient-norm clipping threshold.
        seed : int or None
            Random seed for train/val split.
        verbose : bool
            Print per-epoch summary.

        Returns
        -------
        self
        """
        from ..training.dataset import EMDataset
        from ..training.trainer import EMTrainer

        # ── Load data ────────────────────────────────────────────────────
        ds = self._load_dataset(X, y)

        # ── Build EMDataset (normalise) ──────────────────────────────────
        em_ds = EMDataset(
            ds,
            n_layers=self.n_layers,
            log_thickness=self.log_thickness,
            augment_noise=self.augment_noise,
        )
        train_ds, val_ds = em_ds.split(val_frac=val_frac, seed=seed)
        self._em_dataset = em_ds
        self._n_features = em_ds.n_features
        self._n_out = em_ds.n_params

        # ── Build network ────────────────────────────────────────────────
        self._network = self._build_network()

        # ── Train ────────────────────────────────────────────────────────
        device = self._resolve_device()
        trainer = EMTrainer(
            self._network,
            lr=lr,
            patience=patience,
            batch_size=batch_size,
            device=device,
            grad_clip=grad_clip,
            verbose=verbose,
        )
        trainer.fit(train_ds, val_ds, epochs=epochs)
        self._history = trainer.history
        self._meta["best_epoch"] = trainer.best_epoch
        self._meta["best_val_loss"] = trainer.best_val_loss
        self._meta["n_features"] = self._n_features
        self._meta["n_out"] = self._n_out
        self._is_fitted = True
        return self

    # ─── predict ──────────────────────────────────────────────────────────

    def predict(
        self,
        X,
        *,
        as_log_rho: bool = True,
    ) -> np.ndarray:
        """
        Predict model parameters for new data.

        Parameters
        ----------
        X : ndarray (n_samples, n_features) | Z | list of Z | ForwardResponse
            Input data.
        as_log_rho : bool
            If ``True`` (default), the returned resistivity values are
            in log₁₀(Ω·m) space matching the training targets.
            If ``False``, resistivities are back-transformed to Ω·m.

        Returns
        -------
        y_pred : ndarray, shape (n_samples, n_params)
            Parameter vector: ``[log10(ρ), log10(h)]`` (or linear if
            ``as_log_rho=False``).
        """
        self._check_fitted()
        import torch

        X_arr = self._coerce_input(
            X, include_phase=self.include_phase, log_rho=True
        )
        X_norm = self._em_dataset.x_norm.transform(X_arr)

        device = self._resolve_device()
        net = self._network.to(device).eval()

        with torch.no_grad():
            t = torch.from_numpy(X_norm.astype(np.float32)).to(device)
            y_norm = net(t).cpu().numpy()

        # Inverse-normalise → [log10(ρ), log10(h)] or [log10(ρ), h]
        y_raw = self._em_dataset.y_norm.inverse_transform(y_norm)
        if not as_log_rho:
            # Convert log10(ρ) → ρ for the first n_layers columns
            n = self.n_layers
            y_raw[:, :n] = 10.0 ** y_raw[:, :n]
            if self.log_thickness:
                y_raw[:, n:] = 10.0 ** y_raw[:, n:]
        return y_raw

    def predict_models(self, X) -> list:
        """
        Predict and return a list of :class:`~pycsamt.forward.synthetic.LayeredModel`.

        Resistivity and thickness are back-transformed from log space.
        """
        from pycsamt.forward.synthetic import LayeredModel

        y_pred = self.predict(X, as_log_rho=True)
        models = []
        for row in y_pred:
            n = self.n_layers
            rho = 10.0 ** row[:n]
            thick = 10.0 ** row[n:] if self.log_thickness else row[n:]
            try:
                m = LayeredModel(
                    resistivity=np.maximum(rho, 1e-3),
                    thickness=np.maximum(thick, 1.0),
                )
            except Exception:
                m = None
            models.append(m)
        return models

    def predict_response(self, response) -> "LayeredModel":
        """
        Convenience: invert a single
        :class:`~pycsamt.forward.em1d.ForwardResponse` and return the
        predicted :class:`~pycsamt.forward.synthetic.LayeredModel`.
        """
        X = self._coerce_input(
            response, include_phase=self.include_phase, log_rho=True
        )
        models = self.predict_models(X)
        return models[0]

    # ─── serialisation ────────────────────────────────────────────────────

    def save(self, path: Union[str, Path]) -> None:
        """Save weights + normaliser + hyperparameters to *path*."""
        import torch
        device = self._resolve_device()
        weights = {
            k: v.cpu().numpy()
            for k, v in self._network.to(device).state_dict().items()
        }
        # Store normaliser as JSON-friendly dicts
        meta = dict(self._meta)
        if self._em_dataset is not None:
            meta["x_norm"] = self._em_dataset.x_norm.to_dict()
            meta["y_norm"] = self._em_dataset.y_norm.to_dict()
            meta["n_layers_ds"] = self._em_dataset._n_layers
            meta["log_thickness"] = self._em_dataset._log_thickness

        ckpt = EMCheckpoint(
            params=self._get_params(),
            weights=weights,
            history=self._history,
            meta=meta,
        )
        ckpt.save(path)

    @classmethod
    def load(cls, path: Union[str, Path]) -> "EMInverter1D":
        """Load a saved inverter from *path*."""
        import torch
        from ..training.dataset import Normalizer, EMDataset

        ckpt = EMCheckpoint.load(path)
        p = ckpt.params
        obj = cls(
            arch=p["arch"],
            n_layers=p["n_layers"],
            solver=p["solver"],
            device=p["device"],
            log_thickness=p.get("log_thickness", True),
            include_phase=p.get("include_phase", True),
        )
        obj._n_features = ckpt.meta.get("n_features")
        obj._n_out = ckpt.meta.get("n_out", 2 * obj.n_layers - 1)

        # Restore normalisers
        if "x_norm" in ckpt.meta:
            obj._em_dataset = EMDataset.__new__(EMDataset)
            obj._em_dataset.x_norm = Normalizer.from_dict(ckpt.meta["x_norm"])
            obj._em_dataset.y_norm = Normalizer.from_dict(ckpt.meta["y_norm"])
            obj._em_dataset._n_layers = ckpt.meta.get("n_layers_ds")
            obj._em_dataset._log_thickness = ckpt.meta.get("log_thickness", True)

        # Rebuild and load network
        obj._network = obj._build_network()
        state = {k: torch.from_numpy(v) for k, v in ckpt.weights.items()}
        obj._network.load_state_dict(state)
        obj._history = ckpt.history
        obj._meta = ckpt.meta
        obj._is_fitted = bool(ckpt.weights)
        return obj

    # ─── model zoo ────────────────────────────────────────────────────────

    @classmethod
    def from_pretrained(
        cls,
        name: str,
        *,
        cache_dir: Optional[str] = None,
    ) -> "EMInverter1D":
        """
        Load a pre-trained model from the pycsamt model zoo.

        Pre-trained weights are hosted at
        https://github.com/earthai-tech/pycsamt-models
        and are downloaded to ``~/.pycsamt/model_zoo/`` on first call.

        Parameters
        ----------
        name : str
            Model identifier.  Call
            :func:`~pycsamt.ai._zoo.list_pretrained` to see available
            models.
        cache_dir : str or None
            Override the default download directory.

        Returns
        -------
        EMInverter1D

        Raises
        ------
        KeyError
            If *name* is not in the model zoo registry.
        RuntimeError
            If the download fails (weights not yet publicly available).

        Examples
        --------
        >>> from pycsamt.ai.inversion import EMInverter1D
        >>> inv = EMInverter1D.from_pretrained("mt1d-resnet-5layer-v1")  # doctest: +SKIP
        """
        return _load_pretrained(name, cache_dir=cache_dir)

    # ─── BaseEMNet hooks ──────────────────────────────────────────────────

    def _build_network(self):
        """Instantiate the PyTorch module from the chosen architecture."""
        if self._n_features is None:
            raise RuntimeError(
                "Cannot build network before n_features is known.  "
                "Call fit() first (n_features is inferred from training data)."
            )
        _lazy_register()
        arch = self.arch.lower()
        if arch not in _ARCH_REGISTRY:
            raise ValueError(
                f"Unknown architecture {self.arch!r}.  "
                f"Available: {list(_ARCH_REGISTRY)}"
            )
        net = _ARCH_REGISTRY[arch](
            self._n_features, self._n_out, **self._net_kwargs
        )
        return net

    def _get_params(self) -> dict:
        p = super()._get_params()
        p["log_thickness"] = self.log_thickness
        p["include_phase"] = self.include_phase
        p["augment_noise"] = self.augment_noise
        return p

    # ─── internals ────────────────────────────────────────────────────────

    def _load_dataset(self, X, y):
        """Return a ForwardDataset from X (which may be a path, a ds, or arrays)."""
        from pycsamt.forward.batch import ForwardDataset
        if isinstance(X, ForwardDataset):
            return X
        if isinstance(X, (str, Path)):
            return ForwardDataset.load(str(X))
        if isinstance(X, np.ndarray):
            if y is None:
                raise ValueError(
                    "When X is a numpy array, y must also be provided."
                )
            n = len(X)
            meta = np.zeros(n, dtype=[("n_layers", "i4"), ("noise_level", "f4")])
            meta["n_layers"] = self.n_layers
            return ForwardDataset(
                X=X.astype(np.float32),
                y=y.astype(np.float32),
                meta=meta,
                solver=self.solver,
            )
        raise TypeError(
            f"X must be a ForwardDataset, path, or ndarray; got {type(X)!r}"
        )

    def _check_fitted(self):
        if not self._is_fitted or self._network is None:
            raise RuntimeError(
                "Inverter not fitted.  Call fit() first."
            )
