# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.ai
==========

Artificial intelligence and machine learning for EM data processing
and inversion.

Phase 2 additions
-----------------
* :class:`~pycsamt.ai.inversion.inv1d.EMInverter1D` — full 1-D inversion
  pipeline (data loading → normalisation → training → prediction)
* :mod:`~pycsamt.ai.nets` — CNN1D, ResNet1D, FCN1D architectures
* :mod:`~pycsamt.ai.training` — EMDataset, EMTrainer, metrics
* :mod:`~pycsamt.ai.plot` — EMStyle, plot_compare, plot_convergence

Quick start (requires PyTorch)
------------------------------
>>> from pycsamt.forward.batch import generate_dataset
>>> from pycsamt.ai.inversion import EMInverter1D
>>> ds = generate_dataset(n_samples=2_000, seed=0, n_layers=5)
>>> inv = EMInverter1D(arch="resnet", n_layers=5)
>>> inv.fit(ds, epochs=30)
>>> # Predict on new Z objects or ForwardResponse
>>> y_pred = inv.predict(X_test)
"""

from ._base import BaseEMNet, BaseEMProcessor, EMCheckpoint
from .inversion import EMInverter1D, EMInverter2D, JointInverter
from .nets import CNN1DNet, ResNet1DNet, FCN1DNet, UNet2DNet, DRCNNNet
from .training import (
    Normalizer, EMDataset, EMTrainer,
    rmse, mae, r2, relative_rmse,
    depth_rmse, layer_rmse, masked_mse_loss, summarise,
)
from .plot import (
    EMStyle, EM_COLORS, EM_CMAPS, EM_FIGSIZE,
    em_context, add_colorbar,
    plot_compare, plot_profile_pair,
    plot_convergence, plot_lr_schedule,
    plot_section, plot_section_pair, plot_pseudo_section,
    plot_confusion_matrix, plot_residuals, plot_layer_errors,
    plot_uncertainty_bands, plot_feature_importance,
)
from .processing import (
    EMDenoiser, prepare_z_features,
    EMQCScorer,
    AnomalyDetector,
    DimensionalityClassifier,
)

__all__ = [
    # base
    "BaseEMNet", "BaseEMProcessor", "EMCheckpoint",
    # inversion
    "EMInverter1D", "EMInverter2D", "JointInverter",
    # nets
    "CNN1DNet", "ResNet1DNet", "FCN1DNet", "UNet2DNet", "DRCNNNet",
    # training
    "Normalizer", "EMDataset", "EMTrainer",
    "rmse", "mae", "r2", "relative_rmse",
    "depth_rmse", "layer_rmse", "masked_mse_loss", "summarise",
    # plot
    "EMStyle", "EM_COLORS", "EM_CMAPS", "EM_FIGSIZE",
    "em_context", "add_colorbar",
    "plot_compare", "plot_profile_pair",
    "plot_convergence", "plot_lr_schedule",
    # section & diagnostics (Phase 4)
    "plot_section", "plot_section_pair", "plot_pseudo_section",
    "plot_confusion_matrix", "plot_residuals", "plot_layer_errors",
    "plot_uncertainty_bands", "plot_feature_importance",
    # processing (Phase 3)
    "EMDenoiser", "prepare_z_features",
    "EMQCScorer",
    "AnomalyDetector",
    "DimensionalityClassifier",
]
