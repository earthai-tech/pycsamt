# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Artificial intelligence and machine learning for EM processing and
inversion.

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
from ._zoo import (
    download_checkpoint,
    get_pretrained_info,
    list_pretrained,
)
from .inversion import (
    ConformalPredictor,
    DUHIInverter2D,
    DUHIPreparation,
    EMInverter1D,
    EMInverter2D,
    EnsembleInverter,
    GCNInverter3D,
    JointInverter,
    PosteriorCalibrator,
    combine_observation_reliability,
    dimensionality_reliability,
)
from .nets import (
    CNN1DNet,
    DRCNNNet,
    FCN1DNet,
    GCNNet,
    ResNet1DNet,
    UNet2DNet,
    build_adjacency,
)
from .plot import (
    EM_CMAPS,
    EM_COLORS,
    EM_FIGSIZE,
    EMStyle,
    add_colorbar,
    em_context,
    plot_compare,
    plot_confusion_matrix,
    plot_convergence,
    plot_feature_importance,
    plot_layer_errors,
    plot_lr_schedule,
    plot_profile_pair,
    plot_pseudo_section,
    plot_residuals,
    plot_section,
    plot_section_pair,
    plot_uncertainty_bands,
)
from .processing import (
    AnomalyDetector,
    DimensionalityClassifier,
    EMDenoiser,
    EMQCScorer,
    prepare_z_features,
)
from .training import (
    AugmentFreqDrop,
    AugmentMixup,
    AugmentNoise,
    AugmentStaticShift,
    Compose,
    EMDataset,
    EMTrainer,
    Normalizer,
    RandomApply,
    depth_rmse,
    layer_rmse,
    mae,
    masked_mse_loss,
    r2,
    relative_rmse,
    rmse,
    summarise,
)

__all__ = [
    # base
    "BaseEMNet",
    "BaseEMProcessor",
    "EMCheckpoint",
    # inversion
    "EMInverter1D",
    "EMInverter2D",
    "GCNInverter3D",
    "JointInverter",
    "EnsembleInverter",
    "DUHIInverter2D",
    "DUHIPreparation",
    "dimensionality_reliability",
    "combine_observation_reliability",
    # calibrated UQ
    "ConformalPredictor",
    "PosteriorCalibrator",
    # model zoo
    "list_pretrained",
    "get_pretrained_info",
    "download_checkpoint",
    # nets
    "CNN1DNet",
    "ResNet1DNet",
    "FCN1DNet",
    "UNet2DNet",
    "DRCNNNet",
    "GCNNet",
    "build_adjacency",
    # training
    "Normalizer",
    "EMDataset",
    "EMTrainer",
    "rmse",
    "mae",
    "r2",
    "relative_rmse",
    "depth_rmse",
    "layer_rmse",
    "masked_mse_loss",
    "summarise",
    # augmentation (Phase 5)
    "AugmentNoise",
    "AugmentStaticShift",
    "AugmentFreqDrop",
    "AugmentMixup",
    "Compose",
    "RandomApply",
    # plot
    "EMStyle",
    "EM_COLORS",
    "EM_CMAPS",
    "EM_FIGSIZE",
    "em_context",
    "add_colorbar",
    "plot_compare",
    "plot_profile_pair",
    "plot_convergence",
    "plot_lr_schedule",
    # section & diagnostics (Phase 4)
    "plot_section",
    "plot_section_pair",
    "plot_pseudo_section",
    "plot_confusion_matrix",
    "plot_residuals",
    "plot_layer_errors",
    "plot_uncertainty_bands",
    "plot_feature_importance",
    # processing (Phase 3)
    "EMDenoiser",
    "prepare_z_features",
    "EMQCScorer",
    "AnomalyDetector",
    "DimensionalityClassifier",
]
