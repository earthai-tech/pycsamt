"""
pycsamt.ai.training
===================

Framework-agnostic training utilities.
"""
from .augment import (
    AugmentFreqDrop,
    AugmentMixup,
    AugmentNoise,
    AugmentStaticShift,
    Compose,
    RandomApply,
)
from .dataset import EMDataset, Normalizer
from .metrics import (
    depth_rmse,
    layer_rmse,
    mae,
    masked_mse_loss,
    r2,
    relative_rmse,
    rmse,
    summarise,
)
from .trainer import EMTrainer

__all__ = [
    "Normalizer", "EMDataset",
    "EMTrainer",
    "rmse", "mae", "r2", "relative_rmse",
    "depth_rmse", "layer_rmse",
    "masked_mse_loss", "summarise",
    "AugmentNoise", "AugmentStaticShift", "AugmentFreqDrop",
    "AugmentMixup", "Compose", "RandomApply",
]
