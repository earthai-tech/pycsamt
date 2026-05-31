# -*- coding: utf-8 -*-
"""
pycsamt.ai.training
===================

Framework-agnostic training utilities.
"""
from .dataset import Normalizer, EMDataset
from .trainer import EMTrainer
from .metrics import (
    rmse, mae, r2, relative_rmse,
    depth_rmse, layer_rmse,
    masked_mse_loss, summarise,
)

__all__ = [
    "Normalizer", "EMDataset",
    "EMTrainer",
    "rmse", "mae", "r2", "relative_rmse",
    "depth_rmse", "layer_rmse",
    "masked_mse_loss", "summarise",
]
