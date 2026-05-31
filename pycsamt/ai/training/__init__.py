# -*- coding: utf-8 -*-
"""
pycsamt.ai.training
===================

Framework-agnostic training utilities for EM neural networks.

Planned modules (Phase 2)
--------------------------
dataset.py
    :class:`EMDataset` — PyTorch / TensorFlow Dataset wrapper around
    :class:`~pycsamt.forward.batch.ForwardDataset`.  Handles
    normalisation, on-the-fly augmentation, and balanced sampling.

trainer.py
    :class:`EMTrainer` — training loop with early stopping,
    learning-rate scheduling, gradient clipping, and TensorBoard
    logging.  Works with any :class:`~pycsamt.ai._base.BaseEMNet`.

augment.py
    EM-aware data augmentation strategies:
    * Frequency permutation (shuffle a subset of frequency channels)
    * Noise injection (calls :mod:`pycsamt.forward.noise`)
    * Resistivity rescaling (multiply all resistivities by a constant)
    * Phase jitter (small random perturbation on phase values)

metrics.py
    Evaluation metrics for EM inversion:
    RMSE on log₁₀(ρ), relative RMSE, R², depth-weighted RMSE,
    IoU for binary salt/anomaly delineation, and calibration curves
    for ensemble uncertainty.
"""
__all__: list = []
