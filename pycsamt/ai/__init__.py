# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.ai
==========

Artificial intelligence and machine learning for electromagnetic data
processing and inversion.

This package introduces AI-based workflows as a core component of
pycsamt v2.  It is structured around five capability areas:

``nets/``
    Neural network architectures (CNN, ResNet, U-Net, DRCNN, FCN).

``inversion/``
    High-level 1-D and 2-D inversion workflows that consume
    :class:`~pycsamt.z.z.Z` objects or
    :class:`~pycsamt.forward.em1d.ForwardResponse` and return
    subsurface resistivity models.

``processing/``
    ML-based data preparation: denoising, automated QC scoring,
    anomaly detection, and dimensionality/strike classification.

``training/``
    Framework-agnostic training utilities: :class:`EMDataset`,
    training loop with early stopping, EM-aware augmentation, and
    evaluation metrics.

``plot/``
    Unified, publication-quality visualisation with a shared
    :class:`~pycsamt.ai.plot._style.EMStyle` context manager.

Dependencies
------------
The core abstract classes in this package (:mod:`pycsamt.ai._base`)
require only NumPy.  PyTorch (or TensorFlow) is imported **lazily**
inside concrete network classes and is only needed when a network is
actually instantiated or trained.

References
----------
Papers in ``related-works/ml/`` that informed the architecture design:

* Puzyrev & Swidinsky (2021) — 1D CNN inversion for CSEM and TEM
* Liu et al. (2021) — 18-layer ResNet for AMT inversion
* Moghadas (2020) — FCN for EMI inversion
* Oh et al. (2019, 2020) — CNN/U-Net for CSEM salt delineation
* Guo et al. (2021) — DRCNN for joint AMT+seismic inversion

Quick start (Phase 2+)
----------------------
>>> from pycsamt.ai.inversion import EMInverter1D
>>> inv = EMInverter1D(arch='resnet', n_layers=6, solver='mt1d')
>>> inv.fit('mt1d_train.npz', epochs=50)
>>> model = inv.predict(z_obj)
"""

from ._base import BaseEMNet, BaseEMProcessor, EMCheckpoint

__all__ = [
    "BaseEMNet",
    "BaseEMProcessor",
    "EMCheckpoint",
]
