# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.forward
===============

Physics-based 1-D electromagnetic forward solvers and synthetic dataset
generation for AI/ML training.

This package has **no machine-learning dependencies** — only NumPy and
SciPy are required.  It provides the building blocks consumed by
:mod:`pycsamt.ai` for training data generation.

Submodules
----------
em1d
    Three 1-D solvers: :class:`~pycsamt.forward.em1d.MT1DForward`,
    :class:`~pycsamt.forward.em1d.TEM1DForward`,
    :class:`~pycsamt.forward.em1d.CSAMT1DForward`.

synthetic
    :class:`~pycsamt.forward.synthetic.LayeredModel` — the central data
    object for all 1-D models, with ``random``, ``blocky``, ``smooth``,
    and ``from_geology`` constructors.

noise
    Noise models for synthetic data:
    :class:`~pycsamt.forward.noise.GaussianNoise`,
    :class:`~pycsamt.forward.noise.FieldRealisticNoise`,
    :class:`~pycsamt.forward.noise.MultiplicativeNoise`.

batch
    :func:`~pycsamt.forward.batch.generate_dataset` — parallelised
    (model, response) pair generation for large training sets.

Quick start
-----------
>>> import numpy as np
>>> from pycsamt.forward import MT1DForward, LayeredModel
>>> model = LayeredModel([100, 10, 500], [300, 800])
>>> resp = MT1DForward(np.logspace(-3, 4, 30)).run(model)
>>> resp.rho_a.shape
(30,)
"""

from .em1d import (
    MT1DForward,
    TEM1DForward,
    CSAMT1DForward,
    ForwardResponse,
)
from .synthetic import LayeredModel, GEOLOGY_PRIORS
from .noise import (
    GaussianNoise,
    FieldRealisticNoise,
    MultiplicativeNoise,
    add_gaussian_noise,
    add_noise,
)
from .batch import generate_dataset, ForwardDataset

__all__ = [
    # solvers
    "MT1DForward",
    "TEM1DForward",
    "CSAMT1DForward",
    "ForwardResponse",
    # models
    "LayeredModel",
    "GEOLOGY_PRIORS",
    # noise
    "GaussianNoise",
    "FieldRealisticNoise",
    "MultiplicativeNoise",
    "add_gaussian_noise",
    "add_noise",
    # batch
    "generate_dataset",
    "ForwardDataset",
]
