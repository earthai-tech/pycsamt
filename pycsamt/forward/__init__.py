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
from .batch import generate_dataset, ForwardDataset, generate_dataset_3d, SurveyDataset3D
from .config import ForwardConfig
from .grid2d import Grid2D, make_padding
from .em2d import MT2DForward, ForwardResponse2D
from .config2d import ForwardConfig2D
from .plot import (
    plot_response_1d,
    plot_model_1d,
    plot_response_and_model_1d,
    plot_model_2d,
    plot_pseudosection_2d,
    plot_response_profiles,
)

__all__ = [
    # solvers 1-D
    "MT1DForward",
    "TEM1DForward",
    "CSAMT1DForward",
    "ForwardResponse",
    # solvers 2-D
    "MT2DForward",
    "ForwardResponse2D",
    # models
    "LayeredModel",
    "GEOLOGY_PRIORS",
    "Grid2D",
    "make_padding",
    # noise
    "GaussianNoise",
    "FieldRealisticNoise",
    "MultiplicativeNoise",
    "add_gaussian_noise",
    "add_noise",
    # batch 1-D
    "generate_dataset",
    "ForwardDataset",
    # batch 3-D
    "generate_dataset_3d",
    "SurveyDataset3D",
    # configuration
    "ForwardConfig",
    "ForwardConfig2D",
    # plotting
    "plot_response_1d",
    "plot_model_1d",
    "plot_response_and_model_1d",
    "plot_model_2d",
    "plot_pseudosection_2d",
    "plot_response_profiles",
]
