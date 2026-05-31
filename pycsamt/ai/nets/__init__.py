# -*- coding: utf-8 -*-
"""
pycsamt.ai.nets
===============

Neural network architecture implementations.

All classes are framework-agnostic factory wrappers that return a
``torch.nn.Module`` on ``.build()``.  PyTorch is imported lazily
inside each ``build()`` call.
"""
from .cnn1d import CNN1DNet
from .resnet import ResNet1DNet
from .fcn import FCN1DNet

__all__ = ["CNN1DNet", "ResNet1DNet", "FCN1DNet"]
