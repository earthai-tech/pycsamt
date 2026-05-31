# -*- coding: utf-8 -*-
"""
pycsamt.ai.nets
===============

Neural network architecture implementations for EM inversion.

Planned modules (Phase 2)
--------------------------
cnn1d.py
    1-D CNN with encoder–decoder structure (Puzyrev 2019/2021 style).
    Suitable for MT/CSAMT/TEM 1-D inversion.

resnet.py
    Deep residual convolutional network — 18-layer I8RFCN
    (Liu et al. 2021 style).  Improved gradient flow for deep models.

fcn.py
    Fully convolutional network (Moghadas 2020 style).  No fully-
    connected head; works on variable-length input.

unet.py
    U-Net with skip connections for 2-D section inversion
    (Oh et al. 2020 style).  Used for image-to-image tasks.

drcnn.py
    Dual-path deep residual CNN for joint EM + seismic inversion
    (Guo et al. 2021 style).  Takes two input branches.

All architectures in this subpackage require PyTorch.
"""
__all__: list = []
