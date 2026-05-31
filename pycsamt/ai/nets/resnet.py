# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
1-D Residual CNN for EM inversion — Liu et al. (2021) I8RFCN style.

Skip (residual) connections allow training very deep 1-D CNNs without
vanishing gradients, enabling the network to learn both coarse resistivity
contrasts (long-range dependencies) and fine-grained frequency patterns
(short-range) simultaneously.

Three residual stages progressively double the channel width while
halving the spatial resolution.  A global average-pool collapses the
spatial dimension before the final linear output head.

References
----------
Liu, W. et al. (2021). Deep learning AMT inversion using residual-based
deep convolutional neural network. *Journal of Geophysics and Engineering*,
18(6), 876–888.
"""
from __future__ import annotations

from typing import Sequence

__all__ = ["ResNet1DNet"]


def _require_torch():
    try:
        import torch
        import torch.nn as nn
        return torch, nn
    except ImportError:
        raise ImportError(
            "PyTorch is required for pycsamt.ai.nets.  "
            "Install with: pip install torch"
        )


class ResNet1DNet:
    """
    Factory wrapper for the 1-D residual CNN.

    Parameters
    ----------
    n_features : int
        Input feature-vector length.
    n_out : int
        Output parameter vector length (``2*n_layers - 1``).
    channels : sequence of int
        Filter counts for the three residual stages.
        Default ``(64, 128, 256)`` replicates the I8RFCN paper.
    dropout : float
        Dropout before the final linear layer.
    n_blocks : int
        Number of residual blocks per stage (default 2).
    """

    def __init__(
        self,
        n_features: int,
        n_out: int,
        *,
        channels: Sequence[int] = (64, 128, 256),
        dropout: float = 0.3,
        n_blocks: int = 2,
    ):
        self.n_features = n_features
        self.n_out = n_out
        self.channels = tuple(channels)
        self.dropout = dropout
        self.n_blocks = n_blocks

    def build(self):
        """Return the ``nn.Module``."""
        return _build_resnet1d(
            self.n_features, self.n_out,
            channels=self.channels,
            dropout=self.dropout,
            n_blocks=self.n_blocks,
        )

    def __repr__(self):
        return (
            f"ResNet1DNet(n_features={self.n_features}, n_out={self.n_out}, "
            f"channels={self.channels}, n_blocks={self.n_blocks})"
        )


# ─── Internal build ────────────────────────────────────────────────────────

def _build_resnet1d(n_features, n_out, channels, dropout, n_blocks):
    torch, nn = _require_torch()
    import torch.nn.functional as F

    class _ResBlock(nn.Module):
        """Basic 1-D residual block with optional downsampling."""
        def __init__(self, in_ch, out_ch, stride=1):
            super().__init__()
            self.conv1 = nn.Conv1d(
                in_ch, out_ch, 3, stride=stride, padding=1, bias=False
            )
            self.bn1 = nn.BatchNorm1d(out_ch)
            self.conv2 = nn.Conv1d(
                out_ch, out_ch, 3, padding=1, bias=False
            )
            self.bn2 = nn.BatchNorm1d(out_ch)
            # Projection shortcut when dimensions change
            if stride != 1 or in_ch != out_ch:
                self.shortcut = nn.Sequential(
                    nn.Conv1d(in_ch, out_ch, 1, stride=stride, bias=False),
                    nn.BatchNorm1d(out_ch),
                )
            else:
                self.shortcut = nn.Identity()

        def forward(self, x):
            out = F.relu(self.bn1(self.conv1(x)), inplace=True)
            out = self.bn2(self.conv2(out))
            return F.relu(out + self.shortcut(x), inplace=True)

    class _ResNet1DModule(nn.Module):
        def __init__(self):
            super().__init__()
            # Stem: initial wide conv
            self.stem = nn.Sequential(
                nn.Conv1d(1, channels[0], 7, padding=3, bias=False),
                nn.BatchNorm1d(channels[0]),
                nn.ReLU(inplace=True),
            )
            # Residual stages
            stage_list = []
            in_ch = channels[0]
            for stage_i, out_ch in enumerate(channels):
                for block_i in range(n_blocks):
                    stride = 2 if (stage_i > 0 and block_i == 0) else 1
                    stage_list.append(_ResBlock(in_ch, out_ch, stride))
                    in_ch = out_ch
            self.stages = nn.Sequential(*stage_list)
            self.pool = nn.AdaptiveAvgPool1d(1)   # global average pool
            self.dropout = nn.Dropout(dropout)
            self.fc = nn.Linear(channels[-1], n_out)

        def forward(self, x):
            # x: (batch, n_features)
            x = x.unsqueeze(1)              # (batch, 1, n_features)
            x = self.stem(x)
            x = self.stages(x)
            x = self.pool(x).squeeze(-1)    # (batch, channels[-1])
            x = self.dropout(x)
            return self.fc(x)

    return _ResNet1DModule()
