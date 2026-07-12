# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
1-D CNN for EM inversion — Puzyrev (2019/2021) architecture.

The input feature vector (log-scaled apparent resistivity + phase, or
log-scaled TEM dBz/dt) is treated as a 1-D sequence over
frequency/time channels.  Three convolutional blocks encode local
frequency dependencies; a small fully-connected head maps to the
output model parameter vector.

References
----------
Puzyrev, V. et al. (2019). Deep CNNs for 1D inversion of EM data.
*EAGE Conference 2019*.

Puzyrev, V. & Swidinsky, A. (2021). Inversion of 1D frequency- and
time-domain EM data with CNNs. *Computers & Geosciences*, 149, 104681.
"""
from __future__ import annotations

from collections.abc import Sequence

__all__ = ["CNN1DNet"]


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


class CNN1DNet:
    """
    Factory wrapper — call :func:`build` to get an ``nn.Module``.

    Use :class:`~pycsamt.ai.inversion.inv1d.EMInverter1D` instead of
    instantiating this class directly.

    Parameters
    ----------
    n_features : int
        Length of the input feature vector.
    n_out : int
        Length of the output parameter vector (``2*n_layers - 1``).
    channels : sequence of int
        Number of filters in each convolutional block.
    kernel_size : int
        Convolutional kernel width (same padding applied).
    dropout : float
        Dropout probability in the FC head.
    """

    def __init__(
        self,
        n_features: int,
        n_out: int,
        *,
        channels: Sequence[int] = (32, 64, 128),
        kernel_size: int = 5,
        dropout: float = 0.3,
    ):
        self.n_features = n_features
        self.n_out = n_out
        self.channels = tuple(channels)
        self.kernel_size = kernel_size
        self.dropout = dropout

    def build(self):
        """Return the ``nn.Module``."""
        torch, nn = _require_torch()
        return _CNN1DModule(
            self.n_features, self.n_out,
            channels=self.channels,
            kernel_size=self.kernel_size,
            dropout=self.dropout,
        )

    def __repr__(self):
        return (
            f"CNN1DNet(n_features={self.n_features}, n_out={self.n_out}, "
            f"channels={self.channels})"
        )


# ─── Internal nn.Module ───────────────────────────────────────────────────────

def _build_cnn1d(n_features, n_out, channels, kernel_size, dropout):
    """Build and return a CNN1D nn.Module (called lazily)."""
    torch, nn = _require_torch()

    class _CNN1DModule(nn.Module):
        def __init__(self):
            super().__init__()
            pad = kernel_size // 2
            blocks = []
            in_ch = 1
            length = int(n_features)
            for out_ch in channels:
                blocks += [
                    nn.Conv1d(in_ch, out_ch, kernel_size, padding=pad),
                    nn.BatchNorm1d(out_ch),
                    nn.ReLU(inplace=True),
                ]
                # pooling a length-1 sequence would produce size 0;
                # skip it once the sequence cannot be halved
                if length >= 2:
                    blocks.append(nn.MaxPool1d(2))
                    length //= 2
                in_ch = out_ch
            self.encoder = nn.Sequential(*blocks)

            # Compute flattened size after encoder. The probe must
            # run in eval mode: BatchNorm cannot normalise a
            # single-sample batch in training mode.
            with torch.no_grad():
                self.encoder.eval()
                dummy = torch.zeros(1, 1, n_features)
                flat = self.encoder(dummy).view(1, -1).shape[1]
                self.encoder.train()

            self.head = nn.Sequential(
                nn.Linear(flat, 256),
                nn.ReLU(inplace=True),
                nn.Dropout(dropout),
                nn.Linear(256, 128),
                nn.ReLU(inplace=True),
                nn.Linear(128, n_out),
            )

        def forward(self, x):
            # x: (batch, n_features)
            x = x.unsqueeze(1)               # (batch, 1, n_features)
            x = self.encoder(x)
            x = x.reshape(x.size(0), -1)
            return self.head(x)

    return _CNN1DModule()


# Monkey-patch build to use the lazy factory
def _cnn1d_build(self):
    return _build_cnn1d(
        self.n_features, self.n_out,
        self.channels, self.kernel_size, self.dropout
    )

CNN1DNet.build = _cnn1d_build
