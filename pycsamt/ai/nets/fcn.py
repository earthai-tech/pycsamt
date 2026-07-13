# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Fully Convolutional Network for EM inversion — Moghadas (2020) style.

Contains no fully-connected layers.  All weights are convolutional,
including the final output projection (1x1 convolution).  This makes
the network invariant to input length, which is useful when different
survey setups use different numbers of frequencies.

The global average pool collapses the spatial dimension to a single
vector, and a final 1x1 convolution maps to the output size.

References
----------
Moghadas, D. (2020). One-dimensional deep learning inversion of EM
induction data using CNN. *Geophysical Journal International*,
223(1), 627-641.
"""

from __future__ import annotations

from collections.abc import Sequence

__all__ = ["FCN1DNet"]


def _require_torch():
    try:
        import torch
        import torch.nn as nn

        return torch, nn
    except ImportError:
        raise ImportError(
            "PyTorch is required for pycsamt.ai.nets. "
            "Install with: pip install torch"
        )


# ── Module-level class definition ────────────────
# At module scope so pickle can locate it when
# coordinator.py checkpoints AgentResult dicts.

try:
    import torch.nn as nn

    class _FCN1DModule(nn.Module):
        def __init__(self, channels, dropout, n_out):
            super().__init__()
            encoder = []
            in_ch = 1
            for out_ch in channels[:-1]:
                encoder += [
                    nn.Conv1d(in_ch, out_ch, 3, padding=1),
                    nn.BatchNorm1d(out_ch),
                    nn.ReLU(inplace=True),
                    nn.MaxPool1d(2),
                ]
                in_ch = out_ch
            self.encoder = nn.Sequential(*encoder)
            self.bottleneck = nn.Sequential(
                nn.Conv1d(in_ch, channels[-1], 1),
                nn.BatchNorm1d(channels[-1]),
                nn.ReLU(inplace=True),
                nn.Dropout(dropout),
            )
            self.gap = nn.AdaptiveAvgPool1d(1)
            self.out_proj = nn.Conv1d(channels[-1], n_out, 1)

        def forward(self, x):
            x = x.unsqueeze(1)
            x = self.encoder(x)
            x = self.bottleneck(x)
            x = self.gap(x)
            x = self.out_proj(x).squeeze(-1)
            return x

except ImportError:
    pass  # torch not installed


# ── Factory class ─────────────────────────────────


class FCN1DNet:
    r"""
    Factory wrapper for the fully-convolutional
    1-D EM network.

    Parameters
    ----------
    n_features : int
        Input feature vector length.
    n_out : int
        Output parameter vector length
        (``2*n_layers - 1``).
    channels : sequence of int
        Channel widths for each convolutional block.
        The last entry is used for the bottleneck
        and global pool.
    dropout : float
        Spatial dropout rate.
    """

    def __init__(
        self,
        n_features: int,
        n_out: int,
        *,
        channels: Sequence[int] = (32, 64, 128, 64),
        dropout: float = 0.2,
    ):
        self.n_features = n_features
        self.n_out = n_out
        self.channels = tuple(channels)
        self.dropout = dropout

    def build(self):
        """Return the ``nn.Module``."""
        return _build_fcn1d(
            self.n_features,
            self.n_out,
            channels=self.channels,
            dropout=self.dropout,
        )

    def __repr__(self):
        return (
            f"FCN1DNet("
            f"n_features={self.n_features}, "
            f"n_out={self.n_out}, "
            f"channels={self.channels})"
        )


# ── Internal build ────────────────────────────────


def _build_fcn1d(n_features, n_out, channels, dropout):
    _require_torch()
    return _FCN1DModule(tuple(channels), dropout, n_out)
