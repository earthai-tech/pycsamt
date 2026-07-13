# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
DRCNNNet — Dense Residual CNN for joint / multi-modal EM inversion.

Adapts the Dense-Residual Convolutional Neural Network architecture of
Guo et al. (2021) to the general case of fusing two or more 1-D EM
datasets (e.g., AMT + TEM, MT + gravity) into a single subsurface
model prediction.

Architecture
------------
Each modality is encoded by an independent 1-D dense block.  Dense
blocks use DenseNet-style growth: every sub-layer receives the
concatenation of all previous sub-layer outputs plus the original
input.  A residual shortcut connects the block input to its output so
the decoder sees both a bottleneck and the original features.

After encoding, features from all modalities are concatenated and
processed by a shared fusion dense block, followed by a linear output
head.

.. math::

    \\mathbf{h}_k^{(l)} = \\sigma\\bigl(
        W_k^{(l)} [\\mathbf{x}_k; \\mathbf{h}_k^{(1)};
                    \\ldots;
                    \\mathbf{h}_k^{(l-1)}]
    \\bigr)

where :math:`k` indexes the modality and :math:`l` indexes the
dense-block sub-layer.

References
----------
Guo R. et al. (2021) *IEEE TGRS* — DRCNN for joint AMT+seismic.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

__all__ = ["DRCNNNet"]


# ── Module-level class definitions ───────────────
# At module scope so pickle can locate them when
# coordinator.py checkpoints AgentResult dicts.

try:
    import torch
    import torch.nn as nn

    class _DenseBlock(nn.Module):
        def __init__(
            self,
            in_features: int,
            out_features: int,
            growth_rate: int,
            n_layers: int,
            dropout: float,
        ) -> None:
            super().__init__()
            dims: list[int] = [in_features]
            self.sub_layers = nn.ModuleList()
            for _ in range(n_layers):
                d_in = sum(dims)
                self.sub_layers.append(
                    nn.Sequential(
                        nn.Linear(d_in, growth_rate),
                        nn.BatchNorm1d(growth_rate),
                        nn.ReLU(),
                        nn.Dropout(dropout),
                    )
                )
                dims.append(growth_rate)
            total_dim = sum(dims)
            self.transition = nn.Sequential(
                nn.Linear(total_dim, out_features),
                nn.BatchNorm1d(out_features),
                nn.ReLU(),
            )
            if in_features != out_features:
                self.shortcut = nn.Linear(
                    in_features,
                    out_features,
                    bias=False,
                )
            else:
                self.shortcut = nn.Identity()

        def forward(self, x):
            outputs = [x]
            for layer in self.sub_layers:
                cat = torch.cat(outputs, dim=-1)
                h = layer(cat)
                outputs.append(h)
            all_cat = torch.cat(outputs, dim=-1)
            out = self.transition(all_cat)
            return out + self.shortcut(x)

    class _DRCNN(nn.Module):
        def __init__(
            self,
            n_features_list: tuple[int, ...],
            n_out: int,
            growth_rate: int,
            n_layers: int,
            hidden_dim: int,
            dropout: float,
        ) -> None:
            super().__init__()
            self.encoders = nn.ModuleList(
                [
                    _DenseBlock(
                        nf,
                        hidden_dim,
                        growth_rate,
                        n_layers,
                        dropout,
                    )
                    for nf in n_features_list
                ]
            )
            fused_dim = hidden_dim * len(n_features_list)
            self.fusion = _DenseBlock(
                fused_dim,
                hidden_dim,
                growth_rate,
                n_layers,
                dropout,
            )
            self.head = nn.Sequential(
                nn.Linear(hidden_dim, hidden_dim // 2),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_dim // 2, n_out),
            )

        def forward(self, *inputs):
            if len(inputs) != len(self.encoders):
                raise ValueError(
                    f"Expected {len(self.encoders)} inputs, got {len(inputs)}"
                )
            encoded = [enc(x) for enc, x in zip(self.encoders, inputs)]
            fused = torch.cat(encoded, dim=-1)
            out = self.fusion(fused)
            return self.head(out)

except ImportError:
    pass  # torch not installed


# ── Factory class ─────────────────────────────────


class DRCNNNet:
    r"""
    Factory wrapper for the Dense Residual CNN.

    Parameters
    ----------
    n_features_list : sequence of int
        Feature vector length for each input modality.
        E.g. ``(120, 48)`` for 120 MT features and
        48 seismic features.
    n_out : int
        Output dimension (number of subsurface
        parameters to predict).
    growth_rate : int, default 32
        New channels added by each sub-layer in a
        dense block.
    n_layers : int, default 6
        Number of sub-layers per dense block.
    hidden_dim : int, default 256
        Dimension of the encoded representation from
        each modality and the fusion block output.
    dropout : float, default 0.2

    Examples
    --------
    >>> # MT (120 features) + TEM (48 features)
    >>> drcnn = DRCNNNet((120, 48), n_out=9)
    >>> model = drcnn.build()   # doctest: +SKIP
    """

    def __init__(
        self,
        n_features_list: Sequence[int],
        n_out: int,
        *,
        growth_rate: int = 32,
        n_layers: int = 6,
        hidden_dim: int = 256,
        dropout: float = 0.2,
    ) -> None:
        self.n_features_list = tuple(int(n) for n in n_features_list)
        self.n_out = int(n_out)
        self.growth_rate = int(growth_rate)
        self.n_layers = int(n_layers)
        self.hidden_dim = int(hidden_dim)
        self.dropout = float(dropout)

    def build(self) -> Any:
        """Instantiate and return the ``nn.Module``."""
        return _build_drcnn(
            self.n_features_list,
            self.n_out,
            self.growth_rate,
            self.n_layers,
            self.hidden_dim,
            self.dropout,
        )

    def __repr__(self) -> str:
        return (
            f"DRCNNNet("
            f"n_features_list={self.n_features_list}, "
            f"n_out={self.n_out}, "
            f"hidden_dim={self.hidden_dim})"
        )


# ── Internal build helpers ────────────────────────


def _dense_block_1d(
    in_features: int,
    out_features: int,
    growth_rate: int,
    n_layers: int,
    dropout: float,
) -> Any:
    """Return a 1-D dense block as an ``nn.Module``."""
    if "_DenseBlock" not in globals():
        raise ImportError(
            "PyTorch is required for DRCNNNet. "
            "Install with: pip install torch"
        )
    return _DenseBlock(
        in_features,
        out_features,
        growth_rate,
        n_layers,
        dropout,
    )


def _build_drcnn(
    n_features_list: tuple[int, ...],
    n_out: int,
    growth_rate: int,
    n_layers: int,
    hidden_dim: int,
    dropout: float,
) -> Any:
    """Build the full multi-modal DRCNN."""
    if "_DRCNN" not in globals():
        raise ImportError(
            "PyTorch is required for DRCNNNet. "
            "Install with: pip install torch"
        )
    return _DRCNN(
        tuple(n_features_list),
        n_out,
        growth_rate,
        n_layers,
        hidden_dim,
        dropout,
    )
