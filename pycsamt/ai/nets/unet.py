# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
UNet2DNet — 2-D U-Net architecture for EM section inversion.

Based on the encoder-decoder with skip connections introduced by
Ronneberger et al. (2015) and adapted for 2-D CSEM / MT inversion by
Oh et al. (2019, 2020).

The network maps a 2-D observed data panel

.. math::

    \\mathbf{D} \\in \\mathbb{R}^{C_{\\text{in}} \\times N_f \\times N_s}

(components × frequencies × stations) to a 2-D subsurface resistivity
section

.. math::

    \\hat{\\boldsymbol{\\rho}} \\in \\mathbb{R}^{1 \\times N_z \\times N_s}

(depth × stations).  Both spatial dimensions may differ between input
and output; bilinear upsampling handles the mismatch.

References
----------
Oh S. et al. (2019) *Geophysics* — 2D CSEM inversion with a U-Net.
Oh S. et al. (2020) *JGR Solid Earth* — generalisation to salt imaging.
"""
from __future__ import annotations

from typing import Any, Dict, Optional, Tuple

__all__ = ["UNet2DNet"]


# ─────────────────────────────────────────────────────────────────────────────
# Internal network builder (lazy torch import)
# ─────────────────────────────────────────────────────────────────────────────

def _conv_block(in_ch: int, out_ch: int, dropout: float) -> Any:
    import torch.nn as nn
    return nn.Sequential(
        nn.Conv2d(in_ch, out_ch, 3, padding=1, bias=False),
        nn.BatchNorm2d(out_ch),
        nn.ReLU(inplace=True),
        nn.Conv2d(out_ch, out_ch, 3, padding=1, bias=False),
        nn.BatchNorm2d(out_ch),
        nn.ReLU(inplace=True),
        nn.Dropout2d(dropout),
    )


def _build_unet2d(
    n_in: int,
    n_out: int,
    channels: Tuple[int, ...],
    dropout: float,
) -> Any:
    """Return an ``nn.Module`` implementing 2-D U-Net."""
    try:
        import torch
        import torch.nn as nn
        import torch.nn.functional as F
    except ImportError as exc:
        raise ImportError("PyTorch is required for UNet2DNet") from exc

    ch = list(channels)
    n_stages = len(ch) - 1  # last entry is the bridge width

    class _UNet2D(nn.Module):
        def __init__(self) -> None:
            super().__init__()
            # Encoder stages
            in_chs = [n_in] + ch[:-1]
            self.encoders = nn.ModuleList(
                [_conv_block(in_chs[i], ch[i], dropout) for i in range(n_stages)]
            )
            self.pools = nn.ModuleList(
                [nn.MaxPool2d(2) for _ in range(n_stages)]
            )
            # Bridge
            self.bridge = _conv_block(ch[n_stages - 1], ch[-1], dropout)
            # Decoder stages (reverse)
            dec_in = [ch[-1] + ch[i] for i in range(n_stages - 1, -1, -1)]
            dec_out = [ch[i] for i in range(n_stages - 1, -1, -1)]
            self.decoders = nn.ModuleList(
                [_conv_block(dec_in[i], dec_out[i], dropout)
                 for i in range(n_stages)]
            )
            self.out_conv = nn.Conv2d(ch[0], n_out, 1)

        def forward(self, x):
            skips = []
            for enc, pool in zip(self.encoders, self.pools):
                x = enc(x)
                skips.append(x)
                x = pool(x)

            x = self.bridge(x)

            for dec, skip in zip(self.decoders, reversed(skips)):
                x = F.interpolate(
                    x, size=skip.shape[-2:], mode="bilinear", align_corners=False
                )
                x = torch.cat([x, skip], dim=1)
                x = dec(x)

            return self.out_conv(x)

    return _UNet2D()


# ─────────────────────────────────────────────────────────────────────────────
# UNet2DNet factory wrapper
# ─────────────────────────────────────────────────────────────────────────────

class UNet2DNet:
    """
    Factory wrapper for the 2-D U-Net.

    All heavy imports are deferred to :meth:`build`.

    Parameters
    ----------
    n_in : int
        Input channels — typically ``n_components`` (e.g. 4 for
        off-diagonal MT impedance: ``log|Zxy|``, ``φ_xy``,
        ``log|Zyx|``, ``φ_yx``).
    n_out : int
        Output channels — 1 for a single log₁₀(ρ) section.
    channels : tuple of int, default (32, 64, 128, 256, 512)
        Channel widths at each encoder stage plus the bridge.
        ``len(channels) - 1`` determines the number of
        pooling/upsampling stages.
    dropout : float, default 0.2
        2-D spatial dropout probability in each convolutional block.

    Examples
    --------
    >>> net_factory = UNet2DNet(n_in=4, n_out=1)
    >>> model = net_factory.build()   # returns nn.Module  # doctest: +SKIP
    """

    def __init__(
        self,
        n_in: int,
        n_out: int = 1,
        *,
        channels: Tuple[int, ...] = (32, 64, 128, 256, 512),
        dropout: float = 0.2,
    ) -> None:
        self.n_in = int(n_in)
        self.n_out = int(n_out)
        self.channels = tuple(channels)
        self.dropout = float(dropout)

    def build(self) -> Any:
        """Instantiate and return the ``torch.nn.Module``."""
        return _build_unet2d(self.n_in, self.n_out, self.channels, self.dropout)

    def __repr__(self) -> str:
        return (
            f"UNet2DNet(n_in={self.n_in}, n_out={self.n_out}, "
            f"channels={self.channels})"
        )
