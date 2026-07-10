# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
UNet2DNet — 2-D U-Net architecture for EM section inversion.

Based on the encoder-decoder with skip connections introduced by
Ronneberger et al. (2015) and adapted for 2-D CSEM / MT inversion
by Oh et al. (2019, 2020).

The network maps a 2-D observed data panel

.. math::

    \\mathbf{D} \\in
    \\mathbb{R}^{C_{\\text{in}} \\times N_f \\times N_s}

(components x frequencies x stations) to a 2-D subsurface
resistivity section

.. math::

    \\hat{\\boldsymbol{\\rho}} \\in
    \\mathbb{R}^{1 \\times N_z \\times N_s}

(depth x stations).  Both spatial dimensions may differ between
input and output; bilinear upsampling handles the mismatch.

References
----------
Oh S. et al. (2019) *Geophysics* — 2D CSEM inversion with U-Net.
Oh S. et al. (2020) *JGR Solid Earth* — generalisation to salt.
"""
from __future__ import annotations

from typing import Any

__all__ = ["UNet2DNet"]


# ── Module-level class definitions ───────────────
# At module scope so pickle can locate them when
# coordinator.py checkpoints AgentResult dicts.

try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F

    def _conv_block(
        in_ch: int, out_ch: int, dropout: float
    ) -> nn.Sequential:
        return nn.Sequential(
            nn.Conv2d(
                in_ch, out_ch, 3,
                padding=1, bias=False,
            ),
            nn.BatchNorm2d(out_ch),
            nn.ReLU(inplace=True),
            nn.Conv2d(
                out_ch, out_ch, 3,
                padding=1, bias=False,
            ),
            nn.BatchNorm2d(out_ch),
            nn.ReLU(inplace=True),
            nn.Dropout2d(dropout),
        )

    class _UNet2D(nn.Module):
        def __init__(
            self,
            n_in: int,
            n_out: int,
            channels: tuple,
            dropout: float,
        ) -> None:
            super().__init__()
            ch = list(channels)
            n_stages = len(ch) - 1
            in_chs = [n_in] + ch[:-1]
            self.encoders = nn.ModuleList(
                [
                    _conv_block(
                        in_chs[i], ch[i], dropout
                    )
                    for i in range(n_stages)
                ]
            )
            self.pools = nn.ModuleList(
                [
                    nn.MaxPool2d(2)
                    for _ in range(n_stages)
                ]
            )
            self.bridge = _conv_block(
                ch[n_stages - 1], ch[-1], dropout
            )
            dec_out = [
                ch[i]
                for i in range(n_stages - 1, -1, -1)
            ]
            prev_ch = [ch[-1]] + dec_out[:-1]
            skip_ch = dec_out
            dec_in = [
                prev_ch[j] + skip_ch[j]
                for j in range(n_stages)
            ]
            self.decoders = nn.ModuleList(
                [
                    _conv_block(
                        dec_in[i], dec_out[i], dropout
                    )
                    for i in range(n_stages)
                ]
            )
            self.out_conv = nn.Conv2d(ch[0], n_out, 1)

        def forward(self, x):
            skips = []
            for enc, pool in zip(
                self.encoders, self.pools
            ):
                x = enc(x)
                skips.append(x)
                x = pool(x)
            x = self.bridge(x)
            for dec, skip in zip(
                self.decoders, reversed(skips)
            ):
                x = F.interpolate(
                    x,
                    size=skip.shape[-2:],
                    mode="bilinear",
                    align_corners=False,
                )
                x = torch.cat([x, skip], dim=1)
                x = dec(x)
            return self.out_conv(x)

except ImportError:
    pass  # torch not installed


# ── Factory class ─────────────────────────────────

class UNet2DNet:
    r"""
    Factory wrapper for the 2-D U-Net.

    All heavy imports are deferred to :meth:`build`.

    Parameters
    ----------
    n_in : int
        Input channels — typically ``n_components``
        (e.g. 4 for off-diagonal MT impedance:
        ``log|Zxy|``, ``phi_xy``,
        ``log|Zyx|``, ``phi_yx``).
    n_out : int
        Output channels — 1 for a single
        :math:`\\log_{10}(\\rho)` section.
    channels : tuple of int, default
        ``(32, 64, 128, 256, 512)``
        Channel widths at each encoder stage plus
        the bridge. ``len(channels) - 1`` determines
        the number of pooling/upsampling stages.
    dropout : float, default 0.2
        2-D spatial dropout probability in each
        convolutional block.

    Examples
    --------
    >>> net = UNet2DNet(n_in=4, n_out=1)
    >>> model = net.build()   # doctest: +SKIP
    """

    def __init__(
        self,
        n_in: int,
        n_out: int = 1,
        *,
        channels: tuple[int, ...] = (
            32, 64, 128, 256, 512
        ),
        dropout: float = 0.2,
    ) -> None:
        self.n_in = int(n_in)
        self.n_out = int(n_out)
        self.channels = tuple(channels)
        self.dropout = float(dropout)

    def build(self) -> Any:
        """Instantiate and return the ``nn.Module``."""
        return _build_unet2d(
            self.n_in,
            self.n_out,
            self.channels,
            self.dropout,
        )

    def __repr__(self) -> str:
        return (
            f"UNet2DNet("
            f"n_in={self.n_in}, "
            f"n_out={self.n_out}, "
            f"channels={self.channels})"
        )


# ── Internal build ────────────────────────────────

def _build_unet2d(
    n_in: int,
    n_out: int,
    channels: tuple[int, ...],
    dropout: float,
) -> Any:
    """Return an ``nn.Module`` implementing 2-D U-Net."""
    if "_UNet2D" not in globals():
        raise ImportError(
            "PyTorch is required for UNet2DNet. "
            "Install with: pip install torch"
        )
    return _UNet2D(n_in, n_out, tuple(channels), dropout)
