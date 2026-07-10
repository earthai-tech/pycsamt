# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
A convenient aggregator for all Quality-Control (QC)
components.

This module re-exports the individual QC metric classes from
the `var_pc` (percent-variation) and `var_std` (phase standard
deviation) modules into a single, easy-to-import namespace.
"""
from __future__ import annotations

from .var_pc import (
    EmagPctErr,
    HmagPctErr,
    PcEmag,
    PcHmag,
    PcRho,
    RhoPctErr,
)
from .var_std import PhaseSigma, SEphz, SHphz, SPhz

__all__ = [
    "PcEmag",
    "PcHmag",
    "PcRho",
    "EmagPctErr",
    "HmagPctErr",
    "RhoPctErr",
    "SEphz",
    "SHphz",
    "SPhz",
    "PhaseSigma",
]
