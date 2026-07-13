# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
pycsamt.zonge.components
------------------------

Top-level module that aggregates all Zonge data components into a
single namespace for convenient importing.

This module re-exports all primary component classes for handling
various aspects of Zonge AVG data, including measurements, survey
geometry, core scientific values (resistivity, phase, impedance),
and quality-control metrics.
"""

from __future__ import annotations

from .meas import Amps, CompMeas, Frequency
from .qc import (
    PcEmag,
    PcHmag,
    PcRho,
    SEphz,
    SHphz,
    SPhz,
)
from .resphase import Phase, Resistivity
from .survey import Station
from .z import Z

__all__ = [
    # Core measurements
    "Frequency",
    "Amps",
    "CompMeas",
    # Survey geometry
    "Station",
    # Scientific outputs
    "Resistivity",
    "Phase",
    "Z",
    # Quality Control (QC) metrics
    "PcEmag",
    "PcHmag",
    "PcRho",
    "SEphz",
    "SHphz",
    "SPhz",
]
