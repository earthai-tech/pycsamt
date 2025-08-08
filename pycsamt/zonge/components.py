# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
One-stop public façade that re-exports all high-level containers
(metadata, measurements, QC helpers, tensor stacks…) so clients can
simply write ::

    from pycsamt.zonge.components import Z, Phase, Station, …

without digging into the internal sub-packages.
"""

from .heads import Head
from .properties import (
    Hardware,
    SurveyAnnotation,
    SurveyConfiguration,
    Receiver,
    Transmitter,
    SkipFlag,
)
from .meas import (
    CompMeas,
    Amps,
    Frequency,
)
from .resphase import Resistivity, Phase
from .survey   import Station
from .var      import (
    PcEmag,
    PcHmag,
    PcRho,
    SPhz,
    SHphz,
    SEphz,
)
from .z import Z
from .info import DataInfo 

__all__ = [
    # header & metadata
    "Head",
    "Hardware",
    "SurveyAnnotation",
    "SurveyConfiguration",
    "Receiver",
    "Transmitter",
    "SkipFlag",

    # measurements / site-level info
    "CompMeas",
    "Amps",
    "Frequency",
    "Station",

    # data-quality / variation metrics
    "PcEmag",
    "PcHmag",
    "PcRho",
    "SPhz",
    "SHphz",
    "SEphz",

    # apparent-field estimates
    "Resistivity",
    "Phase",

    # full impedance tensor stack
    "Z",
    
    # Info 
    "DataInfo"
]
