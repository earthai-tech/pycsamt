
"""
Aggregate and expose SEG-EDI component classes for convenient imports.

Example
-------
from pycsamt.seg.components import (
    Head, Info, DefineMeasurement,
    Hmeasurement, Emeasurement,
    MTEMAP, Spectra, TimeSeries,
    References, Copyright, Source,
    Processing, Software, Person,
    IsEdi, EdiComponentBase,
)
"""
from __future__ import annotations

# Base
from .base import EdiComponentBase

# Properties (metadata and validators)
from .properties import (
    References,
    Copyright,
    Source,
    Processing,
    Software,
    Person,
    IsEdi,
)

# Head & Info blocks
from .heads import Head, Info

# DefineMeasurement section
from .definemeas import DefineMeasurement

# HMEAS / EMEAS
from .he_meas import Hmeasurement, Emeasurement

# MT/EMAP section
from .mtemap import MTEMAP

# Spectra section
from .spectra import Spectra

# Time series section
from .time_series import TimeSeries


__all__ = [
    # base
    "EdiComponentBase",

    # properties
    "References",
    "Copyright",
    "Source",
    "Processing",
    "Software",
    "Person",
    "IsEdi",

    # heads
    "Head",
    "Info",

    # definemeas
    "DefineMeasurement",

    # H/E measurements
    "Hmeasurement",
    "Emeasurement",

    # mtemap
    "MTEMAP",

    # spectra
    "Spectra",

    # time series
    "TimeSeries",
]
