# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
pycsamt.tdem — Time-domain EM (TEM) processing and MT conversion
=================================================================

This subpackage converts field TEM data to frequency-domain
apparent impedance, producing :class:`~pycsamt.seg.collection.EDICollection`
objects that feed directly into the rest of the pyCSAMT pipeline
(processing, inversion, interpretation).

Workflow
--------
1. Load raw TEM data with a reader from :mod:`pycsamt.tdem.io`.
2. (Optional) attach a transmitter waveform from
   :mod:`pycsamt.tdem.waveform`.
3. Convert to frequency domain with
   :class:`~pycsamt.tdem.transform.TEMtoEDI`.
4. Inspect the result with :mod:`pycsamt.tdem.plot`.

Quick example
-------------
>>> import numpy as np
>>> from pycsamt.tdem import TEMSounding, TEMtoEDI
>>> t = np.logspace(-5, -2, 30)
>>> dBdt = 5e-5 * t ** (-5.0 / 2.0)           # synthetic decay
>>> snd = TEMSounding(
...     t, dBdt,
...     current=8.0, tx_area=100.0 ** 2,
...     data_type="dBdt",
...     station_name="S01", x=1000.0, y=500.0,
... )
>>> conv = TEMtoEDI(method="late_time", phase_mode="weidelt")
>>> coll = conv.transform(snd)                 # → EDICollection (1 site)
"""

from ._base import TEMSounding
from .transform import LateTimeTransform, FourierTransform, TEMtoEDI
from .waveform import SquareWaveform, RampWaveform, HalfSineWaveform, CustomWaveform
from . import io
from . import plot

__all__ = [
    # data model
    "TEMSounding",
    # transforms
    "LateTimeTransform",
    "FourierTransform",
    "TEMtoEDI",
    # waveforms
    "SquareWaveform",
    "RampWaveform",
    "HalfSineWaveform",
    "CustomWaveform",
    # sub-modules
    "io",
    "plot",
]
