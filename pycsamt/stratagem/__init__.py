# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.stratagem
=================

Tools for ingesting and pre-processing AMT surveys collected with
the Geometrics/EMI Stratagem hardware.

Workflow
--------
Raw Stratagem files (.HXX / .MDX)
  → WinGLink EDI export  (external tool, required once)
    → :class:`EDIBatch`              load the EDI directory
      → :class:`CoordinateInjector` inject GPS coords (CSV / XLS / XLSX)
        → export to corrected EDIs   ready for emtools processing

Classes
-------
StratagemRawReader
    Parse raw 19-column Stratagem ASCII files for QC diagnostics and
    SNR masks.  Does not produce EDI output.
EDIBatch
    Load a WinGLink-exported EDI directory as a list of
    :class:`~pycsamt.seg.edi.EDIFile` objects with natural-sort ordering.
CoordinateInjector
    Read a GPS coordinate table (CSV / XLS / XLSX), convert projected
    coordinates to WGS84, and inject them into EDI ``>HEAD`` sections.
StationLocator
    Detect or specify the mapping between EDI acquisition order and GPS
    table row order (forward / reversed / custom index).
"""

from .gis_correct import CoordinateInjector, StationLocator
from .io import EDIBatch, StratagemRawReader
from .process import NoiseRemover, StaticShiftCorrector
from .qc import FrequencyFilter, QualityController
from .rename import EDIRenamer, EDIWriter
from .survey import StratagemSurvey

__all__ = [
    # Phase 1 – I/O & coordinate injection
    "EDIBatch",
    "StratagemRawReader",
    "CoordinateInjector",
    "StationLocator",
    # Phase 2 – QC & processing
    "FrequencyFilter",
    "QualityController",
    "NoiseRemover",
    "StaticShiftCorrector",
    # Phase 3 – rename & orchestration
    "EDIRenamer",
    "EDIWriter",
    "StratagemSurvey",
]
