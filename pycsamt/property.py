# -*- coding: utf-8 -*-
# Author: Kouadio Laurent alias Daniel <etanoyau@gmail.com>
# License: GPL-3.0
"""
pycsamt.property

Central definitions and file-type recognition utilities for pycsamt.

Contains:
 - TermDefinitions: human-readable definitions of key AMT/CSAMT concepts
 - FieldAliases: standardized lists of common column/field name variants
 - FileRecognizer: logic to infer file type based on path or content
"""
# import warnings 
from __future__ import annotations 
import os
import warnings
from abc import ( 
    ABC, 
    abstractmethod, 
    )

from typing import List, Optional, Tuple, Union

import numpy as np
import pandas as pd 

from pycsamt.log.logger import get_logger
from pycsamt.exceptions import FileHandlingError, EdIDataError

logger = get_logger(__name__)

__all__ = ['TermDefinitions', 'FieldAliases', 'FileRecognizer']

# UTM zone designator ranges: letter -> [min_lon, max_lon]
UTM_ZONE_DESIGNATOR = {
    'X': [72, 84],
    'W': [64, 72],
    'V': [56, 64],
    'U': [48, 56],
    'T': [40, 48],
    'S': [32, 40],
    'R': [24, 32],
    'Q': [16, 24],
    'P': [8, 16],
    'N': [0, 8],
    'M': [-8, 0],
    'L': [-16, -8],
    'K': [-24, -16],
    'J': [-32, -24],
    'H': [-40, -32],
    'G': [-48, -40],
    'F': [-56, -48],
    'E': [-64, -56],
    'D': [-72, -64],
    'C': [-80, -72],
    'Z': [-80, 84],
}
 

class TermDefinitions:
    """
    Definitions of technical terms used in AMT/CSAMT
    processing. Access attributes to retrieve a description
    string.
    """
    reference_frequency = (
        "Reference frequency is the highest frequency with clean "
        "data. See Weik (2001) Computer Science and Communications "
        "Dictionary for details."
    )

    trimmed_moving_average = (
        "Trimmed Moving Average filter: removes single-station "
        "offsets while preserving broader trends to estimate "
        "apparent resistivities."
    )

    avg = (
        "AVG files from Zonge Engineering: average raw CSAMT/CSHA "
        "data, including multiple auxiliary files (.LOG, .AL, .AD, "
        ".AVG, .Z, .Xnn)."
    )

    fixed_moving_average = (
        "Fixed-length Moving Average (FMLA) filter: computes static-"
        "corrected apparent resistivities at a reference frequency "
        "along a profile."
    )

    hanning_window = (
        "Hanning window: a weighted cosine taper (aka Hann window), "
        "named after Julius von Hann."
    )

    weighted_beta = (
        "Beta coefficients weight dipole responses to emulate a "
        "Hanning window via least-squares fitting."
    )

    j_format = (
        "J-format: MT data file format by A.G. Jones. "
        "E-polarization is RXY, B-polarization is RYX. See "
        "http://mtnet.info/docs/jformat.txt for details."
    )


class FieldAliases:
    """
    Common field/column name variants for ease of lookup.
    """
    missing_values: List = [' ', 'nan', np.nan, '*', 'NaN', 'none', None]

    longitude: List[str] = ['lon', 'longitude', 'LONG', 'LON']
    latitude:  List[str] = ['lat', 'latitude', 'LAT', 'LATITUDE']
    easting:   List[str] = ['e', 'east', 'easting', 'EASTING']
    northing:  List[str] = ['n', 'north', 'northing', 'NORTHING']
    station:   List[str] = ['sta', 'station', 'stn']
    elevation: List[str] = ['elev', 'elevation', 'ELEV', 'ELEVATION']
    azimuth:   List[str] = ['azim', 'azimuth']


class FileRecognizer:
    """
    Infer file type from file path or content.

    Supported types: 'avg', 'j', 'edi', 'stn', 'resp', 'mesh', 'model',
    'startup', 'iter', 'logfile'.
    """
    _type_signatures = {
        'edi': ['>HEAD', '>END'],
        'j':   ['>AZIMUTH', 'RXX', 'RXY'],
        'avg': ['AMTAVG', 'skp', '%Rho'],
        'stn': ['Station', 'Freq'],
        'resp': None,  # numeric-only rows
        'mesh': ['ZZZZZZZZZZZZ', '????????????'],
        'model': ['MODEL NAME', 'NUM LAYERS'],
        'startup': ['STARTUP', 'ITERATION'],
        'iter': ['Iteration', 'Misfit'],
        'logfile': ['MISFIT', 'ROUGHNESS'],
    }

    @classmethod
    def recognize(cls, filepath: str, deep: bool = True) -> str:
        """
        Recognize file type by extension (shallow) or content (deep).

        Parameters
        ----------
        filepath : str
            Path to the file to analyze.
        deep : bool
            If False, use extension only; if True, inspect file contents.

        Returns
        -------
        str
            One of the supported file types.

        Raises
        ------
        FileHandlingError
            If file is unreadable or type cannot be determined.
        """
        if not filepath or not os.path.isfile(filepath):
            raise FileHandlingError(f"Invalid file path: {filepath}")

        ext = os.path.splitext(filepath)[1].lstrip('.').lower()
        if not deep and ext in cls._type_signatures:
            return ext

        # Read lines
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                lines = f.read().splitlines()
        except Exception as e:
            raise FileHandlingError(f"Cannot read file: {e}")

        # Try signature matching
        for ftype, sigs in cls._type_signatures.items():
            if sigs is None and cls._is_numeric_file(lines):
                return ftype
            if sigs and any(s in line for s in sigs for line in lines[:10]):
                return ftype

        raise FileHandlingError(f"Unrecognized file format: {filepath}")

    @staticmethod
    def _is_numeric_file(lines: List[str]) -> bool:
        """
        Heuristic: all tokens in first and last lines are numeric.
        """
        for line in (lines[0], lines[-1]):
            tokens = line.split()
            if not tokens:
                return False
            try:
                _ = [float(t) for t in tokens]
            except ValueError:
                return False
        return True



