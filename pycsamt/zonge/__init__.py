# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
The `zonge` subpackage provides a comprehensive toolkit for
reading, processing, and analyzing Zonge AVG/AMTAVG data files.
"""

from .avg import AMTAVG, AVG, BaseAVG
from .base import (
    AVGComponentBase,
    AVGFrame,
    AvgRow,
    FieldAliases,
)
from .config import Zonge
from .info import DataInfo
from .meas import Amps, CompMeas, Frequency
from .proc_utils import (
    ama,
    flma,
    get_reference_frequency,
    get_skew,
    get_strike,
    interpolate_to_log_space,
    smooth_rho_from_phase,
    tma,
)
from .processing import ASTATIC
from .property import (
    Hardware,
    Receiver,
    SkipFlag,
    SurveyAnnotation,
    SurveyConfiguration,
    Transmitter,
)
from .resphase import Phase, Resistivity
from .schema import _CANON_TO_LEGACY as CANON_TO_LEGACY_MAP
from .schema import _CANON_TO_MODERN as CANON_TO_MODERN_MAP
from .schema import _CANONICAL_MAP as CANONICAL_MAP
from .schema import ALL_ALIASES, QC_ALIASES, get_aliases
from .survey import Station
from .tensor import TensorBase
from .tipper import Tipper
from .utils import load_avg, write_avg

__all__ = [
    "AVG",
    "AMTAVG",
    "BaseAVG",
    "Zonge",
    "AVGFrame",
    "AvgRow",
    "FieldAliases",
    "AVGComponentBase",
    "DataInfo",
    "Amps",
    "CompMeas",
    "Frequency",
    "ASTATIC",
    "Hardware",
    "SurveyAnnotation",
    "SurveyConfiguration",
    "Receiver",
    "Transmitter",
    "SkipFlag",
    "Resistivity",
    "Phase",
    "CANONICAL_MAP",
    "CANON_TO_MODERN_MAP",
    "CANON_TO_LEGACY_MAP",
    "ALL_ALIASES",
    "QC_ALIASES",
    "get_aliases",
    "Station",
    "TensorBase",
    "Tipper",
    "load_avg",
    "write_avg",
    "tma",
    "flma",
    "ama",
    "interpolate_to_log_space",
    "smooth_rho_from_phase",
    "get_reference_frequency",
    "get_skew",
    "get_strike",
]
