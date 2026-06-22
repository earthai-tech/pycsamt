# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
pycsamt.tdem - Time-domain EM (TEM) processing and MT conversion
=================================================================

This subpackage converts field TEM data to frequency-domain
apparent impedance, producing
:class:`~pycsamt.seg.collection.EDICollection` objects that feed
directly into the rest of the pyCSAMT pipeline.

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
>>> coll = conv.transform(snd)                 # EDICollection (1 site)
"""

from . import io, plot
from ._base import TEMSounding
from .reader import TEMReader
from .avg import TEMAVG, TEMAVGRecord, is_temavg_file
from .coordinates import (
    TEMCoordinate,
    TEMCoordinateTable,
    read_tem_coordinates,
)
from .log import TEMLog, TEMLogRecord, is_tem_log_file
from .plot import (
    PlotDecayCurve,
    PlotElevationProfile,
    PlotGateProfile,
    PlotSurveyMap,
    PlotSurveyOverview,
    PlotTEMAVGSection,
    PlotTEMDashboard,
    PlotTEMZSection,
    PlotTransformedRho,
    StationTickConfig,
    TDEMPlotStyle,
    plot_decay,
    plot_elevation_profile,
    plot_gate_profile,
    plot_survey_map,
    plot_survey_overview,
    plot_tem_dashboard,
    plot_tem_z_section,
    plot_temavg_section,
    plot_transformed_rho,
)
from .survey import TEMSurvey, read_temavg_survey
from .transform import FourierTransform, LateTimeTransform, TEMtoEDI
from .waveform import (
    CustomWaveform,
    HalfSineWaveform,
    RampWaveform,
    SquareWaveform,
)
from .workflow import (
    TEMAVGConversion,
    read_temavg_soundings,
    transform_temavg_survey,
)
from .zplot import TEMZPlot, TEMZRecord, is_tem_z_file

__all__ = [
    # reader
    "TEMReader",
    # data model
    "TEMSounding",
    "TEMAVG",
    "TEMAVGRecord",
    "TEMCoordinate",
    "TEMCoordinateTable",
    "TEMLog",
    "TEMLogRecord",
    "TEMSurvey",
    "TEMAVGConversion",
    "TEMZPlot",
    "TEMZRecord",
    "TDEMPlotStyle",
    "StationTickConfig",
    "PlotDecayCurve",
    "PlotElevationProfile",
    "PlotTransformedRho",
    "PlotTEMAVGSection",
    "PlotTEMZSection",
    "PlotSurveyMap",
    "PlotSurveyOverview",
    "PlotGateProfile",
    "PlotTEMDashboard",
    "is_temavg_file",
    "is_tem_z_file",
    "is_tem_log_file",
    "read_temavg_survey",
    "read_tem_coordinates",
    "read_temavg_soundings",
    "transform_temavg_survey",
    "plot_decay",
    "plot_elevation_profile",
    "plot_transformed_rho",
    "plot_temavg_section",
    "plot_tem_z_section",
    "plot_survey_map",
    "plot_survey_overview",
    "plot_gate_profile",
    "plot_tem_dashboard",
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
