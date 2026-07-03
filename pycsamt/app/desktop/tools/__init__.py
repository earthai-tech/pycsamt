# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt desktop Tools — one-shot utility dialogs.

Each tool is a self-contained QDialog imported lazily from MainWindow so
the app starts fast even when heavy scientific back-ends are not yet loaded.

Original 5 tools
----------------
StrikeAnalyzerDialog    — regional-strike rose diagram + per-station table
EDIValidatorDialog      — per-station data-quality checklist
FormatConverterDialog   — source-folder → EDI / CSV / JSON re-export
BatchExportDialog       — export every open canvas to PNG / PDF / SVG
CoordTransformDialog    — UTM ↔ Geographic (Lat / Lon) converter

Extended 8 tools
----------------
StationResponseDialog     — full impedance tensor + tipper for one station
StrikeProfileDialog       — strike angle vs. station-position line plot
PhaseTensorMapDialog      — geographic ellipse map at a chosen period
PhaseTensorStripGridDialog — multi-profile phase-tensor ellipse-strip grid
DimensionalityDialog      — 1D / 2D / 3D classification table + bar chart
FrequencyEditorDialog     — confidence-based frequency QC workflow
LayeredModelDialog        — interactive 1-D earth model builder
ElevationEnrichDialog     — fetch elevation for all stations via open API
"""
from .strike_tool           import StrikeAnalyzerDialog
from .validator_tool        import EDIValidatorDialog
from .converter_tool        import FormatConverterDialog
from .batch_export_tool     import BatchExportDialog
from .coord_tool            import CoordTransformDialog
from .station_response_tool import StationResponseDialog
from .strike_profile_tool   import StrikeProfileDialog
from .phase_tensor_map_tool import PhaseTensorMapDialog
from .phase_tensor_strip_grid_tool import PhaseTensorStripGridDialog
from .dimensionality_tool   import DimensionalityDialog
from .frequency_editor_tool import FrequencyEditorDialog
from .layered_model_tool    import LayeredModelDialog
from .elevation_tool        import ElevationEnrichDialog

__all__ = [
    "StrikeAnalyzerDialog",
    "EDIValidatorDialog",
    "FormatConverterDialog",
    "BatchExportDialog",
    "CoordTransformDialog",
    "StationResponseDialog",
    "StrikeProfileDialog",
    "PhaseTensorMapDialog",
    "PhaseTensorStripGridDialog",
    "DimensionalityDialog",
    "FrequencyEditorDialog",
    "LayeredModelDialog",
    "ElevationEnrichDialog",
]
