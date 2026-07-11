"""EDI section models for metadata, spectra, time series, and impedance data."""

from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING

__all__ = [
    "EDIMixin", "EDIOMixin", "EDIFile",
    "Spectra", "SpectraSECT", "SpectraIO",
    "TimeSeries", "TSect", "TSIO",
    "OtherSECT", "OtherIO",
    "IsEdi",
    "SurveyBase", "EDIProfile", "Stations", "Topography",
    "XAMixin", "build_dataset", "EDIAcc",
]

# Map exported names to (module, attribute)
_EXPORTS = {
    "EDIMixin": ("edi", "EDIMixin"),
    "EDIOMixin": ("edi", "EDIOMixin"),
    "EDIFile": ("edi", "EDIFile"),
    "Spectra": ("spectra", "Spectra"),
    "SpectraSECT": ("spectra", "SpectraSECT"),
    "SpectraIO": ("spectra", "SpectraIO"),
    "TimeSeries": ("time_series", "TimeSeries"),
    "TSect": ("time_series", "TSect"),
    "TSIO": ("time_series", "TSIO"),
    "OtherSECT": ("other", "OtherSECT"),
    "OtherIO": ("other", "OtherIO"),
    "IsEdi": ("validation", "IsEdi"),
    "SurveyBase": ("survey", "SurveyBase"),
    "EDIProfile": ("survey", "EDIProfile"),
    "Stations": ("survey", "Stations"),
    "Topography": ("survey", "Topography"),
    "XAMixin": ("xa", "XAMixin"),
    "build_dataset": ("xa", "build_dataset"),
    "EDIAcc": ("xa", "EDIAcc"),
}

def __getattr__(name: str):
    if name in _EXPORTS:
        mod, attr = _EXPORTS[name]
        try:
            m = import_module(f".{mod}", package=__name__)
            val = getattr(m, attr)
        except Exception as exc:  # give a helpful message
            raise AttributeError(
                f"{name} is unavailable because importing '{mod}' failed. "
                f"Install optional deps "
                "(e.g. xarray/matplotlib) or check circular imports."
            ) from exc
        globals()[name] = val  # cache after first access
        return val
    raise AttributeError(name)

# Type checkers/IDE intellisense still see the symbols:
if TYPE_CHECKING:
    from .edi import EDIFile, EDIMixin, EDIOMixin
    from .other import OtherIO, OtherSECT
    from .spectra import Spectra, SpectraIO, SpectraSECT
    from .survey import (
        EDIProfile,
        Stations,
        SurveyBase,
        Topography,
    )
    from .time_series import TSIO, TimeSeries, TSect
    from .validation import IsEdi
    from .xa import EDIAcc, XAMixin, build_dataset
