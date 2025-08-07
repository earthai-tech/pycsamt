# -*- coding: utf-8 -*-
# Author: Kouadio Laurent alias Daniel <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.exceptions

Custom exception hierarchy for pycsamt (v2.0).
"""

class PycsamtError(Exception):
    """Base exception for all pycsamt errors."""
    pass

class NotFittedError(PycsamtError):
    """
    Raised when an estimator or reader method is invoked
    before `fit()` has been called to load the underlying data.
    """
    pass

class GisError(PycsamtError):
    """Errors in GIS coordinate transformations and lookups."""
    pass

class FileIOError(PycsamtError):
    """General file I/O errors in pycsamt readers/writers."""
    pass

class ValidationError(PycsamtError):
    """Raised when input data fails validation checks."""
    pass

class Occam2DPlotError(PycsamtError):
    """Error plotting Occam2D inversion results."""
    pass

class Occam2DError(PycsamtError):
    """Error during Occam2D inversion process."""
    pass

class Occam2DIterToDataError(PycsamtError):
    """Error converting Occam2D iteration output to data format."""
    pass

class ProfileError(PycsamtError):
    """Error processing profile data."""
    pass

class JFileError(PycsamtError):
    """Error handling J-format (joint) data files."""
    pass

class SiteError(PycsamtError):
    """Error with site metadata or operations."""
    pass

class EdIDataError(PycsamtError):
    """Error reading or writing EDI files."""
    pass

class LocationError(PycsamtError):
    """Error with geographic location data."""
    pass

class WellBuildingError(PycsamtError):
    """Error constructing well-building models."""
    pass

class FloatConversionError(PycsamtError):
    """Error converting values to float."""
    pass

class FrequencyError(PycsamtError):
    """Error with frequency values or units."""
    pass

class ArgumentError(PycsamtError):
    """Invalid function or method arguments."""
    pass

class InputError(PycsamtError):
    """Invalid Inputs arguments."""
    pass


class HeaderError(PycsamtError):
    """Error processing file headers."""
    pass

class ParameterValueError(PycsamtError):
    """Invalid parameter values provided by user or file."""
    pass

class AzimuthError(PycsamtError):
    """Error in azimuth calculations or input."""
    pass

class TimeSeriesDataError(PycsamtError):
    """Error processing time-series data."""
    pass

class ConfigFileError(PycsamtError):
    """Error loading or parsing the configuration file."""
    pass

class FileHandlingError(PycsamtError):
    """General file I/O error."""
    pass

class AvgFileError(PycsamtError):
    """Error reading or writing average data files."""
    pass

class AvgDataError(PycsamtError):
    """Error processing average data content."""
    pass

class PlotError(PycsamtError):
    """General plotting error."""
    pass

class EmagError(PycsamtError):
    """Error in electric field magnitude (E_mag) processing."""
    pass

class HmagError(PycsamtError):
    """Error in magnetic field magnitude (H_mag) processing."""
    pass

class EphzError(PycsamtError):
    """Error in phase (phi-z) calculations."""
    pass

class ParameterCountError(PycsamtError):
    """Incorrect number of parameters provided."""
    pass

class ProcessingError(PycsamtError):
    """General processing or algorithm error."""
    pass

class ModuleImportError(PycsamtError):
    """Error importing a required module."""
    pass

class PhaseError(PycsamtError):
    """Error in phase-specific calculations."""
    pass

class ResistivityError(PycsamtError):
    """Error calculating resistivity (rho)."""
    pass

class StationError(PycsamtError):
    """Error with station metadata or operations."""
    pass

class StructuralError(PycsamtError):
    """Error with structural model data."""
    pass

class StrataError(PycsamtError):
    """Error with geological strata data."""
    pass

class MemoryLimitError(PycsamtError):
    """Memory-related error; insufficient resources."""
    pass

class GeoPlotArgumentError(PycsamtError):
    """Invalid arguments for geophysical plotting functions."""
    pass

class EMError(PycsamtError):
    """Base error for EM utilities."""
    pass

class GeodrillInputError(PycsamtError):
    """Error with input arguments for the geodrill module."""
    pass

class SQLError(PycsamtError):
    """General SQL database error."""
    pass

class ZError(PycsamtError):
    """General ZError."""
    pass

class ZDataError(PycsamtError):
    """General Z data Error."""
    pass

class SQLManagerError(SQLError):
    """Error managing the SQL database operations."""
    pass

class SQLInterfaceError(SQLError):
    """Error interfacing with the SQL database."""
    pass

class SQLGeoDatabaseError(SQLError):
    """Error with geo-database operations in SQL."""
    pass

class SQLUpdateGeoInfoError(SQLError):
    """Error updating geographic information in SQL database."""
    pass

class TopographyError (ProfileError): 
    pass 

class StatsError(Exception):
    """Base exception for stats module."""
    pass


__all__ = [
    "PycsamtError",
    "Occam2DPlotError",
    "Occam2DError",
    "Occam2DIterToDataError",
    "ProfileError",
    "JFileError",
    "SiteError",
    "EdIDataError",
    "LocationError",
    "WellBuildingError",
    "FloatConversionError",
    "FrequencyError",
    "ArgumentError",
    "HeaderError",
    "ParameterValueError",
    "AzimuthError",
    "TimeSeriesDataError",
    "ConfigFileError",
    "FileHandlingError",
    "AvgFileError",
    "AvgDataError",
    "PlotError",
    "EmagError",
    "HmagError",
    "EphzError",
    "ParameterCountError",
    "ProcessingError",
    "ModuleImportError",
    "PhaseError",
    "ResistivityError",
    "StationError",
    "StructuralError",
    "StrataError",
    "MemoryLimitError",
    "GeoPlotArgumentError",
    "GeodrillInputError",
    "SQLError",
    "SQLManagerError",
    "SQLInterfaceError",
    "SQLGeoDatabaseError",
    "SQLUpdateGeoInfoError",
    'NotFittedError',
    'GisError',
    'FileIOError',
    'ValidationError',
    'TopographyError', 
]
