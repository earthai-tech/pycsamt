"""Survey metadata, frequency bands, references, and software provenance."""

from ._property import (
    CopyrightInfo,
    Person,
    Reference,
    Software,
)
from .frequency import (
    MT_BANDS,
    REGISTRY,
    FrequencyBand,
    band_for_frequency,
    doi_estimate,
    frequency_range,
    register_band,
)
from .geology import (
    CATALOG,
    Formation,
    GeologyCatalog,
    geology_prior,
)
from .instrument import (
    KNOWN_SYSTEMS,
    InstrumentMeta,
    SensorSpec,
    known_system,
    list_presets,
)
from .quality import (
    ComponentQuality,
    DataQuality,
    QualityFlag,
    assess_collection,
    quality_dataframe,
)
from .rocks import (
    GEO_ROCK_RESISTIVITY,
    ROCK_HATCH_PATTERNS,
    RockProperties,
)
from .survey import BBox, SurveyMeta

__all__ = [
    # legacy / bibliographic
    "Reference",
    "CopyrightInfo",
    "Person",
    "Software",
    # rocks
    "GEO_ROCK_RESISTIVITY",
    "ROCK_HATCH_PATTERNS",
    "RockProperties",
    # survey
    "BBox",
    "SurveyMeta",
    # geology
    "Formation",
    "GeologyCatalog",
    "CATALOG",
    "geology_prior",
    # quality
    "QualityFlag",
    "ComponentQuality",
    "DataQuality",
    "assess_collection",
    "quality_dataframe",
    # frequency
    "FrequencyBand",
    "MT_BANDS",
    "REGISTRY",
    "band_for_frequency",
    "frequency_range",
    "doi_estimate",
    "register_band",
    # instrument
    "SensorSpec",
    "InstrumentMeta",
    "KNOWN_SYSTEMS",
    "known_system",
    "list_presets",
]
