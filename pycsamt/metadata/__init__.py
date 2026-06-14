from ._property import (
    Reference,
    CopyrightInfo,
    Person,
    Software,
)
from .rocks import (
    GEO_ROCK_RESISTIVITY,
    ROCK_HATCH_PATTERNS,
    RockProperties,
)
from .survey import BBox, SurveyMeta
from .geology import Formation, GeologyCatalog, CATALOG, geology_prior
from .quality import (
    QualityFlag,
    ComponentQuality,
    DataQuality,
    assess_collection,
    quality_dataframe,
)
from .frequency import (
    FrequencyBand,
    MT_BANDS,
    REGISTRY,
    band_for_frequency,
    frequency_range,
    doi_estimate,
    register_band,
)
from .instrument import (
    SensorSpec,
    InstrumentMeta,
    KNOWN_SYSTEMS,
    known_system,
    list_presets,
)

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
