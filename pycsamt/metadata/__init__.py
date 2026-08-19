"""Survey metadata, frequency bands, references, and software provenance."""

from ._property import (
    CopyrightInfo,
    MetadataError,
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
from .channels import ChannelMeta, SiteLayout
from .instrument import (
    KNOWN_SYSTEMS,
    InstrumentMeta,
    SensorSpec,
    known_system,
    list_presets,
)
from .processing import (
    ProcessingMeta,
    RemoteReferenceMeta,
    normalize_sign_convention,
)
from .provenance import ProvenanceMeta
from .quality import (
    ComponentQuality,
    DataQuality,
    QualityComment,
    QualityFlag,
    TransferFunctionQuality,
    assess_collection,
    quality_dataframe,
)
from .rocks import (
    GEO_ROCK_RESISTIVITY,
    ROCK_HATCH_PATTERNS,
    RockProperties,
)
from .orientation import OrientationMeta
from .site import LocationMeta, SiteMeta
from .survey import BBox, SurveyMeta

__all__ = [
    # legacy / bibliographic
    "Reference",
    "CopyrightInfo",
    "Person",
    "Software",
    "MetadataError",
    # transfer-function metadata
    "ProvenanceMeta",
    "LocationMeta",
    "SiteMeta",
    "OrientationMeta",
    "ChannelMeta",
    "SiteLayout",
    "RemoteReferenceMeta",
    "ProcessingMeta",
    "normalize_sign_convention",
    "QualityComment",
    "TransferFunctionQuality",
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
