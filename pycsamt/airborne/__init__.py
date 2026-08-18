"""Format-neutral airborne electromagnetic survey containers.

The package provides common survey/navigation objects, a technology registry,
native-I/O extension points, and structural QC. Technology subpackages map
decoded scientific arrays into this model; native vendor readers remain
real-data driven and are not guessed from literature descriptions.
"""

from .base import AirborneEMDataset, AirborneEMLine, AirborneEMRecord
from .io import (
    AirborneIOError,
    available_airborne_readers,
    available_airborne_writers,
    read_airborne,
    write_airborne,
)
from .navigation import NavigationTrack
from .qc import (
    AirborneInspection,
    AirborneQCIssue,
    AirborneQCReport,
    assess_airborne_qc,
    inspect_airborne,
)
from .registry import (
    AirborneFormatDefinition,
    AirborneFormatDetectionError,
    AirborneRegistryError,
    AirborneTechnologyAmbiguityError,
    AirborneTechnologyDefinition,
    detect_airborne_format,
    detect_airborne_technology,
    get_airborne_format,
    get_airborne_technology,
    identify_airborne_technologies,
    list_airborne_formats,
    list_airborne_technologies,
    register_airborne_format,
    register_airborne_technology,
)

__all__ = [
    "NavigationTrack",
    "AirborneEMRecord",
    "AirborneEMLine",
    "AirborneEMDataset",
    "AirborneTechnologyDefinition",
    "AirborneFormatDefinition",
    "AirborneRegistryError",
    "AirborneTechnologyAmbiguityError",
    "AirborneFormatDetectionError",
    "register_airborne_technology",
    "get_airborne_technology",
    "list_airborne_technologies",
    "identify_airborne_technologies",
    "detect_airborne_technology",
    "register_airborne_format",
    "get_airborne_format",
    "list_airborne_formats",
    "detect_airborne_format",
    "AirborneIOError",
    "read_airborne",
    "write_airborne",
    "available_airborne_readers",
    "available_airborne_writers",
    "AirborneQCIssue",
    "AirborneInspection",
    "AirborneQCReport",
    "inspect_airborne",
    "assess_airborne_qc",
]
