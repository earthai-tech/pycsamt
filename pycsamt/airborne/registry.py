# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Technology and native-format registries for airborne EM data.

The registry separates two concepts deliberately:

* a *technology* describes scientific semantics (MobileMT, ZTEM, AFMAG);
* a *format* describes a concrete native delivery that can be read/written.

The built-in technologies are registered immediately because their scientific
contracts are stable. No native vendor format is registered until a genuine
sample or authoritative format specification is available.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

from ..core.base import CoreObject
from ..core._named_registry import (
    NamedRegistry,
    normalize_extension as _normal_extension,
    normalize_key as _normal_key,
)

__all__ = [
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
]


class AirborneRegistryError(ValueError):
    """Base exception for airborne technology/format registry errors."""


class AirborneTechnologyAmbiguityError(AirborneRegistryError):
    """Raised when an object contains more than one airborne technology."""


class AirborneFormatDetectionError(AirborneRegistryError):
    """Raised when native airborne format detection is ambiguous."""


@dataclass(frozen=True)
class AirborneTechnologyDefinition:
    """Describe one scientific airborne-EM technology contract."""

    name: str
    label: str
    family: str
    aliases: tuple[str, ...] = field(default_factory=tuple)
    primary_tf_names: tuple[str, ...] = field(default_factory=tuple)
    reference_required: bool = False
    infer_from_tf: bool = False
    description: str = ""

    def __post_init__(self) -> None:
        name = _normal_key(self.name)
        label = str(self.label).strip()
        family = _normal_key(self.family)
        aliases = tuple(_normal_key(value) for value in self.aliases)
        primary = tuple(_normal_key(value) for value in self.primary_tf_names)
        if not name or not label or not family:
            raise ValueError("technology name, label, and family are required")
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "label", label)
        object.__setattr__(self, "family", family)
        object.__setattr__(self, "aliases", aliases)
        object.__setattr__(self, "primary_tf_names", primary)
        object.__setattr__(
            self,
            "description",
            str(self.description).strip(),
        )


@dataclass(frozen=True)
class AirborneFormatDefinition:
    """Describe one concrete native airborne delivery format."""

    name: str
    technology: str
    reader: Callable[..., Any] | None = None
    writer: Callable[..., Any] | None = None
    detector: Callable[[Any], bool] | None = None
    extensions: tuple[str, ...] = field(default_factory=tuple)
    aliases: tuple[str, ...] = field(default_factory=tuple)
    description: str = ""

    def __post_init__(self) -> None:
        name = _normal_key(self.name)
        technology = _normal_key(self.technology)
        aliases = tuple(_normal_key(value) for value in self.aliases)
        extensions = tuple(
            _normal_extension(value) for value in self.extensions
        )
        if not name or not technology:
            raise ValueError("format name and technology are required")
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "technology", technology)
        object.__setattr__(self, "aliases", aliases)
        object.__setattr__(self, "extensions", extensions)
        object.__setattr__(
            self,
            "description",
            str(self.description).strip(),
        )

    @property
    def readable(self) -> bool:
        """Whether a native reader has been registered."""
        return self.reader is not None

    @property
    def writable(self) -> bool:
        """Whether a native writer has been registered."""
        return self.writer is not None


_tech_registry: NamedRegistry[AirborneTechnologyDefinition] = NamedRegistry(
    kind="airborne technology",
    error_cls=AirborneRegistryError,
)
_format_registry: NamedRegistry[AirborneFormatDefinition] = NamedRegistry(
    kind="airborne format",
    error_cls=AirborneRegistryError,
)


def register_airborne_technology(
    definition: AirborneTechnologyDefinition,
    *,
    replace: bool = False,
) -> AirborneTechnologyDefinition:
    """Register one scientific technology definition."""
    if not isinstance(definition, AirborneTechnologyDefinition):
        raise TypeError("definition must be AirborneTechnologyDefinition")
    _tech_registry.register(
        definition.name,
        definition,
        aliases=definition.aliases,
        replace=replace,
    )
    return definition


def get_airborne_technology(
    name: str,
) -> AirborneTechnologyDefinition | None:
    """Return a technology definition by canonical name or alias."""
    return _tech_registry.get(name)


def list_airborne_technologies() -> tuple[AirborneTechnologyDefinition, ...]:
    """Return registered technologies in registration order."""
    return _tech_registry.all()


def _technology_from_text(value: Any) -> str | None:
    if value is None:
        return None
    definition = get_airborne_technology(str(value))
    return None if definition is None else definition.name


def _collect_object_technologies(obj: Any, out: set[str]) -> None:
    attrs = getattr(obj, "attrs", None)
    if isinstance(attrs, dict):
        technology = _technology_from_text(attrs.get("technology"))
        if technology is not None:
            out.add(technology)

    subtype = getattr(obj, "subtype", None)
    subtype_key = _normal_key(subtype) if subtype is not None else ""
    subtype_map = {
        "mobilemt": "mobilemt",
        "ztem": "ztem",
        "afmag_original": "afmag",
        "afmag_airmt": "airmt",
    }
    if subtype_key in subtype_map:
        out.add(subtype_map[subtype_key])

    transfer_functions = getattr(obj, "transfer_functions", None)
    if isinstance(transfer_functions, dict):
        for tf in transfer_functions.values():
            tf_attrs = getattr(tf, "attrs", None)
            if isinstance(tf_attrs, dict):
                technology = _technology_from_text(
                    tf_attrs.get("technology")
                )
                if technology is not None:
                    out.add(technology)
            tf_name = _normal_key(getattr(tf, "name", ""))
            for definition in _tech_registry.all():
                if (
                    definition.infer_from_tf
                    and tf_name in definition.primary_tf_names
                ):
                    out.add(definition.name)

    emtf = getattr(obj, "emtf", None)
    if emtf is not None:
        _collect_object_technologies(emtf, out)

    records = getattr(obj, "records", None)
    if isinstance(records, dict):
        for record in records.values():
            _collect_object_technologies(record, out)

    lines = getattr(obj, "lines", None)
    if isinstance(lines, dict):
        for line in lines.values():
            _collect_object_technologies(line, out)


def identify_airborne_technologies(obj: Any) -> tuple[str, ...]:
    """Return canonical technologies explicitly or safely identified.

    Only response types that are unique to one technology are inferred from
    transfer-function names. Standard tipper ``T`` and interstation ``TI``
    are intentionally not enough by themselves to identify ZTEM or AirMt.
    """
    out: set[str] = set()
    _collect_object_technologies(obj, out)
    order = {
        definition.name: index
        for index, definition in enumerate(_tech_registry.all())
    }
    return tuple(sorted(out, key=lambda name: order.get(name, 10**9)))


def detect_airborne_technology(
    obj: Any,
    *,
    strict: bool = True,
) -> str | None:
    """Return one canonical technology or report ambiguous mixed content."""
    technologies = identify_airborne_technologies(obj)
    if len(technologies) == 1:
        return technologies[0]
    if not technologies:
        return None
    if strict:
        raise AirborneTechnologyAmbiguityError(
            "multiple airborne technologies are present: "
            + ", ".join(technologies)
        )
    return None


def register_airborne_format(
    definition: AirborneFormatDefinition,
    *,
    replace: bool = False,
) -> AirborneFormatDefinition:
    """Register a concrete native airborne file/delivery format."""
    if not isinstance(definition, AirborneFormatDefinition):
        raise TypeError("definition must be AirborneFormatDefinition")
    technology = get_airborne_technology(definition.technology)
    if technology is None:
        raise AirborneRegistryError(
            "format technology is not registered: "
            f"{definition.technology!r}"
        )
    if definition.technology != technology.name:
        definition = AirborneFormatDefinition(
            name=definition.name,
            technology=technology.name,
            reader=definition.reader,
            writer=definition.writer,
            detector=definition.detector,
            extensions=definition.extensions,
            aliases=definition.aliases,
            description=definition.description,
        )
    _format_registry.register(
        definition.name,
        definition,
        aliases=definition.aliases,
        replace=replace,
    )
    return definition


def get_airborne_format(name: str) -> AirborneFormatDefinition | None:
    """Return a native format definition by name or alias."""
    return _format_registry.get(name)


def list_airborne_formats(
    *,
    technology: str | None = None,
) -> tuple[AirborneFormatDefinition, ...]:
    """Return registered native formats, optionally for one technology."""
    if technology is None:
        return _format_registry.all()
    definition = get_airborne_technology(technology)
    if definition is None:
        raise AirborneRegistryError(f"unknown technology: {technology!r}")
    return tuple(
        fmt for fmt in _format_registry.all()
        if fmt.technology == definition.name
    )


def _extension_of(source: Any) -> str | None:
    if isinstance(source, (str, Path)):
        try:
            return Path(source).suffix.lower() or None
        except (TypeError, ValueError):
            return None
    return None


def detect_airborne_format(source: Any) -> str | None:
    """Detect a registered native format using detectors, then extensions.

    Detectors have priority. Extensions are only hints and are used when they
    map to exactly one registered format. No built-in vendor formats are
    registered merely from published system descriptions.
    """
    matches = _format_registry.match_by_detector(source)
    if len(matches) > 1:
        raise AirborneFormatDetectionError(
            "multiple airborne formats matched source: "
            + ", ".join(matches)
        )
    if matches:
        return matches[0]

    extension = _extension_of(source)
    if extension is None:
        return None
    extension_matches = _format_registry.match_by_extension(extension)
    if len(extension_matches) > 1:
        raise AirborneFormatDetectionError(
            "file extension is ambiguous across airborne formats: "
            + ", ".join(extension_matches)
        )
    return extension_matches[0] if extension_matches else None


def _register_builtin_technologies() -> None:
    builtins = (
        AirborneTechnologyDefinition(
            name="mobilemt",
            label="MobileMT",
            family="natural_field_airborne_em",
            aliases=("mobile_mt",),
            primary_tf_names=("mobilemt_admittance",),
            reference_required=True,
            infer_from_tf=True,
            description="Three-output/two-input magnetic/electric admittance.",
        ),
        AirborneTechnologyDefinition(
            name="ztem",
            label="ZTEM",
            family="natural_field_airborne_em",
            aliases=("z_tem",),
            primary_tf_names=("tipper",),
            reference_required=True,
            infer_from_tf=False,
            description="Airborne Hz to fixed-ground Hx/Hy tipper.",
        ),
        AirborneTechnologyDefinition(
            name="afmag",
            label="AFMAG (original comparator)",
            family="afmag",
            aliases=("original_afmag", "comparator_afmag"),
            primary_tf_names=("afmag_tilt",),
            reference_required=False,
            infer_from_tf=True,
            description="Historical scalar comparator/polarization tilt.",
        ),
        AirborneTechnologyDefinition(
            name="airmt",
            label="AirMt / tensor AFMAG",
            family="afmag",
            aliases=("tensor_afmag", "afmag_tensor"),
            primary_tf_names=("interstation_transfer_functions",),
            reference_required=True,
            infer_from_tf=False,
            description="Three-output/two-input magnetic interstation TF.",
        ),
    )
    for definition in builtins:
        if get_airborne_technology(definition.name) is None:
            register_airborne_technology(definition)


_register_builtin_technologies()
