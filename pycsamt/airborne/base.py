# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Format-neutral scientific containers for airborne EM surveys."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Iterator

from ..core.base import CoreObject, MTBase
from ..metadata import BBox, InstrumentMeta, SurveyMeta
from .navigation import NavigationTrack

if TYPE_CHECKING:
    from ..emtf.document import EMTF


def _emtf_class():
    from ..emtf.document import EMTF

    return EMTF

__all__ = [
    "AirborneEMRecord",
    "AirborneEMLine",
    "AirborneEMDataset",
]


def _clean_identifier(value: Any, *, name: str) -> str:
    text = str(value).strip()
    if not text:
        raise ValueError(f"{name} must be non-empty")
    return text


@dataclass(repr=False)
class AirborneEMRecord(CoreObject):
    """One sample-aligned airborne EM scientific record.

    ``EMTF`` is reused as the transfer-function payload so MobileMT, ZTEM,
    AFMAG, and future passive systems do not need parallel matrix classes.
    Raw or auxiliary field values may be retained in :attr:`fields` until a
    technology-specific adapter establishes stronger semantics.
    """

    sample_id: str
    emtf: "EMTF | None" = None
    fields: dict[str, Any] = field(default_factory=dict)
    quality: dict[str, Any] = field(default_factory=dict)
    attrs: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        self.sample_id = _clean_identifier(
            self.sample_id,
            name="sample_id",
        )
        if self.emtf is not None and not isinstance(self.emtf, _emtf_class()):
            raise TypeError("emtf must be an EMTF instance or None")
        self.fields = dict(self.fields or {})
        self.quality = dict(self.quality or {})
        self.attrs = dict(self.attrs or {})

    @property
    def transfer_function_names(self) -> tuple[str, ...]:
        """Transfer-function names available for this sample."""
        if self.emtf is None:
            return ()
        return tuple(self.emtf.transfer_functions)


@dataclass(repr=False)
class AirborneEMLine(CoreObject):
    """One airborne flight line with navigation and sparse EM records.

    Records are keyed by navigation ``sample_id`` and may be sparse.  This is
    deliberate: missing or rejected EM samples should not require deleting the
    corresponding navigation point or fabricating a transfer function.
    """

    line_id: str
    navigation: NavigationTrack
    records: dict[str, AirborneEMRecord] = field(default_factory=dict)
    attrs: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        self.line_id = _clean_identifier(self.line_id, name="line_id")
        if not isinstance(self.navigation, NavigationTrack):
            raise TypeError("navigation must be a NavigationTrack")
        incoming = dict(self.records or {})
        self.records = {}
        for key, record in incoming.items():
            if not isinstance(record, AirborneEMRecord):
                raise TypeError("records must contain AirborneEMRecord values")
            if str(key).strip() != record.sample_id:
                raise ValueError(
                    "record mapping key must match record.sample_id"
                )
            self.add_record(record)
        self.attrs = dict(self.attrs or {})

    @property
    def n_samples(self) -> int:
        """Number of navigation samples on the line."""
        return self.navigation.n_samples

    @property
    def n_records(self) -> int:
        """Number of EM records currently attached to the line."""
        return len(self.records)

    @property
    def bbox(self) -> BBox | None:
        """Geographic bounding box when navigation coordinates exist."""
        return self.navigation.bbox

    @property
    def missing_sample_ids(self) -> tuple[str, ...]:
        """Navigation samples that currently have no EM record."""
        return tuple(
            sample_id
            for sample_id in self.navigation.sample_ids
            if sample_id not in self.records
        )

    @property
    def transfer_function_names(self) -> tuple[str, ...]:
        """Sorted union of transfer-function names present on this line."""
        names: set[str] = set()
        for record in self.records.values():
            names.update(record.transfer_function_names)
        return tuple(sorted(names))

    def add_record(
        self,
        record: AirborneEMRecord,
        *,
        replace: bool = False,
    ) -> "AirborneEMLine":
        """Attach one record after verifying navigation alignment."""
        if not isinstance(record, AirborneEMRecord):
            raise TypeError("record must be an AirborneEMRecord")
        self.navigation.index_of(record.sample_id)
        if record.sample_id in self.records and not replace:
            raise ValueError(
                f"airborne record already exists: {record.sample_id}"
            )
        self.records[record.sample_id] = record
        return self

    def add_emtf(
        self,
        sample_id: str,
        emtf: "EMTF",
        *,
        fields: dict[str, Any] | None = None,
        quality: dict[str, Any] | None = None,
        attrs: dict[str, Any] | None = None,
        replace: bool = False,
    ) -> AirborneEMRecord:
        """Convenience helper for attaching an :class:`EMTF` sample."""
        record = AirborneEMRecord(
            sample_id=str(sample_id),
            emtf=emtf,
            fields={} if fields is None else fields,
            quality={} if quality is None else quality,
            attrs={} if attrs is None else attrs,
        )
        self.add_record(record, replace=replace)
        return record

    def get_record(self, sample_id: str) -> AirborneEMRecord | None:
        """Return a record by sample identifier, or ``None`` when absent."""
        key = str(sample_id).strip()
        self.navigation.index_of(key)
        return self.records.get(key)

    def record_at(self, index: int) -> AirborneEMRecord | None:
        """Return the record aligned with a navigation index."""
        sample_id = self.navigation.sample_ids[index]
        return self.records.get(sample_id)

    def iter_records(self) -> Iterator[AirborneEMRecord]:
        """Iterate records in navigation order, skipping missing samples."""
        for sample_id in self.navigation.sample_ids:
            record = self.records.get(sample_id)
            if record is not None:
                yield record


@dataclass(repr=False)
class AirborneEMDataset(MTBase):
    """Format-neutral collection of airborne EM flight lines.

    The dataset is intentionally an organisational layer above :class:`EMTF`.
    It does not define a MobileMT, ZTEM, or AFMAG file schema.  Technology
    adapters added later should populate this object rather than introducing
    separate transfer-function mathematics.
    """

    name: str
    lines: dict[str, AirborneEMLine] = field(default_factory=dict)
    survey: SurveyMeta | None = None
    instrument: InstrumentMeta | None = None
    method: str = "AEM"
    attrs: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        self.name = _clean_identifier(self.name, name="name")
        method = str(self.method).strip().upper()
        if not method:
            raise ValueError("method must be non-empty")
        self.method = method
        if self.survey is not None and not isinstance(self.survey, SurveyMeta):
            raise TypeError("survey must be a SurveyMeta or None")
        if self.instrument is not None and not isinstance(
            self.instrument,
            InstrumentMeta,
        ):
            raise TypeError("instrument must be an InstrumentMeta or None")

        incoming = dict(self.lines or {})
        self.lines = {}
        for key, line in incoming.items():
            if not isinstance(line, AirborneEMLine):
                raise TypeError("lines must contain AirborneEMLine values")
            if str(key).strip() != line.line_id:
                raise ValueError("line mapping key must match line.line_id")
            self.add_line(line)
        self.attrs = dict(self.attrs or {})

    @property
    def line_ids(self) -> tuple[str, ...]:
        """Flight-line identifiers in insertion order."""
        return tuple(self.lines)

    @property
    def n_lines(self) -> int:
        """Number of flight lines."""
        return len(self.lines)

    @property
    def n_samples(self) -> int:
        """Total number of navigation samples across all lines."""
        return sum(line.n_samples for line in self.lines.values())

    @property
    def n_records(self) -> int:
        """Total number of attached EM records across all lines."""
        return sum(line.n_records for line in self.lines.values())

    @property
    def transfer_function_names(self) -> tuple[str, ...]:
        """Sorted union of transfer-function types in the dataset."""
        names: set[str] = set()
        for line in self.lines.values():
            names.update(line.transfer_function_names)
        return tuple(sorted(names))

    @property
    def bbox(self) -> BBox | None:
        """Geographic bounding box over all lines with finite coordinates."""
        boxes = [line.bbox for line in self.lines.values()]
        boxes = [box for box in boxes if box is not None]
        if not boxes:
            return None
        return BBox(
            lat_min=min(box.lat_min for box in boxes),
            lat_max=max(box.lat_max for box in boxes),
            lon_min=min(box.lon_min for box in boxes),
            lon_max=max(box.lon_max for box in boxes),
        )

    def add_line(
        self,
        line: AirborneEMLine,
        *,
        replace: bool = False,
    ) -> "AirborneEMDataset":
        """Attach one flight line."""
        if not isinstance(line, AirborneEMLine):
            raise TypeError("line must be an AirborneEMLine")
        if line.line_id in self.lines and not replace:
            raise ValueError(f"airborne line already exists: {line.line_id}")
        self.lines[line.line_id] = line
        return self

    def get_line(self, line_id: str) -> AirborneEMLine | None:
        """Return a line by identifier."""
        return self.lines.get(str(line_id).strip())

    def iter_lines(self) -> Iterator[AirborneEMLine]:
        """Iterate flight lines in insertion order."""
        yield from self.lines.values()

    def iter_records(
        self,
    ) -> Iterator[tuple[str, AirborneEMRecord]]:
        """Iterate ``(line_id, record)`` pairs in navigation order."""
        for line in self.lines.values():
            for record in line.iter_records():
                yield line.line_id, record

    def emtf_records(self) -> dict[tuple[str, str], "EMTF"]:
        """Return all non-empty EMTF records keyed by line/sample ID."""
        out: dict[tuple[str, str], "EMTF"] = {}
        for line_id, record in self.iter_records():
            if record.emtf is not None:
                out[(line_id, record.sample_id)] = record.emtf
        return out

    def inspect(self):
        """Return the common airborne inspection summary lazily."""
        from .qc import inspect_airborne

        return inspect_airborne(self)

    def qc(self):
        """Return the common structural airborne QC report lazily."""
        from .qc import assess_airborne_qc

        return assess_airborne_qc(self)
