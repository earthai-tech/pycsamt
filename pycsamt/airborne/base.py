# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Format-neutral scientific containers for airborne EM surveys.

Three dataclasses form the common in-memory model that every technology
adapter (:mod:`pycsamt.airborne.mobilemt`, :mod:`pycsamt.airborne.ztem`,
:mod:`pycsamt.airborne.afmag`) populates instead of inventing its own
survey/line/sample containers:

* :class:`AirborneEMRecord` -- one sample-aligned EM response;
* :class:`AirborneEMLine` -- one flight line of navigation plus sparse
  records;
* :class:`AirborneEMDataset` -- the survey-level collection of lines.

All three inherit :class:`~pycsamt.core.base.CoreObject` rather than
:class:`~pycsamt.core.base.MTBase`: they organize and index scientific
content but do not themselves perform electromagnetic arithmetic. That
arithmetic lives in :class:`~pycsamt.emtf.document.EMTF` and
:class:`~pycsamt.emtf.transfer.TransferFunction`, which these containers
hold rather than duplicate.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

from ..core.base import CoreObject
from ..metadata import BBox, InstrumentMeta, SurveyMeta
from .navigation import NavigationTrack
from .validation import clean_identifier, emtf_class

if TYPE_CHECKING:
    from ..emtf.document import EMTF

__all__ = [
    "AirborneEMRecord",
    "AirborneEMLine",
    "AirborneEMDataset",
]


@dataclass(repr=False)
class AirborneEMRecord(CoreObject):
    """One sample-aligned airborne EM scientific record.

    Parameters
    ----------
    sample_id : str
        Identifier matching one entry of the owning line's
        ``navigation.sample_ids``. Stripped and required to be
        non-empty.
    emtf : EMTF, optional
        Transfer-function payload for this sample. ``EMTF`` is reused
        directly rather than introduced as a parallel matrix class, so
        MobileMT, ZTEM, AFMAG, and future passive systems share one
        scientific representation instead of duplicating it.
    fields : dict, optional
        Auxiliary decoded scalar/array fields with no stronger
        scientific type yet, for example a processed apparent
        conductivity vector. Content is technology-defined.
    quality : dict, optional
        Sample-level quality flags or scores. Content is
        technology-defined.
    attrs : dict, optional
        Free-form extension metadata.

    Raises
    ------
    ValueError
        If ``sample_id`` is empty after stripping.
    TypeError
        If ``emtf`` is supplied and is not an
        :class:`~pycsamt.emtf.EMTF` instance.

    Examples
    --------
    >>> from pycsamt.airborne import AirborneEMRecord
    >>> record = AirborneEMRecord(sample_id=" S001 ")
    >>> record.sample_id
    'S001'
    >>> record.transfer_function_names
    ()
    """

    sample_id: str
    emtf: EMTF | None = None
    fields: dict[str, Any] = field(default_factory=dict)
    quality: dict[str, Any] = field(default_factory=dict)
    attrs: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        """Normalize identifier/dict fields and check the EMTF type."""
        self.sample_id = clean_identifier(
            self.sample_id,
            name="sample_id",
        )
        if self.emtf is not None and not isinstance(self.emtf, emtf_class()):
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

    Parameters
    ----------
    line_id : str
        Flight-line identifier. Stripped and required to be non-empty.
    navigation : NavigationTrack
        Sample-aligned navigation/attitude track defining the line's
        common sample axis. Every record's ``sample_id`` must appear
        in ``navigation.sample_ids``.
    records : dict of str to AirborneEMRecord, optional
        Records keyed by ``sample_id``. The mapping key must equal
        ``record.sample_id`` for every entry.
    attrs : dict, optional
        Free-form extension metadata.

    Raises
    ------
    ValueError
        If ``line_id`` is empty, a record key does not match
        ``record.sample_id``, or a record's ``sample_id`` is not a
        known navigation sample.
    TypeError
        If ``navigation`` is not a :class:`NavigationTrack`, or
        ``records`` contains a non-:class:`AirborneEMRecord` value.

    Notes
    -----
    Records are keyed by navigation ``sample_id`` and may be sparse.
    This is deliberate: a missing or rejected EM sample should not
    require deleting the corresponding navigation point, nor
    fabricating a transfer function to fill the gap.

    Examples
    --------
    >>> from pycsamt.airborne import AirborneEMLine, NavigationTrack
    >>> nav = NavigationTrack(sample_ids=("S1", "S2"))
    >>> line = AirborneEMLine(line_id="L001", navigation=nav)
    >>> line.n_samples, line.n_records
    (2, 0)
    >>> line.missing_sample_ids
    ('S1', 'S2')
    """

    line_id: str
    navigation: NavigationTrack
    records: dict[str, AirborneEMRecord] = field(default_factory=dict)
    attrs: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        """Normalize the identifier and re-attach incoming records."""
        self.line_id = clean_identifier(self.line_id, name="line_id")
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
    ) -> AirborneEMLine:
        """Attach one record after verifying navigation alignment.

        Parameters
        ----------
        record : AirborneEMRecord
            Record whose ``sample_id`` must already exist on
            ``navigation``.
        replace : bool, default False
            Whether to overwrite an existing record for the same
            sample instead of raising.

        Returns
        -------
        AirborneEMLine
            ``self``, to support call chaining.

        Raises
        ------
        TypeError
            If ``record`` is not an :class:`AirborneEMRecord`.
        KeyError
            If ``record.sample_id`` is not a known navigation sample.
        ValueError
            If a record already exists for that sample and
            ``replace`` is ``False``.
        """
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
        emtf: EMTF,
        *,
        fields: dict[str, Any] | None = None,
        quality: dict[str, Any] | None = None,
        attrs: dict[str, Any] | None = None,
        replace: bool = False,
    ) -> AirborneEMRecord:
        """Build and attach one :class:`AirborneEMRecord` from an EMTF.

        Convenience wrapper around :meth:`add_record` for the common
        case of attaching a decoded EMTF response without constructing
        the record explicitly.

        Parameters
        ----------
        sample_id : str
            Navigation sample identifier for the new record.
        emtf : EMTF
            Transfer-function payload for the sample.
        fields, quality, attrs : dict, optional
            Forwarded to :class:`AirborneEMRecord`.
        replace : bool, default False
            Forwarded to :meth:`add_record`.

        Returns
        -------
        AirborneEMRecord
            The record that was attached.
        """
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
        """Return a record by sample identifier, or ``None`` when absent.

        Raises
        ------
        KeyError
            If ``sample_id`` is not a known navigation sample.
        """
        key = str(sample_id).strip()
        self.navigation.index_of(key)
        return self.records.get(key)

    def record_at(self, index: int) -> AirborneEMRecord | None:
        """Return the record aligned with navigation index *index*."""
        sample_id = self.navigation.sample_ids[index]
        return self.records.get(sample_id)

    def iter_records(self) -> Iterator[AirborneEMRecord]:
        """Iterate records in navigation order, skipping missing samples."""
        for sample_id in self.navigation.sample_ids:
            record = self.records.get(sample_id)
            if record is not None:
                yield record


@dataclass(repr=False)
class AirborneEMDataset(CoreObject):
    """Format-neutral collection of airborne EM flight lines.

    Parameters
    ----------
    name : str
        Survey/dataset name. Stripped and required to be non-empty.
    lines : dict of str to AirborneEMLine, optional
        Flight lines keyed by ``line_id``. The mapping key must equal
        ``line.line_id`` for every entry.
    survey : SurveyMeta, optional
        Survey-level metadata.
    instrument : InstrumentMeta, optional
        System/instrument metadata.
    method : str, default "AEM"
        Survey method label, upper-cased on construction (for example
        ``"AEM"``).
    attrs : dict, optional
        Free-form extension metadata.

    Raises
    ------
    ValueError
        If ``name`` or ``method`` is empty, or a line mapping key does
        not match ``line.line_id``.
    TypeError
        If ``survey``, ``instrument``, or an entry of ``lines`` has
        the wrong type.

    Notes
    -----
    The dataset is intentionally an organisational layer above
    :class:`~pycsamt.emtf.EMTF`. It does not define a MobileMT, ZTEM,
    or AFMAG file schema, and it inherits
    :class:`~pycsamt.core.base.CoreObject` rather than
    :class:`~pycsamt.core.base.MTBase`: aggregating flight lines is
    not itself electromagnetic arithmetic, so this class should not
    carry :class:`~pycsamt.core.base.MTBase`'s numeric EM utilities
    (those belong to the :class:`~pycsamt.emtf.EMTF`/
    :class:`~pycsamt.emtf.TransferFunction` objects it holds).
    Technology adapters populate this object rather than introducing
    separate transfer-function mathematics.

    Examples
    --------
    >>> from pycsamt.airborne import AirborneEMDataset
    >>> dataset = AirborneEMDataset(name="survey-001")
    >>> dataset.method, dataset.n_lines
    ('AEM', 0)
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
        """Normalize identifier/method fields and re-attach lines."""
        self.name = clean_identifier(self.name, name="name")
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
    ) -> AirborneEMDataset:
        """Attach one flight line.

        Parameters
        ----------
        line : AirborneEMLine
            Flight line to attach.
        replace : bool, default False
            Whether to overwrite an existing line with the same
            ``line_id`` instead of raising.

        Returns
        -------
        AirborneEMDataset
            ``self``, to support call chaining.

        Raises
        ------
        TypeError
            If ``line`` is not an :class:`AirborneEMLine`.
        ValueError
            If a line already exists for that ``line_id`` and
            ``replace`` is ``False``.
        """
        if not isinstance(line, AirborneEMLine):
            raise TypeError("line must be an AirborneEMLine")
        if line.line_id in self.lines and not replace:
            raise ValueError(f"airborne line already exists: {line.line_id}")
        self.lines[line.line_id] = line
        return self

    def get_line(self, line_id: str) -> AirborneEMLine | None:
        """Return a line by identifier, or ``None`` when absent."""
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

    def emtf_records(self) -> dict[tuple[str, str], EMTF]:
        """Return all non-empty EMTF records keyed by line/sample ID.

        Returns
        -------
        dict of (str, str) to EMTF
            Mapping from ``(line_id, sample_id)`` to the attached
            :class:`~pycsamt.emtf.EMTF`. Records with no EMTF payload
            are omitted rather than represented with a placeholder.
        """
        out: dict[tuple[str, str], EMTF] = {}
        for line_id, record in self.iter_records():
            if record.emtf is not None:
                out[(line_id, record.sample_id)] = record.emtf
        return out

    def inspect(self):
        """Return the common airborne inspection summary lazily.

        Returns
        -------
        AirborneInspection
            Compact scientific inventory; see
            :func:`pycsamt.airborne.qc.inspect_airborne`.

        Notes
        -----
        The import is deferred to avoid a hard import-time dependency
        between :mod:`pycsamt.airborne.base` and
        :mod:`pycsamt.airborne.qc`, which itself imports this module.
        """
        from .qc import inspect_airborne

        return inspect_airborne(self)

    def qc(self):
        """Return the common structural airborne QC report lazily.

        Returns
        -------
        AirborneQCReport
            Structural/metadata completeness report; see
            :func:`pycsamt.airborne.qc.assess_airborne_qc`.
        """
        from .qc import assess_airborne_qc

        return assess_airborne_qc(self)
