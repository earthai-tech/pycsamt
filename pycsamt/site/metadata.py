"""Declarative, auditable editing of EDI station metadata."""

from __future__ import annotations

import copy
import inspect
import math
import warnings
from collections.abc import Callable, Mapping, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Literal

import pandas as pd

from ..api.view import maybe_wrap_frame
from ..seg.heads import Info
from .base import Site, Sites
from .export import write_sites
from .utils import _ensure_head, get_coords, set_coords, station_name

__all__ = [
    "MetadataChange",
    "SiteMetadataEditor",
    "rename_sites",
    "update_metadata",
    "update_metadata_all",
]

MissingPolicy = Literal["raise", "warn", "ignore"]
ErrorPolicy = Literal["raise", "warn", "ignore"]
_SPEC_KEYS = {
    "name",
    "station",
    "lat",
    "lon",
    "long",
    "elev",
    "coords",
    "head",
    "info",
    "sections",
    "set",
    "unset",
    "transform",
}
_STATION_COLUMNS = ("station", "name", "site", "dataid", "id")


@dataclass(frozen=True)
class MetadataChange:
    r"""Describe the metadata changes attempted for one station.

    Parameters
    ----------
    index : int
        Zero-based position of the station in the input collection.
    old_name : str
        Station identity before editing.
    new_name : str
        Station identity after editing. For a failed operation this remains
        equal to ``old_name``.
    changed_fields : tuple of str
        Canonical paths of fields whose values changed.
    status : {'updated', 'unchanged', 'error'}, default='updated'
        Outcome of the station-level operation.
    error : str or None, default=None
        Error message captured when ``status='error'``.
    requested_fields : tuple of str, default=()
        Canonical paths requested by the metadata specification.
    before, after : mapping or None, default=None
        Audit snapshots surrounding the operation.

    Attributes
    ----------
    index : int
    old_name, new_name : str
    changed_fields, requested_fields : tuple of str
    status : str
    error : str or None
    before, after : mapping or None

    Notes
    -----
    Instances are immutable. Use :meth:`to_dict` when serializing an audit or
    constructing a tabular report.

    Examples
    --------
    >>> from pycsamt.site import MetadataChange
    >>> change = MetadataChange(
    ...     0, "18-012A", "L01_012", ("name", "head.project")
    ... )
    >>> change.status
    'updated'
    >>> change.to_dict()["changed_fields"]
    ['name', 'head.project']

    See Also
    --------
    SiteMetadataEditor.audit : Return all station records as a DataFrame.
    """

    index: int
    old_name: str
    new_name: str
    changed_fields: tuple[str, ...]
    status: str = "updated"
    error: str | None = None
    requested_fields: tuple[str, ...] = ()
    before: Mapping[str, Any] | None = None
    after: Mapping[str, Any] | None = None

    def to_dict(self) -> dict[str, Any]:
        r"""Return the record as a serialization-friendly dictionary.

        Returns
        -------
        dict of str to Any
            Dataclass fields with tuple-valued field lists converted to
            ordinary lists.

        Examples
        --------
        >>> from pycsamt.site import MetadataChange
        >>> record = MetadataChange(0, "A01", "B01", ("name",))
        >>> record.to_dict()["changed_fields"]
        ['name']

        See Also
        --------
        SiteMetadataEditor.audit : Build a DataFrame from change records.
        """
        data = asdict(self)
        data["changed_fields"] = list(self.changed_fields)
        data["requested_fields"] = list(self.requested_fields)
        return data


class SiteMetadataEditor:
    r"""Apply declarative, validated, and auditable EDI metadata changes.

    Parameters
    ----------
    updates : mapping, sequence, callable, pandas.DataFrame, or path-like
        Metadata source. A mapping may be keyed by current station identity or
        may be one specification applied to every station. A sequence is
        aligned with input order. A callable receives an EDI object and,
        optionally, its zero-based index. DataFrames and CSV files require a
        station column named ``station``, ``name``, ``site``, ``dataid``, or
        ``id``.
    missing : {'raise', 'warn', 'ignore'}, default='raise'
        Policy for source keys that match no input station.
    allow_duplicates : bool, default=False
        Permit duplicate final station identities. Keeping the default avoids
        ambiguous selection and export filenames.
    on_error : {'raise', 'warn', 'ignore'}, default='raise'
        Station-level failure policy. ``raise`` preserves batch atomicity;
        ``warn`` and ``ignore`` retain failed stations unchanged and commit
        successful stations.
    validate_coordinates : bool, default=True
        Validate finite latitude, longitude, and elevation values and enforce
        geographic latitude/longitude bounds.
    allow_empty_names : bool, default=False
        Permit an empty final station identity.
    validators : sequence of callable or None, default=None
        Additional validators called after each staged update. A validator
        receives the staged EDI and, optionally, its index. Returning ``False``
        rejects the station; raising an exception records that exception.

    Attributes
    ----------
    updates : Any
        Original metadata source.
    missing, on_error : str
        Configured unmatched-key and station-error policies.
    allow_duplicates, validate_coordinates, allow_empty_names : bool
        Validation switches.
    validators : tuple of callable
        Custom validators in execution order.
    records_ : list of MetadataChange
        Audit records from the latest :meth:`apply` or :meth:`plan` call.
    output_paths_ : list of pathlib.Path
        Paths written by the latest :meth:`apply_and_write` call.

    Notes
    -----
    A station specification recognizes ``name``/``station``, ``lat``, ``lon``,
    ``elev``, ``coords``, ``head``, ``info``, ``sections``, ``set``, ``unset``,
    and ``transform``. Generic paths default to ``HEAD``. Explicit path forms
    are ``head.<field>``, ``info.<field>``, ``edi.<field>``, and
    ``section.<name>.<field>``.

    All changes are staged on deep copies. With ``on_error='raise'``, an
    in-place batch is committed only after every station and final identity
    constraint passes validation.

    Examples
    --------
    Rename stations and update acquisition metadata:

    >>> from pycsamt.site import SiteMetadataEditor
    >>> editor = SiteMetadataEditor(
    ...     {
    ...         "18-012A": {
    ...             "name": "L01_012",
    ...             "coords": (5.25, -3.75, 120.0),
    ...             "head": {"project": "LINE_01"},
    ...         }
    ...     }
    ... )
    >>> # updated = editor.apply(sites)
    >>> # editor.audit()[["old_name", "new_name", "status"]]

    Generic actions can address nested fields:

    >>> editor = SiteMetadataEditor(
    ...     {
    ...         "18-012A": {
    ...             "set": {"info.processingtag": "reviewed"},
    ...             "transform": {"head.elev": lambda value: value + 1.5},
    ...             "unset": ["head.county"],
    ...         }
    ...     }
    ... )

    See Also
    --------
    update_metadata : Update one site or EDI-like object.
    update_metadata_all : Update a station collection.
    rename_sites : Rename from a mapping, sequence, or callable.
    pycsamt.site.export.write_sites : Export edited stations separately.
    """

    def __init__(
        self,
        updates: Any,
        *,
        missing: MissingPolicy = "raise",
        allow_duplicates: bool = False,
        on_error: ErrorPolicy = "raise",
        validate_coordinates: bool = True,
        allow_empty_names: bool = False,
        validators: Sequence[Callable[..., Any]] | None = None,
    ) -> None:
        if missing not in {"raise", "warn", "ignore"}:
            raise ValueError("missing must be 'raise', 'warn', or 'ignore'")
        if on_error not in {"raise", "warn", "ignore"}:
            raise ValueError("on_error must be 'raise', 'warn', or 'ignore'")
        self.updates = updates
        self.missing = missing
        self.allow_duplicates = bool(allow_duplicates)
        self.on_error = on_error
        self.validate_coordinates = bool(validate_coordinates)
        self.allow_empty_names = bool(allow_empty_names)
        self.validators = tuple(validators or ())
        self.records_: list[MetadataChange] = []
        self.output_paths_: list[Path] = []

    def apply(self, source: Any, *, inplace: bool = False) -> Any:
        r"""Apply the configured metadata updates.

        Parameters
        ----------
        source : Site, Sites, EDI-like object, or iterable of EDI-like objects
            Station data to edit.
        inplace : bool, default=False
            Commit successful staged objects back into ``source``. With the
            default, return independently copied objects.

        Returns
        -------
        Site, Sites, or EDI-like object
            One-site inputs retain their logical type. Collection inputs return
            :class:`~pycsamt.site.base.Sites` unless the input is already a
            ``Sites`` object edited in place.

        Raises
        ------
        KeyError
            If metadata keys are unmatched and ``missing='raise'``.
        ValueError
            If names, coordinates, actions, or validators fail validation.
        TypeError
            If the source cannot be staged safely or an action has an invalid
            type.

        Notes
        -----
        All edits are first performed on private copies. With ``inplace=True``
        and ``on_error='raise'``, the original object is updated only after the
        complete batch passes validation. ``warn`` and ``ignore`` deliberately
        commit successful stations while retaining failed stations unchanged.

        Examples
        --------
        >>> from pycsamt.site import SiteMetadataEditor
        >>> editor = SiteMetadataEditor({"A01": {"name": "L01_001"}})
        >>> # renamed = editor.apply(sites)
        >>> # renamed["L01_001"].name

        See Also
        --------
        plan : Preview and validate without modifying the source.
        apply_and_write : Apply and export in one operation.
        audit : Return records from the latest operation.
        """
        self.output_paths_ = []
        items, kind = _unpack(source)
        if not items:
            self.records_ = []
            return _repack(source, [], kind, inplace=inplace)

        original_names = [station_name(ed) for ed in items]
        specs, unused = _resolve_specs(self.updates, items, original_names)
        specs = [
            _materialize_identity(spec, ed, index, old)
            for index, (spec, ed, old) in enumerate(zip(specs, items, original_names))
        ]
        if unused:
            self._handle_missing(
                "metadata keys did not match any station: " + ", ".join(sorted(unused))
            )

        desired = [
            _desired_name(spec, ed, i, old)
            for i, (spec, ed, old) in enumerate(zip(specs, items, original_names))
        ]
        if not self.allow_empty_names:
            empty = [str(i) for i, name in enumerate(desired) if not name.strip()]
            if empty:
                raise ValueError(
                    "metadata update would create empty station names at indices: "
                    + ", ".join(empty)
                )
        if not self.allow_duplicates:
            duplicates = _duplicates(desired)
            if duplicates:
                raise ValueError(
                    "metadata update would create duplicate station names: "
                    + ", ".join(sorted(duplicates))
                )

        # Always stage on clones; this gives both copy mode and in-place mode
        # the same atomic failure semantics.
        targets = [_clone(ed) for ed in items]
        records: list[MetadataChange] = []
        for index, (ed, old, spec) in enumerate(zip(targets, original_names, specs)):
            if spec is None:
                snapshot = _snapshot(ed)
                records.append(
                    MetadataChange(
                        index,
                        old,
                        old,
                        (),
                        status="unchanged",
                        before=snapshot,
                        after=snapshot,
                    )
                )
                continue
            requested = tuple(_requested_fields(spec))
            before = _snapshot(ed, spec)
            try:
                changed = _apply_spec(
                    ed,
                    spec,
                    index=index,
                    validate_coordinates=self.validate_coordinates,
                )
                for validator in self.validators:
                    verdict = _call_validator(validator, ed, index)
                    if verdict is False:
                        raise ValueError(f"metadata validator rejected station {old!r}")
                new = station_name(ed)
                after = _snapshot(ed, spec)
                records.append(
                    MetadataChange(
                        index,
                        old,
                        new,
                        tuple(changed),
                        status="updated" if changed else "unchanged",
                        requested_fields=requested,
                        before=before,
                        after=after,
                    )
                )
            except Exception as exc:
                records.append(
                    MetadataChange(
                        index,
                        old,
                        old,
                        (),
                        status="error",
                        error=str(exc),
                        requested_fields=requested,
                        before=before,
                        after=before,
                    )
                )
                # Discard this station's failed working copy.
                targets[index] = _clone(items[index])
                if self.on_error == "raise":
                    self.records_ = records
                    raise
                if self.on_error == "warn":
                    warnings.warn(str(exc), UserWarning, stacklevel=2)

        self.records_ = records
        final_names = [station_name(ed) for ed in targets]
        if not self.allow_empty_names and any(not name.strip() for name in final_names):
            raise ValueError("metadata actions produced an empty station name")
        if not self.allow_duplicates:
            duplicates = _duplicates(final_names)
            if duplicates:
                raise ValueError(
                    "metadata actions produced duplicate station names: "
                    + ", ".join(sorted(duplicates))
                )

        return _repack(source, targets, kind, inplace=inplace)

    def plan(self, source: Any, *, api: bool | None = False) -> Any:
        r"""Validate and preview changes without modifying the source.

        Parameters
        ----------
        source : Site, Sites, EDI-like object, or iterable of EDI-like objects
            Station data used to evaluate the configured changes.
        api : bool or None, default=False
            Passed to the API-view wrapper. ``False`` returns a pandas
            DataFrame, ``True`` forces an API frame, and ``None`` defers to the
            global API-view configuration.

        Returns
        -------
        pandas.DataFrame or APIFrame
            Audit preview with one row per input station.

        Raises
        ------
        KeyError, ValueError, TypeError
            Propagated from staged resolution and validation, according to the
            configured policies.

        Notes
        -----
        Actions run only on staged copies, so callable transformations and
        validators are evaluated realistically. A later :meth:`apply` invokes
        callables again; stateful callables should therefore be avoided.

        Examples
        --------
        >>> from pycsamt.site import SiteMetadataEditor
        >>> editor = SiteMetadataEditor({"A01": {"elev": 125.0}})
        >>> # preview = editor.plan(sites)
        >>> # preview[["old_name", "changed_fields", "status"]]

        See Also
        --------
        apply : Apply the configured updates.
        audit : Return the most recently generated records.
        """
        self.apply(source, inplace=False)
        return self.audit(api=api)

    def apply_and_write(
        self,
        source: Any,
        outdir: str | Path,
        *,
        inplace: bool = False,
        template: str = "{station}.edi",
        exist_ok: bool = False,
        manifest_csv: str | Path | None = None,
    ) -> Any:
        r"""Apply metadata changes and export the resulting stations.

        Parameters
        ----------
        source : Site, Sites, EDI-like object, or iterable of EDI-like objects
            Station data to edit and export.
        outdir : path-like
            Destination directory.
        inplace : bool, default=False
            Commit staged metadata changes back into ``source``.
        template : str, default='{station}.edi'
            Export filename template accepted by
            :func:`pycsamt.site.export.write_sites`.
        exist_ok : bool, default=False
            Permit destinations that already exist.
        manifest_csv : path-like or None, default=None
            Optional manifest CSV destination.

        Returns
        -------
        Site, Sites, or EDI-like object
            Edited result. Written paths are available in ``output_paths_``.

        Raises
        ------
        KeyError, ValueError, TypeError
            Propagated from metadata resolution and validation.
        FileExistsError
            If an export destination exists and ``exist_ok=False``.
        RuntimeError
            If an EDI backend cannot write a station.

        Notes
        -----
        Editing and persistence remain separate internally: :meth:`apply` is
        completed before :func:`pycsamt.site.export.write_sites` is called.

        Examples
        --------
        >>> from pycsamt.site import SiteMetadataEditor
        >>> editor = SiteMetadataEditor({"A01": {"name": "L01_001"}})
        >>> # result = editor.apply_and_write(sites, "renamed_edi")
        >>> # [path.name for path in editor.output_paths_]

        See Also
        --------
        apply : Apply without writing files.
        pycsamt.site.export.write_sites : Export an existing collection.
        """
        result = self.apply(source, inplace=inplace)
        self.output_paths_ = write_sites(
            result,
            outdir,
            template=template,
            exist_ok=exist_ok,
            manifest_csv=manifest_csv,
        )
        return result

    def audit(self, *, api: bool | None = False) -> Any:
        r"""Return records from the latest operation as a table.

        Parameters
        ----------
        api : bool or None, default=False
            ``False`` returns a pandas DataFrame, ``True`` forces an API frame,
            and ``None`` defers to the global API-view configuration.

        Returns
        -------
        pandas.DataFrame or APIFrame
            Columns correspond to :class:`MetadataChange` fields. Before any
            operation, an empty table with the stable audit schema is returned.

        Examples
        --------
        >>> from pycsamt.site import SiteMetadataEditor
        >>> editor = SiteMetadataEditor({"A01": {"name": "B01"}})
        >>> list(editor.audit(api=False).columns[:4])
        ['index', 'old_name', 'new_name', 'changed_fields']

        See Also
        --------
        MetadataChange : Station-level audit record.
        plan : Populate the audit through a non-mutating preview.
        apply : Populate the audit while applying updates.
        """
        columns = [field.name for field in MetadataChange.__dataclass_fields__.values()]
        frame = pd.DataFrame(
            [record.to_dict() for record in self.records_], columns=columns
        )
        return maybe_wrap_frame(
            frame,
            api=api,
            name="site_metadata_audit",
            kind="metadata",
            source="SiteMetadataEditor",
        )

    def _handle_missing(self, message: str) -> None:
        if self.missing == "raise":
            raise KeyError(message)
        if self.missing == "warn":
            warnings.warn(message, UserWarning, stacklevel=3)


def update_metadata(
    site: Any,
    update: Mapping[str, Any],
    *,
    inplace: bool = False,
    validate_coordinates: bool = True,
    validators: Sequence[Callable[..., Any]] | None = None,
) -> Any:
    r"""Update metadata for one site or EDI-like object.

    Parameters
    ----------
    site : Site or EDI-like object
        Object to update.
    update : mapping
        One metadata specification. Supported keys are ``name``, ``station``,
        ``lat``, ``lon``, ``long``, ``elev``, ``coords``, ``head``, ``info``,
        ``sections``, ``set``, ``unset``, and ``transform``.
    inplace : bool, default=False
        Commit the staged state into ``site`` rather than returning an
        independent copy.
    validate_coordinates : bool, default=True
        Enforce finite geographic coordinate values and valid latitude and
        longitude ranges.
    validators : sequence of callable or None, default=None
        Additional staged-object validators.

    Returns
    -------
    Site or EDI-like object
        Updated object with the same logical single-site form as the input.

    Raises
    ------
    ValueError
        If a field, coordinate, station identity, or validator is invalid.
    TypeError
        If the update specification or input cannot be handled safely.

    Examples
    --------
    >>> from pycsamt.site import update_metadata
    >>> update = {
    ...     "name": "L01_012",
    ...     "coords": (5.25, -3.75, 120.0),
    ...     "info": {"processingtag": "reviewed"},
    ... }
    >>> # reviewed = update_metadata(site, update)

    See Also
    --------
    update_metadata_all : Apply station-specific updates to a collection.
    SiteMetadataEditor : Configure validation, planning, and audit behavior.
    rename_sites : Rename one or many stations using a compact interface.
    """
    return SiteMetadataEditor(
        dict(update),
        validate_coordinates=validate_coordinates,
        validators=validators,
    ).apply(site, inplace=inplace)


def update_metadata_all(
    sites: Any,
    updates: Any,
    *,
    inplace: bool = False,
    missing: MissingPolicy = "raise",
    allow_duplicates: bool = False,
    on_error: ErrorPolicy = "raise",
    validate_coordinates: bool = True,
    allow_empty_names: bool = False,
    validators: Sequence[Callable[..., Any]] | None = None,
) -> Any:
    r"""Apply metadata specifications to a station collection.

    Parameters
    ----------
    sites : Sites, iterable of EDI-like objects, Site, or EDI-like object
        Input station data.
    updates : mapping, sequence, callable, pandas.DataFrame, or path-like
        Metadata source accepted by :class:`SiteMetadataEditor`.
    inplace : bool, default=False
        Commit staged objects back into the supplied input.
    missing : {'raise', 'warn', 'ignore'}, default='raise'
        Policy for update keys that match no station.
    allow_duplicates : bool, default=False
        Permit duplicate final station identities.
    on_error : {'raise', 'warn', 'ignore'}, default='raise'
        Station-level failure policy.
    validate_coordinates : bool, default=True
        Validate geographic coordinates before committing.
    allow_empty_names : bool, default=False
        Permit empty final station identities.
    validators : sequence of callable or None, default=None
        Additional validators applied to each staged station.

    Returns
    -------
    Sites, Site, or EDI-like object
        Updated data. Collection-like inputs normally return
        :class:`~pycsamt.site.base.Sites`.

    Raises
    ------
    KeyError
        If station-keyed metadata contains unmatched keys and
        ``missing='raise'``.
    ValueError
        If the batch violates metadata or identity constraints.
    TypeError
        If the source, metadata source, or action is unsupported.

    Examples
    --------
    Use an explicit station mapping:

    >>> from pycsamt.site import update_metadata_all
    >>> updates = {
    ...     "A01": {"name": "L01_001", "head": {"project": "L01"}},
    ...     "A02": {"name": "L01_002", "elev": 121.0},
    ... }
    >>> # updated = update_metadata_all(sites, updates)

    A DataFrame or CSV review table can use columns such as ``station``,
    ``new_name``, ``latitude``, and ``head.project``.

    See Also
    --------
    update_metadata : Update one site.
    SiteMetadataEditor.apply : Apply with a reusable configured editor.
    SiteMetadataEditor.plan : Preview a batch before committing.
    rename_sites : Rename a collection without a full metadata specification.
    """
    return SiteMetadataEditor(
        updates,
        missing=missing,
        allow_duplicates=allow_duplicates,
        on_error=on_error,
        validate_coordinates=validate_coordinates,
        allow_empty_names=allow_empty_names,
        validators=validators,
    ).apply(sites, inplace=inplace)


def rename_sites(
    sites: Any,
    names: Mapping[str, str] | Sequence[str] | Callable[..., str],
    *,
    inplace: bool = False,
    missing: MissingPolicy = "raise",
    allow_duplicates: bool = False,
    allow_empty_names: bool = False,
) -> Any:
    r"""Rename stations from a mapping, aligned sequence, or callable.

    Parameters
    ----------
    sites : Sites, iterable of EDI-like objects, Site, or EDI-like object
        Stations to rename.
    names : mapping, sequence of str, or callable
        A mapping relates current names to new names. A sequence is aligned
        with input order. A callable receives an EDI object and, optionally,
        its zero-based index, and returns the new name.
    inplace : bool, default=False
        Commit synchronized identities back into the input.
    missing : {'raise', 'warn', 'ignore'}, default='raise'
        Policy for mapping keys that match no station.
    allow_duplicates : bool, default=False
        Permit duplicate final identities.
    allow_empty_names : bool, default=False
        Permit empty final identities.

    Returns
    -------
    Sites, Site, or EDI-like object
        Renamed station data.

    Raises
    ------
    KeyError
        If a mapping key is unmatched and ``missing='raise'``.
    ValueError
        If names are duplicate or empty under the configured policy, or a
        sequence length differs from the number of stations.
    TypeError
        If ``names`` or the station source is unsupported.

    Notes
    -----
    Renaming synchronizes object-level identity, common ``HEAD`` aliases, and
    linked ``SECTID`` values. It does not rename an existing source file;
    exporting with ``template='{station}.edi'`` uses the new identity.

    Examples
    --------
    >>> from pycsamt.site import rename_sites
    >>> mapping = {"18-012A": "L01_012", "18-013A": "L01_013"}
    >>> # renamed = rename_sites(sites, mapping)

    Generate names from input order:

    >>> # renamed = rename_sites(
    >>> #     sites, lambda _edi, index: f"L22_{index + 1:03d}"
    >>> # )

    See Also
    --------
    update_metadata : Update one station and its metadata.
    update_metadata_all : Apply richer station-specific specifications.
    SiteMetadataEditor : Reusable editor with planning and audit records.
    pycsamt.site.export.write_sites : Export using updated station names.
    """
    if isinstance(names, Mapping):
        updates = {str(old): {"name": str(new)} for old, new in names.items()}
    elif callable(names):

        def updates(site: Any, index: int) -> dict[str, str]:
            return {"name": str(_call(names, site, index))}

    else:
        updates = [{"name": str(name)} for name in names]
    return update_metadata_all(
        sites,
        updates,
        inplace=inplace,
        missing=missing,
        allow_duplicates=allow_duplicates,
        allow_empty_names=allow_empty_names,
    )


def _resolve_specs(
    updates: Any, items: list[Any], names: list[str]
) -> tuple[list[dict[str, Any] | None], set[str]]:
    updates = _coerce_update_source(updates)
    if callable(updates):
        return [
            _normalize_spec(_call(updates, ed, i)) for i, ed in enumerate(items)
        ], set()

    if isinstance(updates, Mapping):
        if _is_spec(updates):
            return [dict(updates) for _ in items], set()
        lookup: dict[str, tuple[str, Any]] = {}
        for key, value in updates.items():
            original = str(key)
            folded = original.casefold()
            if folded in lookup:
                raise ValueError(
                    f"metadata mapping contains duplicate case-insensitive key: {original!r}"
                )
            lookup[folded] = (original, value)
        used: set[str] = set()
        specs: list[dict[str, Any] | None] = []
        for name in names:
            key = name.casefold()
            if key in lookup:
                used.add(key)
                specs.append(_normalize_spec(lookup[key][1]))
            else:
                specs.append(None)
        return specs, {lookup[key][0] for key in set(lookup) - used}

    if isinstance(updates, Sequence) and not isinstance(updates, (str, bytes)):
        if len(updates) != len(items):
            raise ValueError(
                f"expected {len(items)} metadata specifications, got {len(updates)}"
            )
        return [_normalize_spec(value) for value in updates], set()
    raise TypeError("updates must be a mapping, sequence, or callable")


def _coerce_update_source(updates: Any) -> Any:
    """Normalize CSV/DataFrame-like update tables to station mappings."""
    if isinstance(updates, (str, Path)):
        path = Path(updates)
        if not path.exists():
            raise FileNotFoundError(path)
        updates = pd.read_csv(path)

    if hasattr(updates, "to_dict") and hasattr(updates, "columns"):
        columns = {str(column).casefold(): column for column in updates.columns}
        station_column = next(
            (columns[name] for name in _STATION_COLUMNS if name in columns),
            None,
        )
        if station_column is None:
            raise ValueError(
                "metadata table needs a station column; accepted names are: "
                + ", ".join(_STATION_COLUMNS)
            )
        result: dict[str, dict[str, Any]] = {}
        seen: set[str] = set()
        for row in updates.to_dict(orient="records"):
            station = row.pop(station_column)
            if _is_missing(station):
                raise ValueError("metadata table contains an empty station key")
            spec: dict[str, Any] = {}
            for key, value in row.items():
                if _is_missing(value):
                    continue
                _assign_table_field(spec, str(key), value)
            name = str(station)
            folded_name = name.casefold()
            if folded_name in seen:
                raise ValueError(
                    f"metadata table contains duplicate station row: {name!r}"
                )
            seen.add(folded_name)
            result[name] = spec
        return result
    return updates


def _assign_table_field(spec: dict[str, Any], key: str, value: Any) -> None:
    normalized = key.strip()
    folded = normalized.casefold()
    aliases = {
        "new_name": "name",
        "latitude": "lat",
        "longitude": "lon",
        "elevation": "elev",
    }
    if folded in aliases:
        spec[aliases[folded]] = value
    elif folded in _SPEC_KEYS:
        spec[folded] = value
    elif "." in normalized and folded.split(".", 1)[0] in {"head", "info"}:
        section, path = normalized.split(".", 1)
        spec.setdefault(section.casefold(), {})[path] = value
    elif folded.startswith("section.") and normalized.count(".") >= 2:
        _, section, path = normalized.split(".", 2)
        spec.setdefault("sections", {}).setdefault(section, {})[path] = value
    else:
        spec.setdefault("set", {})[normalized] = value


def _is_missing(value: Any) -> bool:
    if value is None:
        return True
    try:
        result = math.isnan(value)
        return bool(result)
    except (TypeError, ValueError):
        return type(value).__name__ in {"NAType", "NaTType"}


def _values_equal(left: Any, right: Any) -> bool:
    if _is_missing(left) and _is_missing(right):
        return True
    try:
        result = left == right
        return bool(result) if not hasattr(result, "all") else bool(result.all())
    except Exception:
        return False


def _normalize_spec(value: Any) -> dict[str, Any]:
    if isinstance(value, str):
        return {"name": value}
    if not isinstance(value, Mapping):
        raise TypeError("each metadata specification must be a mapping or name")
    unknown = set(value) - _SPEC_KEYS
    if unknown:
        raise ValueError("unknown metadata fields: " + ", ".join(sorted(unknown)))
    return dict(value)


def _is_spec(value: Mapping[Any, Any]) -> bool:
    if not value or not set(value).issubset(_SPEC_KEYS):
        return False
    # A real station may itself be named "station" or "name". A nested
    # mapping at either key is therefore a keyed batch, not one specification.
    return not any(
        key in value and isinstance(value[key], Mapping) for key in ("station", "name")
    )


def _apply_spec(
    ed: Any,
    spec: Mapping[str, Any],
    *,
    index: int,
    validate_coordinates: bool,
) -> list[str]:
    changed: list[str] = []
    old_name = station_name(ed)
    new_name = spec.get("name", spec.get("station"))
    new_name = _resolve_value(new_name, ed, index, old_name)
    if new_name is not None and str(new_name) != old_name:
        _set_identity(ed, str(new_name))
        changed.append("name")

    before = get_coords(ed)
    coords = spec.get("coords")
    coord_values = _normalize_coords(coords, ed, index) if coords is not None else {}
    lat = spec.get("lat", coord_values.get("lat"))
    lon = spec.get("lon", spec.get("long", coord_values.get("lon")))
    elev = spec.get("elev", coord_values.get("elev"))
    lat = _resolve_value(lat, ed, index, before.lat)
    lon = _resolve_value(lon, ed, index, before.lon)
    elev = _resolve_value(elev, ed, index, before.elev)
    if any(value is not None for value in (lat, lon, elev)):
        if validate_coordinates:
            _validate_coords(lat=lat, lon=lon, elev=elev)
        set_coords(ed, lat=lat, lon=lon, elev=elev, inplace=True)
        after = get_coords(ed)
        for field in ("lat", "lon", "elev"):
            if not _values_equal(getattr(before, field), getattr(after, field)):
                changed.append(field)

    head_values = spec.get("head")
    if head_values is not None:
        if not isinstance(head_values, Mapping):
            raise TypeError("head metadata must be a mapping")
        head = _ensure_head(ed)
        for key, value in head_values.items():
            value = _resolve_value(value, ed, index, _get_path(head, str(key)))
            if _set_path(head, str(key), value):
                changed.append(f"head.{key}")

    info_values = spec.get("info")
    if info_values is not None:
        if not isinstance(info_values, Mapping):
            raise TypeError("info metadata must be a mapping")
        info = _ensure_section(ed, "info")
        for key, value in info_values.items():
            key = str(key)
            value = _resolve_value(value, ed, index, _get_info_path(info, key))
            updater = getattr(info, "update", None)
            if "." not in key and callable(updater):
                before = _info_value(info, key)
                updater(**{key: value})
                known = {
                    str(item).replace("_", "").casefold()
                    for item in getattr(info, "infokeys", ())
                }
                if value is not None and key.replace("_", "").casefold() not in known:
                    _upsert_info_text(info, key, value)
                elif value is None:
                    _remove_info_text(info, key)
                did_change = before != _info_value(info, key)
            else:
                did_change = _set_path(info, key, value)
            if did_change:
                changed.append(f"info.{key}")

    section_values = spec.get("sections")
    if section_values is not None:
        if not isinstance(section_values, Mapping):
            raise TypeError("sections action must map section names to fields")
        for section_name, values in section_values.items():
            if not isinstance(values, Mapping):
                raise TypeError(
                    f"metadata for section {section_name!r} must be a mapping"
                )
            section = _ensure_section(ed, str(section_name))
            for key, value in values.items():
                key = str(key)
                value = _resolve_value(value, ed, index, _get_path(section, key))
                if _set_path(section, key, value):
                    changed.append(f"section.{section_name}.{key}")

    set_values = spec.get("set")
    if set_values is not None:
        if not isinstance(set_values, Mapping):
            raise TypeError("set action must be a mapping of path to value")
        for path, value in set_values.items():
            path = str(path)
            current = _get_field(ed, path)
            value = _resolve_value(value, ed, index, current)
            if _set_field(ed, path, value):
                changed.append(_canonical_path(path))

    transforms = spec.get("transform")
    if transforms is not None:
        if not isinstance(transforms, Mapping):
            raise TypeError("transform action must be a mapping of path to callable")
        for path, transform in transforms.items():
            if not callable(transform):
                raise TypeError(f"transform for {path!r} must be callable")
            path = str(path)
            current = _get_field(ed, path)
            value = _call_transform(transform, current, ed, index)
            if _set_field(ed, path, value):
                changed.append(_canonical_path(path))

    unset = spec.get("unset")
    if unset is not None:
        paths = [unset] if isinstance(unset, str) else list(unset)
        for path in paths:
            path = str(path)
            if _unset_field(ed, path):
                changed.append(_canonical_path(path))
    return changed


def _desired_name(
    spec: Mapping[str, Any] | None,
    ed: Any,
    index: int,
    old: str,
) -> str:
    if spec is None:
        return old
    value = spec.get("name", spec.get("station", old))
    resolved = _resolve_value(value, ed, index, old)
    return old if resolved is None else str(resolved)


def _materialize_identity(
    spec: dict[str, Any] | None,
    ed: Any,
    index: int,
    old: str,
) -> dict[str, Any] | None:
    if spec is None:
        return None
    result = dict(spec)
    key = "name" if "name" in result else "station" if "station" in result else None
    if key is not None and callable(result[key]):
        result[key] = _resolve_value(result[key], ed, index, old)
    return result


def _normalize_coords(coords: Any, ed: Any, index: int) -> dict[str, Any]:
    coords = _resolve_value(coords, ed, index, get_coords(ed))
    if isinstance(coords, Mapping):
        return {
            "lat": coords.get("lat", coords.get("latitude")),
            "lon": coords.get("lon", coords.get("long", coords.get("longitude"))),
            "elev": coords.get("elev", coords.get("elevation")),
        }
    if isinstance(coords, Sequence) and not isinstance(coords, (str, bytes)):
        if len(coords) not in {2, 3}:
            raise ValueError("coords must contain (lat, lon) or (lat, lon, elev)")
        return {
            "lat": coords[0],
            "lon": coords[1],
            "elev": coords[2] if len(coords) == 3 else None,
        }
    raise TypeError("coords must be a mapping, sequence, or callable")


def _validate_coords(*, lat: Any, lon: Any, elev: Any) -> None:
    for label, value, lower, upper in (
        ("latitude", lat, -90.0, 90.0),
        ("longitude", lon, -180.0, 180.0),
    ):
        if value is None:
            continue
        number = float(value)
        if not math.isfinite(number) or not lower <= number <= upper:
            raise ValueError(f"{label} must be finite and in [{lower}, {upper}]")
    if elev is not None and not math.isfinite(float(elev)):
        raise ValueError("elevation must be finite when supplied")


def _resolve_value(value: Any, ed: Any, index: int, current: Any) -> Any:
    if not callable(value):
        return value
    return _call_value(value, current, ed, index)


def _call_value(fn: Callable[..., Any], current: Any, ed: Any, index: int) -> Any:
    count = _parameter_count(fn)
    if count >= 3:
        return fn(current, ed, index)
    if count == 2:
        return fn(current, ed)
    if count == 1:
        return fn(current)
    return fn()


def _call_transform(fn: Callable[..., Any], current: Any, ed: Any, index: int) -> Any:
    return _call_value(fn, current, ed, index)


def _canonical_path(path: str) -> str:
    prefix = path.split(".", 1)[0].casefold()
    return path if prefix in {"head", "info", "edi", "section"} else f"head.{path}"


def _split_field(ed: Any, path: str) -> tuple[Any, str, str]:
    text = path.strip()
    if not text:
        raise ValueError("metadata field path cannot be empty")
    if "." in text:
        prefix, remainder = text.split(".", 1)
        prefix = prefix.casefold()
    else:
        prefix, remainder = "head", text
    if prefix == "head":
        return _ensure_head(ed), remainder, "head"
    if prefix == "info":
        return _ensure_section(ed, "info"), remainder, "info"
    if prefix == "edi":
        return ed, remainder, "edi"
    if prefix == "section":
        if "." not in remainder:
            raise ValueError("section paths use 'section.<name>.<field>' syntax")
        section_name, field = remainder.split(".", 1)
        return (
            _ensure_section(ed, section_name),
            field,
            f"section.{section_name}",
        )
    # Unrecognized prefixes are ordinary nested HEAD paths.
    return _ensure_head(ed), text, "head"


def _get_field(ed: Any, path: str) -> Any:
    root, remainder, section = _split_field(ed, path)
    if section == "info":
        return _get_info_path(root, remainder)
    return _get_path(root, remainder)


def _set_field(ed: Any, path: str, value: Any) -> bool:
    normalized = _canonical_path(path).casefold()
    if normalized in {
        "edi.name",
        "edi.station",
        "head.dataid",
        "head.station",
        "head.sitename",
        "head.name",
    }:
        before = station_name(ed)
        _set_identity(ed, str(value))
        return before != str(value)
    root, remainder, section = _split_field(ed, path)
    if section == "info" and "." not in remainder:
        before = _info_value(root, remainder)
        updater = getattr(root, "update", None)
        if callable(updater):
            updater(**{remainder: value})
            known = {
                str(item).replace("_", "").casefold()
                for item in getattr(root, "infokeys", ())
            }
            if value is not None and remainder.replace("_", "").casefold() not in known:
                _upsert_info_text(root, remainder, value)
            return before != _info_value(root, remainder)
    return _set_path(root, remainder, value)


def _unset_field(ed: Any, path: str) -> bool:
    if _canonical_path(path).casefold() in {
        "edi.name",
        "edi.station",
        "head.dataid",
        "head.station",
        "head.sitename",
        "head.name",
    }:
        raise ValueError("station identity cannot be unset; rename it explicitly")
    root, remainder, section = _split_field(ed, path)
    current = (
        _get_info_path(root, remainder)
        if section == "info"
        else _get_path(root, remainder)
    )
    if current is None:
        return False
    if section == "info" and "." not in remainder:
        _remove_info_text(root, remainder)
    _set_field(ed, path, None)
    return True


def _get_path(root: Any, path: str) -> Any:
    target = root
    for part in path.split("."):
        if isinstance(target, Mapping):
            actual = _mapping_key(target, part)
            target = target.get(actual) if actual is not None else None
        else:
            actual = _attribute_name(target, part)
            target = getattr(target, actual, None) if actual is not None else None
        if target is None:
            break
    return target


def _get_info_path(info: Any, path: str) -> Any:
    return _info_value(info, path) if "." not in path else _get_path(info, path)


def _set_identity(ed: Any, name: str) -> None:
    try:
        ed.name = name
    except Exception:
        pass
    try:
        ed.station = name
    except Exception:
        pass
    head = _ensure_head(ed)
    for field in ("dataid", "station", "sitename", "name", "STATION"):
        try:
            setattr(head, field, name)
        except Exception:
            pass
    candidates = list(getattr(ed, "__dict__", {}).values())
    sections = getattr(ed, "sections", {})
    if isinstance(sections, Mapping):
        candidates.extend(sections.values())
    seen: set[int] = set()
    for section in candidates:
        if id(section) in seen:
            continue
        seen.add(id(section))
        if hasattr(section, "sectid"):
            try:
                section.sectid = name
            except Exception:
                pass


def _ensure_section(ed: Any, key: str) -> Any:
    section = None
    getter = getattr(ed, "get_section", None)
    if callable(getter):
        section = getter(key)
    if section is None:
        attribute = _attribute_name(ed, key)
        section = getattr(ed, attribute, None) if attribute is not None else None
    if section is None:
        if key == "info":
            section = Info()
        else:
            section = SimpleNamespace()
        setter = getattr(ed, "set_section", None)
        if callable(setter):
            setter(key, section)
        sections = getattr(ed, "sections", None)
        if isinstance(sections, dict):
            sections[str(key).lower()] = section
        try:
            setattr(ed, key.capitalize(), section)
        except Exception:
            pass
    return section


def _set_path(root: Any, path: str, value: Any) -> bool:
    parts = path.split(".")
    target = root
    for part in parts[:-1]:
        if isinstance(target, Mapping):
            actual = _mapping_key(target, part)
            child = target.get(actual) if actual is not None else None
        else:
            actual = _attribute_name(target, part)
            child = getattr(target, actual, None) if actual is not None else None
        if child is None:
            child = SimpleNamespace()
            if isinstance(target, dict):
                target[part] = child
            else:
                setattr(target, part, child)
        target = child
    leaf = parts[-1]
    if isinstance(target, Mapping):
        actual_leaf = _mapping_key(target, leaf) or leaf
        old = target.get(actual_leaf)
    else:
        actual_leaf = _attribute_name(target, leaf) or leaf
        old = getattr(target, actual_leaf, None)
    if isinstance(target, dict):
        target[actual_leaf] = value
    else:
        setattr(target, actual_leaf, value)
    return old != value


def _mapping_key(mapping: Mapping[Any, Any], requested: str) -> Any:
    normalized = requested.replace("_", "").casefold()
    for key in mapping:
        if str(key).replace("_", "").casefold() == normalized:
            return key
    return None


def _attribute_name(obj: Any, requested: str) -> str | None:
    if hasattr(obj, requested):
        return requested
    normalized = requested.replace("_", "").casefold()
    for candidate in dir(obj):
        if candidate.replace("_", "").casefold() == normalized:
            return candidate
    return None


def _info_value(info: Any, key: str) -> Any:
    as_dict = getattr(info, "as_dict", None)
    if callable(as_dict):
        values = as_dict()
        normalized = key.replace("_", "").casefold()
        for candidate, value in values.items():
            if str(candidate).replace("_", "").casefold() == normalized:
                return value
    actual = _attribute_name(info, key)
    return getattr(info, actual, None) if actual is not None else None


def _upsert_info_text(info: Any, key: str, value: Any) -> None:
    lines = list(getattr(info, "info_text", []) or [])
    prefix = key.strip().upper() + "="
    replacement = f"{prefix}{value}"
    for index, line in enumerate(lines):
        if str(line).strip().upper().startswith(prefix):
            lines[index] = replacement
            break
    else:
        lines.append(replacement)
    info.info_text = lines


def _remove_info_text(info: Any, key: str) -> None:
    prefix = key.strip().upper() + "="
    info.info_text = [
        line
        for line in list(getattr(info, "info_text", []) or [])
        if not str(line).strip().upper().startswith(prefix)
    ]


def _requested_fields(spec: Mapping[str, Any]) -> list[str]:
    fields: list[str] = []
    for key, value in spec.items():
        if key in {"head", "info", "set", "transform"} and isinstance(value, Mapping):
            fields.extend(
                f"{key}.{field}"
                if key in {"head", "info"}
                else _canonical_path(str(field))
                for field in value
            )
        elif key == "sections" and isinstance(value, Mapping):
            for section, section_fields in value.items():
                if isinstance(section_fields, Mapping):
                    fields.extend(
                        f"section.{section}.{field}" for field in section_fields
                    )
        elif key == "unset":
            values = [value] if isinstance(value, str) else list(value)
            fields.extend(_canonical_path(str(field)) for field in values)
        elif key == "coords":
            fields.extend(("lat", "lon", "elev"))
        else:
            fields.append(str(key))
    return list(dict.fromkeys(fields))


def _snapshot(ed: Any, spec: Mapping[str, Any] | None = None) -> dict[str, Any]:
    coords = get_coords(ed)
    snapshot: dict[str, Any] = {
        "name": station_name(ed),
        "lat": coords.lat,
        "lon": coords.lon,
        "elev": coords.elev,
    }
    if spec is None:
        return snapshot
    for path in _requested_fields(spec):
        if path in snapshot or path in {"station", "coords"}:
            continue
        try:
            snapshot[path] = _safe_audit_value(_get_field(ed, path))
        except Exception:
            snapshot[path] = None
    return snapshot


def _safe_audit_value(value: Any) -> Any:
    if isinstance(value, (str, int, float, bool, type(None))):
        return value
    if isinstance(value, Mapping):
        return {str(key): _safe_audit_value(item) for key, item in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return [_safe_audit_value(item) for item in value]
    return repr(value)


def _call(fn: Callable[..., Any], site: Any, index: int) -> Any:
    nparams = _parameter_count(fn, default=2)
    return fn(site, index) if nparams >= 2 else fn(site)


def _parameter_count(fn: Callable[..., Any], *, default: int = 1) -> int:
    try:
        return len(inspect.signature(fn).parameters)
    except (TypeError, ValueError):
        return default


def _call_validator(fn: Callable[..., Any], ed: Any, index: int) -> Any:
    count = _parameter_count(fn)
    return fn(ed, index) if count >= 2 else fn(ed)


def _clone(ed: Any) -> Any:
    try:
        return copy.deepcopy(ed)
    except Exception as exc:
        raise TypeError(
            f"cannot stage metadata safely for {type(ed).__name__}; "
            "the object must support deepcopy"
        ) from exc


def _unpack(source: Any) -> tuple[list[Any], str]:
    if isinstance(source, Site):
        return [source.edi], "site"
    if isinstance(source, Sites):
        return [site.edi for site in source], "sites"
    if hasattr(source, "get_section") or hasattr(source, "Z"):
        return [source], "edi"
    values = list(source)
    return [value.edi if isinstance(value, Site) else value for value in values], "many"


def _repack(source: Any, items: list[Any], kind: str, *, inplace: bool) -> Any:
    if kind == "site":
        if inplace:
            _replace_state(source.edi, items[0])
            return source
        return Site(items[0])
    if kind == "edi":
        if inplace:
            _replace_state(source, items[0])
            return source
        return items[0]
    if kind == "sites":
        if inplace:
            for site, staged in zip(source._items, items):
                _replace_state(site.edi, staged)
            return source
        return Sites(items)
    if inplace:
        if not isinstance(source, Sequence):
            raise TypeError(
                "inplace=True requires a reusable sequence, Site, Sites, "
                "or EDI object; one-shot iterables cannot be committed safely"
            )
        originals = list(source)
        raw = [item.edi if isinstance(item, Site) else item for item in originals]
        for original, staged in zip(raw, items):
            _replace_state(original, staged)
        return Sites(raw)
    return Sites(items)


def _replace_state(target: Any, staged: Any) -> None:
    """Commit a staged object while preserving *target* identity."""
    target_state = getattr(target, "__dict__", None)
    staged_state = getattr(staged, "__dict__", None)
    if isinstance(target_state, dict) and isinstance(staged_state, dict):
        replacement = copy.deepcopy(staged_state)
        target_state.clear()
        target_state.update(replacement)
        return
    raise TypeError(
        f"cannot commit metadata in place for {type(target).__name__}; use inplace=False"
    )


def _duplicates(values: Sequence[str]) -> set[str]:
    seen: set[str] = set()
    duplicates: set[str] = set()
    for value in values:
        folded = value.casefold()
        if folded in seen:
            duplicates.add(value)
        seen.add(folded)
    return duplicates
