# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

import csv
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass, field
from os import PathLike
from os.path import commonpath
from pathlib import Path
from typing import (
    Any,
    Callable,
)

import numpy as np

from ..api.view import progress_enabled
from ..seg.edi import EDIFile
from .base import Sites, to_edis
from .edit import (
    _get_attr_any,
    _rotm,
    _set_attr_first,
    fill_missing,
    recompute_res_phase,
    rename,
    select_freq,
)
from .edit import (
    rotate as _rotate_site,
)
from .utils import maybe_copy, station_name

__all__ = [
    "EDIRecomputeRecord",
    "EDIRecomputeResult",
    "EDIRecomputer",
    "recompute_edi",
    "recompute_edis",
]


@dataclass
class EDIRecomputeRecord:
    r"""
    Per-station outcome for an EDI recomputation workflow.

    Attributes
    ----------
    source : pathlib.Path or None
        Source EDI path when known.
    output : pathlib.Path or None
        Written EDI path when ``write=True`` and the station
        was written successfully.
    line : str or None
        Line/group name inferred from the source directory.
    station : str
        Station name after recomputation and optional renaming.
    status : str
        ``"ok"`` for success or ``"failed"`` when processing
        continued after an error.
    message : str, default ""
        Optional diagnostic message.
    """

    source: Path | None
    output: Path | None
    line: str | None
    station: str
    status: str
    message: str = ""


@dataclass
class EDIRecomputeResult:
    r"""
    Result returned by :class:`EDIRecomputer`.

    Attributes
    ----------
    sites : pycsamt.site.base.Sites
        Recomputed EDI objects wrapped as sites.
    records : list of EDIRecomputeRecord
        Per-station processing and writing outcomes.
    output_root : pathlib.Path or None
        Root directory used for exported EDI files.
    """

    sites: Sites
    records: list[EDIRecomputeRecord]
    output_root: Path | None = None
    items: list[EDIFile] = field(default_factory=list, repr=False)

    @property
    def edis(self) -> list[EDIFile]:
        """Return recomputed raw EDI objects."""

        return list(self.items) if self.items else self.sites.to_edis()

    @property
    def paths(self) -> list[Path]:
        """Return successfully written output paths."""

        return [
            rec.output
            for rec in self.records
            if rec.output is not None and rec.status == "ok"
        ]

    @property
    def failed(self) -> list[EDIRecomputeRecord]:
        """Return failed records."""

        return [rec for rec in self.records if rec.status != "ok"]

    def to_manifest(self, path: str | Path) -> Path:
        """Write the records to a CSV manifest."""

        return _write_manifest(self.records, Path(path))


@dataclass
class _LineGroup:
    name: str | None
    sources: list[Path] | list[EDIFile]
    source_root: Path | None = None


@dataclass
class EDIRecomputer:
    r"""
    Recompute and rewrite EDI files using pyCSAMT conventions.

    This class is a workflow layer over the lower-level site
    helpers. It accepts EDI objects, EDI collections, files,
    line folders, or survey folders, applies a sequence of
    optional normalizations, and can write pyCSAMT-generated
    EDI files to a new output tree.

    Parameters
    ----------
    output_root : str or pathlib.Path, optional
        Root directory for recomputed EDI files. If omitted,
        a ``recomputed_edis`` directory is created next to a
        selected line folder, or inside a selected survey folder.
    output_name : str, default "recomputed_edis"
        Directory name used when ``output_root`` is not provided.
    preserve_line_dirs : bool, default True
        If ``True``, write each inferred line into its own
        subdirectory under ``output_root``. If ``False``, write
        all recomputed EDI files directly under ``output_root``.
    template : str, default "{station}.edi"
        Filename template. Available keys are ``station``,
        ``index``, ``line``, and ``source_stem``.
    overwrite : bool, default False
        Allow replacing existing output files.
    write : bool, default True
        If ``False``, only return recomputed in-memory objects.
    manifest_csv : bool or str or pathlib.Path, default True
        Write a CSV manifest. ``True`` writes
        ``<output_root>/manifest.csv``. A path writes there.
    rotate_angle : float, optional
        Rotation angle in degrees. If omitted, no rotation is
        applied.
    rotate_components : iterable of str, default ("Z", "Tip")
        Components to rotate. Accepts values such as ``"Z"``,
        ``"impedance"``, ``"Tip"``, and ``"tipper"``.
    fmin, fmax : float, optional
        Frequency range to keep before recomputation.
    keep_freq : iterable of int, optional
        Explicit frequency indices or mask passed to
        :func:`pycsamt.site.edit.select_freq`.
    fill_missing_values : {"zero", "nan"}, optional
        Missing-value policy passed to
        :func:`pycsamt.site.edit.fill_missing`.
    recompute_resphase : bool, default True
        Recompute apparent resistivity and phase from impedance.
    rename_policy : callable, optional
        Function mapping the current station name to a new name.
    datatype : str, optional
        EDI writer datatype override, for example ``"mt"`` or
        ``"emap"``.
    synthesize_spectra : bool, default False
        Ask the pyCSAMT EDI writer to synthesize spectra when
        possible and missing.
    recursive : bool, default True
        Recurse while discovering EDI files inside line folders.
    strict : bool, default False
        Raise on read/process/write errors instead of recording
        failed manifest rows.
    on_dup : {"replace", "keep"}, default "replace"
        Duplicate station policy during loading.
    copy : bool, default True
        Recompute copies of loaded EDI objects. Keep this enabled
        when the original objects should remain unchanged.
    progress : bool or {"auto"}, default False
        Show progress while recomputing.
    verbose : int, default 0
        Verbosity for loading, processing, and writing.
    """

    output_root: str | Path | None = None
    output_name: str = "recomputed_edis"
    preserve_line_dirs: bool = True
    template: str = "{station}.edi"
    overwrite: bool = False
    write: bool = True
    manifest_csv: bool | str | Path = True
    rotate_angle: float | None = None
    rotate_components: Iterable[str] = field(
        default_factory=lambda: ("Z", "Tip")
    )
    fmin: float | None = None
    fmax: float | None = None
    keep_freq: Iterable[int] | None = None
    fill_missing_values: str | None = None
    recompute_resphase: bool = True
    rename_policy: Callable[[str], str] | None = None
    datatype: str | None = None
    synthesize_spectra: bool = False
    recursive: bool = True
    strict: bool = False
    on_dup: str = "replace"
    copy: bool = True
    progress: bool | str = False
    verbose: int = 0
    progress_callback: Callable[[int, int, str, str, str], None] | None = (
        field(
            default=None,
            repr=False,
        )
    )

    def run(self, source: Any) -> EDIRecomputeResult:
        """Run the recomputation workflow."""

        groups = _discover_groups(
            source,
            recursive=self.recursive,
            strict=self.strict,
            on_dup=self.on_dup,
            verbose=self.verbose,
        )
        output_root = self._resolve_output_root(source, groups)
        total = sum(len(group.sources) for group in groups)

        edis: list[EDIFile] = []
        records: list[EDIRecomputeRecord] = []
        counter = 0
        iterator = _iter_progress_groups(
            groups,
            enabled=self.progress,
            total=total,
        )
        for group, item in iterator:
            counter += 1
            name = _safe_station(item) or _source_stem(item)
            try:
                ed = _read_group_item(
                    item,
                    strict=self.strict,
                    verbose=self.verbose,
                )
                if ed is None:
                    self._notify_progress(counter, total, name, "skipped", "")
                    continue
                out = self.recompute_one(ed)
                edis.append(out)
                dest = None
                if self.write:
                    dest = self._destination_for(
                        out,
                        source=ed,
                        line=group.name,
                        index=counter - 1,
                        output_root=output_root,
                    )
                    if dest.exists() and not self.overwrite:
                        raise FileExistsError(dest)
                    self._write_one(out, dest)
                rec = EDIRecomputeRecord(
                    source=_source_path(ed),
                    output=dest,
                    line=group.name,
                    station=station_name(out),
                    status="ok",
                )
                records.append(rec)
                self._notify_progress(counter, total, rec.station, "ok", "")
                if self.verbose:
                    msg = f"recomputed {rec.station}"
                    if rec.output is not None:
                        msg += f" -> {rec.output}"
                    print(msg)
            except Exception as exc:
                if self.strict:
                    raise
                message = f"{type(exc).__name__}: {exc}"
                rec = EDIRecomputeRecord(
                    source=_source_path(item),
                    output=None,
                    line=group.name,
                    station=name,
                    status="failed",
                    message=message,
                )
                records.append(rec)
                self._notify_progress(counter, total, name, "failed", message)
                if self.verbose:
                    print(f"failed {rec.station}: {rec.message}")

        sites = Sites(edis)
        result = EDIRecomputeResult(
            sites=sites,
            records=records,
            output_root=output_root if self.write else None,
            items=edis,
        )
        if self.write and self.manifest_csv:
            manifest = (
                output_root / "manifest.csv"
                if self.manifest_csv is True
                else Path(self.manifest_csv)
            )
            result.to_manifest(manifest)
        return result

    def _notify_progress(
        self,
        done: int,
        total: int,
        station: str,
        status: str,
        message: str = "",
    ) -> None:
        if self.progress_callback is not None:
            self.progress_callback(done, total, station, status, message)

    def recompute_one(self, edi: Any) -> EDIFile:
        """Recompute one EDI object and return the result."""

        return recompute_edi(
            edi,
            rotate_angle=self.rotate_angle,
            rotate_components=self.rotate_components,
            fmin=self.fmin,
            fmax=self.fmax,
            keep_freq=self.keep_freq,
            fill_missing_values=self.fill_missing_values,
            recompute_resphase=self.recompute_resphase,
            rename_policy=self.rename_policy,
            copy=self.copy,
        )

    def _resolve_output_root(
        self,
        source: Any,
        groups: Sequence[_LineGroup],
    ) -> Path:
        if self.output_root is not None:
            return Path(self.output_root).expanduser().resolve()

        if groups:
            root = groups[0].source_root
            if root is not None:
                return (root / self.output_name).resolve()

        if _is_pathlike(source):
            p = Path(source).expanduser().resolve()
            base = p.parent if p.is_file() else p
            return (base / self.output_name).resolve()

        return (Path.cwd() / self.output_name).resolve()

    def _destination_for(
        self,
        edi: EDIFile,
        *,
        source: Any,
        line: str | None,
        index: int,
        output_root: Path,
    ) -> Path:
        ctx = {
            "station": station_name(edi) or f"site_{index:04d}",
            "index": index,
            "line": line or "",
            "source_stem": _source_stem(source),
        }
        name = _render_name(self.template, ctx) or ctx["station"]
        if not name.lower().endswith(".edi"):
            name += ".edi"
        outdir = output_root
        if self.preserve_line_dirs and line:
            outdir = outdir / line
        return outdir / name

    def _write_one(self, edi: EDIFile, dest: Path) -> None:
        dest.parent.mkdir(parents=True, exist_ok=True)
        try:
            edi.write(
                new_edifn=dest.name,
                savepath=dest.parent,
                datatype=self.datatype,
                synthesize_spectra=self.synthesize_spectra,
                verbose=self.verbose,
            )
        except TypeError:
            edi.write(new_edifn=str(dest))


def recompute_edi(
    edi: Any,
    *,
    rotate_angle: float | None = None,
    rotate_components: Iterable[str] = ("Z", "Tip"),
    fmin: float | None = None,
    fmax: float | None = None,
    keep_freq: Iterable[int] | None = None,
    fill_missing_values: str | None = None,
    recompute_resphase: bool = True,
    rename_policy: Callable[[str], str] | None = None,
    copy: bool = True,
) -> EDIFile:
    r"""
    Recompute one EDI object.

    The operation is copy-returning by default. It can rotate
    transfer functions, subset frequencies, fill missing values,
    recompute apparent resistivity/phase, and rename station ids.
    """

    out = maybe_copy(edi) if copy else edi
    if rotate_angle is not None:
        out = _rotate_selected(
            out,
            float(rotate_angle),
            components=rotate_components,
        )
    if fmin is not None or fmax is not None or keep_freq is not None:
        out = select_freq(
            out,
            fmin=fmin,
            fmax=fmax,
            keep=keep_freq,
            inplace=True,
        )
    if fill_missing_values is not None:
        out = fill_missing(out, how=fill_missing_values, inplace=True)
    if recompute_resphase:
        out = recompute_res_phase(out, inplace=True)
    if rename_policy is not None:
        out = rename(out, policy=rename_policy, inplace=True)
    return out


def recompute_edis(source: Any, **kwargs: Any) -> EDIRecomputeResult:
    r"""
    Convenience function for :class:`EDIRecomputer`.

    Examples
    --------
    >>> result = recompute_edis("WILLY_DATA", rotate_angle=30.0)
    >>> result.paths
    [...]
    """

    return EDIRecomputer(**kwargs).run(source)


def _rotate_selected(
    edi: Any,
    angle_deg: float,
    *,
    components: Iterable[str],
) -> Any:
    comps = {str(c).strip().lower() for c in components}
    do_z = bool(comps & {"z", "imp", "impedance"})
    do_t = bool(comps & {"t", "tip", "tipper"})
    if not comps or (do_z and do_t):
        return _rotate_site(edi, angle_deg, inplace=True)

    r, rinv = _rotm(angle_deg)
    zobj = getattr(edi, "Z", None)

    if do_z and zobj is not None:
        zz = _get_attr_any(zobj, "z", "impedance", "_z")
        if zz is not None:
            arr = np.asarray(zz, complex)
            if arr.ndim == 3 and arr.shape[-2:] == (2, 2):
                _set_attr_first(
                    zobj,
                    ("z", "impedance", "_z"),
                    r[None, :, :] @ arr @ rinv[None, :, :],
                )
        ze = _get_attr_any(
            zobj, "z_error", "z_err", "impedance_err", "_z_err"
        )
        if ze is not None:
            err = np.asarray(ze, float)
            if err.ndim == 3 and err.shape[-2:] == (2, 2):
                ar = np.abs(r)
                ari = np.abs(rinv)
                _set_attr_first(
                    zobj,
                    ("z_error", "z_err", "impedance_err", "_z_err"),
                    ar[None, :, :] * err * ari[None, :, :],
                )

    if do_t:
        tip_obj = None
        for cand in ("T", "TIP", "Tip"):
            tip_obj = getattr(edi, cand, None)
            if tip_obj is not None:
                break
        if tip_obj is not None:
            tip = _get_attr_any(tip_obj, "tipper", "_tipper")
            if tip is not None:
                arr = np.asarray(tip, complex)
                if arr.ndim == 2 and arr.shape[1] == 2:
                    _set_attr_first(
                        tip_obj,
                        ("tipper", "_tipper"),
                        arr @ rinv.T,
                    )
                elif arr.ndim == 3 and arr.shape[-1] == 2:
                    _set_attr_first(
                        tip_obj,
                        ("tipper", "_tipper"),
                        arr @ rinv.T,
                    )
        elif zobj is not None:
            tip = _get_attr_any(zobj, "tipper", "tip", "_tipper")
            if tip is not None:
                arr = np.asarray(tip, complex)
                if arr.ndim == 2 and arr.shape[1] == 2:
                    _set_attr_first(
                        zobj,
                        ("tipper", "tip", "_tipper"),
                        arr @ rinv.T,
                    )
    return edi


def _discover_groups(
    source: Any,
    *,
    recursive: bool,
    strict: bool,
    on_dup: str,
    verbose: int,
) -> list[_LineGroup]:
    if _is_pathlike(source):
        return _discover_path_groups(Path(source), recursive=recursive)
    if _is_seq_of_pathlike(source):
        return _discover_many_path_groups(source, recursive=recursive)

    edis = _object_source_to_edis(source, strict=strict, verbose=verbose)
    return _discover_edi_object_groups(edis)


def _object_source_to_edis(
    source: Any,
    *,
    strict: bool,
    verbose: int,
) -> list[Any]:
    if hasattr(source, "get_section") and hasattr(source, "Z"):
        return [source]
    try:
        items = list(source)
    except Exception:
        items = [source]

    edis: list[Any] = []
    for item in items:
        raw = to_edis(
            item,
            copy=False,
            strict=strict,
            on_dup="keep",
            verbose=verbose,
        )
        if isinstance(raw, list):
            edis.extend(raw)
        elif raw is not None:
            edis.append(raw)
    return edis


def _discover_edi_object_groups(edis: Sequence[Any]) -> list[_LineGroup]:
    edis = list(edis)
    if not edis:
        return [_LineGroup(name=None, sources=[], source_root=None)]

    paths = [_source_path(ed) for ed in edis]
    valid_paths = [p for p in paths if p is not None]
    if not valid_paths:
        return [_LineGroup(name=None, sources=edis, source_root=None)]

    parents = [p.parent for p in valid_paths]
    unique_parents = sorted(set(parents))
    if len(unique_parents) == 1:
        parent = unique_parents[0]
        return [
            _LineGroup(
                name=parent.name,
                sources=edis,
                source_root=parent.parent,
            )
        ]

    try:
        root = Path(commonpath([str(parent) for parent in parents]))
    except Exception:
        root = unique_parents[0].parent

    grouped: dict[str | None, list[Any]] = {}
    for ed, path in zip(edis, paths):
        if path is None:
            grouped.setdefault(None, []).append(ed)
            continue
        try:
            rel = path.parent.relative_to(root)
            line = rel.parts[0] if rel.parts else path.parent.name
        except Exception:
            line = path.parent.name
        grouped.setdefault(line, []).append(ed)

    return [
        _LineGroup(name=line, sources=items, source_root=root)
        for line, items in grouped.items()
    ]


def _discover_path_groups(
    source: Path,
    *,
    recursive: bool,
) -> list[_LineGroup]:
    p = source.expanduser().resolve()
    if p.is_file():
        return [
            _LineGroup(
                name=p.parent.name,
                sources=[p],
                source_root=p.parent,
            )
        ]
    if not p.is_dir():
        return [_LineGroup(name=None, sources=[p], source_root=p.parent)]

    direct = sorted(p.glob("*.edi"))
    child_groups: list[_LineGroup] = []
    for child in sorted(p.iterdir()):
        if not child.is_dir():
            continue
        files = sorted(
            child.rglob("*.edi") if recursive else child.glob("*.edi")
        )
        if files:
            child_groups.append(
                _LineGroup(
                    name=child.name,
                    sources=files,
                    source_root=p,
                )
            )

    groups: list[_LineGroup] = []
    if direct:
        groups.append(
            _LineGroup(
                name=p.name,
                sources=direct,
                source_root=p.parent,
            )
        )
    if child_groups:
        groups.extend(child_groups)
    if groups:
        return groups

    files = sorted(p.rglob("*.edi") if recursive else p.glob("*.edi"))
    if files:
        return [
            _LineGroup(
                name=p.name,
                sources=files,
                source_root=p.parent,
            )
        ]
    return [_LineGroup(name=p.name, sources=[], source_root=p.parent)]


def _discover_many_path_groups(
    sources: Any,
    *,
    recursive: bool,
) -> list[_LineGroup]:
    groups: list[_LineGroup] = []
    for src in sources:
        groups.extend(_discover_path_groups(Path(src), recursive=recursive))
    return groups


def _read_group_item(
    item: Any,
    *,
    strict: bool,
    verbose: int,
) -> EDIFile | None:
    if isinstance(item, EDIFile):
        return item
    edi = getattr(item, "edi", None)
    if edi is not None:
        if hasattr(edi, "get_section") and hasattr(edi, "Z"):
            return edi
        if _is_pathlike(edi):
            return EDIFile(Path(edi), verbose=verbose)
    if _is_pathlike(item):
        return EDIFile(Path(item), verbose=verbose)
    if hasattr(item, "get_section") and hasattr(item, "Z"):
        return item
    if strict:
        raise TypeError(f"Not an EDI-like item: {type(item).__name__}")
    return None


def _iter_progress_groups(
    groups: Sequence[_LineGroup],
    *,
    enabled: bool | str,
    total: int,
) -> Iterator[tuple[_LineGroup, Any]]:
    items = [(group, item) for group in groups for item in group.sources]
    if not progress_enabled(enabled):
        yield from items
        return
    try:
        from tqdm import tqdm
    except Exception:
        yield from items
        return
    bar = tqdm(items, total=total, unit="edi", leave=False)
    for group, item in bar:
        name = _safe_station(item) or _source_stem(item)
        prefix = f"{group.name}/" if group.name else ""
        bar.set_description(f"Recomputing {prefix}{name}")
        yield group, item


def _write_manifest(
    records: Sequence[EDIRecomputeRecord],
    path: Path,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "source",
        "output",
        "line",
        "station",
        "status",
        "message",
    ]
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for rec in records:
            writer.writerow(
                {
                    "source": "" if rec.source is None else str(rec.source),
                    "output": "" if rec.output is None else str(rec.output),
                    "line": rec.line or "",
                    "station": rec.station,
                    "status": rec.status,
                    "message": rec.message,
                }
            )
    return path


class _SafeDict(dict):
    def __missing__(self, key: str) -> str:
        return ""


def _render_name(template: str, ctx: dict[str, Any]) -> str:
    return str(template).format_map(_SafeDict(ctx))


def _is_pathlike(obj: Any) -> bool:
    return isinstance(obj, (str, bytes, Path, PathLike))


def _is_seq_of_pathlike(x: Any) -> bool:
    if _is_pathlike(x):
        return False
    try:
        xs = list(x)
    except Exception:
        return False
    return bool(xs) and all(_is_pathlike(item) for item in xs)


def _source_path(obj: Any) -> Path | None:
    if _is_pathlike(obj):
        return Path(obj).expanduser().resolve()
    path = getattr(obj, "path", None)
    if path:
        try:
            return Path(path).expanduser().resolve()
        except Exception:
            return None
    return None


def _source_stem(obj: Any) -> str:
    path = _source_path(obj)
    if path is not None:
        return path.stem
    return _safe_station(obj)


def _safe_station(obj: Any) -> str:
    try:
        if _is_pathlike(obj):
            return Path(obj).stem
        return station_name(obj) or ""
    except Exception:
        return ""
