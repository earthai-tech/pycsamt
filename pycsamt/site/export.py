# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

import csv
import zipfile
from collections.abc import Iterable
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

from .utils import (
    _ensure_head,
    station_name,
)

__all__ = ["write_site", "write_sites", "pack_zip"]


def write_site(site: Any, path: str | Path) -> Path:
    r"""
    Write a single site (EDI) to a target path.

    This is a thin, best-effort adapter around several common EDI
    writer spellings. The function will create parent directories
    as needed and then try, in order, the following methods on
    ``site`` until one succeeds:

    1) ``write(new_edifn=path)``
    2) ``write(path)``
    3) ``to_file(path)``
    4) ``save(path)``

    If none of these exist or succeed, a ``RuntimeError`` is raised.

    Parameters
    ----------
    site : Any
        An EDI-like object. It can be a ``pycsamt.seg.edi.EDIFile``
        or any object exposing one of the writer methods listed
        above.
    path : str or pathlib.Path
        Destination file path. Parent directories are created
        if they do not exist.

    Returns
    -------
    pathlib.Path
        The resolved output path.

    Notes
    -----
    This function does not enforce overwrite policy. Whether an
    existing file is replaced depends on the underlying writer
    implementation of the provided ``site`` object.

    Examples
    --------
    >>> from pathlib import Path
    >>> from pycsamt.site.export import write_site
    >>> class Dummy:
    ...     def to_file(self, p):  # minimal writer
    ...         Path(p).write_text("# dummy edi\\n", encoding="utf-8")
    >>> out = write_site(Dummy(), Path("out") / "S01.edi")
    >>> out.name
    'S01.edi'
    >>> out.exists()
    True

    See Also
    --------
    pycsamt.site.export.write_sites : Batch writing with templated
        filenames and optional manifest.
    pycsamt.site.export.pack_zip : Archive a set of sites into a zip.

    References
    ----------
    .. [1] Python Software Foundation. "pathlib" and "io" modules.
    """

    out = Path(path)
    _ensure_parent(out)
    _write_via_backend(site, out)
    return out


def write_sites(
    sites: Any,
    outdir: str | Path,
    *,
    template: str = "{station}.edi",
    exist_ok: bool = False,
    manifest_csv: str | Path | None = None,
) -> list[Path]:
    r"""
    Write a collection of sites to a directory using a filename
    template.

    The function accepts many input forms (``Sites``, an
    ``EDICollection``, any iterable of EDI-like objects, or a single
    object) and writes each item to ``outdir``. Filenames are
    rendered from a context and the ``template`` string.

    Supported template keys (filled via a safe formatter):

    - ``{station}`` : current station name
    - ``{index}`` : zero-based index in the iteration order
    - ``{lat}``, ``{lon}``, ``{elev}`` : header coordinates, or NaN
    - ``{chainage}`` : optional header chainage, or NaN

    If the rendered name does not end with ``.edi``, the extension
    is appended automatically.

    Parameters
    ----------
    sites : Any
        A ``Sites`` instance, an ``EDICollection``, any iterable of
        EDI-like objects, or a single EDI-like object. EDI-like means
        it implements one of: ``write(new_edifn=...)``, ``write(...)``,
        ``to_file(...)``, or ``save(...)``.
    outdir : str or pathlib.Path
        Output directory. It is created if it does not exist.
    template : str, optional
        Filename template. Defaults to ``"{station}.edi"``.
    exist_ok : bool, optional
        If ``False`` (default), raise ``FileExistsError`` on the first
        name collision inside ``outdir``. If ``True``, allow
        overwriting subject to the writer behavior.
    manifest_csv : str or pathlib.Path or None, optional
        If provided, write a CSV manifest with one row per written
        site. Columns are:
        ``index, station, lat, lon, elev, chainage, filename, path``.

    Returns
    -------
    list of pathlib.Path
        Paths to the files written, in the same order as the input
        iteration.

    Notes
    -----
    The ``index`` used in templating and in the manifest is the
    zero-based position in the input order. Coordinate fields come
    from the EDI header when available; missing values are written
    as NaN.

    Examples
    --------
    >>> from pathlib import Path
    >>> from pycsamt.site.export import write_sites
    >>> class EdiToFile:
    ...     def __init__(self, name):
    ...         self._n = name
    ...
    ...     def to_file(self, p):
    ...         Path(p).write_text(f"# {self._n}\\n", encoding="utf-8")
    ...
    ...     # station name is taken from header helpers when present,
    ...     # but the template can still use {index}.
    >>> outdir = Path("eds_out")
    >>> paths = write_sites(
    ...     [EdiToFile("S01"), EdiToFile("S02")],
    ...     outdir,
    ...     template="{index:03d}_{station}",
    ... )
    >>> [p.exists() for p in paths]
    [True, True]

    >>> # Write with a manifest
    >>> mpath = Path("eds_out") / "manifest.csv"
    >>> _ = write_sites(
    ...     [EdiToFile("S01"), EdiToFile("S02")],
    ...     outdir,
    ...     template="{station}",
    ...     exist_ok=True,
    ...     manifest_csv=mpath,
    ... )
    >>> mpath.exists()
    True

    See Also
    --------
    pycsamt.site.base.Sites.write : Higher-level convenience bound to
        a ``Sites`` collection.
    pycsamt.site.export.pack_zip : Create a zip archive instead of a
        directory tree.

    References
    ----------
    .. [1] Python Software Foundation. "csv" module.
    """

    out_root = Path(outdir)
    out_root.mkdir(parents=True, exist_ok=True)

    written: list[tuple[int, Any, Path]] = []
    used: set[str] = set()
    for i, ed in enumerate(_iter_any_sites(sites)):
        ctx = _context_for(ed, i)
        name = _render_name(template, ctx) or ""
        # Derive a non-empty stem. An empty station name would render to
        # ".edi" (a stem-less, often hidden file) — ``ensure_sites`` then
        # ignores it, so the export silently loads back as zero stations
        # ("No stations with valid impedance data found"). Fall back to a
        # positional name so every site lands in a real file.
        stem = name[:-4] if name.lower().endswith(".edi") else name
        stem = stem.strip()
        if not stem:
            stem = f"site_{i:03d}"
        name = f"{stem}.edi"
        # Disambiguate collisions. Two sites sharing a station name would
        # otherwise overwrite each other on disk (silent data loss while
        # still returning one path per input), corrupting the dataset.
        if name in used:
            name = f"{stem}_{i:03d}.edi"
        used.add(name)
        dest = out_root / name
        if dest.exists() and not exist_ok:
            raise FileExistsError(dest)
        _write_via_backend(ed, dest)
        written.append((i, ed, dest))

    if manifest_csv:
        rows = _rows_for_manifest(written)
        _write_manifest_csv(rows, Path(manifest_csv))

    return [p for _, _, p in written]


def pack_zip(
    sites: Any,
    out_zip: str | Path,
    *,
    template: str = "{station}.edi",
    manifest_csv: str | Path | None = None,
) -> Path:
    r"""
    Pack a set of sites into a zip archive using a filename template.

    Each input item is written to a temporary directory first, then
    added to the ``out_zip`` archive using ``ZIP_DEFLATED``. Filenames
    inside the archive are rendered from the same context as in
    :func:`write_sites`. If a name lacks the ``.edi`` suffix, it is
    appended automatically.

    Optionally, a CSV manifest can be written alongside the archive.

    Parameters
    ----------
    sites : Any
        A ``Sites`` instance, an ``EDICollection``, any iterable of
        EDI-like objects, or a single EDI-like object.
    out_zip : str or pathlib.Path
        Destination zip file path. Parent directories are created as
        needed.
    template : str, optional
        Filename template for entries stored in the archive.
        Defaults to ``"{station}.edi"``.
    manifest_csv : str or pathlib.Path or None, optional
        If provided, write a CSV manifest next to the zip. Columns:
        ``index, station, lat, lon, elev, chainage, filename, path``.

    Returns
    -------
    pathlib.Path
        The path to the created zip archive.

    Notes
    -----
    Files are staged in a temporary directory and then compressed
    with ``zipfile.ZIP_DEFLATED``. The ``index`` used in templating
    and the manifest corresponds to the input iteration order. This
    function does not delete or modify any original EDI sources.

    Examples
    --------
    >>> from pathlib import Path
    >>> from zipfile import ZipFile
    >>> from pycsamt.site.export import pack_zip
    >>> class EdiSave:
    ...     def __init__(self, name):
    ...         self._n = name
    ...
    ...     def to_file(self, p):
    ...         Path(p).write_text(f"# {self._n}\\n", encoding="utf-8")
    >>> zpath = Path("out_bundle") / "sites.zip"
    >>> out = pack_zip(
    ...     [EdiSave("A01"), EdiSave("A02")],
    ...     zpath,
    ...     template="{station}.edi",
    ...     manifest_csv=Path("out_bundle") / "manifest.csv",
    ... )
    >>> out == zpath, zpath.exists()
    (True, True)
    >>> with ZipFile(zpath, "r") as zf:
    ...     sorted(zf.namelist())
    ['A01.edi', 'A02.edi']

    See Also
    --------
    pycsamt.site.export.write_sites : Write files to a directory
        instead of an archive.

    References
    ----------
    .. [1] Python Software Foundation. "zipfile" module.
    """

    out_zip = Path(out_zip)
    _ensure_parent(out_zip)

    rows: list[dict[str, Any]] = []
    with TemporaryDirectory() as td:
        tmp = Path(td)
        file_map: list[tuple[Path, str]] = []
        rows: list[dict] = []

        for i, ed in enumerate(_iter_any_sites(sites)):
            ctx = _context_for(ed, i)
            name = _render_name(template, ctx) or "site.edi"
            if not name.lower().endswith(".edi"):
                name = f"{name}.edi"
            dst = tmp / name
            _ensure_parent(dst)
            _write_via_backend(ed, dst)
            file_map.append((dst, name))
            rows.append(
                {
                    "index": ctx["index"],
                    "station": ctx["station"],
                    "lat": ctx["lat"],
                    "lon": ctx["lon"],
                    "elev": ctx["elev"],
                    "chainage": ctx["chainage"],
                    "filename": name,
                    "path": str(out_zip),
                }
            )

        with zipfile.ZipFile(out_zip, "w", compression=zipfile.ZIP_DEFLATED) as zf:
            for fp, arcname in file_map:
                zf.write(fp, arcname)

    if manifest_csv:
        _write_manifest_csv(rows, Path(manifest_csv))

    return out_zip


# ----- helpers ------------


class _SafeDict(dict):
    def __missing__(self, key):
        return ""


def _ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def _get_number(v: Any) -> float:
    try:
        if v is None:
            return float("nan")
        return float(v)
    except Exception:
        return float("nan")


def _ctx_coords(h: Any) -> tuple[float, float, float]:
    # be robust to different backends / spellings
    la = getattr(h, "lat", None)
    if la is None:
        la = getattr(h, "latitude", None)

    lo = getattr(h, "lon", None)
    if lo is None:
        lo = getattr(h, "long", None)
    if lo is None:
        lo = getattr(h, "longitude", None)

    ev = getattr(h, "elev", None)
    if ev is None:
        ev = getattr(h, "elevation", None)

    return _get_number(la), _get_number(lo), _get_number(ev)


def _context_for(ed: Any, index: int) -> dict[str, Any]:
    h = _ensure_head(ed)
    nm = station_name(ed)
    lat, lon, elev = _ctx_coords(h)
    ch = getattr(ed, "chainage", None)
    return {
        "station": str(nm or ""),
        "index": int(index),
        "lat": float(lat),
        "lon": float(lon),
        "elev": float(elev),
        "chainage": _get_number(ch),
    }


def _render_name(template: str, ctx: dict[str, Any]) -> str:
    return str(template).format_map(_SafeDict(ctx))


def _write_via_backend(ed: Any, path: Path) -> None:
    """
    Best-effort write using common backend spellings.
    """
    # pycsamt EDI backends usually accept `write(new_edifn=...)`
    try:
        ed.write(new_edifn=path.name, savepath=path.parent)
        if path.exists():
            return
    except Exception:
        pass

    # Some backends accept a complete destination as `new_edifn`.
    try:
        ed.write(new_edifn=str(path))
        if path.exists():
            return
    except Exception:
        pass

    # or `write(path)`
    try:
        ed.write(str(path))
        if path.exists():
            return
    except Exception:
        pass

    # or `to_file(path)`
    try:
        ed.to_file(str(path))
        if path.exists():
            return
    except Exception:
        pass

    # or `save(path)`
    try:
        ed.save(str(path))
        if path.exists():
            return
    except Exception:
        pass

    raise RuntimeError("Cannot write EDI file for site")


def _rows_for_manifest(
    items: Iterable[tuple[int, Any, Path]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for idx, ed, fp in items:
        c = _context_for(ed, idx)
        rows.append(
            {
                "index": c["index"],
                "station": c["station"],
                "lat": c["lat"],
                "lon": c["lon"],
                "elev": c["elev"],
                "chainage": c["chainage"],
                "filename": str(fp.name),
                "path": str(fp),
            }
        )
    return rows


def _write_manifest_csv(rows: list[dict[str, Any]], csv_path: Path) -> None:
    if not rows:
        return
    _ensure_parent(csv_path)
    fields = [
        "index",
        "station",
        "lat",
        "lon",
        "elev",
        "chainage",
        "filename",
        "path",
    ]
    with csv_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def _iter_any_sites(sites: Any):
    # Accept Sites wrapper
    try:
        from .base import Sites

        if isinstance(sites, Sites):
            # Sites exposes .as_list() of EDI-like items
            for ed in sites.as_list():
                yield ed
            return
    except Exception:
        pass

    # Accept EDICollection (if present)
    from .utils import as_edicollection

    coll = as_edicollection(sites)
    if coll is not None:
        for ed in coll:
            yield ed
        return

    # Accept generic iterables
    if hasattr(sites, "__iter__") and not isinstance(sites, (str, bytes)):
        for ed in sites:
            yield ed
        return

    # Fallback: single item
    yield sites
