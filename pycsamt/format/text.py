# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""PCSM — ASCII sibling of PCSF (pyCSAMT Common Subsurface Markup).

PCSF (:mod:`pycsamt.format.io`, HDF5, ``.pcsf``) is the canonical,
machine-optimal container. PCSM (``.pcsm``) is a lossless, hand-editable
*text projection* of the exact same :class:`~pycsamt.format.schema.PCSFModel`
-- the same relationship netCDF (binary) has to CDL, its ``ncdump``/
``ncgen`` ASCII text form: PCSF stays the source of truth; PCSM is a
derived, round-trippable serialisation for scripting, teaching, diffing,
and hand inspection, not a second schema to keep in sync.

Grammar (deliberately not YAML/JSON, so a reader can be hand-written in
any language with no library, the same way Occam2D/ModEM/MARE2DEM's own
native ASCII files already are):

* ``KEYWORD value`` header lines for scalars.
* ``KEYWORD ... END_KEYWORD`` blocks for arrays -- values may be spread
  across any number of lines; the terminator, not a declared count,
  ends the block, so a file remains easy to hand-edit (reflow, add
  blank lines, add/remove values) and errors are caught by comparing
  the collected count against the count declared earlier (``NX``,
  ``N_STATIONS``, ...) rather than silently misaligning.
* ``#`` starts a comment: a whole line, or the remainder of a line
  after data, is ignored -- something none of the three native solver
  formats PCSM interoperates with support. The handful of genuinely
  free-text fields (``DESCRIPTION``, ``CREATED_BY``, ``CREATED_AT``,
  ``CRS``, ``SURVEY_JSON``, ``METADATA_JSON``) are the one exception:
  their value is everything after the keyword to end of line, taken
  verbatim, so a value that happens to contain ``#`` is not truncated.
* Floats are written with Python's ``repr()`` (the shortest string that
  round-trips exactly) -- unlike ModEM's own ASCII ``.rho`` format
  (~5 significant figures), a ``.pcsm`` round-trips a resistivity
  volume bit-exactly back through ``.pcsf``.

A large volume (a native ``grid3d``, a dense ``mesh_unstructured``
mesh) produces a large text file -- the same trade-off ModEM's own
ASCII format already accepts, inherent to any plain-text encoding of
bulk numeric data, not fixable while staying hand-editable. A path
ending in ``.gz`` is written/read gzip-compressed as a mitigation
(5.0x smaller on the real ModEM case in
``examples/pcsm_conversion_demo``) at the cost of no longer being
casually text-editor-readable, so it is opt-in rather than default.
"""

from __future__ import annotations

import gzip
import json
import re
import warnings
from datetime import datetime, timezone
from os import PathLike
from pathlib import Path
from typing import Any

import numpy as np

from ._version import check_pcsf_version
from .schema import (
    PCSF_VERSION,
    RESISTIVITY_UNIT,
    DerivedVolume,
    Grid2DGeometry,
    Grid3DGeometry,
    LineEntry,
    MultilineGeometry,
    PCSFModel,
    StationTable,
    TopographyPerStation,
    TopographyRaster,
    UnstructuredMeshGeometry,
)

__all__ = [
    "write_pcsm",
    "read_pcsm",
    "pcsf_to_pcsm",
    "pcsm_to_pcsf",
    "read_pcsf_or_pcsm",
    "peek_kind",
]

_LINE_ID_NA = "-"
_DEFAULT_ROW_WIDTH = 8


def _capped_row(n: int | None) -> int | None:
    """Cap a "one text line per geometric row" width at
    ``_DEFAULT_ROW_WIDTH``. A small grid still reads as one line per
    z-layer (nice, and still readable); a real, large grid (NX in the
    hundreds) does not turn into a single unreadable multi-thousand-
    character line just because that happens to match the grid's own
    row length -- readability wins over preserving row structure once
    a row is long enough that the structure isn't visible anyway.
    """
    return None if not n else min(n, _DEFAULT_ROW_WIDTH)


def _write_text(path: Path, text: str) -> None:
    if path.suffix == ".gz":
        with gzip.open(path, "wt", encoding="utf-8") as fh:
            fh.write(text)
    else:
        path.write_text(text, encoding="utf-8")


def _read_text(path: Path) -> str:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as fh:
            return fh.read()
    return path.read_text(encoding="utf-8")


# ---------------------------------------------------------------------
# Writer: tiny line-buffer helper
# ---------------------------------------------------------------------


class _Writer:
    def __init__(self) -> None:
        self.lines: list[str] = []

    def line(self, *tokens: Any, comment: str | None = None) -> None:
        text = " ".join(str(t) for t in tokens)
        if comment:
            text = f"{text}  # {comment}"
        self.lines.append(text)

    def raw(self, keyword: str, value: str) -> None:
        # Verbatim, no comment-stripping on read -- see module docstring.
        self.lines.append(f"{keyword} {value}")

    def block(
        self,
        keyword: str,
        values: Any,
        *,
        per_row: int | None = None,
        comment: str | None = None,
    ) -> None:
        header = keyword if not comment else f"{keyword}  # {comment}"
        self.lines.append(header)
        flat = np.asarray(values).reshape(-1)
        if per_row is None:
            # Wrap 1-D arrays at a fixed width rather than one giant
            # line -- readability matters most for exactly these
            # (coordinates, node ids), not just the explicitly
            # row-shaped 2-D/3-D blocks that pass per_row themselves.
            per_row = min(flat.size, _DEFAULT_ROW_WIDTH) or 1
        # Right-justify every value to the block's own widest formatted
        # value: same repr() precision (so round-trip fidelity is
        # unaffected), just padded so columns line up when the file is
        # opened in an editor -- the whole point of a *hand-editable*
        # sibling format.
        formatted = [_fmt(v) for v in flat]
        width = max((len(s) for s in formatted), default=0)
        for start in range(0, len(formatted), per_row):
            row = formatted[start : start + per_row]
            self.lines.append(" ".join(s.rjust(width) for s in row))
        self.lines.append(f"END_{keyword}")

    def aligned_rows(self, rows: list[list[str]], *, numeric: list[bool]) -> None:
        """Append *rows* (already-tokenised strings) with each column
        padded to that column's own widest value -- numeric columns
        right-justified, text columns (station names, ids) left-justified.
        """
        if not rows:
            return
        widths = [
            max(len(row[i]) for row in rows) for i in range(len(rows[0]))
        ]
        for row in rows:
            cells = [
                cell.rjust(widths[i]) if numeric[i] else cell.ljust(widths[i])
                for i, cell in enumerate(row)
            ]
            self.lines.append(" ".join(cells).rstrip())

    def text(self) -> str:
        return "\n".join(self.lines) + "\n"


def _fmt(value: Any) -> str:
    if isinstance(value, (np.floating, float)):
        return repr(float(value))
    return str(value)


# ---------------------------------------------------------------------
# Resistivity encoding labels -- always visible next to the data itself
# (not only in a separate header line), per the request that a reader
# must be able to tell linear from log10 at the point they see numbers.
# ---------------------------------------------------------------------

_LINEAR_LABEL = "linear ohm.m (canonical, see model/resistivity in SPEC.md S4)"
_LOG10_VIEW_LABEL = (
    "log10(ohm.m) -- derived view for human inspection only; "
    "read_pcsm() discards this, it is never part of the model"
)


def _native_label(encoding: str | None) -> str:
    return f"source-native encoding: {encoding}"


def _write_log10_view(
    w: _Writer, base_keyword: str, resistivity: np.ndarray, *, per_row: int | None
) -> None:
    with warnings.catch_warnings():
        # log10 of a non-positive resistivity (masked/placeholder cells)
        # legitimately produces -inf/nan; this view is write-only and
        # never parsed back, so that is fine -- suppress the runtime
        # warning rather than let it leak to the caller's console.
        warnings.simplefilter("ignore", category=RuntimeWarning)
        log10_values = np.log10(resistivity)
    w.block(
        f"{base_keyword}_LOG10", log10_values, per_row=per_row,
        comment=_LOG10_VIEW_LABEL,
    )


# ---------------------------------------------------------------------
# Reader: tokenised, comment-aware line cursor
# ---------------------------------------------------------------------

_VERBATIM_KEYWORDS = {
    "DESCRIPTION",
    "CREATED_BY",
    "CREATED_AT",
    "CRS",
    "SURVEY_JSON",
    "METADATA_JSON",
}


class _Cursor:
    """Walks a PCSM file's lines, stripping comments except on
    :data:`_VERBATIM_KEYWORDS` lines, and skipping blank lines.
    """

    def __init__(self, text: str) -> None:
        self._raw = text.splitlines()
        self._i = 0

    def _strip(self, raw: str) -> str:
        stripped = raw.lstrip()
        first_token = stripped.split(None, 1)[0] if stripped else ""
        if first_token in _VERBATIM_KEYWORDS:
            return raw.rstrip("\n")
        if "#" in raw:
            raw = raw[: raw.index("#")]
        return raw.rstrip()

    def peek(self) -> str | None:
        i = self._i
        while i < len(self._raw):
            line = self._strip(self._raw[i])
            if line.strip():
                return line
            i += 1
        return None

    def next(self) -> str | None:
        while self._i < len(self._raw):
            line = self._strip(self._raw[self._i])
            self._i += 1
            if line.strip():
                return line
        return None

    def expect_keyword(self, keyword: str) -> str:
        line = self.next()
        if line is None or not line.split(None, 1)[0] == keyword:
            raise ValueError(
                f"PCSM parse error: expected {keyword!r}, got {line!r}"
            )
        rest = line[len(keyword) :].strip()
        return rest

    def read_block(self, keyword: str) -> list[str]:
        header = self.next()
        if header is None or header.split(None, 1)[0] != keyword:
            raise ValueError(
                f"PCSM parse error: expected block {keyword!r}, got {header!r}"
            )
        end = f"END_{keyword}"
        tokens: list[str] = []
        while True:
            line = self.next()
            if line is None:
                raise ValueError(f"PCSM parse error: missing {end!r}")
            if line.strip() == end:
                return tokens
            tokens.extend(line.split())

    def try_peek_keyword(self) -> str | None:
        line = self.peek()
        if line is None:
            return None
        return line.split(None, 1)[0]


def _read_floats(cursor: _Cursor, keyword: str) -> np.ndarray:
    return np.asarray(cursor.read_block(keyword), dtype=float)


def _read_ints(cursor: _Cursor, keyword: str) -> np.ndarray:
    return np.asarray(cursor.read_block(keyword), dtype=np.int64)


def _read_scalar_line(cursor: _Cursor, keyword: str) -> str | None:
    if cursor.try_peek_keyword() != keyword:
        return None
    return cursor.expect_keyword(keyword)


# ---------------------------------------------------------------------
# Geometry (de)serialisation
# ---------------------------------------------------------------------


def _write_grid2d(w: _Writer, geo: Grid2DGeometry) -> None:
    w.line("NX", geo.x.shape[0])
    w.line("NZ", geo.z.shape[0])
    w.block("X_COORDS", geo.x)
    w.block("Z_COORDS", geo.z)
    if geo.x_nodes is not None:
        w.block("X_NODES", geo.x_nodes)
    if geo.z_nodes is not None:
        w.block("Z_NODES", geo.z_nodes)
    if geo.origin is not None:
        w.line("ORIGIN", *[_fmt(v) for v in geo.origin])
    if geo.azimuth_deg is not None:
        w.line("AZIMUTH_DEG", _fmt(geo.azimuth_deg))


def _read_grid2d(cursor: _Cursor) -> Grid2DGeometry:
    n_x = int(cursor.expect_keyword("NX"))
    n_z = int(cursor.expect_keyword("NZ"))
    x = _read_floats(cursor, "X_COORDS")
    z = _read_floats(cursor, "Z_COORDS")
    if x.size != n_x:
        raise ValueError(f"PCSM: X_COORDS has {x.size} values, NX declared {n_x}")
    if z.size != n_z:
        raise ValueError(f"PCSM: Z_COORDS has {z.size} values, NZ declared {n_z}")
    x_nodes = _read_floats(cursor, "X_NODES") if cursor.try_peek_keyword() == "X_NODES" else None
    z_nodes = _read_floats(cursor, "Z_NODES") if cursor.try_peek_keyword() == "Z_NODES" else None
    origin = None
    if cursor.try_peek_keyword() == "ORIGIN":
        origin = np.asarray(
            [float(v) for v in cursor.expect_keyword("ORIGIN").split()]
        )
    azimuth_deg = None
    raw_az = _read_scalar_line(cursor, "AZIMUTH_DEG")
    if raw_az is not None:
        azimuth_deg = float(raw_az)
    return Grid2DGeometry(
        x=x, z=z, x_nodes=x_nodes, z_nodes=z_nodes, origin=origin,
        azimuth_deg=azimuth_deg,
    )


def _write_grid3d(w: _Writer, geo: Grid3DGeometry) -> None:
    w.line("NX", geo.x.shape[0])
    w.line("NY", geo.y.shape[0])
    w.line("NZ", geo.z.shape[0])
    w.block("X_COORDS", geo.x)
    w.block("Y_COORDS", geo.y)
    w.block("Z_COORDS", geo.z)
    if geo.x_nodes is not None:
        w.block("X_NODES", geo.x_nodes)
    if geo.y_nodes is not None:
        w.block("Y_NODES", geo.y_nodes)
    if geo.z_nodes is not None:
        w.block("Z_NODES", geo.z_nodes)
    if geo.origin is not None:
        w.line("ORIGIN", *[_fmt(v) for v in geo.origin])
    w.line("ROTATION_DEG", _fmt(geo.rotation_deg))
    w.line("N_AIR", geo.n_air)


def _read_grid3d(cursor: _Cursor) -> Grid3DGeometry:
    n_x = int(cursor.expect_keyword("NX"))
    n_y = int(cursor.expect_keyword("NY"))
    n_z = int(cursor.expect_keyword("NZ"))
    x = _read_floats(cursor, "X_COORDS")
    y = _read_floats(cursor, "Y_COORDS")
    z = _read_floats(cursor, "Z_COORDS")
    for name, arr, n in (("X_COORDS", x, n_x), ("Y_COORDS", y, n_y), ("Z_COORDS", z, n_z)):
        if arr.size != n:
            raise ValueError(f"PCSM: {name} has {arr.size} values, declared {n}")
    x_nodes = _read_floats(cursor, "X_NODES") if cursor.try_peek_keyword() == "X_NODES" else None
    y_nodes = _read_floats(cursor, "Y_NODES") if cursor.try_peek_keyword() == "Y_NODES" else None
    z_nodes = _read_floats(cursor, "Z_NODES") if cursor.try_peek_keyword() == "Z_NODES" else None
    origin = None
    if cursor.try_peek_keyword() == "ORIGIN":
        origin = np.asarray(
            [float(v) for v in cursor.expect_keyword("ORIGIN").split()]
        )
    rotation_deg = float(cursor.expect_keyword("ROTATION_DEG"))
    n_air = int(cursor.expect_keyword("N_AIR"))
    return Grid3DGeometry(
        x=x, y=y, z=z, x_nodes=x_nodes, y_nodes=y_nodes, z_nodes=z_nodes,
        origin=origin, rotation_deg=rotation_deg, n_air=n_air,
    )


def _write_mesh(w: _Writer, geo: UnstructuredMeshGeometry) -> None:
    w.line("N_NODES", geo.nodes.shape[0])
    w.line("N_TRIANGLES", geo.connectivity.shape[0])
    w.line("PLANE", geo.plane)
    w.block("NODES", geo.nodes, per_row=geo.nodes.shape[1])
    w.block("CONNECTIVITY", geo.connectivity, per_row=3)
    w.block("REGION_IDS", geo.region_ids)


def _read_mesh(cursor: _Cursor) -> UnstructuredMeshGeometry:
    n_nodes = int(cursor.expect_keyword("N_NODES"))
    n_tri = int(cursor.expect_keyword("N_TRIANGLES"))
    plane = cursor.expect_keyword("PLANE")
    node_tokens = cursor.read_block("NODES")
    ndim = len(node_tokens) // n_nodes if n_nodes else 2
    nodes = np.asarray(node_tokens, dtype=float).reshape(n_nodes, ndim)
    conn_tokens = cursor.read_block("CONNECTIVITY")
    connectivity = np.asarray(conn_tokens, dtype=np.int64).reshape(n_tri, 3)
    region_ids = _read_ints(cursor, "REGION_IDS")
    return UnstructuredMeshGeometry(
        nodes=nodes, connectivity=connectivity, region_ids=region_ids, plane=plane,
    )


def _write_line_entry(w: _Writer, line: LineEntry, *, log10_view: bool = False) -> None:
    w.line("LINE_BEGIN", line.line_id)
    _write_grid2d(w, line.geometry)
    row = _capped_row(line.geometry.x.shape[0])
    w.block("RESISTIVITY", line.resistivity, per_row=row, comment=_LINEAR_LABEL)
    if log10_view:
        _write_log10_view(w, "RESISTIVITY", line.resistivity, per_row=row)
    w.line("OFFSET_Y", _fmt(line.offset_y))
    w.line("OFFSET_KIND", line.offset_kind)
    if line.azimuth_deg is not None:
        w.line("LINE_AZIMUTH_DEG", _fmt(line.azimuth_deg))
    w.line("LINE_END", line.line_id)


def _read_line_entry(cursor: _Cursor, line_id: str) -> LineEntry:
    got_id = cursor.expect_keyword("LINE_BEGIN")
    if got_id != line_id:
        raise ValueError(f"PCSM: expected LINE_BEGIN {line_id!r}, got {got_id!r}")
    geometry = _read_grid2d(cursor)
    resistivity = _read_floats(cursor, "RESISTIVITY").reshape(
        geometry.resistivity_shape
    )
    if cursor.try_peek_keyword() == "RESISTIVITY_LOG10":
        cursor.read_block("RESISTIVITY_LOG10")  # write-only view, discarded
    offset_y = float(cursor.expect_keyword("OFFSET_Y"))
    offset_kind = cursor.expect_keyword("OFFSET_KIND")
    azimuth_deg = None
    raw_az = _read_scalar_line(cursor, "LINE_AZIMUTH_DEG")
    if raw_az is not None:
        azimuth_deg = float(raw_az)
    end_id = cursor.expect_keyword("LINE_END")
    if end_id != line_id:
        raise ValueError(f"PCSM: expected LINE_END {line_id!r}, got {end_id!r}")
    return LineEntry(
        line_id=line_id, geometry=geometry, resistivity=resistivity,
        offset_y=offset_y, offset_kind=offset_kind, azimuth_deg=azimuth_deg,
    )


def _write_multiline(w: _Writer, geo: MultilineGeometry, *, log10_view: bool = False) -> None:
    w.line("LINE_ORDER", *[line.line_id for line in geo.lines])
    for line in geo.lines:
        _write_line_entry(w, line, log10_view=log10_view)
    if geo.derived_volume is not None:
        dv = geo.derived_volume
        w.line("DERIVED_VOLUME_BEGIN")
        _write_grid3d(w, dv.grid)
        w.block("DV_RESISTIVITY", dv.resistivity, per_row=_capped_row(dv.grid.x.shape[0]))
        w.line("DERIVATION_METHOD", dv.derivation_method)
        w.line("SYNTHESIZED", "true" if dv.synthesized else "false")
        w.line("DERIVED_FROM", *dv.derived_from)
        w.line("DERIVED_VOLUME_END")


def _read_multiline(cursor: _Cursor) -> MultilineGeometry:
    line_ids = cursor.expect_keyword("LINE_ORDER").split()
    lines = [_read_line_entry(cursor, line_id) for line_id in line_ids]
    derived_volume = None
    if cursor.try_peek_keyword() == "DERIVED_VOLUME_BEGIN":
        cursor.expect_keyword("DERIVED_VOLUME_BEGIN")
        grid = _read_grid3d(cursor)
        resistivity = _read_floats(cursor, "DV_RESISTIVITY").reshape(
            grid.resistivity_shape
        )
        derivation_method = cursor.expect_keyword("DERIVATION_METHOD")
        synthesized = cursor.expect_keyword("SYNTHESIZED").lower() == "true"
        derived_from = cursor.expect_keyword("DERIVED_FROM").split()
        cursor.expect_keyword("DERIVED_VOLUME_END")
        derived_volume = DerivedVolume(
            grid=grid, resistivity=resistivity, derivation_method=derivation_method,
            derived_from=derived_from, synthesized=synthesized,
        )
    return MultilineGeometry(lines=lines, derived_volume=derived_volume)


_GEOMETRY_WRITERS = {
    "grid2d": _write_grid2d,
    "grid3d": _write_grid3d,
    "mesh_unstructured": _write_mesh,
    "multiline": _write_multiline,
}
_GEOMETRY_READERS = {
    "grid2d": _read_grid2d,
    "grid3d": _read_grid3d,
    "mesh_unstructured": _read_mesh,
    "multiline": _read_multiline,
}


# ---------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------


def write_pcsm(
    model: PCSFModel, path: str | PathLike, *, log10_view: bool = False
) -> Path:
    """Write a :class:`PCSFModel` to a ``.pcsm`` (ASCII) file.

    The same model :func:`~pycsamt.format.io.write_pcsf` would encode
    as HDF5, encoded instead as human-readable, hand-editable text.
    See the module docstring for the grammar.

    Every resistivity block carries its own encoding as an inline
    comment right next to the data, not only in a separate header line
    -- e.g. ``RESISTIVITY  # linear ohm.m (canonical, ...)`` -- so a
    reader never has to guess or scroll elsewhere to tell linear from
    log10. The canonical ``RESISTIVITY`` field is **always** linear
    ohm.m (PCSF's own binding rule, SPEC.md section 2): this cannot be
    changed by an option, since a reader/consumer relies on that
    invariant unconditionally, and reversing a ``10**log10(x)``
    round-trip is not guaranteed bit-exact the way the rest of this
    format's ``repr()``-based float round-trip is.

    Pass ``log10_view=True`` to additionally write a
    ``RESISTIVITY_LOG10`` block (and, for a ``multiline`` file, one per
    line) alongside it -- a convenience for reading resistivity in the
    log space many EM inversions (Occam2D among them) actually work
    in. This is a write-only, clearly-labelled derived view:
    :func:`read_pcsm` discards it, it never becomes part of the
    returned :class:`PCSFModel`, and it never influences the canonical
    ``RESISTIVITY`` block. (An Occam2D-sourced model already carries
    its own real log10 array as ``RESISTIVITY_NATIVE`` whenever
    ``resistivity_native`` is set -- that block is unaffected by this
    option and is written either way.)

    A path ending in ``.gz`` (e.g. ``model.pcsm.gz``) is written
    gzip-compressed. Large volumes (a native ``grid3d`` or a dense
    ``mesh_unstructured`` mesh) produce large plain-text files -- the
    same trade-off ModEM's own ASCII ``.rho`` format already accepts --
    and gzip substantially reduces that in practice (measured on the
    real, bundled ``willy_27freq_watex_line02_sample`` 590,400-cell
    ModEM volume, see ``examples/pcsm_conversion_demo``: 19.7 MB
    plain, 3.9 MB gzipped -- 5.0x smaller, real inverted resistivity
    compressing far better than synthetic/random test data would --
    vs. 4.1 MB for the equivalent ``.pcsf``). A gzipped file is no
    longer something to casually open in a text editor, so this stays
    opt-in rather than the default: PCSM's purpose is hand-editability,
    which a compressed file gives up.

    Parameters
    ----------
    model : PCSFModel
        The model to serialize. Validated before anything is written.
    path : path-like
        Destination file. Parent directories are created if missing.
        A ``.gz`` suffix writes a gzip-compressed file.
    log10_view : bool, default False
        Also write a ``RESISTIVITY_LOG10`` convenience block (see
        above). Off by default so existing files/output size are
        unaffected unless explicitly requested.

    Returns
    -------
    pathlib.Path
        The path written to.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.format import Grid2DGeometry, PCSFModel
    >>> from pycsamt.format.text import write_pcsm, read_pcsm
    >>> geometry = Grid2DGeometry(x=np.array([0.0, 100.0]), z=np.array([10.0, 50.0]))
    >>> model = PCSFModel(
    ...     geometry=geometry,
    ...     resistivity=np.array([[100.0, 120.0], [50.0, 60.0]]),
    ...     source_backend="occam2d",
    ... )
    >>> path = write_pcsm(model, "example.pcsm")  # doctest: +SKIP
    >>> round_tripped = read_pcsm(path)  # doctest: +SKIP
    >>> gz_path = write_pcsm(model, "example.pcsm.gz")  # doctest: +SKIP
    >>> with_log10 = write_pcsm(model, "example_log10.pcsm", log10_view=True)  # doctest: +SKIP
    """
    if not isinstance(model, PCSFModel):
        raise TypeError(f"model must be a PCSFModel, got {type(model)!r}")
    model.validate()

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    w = _Writer()
    w.lines.append(
        "# PCSM -- pyCSAMT Common Subsurface Markup (ASCII sibling of a .pcsf file)"
    )
    w.line("PCSM_VERSION", PCSF_VERSION)
    w.line("SOURCE_BACKEND", model.source_backend)
    w.raw("CREATED_BY", model.created_by or "")
    w.raw("CREATED_AT", model.created_at or datetime.now(timezone.utc).isoformat())
    w.line("RESISTIVITY_UNIT", RESISTIVITY_UNIT)
    if model.crs:
        w.raw("CRS", model.crs)
    if model.description:
        w.raw("DESCRIPTION", model.description)
    w.lines.append("")

    w.line("GEOMETRY_KIND", model.geometry.kind)
    if model.geometry.kind == "multiline":
        _write_multiline(w, model.geometry, log10_view=log10_view)
    else:
        _GEOMETRY_WRITERS[model.geometry.kind](w, model.geometry)
    w.lines.append("")

    if model.geometry.kind != "multiline":
        if model.resistivity_native_encoding is not None:
            w.line("RESISTIVITY_NATIVE_ENCODING", model.resistivity_native_encoding)
        row = _capped_row(
            model.geometry.x.shape[0] if hasattr(model.geometry, "x") else None
        )
        w.block("RESISTIVITY", model.resistivity, per_row=row, comment=_LINEAR_LABEL)
        if model.resistivity_native is not None:
            w.block(
                "RESISTIVITY_NATIVE", model.resistivity_native, per_row=row,
                comment=_native_label(model.resistivity_native_encoding),
            )
        if log10_view:
            _write_log10_view(w, "RESISTIVITY", model.resistivity, per_row=row)
        if model.uncertainty is not None:
            w.block("UNCERTAINTY", model.uncertainty, per_row=row)
        if model.sensitivity is not None:
            w.block("SENSITIVITY", model.sensitivity, per_row=row)
        if model.resistivity_by_region is not None:
            w.block("RESISTIVITY_BY_REGION", model.resistivity_by_region)
        if model.resistivity_by_node is not None:
            w.block("RESISTIVITY_BY_NODE", model.resistivity_by_node)
        w.lines.append("")

    if model.stations is not None:
        n = len(model.stations.name)
        w.line(
            "STATIONS",
            comment="name x y z line_id lon lat (lon/lat: nan if unset)",
        )
        line_ids = model.stations.line_id or [_LINE_ID_NA] * n
        lons = (
            model.stations.lon
            if model.stations.lon is not None
            else np.full(n, np.nan)
        )
        lats = (
            model.stations.lat
            if model.stations.lat is not None
            else np.full(n, np.nan)
        )
        rows = [
            [str(name), _fmt(x), _fmt(y), _fmt(z), str(lid), _fmt(lon), _fmt(lat)]
            for name, x, y, z, lid, lon, lat in zip(
                model.stations.name, model.stations.x, model.stations.y,
                model.stations.z, line_ids, lons, lats,
            )
        ]
        w.aligned_rows(
            rows, numeric=[False, True, True, True, False, True, True]
        )
        w.lines.append("END_STATIONS")
        w.lines.append("")

    if model.topography is not None:
        w.line("TOPOGRAPHY_KIND", model.topography.kind)
        if model.topography.kind == "raster":
            w.block("TOPOGRAPHY_X", model.topography.x)
            w.block("TOPOGRAPHY_Y", model.topography.y)
            w.block(
                "TOPOGRAPHY_ELEVATION",
                model.topography.elevation,
                per_row=_capped_row(model.topography.x.shape[0]),
            )
        else:
            w.line("TOPOGRAPHY", comment="station_id elevation")
            rows = [
                [str(sid), _fmt(elev)]
                for sid, elev in zip(
                    model.topography.station_id, model.topography.elevation
                )
            ]
            w.aligned_rows(rows, numeric=[False, True])
            w.lines.append("END_TOPOGRAPHY")
        w.lines.append("")

    if model.survey:
        w.raw("SURVEY_JSON", json.dumps(model.survey, default=_json_default))
    if model.metadata:
        w.raw("METADATA_JSON", json.dumps(model.metadata, default=_json_default))
    if model.survey or model.metadata:
        w.lines.append("")

    if model.history:
        w.line("HISTORY_BEGIN")
        for key, values in model.history.items():
            w.lines.append(f"{key} " + " ".join(_fmt(v) for v in values))
        w.line("HISTORY_END")

    _write_text(path, w.text())
    return path


def read_pcsm(path: str | PathLike) -> PCSFModel:
    """Read a :class:`PCSFModel` back from a ``.pcsm`` (ASCII) file.

    Parameters
    ----------
    path : path-like
        Source file.

    Returns
    -------
    PCSFModel
        Fully reconstructed and re-validated model.

    Raises
    ------
    ValueError
        If ``PCSM_VERSION`` is missing, malformed, or names an
        unrecognised MAJOR version (see ``pycsamt/format/SPEC.md``
        section 5), if the file is otherwise malformed, or if the
        reconstructed model fails :meth:`PCSFModel.validate`.
    """
    path = Path(path)
    cursor = _Cursor(_read_text(path))

    version = _read_scalar_line(cursor, "PCSM_VERSION")
    if version is None:
        raise ValueError(f"{path}: missing required 'PCSM_VERSION' header")
    check_pcsf_version(version, PCSF_VERSION)

    source_backend = cursor.expect_keyword("SOURCE_BACKEND")
    created_by = cursor.expect_keyword("CREATED_BY")
    created_at = cursor.expect_keyword("CREATED_AT")
    cursor.expect_keyword("RESISTIVITY_UNIT")
    crs = _read_scalar_line(cursor, "CRS")
    description = _read_scalar_line(cursor, "DESCRIPTION") or ""

    kind = cursor.expect_keyword("GEOMETRY_KIND")
    reader = _GEOMETRY_READERS.get(kind)
    if reader is None:
        raise ValueError(
            f"{path}: unknown geometry kind {kind!r}; expected one of "
            f"{tuple(_GEOMETRY_READERS)}"
        )
    geometry = reader(cursor)

    resistivity = None
    resistivity_native = None
    resistivity_native_encoding = None
    uncertainty = None
    sensitivity = None
    resistivity_by_region = None
    resistivity_by_node = None
    if kind != "multiline":
        if cursor.try_peek_keyword() == "RESISTIVITY_NATIVE_ENCODING":
            resistivity_native_encoding = cursor.expect_keyword(
                "RESISTIVITY_NATIVE_ENCODING"
            )
        resistivity = _read_shaped(cursor, "RESISTIVITY", geometry, kind)
        if cursor.try_peek_keyword() == "RESISTIVITY_NATIVE":
            resistivity_native = _read_shaped(
                cursor, "RESISTIVITY_NATIVE", geometry, kind
            )
        if cursor.try_peek_keyword() == "RESISTIVITY_LOG10":
            cursor.read_block("RESISTIVITY_LOG10")  # write-only view, discarded
        if cursor.try_peek_keyword() == "UNCERTAINTY":
            uncertainty = _read_shaped(cursor, "UNCERTAINTY", geometry, kind)
        if cursor.try_peek_keyword() == "SENSITIVITY":
            sensitivity = _read_shaped(cursor, "SENSITIVITY", geometry, kind)
        if cursor.try_peek_keyword() == "RESISTIVITY_BY_REGION":
            resistivity_by_region = _read_floats(cursor, "RESISTIVITY_BY_REGION")
        if cursor.try_peek_keyword() == "RESISTIVITY_BY_NODE":
            resistivity_by_node = _read_floats(cursor, "RESISTIVITY_BY_NODE")

    stations = None
    if cursor.try_peek_keyword() == "STATIONS":
        cursor.next()  # header line (may carry a trailing comment, discarded)
        names: list[str] = []
        xs: list[float] = []
        ys: list[float] = []
        zs: list[float] = []
        lids: list[str] = []
        lons: list[float] = []
        lats: list[float] = []
        while True:
            line = cursor.next()
            if line is None:
                raise ValueError("PCSM: missing END_STATIONS")
            if line.strip() == "END_STATIONS":
                break
            parts = line.split()
            if len(parts) == 7:
                name, x, y, z, lid, lon, lat = parts
            elif len(parts) == 5:
                # A file written before lon/lat columns existed --
                # still readable, just carries no real-world position.
                name, x, y, z, lid = parts
                lon, lat = "nan", "nan"
            else:
                raise ValueError(
                    "PCSM: malformed STATIONS row (expected 5 or 7 "
                    f"fields, got {len(parts)}): {line!r}"
                )
            names.append(name)
            xs.append(float(x))
            ys.append(float(y))
            zs.append(float(z))
            lids.append(lid)
            lons.append(float(lon))
            lats.append(float(lat))
        line_id = None if all(v == _LINE_ID_NA for v in lids) else lids
        lon_arr = np.asarray(lons)
        lat_arr = np.asarray(lats)
        has_lonlat = not np.all(np.isnan(lon_arr))
        stations = StationTable(
            name=names, x=np.asarray(xs), y=np.asarray(ys), z=np.asarray(zs),
            line_id=line_id,
            lon=lon_arr if has_lonlat else None,
            lat=lat_arr if has_lonlat else None,
        )

    topography = None
    if cursor.try_peek_keyword() == "TOPOGRAPHY_KIND":
        topo_kind = cursor.expect_keyword("TOPOGRAPHY_KIND")
        if topo_kind == "raster":
            tx = _read_floats(cursor, "TOPOGRAPHY_X")
            ty = _read_floats(cursor, "TOPOGRAPHY_Y")
            telev = _read_floats(cursor, "TOPOGRAPHY_ELEVATION").reshape(
                ty.size, tx.size
            )
            topography = TopographyRaster(x=tx, y=ty, elevation=telev)
        else:
            cursor.next()  # TOPOGRAPHY header line
            sids: list[str] = []
            elevs: list[float] = []
            while True:
                line = cursor.next()
                if line is None:
                    raise ValueError("PCSM: missing END_TOPOGRAPHY")
                if line.strip() == "END_TOPOGRAPHY":
                    break
                sid, elev = line.split()
                sids.append(sid)
                elevs.append(float(elev))
            topography = TopographyPerStation(
                station_id=sids, elevation=np.asarray(elevs)
            )

    survey: dict[str, Any] = {}
    if cursor.try_peek_keyword() == "SURVEY_JSON":
        survey = json.loads(cursor.expect_keyword("SURVEY_JSON"))
    metadata: dict[str, Any] = {}
    if cursor.try_peek_keyword() == "METADATA_JSON":
        metadata = json.loads(cursor.expect_keyword("METADATA_JSON"))

    history: dict[str, np.ndarray] = {}
    if cursor.try_peek_keyword() == "HISTORY_BEGIN":
        cursor.next()
        while True:
            line = cursor.next()
            if line is None:
                raise ValueError("PCSM: missing HISTORY_END")
            if line.strip() == "HISTORY_END":
                break
            key, *values = line.split()
            history[key] = np.asarray(values, dtype=float)

    model = PCSFModel(
        geometry=geometry,
        resistivity=resistivity,
        resistivity_native=resistivity_native,
        resistivity_native_encoding=resistivity_native_encoding,
        uncertainty=uncertainty,
        sensitivity=sensitivity,
        resistivity_by_region=resistivity_by_region,
        resistivity_by_node=resistivity_by_node,
        stations=stations,
        topography=topography,
        survey=survey,
        history=history,
        source_backend=source_backend,
        created_by=created_by,
        created_at=created_at,
        crs=crs,
        description=description,
        metadata=metadata,
    )
    model.validate()
    return model


def _read_shaped(cursor: _Cursor, keyword: str, geometry: Any, kind: str) -> np.ndarray:
    flat = _read_floats(cursor, keyword)
    if kind == "mesh_unstructured":
        return flat
    return flat.reshape(geometry.resistivity_shape)


def _json_default(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    raise TypeError(f"object of type {type(value)!r} is not JSON serializable")


def pcsf_to_pcsm(
    pcsf_path: str | PathLike, pcsm_path: str | PathLike, *, log10_view: bool = False
) -> Path:
    """Convert a ``.pcsf`` (HDF5) file to a ``.pcsm`` (ASCII) file.

    ``log10_view`` is passed through to :func:`write_pcsm` — see there.
    """
    from .io import read_pcsf

    return write_pcsm(read_pcsf(pcsf_path), pcsm_path, log10_view=log10_view)


def pcsm_to_pcsf(pcsm_path: str | PathLike, pcsf_path: str | PathLike) -> Path:
    """Convert a ``.pcsm`` (ASCII) file to a ``.pcsf`` (HDF5) file."""
    from .io import write_pcsf

    return write_pcsf(read_pcsm(pcsm_path), pcsf_path)


def read_pcsf_or_pcsm(path: str | PathLike) -> PCSFModel:
    """Read a :class:`PCSFModel` from either encoding, by extension.

    A single entry point for a consumer (e.g. ``app/mapview``'s
    inversion-result importer) that wants to accept whichever of the
    two lossless PCSF encodings (see ``SPEC.md`` S9) a user hands it,
    without duplicating the ``.pcsf`` vs. ``.pcsm``/``.pcsm.gz``
    dispatch itself. Extension is the only signal used -- content is
    never sniffed -- matching every other reader in this package.

    Parameters
    ----------
    path : path-like
        A ``.pcsf`` (binary HDF5), ``.pcsm`` (ASCII), or ``.pcsm.gz``
        (gzip-compressed ASCII) file.

    Raises
    ------
    ValueError
        If *path*'s extension is none of the three recognised forms.
    """
    name = str(path).lower()
    if name.endswith(".pcsm") or name.endswith(".pcsm.gz"):
        return read_pcsm(path)
    if name.endswith(".pcsf"):
        from .io import read_pcsf

        return read_pcsf(path)
    raise ValueError(
        f"{path}: unrecognised PCSF/PCSM extension -- expected one of "
        "'.pcsf', '.pcsm', '.pcsm.gz'."
    )


def peek_kind(path: str | PathLike) -> str:
    """Return a PCSF/PCSM file's ``geometry.kind`` without loading any
    array data -- just the root/geometry attributes for ``.pcsf``, or
    the ``GEOMETRY_KIND`` header line for ``.pcsm``/``.pcsm.gz``.

    Meant for a caller that wants to label several candidate files
    (e.g. a file picker UI, deciding which are ``grid2d``/``multiline``
    and therefore importable by a given consumer) before committing to
    a full :func:`read_pcsf_or_pcsm` load of any one of them.

    Raises
    ------
    ValueError
        If *path*'s extension is not recognised, or the kind
        attribute/header cannot be found (e.g. a truncated file).
    """
    name = str(path).lower()
    if name.endswith(".pcsm") or name.endswith(".pcsm.gz"):
        text = _read_text(Path(path))
        match = re.search(r"^GEOMETRY_KIND\s+(\S+)", text, re.MULTILINE)
        if not match:
            raise ValueError(f"{path}: no GEOMETRY_KIND header found")
        return match.group(1)
    if name.endswith(".pcsf"):
        import h5py

        with h5py.File(path, "r") as fh:
            kind = fh["geometry"].attrs.get("kind")
        if kind is None:
            raise ValueError(f"{path}: geometry group has no 'kind' attribute")
        return str(kind)
    raise ValueError(
        f"{path}: unrecognised PCSF/PCSM extension -- expected one of "
        "'.pcsf', '.pcsm', '.pcsm.gz'."
    )
