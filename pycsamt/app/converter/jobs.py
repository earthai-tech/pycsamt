# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.app.converter.jobs
============================

Plain-Python, Qt-free conversion jobs behind every page of the
converter app. Each function takes plain paths/options and returns a
JSON-friendly ``dict`` (or raises a plain exception on failure) -- no
widget, dialog, or Qt import appears anywhere in this module, so every
function here is directly unit-testable with ``pytest`` and no display.

Widgets only collect input and hand it to one of these functions inside
:class:`pycsamt.app.converter.workers.ConversionWorker`.

Every function reuses the same engine the CLI drives
(:mod:`pycsamt.format.convert_engine`, :mod:`pycsamt.emtf`,
:mod:`pycsamt.format.borehole`, :mod:`pycsamt.format.geology`,
:mod:`pycsamt.format.structure`, :mod:`pycsamt.format.pointset`) --
this module adds no parallel conversion logic of its own.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

__all__ = [
    "detect_source_job",
    "classify_batch_item",
    "run_batch_job",
    "convert_to_pcsf_job",
    "transcode_pcsf_pcsm_job",
    "validate_pcsf_job",
    "info_pcsf_job",
    "edi_to_xml_job",
    "xml_to_edi_job",
    "build_pcbh_job",
    "build_pcgl_job",
    "build_pcgs_job",
    "build_pcpt_job",
]


# ---------------------------------------------------------------------------
# Detection
# ---------------------------------------------------------------------------


def detect_source_job(path: Path, solver: str | None = None) -> dict[str, Any]:
    """Classify *path* -- returns :meth:`SourceKind.to_dict`."""
    from pycsamt.format import convert_engine

    sk = convert_engine.detect(Path(path), solver)
    return sk.to_dict()


# ---------------------------------------------------------------------------
# Batch queue
# ---------------------------------------------------------------------------


def classify_batch_item(path: Path) -> str:
    """Cheap kind guess for the batch queue: ``"edi"``, ``"xml"``,
    ``"inversion"`` (any PCSF-convertible source), or ``"unknown"``.
    """
    from pycsamt.format import convert_engine

    path = Path(path)
    if path.is_file() and path.suffix.lower() == ".edi":
        return "edi"
    if path.is_file() and path.suffix.lower() == ".xml":
        return "xml"
    try:
        sk = convert_engine.detect(path)
    except (FileNotFoundError, ValueError):
        return "unknown"
    if sk.category in {"solver", "ai_arrays", "pcsf", "pcsm"} and sk.convertible:
        return "inversion"
    return "unknown"


def run_batch_job(
    items: list[dict[str, Any]],
    output_dir: Path,
    *,
    to_format: str = "pcsf",
    overwrite: bool = False,
    prefer_spectra: bool = True,
    on_loss: str = "warn",
    xml_strict: bool = True,
    progress: Any = None,
) -> list[dict[str, Any]]:
    """Run a mixed queue of items sequentially.

    ``items`` is a list of ``{"path": ..., "kind": "edi"|"xml"|"inversion"}``
    (see :func:`classify_batch_item`). Each item is converted independently
    -- one failure is recorded and does not stop the rest of the queue.
    """
    output_dir = Path(output_dir)
    results: list[dict[str, Any]] = []
    total = len(items)

    for idx, item in enumerate(items, start=1):
        path = Path(item["path"])
        kind = item.get("kind") or classify_batch_item(path)
        entry: dict[str, Any] = {"path": str(path), "kind": kind}
        try:
            if kind == "inversion":
                report = convert_to_pcsf_job(
                    path, None, to_format, output_dir, overwrite=overwrite
                )
                entry.update(status="done", detail=report["file"])
            elif kind == "edi":
                written = edi_to_xml_job(
                    path, output_dir, prefer_spectra=prefer_spectra, overwrite=overwrite
                )
                entry.update(status="done", detail=written[0]["target"])
            elif kind == "xml":
                written = xml_to_edi_job(
                    path, output_dir, on_loss=on_loss, strict=xml_strict, overwrite=overwrite
                )
                entry.update(status="done", detail=written[0]["target"])
            else:
                entry.update(status="skipped", detail="Unrecognized source kind")
        except Exception as exc:  # noqa: BLE001 - one bad row must not abort the queue
            entry.update(status="error", detail=str(exc))
        results.append(entry)
        if progress is not None:
            progress(idx, total, path.name)

    return results


# ---------------------------------------------------------------------------
# Inversion / AI arrays / PCSF / PCSM -> PCSF or PCSM
# ---------------------------------------------------------------------------


def convert_to_pcsf_job(
    source: Path,
    target: Path | None,
    to_format: str | None,
    output_dir: Path,
    *,
    overwrite: bool = False,
    solver: str | None = None,
    iteration: int | None = None,
    topo: Path | None = None,
    epsg: int | None = None,
    utm_zone: str | None = None,
    encoding: str | None = None,
    station_z_convention: str = "auto",
    air_threshold_ohm_m: float | None = 1e8,
    poly: Path | None = None,
    origin: tuple[float, ...] | None = None,
    azimuth_deg: float | None = None,
    created_by: str = "pyCSAMT Format Studio",
    description: str = "",
    log10_view: bool = False,
) -> dict[str, Any]:
    """Any solver result / AI array bundle / .pcsf / .pcsm -> .pcsf or .pcsm.

    Mirrors ``pycsamt format convert``: auto-detects *source*, then
    either transcodes (PCSF<->PCSM) or rebuilds a
    :class:`~pycsamt.format.schema.PCSFModel` via
    :mod:`pycsamt.format.convert_engine` and writes it out.
    """
    from pycsamt.format import convert_engine

    source = Path(source)
    sk = convert_engine.detect(source, solver)
    if not sk.convertible:
        raise ValueError(f"Nothing to convert — {source} looks like: {sk.detail}")

    dst, canonical = convert_engine.resolve_target(
        source, Path(target) if target else None, to_format, Path(output_dir)
    )

    same_format = (
        sk.category in {"pcsf", "pcsm"}
        and canonical.split(".")[0] == sk.category
    )
    if same_format and (target is None or source.resolve() == dst.resolve()):
        raise ValueError(
            f"{source.name} is already {sk.category.upper()}; choose a "
            "different target or output format."
        )

    if dst.exists() and not overwrite:
        raise FileExistsError(f"{dst} exists — enable overwrite to replace it.")

    if sk.category in {"pcsf", "pcsm"}:
        report = _transcode_pcsf_pcsm(source, dst, sk.category, canonical, log10_view)
        report["transcoded"] = True
        report["detected"] = sk.to_dict()
        return report

    opts = convert_engine.ConvertOptions(
        iteration=iteration,
        topo=Path(topo) if topo else None,
        epsg=epsg,
        utm_zone=utm_zone,
        encoding=encoding,
        station_z_convention=station_z_convention,
        air_threshold_ohm_m=air_threshold_ohm_m,
        poly=Path(poly) if poly else None,
        origin=origin,
        azimuth_deg=azimuth_deg,
        created_by=created_by,
        description=description,
    )
    model = convert_engine.build_model(sk, opts)
    written = convert_engine.write_model(model, dst, canonical, log10_view)
    report = convert_engine.model_report(model, Path(written))
    report["transcoded"] = False
    report["detected"] = sk.to_dict()
    return report


def _transcode_pcsf_pcsm(
    source: Path,
    dst: Path,
    src_category: str,
    canonical: str,
    log10_view: bool,
) -> dict[str, Any]:
    from pycsamt.format import convert_engine
    from pycsamt.format.text import (
        pcsf_to_pcsm,
        pcsm_to_pcsf,
        read_pcsf_or_pcsm,
    )

    dst.parent.mkdir(parents=True, exist_ok=True)
    target_base = canonical.split(".")[0]

    if src_category == "pcsf" and target_base == "pcsm":
        pcsf_to_pcsm(source, dst, log10_view=log10_view)
    elif src_category == "pcsm" and target_base == "pcsf":
        pcsm_to_pcsf(source, dst)
    else:
        model = read_pcsf_or_pcsm(source)
        convert_engine.write_model(model, dst, target_base, log10_view)

    model = read_pcsf_or_pcsm(dst)
    return convert_engine.model_report(model, dst)


def transcode_pcsf_pcsm_job(
    source: Path,
    dst: Path,
    target_format: str,
    *,
    log10_view: bool = False,
    overwrite: bool = False,
) -> dict[str, Any]:
    """Direct PCSF<->PCSM transcode without touching solver/AI sources."""
    from pycsamt.format.text import peek_kind  # noqa: F401 (validates readability)

    source = Path(source)
    dst = Path(dst)
    if dst.exists() and not overwrite:
        raise FileExistsError(f"{dst} exists — enable overwrite to replace it.")
    src_category = "pcsm" if source.name.lower().endswith((".pcsm", ".pcsm.gz")) else "pcsf"
    return _transcode_pcsf_pcsm(source, dst, src_category, target_format, log10_view)


# ---------------------------------------------------------------------------
# Validate / Info
# ---------------------------------------------------------------------------


def validate_pcsf_job(file: Path, *, roundtrip: bool = True) -> dict[str, Any]:
    """Structural + round-trip check -- same steps as ``pycsamt format validate``."""
    import tempfile

    import numpy as np

    from pycsamt.format.io import write_pcsf
    from pycsamt.format.text import peek_kind, read_pcsf_or_pcsm, write_pcsm

    file = Path(file)
    name = file.name.lower()
    is_pcsm = name.endswith((".pcsm", ".pcsm.gz"))
    if not (is_pcsm or name.endswith(".pcsf")):
        raise ValueError(f"{file.name} is not a .pcsf / .pcsm / .pcsm.gz file.")

    checks: list[dict[str, Any]] = []

    def _record(step: str, ok: bool, detail: str = "") -> bool:
        checks.append({"step": step, "ok": ok, "detail": detail})
        return ok

    kind = None
    try:
        kind = peek_kind(file)
        _record("header", True, f"geometry={kind}")
    except Exception as exc:  # noqa: BLE001
        _record("header", False, str(exc))

    model = None
    if kind is not None:
        try:
            model = read_pcsf_or_pcsm(file)
            _record("load", True, f"kind={model.kind}")
        except Exception as exc:  # noqa: BLE001
            _record("load", False, str(exc))

    if model is not None:
        try:
            model.validate()
            _record("schema", True)
        except Exception as exc:  # noqa: BLE001
            _record("schema", False, str(exc))

    if model is not None and roundtrip:
        try:
            with tempfile.TemporaryDirectory() as tmp:
                if is_pcsm:
                    rt = Path(tmp) / "rt.pcsm"
                    write_pcsm(model, rt)
                else:
                    rt = Path(tmp) / "rt.pcsf"
                    write_pcsf(model, rt)
                back = read_pcsf_or_pcsm(rt)
                a, b = model.resistivity, back.resistivity
                if a is None and b is None:
                    _record("roundtrip", True, "no resistivity array")
                elif a is None or b is None:
                    _record("roundtrip", False, "resistivity presence changed")
                else:
                    np.testing.assert_allclose(
                        np.asarray(a), np.asarray(b),
                        rtol=1e-9, atol=0, equal_nan=True,
                    )
                    _record("roundtrip", True, f"resistivity {list(a.shape)} matches")
        except Exception as exc:  # noqa: BLE001
            _record("roundtrip", False, str(exc))

    ok = all(c["ok"] for c in checks) and len(checks) > 0
    return {"file": str(file), "valid": ok, "checks": checks}


def info_pcsf_job(file: Path) -> dict[str, Any]:
    """Full summary of a .pcsf / .pcsm file -- same fields as ``format info``."""
    import numpy as np

    from pycsamt.format.text import peek_kind, read_pcsf_or_pcsm

    file = Path(file)
    name = file.name.lower()
    if not name.endswith((".pcsf", ".pcsm", ".pcsm.gz")):
        raise ValueError(f"{file.name} is not a .pcsf / .pcsm / .pcsm.gz file.")

    kind = peek_kind(file)
    model = read_pcsf_or_pcsm(file)
    geom = model.geometry

    report: dict[str, Any] = {
        "file": str(file),
        "size_bytes": file.stat().st_size,
        "container": "pcsm" if name.endswith((".pcsm", ".pcsm.gz")) else "pcsf",
        "geometry_kind": kind,
        "source_backend": model.source_backend,
        "created_by": model.created_by,
        "created_at": model.created_at,
        "description": model.description,
        "crs": model.crs,
        "native_encoding": model.resistivity_native_encoding,
    }

    for label, arr in (
        ("resistivity", model.resistivity),
        ("uncertainty", model.uncertainty),
        ("sensitivity", model.sensitivity),
    ):
        if arr is None:
            continue
        arr = np.asarray(arr, dtype=float)
        finite = arr[np.isfinite(arr)]
        entry: dict[str, Any] = {"shape": list(arr.shape), "n_nan": int(arr.size - finite.size)}
        if finite.size:
            entry.update(
                min=float(finite.min()), max=float(finite.max()),
                median=float(np.median(finite)),
            )
        report[label] = entry

    if model.stations is not None:
        st = model.stations
        report["stations"] = {
            "n": len(st.name),
            "has_lonlat": st.lon is not None and st.lat is not None,
        }
    if model.topography is not None:
        report["has_topography"] = True

    del geom
    return report


# ---------------------------------------------------------------------------
# EDI <-> EMTF-XML
# ---------------------------------------------------------------------------


def _iter_matching(source: Path, patterns: tuple[str, ...]) -> list[Path]:
    if source.is_file():
        return [source]
    found: dict[Path, None] = {}
    for pattern in patterns:
        for path in source.glob(pattern):
            found.setdefault(path.resolve(), None)
    return sorted(found)


def edi_to_xml_job(
    source: Path,
    output_dir: Path,
    *,
    prefer_spectra: bool = True,
    overwrite: bool = False,
    progress: Any = None,
) -> list[dict[str, Any]]:
    """Convert one .edi file or a directory of them to EMTF-XML."""
    from pycsamt.emtf import edi_to_emtf, write_emtf_xml

    output_dir = Path(output_dir)
    sources = _iter_matching(Path(source), ("*.edi", "*.EDI"))
    if not sources:
        raise FileNotFoundError(f"No .edi files found under {source}.")

    results = []
    for idx, path in enumerate(sources, start=1):
        document = edi_to_emtf(path, prefer_spectra=prefer_spectra)
        dst = output_dir / f"{path.stem}.xml"
        if dst.exists() and not overwrite:
            raise FileExistsError(f"{dst} exists — enable overwrite to replace it.")
        write_emtf_xml(document, dst)
        results.append({"source": str(path), "target": str(dst)})
        if progress is not None:
            progress(idx, len(sources), path.name)
    return results


def xml_to_edi_job(
    source: Path,
    output_dir: Path,
    *,
    on_loss: str = "warn",
    strict: bool = True,
    overwrite: bool = False,
    progress: Any = None,
) -> list[dict[str, Any]]:
    """Convert one EMTF-XML file or a directory of them to EDI."""
    from pycsamt.emtf import EMTFXMLReader, write_edi

    output_dir = Path(output_dir)
    sources = _iter_matching(Path(source), ("*.xml", "*.XML"))
    if not sources:
        raise FileNotFoundError(f"No .xml files found under {source}.")

    reader = EMTFXMLReader(strict=strict)
    results = []
    for idx, path in enumerate(sources, start=1):
        document = reader.read(path)
        name = document.station or path.stem
        dst = output_dir / f"{name}.edi"
        if dst.exists() and not overwrite:
            raise FileExistsError(f"{dst} exists — enable overwrite to replace it.")
        write_edi(document, dst, on_loss=on_loss)
        results.append({"source": str(path), "target": str(dst)})
        if progress is not None:
            progress(idx, len(sources), path.name)
    return results


# ---------------------------------------------------------------------------
# PCBH builder
# ---------------------------------------------------------------------------


def build_pcbh_job(
    source: Path,
    output: Path,
    *,
    source_kind: str = "auto",
    collar_id: str | None = None,
    x: float | None = None,
    y: float | None = None,
    z: float | None = None,
    crs_horizontal: str | None = None,
    document_id: str | None = None,
    created_by: str = "pyCSAMT Format Studio",
    overwrite: bool = False,
) -> dict[str, Any]:
    """Build a PCBH document from a combined CSV, a relational CSV
    directory, an XLSX workbook, or a single LAS log (with an explicit
    collar). Same dispatch as ``pycsamt format build-pcbh``."""
    from pycsamt.format.borehole.jsonio import write_pcbh

    source = Path(source)
    output = Path(output)
    if output.exists() and not overwrite:
        raise FileExistsError(f"{output} exists — enable overwrite to replace it.")

    kind = source_kind if source_kind != "auto" else _infer_pcbh_kind(source)
    report = None

    if kind == "csv":
        from pycsamt.format.borehole.csvio import boreholes_from_csv

        document, report = boreholes_from_csv(
            source, document_id=document_id, created_by=created_by
        )
    elif kind == "csv-dir":
        from pycsamt.format.borehole.relational import boreholes_from_csv_directory

        document, report = boreholes_from_csv_directory(
            source, document_id=document_id, created_by=created_by
        )
    elif kind == "xlsx":
        from pycsamt.format.borehole.xlsxio import boreholes_from_xlsx

        document, report = boreholes_from_xlsx(
            source, document_id=document_id, created_by=created_by
        )
    elif kind == "las":
        if None in (collar_id, x, y, z, crs_horizontal):
            raise ValueError(
                "LAS import needs collar_id, x, y, z and crs_horizontal."
            )
        from pycsamt.format.borehole.lasio import borehole_from_las
        from pycsamt.format.borehole.schema import Collar

        collar = Collar(x=x, y=y, z=z)
        document, report = borehole_from_las(
            source, collar=collar, crs_horizontal=crs_horizontal
        )
        document.document_id = document_id or document.document_id
        document.created_by = created_by
    else:
        raise ValueError(f"Unknown PCBH source kind: {kind!r}")

    written = write_pcbh(document, output)
    result: dict[str, Any] = {
        "file": str(written),
        "document_id": document.document_id,
        "n_boreholes": len(document.boreholes),
    }
    if report is not None:
        result["rows_read"] = report.rows_read
        result["rows_accepted"] = report.rows_accepted
        result["rows_rejected"] = report.rows_rejected
        result["issues"] = [
            {
                "severity": i.severity,
                "row": i.row,
                "column": i.column,
                "message": i.message,
            }
            for i in report.issues
        ]
    return result


def _infer_pcbh_kind(source: Path) -> str:
    if source.is_dir():
        return "csv-dir"
    suffix = source.suffix.lower()
    if suffix in (".xlsx", ".xlsm"):
        return "xlsx"
    if suffix == ".las":
        return "las"
    return "csv"


# ---------------------------------------------------------------------------
# PCGL builder
# ---------------------------------------------------------------------------


def build_pcgl_job(
    source: Path,
    output: Path,
    *,
    title: str = "",
    document_id: str | None = None,
    created_by: str = "pyCSAMT Format Studio",
    overwrite: bool = False,
) -> dict[str, Any]:
    """Build a PCGL geology legend from a CSV -- see
    :class:`pycsamt.format.geology.GeologyLegend`."""
    from pycsamt.format.geology import legend_from_csv, write_legend

    output = Path(output)
    if output.exists() and not overwrite:
        raise FileExistsError(f"{output} exists — enable overwrite to replace it.")

    legend = legend_from_csv(
        Path(source), document_id=document_id, created_by=created_by, title=title
    )
    written = write_legend(legend, output)
    return {
        "file": str(written),
        "document_id": legend.document_id,
        "title": legend.title,
        "n_entries": len(legend.entries),
    }


# ---------------------------------------------------------------------------
# PCGS builder
# ---------------------------------------------------------------------------


def build_pcgs_job(
    output: Path,
    *,
    planar_path: Path | None = None,
    linear_path: Path | None = None,
    faults_path: Path | None = None,
    title: str = "",
    document_id: str | None = None,
    created_by: str = "pyCSAMT Format Studio",
    overwrite: bool = False,
) -> dict[str, Any]:
    """Build a PCGS structural-evidence document from up to three CSVs --
    see :class:`pycsamt.format.structure.StructModel`."""
    from pycsamt.format.structure import structure_from_csv, write_structure

    if not (planar_path or linear_path or faults_path):
        raise ValueError("Pass at least one of planar_path / linear_path / faults_path.")

    output = Path(output)
    if output.exists() and not overwrite:
        raise FileExistsError(f"{output} exists — enable overwrite to replace it.")

    model = structure_from_csv(
        planar_path=Path(planar_path) if planar_path else None,
        linear_path=Path(linear_path) if linear_path else None,
        faults_path=Path(faults_path) if faults_path else None,
        document_id=document_id,
        created_by=created_by,
        title=title,
    )
    written = write_structure(model, output)
    return {
        "file": str(written),
        "document_id": model.document_id,
        "n_planar": len(model.model.planar),
        "n_linear": len(model.model.linear),
        "n_faults": len(model.model.faults),
    }


# ---------------------------------------------------------------------------
# PCPT builder
# ---------------------------------------------------------------------------


def build_pcpt_job(
    source: Path,
    output: Path,
    *,
    sheet: str | int | None = None,
    header_row: int | None = None,
    crs: str | None = None,
    document_id: str | None = None,
    overwrite: bool = False,
) -> dict[str, Any]:
    """Build a PCPT points-of-interest document from a CSV or XLSX -- see
    :class:`pycsamt.format.pointset.PointSet`."""
    from pycsamt.format.pointset import points_from_csv, points_from_xlsx, write_points

    source = Path(source)
    output = Path(output)
    if output.exists() and not overwrite:
        raise FileExistsError(f"{output} exists — enable overwrite to replace it.")

    if source.suffix.lower() in (".xlsx", ".xlsm"):
        point_set = points_from_xlsx(
            source, sheet=sheet, header_row=header_row, crs=crs,
            document_id=document_id,
        )
    else:
        point_set = points_from_csv(source, crs=crs, document_id=document_id)

    written = write_points(point_set, output)
    return {
        "file": str(written),
        "document_id": point_set.document_id,
        "n_points": len(point_set.points),
        "crs": point_set.crs,
    }
