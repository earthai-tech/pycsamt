# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Loss-explicit LAS 2.0 subset import and export for PCBH."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from ...api.property import PyCSAMTObject
from .mapping import ImportReport
from .schema import (
    Collar,
    CoordinateReferenceSystem,
    LogInterval,
    PCBHBorehole,
    PCBHDocument,
    Trajectory,
    VocabularyEntry,
)

__all__ = ["LASExportReport", "borehole_from_las", "write_las_subset"]


@dataclass(repr=False)
class LASExportReport(PyCSAMTObject):
    """Preservation and loss summary for a LAS subset export."""

    path: str
    curves_written: list[str] = field(default_factory=list)
    losses: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


def borehole_from_las(
    path: str | Path,
    *,
    collar: Collar,
    crs_horizontal: str,
    kind: str = "unknown",
    status: str = "unknown",
    depth_curve: str = "DEPT",
    resistivity_curve: str = "RESD",
    lithology_curve: str | None = "LITH",
    max_samples: int = 250_000,
) -> tuple[PCBHDocument, ImportReport]:
    """Import a LAS 2.0 curve subset while retaining curve metadata."""
    if not isinstance(collar, Collar):
        raise TypeError("collar must be a PCBH Collar")
    collar.validate()
    if not isinstance(crs_horizontal, str) or not crs_horizontal.strip():
        raise ValueError("crs_horizontal must be a non-empty string")
    if isinstance(max_samples, bool) or not isinstance(max_samples, int):
        raise TypeError("max_samples must be an integer")
    if max_samples <= 0:
        raise ValueError("max_samples must be greater than zero")
    source = Path(path)
    sections, curves, rows = _read_las(source, max_samples=max_samples)
    if _well_text(sections, "WRAP").upper() not in {"", "NO"}:
        raise ValueError("wrapped LAS data is not supported; require WRAP.NO")
    depth_key = depth_curve.upper()
    if depth_key not in curves:
        raise ValueError(f"depth curve {depth_curve!r} is absent")
    depth_unit = curves[depth_key]["unit"].upper()
    depth_factor = {"M": 1.0, "METRE": 1.0, "METER": 1.0, "FT": 0.3048}.get(
        depth_unit
    )
    if depth_factor is None:
        raise ValueError(f"unsupported LAS depth unit {depth_unit!r}")
    index = [value * depth_factor for value in curves[depth_key]["index"]]
    null = _well_float(sections, "NULL", -9999.25)
    samples = []
    for mnemonic, curve in curves.items():
        if mnemonic == depth_key:
            continue
        values = []
        for depth, row in zip(index, rows):
            raw = row[curve["column"]]
            value = None if math.isclose(raw, null, abs_tol=1e-12) else raw
            values.append([depth, value])
        samples.append(
            {
                "mnemonic": mnemonic,
                "name": curve["description"],
                "unit": curve["unit"],
                "index": "md",
                "samples": values,
                "null": None,
            }
        )
    raw_step = _well_float(
        sections, "STEP", _median_step(index) / depth_factor
    )
    step = abs(raw_step * depth_factor)
    if any(second <= first for first, second in zip(index, index[1:])):
        raise ValueError("LAS depth curve must increase strictly")
    total_depth = max(index) + step if index else 0.0
    if total_depth <= 0:
        raise ValueError(
            "LAS depth curve does not define positive total depth"
        )
    well_name = _well_text(sections, "WELL") or source.stem
    intervals, vocabulary = _las_intervals(
        index,
        rows,
        curves,
        lithology_curve,
        resistivity_curve,
        null,
        step,
    )
    hole = PCBHBorehole(
        id=well_name,
        name=well_name,
        kind=kind,
        status=status,
        collar=collar,
        total_depth_md=total_depth,
        trajectory=Trajectory(method="vertical"),
        interval_logs={"lithology": intervals} if intervals else {},
        extensions={
            "pcbh:continuous_curves": samples,
            "pcbh:las_metadata": {
                "version": _well_text(sections, "VERS"),
                "wrap": _well_text(sections, "WRAP"),
                "null_value": null,
                "well": dict(sections.get("W", {})),
            },
        },
    )
    document = PCBHDocument(
        document_id=f"las:{source.stem}",
        created_at=datetime.now(timezone.utc)
        .isoformat()
        .replace("+00:00", "Z"),
        created_by="pycsamt LAS importer",
        crs=CoordinateReferenceSystem(crs_horizontal),
        boreholes=[hole],
        lithologies=vocabulary,
    )
    document.validate()
    report = ImportReport(
        source=str(source),
        source_sha256=__import__("hashlib")
        .sha256(source.read_bytes())
        .hexdigest(),
        delimiter="whitespace",
        strict=True,
        rows_read=len(rows),
        rows_accepted=len(rows),
    )
    report.source_files[source.name] = report.source_sha256
    report.inferred_values.append(
        f"total_depth_md={total_depth!r} from last depth plus STEP"
    )
    if depth_factor != 1.0:
        report.unit_conversions.append(
            f"LAS depth {depth_unit} converted to PCBH metres"
        )
    return document, report


def write_las_subset(
    borehole: PCBHBorehole,
    path: str | Path,
    *,
    null_value: float = -9999.25,
    company: str = "pycsamt",
) -> tuple[Path, LASExportReport]:
    """Write inline PCBH curves, or an interval-derived LAS subset."""
    borehole.validate()
    if not math.isfinite(null_value):
        raise ValueError("null_value must be finite")
    if not isinstance(company, str) or not company.strip():
        raise ValueError("company must be a non-empty string")
    curves = borehole.extensions.get("pcbh:continuous_curves", [])
    losses = [
        "CRS and absolute collar coordinates",
        "structures and non-lithology interval families",
        "samples, assays, construction, and PCBH provenance",
        "original LAS headers not represented by the subset writer",
    ]
    if borehole.trajectory.method != "vertical":
        losses.append("deviated trajectory survey")
    if curves:
        depth, definitions, columns = _extension_columns(curves, null_value)
    else:
        depth, definitions, columns = _interval_columns(borehole, null_value)
        losses.append(
            "continuous curves unavailable; values derived from intervals"
        )
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    step = _median_step(depth)
    lines = [
        "~VERSION INFORMATION",
        " VERS. 2.0: CWLS LOG ASCII STANDARD",
        " WRAP. NO: ONE LINE PER DEPTH STEP",
        "~WELL INFORMATION",
        f" STRT.M {depth[0]:.6g}: START DEPTH",
        f" STOP.M {depth[-1]:.6g}: STOP DEPTH",
        f" STEP.M {step:.6g}: STEP",
        f" NULL. {null_value:.12g}: NULL VALUE",
        f" COMP. {company}: COMPANY",
        f" WELL. {borehole.name}: WELL",
        "~CURVE INFORMATION",
    ]
    lines.extend(
        f" {name}.{unit}: {description}"
        for name, unit, description in definitions
    )
    lines.append("~A " + " ".join(name for name, _, _ in definitions))
    for index in range(len(depth)):
        lines.append(" ".join(f"{column[index]:.12g}" for column in columns))
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output, LASExportReport(
        path=str(output),
        curves_written=[item[0] for item in definitions],
        losses=losses,
    )


def _read_las(path: Path, *, max_samples: int):
    sections: dict[str, dict[str, str]] = {}
    curve_defs, data = [], []
    section = ""
    for raw in path.read_text(encoding="utf-8-sig").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("~"):
            section = line[1:2].upper()
            continue
        if section == "A":
            values = [float(value) for value in line.split()]
            if not all(math.isfinite(value) for value in values):
                raise ValueError("LAS ASCII data values must be finite")
            data.append(values)
            if len(data) > max_samples:
                raise ValueError("LAS sample limit exceeded")
        elif section == "C":
            mnemonic, unit, value, description = _las_item(line)
            curve_defs.append((mnemonic, unit, description))
        else:
            mnemonic, _, value, _ = _las_item(line)
            sections.setdefault(section, {})[mnemonic] = value
    if not curve_defs or not data:
        raise ValueError("LAS file requires curve definitions and ASCII data")
    if any(len(row) != len(curve_defs) for row in data):
        raise ValueError("LAS data width does not match curve definitions")
    if len({item[0] for item in curve_defs}) != len(curve_defs):
        raise ValueError("LAS curve mnemonics must be unique")
    curves = {
        name: {
            "unit": unit,
            "description": description,
            "column": index,
            "index": [row[index] for row in data],
        }
        for index, (name, unit, description) in enumerate(curve_defs)
    }
    return sections, curves, data


def _las_item(line: str):
    left, _, description = line.partition(":")
    mnemonic_unit, _, value = left.partition(" ")
    mnemonic, _, unit = mnemonic_unit.partition(".")
    return (
        mnemonic.strip().upper(),
        unit.strip(),
        value.strip(),
        description.strip(),
    )


def _well_text(sections, key):
    for section in ("W", "V"):
        if key in sections.get(section, {}):
            return sections[section][key]
    return ""


def _well_float(sections, key, default):
    try:
        return float(_well_text(sections, key))
    except (TypeError, ValueError):
        return default


def _median_step(depth):
    if len(depth) < 2:
        return 1.0
    differences = sorted(abs(b - a) for a, b in zip(depth, depth[1:]))
    return differences[len(differences) // 2]


def _las_intervals(index, rows, curves, lithology, resistivity, null, step):
    if not lithology or lithology.upper() not in curves:
        return [], []
    lith_col = curves[lithology.upper()]["column"]
    resistivity_definition = curves.get(resistivity.upper(), {})
    resistivity_unit = resistivity_definition.get("unit", "").upper()
    res_col = (
        resistivity_definition.get("column")
        if resistivity_unit in {"OHMM", "OHM.M", "OHM-M"}
        else None
    )
    intervals, vocabulary = [], {}
    start = 0
    for end in range(1, len(rows) + 1):
        if end < len(rows) and rows[end][lith_col] == rows[start][lith_col]:
            continue
        code = f"{rows[start][lith_col]:g}"
        values = (
            []
            if res_col is None
            else [
                row[res_col]
                for row in rows[start:end]
                if not math.isclose(row[res_col], null)
            ]
        )
        intervals.append(
            LogInterval(
                index[start],
                index[end - 1] + step,
                code=code,
                label=code,
                resistivity_ohm_m=sum(values) / len(values)
                if values
                else None,
                data_nature="observed",
            )
        )
        vocabulary.setdefault(code, VocabularyEntry(code, code))
        start = end
    return intervals, list(vocabulary.values())


def _extension_columns(curves, null):
    depth = [float(item[0]) for item in curves[0]["samples"]]
    definitions = [("DEPT", "M", "MEASURED DEPTH")]
    columns = [depth]
    for curve in curves:
        if [float(item[0]) for item in curve["samples"]] != depth:
            raise ValueError(
                "LAS export requires curves on one shared depth index"
            )
        definitions.append(
            (curve["mnemonic"], curve.get("unit", ""), curve.get("name", ""))
        )
        columns.append(
            [
                null if item[1] is None else float(item[1])
                for item in curve["samples"]
            ]
        )
    return depth, definitions, columns


def _interval_columns(borehole, null):
    intervals = borehole.interval_logs.get("lithology", [])
    if not intervals:
        raise ValueError(
            "LAS export needs inline curves or lithology intervals"
        )
    depth = [item.from_md for item in intervals]
    definitions = [
        ("DEPT", "M", "MEASURED DEPTH"),
        ("RESD", "OHMM", "RESISTIVITY"),
        ("LITH", "", "LITHOLOGY CODE"),
    ]
    codes = {
        item.code or item.label: index + 1
        for index, item in enumerate(intervals)
    }
    return (
        depth,
        definitions,
        [
            depth,
            [item.resistivity_ohm_m or null for item in intervals],
            [codes[item.code or item.label] for item in intervals],
        ],
    )
