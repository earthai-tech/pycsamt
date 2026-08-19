# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt info
============

Inspect EDI and EMTF XML transfer-function metadata.

Usage
-----
::

    pycsamt info station.edi
    pycsamt info station.xml
    pycsamt info survey/
    pycsamt info survey/ --format json
    pycsamt info survey/ --stations S01,S02 -v
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import click

from ...api.cli.config import configure_cli
from ...api.cli.options import (
    format_option,
    fresh_option,
    no_color_option,
    output_dir_option,
    survey_option,
    verbose_option,
)
from ...api.cli.params import StationList
from ..survey import SurveyContext, resolve_survey

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _collect_edi_paths(source: Path | None, sites: Any = None) -> list[Path]:
    """Return a sorted list of .edi paths from a path or a Sites object."""
    if source is not None:
        if source.is_dir():
            return sorted(source.glob("*.edi"))
        if source.suffix.lower() == ".edi":
            return [source]
        return []
    if sites is not None:
        paths = []
        for site in sites:
            try:
                p = Path(site.edi.path)
                if p.exists():
                    paths.append(p)
            except AttributeError:
                pass
        return sorted(paths)
    return []


def _edi_info(path: Path) -> dict[str, Any]:
    """Extract metadata from a single EDI file without external imports.

    Returns a plain dict so callers can format freely (text / json / csv).
    Falls back gracefully when the jones parser is unavailable.
    """
    record: dict[str, Any] = {
        "file": str(path),
        "station": None,
        "latitude": None,
        "longitude": None,
        "elevation": None,
        "n_frequencies": None,
        "freq_min_hz": None,
        "freq_max_hz": None,
        "components": [],
        "quality": {},
        "error": None,
    }

    try:
        # Try the jones fast-path for J files accidentally renamed to .edi;
        # for genuine EDI files use the lightweight header-only parser.
        _parse_edi_header(path, record)
    except Exception:  # noqa: BLE001
        try:
            _parse_edi_header(path, record)
        except Exception as exc2:  # noqa: BLE001
            record["error"] = str(exc2)

    return record


def _parse_edi_header(path: Path, record: dict[str, Any]) -> None:
    """Populate *record* from the [HEAD] and [FREQ] sections of an EDI file."""
    import numpy as np  # noqa: PLC0415

    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()

    # --- HEAD block ---
    in_head = False
    freq_lines: list[str] = []
    in_freq = False
    components: set[str] = set()

    for line in lines:
        stripped = line.strip()
        upper = stripped.upper()

        if upper.startswith(">HEAD"):
            in_head = True
            in_freq = False
            continue
        if upper.startswith(">FREQ"):
            in_head = False
            in_freq = True
            continue
        if upper.startswith(">") and not upper.startswith(">!"):
            in_head = False
            # detect component blocks: >ZXX, >ZXYR, >TIPPER, etc.
            token = upper.lstrip(">").split()[0] if upper.lstrip(">") else ""
            for comp in (
                "ZXX",
                "ZXY",
                "ZYX",
                "ZYY",
                "TXR",
                "TXI",
                "TYR",
                "TYI",
            ):
                if token.startswith(comp):
                    components.add(comp[:3])
            in_freq = False
            continue

        if in_head and "=" in stripped:
            key, _, val = stripped.partition("=")
            key = key.strip().upper()
            val = val.strip().strip('"').strip("'")
            if key == "DATAID":
                record["station"] = val
            elif key == "LAT":
                try:
                    record["latitude"] = float(val)
                except ValueError:
                    pass
            elif key == "LONG":
                try:
                    record["longitude"] = float(val)
                except ValueError:
                    pass
            elif key == "ELEV":
                try:
                    record["elevation"] = float(val)
                except ValueError:
                    pass

        if (
            in_freq
            and stripped
            and not stripped.startswith(">")
            and not stripped.startswith("!")
        ):
            freq_lines.append(stripped)

    # --- parse frequencies ---
    freqs: list[float] = []
    for fl in freq_lines:
        for token in fl.split():
            try:
                freqs.append(float(token))
            except ValueError:
                pass

    if freqs:
        arr = np.asarray(freqs)
        record["n_frequencies"] = int(arr.size)
        record["freq_min_hz"] = float(arr.min())
        record["freq_max_hz"] = float(arr.max())

    record["components"] = sorted(components)




def _emtf_info(path: Path) -> dict[str, Any]:
    """Extract a format-neutral scientific summary from EDI or EMTF XML."""
    import numpy as np  # noqa: PLC0415

    from ...emtf import EMTF
    from ...io import detect_tf_format

    fmt = detect_tf_format(path)
    doc = EMTF.from_edi(path) if fmt == "edi" else EMTF.from_xml(path)
    freq = (
        np.asarray(doc.frequency, dtype=float)
        if doc.frequency is not None
        else np.array([])
    )
    location = doc.site.location if doc.site is not None else None
    orientation = doc.orientation
    components = []
    covariance: dict[str, list[str]] = {}
    for key, tf in doc.transfer_functions.items():
        components.append(key)
        covariance[key] = sorted(tf.estimates)

    return {
        "file": str(path),
        "format": fmt,
        "station": doc.station,
        "latitude": None if location is None else location.latitude,
        "longitude": None if location is None else location.longitude,
        "elevation": None if location is None else location.elevation,
        "n_frequencies": int(freq.size),
        "freq_min_hz": float(freq.min()) if freq.size else None,
        "freq_max_hz": float(freq.max()) if freq.size else None,
        "components": sorted(components),
        "tags": list(doc.tags),
        "orientation": (
            None
            if orientation is None
            else {
                "mode": orientation.mode,
                "angle_to_geographic_north": (
                    orientation.angle_to_geographic_north
                ),
            }
        ),
        "covariance": covariance,
        "quality": {},
        "error": None,
    }


def _collect_tf_paths(source: Path) -> list[Path]:
    """Collect supported TF paths for explicit file/directory inspection."""
    if source.is_file():
        return [source] if source.suffix.lower() in {".edi", ".xml"} else []
    if source.is_dir():
        return sorted(
            p
            for p in source.iterdir()
            if p.is_file() and p.suffix.lower() in {".edi", ".xml"}
        )
    return []

# ---------------------------------------------------------------------------
# Formatters
# ---------------------------------------------------------------------------


def _fmt_text(records: list[dict[str, Any]], verbose: int) -> str:
    lines: list[str] = []
    for r in records:
        lines.append(f"Station : {r.get('station') or Path(r['file']).stem}")
        lines.append(f"  File       : {r['file']}")
        if r.get("format"):
            lines.append(f"  Format     : {r['format']}")
        if r.get("error"):
            lines.append(f"  ERROR      : {r['error']}")
            lines.append("")
            continue
        lines.append(
            f"  Lat / Lon  : {r.get('latitude')}, {r.get('longitude')}"
        )
        lines.append(f"  Elevation  : {r.get('elevation')} m")
        nfreq = r.get("n_frequencies")
        fmin = r.get("freq_min_hz")
        fmax = r.get("freq_max_hz")
        if nfreq and fmin is not None and fmax is not None:
            lines.append(
                f"  Frequencies: {nfreq}  [{fmin:.4g} – {fmax:.4g} Hz]"
            )
        else:
            lines.append(f"  Frequencies: {nfreq or 0}")
        lines.append(
            f"  Components : {', '.join(r.get('components', [])) or '—'}"
        )
        if r.get("orientation"):
            orient = r["orientation"]
            lines.append(
                "  Orientation: "
                f"{orient.get('mode')} @ "
                f"{orient.get('angle_to_geographic_north')}°"
            )
        if verbose >= 1 and r.get("covariance"):
            lines.append(f"  Estimates  : {r['covariance']}")
        if verbose >= 1 and r.get("quality"):
            lines.append(f"  Quality    : {r['quality']}")
        lines.append("")
    return "\n".join(lines)


def _fmt_json(records: list[dict[str, Any]]) -> str:
    return json.dumps(records, indent=2, default=str)


def _fmt_csv(records: list[dict[str, Any]]) -> str:
    if not records:
        return ""
    header = ",".join(records[0].keys())
    rows = [header]
    for r in records:
        rows.append(",".join(str(v) for v in r.values()))
    return "\n".join(rows)


# ---------------------------------------------------------------------------
# Command
# ---------------------------------------------------------------------------


@click.command("info")
@click.argument(
    "source",
    type=click.Path(exists=True, path_type=Path),
    metavar="EDI_FILE_OR_DIR",
    required=False,
    default=None,
)
@survey_option
@fresh_option
@click.option(
    "--stations",
    type=StationList(),
    default=None,
    metavar="S01,S02,…",
    help="Filter output to specific station IDs (comma-separated).",
)
@verbose_option
@no_color_option
@format_option
@output_dir_option
@click.pass_context
def info(
    ctx: click.Context,
    source: Path | None,
    survey_path: Path | None,
    fresh: bool,
    stations: list[str] | None,
    verbose: int,
    no_color: bool,
    output_format: str,
    output_dir: Path,
) -> None:
    """Inspect metadata for EDI or EMTF XML transfer functions.

    SOURCE is optional when an active survey is set with
    ``pycsamt survey set``.  Priority: SOURCE > --survey > active context.

    \b
    Examples:
      # with an active survey set:
      pycsamt info
      pycsamt info --stations S01,S02

      # explicit path:
      pycsamt info station.edi
      pycsamt info station.xml
      pycsamt info survey/
      pycsamt info survey/ --format json -v
    """
    configure_cli(
        log__level=verbose,
        log__color=not no_color,
        output__format=output_format,
        output__dir=output_dir,
    )

    explicit = source or survey_path

    # Explicit paths can now contain either EDI or EMTF XML. Active survey
    # context remains EDI-oriented for backward compatibility.
    if explicit is not None and not fresh:
        edi_paths = _collect_tf_paths(explicit)
    else:
        sites = resolve_survey(explicit, fresh=fresh, verbose=verbose)
        edi_paths = _collect_edi_paths(None, sites)
        if not edi_paths and explicit is not None:
            edi_paths = _collect_tf_paths(explicit)
        elif not edi_paths:
            ctx_obj = SurveyContext.load()
            if ctx_obj is not None:
                edi_paths = _collect_edi_paths(ctx_obj.survey_path)

    if not edi_paths:
        click.echo(
            "No supported EDI/XML transfer-function files found.", err=True
        )
        sys.exit(1)

    if stations:
        station_set = {s.upper() for s in stations}
        edi_paths = [p for p in edi_paths if p.stem.upper() in station_set]
        if not edi_paths:
            click.echo(
                f"No files matched stations: {', '.join(stations)}", err=True
            )
            sys.exit(1)

    if verbose >= 1:
        click.echo(f"Scanning {len(edi_paths)} TF file(s)…", err=True)

    records = [
        _emtf_info(p) if p.suffix.lower() == ".xml" else _edi_info(p)
        for p in edi_paths
    ]

    if output_format == "json":
        output = _fmt_json(records)
    elif output_format == "csv":
        output = _fmt_csv(records)
    else:
        output = _fmt_text(records, verbose)

    click.echo(output)
