# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt map - root Click group and shared helpers.

The map command group collects station-location and map-display workflows
that are useful from the terminal.  Submodules import ``map`` from here
and register commands via ``@map.command(...)``.

Shared utilities
----------------
_get_sites        Resolve positional path / --survey / active context.
_station_rows     Extract station coordinate rows from a Sites object.
_format_rows      Render coordinate rows as text, JSON, or CSV.
"""

from __future__ import annotations

import csv
import io
import json
import math
from pathlib import Path
from typing import Any

import click


def _get_sites(
    explicit: Path | None,
    survey_path: Path | None,
    fresh: bool,
    verbose: int,
) -> Any:
    """Resolve an EDI source to a Sites object.

    Priority is explicit path, then ``--survey``, then active context.
    This mirrors ``site`` and ``invert build`` so map commands behave
    consistently with the rest of the CLI.
    """
    from pycsamt.cli.survey import (
        resolve_survey,  # noqa: PLC0415
    )

    return resolve_survey(
        explicit or survey_path, fresh=fresh, verbose=verbose
    )


def _finite_float(value: Any) -> float | None:
    """Return a finite float or None for missing/non-numeric values."""
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def _site_summary(site: Any) -> dict[str, Any]:
    """Return a tolerant one-site summary dictionary."""
    try:
        summary = site.summary()
        if isinstance(summary, dict):
            return summary
    except Exception:  # noqa: BLE001
        pass

    lat = lon = elev = None
    try:
        coords = site.coords
        lat, lon, elev = coords[0], coords[1], coords[2]
    except Exception:  # noqa: BLE001
        pass

    return {
        "name": getattr(site, "name", None)
        or getattr(site, "station", "site"),
        "lat": lat,
        "lon": lon,
        "elev": elev,
        "nfreq": getattr(site, "nfreq", None),
    }


def _station_rows(
    sites: Any, *, drop_missing: bool = False
) -> list[dict[str, Any]]:
    """Extract station coordinates from *sites* as JSON-ready rows."""
    api_rows = _station_rows_from_map_api(sites, drop_missing=drop_missing)
    if api_rows:
        return api_rows
    if api_rows == [] and _is_empty(sites):
        return []

    rows: list[dict[str, Any]] = []
    for idx, site in enumerate(sites, start=1):
        summary = _site_summary(site)
        lat = _finite_float(summary.get("lat"))
        lon = _finite_float(summary.get("lon"))
        elev = _finite_float(summary.get("elev"))
        if drop_missing and (lat is None or lon is None):
            continue
        rows.append(
            {
                "index": idx,
                "station": str(summary.get("name") or f"S{idx:03d}"),
                "lat": lat,
                "lon": lon,
                "elev": elev,
                "nfreq": summary.get("nfreq"),
            }
        )
    return rows


def _is_empty(value: Any) -> bool:
    """Return True when *value* is a sized empty container."""
    try:
        return len(value) == 0
    except Exception:  # noqa: BLE001
        return False


def _station_rows_from_map_api(
    sites: Any,
    *,
    drop_missing: bool,
) -> list[dict[str, Any]] | None:
    """Return station rows through ``pycsamt.map`` when possible."""
    try:
        from pycsamt.map import (
            ensure_map_data,  # noqa: PLC0415
        )
        from pycsamt.map._core import (
            station_dataframe,  # noqa: PLC0415
        )

        data = ensure_map_data(sites)
        df = station_dataframe(data)
    except Exception:  # noqa: BLE001
        return None

    rows: list[dict[str, Any]] = []
    for idx, record in enumerate(df.to_dict(orient="records"), start=1):
        lat = _finite_float(record.get("Latitude"))
        lon = _finite_float(record.get("Longitude"))
        elev = _finite_float(record.get("Elevation"))
        if drop_missing and (lat is None or lon is None):
            continue
        rows.append(
            {
                "index": int(record.get("Index", idx)),
                "station": str(record.get("ID") or f"S{idx:03d}"),
                "lat": lat,
                "lon": lon,
                "elev": elev,
                "nfreq": record.get("nfreq"),
            }
        )
    return rows


def _format_rows(rows: list[dict[str, Any]], output_format: str) -> str:
    """Render station rows in text, JSON, or CSV format."""
    fmt = output_format.lower()
    if fmt == "json":
        return json.dumps(rows, indent=2, default=str)

    if fmt == "csv":
        buf = io.StringIO()
        writer = csv.DictWriter(
            buf,
            fieldnames=["index", "station", "lat", "lon", "elev", "nfreq"],
        )
        writer.writeheader()
        writer.writerows(rows)
        return buf.getvalue().rstrip()

    if not rows:
        return "No stations found."

    header = f"{'idx':>4}  {'station':<18}  {'lat':>12}  {'lon':>12}  {'elev':>10}  {'nfreq':>7}"
    sep = "-" * len(header)
    lines = [header, sep]
    for row in rows:
        lat = _format_number(row["lat"], 6)
        lon = _format_number(row["lon"], 6)
        elev = _format_number(row["elev"], 2)
        nfreq = "" if row["nfreq"] is None else str(row["nfreq"])
        lines.append(
            f"{row['index']:>4}  {row['station']:<18.18}  "
            f"{lat:>12}  {lon:>12}  {elev:>10}  {nfreq:>7}"
        )
    return "\n".join(lines)


def _format_number(value: Any, ndigits: int) -> str:
    """Format a nullable number for text tables."""
    if value is None:
        return ""
    try:
        return f"{float(value):.{ndigits}f}"
    except (TypeError, ValueError):
        return str(value)


@click.group("map")
@click.pass_context
def map(ctx: click.Context) -> None:  # noqa: A001
    """Inspect and plot station map geometry.

    \b
    Typical workflow:
      pycsamt survey set data/AMT/WILLY_DATA/L18PLT
      pycsamt map stations
      pycsamt map stations --format csv --output station_coords.csv
      pycsamt map plot --output station_map.png
    """
    ctx.ensure_object(dict)
