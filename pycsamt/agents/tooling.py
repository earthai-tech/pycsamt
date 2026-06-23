# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.agents.tooling
======================

:class:`ToolAgent` — analysis and data/IO utilities from the web app's tools
menu, exposed as Agent-Master tasks that return a compact **table** (rendered
as a monospaced block in chat) plus optional **figures**:

* ``strike``          — geoelectric strike per station + rose/analysis figure,
* ``dimensionality``  — 1-D / 2-D / 3-D classification per station × period,
* ``validator``       — per-station EDI quality checklist,
* ``coords``          — transform station lat/lon to UTM easting/northing,
* ``elevation``       — enrich stations with elevation from an open web API,
* ``converter``       — re-export the survey to CSV / JSON / EDI on disk,
* ``batch_export``    — render a bundle of standard plots and save them.

The ``coords`` tool is read-only. The ``elevation`` tool reaches an external
network service, and ``converter`` / ``batch_export`` write files to a
user-supplied folder — these are gated by the parameter modal (the user picks
the API / output folder and clicks *Run* before anything outward-facing
happens) and report exactly what they did.

Like :class:`~pycsamt.agents.plotting.PlotAgent` it never calls an LLM; it
turns a validated parameter set into a result the Agent Master can display.
"""

from __future__ import annotations

import os
import time
from typing import Any

from ._base import AgentResult
from .plotting import _as_list, _filter_sites, _period_range, _truthy

__all__ = ["ToolAgent", "TOOL_KINDS"]

TOOL_KINDS = (
    "strike", "dimensionality", "validator",
    "coords", "elevation", "converter", "batch_export",
)

# Plot bundles offered by the batch_export tool. Values are PlotAgent kinds.
_EXPORT_BUNDLES: dict[str, tuple[str, ...]] = {
    "overview":     ("rhophi", "phase_psection", "pt_psection"),
    "phase_tensor": ("pt_psection", "pt_map"),
    "all":          ("rhophi", "phase_psection", "pt_psection", "pt_map"),
    "rhophi":       ("rhophi",),
    "phase_psection": ("phase_psection",),
    "pt_psection":  ("pt_psection",),
    "pt_map":       ("pt_map",),
}


def _as_fig(obj):
    """Return the matplotlib Figure for a function that returned a Figure
    or an Axes."""
    if obj is None:
        return None
    if hasattr(obj, "savefig"):          # already a Figure
        return obj
    getf = getattr(obj, "get_figure", None)
    if callable(getf):
        return getf()
    return getattr(obj, "figure", None)


def _df_to_text(df, columns=None, max_rows: int = 30, ndigits: int = 2) -> str:
    """Compact fixed-width rendering of a DataFrame for a chat code block."""
    d = df
    if columns:
        keep = [c for c in columns if c in d.columns]
        if keep:
            d = d[keep]
    d = d.copy()
    for c in d.columns:
        try:
            if d[c].dtype.kind == "f":
                d[c] = d[c].round(ndigits)
        except Exception:  # noqa: BLE001
            pass
    n = len(d)
    txt = d.head(max_rows).to_string(index=False)
    if n > max_rows:
        txt += f"\n… ({n - max_rows} more row(s))"
    return txt


def _circular_strike_mean(ang_deg) -> float:
    """Mean of strike angles modulo 180° (geoelectric strike ambiguity)."""
    import numpy as np
    a = np.asarray(ang_deg, float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return float("nan")
    m = np.angle(np.nanmean(np.exp(1j * np.deg2rad(a * 2.0))))
    return float(np.rad2deg(m) / 2.0)


def _get_latlon(ed) -> tuple:
    """Return ``(lat, lon)`` floats for an EDI-like object, or ``(None, None)``.

    Checks the object and its ``.edi`` wrapper for ``lat``/``latitude`` and
    ``lon``/``longitude`` attributes (the two naming conventions used across
    the EDI readers)."""
    for obj in (ed, getattr(ed, "edi", None)):
        if obj is None:
            continue
        lat = getattr(obj, "lat", None)
        if lat is None:
            lat = getattr(obj, "latitude", None)
        lon = getattr(obj, "lon", None)
        if lon is None:
            lon = getattr(obj, "longitude", None)
        if lat is not None and lon is not None:
            try:
                return float(lat), float(lon)
            except (TypeError, ValueError):
                continue
    return None, None


def _station_coords(sites) -> list:
    """Return ``[(name, lat, lon), …]`` for every station in *sites*.

    ``lat`` / ``lon`` are ``None`` when a station carries no coordinates."""
    from ..emtools._core import _iter_items, _name, _unwrap
    out = []
    for i, ed in enumerate(_iter_items(sites)):
        name = _name(ed, i)
        try:
            raw = _unwrap(ed)
        except Exception:  # noqa: BLE001
            raw = ed
        lat, lon = _get_latlon(raw)
        if lat is None and lon is None:
            lat, lon = _get_latlon(ed)
        out.append((name, lat, lon))
    return out


def _ll_to_utm(lat: float, lon: float, zone, hem: str, datum: str):
    """Return ``(easting, northing, zone)`` via pyproj, with a pure-Python
    fallback to :func:`pycsamt.gis.utils.ll_to_utm`."""
    try:
        from pyproj import Proj
        if zone is None:
            zone = int((lon + 180) / 6) + 1
        proj = Proj(proj="utm", zone=zone, datum=datum,
                    south=(hem == "S"), ellps=datum)
        e, n = proj(lon, lat)
        return e, n, zone
    except Exception:  # noqa: BLE001
        from ..gis.utils import ll_to_utm
        res = ll_to_utm(lat, lon)
        return (res["easting"], res["northing"],
                res.get("zone_number", zone or 0))


def _default_out_dir(name: str) -> str:
    """Default output folder under the user's home directory."""
    return os.path.join(os.path.expanduser("~"), name)


def _safe_filename(label: str) -> str:
    """Filesystem-safe figure label."""
    safe = "".join(c if c.isalnum() or c in "-_ " else "_" for c in label)
    return safe.strip().replace(" ", "_") or "figure"


class ToolAgent:
    """Strike / dimensionality / validator analysis tasks."""

    def __init__(self, **_: Any) -> None:
        self._last_cost = 0.0

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        t0 = time.time()
        self._last_cost = 0.0

        kind = str(input_data.get("kind", "")).strip()
        if kind not in TOOL_KINDS:
            return AgentResult.failed(
                f"Unknown tool {kind!r}. Expected one of {TOOL_KINDS}.",
                elapsed=time.time() - t0,
            )
        src = (
            input_data.get("sites")
            or input_data.get("path")
            or input_data.get("data_path")
        )
        if src is None:
            return AgentResult.failed(
                "No data supplied. Load an EDI dataset first.",
                elapsed=time.time() - t0,
            )

        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # noqa: F401

        from ..emtools._core import ensure_sites
        try:
            sites = ensure_sites(src, recursive=True, strict=False, verbose=0)
        except Exception as exc:  # noqa: BLE001
            return AgentResult.failed(
                f"Could not load EDI data: {exc}",
                elapsed=time.time() - t0,
            )

        stations = _as_list(input_data.get("stations"))
        sub = _filter_sites(sites, stations) if stations else sites
        warnings: list[str] = []

        try:
            if kind == "strike":
                summary, table, figs = self._strike(sub, input_data, warnings)
            elif kind == "dimensionality":
                summary, table, figs = self._dimensionality(
                    sub, input_data, warnings
                )
            elif kind == "coords":
                summary, table, figs = self._coords(
                    sub, input_data, warnings
                )
            elif kind == "elevation":
                summary, table, figs = self._elevation(
                    sub, input_data, warnings
                )
            elif kind == "converter":
                summary, table, figs = self._converter(
                    sub, input_data, warnings
                )
            elif kind == "batch_export":
                summary, table, figs = self._batch_export(
                    sub, input_data, warnings
                )
            else:  # validator
                summary, table, figs = self._validator(
                    sub, input_data, warnings
                )
        except Exception as exc:  # noqa: BLE001
            return AgentResult.failed(
                f"{kind} analysis failed: {exc}",
                elapsed=time.time() - t0,
            )

        return AgentResult(
            status="success",
            summary=summary,
            data={"table_text": table, "figures": figs, "tool_kind": kind},
            warnings=warnings,
            elapsed_seconds=time.time() - t0,
            cost_estimate_usd=0.0,
        )

    # ── tools ────────────────────────────────────────────────────────────────
    def _strike(self, sites, d, warnings):
        from ..emtools.strike import (
            estimate_strike_sweep,
            estimate_strike_phase_tensor,
            estimate_strike_consensus,
            plot_strike_analysis,
        )
        method = str(d.get("method", "consensus") or "consensus").lower()
        band = _period_range(d)
        fn = {
            "sweep": estimate_strike_sweep,
            "pt": estimate_strike_phase_tensor,
            "consensus": estimate_strike_consensus,
        }.get(method, estimate_strike_consensus)
        res = fn(sites, band=band, verbose=0)
        df = res.frame if hasattr(res, "frame") else res

        regional = _circular_strike_mean(df["ang"]) if "ang" in df else float("nan")
        summary = (
            f"**Geoelectric strike ({method})** — regional ≈ "
            f"{regional:.1f}° (N{regional:+.0f}°E) across {len(df)} station(s). "
            "Note the inherent 90° ambiguity."
        )
        table = _df_to_text(df, columns=["station", "ang", "iqr", "n"])

        figs = {}
        try:
            fmethod = "pt" if method == "pt" else "sweep"
            fig = _as_fig(
                plot_strike_analysis(sites, method=fmethod, band=band, verbose=0)
            )
            if fig is not None:
                figs["Strike analysis"] = fig
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"strike figure skipped: {exc}")
        return summary, table, figs

    def _dimensionality(self, sites, d, warnings):
        from ..emtools.dimensionality import classify_dimensionality
        try:
            skew_th = float(d.get("skew_th", 3.0) or 3.0)
        except (TypeError, ValueError):
            skew_th = 3.0
        try:
            ellipt_th = float(d.get("ellipt_th", 0.2) or 0.2)
        except (TypeError, ValueError):
            ellipt_th = 0.2

        res = classify_dimensionality(
            sites, skew_th=skew_th, ellipt_th=ellipt_th, verbose=0
        )
        df = res.frame if hasattr(res, "frame") else res
        import numpy as np
        import pandas as pd

        _LBL = {0: "indet.", 1: "1-D", 2: "2-D", 3: "3-D"}
        counts = df["dim"].value_counts().to_dict() if "dim" in df else {}
        total = int(sum(counts.values())) or 1
        dist = ", ".join(
            f"{_LBL.get(k, k)}: {v} ({100 * v / total:.0f}%)"
            for k, v in sorted(counts.items())
        )
        summary = (
            f"**Dimensionality** (skew≤{skew_th:g}, ellipticity≤{ellipt_th:g}) "
            f"over {total} station×period cells — {dist}."
        )

        # per-station dominant class
        if "station" in df and "dim" in df:
            agg = (
                df.groupby("station")["dim"]
                .agg(lambda s: int(s.mode().iloc[0]) if not s.mode().empty else 0)
                .reset_index()
            )
            agg["class"] = agg["dim"].map(_LBL)
            table = _df_to_text(agg, columns=["station", "class"], max_rows=40)
        else:
            table = _df_to_text(df, max_rows=20)

        figs = {}
        try:
            from ..emtools.tensor import plot_dimensionality_psection
            fig = _as_fig(plot_dimensionality_psection(sites, verbose=0))
            if fig is not None:
                figs["Dimensionality pseudo-section"] = fig
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"dimensionality figure skipped: {exc}")
        return summary, table, figs

    def _validator(self, sites, d, warnings):
        from ..agents.loader import _quality_scan
        import pandas as pd
        rows, scan_warn = _quality_scan(sites)
        if not rows:
            return (
                "No stations with valid impedance data found.",
                "(empty)",
                {},
            )
        df = pd.DataFrame(rows)
        # a station is "flagged" if it lacks Z / coords or scores low
        def _flag(r):
            issues = []
            if not r.get("has_z"):
                issues.append("no-Z")
            if not r.get("has_coords"):
                issues.append("no-coords")
            if (r.get("qc_score") or 0) < 50:
                issues.append("low-QC")
            return ",".join(issues) or "ok"
        df["flags"] = df.apply(_flag, axis=1)
        n_flag = int((df["flags"] != "ok").sum())
        summary = (
            f"**EDI validation** — {len(df)} station(s), {n_flag} flagged"
            f"{(' (' + str(len(scan_warn)) + ' warning(s))') if scan_warn else ''}."
        )
        table = _df_to_text(
            df,
            columns=["station", "has_z", "has_coords", "n_freq",
                     "qc_score", "flags"],
            max_rows=40,
        )
        return summary, table, {}

    # ── data / IO tools (Wave C) ──────────────────────────────────────────────
    def _coords(self, sites, d, warnings):
        """Transform every station's lat/lon to UTM easting/northing."""
        import numpy as np
        import pandas as pd

        datum = str(d.get("datum", "WGS84") or "WGS84")
        try:
            zone_in = int(float(d.get("zone", 0) or 0))
        except (TypeError, ValueError):
            zone_in = 0
        zone = zone_in if zone_in > 0 else None
        hem = (str(d.get("hemisphere", "N") or "N").upper() + "N")[0]
        hem = "S" if hem == "S" else "N"

        rows = []
        n_ok = 0
        for name, lat, lon in _station_coords(sites):
            rec = {
                "station": name, "lat": lat, "lon": lon,
                "easting": None, "northing": None, "zone": None,
            }
            if lat is not None and lon is not None:
                try:
                    e, n, z = _ll_to_utm(lat, lon, zone, hem, datum)
                    rec["easting"] = round(float(e), 1)
                    rec["northing"] = round(float(n), 1)
                    rec["zone"] = f"{int(z)}{hem}"
                    n_ok += 1
                except Exception as exc:  # noqa: BLE001
                    warnings.append(f"{name}: transform failed ({exc})")
            rows.append(rec)

        if not rows:
            return "No stations found.", "(empty)", {}
        df = pd.DataFrame(rows)
        n_total = len(df)
        n_missing = n_total - n_ok
        zlabel = f"zone {zone}{hem}" if zone else "auto zone"
        summary = (
            f"**Coordinate transform** (lat/lon → UTM, {datum}, {zlabel}) — "
            f"{n_ok}/{n_total} station(s) converted"
            + (f"; {n_missing} without coordinates." if n_missing else ".")
        )
        table = _df_to_text(
            df,
            columns=["station", "lat", "lon", "easting", "northing", "zone"],
            max_rows=60, ndigits=6,
        )
        return summary, table, {}

    def _elevation(self, sites, d, warnings):
        """Fetch elevation for stations with coordinates via an open web API."""
        import numpy as np
        import pandas as pd

        api = str(d.get("api", "open_meteo") or "open_meteo")
        coords = _station_coords(sites)
        with_coords = [
            (n, la, lo) for n, la, lo in coords
            if la is not None and lo is not None
        ]
        if not with_coords:
            return (
                "No stations carry coordinates — cannot fetch elevation.",
                "(empty)", {},
            )

        # One batched request (the API accepts arrays); fall back to NaN on
        # any network/library failure so the table still renders.
        elev_map: dict[str, float] = {}
        try:
            from ..gis.utils import get_elevation_from_api
            lats = np.array([la for _, la, lo in with_coords], dtype=float)
            lons = np.array([lo for _, la, lo in with_coords], dtype=float)
            res = get_elevation_from_api(lats, lons, api_name=api)
            arr = np.atleast_1d(np.asarray(res, dtype=float))
            for (n, _, _), ev in zip(with_coords, arr):
                elev_map[n] = float(ev)
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"elevation API '{api}' failed: {exc}")

        rows = []
        for name, lat, lon in coords:
            ev = elev_map.get(name, float("nan"))
            rows.append({
                "station": name, "lat": lat, "lon": lon,
                "elevation_m": round(ev, 1) if np.isfinite(ev) else None,
            })
        df = pd.DataFrame(rows)
        n_ok = int(df["elevation_m"].notna().sum())
        summary = (
            f"**Elevation enrichment** via external API ({api}) — "
            f"{n_ok}/{len(df)} station(s) resolved. "
            "Queried an open elevation web service over the network."
        )
        table = _df_to_text(
            df, columns=["station", "lat", "lon", "elevation_m"],
            max_rows=60, ndigits=6,
        )
        return summary, table, {}

    def _converter(self, sites, d, warnings):
        """Re-export the survey metadata (and optionally EDIs) to a folder."""
        import numpy as np
        import pandas as pd
        from ..emtools._core import _get_z_block, _iter_items, _name, _unwrap

        fmt = str(d.get("format", "csv") or "csv").lower()
        if fmt not in ("csv", "json", "edi"):
            fmt = "csv"
        out_dir = (str(d.get("output_dir") or "").strip()
                   or _default_out_dir("pycsamt_export"))

        rows = []
        n_edi = 0
        for i, ed in enumerate(_iter_items(sites)):
            try:
                raw = _unwrap(ed)
            except Exception:  # noqa: BLE001
                raw = ed
            name = _name(ed, i)
            lat, lon = _get_latlon(raw)
            if lat is None and lon is None:
                lat, lon = _get_latlon(ed)
            n_freq = 0
            t_min = t_max = float("nan")
            has_err = False
            try:
                _, z, fr, ze = _get_z_block(raw, with_errors=True)
                if fr is not None and len(fr) > 0:
                    fa = np.asarray(fr, dtype=float)
                    fa = fa[fa > 0]
                    n_freq = int(fa.size)
                    if n_freq:
                        periods = 1.0 / fa
                        t_min, t_max = float(periods.min()), float(periods.max())
                if ze is not None:
                    has_err = bool(np.any(np.isfinite(np.asarray(ze))))
            except Exception:  # noqa: BLE001
                pass
            rows.append({
                "station": name,
                "lat": lat, "lon": lon,
                "n_freq": n_freq,
                "t_min": round(t_min, 6) if np.isfinite(t_min) else None,
                "t_max": round(t_max, 6) if np.isfinite(t_max) else None,
                "has_z_err": has_err,
            })

        if not rows:
            return "No stations to convert.", "(empty)", {}

        os.makedirs(out_dir, exist_ok=True)
        written = []
        try:
            if fmt == "csv":
                import csv
                path = os.path.join(out_dir, "survey_stations.csv")
                with open(path, "w", newline="", encoding="utf-8") as fh:
                    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
                    w.writeheader()
                    w.writerows(rows)
                written.append(path)
            elif fmt == "json":
                import json
                path = os.path.join(out_dir, "survey_stations.json")
                with open(path, "w", encoding="utf-8") as fh:
                    json.dump(rows, fh, indent=2)
                written.append(path)
            else:  # edi re-export, best-effort per station
                for ed in _iter_items(sites):
                    try:
                        edi_obj = getattr(ed, "edi", None) or _unwrap(ed)
                        write_fn = getattr(edi_obj, "write_edi_file", None)
                        if write_fn is None:
                            continue
                        sname = _name(ed, n_edi)
                        epath = os.path.join(out_dir, f"{sname}.edi")
                        write_fn(epath)
                        written.append(epath)
                        n_edi += 1
                    except Exception as exc:  # noqa: BLE001
                        warnings.append(f"EDI write skipped: {exc}")
        except Exception as exc:  # noqa: BLE001
            raise RuntimeError(f"write failed: {exc}") from exc

        if fmt == "edi" and not written:
            warnings.append(
                "no EDI writer available on the loaded objects; "
                "wrote nothing (try CSV/JSON)."
            )
        df = pd.DataFrame(rows)
        what = (f"{len(written)} EDI file(s)" if fmt == "edi"
                else f"{fmt.upper()} ({len(rows)} stations)")
        summary = (
            f"**Format conversion** → wrote {what} to `{out_dir}`."
        )
        table = _df_to_text(
            df,
            columns=["station", "lat", "lon", "n_freq",
                     "t_min", "t_max", "has_z_err"],
            max_rows=40, ndigits=6,
        )
        return summary, table, {}

    def _batch_export(self, sites, d, warnings):
        """Render a bundle of standard plots and save them to a folder."""
        import matplotlib.pyplot as plt
        import pandas as pd
        from .plotting import PlotAgent

        bundle = str(d.get("plots", "overview") or "overview").lower()
        kinds = _EXPORT_BUNDLES.get(bundle, _EXPORT_BUNDLES["overview"])
        fmt = str(d.get("format", "png") or "png").lower().lstrip(".")
        try:
            dpi = int(float(d.get("dpi", 150) or 150))
        except (TypeError, ValueError):
            dpi = 150
        dpi = max(72, min(600, dpi))
        out_dir = (str(d.get("output_dir") or "").strip()
                   or _default_out_dir("pycsamt_figures"))
        os.makedirs(out_dir, exist_ok=True)

        figs_out: dict = {}
        rows = []
        for kind in kinds:
            try:
                res = PlotAgent().execute({
                    "sites": sites, "kind": kind, "publication": "on",
                })
            except Exception as exc:  # noqa: BLE001
                warnings.append(f"{kind}: render failed ({exc})")
                continue
            if res.status != "success":
                warnings.append(f"{kind}: {res.summary}")
                continue
            for title, fig in (res.data.get("figures") or {}).items():
                if not hasattr(fig, "savefig"):
                    continue
                fname = f"{_safe_filename(title)}.{fmt}"
                path = os.path.join(out_dir, fname)
                try:
                    fig.savefig(path, dpi=dpi, bbox_inches="tight")
                    figs_out[title] = fig
                    rows.append({"plot": kind, "file": fname, "saved": "ok"})
                except Exception as exc:  # noqa: BLE001
                    warnings.append(f"{title}: save failed ({exc})")
                    plt.close(fig)

        if not rows:
            return (
                "No figures could be rendered for the selected bundle.",
                "(empty)", {},
            )
        df = pd.DataFrame(rows)
        summary = (
            f"**Batch plot export** ({bundle}) — saved {len(rows)} figure(s) "
            f"as {fmt.upper()} @ {dpi} dpi to `{out_dir}`."
        )
        table = _df_to_text(df, columns=["plot", "file", "saved"], max_rows=40)
        return summary, table, figs_out
