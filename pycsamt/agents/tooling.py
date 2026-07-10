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
* ``batch_export``    — render a bundle of standard plots and save them,
* ``freq_editor``     — confidence-based frequency QC (drop / mask / recover),
* ``layered_model``   — build & preview a synthetic 1-D resistivity model.

``freq_editor`` mutates the survey: it runs out-of-place and hands the edited
``Sites`` back through ``AgentResult.data["corrected_sites"]`` so the chat's
post-processing modal can apply it to the session or export it — the same
pathway used by the static-shift / denoise correction workflows.
``layered_model`` is synthetic and needs no loaded data.

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
from .plotting import (
    _as_list,
    _filter_sites,
    _period_range,
)

__all__ = ["ToolAgent", "TOOL_KINDS"]

TOOL_KINDS = (
    "strike", "dimensionality", "validator",
    "coords", "elevation", "converter", "batch_export",
    "freq_editor", "layered_model", "correction",
)

# Tools that do not need a loaded EDI dataset.
_DATALESS_KINDS = ("layered_model",)

# Plot bundles offered by the batch_export tool. Values are PlotAgent kinds.
_EXPORT_BUNDLES: dict[str, tuple[str, ...]] = {
    "overview":     ("rhophi", "phase_psection", "pt_psection"),
    "phase_tensor": ("pt_psection", "pt_map", "pt_strip_grid"),
    "all":          ("rhophi", "phase_psection", "pt_psection", "pt_map",
                      "pt_strip_grid"),
    "rhophi":       ("rhophi",),
    "phase_psection": ("phase_psection",),
    "pt_psection":  ("pt_psection",),
    "pt_map":       ("pt_map",),
    "pt_strip_grid": ("pt_strip_grid",),
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


def _corr_coord_figure(raw_sites, corrected_sites, label: str):
    """Scatter of station positions before vs after a coordinate correction.

    ρ_a curves are unchanged by coordinate corrections, so a position map is
    the informative before/after view. Returns a Figure or ``None``.
    """
    import matplotlib.pyplot as plt

    from ..gis.coord_correction import _get_coords_df

    try:
        df0 = _get_coords_df(raw_sites)
        df1 = _get_coords_df(corrected_sites)
    except Exception:  # noqa: BLE001
        return None
    if df0 is None or df1 is None or df0.empty:
        return None

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(df0["lon"], df0["lat"], s=28, facecolors="none",
               edgecolors="#3498db", label="Before", zorder=3)
    ax.scatter(df1["lon"], df1["lat"], s=14, color="#e74c3c",
               label="After", zorder=4)
    # connect each station's old → new position
    for (_, r0), (_, r1) in zip(df0.iterrows(), df1.iterrows()):
        ax.plot([r0["lon"], r1["lon"]], [r0["lat"], r1["lat"]],
                color="#999", lw=0.6, zorder=2)
    ax.set_xlabel("Longitude", fontsize=8)
    ax.set_ylabel("Latitude", fontsize=8)
    ax.set_title(f"{label} — station positions", fontsize=9, fontweight="bold")
    ax.legend(fontsize=8)
    ax.tick_params(labelsize=7)
    fig.tight_layout()
    return fig


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
        self._corrected = None  # edited Sites stashed by mutating tools

        kind = str(input_data.get("kind", "")).strip()
        if kind not in TOOL_KINDS:
            return AgentResult.failed(
                f"Unknown tool {kind!r}. Expected one of {TOOL_KINDS}.",
                elapsed=time.time() - t0,
            )

        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # noqa: F401

        # Synthetic tools need no data — handle before the EDI guard.
        if kind in _DATALESS_KINDS:
            warnings: list[str] = []
            try:
                summary, table, figs = self._layered_model(input_data, warnings)
            except Exception as exc:  # noqa: BLE001
                return AgentResult.failed(
                    f"{kind} failed: {exc}", elapsed=time.time() - t0,
                )
            return AgentResult(
                status="success", summary=summary,
                data={"table_text": table, "figures": figs, "tool_kind": kind},
                warnings=warnings, elapsed_seconds=time.time() - t0,
                cost_estimate_usd=0.0,
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
            elif kind == "freq_editor":
                summary, table, figs = self._freq_editor(
                    sub, input_data, warnings
                )
            elif kind == "correction":
                summary, table, figs = self._correction(
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

        data: dict[str, Any] = {
            "table_text": table, "figures": figs, "tool_kind": kind,
        }
        if self._corrected is not None:
            data["corrected_sites"] = self._corrected
        return AgentResult(
            status="success",
            summary=summary,
            data=data,
            warnings=warnings,
            elapsed_seconds=time.time() - t0,
            cost_estimate_usd=0.0,
        )

    # ── tools ────────────────────────────────────────────────────────────────
    def _strike(self, sites, d, warnings):
        from ..emtools.strike import (
            estimate_strike_consensus,
            estimate_strike_phase_tensor,
            estimate_strike_sweep,
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
        from ..emtools.dimensionality import (
            classify_dimensionality,
        )
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
            from ..emtools.tensor import (
                plot_dimensionality_psection,
            )
            fig = _as_fig(plot_dimensionality_psection(sites, verbose=0))
            if fig is not None:
                figs["Dimensionality pseudo-section"] = fig
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"dimensionality figure skipped: {exc}")
        return summary, table, figs

    def _validator(self, sites, d, warnings):
        import pandas as pd

        from ..agents.loader import _quality_scan
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

        from ..emtools._core import (
            _get_z_block,
            _iter_items,
            _name,
            _unwrap,
        )

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

    # ── stateful tools (Wave D) ───────────────────────────────────────────────
    def _freq_editor(self, sites, d, warnings):
        """Confidence-based frequency QC. Edits out-of-place and stashes the
        edited Sites in ``self._corrected`` for the post-processing modal."""
        from ..emtools.frequency import (
            edit_frequencies_by_confidence,
            plot_frequency_edit_summary,
        )

        mode = str(d.get("mode", "recover") or "recover").lower()
        if mode not in ("recover", "drop", "mask"):
            mode = "recover"
        method = str(d.get("method", "composite") or "composite").lower()
        also = str(d.get("also", "both") or "both").lower()
        reject = str(d.get("reject", "drop") or "drop").lower()

        def _f(key, default):
            try:
                return float(d.get(key, default) or default)
            except (TypeError, ValueError):
                return default
        threshold = _f("threshold", 0.50)
        ci_hi = _f("ci_hi", 0.90)
        ci_lo = _f("ci_lo", 0.50)

        result = edit_frequencies_by_confidence(
            sites, mode=mode, method=method, threshold=threshold,
            ci_hi=ci_hi, ci_lo=ci_lo, interpolation="linear",
            reject=reject, also=also, inplace=False, verbose=0,
        )
        edited = getattr(result, "sites", None)
        self._corrected = edited

        decisions = getattr(result, "decisions", None)
        if hasattr(decisions, "frame"):
            decisions = decisions.frame
        elif hasattr(decisions, "_frame"):
            decisions = decisions._frame

        n_drop = int(getattr(result, "n_dropped", 0) or 0)
        n_mask = int(getattr(result, "n_masked", 0) or 0)
        n_recv = int(getattr(result, "n_recovered", 0) or 0)
        summary = (
            f"**Frequency editor** ({mode}, method {method}, "
            f"threshold {threshold:g}) — dropped {n_drop}, masked {n_mask}, "
            f"recovered {n_recv}. "
            + ("Edited data is ready — choose **apply / export** below."
               if edited is not None else
               "No edited survey was returned.")
        )

        table = "(no per-row decisions returned)"
        if decisions is not None and getattr(decisions, "empty", True) is False:
            cols = [c for c in ("station", "period", "confidence", "action")
                    if c in decisions.columns]
            table = _df_to_text(decisions, columns=cols or None,
                                max_rows=40, ndigits=4)

        figs = {}
        try:
            ax = plot_frequency_edit_summary(
                sites, edited if edited is not None else sites,
                method=method, ci_hi=ci_hi, ci_lo=ci_lo,
            )
            fig = _as_fig(ax)
            if fig is not None:
                figs["Frequency edit summary"] = fig
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"summary figure skipped: {exc}")
        return summary, table, figs

    def _correction(self, sites, d, warnings):
        """Apply any catalogue correction with full parameter control.

        Drives :class:`CorrectionController` (the same non-destructive chain
        the desktop / web *correction section* uses) so the algorithm is never
        re-implemented here. The corrected ``Sites`` is stashed in
        ``self._corrected`` for the post-processing modal (apply / export),
        exactly like :meth:`_freq_editor`.
        """
        import matplotlib.pyplot as plt
        import pandas as pd

        from ..app.desktop.controllers.correction_controller import (
            _COORD_FN_NAMES,
            _STRAT_FN_NAMES,
            CorrectionController,
        )
        from ..emtools._core import _iter_items
        from ._corrections import (
            CORRECTION_METHODS,
            coerce_kwargs,
            fn_for,
        )

        wf_id = str(d.get("corr_wf") or "").strip()
        fn_name = str(d.get("fn_name") or "").strip()
        if wf_id and wf_id in CORRECTION_METHODS:
            fn_name = fn_name or fn_for(wf_id)
            meta = CORRECTION_METHODS[wf_id]
            kwargs = coerce_kwargs(wf_id, d)
            label = meta.get("label", fn_name)
            category = meta.get("category", "Correction")
        elif fn_name:
            meta, kwargs, label = {}, {}, fn_name
            category = "Correction"
        else:
            raise ValueError(
                "no correction selected (missing 'corr_wf' / 'fn_name')."
            )

        ctrl = CorrectionController()
        ctrl.dark = False
        is_coord = fn_name in _COORD_FN_NAMES
        is_strat = fn_name in _STRAT_FN_NAMES
        # QC is diagnostic (no corrected data to apply); other strat steps and
        # all impedance/rotation/coord steps return corrected Sites.
        diagnostic = fn_name == "_strat_qc"

        if is_strat:
            # Stratagem operates natively on an EDI directory (edi_objects_),
            # not the filtered Sites — load the directory first.
            edi_dir = str(d.get("path") or d.get("data_path") or "").strip()
            if not edi_dir or not os.path.isdir(edi_dir):
                raise RuntimeError(
                    "Stratagem corrections require an EDI directory path "
                    "(load an EDI folder first)."
                )
            ctrl.load_edi_dir(edi_dir)
            corrected = getattr(ctrl, fn_name)(**kwargs)
        else:
            ctrl.set_raw_sites(sites)
            # _call_fn is the controller's universal dispatcher: impedance
            # (emtools), static-shift / rotation / coordinate wrappers — all
            # return corrected Sites directly (coord wrappers write the new
            # coordinates back, unlike apply()'s DataFrame preview path).
            corrected = ctrl._call_fn(fn_name, sites, **kwargs)

        if corrected is None:
            raise RuntimeError(
                f"correction '{label}' produced no result "
                "(check the dataset has valid impedance / coordinates)."
            )
        # Diagnostic steps must not be offered for apply/export.
        self._corrected = None if diagnostic else corrected

        n_sta = len(list(_iter_items(corrected)))
        param_txt = ", ".join(f"{k}={v}" for k, v in kwargs.items()) or "defaults"
        if diagnostic:
            summary = (
                f"**{category} — {label}** over {n_sta} station(s) "
                f"({param_txt}). Diagnostic report below — no data changed."
            )
        else:
            summary = (
                f"**{category} — {label}** applied to {n_sta} station(s) "
                f"({param_txt}). Corrected data is ready — choose "
                "**apply / export** below."
            )

        # ── table ──────────────────────────────────────────────────────────
        if fn_name == "_strat_qc" and ctrl._strat_qc_report is not None \
                and not ctrl._strat_qc_report.empty:
            table = _df_to_text(ctrl._strat_qc_report, max_rows=40, ndigits=3)
        elif kwargs:
            tbl_df = pd.DataFrame(
                [{"parameter": k, "value": v} for k, v in kwargs.items()]
            )
            table = _df_to_text(
                tbl_df, columns=["parameter", "value"], max_rows=30,
            )
        else:
            table = "(no parameters)"

        # ── figure ─────────────────────────────────────────────────────────
        figs = {}
        try:
            if fn_name == "_strat_qc":
                fig = plt.figure(figsize=(7, 5))
                ctrl.plot_strat_qc(fig)
                figs["Stratagem QC report"] = fig
            elif fn_name == "_strat_static_shift":
                fig, ax = plt.subplots(figsize=(7, 4))
                ctrl.plot_strat_ss_factors(ax)
                fig.tight_layout()
                figs["Static-shift factors"] = fig
            elif is_coord:
                fig = _corr_coord_figure(sites, corrected, label)
                if fig is not None:
                    figs["Station positions before / after"] = fig
            else:
                base = ctrl.raw_sites if ctrl.raw_sites is not None else sites
                fig, axes = plt.subplots(1, 2, figsize=(9, 4), sharey=True)
                ctrl.plot_rho_curves(base, axes[0], title="Before")
                ctrl.plot_rho_curves(corrected, axes[1], title="After")
                fig.suptitle(f"{label} — ρₐ before / after", fontsize=10)
                fig.tight_layout()
                figs["Correction before / after"] = fig
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"figure skipped: {exc}")

        return summary, table, figs

    def _layered_model(self, d, warnings):
        """Build and preview a synthetic 1-D layered resistivity model."""
        import matplotlib.pyplot as plt
        import pandas as pd

        from ..forward.synthetic import LayeredModel

        preset = str(d.get("preset", "custom") or "custom").lower()
        try:
            n_layers = int(float(d.get("n_layers", 3) or 3))
        except (TypeError, ValueError):
            n_layers = 3
        n_layers = max(2, min(20, n_layers))

        if preset in ("random", "blocky", "smooth"):
            try:
                depth_max = float(d.get("depth_max", 2000.0) or 2000.0)
            except (TypeError, ValueError):
                depth_max = 2000.0
            if preset == "random":
                model = LayeredModel.random(
                    n_layers, depth_max=depth_max, seed=0)
            elif preset == "blocky":
                model = LayeredModel.blocky(n_layers)
            else:
                model = LayeredModel.smooth(n_layers)
        else:
            rhos = [float(x) for x in _as_list(d.get("resistivities"))]
            thicks = [float(x) for x in _as_list(d.get("thicknesses"))]
            if not rhos:
                rhos, thicks = [100.0, 10.0, 500.0], [300.0, 800.0]
            if len(thicks) != len(rhos) - 1:
                # coerce: keep n-1 thicknesses, pad with the last/ a default
                thicks = thicks[: len(rhos) - 1]
                while len(thicks) < len(rhos) - 1:
                    thicks.append(thicks[-1] if thicks else 500.0)
                warnings.append(
                    f"adjusted to {len(thicks)} thickness(es) for "
                    f"{len(rhos)} layers."
                )
            model = LayeredModel(resistivity=rhos, thickness=thicks)

        fig, ax = plt.subplots(figsize=(4, 5))
        try:
            model.plot(ax=ax)
            fig.tight_layout()
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"plot failed: {exc}")
            plt.close(fig)
            fig = None

        rows = []
        thick = list(model.thickness)
        for i, rho in enumerate(model.resistivity):
            rows.append({
                "layer": "halfspace" if i >= len(thick) else f"L{i + 1}",
                "rho_ohm_m": round(float(rho), 2),
                "thickness_m": (round(float(thick[i]), 1)
                                if i < len(thick) else None),
                "top_depth_m": round(float(model.depth[i]), 1),
            })
        df = pd.DataFrame(rows)
        summary = (
            f"**Layered model** ({preset}) — {model.n_layers} layers, "
            f"ρ {model.resistivity.min():.1f}–{model.resistivity.max():.1f} "
            f"Ω·m, total depth {float(model.depth[-1]):.0f} m."
        )
        table = _df_to_text(
            df, columns=["layer", "rho_ohm_m", "thickness_m", "top_depth_m"],
            max_rows=22,
        )
        return summary, table, ({"1-D resistivity model": fig} if fig else {})
