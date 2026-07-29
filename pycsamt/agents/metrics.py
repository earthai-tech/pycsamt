# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.agents.metrics
======================

:class:`MetricsAgent` — answer **conversational questions about computed
values** of a survey line (or all lines), inline, without a modal or a figure.

Where :class:`~pycsamt.agents.tooling.ToolAgent` opens a parameter modal and
returns a table + figure, the metrics agent turns a question like *"what's the
strike of L22PLT?"* or *"azimuth of all lines?"* into a short, friendly text
answer with the actual number(s).

Supported value *kinds* (see :data:`METRIC_KINDS`):

* ``strike``        — regional geoelectric strike (consensus, ±90°),
* ``azimuth``       — profile bearing from station geometry (PCA, 0–180°N),
* ``dimensionality``— dominant 1-D / 2-D / 3-D class and the distribution,
* ``skew``          — mean phase-tensor skew |β|,
* ``stations``      — station count (and names),
* ``periods``       — period range (s),
* ``frequencies``   — frequency range (Hz),
* ``coordinates``   — lat/lon bounding box and approximate profile length,
* ``quality``       — mean QC score, flagged stations, tipper presence,
* ``summary``       — a one-paragraph digest combining the above.

The agent never calls an LLM. :func:`parse_metric_request` and
:func:`looks_like_metric_query` are shared, testable helpers the chat router
uses to detect a value question and pull out which kinds / scope were asked.
"""

from __future__ import annotations

import time
from typing import Any

from ._base import AgentResult
from .plotting import (
    _filter_sites,
    _has_tipper,
)
from .tooling import (
    _circular_strike_mean,
    _ll_to_utm,
    _station_coords,
)

__all__ = [
    "MetricsAgent",
    "METRIC_KINDS",
    "parse_metric_request",
    "looks_like_metric_query",
]

METRIC_KINDS = (
    "strike",
    "azimuth",
    "dimensionality",
    "skew",
    "stations",
    "periods",
    "frequencies",
    "coordinates",
    "quality",
    "summary",
)

# Ordered synonym table: most-specific phrases first. First match per kind.
_SYNONYMS: dict[str, tuple[str, ...]] = {
    "summary": (
        "summary",
        "overview",
        "tell me about",
        "describe the line",
        "describe this line",
        "about this line",
        "about the line",
        "give me a rundown",
        "brief on",
        "characterise",
        "characterize",
    ),
    "strike": (
        "geoelectric strike",
        "strike direction",
        "strike angle",
        "regional strike",
        "strike of",
        "strike for",
        "strike",
    ),
    "azimuth": (
        "azimuth",
        "bearing",
        "profile direction",
        "line direction",
        "profile orientation",
        "line orientation",
        "orientation",
    ),
    "dimensionality": (
        "dimensionality",
        "dimension",
        "1d 2d 3d",
        "1d/2d/3d",
        "how many dimensions",
    ),
    "skew": ("skew", "phase tensor skew", "beta", "asymmetry"),
    "stations": (
        "how many stations",
        "number of stations",
        "station count",
        "count of stations",
        "n stations",
        "stations",
    ),
    "periods": (
        "period range",
        "period band",
        "periods",
        "shortest period",
        "longest period",
    ),
    "frequencies": (
        "frequency range",
        "frequency band",
        "frequencies",
        "frequency",
    ),
    "coordinates": (
        "coordinates",
        "bounding box",
        "profile length",
        "line length",
        "how long is the line",
        "extent",
        "location",
    ),
    "quality": (
        "data quality",
        "quality score",
        "qc score",
        "quality",
    ),
}

# Words that mean "the full survey / every line".
_ALL_SCOPE = (
    "all lines",
    "all the lines",
    "every line",
    "each line",
    "all profiles",
    "across lines",
    "all surveys",
    "per line",
)

# Phrasing that marks a *visual* / workflow request, not a value question.
_VISUAL_OR_WORKFLOW = (
    "plot",
    "rose",
    "map",
    "pseudosection",
    "pseudo-section",
    "section",
    "figure",
    "chart",
    "diagram",
    "analyzer",
    "analyser",
    "analysis",
    "run ",
    "compute",
    "estimate",
    "inspect",
    "screen",
    "flag",
    "render",
    "draw",
    "export",
    "save",
)

# Kinds whose names double as pyCSAMT *concepts* ("what is geoelectric
# strike?" is a QUESTION, "what's the strike of L22PLT?" is a value). These
# require explicit data scope to count as a metric question; the survey-only
# kinds (stations / periods / …) accept plain question phrasing.
_CONCEPT_OVERLAP = frozenset(
    {
        "strike",
        "azimuth",
        "dimensionality",
        "skew",
        "summary",
    }
)


def parse_metric_request(text: str) -> tuple[list[str], bool]:
    """Return ``(kinds, all_lines)`` parsed from *text*.

    ``kinds`` is the ordered, de-duplicated list of metric kinds the user
    asked about; ``all_lines`` is ``True`` when the request spans every line.

    Examples
    --------
    >>> parse_metric_request("what is the strike of L22PLT?")
    (['strike'], False)
    >>> parse_metric_request("give me the azimuth for all lines")
    (['azimuth'], True)
    >>> parse_metric_request("tell me about this line")
    (['summary'], False)
    """
    t = (text or "").lower()
    all_lines = any(p in t for p in _ALL_SCOPE)
    kinds: list[str] = []
    for kind, phrases in _SYNONYMS.items():
        if any(p in t for p in phrases) and kind not in kinds:
            kinds.append(kind)
    # "summary" is greedy — if it co-occurs with specific metrics, prefer the
    # specific ones (the user asked for particular values, not a digest).
    if "summary" in kinds and len(kinds) > 1:
        kinds.remove("summary")
    return kinds, all_lines


def looks_like_metric_query(text: str) -> bool:
    """Heuristic: does *text* ask for a computed value of a line?

    True when the message names a value (strike, azimuth, …) **and** is phrased
    as a question / scoped to a line, and is *not* a plot/workflow request.

    Examples
    --------
    >>> looks_like_metric_query("what's the strike of L22PLT?")
    True
    >>> looks_like_metric_query("plot the strike rose")
    False
    >>> looks_like_metric_query("run phase tensor analysis")
    False
    """
    t = (text or "").lower().strip()
    if not t:
        return False
    if any(w in t for w in _VISUAL_OR_WORKFLOW):
        return False
    kinds, all_lines = parse_metric_request(t)
    if not kinds:
        return False
    # A value word alone isn't enough ("strike" could be a workflow). Require
    # question-style or line-scoped phrasing.
    question_like = t.endswith("?") or t.startswith(
        (
            "what",
            "how",
            "which",
            "give",
            "show",
            "tell",
            "list",
            "whats",
            "what's",
        )
    )
    scoped = (
        all_lines
        or " of " in t
        or " for " in t
        or "this line" in t
        or "these lines" in t
        or "the line" in t
        or "the lines" in t
        or "the survey" in t
        or "this survey" in t
        or "my data" in t
        or "my survey" in t
        or "in the line" in t
    )
    # Concept-overlapping value words ("strike", "summary", …) need explicit
    # data scope so "what is geoelectric strike?" / "tell me about the Sites
    # data model" stay package questions, not value lookups.
    if set(kinds) <= _CONCEPT_OVERLAP:
        return scoped
    return bool(question_like or scoped)


def _fmt_deg(v: float) -> str:
    return f"{v:.0f}°"


class MetricsAgent:
    """Compute and phrase per-line survey values as inline chat answers."""

    def __init__(self, **_: Any) -> None:
        self._last_cost = 0.0

    # ── public API ──────────────────────────────────────────────────────────
    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        t0 = time.time()
        self._last_cost = 0.0

        kinds = input_data.get("kinds") or [input_data.get("kind", "summary")]
        kinds = [str(k).strip() for k in kinds if str(k).strip() in METRIC_KINDS]
        if not kinds:
            kinds = ["summary"]

        src = (
            input_data.get("sites")
            or input_data.get("path")
            or input_data.get("data_path")
        )
        if src is None:
            return AgentResult.failed(
                "No data available. Load an EDI dataset or name a survey line.",
                elapsed=time.time() - t0,
            )

        from ..emtools._core import ensure_sites

        try:
            sites = ensure_sites(src, recursive=True, strict=False, verbose=0)
        except Exception as exc:  # noqa: BLE001
            return AgentResult.failed(
                f"Could not load survey data: {exc}",
                elapsed=time.time() - t0,
            )

        stations = input_data.get("stations")
        if stations:
            from .plotting import _as_list

            sites = _filter_sites(sites, _as_list(stations))

        label = str(input_data.get("label") or "").strip()
        warnings: list[str] = []
        values: dict[str, str] = {}
        for kind in kinds:
            try:
                values[kind] = self._compute(kind, sites, warnings)
            except Exception as exc:  # noqa: BLE001
                values[kind] = "could not be computed"
                warnings.append(f"{kind}: {exc}")

        answer = self._phrase(label, kinds, values)
        return AgentResult(
            status="success",
            summary=answer,
            data={"values": values, "kinds": kinds, "tool_kind": "metrics"},
            warnings=warnings,
            elapsed_seconds=time.time() - t0,
            cost_estimate_usd=0.0,
        )

    # ── phrasing ────────────────────────────────────────────────────────────
    @staticmethod
    def _phrase(label: str, kinds: list[str], values: dict[str, str]) -> str:
        where = f" for **{label}**" if label else ""
        if kinds == ["summary"]:
            return f"Here's the rundown{where}: {values['summary']}"
        if len(kinds) == 1:
            k = kinds[0]
            return f"The {k}{where} — {values[k]}."
        lines = "\n".join(f"- **{k}**: {values[k]}" for k in kinds)
        return f"Here's what I found{where}:\n{lines}"

    # ── dispatch ────────────────────────────────────────────────────────────
    def _compute(self, kind: str, sites, warnings) -> str:
        return getattr(self, f"_m_{kind}")(sites, warnings)

    # ── individual metrics ──────────────────────────────────────────────────
    def _m_strike(self, sites, warnings) -> str:
        from ..emtools.strike import estimate_strike_consensus

        band = None
        res = estimate_strike_consensus(sites, band=band, verbose=0)
        df = res.frame if hasattr(res, "frame") else res
        reg = _circular_strike_mean(df["ang"]) if "ang" in df else float("nan")
        import numpy as np

        if not np.isfinite(reg):
            return "no usable strike could be estimated"
        compass = f"N{reg:.0f}°E" if reg >= 0 else f"N{abs(reg):.0f}°W"
        return (
            f"regional geoelectric strike ≈ {compass}, with the inherent "
            f"90° ambiguity, from {len(df)} station(s)"
        )

    def _m_azimuth(self, sites, warnings) -> str:
        import numpy as np

        from ..gis.coord_correction import _pca_azimuth

        es, ns = [], []
        for _, lat, lon in _station_coords(sites):
            if lat is None or lon is None:
                continue
            try:
                e, n, _ = _ll_to_utm(lat, lon, None, "N", "WGS84")
                es.append(e)
                ns.append(n)
            except Exception:  # noqa: BLE001
                continue
        if len(es) < 2:
            return "not enough located stations to estimate a profile azimuth"
        az = _pca_azimuth(np.asarray(es), np.asarray(ns))
        return (
            f"profile azimuth ≈ {_fmt_deg(az)} from north (0–180°), "
            f"from {len(es)} located station(s)"
        )

    def _m_dimensionality(self, sites, warnings) -> str:
        from ..emtools.dimensionality import (
            classify_dimensionality,
        )

        res = classify_dimensionality(sites, verbose=0)
        df = res.frame if hasattr(res, "frame") else res
        if "dim" not in df or df.empty:
            return "could not be classified"
        lbl = {0: "1-D", 1: "2-D", 2: "3-D", 3: "3-D"}
        counts = df["dim"].value_counts().to_dict()
        total = int(sum(counts.values())) or 1
        dist = ", ".join(
            f"{lbl.get(k, k)} {100 * v / total:.0f}%" for k, v in sorted(counts.items())
        )
        dom = max(counts, key=counts.get)
        return (
            f"dominantly {lbl.get(dom, dom)} "
            f"({100 * counts[dom] / total:.0f}% of cells); distribution {dist}"
        )

    def _m_skew(self, sites, warnings) -> str:
        import numpy as np

        from ..emtools.dimensionality import (
            classify_dimensionality,
        )

        res = classify_dimensionality(sites, verbose=0)
        df = res.frame if hasattr(res, "frame") else res
        if "beta_abs" not in df or df.empty:
            return "skew is unavailable"
        v = float(np.nanmean(df["beta_abs"]))
        tag = "→ 3-D / distorted" if v > 3.0 else "→ 1-D / 2-D regime"
        return f"mean phase-tensor skew |β| ≈ {v:.1f}° ({tag})"

    def _m_stations(self, sites, warnings) -> str:
        rows = self._scan(sites)
        names = [r.get("station", "?") for r in rows]
        head = ", ".join(str(n) for n in names[:6])
        more = f" … (+{len(names) - 6} more)" if len(names) > 6 else ""
        return (
            f"{len(names)} station(s): {head}{more}" if names else "no stations found"
        )

    def _m_periods(self, sites, warnings) -> str:
        rows = self._scan(sites)
        tmins = [r["t_min_s"] for r in rows if r.get("t_min_s")]
        tmaxs = [r["t_max_s"] for r in rows if r.get("t_max_s")]
        nfr = [r["n_freq"] for r in rows if r.get("n_freq")]
        if not tmins or not tmaxs:
            return "period information is unavailable"
        return (
            f"period range {min(tmins):.3g}–{max(tmaxs):.3g} s "
            f"(up to {max(nfr) if nfr else 0} periods per station)"
        )

    def _m_frequencies(self, sites, warnings) -> str:
        rows = self._scan(sites)
        tmins = [r["t_min_s"] for r in rows if r.get("t_min_s")]
        tmaxs = [r["t_max_s"] for r in rows if r.get("t_max_s")]
        if not tmins or not tmaxs:
            return "frequency information is unavailable"
        f_hi = 1.0 / min(tmins)
        f_lo = 1.0 / max(tmaxs)
        return f"frequency range {f_lo:.3g}–{f_hi:.3g} Hz"

    def _m_coordinates(self, sites, warnings) -> str:
        import numpy as np

        coords = [
            (la, lo)
            for _, la, lo in _station_coords(sites)
            if la is not None and lo is not None
        ]
        if not coords:
            return "no station coordinates are available"
        lats = np.array([c[0] for c in coords])
        lons = np.array([c[1] for c in coords])
        # approximate length: bounding-box diagonal in UTM
        length_km = float("nan")
        try:
            es, ns = [], []
            for la, lo in coords:
                e, n, _ = _ll_to_utm(la, lo, None, "N", "WGS84")
                es.append(e)
                ns.append(n)
            es, ns = np.asarray(es), np.asarray(ns)
            length_km = (
                float(np.hypot(es.max() - es.min(), ns.max() - ns.min())) / 1000.0
            )
        except Exception:  # noqa: BLE001
            pass
        span = f", profile span ≈ {length_km:.1f} km" if np.isfinite(length_km) else ""
        return (
            f"spans lat {lats.min():.4f}–{lats.max():.4f}°, "
            f"lon {lons.min():.4f}–{lons.max():.4f}°{span} "
            f"({len(coords)} located station(s))"
        )

    def _m_quality(self, sites, warnings) -> str:
        import numpy as np

        rows = self._scan(sites)
        if not rows:
            return "no stations with valid data found"
        scores = np.array([r.get("qc_score") or 0 for r in rows], float)
        n_flag = sum(
            1
            for r in rows
            if not r.get("has_z")
            or not r.get("has_coords")
            or (r.get("qc_score") or 0) < 50
        )
        tip = "present" if _has_tipper(sites) else "absent"
        return (
            f"mean QC score {np.nanmean(scores):.0f}/100, "
            f"{n_flag} of {len(rows)} station(s) flagged; tipper {tip}"
        )

    def _m_summary(self, sites, warnings) -> str:
        parts = []
        for k in (
            "stations",
            "periods",
            "strike",
            "dimensionality",
            "quality",
        ):
            try:
                parts.append(self._compute(k, sites, warnings))
            except Exception as exc:  # noqa: BLE001
                warnings.append(f"{k}: {exc}")
        return "; ".join(parts) + "."

    # ── shared ──────────────────────────────────────────────────────────────
    @staticmethod
    def _scan(sites) -> list:
        from .loader import _quality_scan

        rows, _ = _quality_scan(sites)
        return rows or []
