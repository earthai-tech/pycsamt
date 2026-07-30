# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.agents.plotting
=======================

:class:`PlotAgent` — produce the common visual products the Agent Master
could not make before:

* ``rhophi``         — apparent-resistivity / phase sounding curves
  (:func:`~pycsamt.emtools.inspect.plot_rhoa_phi`),
* ``phase_psection`` — scalar phase pseudo-section
  (:func:`~pycsamt.emtools.inspect.pseudosection` with ``quantity="phi_*"``),
* ``pt_psection``    — phase-tensor (Φ) ellipse pseudo-section
  (:func:`~pycsamt.emtools.tensor.plot_phase_tensor_psection`),
* ``pt_strip``       — single-station phase-tensor ellipse strip vs period
  (:func:`~pycsamt.emtools.tensor.plot_phase_tensor_strip`),
* ``pt_strip_grid``  — phase-tensor ellipse strips tiled by survey line
  (:func:`~pycsamt.emtools.tensor.plot_phase_tensor_strip_grid`).

The agent never calls an LLM; it just turns a validated parameter set into
matplotlib figures, returned under ``AgentResult.data["figures"]`` keyed by a
human title. The Agent Master collects those figures into the chat + the
Figures panel.
"""

from __future__ import annotations

import time
from typing import Any

from ._base import AgentResult

__all__ = ["PlotAgent", "PLOT_KINDS"]

PLOT_KINDS = (
    "rhophi",
    "phase_psection",
    "pt_psection",
    "tipper",
    "pt_map",
    "station_response",
    "strike_profile",
    "pt_strip",
    "pt_strip_grid",
)

# Publication-quality rcParams — applied in a context so global state is
# never mutated.
_PUB_RC = {
    "font.size": 12,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9,
    "axes.linewidth": 1.1,
    "savefig.dpi": 300,
    "figure.dpi": 150,
    "axes.grid": True,
    "grid.alpha": 0.25,
}


def _as_list(val: Any) -> list[str]:
    """Normalise a comma-string / list / None into a clean list of tokens."""
    if val is None:
        return []
    if isinstance(val, (list, tuple)):
        items = list(val)
    else:
        items = str(val).replace(";", ",").split(",")
    out = [str(x).strip() for x in items if str(x).strip()]
    return out


def _norm_components(val: Any) -> tuple[str, ...]:
    """Resolve a component selection to a tuple of ``xy`` / ``yx`` / ``det``."""
    toks = [t.lower() for t in _as_list(val)]
    if not toks or "both" in toks or "all" in toks:
        return ("xy", "yx")
    out: list[str] = []
    for t in toks:
        t = t.replace("rho_", "").replace("phi_", "").replace("z", "")
        if t in ("xy", "yx", "det") and t not in out:
            out.append(t)
    return tuple(out) or ("xy", "yx")


def _period_range(d: dict[str, Any]):
    lo = d.get("period_min")
    hi = d.get("period_max")
    try:
        lo = float(lo) if lo not in (None, "") else None
        hi = float(hi) if hi not in (None, "") else None
    except (TypeError, ValueError):
        return None
    if lo is None and hi is None:
        return None
    # tolerate a single-sided bound
    lo = lo if lo is not None else 1e-9
    hi = hi if hi is not None else 1e9
    return (lo, hi)


def _truthy(v: Any) -> bool:
    if isinstance(v, bool):
        return v
    return str(v).strip().lower() in ("1", "true", "yes", "on", "y")


def _has_tipper(sites: Any) -> bool:
    """True when at least one station carries a vertical-field tipper (T)."""
    from ..emtools._core import _get_t_block, _iter_items

    for ed in _iter_items(sites):
        try:
            _, t, fr = _get_t_block(ed)
        except Exception:  # noqa: BLE001
            continue
        if t is not None and fr is not None:
            return True
    return False


def _norm_parts(val: Any) -> tuple[str, ...]:
    """Resolve tipper real/imag selection."""
    toks = [t.lower() for t in _as_list(val)]
    if not toks or "both" in toks:
        return ("real", "imag")
    out = [t for t in toks if t in ("real", "imag")]
    return tuple(out) or ("real", "imag")


def _filter_sites(sites: Any, station_names: list[str]):
    """Return a Sites with only *station_names* (case-insensitive); the full
    set when the selection is empty or matches nothing."""
    if not station_names:
        return sites
    from ..emtools._core import (
        _iter_items,
        _name,
        _sites_from_items,
        _unwrap,
    )

    wanted = {s.lower() for s in station_names}
    keep = []
    for i, ed in enumerate(_iter_items(sites)):
        if _name(ed, i).lower() in wanted:
            keep.append(_unwrap(ed))
    if not keep:
        return sites
    return _sites_from_items(keep)


class PlotAgent:
    """Render rho/phi curves and phase / phase-tensor pseudo-sections."""

    def __init__(self, **_: Any) -> None:
        # Accept (and ignore) LLM kwargs so the agent registry can build it
        # the same way as every other agent.
        self._last_cost = 0.0

    # ── public API ──────────────────────────────────────────────────────────
    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        t0 = time.time()
        self._last_cost = 0.0

        kind = str(input_data.get("kind", "")).strip()
        if kind not in PLOT_KINDS:
            return AgentResult.failed(
                f"Unknown plot kind {kind!r}. Expected one of {PLOT_KINDS}.",
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
        import matplotlib.pyplot as plt

        from ..emtools._core import ensure_sites

        try:
            sites = ensure_sites(src, recursive=True, strict=False, verbose=0)
        except Exception as exc:  # noqa: BLE001
            return AgentResult.failed(
                f"Could not load EDI data: {exc}",
                elapsed=time.time() - t0,
            )

        # Tipper requires a vertical-field transfer function — many CSAMT /
        # AMT datasets don't record it. Fail clearly (with a reason the caller
        # can turn into "tipper is not available for <line>") instead of
        # drawing an empty axes.
        if kind == "tipper" and not _has_tipper(sites):
            return AgentResult(
                status="failed",
                summary="No tipper (vertical magnetic transfer function) "
                "data is available in this dataset.",
                data={"reason": "no_tipper"},
                elapsed_seconds=time.time() - t0,
            )

        stations = _as_list(input_data.get("stations"))
        pub = _truthy(input_data.get("publication"))
        warnings: list[str] = []

        _METHODS = {
            "rhophi": self._rhophi,
            "phase_psection": self._phase_psection,
            "tipper": self._tipper,
            "pt_map": self._pt_map,
            "pt_psection": self._pt_psection,
            "station_response": self._station_response,
            "strike_profile": self._strike_profile,
            "pt_strip": self._pt_strip,
            "pt_strip_grid": self._pt_strip_grid,
        }
        rc = _PUB_RC if pub else {}
        try:
            with plt.rc_context(rc):
                out = _METHODS[kind](
                    plt, sites, stations, input_data, warnings
                )
        except Exception as exc:  # noqa: BLE001
            return AgentResult.failed(
                f"Plot failed: {exc}",
                hint="Check the station names, components and period range.",
                elapsed=time.time() - t0,
            )

        # A kind returns either (fig, title) or a {title: fig} dict.
        if isinstance(out, dict):
            figures, title = out, kind.replace("_", " ")
        else:
            fig, title = out
            figures = {title: fig}

        if pub:
            for f in figures.values():
                try:
                    f.tight_layout()
                except Exception:  # noqa: BLE001
                    pass

        n_st = stations and len(stations) or "all"
        summary = (
            f"{title} produced ({n_st} station(s)"
            + (", publication style" if pub else "")
            + ")."
        )
        return AgentResult(
            status="success",
            summary=summary,
            data={"figures": figures, "plot_kind": kind},
            warnings=warnings,
            elapsed_seconds=time.time() - t0,
            cost_estimate_usd=0.0,
        )

    # ── kinds ───────────────────────────────────────────────────────────────
    def _rhophi(self, plt, sites, stations, d, warnings):
        from ..emtools.inspect import plot_rhoa_phi

        comps = _norm_components(d.get("components"))
        errorbar = _truthy(d.get("errorbar", True))
        sub = _filter_sites(sites, stations)
        if stations and sub is sites:
            warnings.append(
                "Requested stations not found; plotting all stations."
            )
        figsize = (8.0, 6.5) if _truthy(d.get("publication")) else (7.5, 6.0)
        fig, (ax_r, ax_p) = plt.subplots(2, 1, figsize=figsize, sharex=True)
        plot_rhoa_phi(
            sub,
            components=comps,
            axis="period",
            errorbar=errorbar,
            ax_r=ax_r,
            ax_p=ax_p,
            verbose=0,
        )
        ax_r.set_title(
            "Apparent resistivity & phase  (" + ", ".join(comps) + ")"
        )
        return fig, "Rho/Phi sounding curves"

    def _phase_psection(self, plt, sites, stations, d, warnings):
        from ..emtools.inspect import pseudosection

        comps = _norm_components(d.get("components"))
        pr = _period_range(d)
        sub = _filter_sites(sites, stations)
        n = len(comps)
        figsize = (
            (9.0, 3.4 * n) if _truthy(d.get("publication")) else (8.0, 3.0 * n)
        )
        fig, axes = plt.subplots(n, 1, figsize=figsize, squeeze=False)
        for ax, comp in zip(axes[:, 0], comps):
            pseudosection(
                sub,
                quantity=f"phi_{comp}",
                axis_x="station",
                axis_y="period",
                period_range=pr,
                ax=ax,
                verbose=0,
            )
            ax.set_title(f"Phase pseudo-section  φ_{comp} (deg)")
        return fig, "Phase pseudo-section"

    def _tipper(self, plt, sites, stations, d, warnings):
        sub = _filter_sites(sites, stations)
        if stations and sub is sites:
            warnings.append(
                "Requested stations not found; using all stations."
            )
        view = str(d.get("view", "components") or "components").lower()
        pub = _truthy(d.get("publication"))

        if view in ("arrows", "arrow", "induction", "map", "profile"):
            from ..emtools.tf import plot_induction_arrows

            conv = str(d.get("convention", "park") or "park").lower()
            per = d.get("period")
            if per in (None, ""):
                per = d.get("period_min")
            try:
                periods = [float(per)] if per not in (None, "") else (1.0,)
            except (TypeError, ValueError):
                periods = (1.0,)
            figsize = (8.0, 3.8) if pub else (7.2, 3.4)
            fig, ax = plt.subplots(figsize=figsize)
            plot_induction_arrows(
                sub, periods=periods, convention=conv, ax=ax, verbose=0
            )
            ax.set_title(
                f"Tipper induction arrows  (T≈{periods[0]:g}s, {conv})"
            )
            return fig, "Tipper induction arrows"

        # default: component curves (Tx/Ty, real/imag)
        from ..emtools.inspect import plot_tipper_components

        parts = _norm_parts(d.get("parts"))
        figsize = (8.0, 4.8) if pub else (7.5, 4.5)
        fig, ax = plt.subplots(figsize=figsize)
        plot_tipper_components(
            sub, kind=parts, axis="period", ax=ax, verbose=0
        )
        ax.set_title("Tipper components  (" + ", ".join(parts) + ")")
        return fig, "Tipper components"

    def _station_response(self, plt, sites, stations, d, warnings):
        """Per-station impedance response (Bode ρa/φ consistency diagram),
        one figure per requested component."""
        from ..emtools._core import _iter_items, _name
        from ..emtools.advanced import plot_rho_phase_bode

        comps = [
            c
            for c in _norm_components(d.get("components"))
            if c in ("xy", "yx")
        ] or ["xy"]
        station = stations[0] if stations else None
        if station is None:
            # default to the first station so the inspector always shows one
            for i, ed in enumerate(_iter_items(sites)):
                station = _name(ed, i)
                break
            warnings.append(
                f"No station given; showing the first station ({station})."
            )
        pr = _period_range(d)
        figures: dict = {}
        for comp in comps:
            fig = plot_rho_phase_bode(
                sites,
                station=station,
                component=comp,
                period_range=pr,
                verbose=0,
            )
            figures[f"Response {comp} — {station}"] = fig
        return figures

    def _strike_profile(self, plt, sites, stations, d, warnings):
        from ..emtools.strike import plot_strike_profile

        method = str(d.get("method", "consensus") or "consensus").lower()
        if method not in ("consensus", "sweep", "pt"):
            method = "consensus"
        pr = _period_range(d)
        sort_by = str(d.get("sort_by", "auto") or "auto").lower()
        pub = _truthy(d.get("publication"))
        figsize = (9.2, 4.2) if pub else (8.6, 3.8)
        fig, ax = plt.subplots(figsize=figsize)
        plot_strike_profile(
            _filter_sites(sites, stations),
            method=method,
            band=pr,
            sort_by=sort_by,
            ax=ax,
            verbose=0,
        )
        ax.set_title(f"Strike profile  ({method})")
        return fig, "Strike profile"

    def _pt_map(self, plt, sites, stations, d, warnings):
        from ..emtools.tensor import plot_phase_tensor_map

        per = d.get("period")
        if per in (None, ""):
            per = d.get("period_min")
        try:
            period = float(per) if per not in (None, "") else 10.0
        except (TypeError, ValueError):
            period = 10.0
        c_by = str(d.get("color_by", "") or "").strip().lower() or None
        kw: dict[str, Any] = {}
        if c_by:
            kw["c_by"] = c_by
        scale = d.get("scale")
        if scale not in (None, ""):
            try:
                kw["scale"] = float(scale)
            except (TypeError, ValueError):
                pass
        sub = _filter_sites(sites, stations)
        pub = _truthy(d.get("publication"))
        figsize = (9.0, 7.5) if pub else (8.0, 6.5)
        fig, ax = plt.subplots(figsize=figsize)
        plot_phase_tensor_map(
            sub,
            period=period,
            show_tipper=_has_tipper(sub),
            ax=ax,
            verbose=0,
            **kw,
        )
        ax.set_title(f"Phase-tensor map  (T≈{period:g} s)")
        return fig, "Phase-tensor map"

    def _pt_psection(self, plt, sites, stations, d, warnings):
        from ..emtools.tensor import (
            plot_phase_tensor_psection,
        )

        pr = _period_range(d)
        c_by = str(d.get("color_by", "") or "").strip().lower() or None
        kw: dict[str, Any] = {}
        if c_by:
            kw["c_by"] = c_by
        scale = d.get("scale")
        if scale not in (None, ""):
            try:
                kw["scale"] = float(scale)
            except (TypeError, ValueError):
                pass
        figsize = (11.0, 6.0) if _truthy(d.get("publication")) else (10.0, 5.5)
        fig, ax = plt.subplots(figsize=figsize)
        plot_phase_tensor_psection(
            sites,
            stations=stations or None,
            period_range=pr,
            ax=ax,
            verbose=0,
            **kw,
        )
        return fig, "Phase-tensor (Φ) pseudo-section"

    def _pt_strip(self, plt, sites, stations, d, warnings):
        from ..emtools._core import _iter_items, _name
        from ..emtools.tensor import plot_phase_tensor_strip

        station = stations[0] if stations else None
        if station is None:
            for i, ed in enumerate(_iter_items(sites)):
                station = _name(ed, i)
                break
            warnings.append(
                f"No station given; showing the first station ({station})."
            )
        pr = _period_range(d)
        pub = _truthy(d.get("publication"))
        figsize = (7.5, 1.8) if pub else (7.0, 1.6)
        fig, ax = plt.subplots(figsize=figsize)
        plot_phase_tensor_strip(
            sites,
            station=station,
            period_range=pr,
            ax=ax,
            verbose=0,
        )
        return fig, f"Phase-tensor ellipse strip — {station}"

    def _pt_strip_grid(self, plt, sites, stations, d, warnings):
        from ..emtools._core import _iter_items, _name
        from ..emtools.tensor import (
            plot_phase_tensor_strip_grid,
        )
        from ..site.lines import (
            detect_lines_from_station_ids,
            pick_representative_stations,
        )

        names = [_name(ed, i) for i, ed in enumerate(_iter_items(sites))]
        lines = detect_lines_from_station_ids(names)
        if not lines:
            raise ValueError("No survey lines detected for the strip grid.")

        try:
            per_line = max(1, int(d.get("per_line") or 4))
        except (TypeError, ValueError):
            per_line = 4

        profiles = {
            ln: pick_representative_stations(sts, per_line)
            for ln, sts in lines.items()
        }
        fig = plot_phase_tensor_strip_grid(sites, profiles=profiles, verbose=0)
        return fig, "Phase-tensor ellipse strip grid (by line)"
