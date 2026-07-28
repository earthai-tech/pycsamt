# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.site.report
===================

Human-friendly display layer for :class:`~pycsamt.site.base.Site` and
:class:`~pycsamt.site.base.Sites` objects.

Two classes are provided:

* :class:`SiteReport`  — statistics and rich display for a **single** site.
* :class:`SitesReport` — statistics and rich display for a **collection**.

Both compute their statistics lazily at construction time so that
repeated calls to :meth:`report` are cheap.

Quick start
-----------
::

    from pycsamt.site.report import SiteReport, SitesReport

    # --- single site ---
    report = SiteReport(site)
    report.report()  # prints to terminal
    d = report.to_dict()  # plain dict
    df = report.to_dataframe()  # pandas DataFrame (resphase)

    # --- collection ---
    rep = SitesReport(sites)
    rep.report()  # full survey summary
    rep.report(top=10)  # first 10 stations only
    df = rep.to_dataframe()  # one row per station
    d = rep.to_dict()  # list of per-station dicts

Notes
-----
The :meth:`report` method uses ``rich`` when available and degrades to
plain-text output when ``rich`` is not installed.

All statistics are computed from the arrays exposed by
:class:`~pycsamt.site.base.SiteMixin`: ``freq``, ``z``, ``rho``,
``phase``, ``tipper``.  Missing arrays are handled gracefully and
reported as ``"—"`` in the output.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from ..api.view import maybe_wrap_frame

# ---------------------------------------------------------------------------
# Optional rich import
# ---------------------------------------------------------------------------
try:
    from rich.console import Console
    from rich.panel import Panel
    from rich.table import Table
    from rich.text import Text

    _RICH = True
except ImportError:
    _RICH = False

__all__ = ["SiteReport", "SitesReport"]

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_COMPONENTS = ("Zxx", "Zxy", "Zyx", "Zyy")
_COMP_LOWER = tuple(c.lower() for c in _COMPONENTS)
_BAR_CHARS = "█▓▒░"  # filled → empty
_BAR_WIDTH = 10


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _bar(fraction: float, width: int = _BAR_WIDTH) -> str:
    """Return a Unicode block-character progress bar string.

    Examples
    --------
    >>> _bar(1.0)  # '██████████'
    >>> _bar(0.75)  # '███████░░░'
    >>> _bar(0.0)  # '░░░░░░░░░░'
    """
    fraction = max(0.0, min(1.0, float(fraction)))
    filled = round(fraction * width)
    empty = width - filled
    return "█" * filled + "░" * empty


def _fmt_rho(mean: float | None, std: float | None) -> str:
    if mean is None:
        return "—"
    s = f"{mean:.0f}"
    if std is not None:
        s += f" ± {std:.0f}"
    return s + " Ω·m"


def _fmt_phi(mean: float | None, std: float | None) -> str:
    if mean is None:
        return "—"
    s = f"{mean:.1f}"
    if std is not None:
        s += f" ± {std:.1f}"
    return s + "°"


def _fmt_freq(f: float | None) -> str:
    if f is None:
        return "—"
    if f >= 1e3:
        return f"{f / 1e3:.3g} kHz"
    if f >= 1.0:
        return f"{f:.4g} Hz"
    if f >= 1e-3:
        return f"{f * 1e3:.4g} mHz"
    return f"{f:.3e} Hz"


def _check(has: bool) -> str:
    return "✓" if has else "✗"


def _safe_arr(arr: Any, *, real_only: bool = True) -> np.ndarray | None:
    """Return a non-empty array, or None.

    When *real_only* is True and the array is complex, only the real
    part is kept (avoids ComplexWarning when computing float statistics).
    """
    if arr is None:
        return None
    a = np.asarray(arr)
    if a.size == 0:
        return None
    if real_only and np.iscomplexobj(a):
        a = a.real
    return a.astype(float, copy=False)


def _rho_stats(
    rho: np.ndarray | None,
    comp_idx: int,
) -> tuple[float | None, float | None]:
    """Return (mean, std) of log10-ρ_a for component *comp_idx*."""
    if rho is None:
        return None, None
    try:
        r = np.asarray(rho, dtype=float)
        if r.ndim == 3:  # (n, 2, 2)
            row, col = divmod(comp_idx, 2)
            col_data = r[:, row, col]
        elif r.ndim == 2:  # (n, 4)  Zxx=0 Zxy=1 Zyx=2 Zyy=3
            col_data = r[:, comp_idx]
        else:
            return None, None
        col_data = col_data[np.isfinite(col_data) & (col_data > 0)]
        if col_data.size == 0:
            return None, None
        return float(np.mean(col_data)), float(np.std(col_data))
    except Exception:  # noqa: BLE001
        return None, None


def _phase_stats(
    phase: np.ndarray | None,
    comp_idx: int,
) -> tuple[float | None, float | None]:
    """Return (mean, std) of phase (degrees) for component *comp_idx*."""
    if phase is None:
        return None, None
    try:
        p = np.asarray(phase, dtype=float)
        if p.ndim == 3:
            row, col = divmod(comp_idx, 2)
            col_data = p[:, row, col]
        elif p.ndim == 2:
            col_data = p[:, comp_idx]
        else:
            return None, None
        col_data = col_data[np.isfinite(col_data)]
        if col_data.size == 0:
            return None, None
        return float(np.mean(col_data)), float(np.std(col_data))
    except Exception:  # noqa: BLE001
        return None, None


def _quality_pct(arr: np.ndarray | None, comp_idx: int) -> float | None:
    """Fraction of finite values for one impedance component.

    Handles both real and complex Z arrays without raising ComplexWarning.
    """
    if arr is None:
        return None
    try:
        a = np.asarray(arr)
        if a.ndim == 3:
            row, col = divmod(comp_idx, 2)
            col_data = a[:, row, col]
        elif a.ndim == 2:
            col_data = a[:, comp_idx]
        else:
            return None
        if col_data.size == 0:
            return None
        if np.iscomplexobj(col_data):
            finite = np.isfinite(col_data.real) & np.isfinite(col_data.imag)
        else:
            finite = np.isfinite(col_data.astype(float, copy=False))
        return float(finite.sum() / col_data.size)
    except Exception:  # noqa: BLE001
        return None


def _has_component(site: Any, name: str) -> bool:
    """Return True if *site* has a non-empty component *name*."""
    try:
        return bool(site.has_component(name))
    except AttributeError:
        pass
    try:
        z = site.z
        if z is not None and np.asarray(z).size > 0:
            return True
    except Exception:  # noqa: BLE001
        pass
    return False


def _console() -> Any:
    if _RICH:
        return Console()
    return None


def _print_plain(lines: list[str]) -> None:
    for ln in lines:
        print(ln)  # noqa: T201


# ---------------------------------------------------------------------------
# SiteReport
# ---------------------------------------------------------------------------


class SiteReport:
    """Statistics and display for a single :class:`~pycsamt.site.base.Site`.

    Parameters
    ----------
    site : Site-like
        Any object that exposes ``name``, ``coords``, ``freq``, ``z``,
        ``rho``, ``phase``, and ``tipper`` as per
        :class:`~pycsamt.site.base.SiteMixin`.

    Examples
    --------
    ::

        from pycsamt.site.report import SiteReport

        rep = SiteReport(site)
        rep.report()  # rich terminal output
        d = rep.to_dict()  # machine-readable dict
    """

    def __init__(self, site: Any) -> None:
        self._site = site
        self._stats = self._compute()

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def report(self, *, detail: bool = False) -> None:
        """Print a rich panel with site statistics.

        Parameters
        ----------
        detail : bool
            If ``True``, include per-frequency Z and ρ–φ tables.
        """
        if _RICH:
            self._rich_report(detail=detail)
        else:
            _print_plain(self._plain_lines(detail=detail))

    def summary(self) -> str:
        """Return a one-line summary string."""
        s = self._stats
        return (
            f"SiteReport({s['name']!r}  "
            f"{s['nfreq']} freq  "
            f"ρ_xy={_fmt_rho(s['rho_xy_mean'], None)}  "
            f"φ_xy={_fmt_phi(s['phi_xy_mean'], None)})"
        )

    def to_dict(self) -> dict[str, Any]:
        """Return a plain dict of all computed statistics."""
        return dict(self._stats)

    def to_dataframe(self, kind: str = "resphase", *, api: bool | None = None) -> Any:
        """Export site arrays to a :class:`pandas.DataFrame`.

        Parameters
        ----------
        kind : str
            Passed to :meth:`~pycsamt.site.base.SiteMixin.to_dataframe`.
        api : bool or None, default None
            ``True`` forces API view wrapping; ``False`` returns a bare pandas
            dataframe; ``None`` (default) defers to the global
            :data:`~pycsamt.api.view.config.PYCSAMT_API_VIEW` setting.
        """
        return self._site.to_dataframe(kind, api=api)

    def __repr__(self) -> str:
        return self.summary()

    # ------------------------------------------------------------------
    # Statistics
    # ------------------------------------------------------------------

    def _compute(self) -> dict[str, Any]:
        site = self._site
        try:
            lat, lon, elev = site.coords
        except Exception:  # noqa: BLE001
            lat, lon, elev = None, None, None

        freq = _safe_arr(getattr(site, "freq", None))
        nfreq = int(freq.size) if freq is not None else 0
        freq_min = float(freq.min()) if freq is not None else None
        freq_max = float(freq.max()) if freq is not None else None
        per_min = 1.0 / freq_max if freq_max else None
        per_max = 1.0 / freq_min if freq_min else None

        rho = getattr(site, "rho", None)
        phase = getattr(site, "phase", None)
        z = getattr(site, "z", None)
        tip = getattr(site, "tipper", None)

        # Z component presence
        comp_present = {}
        for i, c in enumerate(_COMPONENTS):
            try:
                comp_present[c] = _has_component(site, c.lower())
            except Exception:  # noqa: BLE001
                comp_present[c] = False
        has_tip = tip is not None and np.asarray(tip).size > 0

        # Rho/phase stats for Zxy (idx=1) and Zyx (idx=2)
        rho_xy_m, rho_xy_s = _rho_stats(rho, 1)
        rho_yx_m, rho_yx_s = _rho_stats(rho, 2)
        phi_xy_m, phi_xy_s = _phase_stats(phase, 1)
        phi_yx_m, phi_yx_s = _phase_stats(phase, 2)

        # Data quality per component
        quality = {}
        for i, c in enumerate(_COMPONENTS):
            q = _quality_pct(z, i)
            quality[c] = q

        return {
            "name": getattr(site, "name", "?"),
            "lat": lat,
            "lon": lon,
            "elev": elev,
            "nfreq": nfreq,
            "freq_min": freq_min,
            "freq_max": freq_max,
            "per_min": per_min,
            "per_max": per_max,
            "components": comp_present,
            "has_tipper": has_tip,
            "rho_xy_mean": rho_xy_m,
            "rho_xy_std": rho_xy_s,
            "rho_yx_mean": rho_yx_m,
            "rho_yx_std": rho_yx_s,
            "phi_xy_mean": phi_xy_m,
            "phi_xy_std": phi_xy_s,
            "phi_yx_mean": phi_yx_m,
            "phi_yx_std": phi_yx_s,
            "quality": quality,
        }

    # ------------------------------------------------------------------
    # Rich output
    # ------------------------------------------------------------------

    def _rich_report(self, *, detail: bool = False) -> None:
        s = self._stats
        con = Console()
        t = Table.grid(padding=(0, 2))
        t.add_column(style="bold dim", no_wrap=True)
        t.add_column(style="white")

        def row(k: str, v: str) -> None:
            t.add_row(k, v)

        coords = (
            f"{s['lat']:.5f}°N  {s['lon']:.5f}°E  elev {s['elev']:.0f} m"
            if s["lat"] is not None
            else "—"
        )
        freq_str = (
            f"{s['nfreq']}  ·  "
            f"{_fmt_freq(s['freq_min'])} → {_fmt_freq(s['freq_max'])}"
            if s["nfreq"] > 0
            else "—"
        )
        per_str = (
            f"{_fmt_freq(s['per_min'])} → {_fmt_freq(s['per_max'])}"
            if s["per_min"] is not None
            else "—"
        )
        comp_str = "  ".join(
            f"[green]{c} ✓[/green]" if v else f"[red]{c} ✗[/red]"
            for c, v in s["components"].items()
        )
        tip_str = "[green]✓  Tx  Ty[/green]" if s["has_tipper"] else "[dim]—[/dim]"

        row("Coordinates", coords)
        row("Frequencies", freq_str)
        row("Periods", per_str)
        row("Components", comp_str)
        row("Tipper", tip_str)
        t.add_row("", "")

        row("ρ_xy", _fmt_rho(s["rho_xy_mean"], s["rho_xy_std"]))
        row("ρ_yx", _fmt_rho(s["rho_yx_mean"], s["rho_yx_std"]))
        row("φ_xy", _fmt_phi(s["phi_xy_mean"], s["phi_xy_std"]))
        row("φ_yx", _fmt_phi(s["phi_yx_mean"], s["phi_yx_std"]))
        t.add_row("", "")

        for comp, q in s["quality"].items():
            if q is not None:
                bar = _bar(q)
                pct = f"{q * 100:.0f}%"
                row(f"Quality {comp}", f"{bar} {pct}")

        panel = Panel(
            t,
            title=f"[bold cyan]Site: {s['name']}[/bold cyan]",
            border_style="bright_blue",
            expand=False,
        )
        con.print()
        con.print(panel)
        con.print()

    # ------------------------------------------------------------------
    # Plain-text output
    # ------------------------------------------------------------------

    def _plain_lines(self, *, detail: bool = False) -> list[str]:
        s = self._stats
        lines = [
            f"\nSite: {s['name']}",
            "─" * 50,
        ]

        def row(k: str, v: str) -> None:
            lines.append(f"  {k:<18} {v}")

        coords = (
            f"{s['lat']:.5f}°N  {s['lon']:.5f}°E  elev {s['elev']:.0f} m"
            if s["lat"] is not None
            else "—"
        )
        row("Coordinates", coords)
        if s["nfreq"] > 0:
            row(
                "Frequencies",
                f"{s['nfreq']}  ·  "
                f"{_fmt_freq(s['freq_min'])} → {_fmt_freq(s['freq_max'])}",
            )
        comp_str = "  ".join(
            f"{c}{'✓' if v else '✗'}" for c, v in s["components"].items()
        )
        row("Components", comp_str)
        row("Tipper", "✓" if s["has_tipper"] else "—")
        lines.append("")
        row("ρ_xy", _fmt_rho(s["rho_xy_mean"], s["rho_xy_std"]))
        row("ρ_yx", _fmt_rho(s["rho_yx_mean"], s["rho_yx_std"]))
        row("φ_xy", _fmt_phi(s["phi_xy_mean"], s["phi_xy_std"]))
        row("φ_yx", _fmt_phi(s["phi_yx_mean"], s["phi_yx_std"]))
        lines.append("")
        for comp, q in s["quality"].items():
            if q is not None:
                row(f"Quality {comp}", f"{_bar(q)} {q * 100:.0f}%")
        lines.append("")
        return lines


# ---------------------------------------------------------------------------
# SitesReport
# ---------------------------------------------------------------------------


class SitesReport:
    """Statistics and display for a :class:`~pycsamt.site.base.Sites` collection.

    Parameters
    ----------
    sites : Sites-like
        Any iterable of :class:`~pycsamt.site.base.Site`-like objects.

    Examples
    --------
    ::

        from pycsamt.site.report import SitesReport

        rep = SitesReport(sites)
        rep.report()  # full survey panel + per-station table
        rep.report(top=10)  # first 10 stations only
        df = rep.to_dataframe()  # one row per station
    """

    def __init__(self, sites: Any) -> None:
        self._sites = list(sites)
        self._records: list[dict[str, Any]] = [
            SiteReport(s)._compute() for s in self._sites
        ]
        self._survey = self._survey_stats()

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def report(
        self,
        *,
        top: int | None = None,
        detail: bool = False,
    ) -> None:
        """Print a full survey report.

        Parameters
        ----------
        top : int, optional
            Limit the per-station table to the first *top* stations.
        detail : bool
            If ``True``, print additional per-station statistics.
        """
        records = self._records[:top] if top else self._records
        if _RICH:
            self._rich_report(records, detail=detail)
        else:
            _print_plain(self._plain_lines(records, detail=detail))

    def summary(self) -> str:
        sv = self._survey
        return (
            f"SitesReport({sv['n_stations']} stations  "
            f"freq {_fmt_freq(sv['freq_min_common'])} → "
            f"{_fmt_freq(sv['freq_max_common'])})"
        )

    def to_dict(self) -> list[dict[str, Any]]:
        """Return a list of per-station stat dicts."""
        return [dict(r) for r in self._records]

    def to_dataframe(self, *, api: bool | None = None) -> Any:
        """Return a :class:`pandas.DataFrame` with one row per station."""
        import pandas as pd  # noqa: PLC0415

        rows = []
        for r in self._records:
            rows.append(
                {
                    "station": r["name"],
                    "lat": r["lat"],
                    "lon": r["lon"],
                    "elev": r["elev"],
                    "nfreq": r["nfreq"],
                    "freq_min": r["freq_min"],
                    "freq_max": r["freq_max"],
                    "has_Zxx": r["components"].get("Zxx", False),
                    "has_Zxy": r["components"].get("Zxy", False),
                    "has_Zyx": r["components"].get("Zyx", False),
                    "has_Zyy": r["components"].get("Zyy", False),
                    "has_tipper": r["has_tipper"],
                    "rho_xy": r["rho_xy_mean"],
                    "rho_xy_std": r["rho_xy_std"],
                    "rho_yx": r["rho_yx_mean"],
                    "rho_yx_std": r["rho_yx_std"],
                    "phi_xy": r["phi_xy_mean"],
                    "phi_xy_std": r["phi_xy_std"],
                    "phi_yx": r["phi_yx_mean"],
                    "phi_yx_std": r["phi_yx_std"],
                }
            )
        df = pd.DataFrame(rows)

        return maybe_wrap_frame(
            df,
            api=api,
            name="sites_report",
            kind="site.report",
            source=self._sites,
            description="Per-station site report statistics.",
        )

    def __repr__(self) -> str:
        return self.summary()

    # ------------------------------------------------------------------
    # Survey-level statistics
    # ------------------------------------------------------------------

    def _survey_stats(self) -> dict[str, Any]:
        recs = self._records
        if not recs:
            return {"n_stations": 0}

        lats = [r["lat"] for r in recs if r["lat"] is not None]
        lons = [r["lon"] for r in recs if r["lon"] is not None]
        elevs = [r["elev"] for r in recs if r["elev"] is not None]

        fmins = [r["freq_min"] for r in recs if r["freq_min"] is not None]
        fmaxs = [r["freq_max"] for r in recs if r["freq_max"] is not None]

        # Component availability counts
        comp_counts: dict[str, int] = {c: 0 for c in _COMPONENTS}
        tip_count = 0
        nfreq_vals = []
        for r in recs:
            for c in _COMPONENTS:
                if r["components"].get(c, False):
                    comp_counts[c] += 1
            if r["has_tipper"]:
                tip_count += 1
            if r["nfreq"] > 0:
                nfreq_vals.append(r["nfreq"])

        return {
            "n_stations": len(recs),
            "lat_min": min(lats) if lats else None,
            "lat_max": max(lats) if lats else None,
            "lon_min": min(lons) if lons else None,
            "lon_max": max(lons) if lons else None,
            "elev_min": min(elevs) if elevs else None,
            "elev_max": max(elevs) if elevs else None,
            "freq_min_common": min(fmins) if fmins else None,
            "freq_max_common": max(fmaxs) if fmaxs else None,
            "nfreq_min": min(nfreq_vals) if nfreq_vals else None,
            "nfreq_max": max(nfreq_vals) if nfreq_vals else None,
            "comp_counts": comp_counts,
            "tip_count": tip_count,
        }

    # ------------------------------------------------------------------
    # Rich output
    # ------------------------------------------------------------------

    def _rich_report(
        self,
        records: list[dict[str, Any]],
        *,
        detail: bool = False,
    ) -> None:
        sv = self._survey
        con = Console()

        # --- header panel ---
        n = sv["n_stations"]
        bbox = ""
        if sv["lat_min"] is not None:
            bbox = (
                f"Lat {sv['lat_min']:.2f}–{sv['lat_max']:.2f}°N  ·  "
                f"Lon {sv['lon_min']:.2f}–{sv['lon_max']:.2f}°E"
            )
            if sv["elev_min"] is not None:
                bbox += f"  ·  Elev {sv['elev_min']:.0f}–{sv['elev_max']:.0f} m"
        nf_str = ""
        if sv["nfreq_min"] is not None:
            nf_str = (
                f"{sv['nfreq_min']}–{sv['nfreq_max']} freq/station  ·  "
                f"{_fmt_freq(sv['freq_min_common'])} → "
                f"{_fmt_freq(sv['freq_max_common'])}"
            )

        hdr = Table.grid(padding=(0, 2))
        hdr.add_column(style="bold dim")
        hdr.add_column(style="white")
        hdr.add_row("Stations", str(n))
        if bbox:
            hdr.add_row("Coverage", bbox)
        if nf_str:
            hdr.add_row("Frequencies", nf_str)

        con.print()
        con.print(
            Panel(
                hdr,
                title="[bold bright_cyan]Survey Summary[/bold bright_cyan]",
                border_style="bright_cyan",
                expand=False,
            )
        )

        # --- per-station table ---
        tbl = Table(
            title=f"[dim]Stations ({len(records)}"
            + (f" of {n}" if len(records) < n else "")
            + ")[/dim]",
            border_style="blue",
            show_lines=False,
        )
        tbl.add_column("Station", style="bold cyan", no_wrap=True)
        tbl.add_column("Freq", style="dim", justify="right")
        for c in _COMPONENTS:
            tbl.add_column(c, justify="center", width=4)
        tbl.add_column("Tip", justify="center", width=4)
        tbl.add_column("ρ_xy Ω·m", justify="right")
        tbl.add_column("φ_xy °", justify="right")
        tbl.add_column("Cover", no_wrap=True)

        for r in records:
            comp_cells = [
                "[green]✓[/green]" if r["components"].get(c, False) else "[red]✗[/red]"
                for c in _COMPONENTS
            ]
            tip_cell = "[green]✓[/green]" if r["has_tipper"] else "[dim]—[/dim]"

            rho_s = (
                f"{r['rho_xy_mean']:.0f}±{r['rho_xy_std']:.0f}"
                if r["rho_xy_mean"] is not None
                else "—"
            )
            phi_s = (
                f"{r['phi_xy_mean']:.1f}±{r['phi_xy_std']:.1f}"
                if r["phi_xy_mean"] is not None
                else "—"
            )

            # coverage = average quality across present components
            qs = [
                q
                for c, q in r["quality"].items()
                if q is not None and r["components"].get(c, False)
            ]
            cover = _bar(float(np.mean(qs)), width=6) if qs else "——"

            tbl.add_row(
                r["name"],
                str(r["nfreq"]) if r["nfreq"] else "—",
                *comp_cells,
                tip_cell,
                rho_s,
                phi_s,
                cover,
            )

        con.print(tbl)

        # --- component availability ---
        avail = Table.grid(padding=(0, 1))
        avail.add_column(style="bold dim", width=6)
        avail.add_column(no_wrap=True)
        avail.add_column(style="dim", justify="right", width=10)
        avail.add_column(style="dim", justify="right", width=6)

        all_comps = list(_COMPONENTS) + ["Tipper"]
        counts = {c: sv["comp_counts"].get(c, 0) for c in _COMPONENTS}
        counts["Tipper"] = sv["tip_count"]

        for comp in all_comps:
            cnt = counts[comp]
            frac = cnt / n if n > 0 else 0.0
            avail.add_row(
                comp,
                _bar(frac, width=16),
                f"{cnt}/{n}",
                f"{frac * 100:.0f}%",
            )

        con.print()
        con.print(
            Panel(
                avail,
                title="[dim]Component Availability[/dim]",
                border_style="dim",
                expand=False,
            )
        )
        con.print()

    # ------------------------------------------------------------------
    # Plain-text output
    # ------------------------------------------------------------------

    def _plain_lines(
        self,
        records: list[dict[str, Any]],
        *,
        detail: bool = False,
    ) -> list[str]:
        sv = self._survey
        n = sv["n_stations"]
        lines = [
            "",
            f"Survey: {n} station(s)",
            "─" * 72,
        ]
        if sv["lat_min"] is not None:
            lines.append(
                f"  Lat {sv['lat_min']:.2f}–{sv['lat_max']:.2f}°N  "
                f"Lon {sv['lon_min']:.2f}–{sv['lon_max']:.2f}°E"
            )
        if sv["freq_min_common"] is not None:
            lines.append(
                f"  Freq {sv['nfreq_min']}–{sv['nfreq_max']}/station  "
                f"· {_fmt_freq(sv['freq_min_common'])} → "
                f"{_fmt_freq(sv['freq_max_common'])}"
            )
        lines.append("")

        hdr = (
            f"  {'Station':<12} {'Freq':>5}  "
            + "  ".join(f"{c:4}" for c in _COMPONENTS)
            + f"  {'Tip':4}  {'ρ_xy':>12}  {'φ_xy':>10}  Cover"
        )
        lines += [hdr, "  " + "─" * (len(hdr) - 2)]

        for r in records:
            comp_cells = "  ".join(
                f"{'✓' if r['components'].get(c, False) else '✗':4}"
                for c in _COMPONENTS
            )
            tip_cell = f"{'✓' if r['has_tipper'] else '—':4}"
            rho_s = (
                f"{r['rho_xy_mean']:.0f}±{r['rho_xy_std']:.0f} Ω·m"
                if r["rho_xy_mean"] is not None
                else "—"
            )
            phi_s = (
                f"{r['phi_xy_mean']:.1f}±{r['phi_xy_std']:.1f}°"
                if r["phi_xy_mean"] is not None
                else "—"
            )
            qs = [
                q
                for c, q in r["quality"].items()
                if q is not None and r["components"].get(c, False)
            ]
            cover = _bar(float(np.mean(qs)), width=6) if qs else "——"

            lines.append(
                f"  {r['name']:<12} {r['nfreq']:>5}  "
                + comp_cells
                + f"  {tip_cell}  {rho_s:>12}  {phi_s:>10}  {cover}"
            )

        lines += ["", "Component availability:"]
        for comp in list(_COMPONENTS) + ["Tipper"]:
            cnt = sv["comp_counts"].get(
                comp, sv["tip_count"] if comp == "Tipper" else 0
            )
            frac = cnt / n if n > 0 else 0.0
            lines.append(f"  {comp:<6}  {_bar(frac, 16)}  {cnt}/{n}  {frac * 100:.0f}%")
        lines.append("")
        return lines
