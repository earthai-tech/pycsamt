"""
pycsamt.agents.tipper_analysis
================================

:class:`TipperAnalysisAgent` — Dedicated tipper vector analysis.

Computes induction arrows (Wiese or Parkinson convention), tipper magnitude
and phase curves, and produces a map-view arrow plot alongside per-station
tipper vs period profiles.  Tipper vectors are a primary diagnostic for
3-D structure, crustal conductors, and coastal effects.

Wraps
-----
* :meth:`~pycsamt.core.base.MTBase.induction_arrows`
* :meth:`~pycsamt.core.base.MTBase.tipper_amp_phase`
* :class:`~pycsamt.z.tipper.Tipper`
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np

from ._base import AgentResult, BaseAgent

_SYSTEM_PROMPT = """\
You are an expert in magnetotelluric tipper analysis and 3-D structure interpretation.
Given a tipper analysis result, write 4-5 sentences that:
1. Describe the general tipper magnitude pattern across the survey (strong/weak, frequency dependence).
2. Identify stations with anomalously large tipper amplitudes and their likely cause.
3. Interpret the induction arrow direction(s) — do they point toward or away from a conductor?
4. Assess the 3-D character of the survey based on tipper consistency along the profile.
5. Recommend whether a 3-D inversion is warranted or whether 2-D is sufficient.
Reply in plain scientific English.
"""


class TipperAnalysisAgent(BaseAgent):
    """Analyse tipper vectors and plot induction arrows.

    Parameters
    ----------
    api_key, model, llm_provider : str
    convention : {'wiese', 'parkinson'}
        Arrow convention.  Wiese (default) points toward conductors in
        the real-part convention; Parkinson points toward them.
    use_imag : bool
        Use imaginary tipper parts for deeper structure (default False).
    period_ref : float or None
        Reference period (s) for the induction arrow map.
        ``None`` uses the geometric mean of available periods.

    Input keys
    ----------
    ``sites`` / ``path`` : Sites or str
    ``convention`` : str, optional
    ``use_imag`` : bool, optional
    ``period_ref`` : float, optional
    ``output_dir`` : str, optional

    Output data keys
    ----------------
    ``tipper_table``      pandas.DataFrame — per-(station, period)
    ``arrow_table``       pandas.DataFrame — induction arrows at period_ref
    ``period_ref``        float — period used for arrow map
    ``n_stations_with_tipper``  int
    ``figures``           dict
    ``figure_paths``      dict

    Examples
    --------
    >>> agent = TipperAnalysisAgent(convention='wiese')
    >>> r = agent.execute({"path": "/data/WILLY_EDIs"})
    >>> r["n_stations_with_tipper"]
    12
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        convention: str = "wiese",
        use_imag: bool = False,
        period_ref: float | None = None,
    ) -> None:
        super().__init__(
            "TipperAnalysisAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="pseudosection",
        )
        self.convention = convention.lower()
        self.use_imag   = use_imag
        self.period_ref = period_ref

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        try:
            from ..core.base import MTBase
        except ImportError as exc:
            return AgentResult.failed(
                f"pycsamt.core.base not available: {exc}", elapsed=time.time() - t0,
            )

        from ..emtools._core import (
            _get_z_block,
            _iter_items,
            _name,
            ensure_sites,
        )

        sites_raw = input_data.get("sites") or input_data.get("path")
        if sites_raw is None:
            return AgentResult.failed("No 'sites' or 'path'.", elapsed=time.time() - t0)
        try:
            sites = ensure_sites(sites_raw, verbose=0)
        except Exception as exc:
            return AgentResult.failed(str(exc), elapsed=time.time() - t0)

        convention = str(input_data.get("convention", self.convention)).lower()
        use_imag   = bool(input_data.get("use_imag",   self.use_imag))
        output_dir = input_data.get("output_dir")

        mtbase = MTBase()

        # ── collect tipper data ───────────────────────────────────────────
        records: list[dict] = []
        all_periods: list[float] = []

        for i, ed in enumerate(_iter_items(sites)):
            nm = _name(ed, i)
            _, z, fr = _get_z_block(ed)
            if fr is None or fr.size == 0:
                continue

            tip_raw = None
            try:
                tip_obj = getattr(ed, "Tip", None)
                if tip_obj is not None:
                    tip_raw = getattr(tip_obj, "tipper", None)
            except Exception:
                pass

            if tip_raw is None:
                warnings.append(f"{nm}: no tipper data.")
                continue

            try:
                t_arr = np.asarray(tip_raw, dtype=complex)
                if t_arr.ndim == 3 and t_arr.shape[1] == 1:
                    t_arr = t_arr[:, 0, :]  # → (n_freq, 2)
                if t_arr.shape != (len(fr), 2):
                    warnings.append(f"{nm}: unexpected tipper shape {t_arr.shape}.")
                    continue

                per = 1.0 / np.where(fr == 0, np.nan, fr)
                amp, phi = mtbase.tipper_amp_phase(t_arr, phase_unit="deg")
                ax_arr, ay_arr = mtbase.induction_arrows(
                    t_arr, convention=convention, use_imag=use_imag,
                )

                for fi in range(len(fr)):
                    if not np.isfinite(per[fi]):
                        continue
                    records.append({
                        "station":   nm,
                        "period_s":  float(per[fi]),
                        "freq_hz":   float(fr[fi]),
                        "Tx_re":     float(t_arr[fi, 0].real),
                        "Tx_im":     float(t_arr[fi, 0].imag),
                        "Ty_re":     float(t_arr[fi, 1].real),
                        "Ty_im":     float(t_arr[fi, 1].imag),
                        "amplitude": float(amp[fi]) if np.ndim(amp) > 0 else float(amp),
                        "phase_deg": float(phi[fi]) if np.ndim(phi) > 0 else float(phi),
                        "arrow_x":   float(ax_arr[fi]) if np.ndim(ax_arr) > 0 else float(ax_arr),
                        "arrow_y":   float(ay_arr[fi]) if np.ndim(ay_arr) > 0 else float(ay_arr),
                    })
                    all_periods.append(float(per[fi]))

            except Exception as exc:
                warnings.append(f"{nm}: tipper extraction failed: {exc}")

        n_sta_tip = len({r["station"] for r in records})

        if n_sta_tip == 0:
            return AgentResult.failed(
                "No tipper data found in any station.",
                hint="Ensure EDI files contain Tipper sections.",
                elapsed=time.time() - t0,
            )

        try:
            import pandas as pd
            tipper_table = pd.DataFrame(records)
        except ImportError:
            tipper_table = records  # fallback: plain list

        # ── choose reference period for arrow map ─────────────────────────
        period_ref = input_data.get("period_ref") or self.period_ref
        if period_ref is None and all_periods:
            period_ref = float(np.exp(np.nanmean(np.log(np.clip(all_periods, 1e-9, None)))))

        # ── build arrow table at reference period ─────────────────────────
        arrow_records: list[dict] = []
        if period_ref is not None and isinstance(tipper_table, list):
            rows = tipper_table
        elif period_ref is not None:
            rows = tipper_table.to_dict("records")
        else:
            rows = []

        station_names_ordered: list[str] = []
        seen: set[str] = set()
        for r in rows:
            if r["station"] not in seen:
                station_names_ordered.append(r["station"])
                seen.add(r["station"])

        if period_ref is not None:
            for sta in station_names_ordered:
                sta_rows = [r for r in rows if r["station"] == sta]
                if not sta_rows:
                    continue
                closest = min(sta_rows,
                              key=lambda r: abs(np.log10(max(r["period_s"], 1e-9))
                                                - np.log10(max(period_ref, 1e-9))))
                arrow_records.append({
                    "station": sta,
                    "period_s": closest["period_s"],
                    "arrow_x":  closest["arrow_x"],
                    "arrow_y":  closest["arrow_y"],
                    "amplitude":closest["amplitude"],
                })

        try:
            import pandas as pd
            arrow_table = pd.DataFrame(arrow_records)
        except ImportError:
            arrow_table = arrow_records

        # ── figures ───────────────────────────────────────────────────────
        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}

        try:
            fig_amp = _plot_tipper_amplitude(rows, station_names_ordered)
            if fig_amp is not None:
                figures["tipper_amplitude"] = fig_amp
                p = self._save_figure(fig_amp, output_dir, "tipper_amplitude",
                                      warnings_list=warnings)
                if p:
                    fig_paths["tipper_amplitude"] = p
        except Exception as exc:
            warnings.append(f"Tipper amplitude figure: {exc}")

        try:
            fig_arr = _plot_induction_arrows(
                arrow_records, period_ref, convention,
            )
            if fig_arr is not None:
                figures["induction_arrows"] = fig_arr
                p = self._save_figure(fig_arr, output_dir, "induction_arrows",
                                      warnings_list=warnings)
                if p:
                    fig_paths["induction_arrows"] = p
        except Exception as exc:
            warnings.append(f"Induction arrow figure: {exc}")

        try:
            fig_pseudo = _plot_tipper_pseudosection(rows, station_names_ordered)
            if fig_pseudo is not None:
                figures["tipper_pseudosection"] = fig_pseudo
                p = self._save_figure(fig_pseudo, output_dir, "tipper_pseudosection",
                                      warnings_list=warnings)
                if p:
                    fig_paths["tipper_pseudosection"] = p
        except Exception as exc:
            warnings.append(f"Tipper pseudosection: {exc}")

        # ── LLM interpretation ────────────────────────────────────────────
        interp: str | None = None
        if self.api_key and records:
            amps = [r["amplitude"] for r in records if np.isfinite(r["amplitude"])]
            mean_amp = float(np.mean(amps)) if amps else float("nan")
            max_tip  = max(arrow_records, key=lambda r: r["amplitude"],
                           default={}).get("station", "?") if arrow_records else "?"
            prompt = (
                f"Tipper analysis summary:\n"
                f"  Convention: {convention}, use_imag: {use_imag}\n"
                f"  Stations with tipper: {n_sta_tip}\n"
                f"  Reference period: {period_ref:.3f} s\n"
                f"  Mean tipper amplitude: {mean_amp:.3f}\n"
                f"  Station with max tipper: {max_tip}\n"
                f"  Warnings: {warnings[:3] if warnings else 'none'}\n\n"
                "Interpret the induction arrows and assess 3-D structure."
            )
            interp = self.query_llm(prompt, max_tokens=250)

        elapsed = time.time() - t0
        return AgentResult(
            status="success",
            summary=(
                f"Tipper analysis ({convention}): {n_sta_tip} stations. "
                f"Reference period: {period_ref:.3f} s. "
                f"{len(figures)} figures."
            ),
            data={
                "tipper_table":            tipper_table,
                "arrow_table":             arrow_table,
                "period_ref":              period_ref,
                "n_stations_with_tipper":  n_sta_tip,
                "convention":              convention,
                "figures":                 figures,
                "figure_paths":            fig_paths,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


# ── plot helpers ──────────────────────────────────────────────────────────────

def _plot_tipper_amplitude(rows: list[dict], station_order: list[str]) -> Any:
    """Tipper amplitude vs period, one curve per station."""
    import matplotlib.pyplot as plt

    if not rows:
        return None

    fig, ax = plt.subplots(figsize=(10, 5))
    for sta in station_order:
        sta_rows = sorted([r for r in rows if r["station"] == sta],
                          key=lambda r: r["period_s"])
        if not sta_rows:
            continue
        per = [r["period_s"] for r in sta_rows]
        amp = [r["amplitude"] for r in sta_rows]
        ax.semilogx(per, amp, "-o", markersize=3, linewidth=0.9, label=sta)

    ax.set_xlabel("Period (s)", fontsize=9)
    ax.set_ylabel("|T|  (dimensionless)", fontsize=9)
    ax.set_title("Tipper amplitude vs period", fontsize=10, fontweight="bold")
    ax.tick_params(labelsize=8)
    if len(station_order) <= 15:
        ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    return fig


def _plot_induction_arrows(
    arrow_records: list[dict],
    period_ref: float | None,
    convention: str,
) -> Any:
    """Map-view induction arrow plot at the reference period."""
    import matplotlib.pyplot as plt

    if not arrow_records:
        return None

    fig, ax = plt.subplots(figsize=(8, 6))

    xs = [i for i in range(len(arrow_records))]  # use station index as x-position
    [0.0] * len(arrow_records)

    for xi, rec in zip(xs, arrow_records):
        ax_val, ay_val = rec["arrow_x"], rec["arrow_y"]
        scale = 2.0
        ax.annotate(
            "", xy=(xi + ax_val * scale, ay_val * scale),
            xytext=(xi, 0.0),
            arrowprops=dict(arrowstyle="->", color="#c0392b", lw=1.5),
        )
        ax.text(xi, -0.15, rec["station"], ha="center", va="top",
                fontsize=6.5, rotation=90)

    ax.axhline(0, color="k", lw=0.5, ls="--")
    ax.set_xlim(-0.5, len(arrow_records) - 0.5)
    ax.set_ylim(-0.8, 1.0)
    ax.set_xlabel("Station index", fontsize=9)
    ax.set_ylabel("Tipper (normalised)", fontsize=9)
    per_str = f"{period_ref:.3f} s" if period_ref is not None else "?"
    ax.set_title(
        f"Induction arrows — {convention.capitalize()} convention  "
        f"(T = {per_str})",
        fontsize=10, fontweight="bold",
    )
    ax.tick_params(labelsize=8)
    fig.tight_layout()
    return fig


def _plot_tipper_pseudosection(
    rows: list[dict],
    station_order: list[str],
) -> Any:
    """Colour-coded tipper amplitude pseudosection (station × log-period)."""
    import matplotlib.pyplot as plt

    if not rows or not station_order:
        return None

    periods = sorted({r["period_s"] for r in rows if np.isfinite(r["period_s"])})
    if not periods:
        return None

    mat = np.full((len(periods), len(station_order)), np.nan)
    per_idx = {p: i for i, p in enumerate(periods)}
    sta_idx = {s: i for i, s in enumerate(station_order)}

    for r in rows:
        pi = per_idx.get(r["period_s"])
        si = sta_idx.get(r["station"])
        if pi is not None and si is not None and np.isfinite(r["amplitude"]):
            mat[pi, si] = r["amplitude"]

    fig, ax = plt.subplots(figsize=(max(8, len(station_order) * 0.5), 5))
    log_per = np.log10(np.clip(periods, 1e-9, None))

    vv = mat[np.isfinite(mat)]
    vmax = float(np.percentile(vv, 95)) if vv.size else 1.0

    im = ax.imshow(
        mat,
        aspect="auto",
        origin="upper",
        extent=(-0.5, len(station_order) - 0.5, log_per[-1], log_per[0]),
        cmap="hot_r",
        vmin=0.0, vmax=vmax,
        interpolation="bilinear",
    )
    ax.set_xticks(range(len(station_order)))
    ax.set_xticklabels(station_order, rotation=90, fontsize=6.5)
    ax.set_ylabel("log₁₀ Period (s)", fontsize=9)
    ax.set_title("Tipper amplitude pseudosection", fontsize=10, fontweight="bold")
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    div = make_axes_locatable(ax)
    cax = div.append_axes("right", size="2%", pad=0.05)
    plt.colorbar(im, cax=cax, label="|T|")
    fig.tight_layout()
    return fig


__all__ = ["TipperAnalysisAgent"]
