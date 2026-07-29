"""
pycsamt.agents.inversion_eval
==============================

:class:`InversionEvaluationAgent` — Evaluate inversion result quality.

Computes:
* Per-station RMS misfit between observed and model-predicted responses.
* Residual phase tensor section (data PT minus model PT).
* Misfit pseudosection (normalised ρa residuals).

Uses :func:`~pycsamt.emtools.inspect.plot_station_response` for per-station
response overlay and the phase tensor pipeline for residual PT.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np

from ._base import AgentResult, BaseAgent

_SYSTEM_PROMPT = """\
You are an expert MT inversion quality assessor.
Given a misfit summary, write 3–4 sentences that:
1. State whether the inversion converged acceptably (RMS 0.8–1.5 = good).
2. Identify stations or period bands with elevated misfit.
3. Diagnose likely causes (3-D effects, noise, model inadequacy).
4. Recommend whether to re-run with adjusted regularisation.
Reply in plain English.
"""


class InversionEvaluationAgent(BaseAgent):
    """Evaluate inversion quality: RMS, residual PT, misfit section.

    Input keys
    ----------
    ``sites_obs`` / ``path_obs`` : Sites or str — observed data
    ``sites_mod`` / ``path_mod`` : Sites or str — model-predicted responses
    ``output_dir`` : str, optional
    ``component`` : str — default ``"xy"``

    Output data keys
    ----------------
    ``rms_per_station``   dict {station: rms}
    ``rms_global``        float
    ``residual_pt_table`` pandas DataFrame
    ``figures``           dict
    ``figure_paths``      dict
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
    ) -> None:
        super().__init__(
            "InversionEvaluationAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="inversion",
        )

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        from ..emtools._core import (
            _get_z_block,
            _iter_items,
            _name,
            ensure_sites,
        )
        from ..emtools.tensor import build_phase_tensor_table

        obs_raw = input_data.get("sites_obs") or input_data.get("path_obs")
        mod_raw = input_data.get("sites_mod") or input_data.get("path_mod")
        component = str(input_data.get("component", "xy")).lower()
        output_dir = input_data.get("output_dir")
        ri, ci = (0, 1) if component == "xy" else (1, 0)

        if obs_raw is None:
            return AgentResult.failed(
                "No 'sites_obs' or 'path_obs' provided.",
                elapsed=time.time() - t0,
            )
        try:
            sites_obs = ensure_sites(obs_raw, verbose=0)
        except Exception as exc:
            return AgentResult.failed(str(exc), elapsed=time.time() - t0)

        if mod_raw is None:
            warnings.append("No model response provided; RMS cannot be computed.")
            sites_mod = None
        else:
            try:
                sites_mod = ensure_sites(mod_raw, verbose=0)
            except Exception as exc:
                warnings.append(f"Could not load model response: {exc}")
                sites_mod = None

        # ── per-station RMS ───────────────────────────────────────────────────
        rms_per: dict[str, float] = {}
        if sites_mod is not None:
            obs_dict = {_name(ed, i): ed for i, ed in enumerate(_iter_items(sites_obs))}
            for i, ed_m in enumerate(_iter_items(sites_mod)):
                nm = _name(ed_m, i)
                ed_o = obs_dict.get(nm)
                if ed_o is None:
                    continue
                try:
                    _, z_o, fr_o = _get_z_block(ed_o)
                    _, z_m, fr_m = _get_z_block(ed_m)
                    if z_o is None or z_m is None:
                        continue
                    rho_o = (0.2 / np.where(fr_o == 0, np.nan, fr_o)) * np.abs(
                        z_o[:, ri, ci]
                    ) ** 2
                    rho_m = (0.2 / np.where(fr_m == 0, np.nan, fr_m)) * np.abs(
                        z_m[:, ri, ci]
                    ) ** 2
                    n = min(len(rho_o), len(rho_m))
                    mask = (
                        np.isfinite(rho_o[:n])
                        & (rho_o[:n] > 0)
                        & np.isfinite(rho_m[:n])
                        & (rho_m[:n] > 0)
                    )
                    if mask.sum() < 2:
                        continue
                    res = np.log10(rho_o[:n][mask]) - np.log10(rho_m[:n][mask])
                    rms_per[nm] = float(np.sqrt(np.mean(res**2)))
                except Exception as exc:
                    warnings.append(f"RMS for {nm}: {exc}")

        rms_global = float(np.mean(list(rms_per.values()))) if rms_per else np.nan

        # ── residual PT ───────────────────────────────────────────────────────
        residual_pt = None
        if sites_mod is not None:
            try:
                pt_obs = build_phase_tensor_table(sites_obs, verbose=0)
                pt_mod = build_phase_tensor_table(sites_mod, verbose=0)
                if not pt_obs.empty and not pt_mod.empty:
                    merged = pt_obs.merge(
                        pt_mod,
                        on=["station", "period"],
                        suffixes=("_obs", "_mod"),
                    )
                    for col in ("skew", "ellipt", "theta"):
                        if (
                            f"{col}_obs" in merged.columns
                            and f"{col}_mod" in merged.columns
                        ):
                            merged[f"d_{col}"] = (
                                merged[f"{col}_obs"] - merged[f"{col}_mod"]
                            )
                    residual_pt = merged
            except Exception as exc:
                warnings.append(f"Residual PT computation: {exc}")

        # ── figures ───────────────────────────────────────────────────────────
        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}

        if sites_mod is not None:
            try:
                from ..emtools.inspect import (
                    plot_station_response,
                )

                fig_resp = plot_station_response(
                    sites_obs,
                    sites_model=sites_mod,
                )
                if fig_resp is not None:
                    f = (
                        fig_resp
                        if hasattr(fig_resp, "savefig")
                        else (
                            fig_resp.get_figure()
                            if hasattr(fig_resp, "get_figure")
                            else None
                        )
                    )
                    if f is not None:
                        figures["station_response"] = f
                        p = self._save_figure(
                            f,
                            output_dir,
                            "inversion_response",
                            warnings_list=warnings,
                        )
                        if p:
                            fig_paths["station_response"] = p
            except Exception as exc:
                warnings.append(f"plot_station_response: {exc}")

        # ── LLM interpretation ────────────────────────────────────────────────
        interp: str | None = None
        if self.api_key and not np.isnan(rms_global):
            bad = [f"{st}={v:.2f}" for st, v in rms_per.items() if v > 1.5]
            prompt = (
                f"Inversion evaluation:\n"
                f"  Global RMS: {rms_global:.3f}\n"
                f"  Stations with RMS > 1.5: {bad or 'none'}\n"
                f"  Warnings: {warnings[:3] if warnings else 'none'}\n\n"
                "Assess inversion quality and recommend next steps."
            )
            interp = self.query_llm(prompt, max_tokens=200)

        elapsed = time.time() - t0
        return AgentResult(
            status="success" if not np.isnan(rms_global) else "needs_review",
            summary=(
                f"Inversion evaluation: global RMS = {rms_global:.3f}"
                if not np.isnan(rms_global)
                else "Inversion evaluation: no model response provided."
            ),
            data={
                "rms_per_station": rms_per,
                "rms_global": rms_global,
                "residual_pt_table": residual_pt,
                "figures": figures,
                "figure_paths": fig_paths,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


__all__ = ["InversionEvaluationAgent"]
