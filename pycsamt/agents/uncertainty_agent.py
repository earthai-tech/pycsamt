# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.agents.uncertainty_agent
==================================

:class:`UncertaintyCalibrationAgent` — Learned per-cell error-floor
re-estimation for MT impedance data.

Wraps :class:`~pycsamt.ai.processing.uncertainty.UncertaintyCalibrator`:
regresses a field-processing-consistent fractional error from the same
signal-quality features :class:`~pycsamt.ai.processing.qc.EMQCScorer`
extracts, then rescales each station's ``Z.z_err`` toward that
calibrated level -- ``Z.z`` itself is never touched, only its stated
uncertainty.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np

from ._base import AgentResult, BaseAgent

_SYSTEM_PROMPT = """\
You are an expert in MT/AMT impedance error-floor calibration and inversion
weighting.
Given an uncertainty-calibration result summary, write 3-4 sentences that:
1. State how the calibrated fractional error compares to the field-processing
   error already on the data (systematically higher, lower, or mixed).
2. Identify which stations show the largest recalibration and whether that
   points to a genuine data-quality issue or a processing-software artefact.
3. Note the held-out validation RMSE/R^2 and whether the calibration is
   trustworthy enough to use directly in an inversion error floor.
4. Recommend whether to apply the calibrated errors as-is or treat them as
   an additional diagnostic alongside the original error floor.
Reply in plain English.
"""


class UncertaintyCalibrationAgent(BaseAgent):
    """Recalibrate MT impedance error floors from signal-quality features.

    Parameters
    ----------
    api_key, model, llm_provider : str
    hidden : tuple of int
        Hidden-layer widths for the calibration network (default
        ``(32, 16)``).
    epochs : int
        Training epochs (default 150).
    val_frac : float
        Held-out fraction used for the quantitative RMSE/R^2 check
        (default 0.2).

    Input keys
    ----------
    ``sites`` / ``path`` : Sites or str
    ``output_dir`` : str, optional

    Output data keys
    ----------------
    ``calibrated_sites``  Sites with Z.z_err rescaled (Z.z untouched)
    ``table``             pandas DataFrame — full per-(station, freq)
                           z_err_frac / z_err_frac_calibrated
    ``rmse``, ``r2``       float — held-out validation metrics
    ``figures``            dict
    ``figure_paths``       dict
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        hidden: tuple[int, ...] = (32, 16),
        epochs: int = 150,
        val_frac: float = 0.2,
    ) -> None:
        super().__init__(
            "UncertaintyCalibrationAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="pseudosection",
        )
        self.hidden = hidden
        self.epochs = epochs
        self.val_frac = val_frac

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        try:
            from ..ai.processing import (
                UncertaintyCalibrator,
                build_uncertainty_features_table,
            )
            from ..backends import get_backend_instance

            if get_backend_instance() is None:
                raise ImportError("No DL backend.")
        except ImportError as exc:
            return AgentResult.failed(
                f"UncertaintyCalibrationAgent requires PyTorch or TensorFlow: {exc}",
                hint="pip install torch  or  pip install tensorflow",
                elapsed=time.time() - t0,
            )

        from ..emtools._core import ensure_sites

        sites_raw = input_data.get("sites") or input_data.get("path")
        if sites_raw is None:
            return AgentResult.failed(
                "No 'sites' or 'path'.", elapsed=time.time() - t0
            )
        try:
            sites = ensure_sites(sites_raw, verbose=0)
        except Exception as exc:
            return AgentResult.failed(str(exc), elapsed=time.time() - t0)

        output_dir = input_data.get("output_dir")

        try:
            feats = build_uncertainty_features_table(sites)
        except Exception as exc:
            return AgentResult.failed(
                f"build_uncertainty_features_table failed: {exc}",
                elapsed=time.time() - t0,
            )
        if len(feats) < 5:
            return AgentResult.failed(
                "Need >= 5 (station, freq) rows to calibrate.",
                elapsed=time.time() - t0,
            )

        # UncertaintyCalibrator.fit() silently drops non-finite/non-positive
        # target rows internally (see its own docstring), so training is
        # robust either way -- but a NaN row landing in *our* held-out split
        # propagates straight into the RMSE/R^2 below unless dropped here
        # first, exactly like fit()'s own validity criterion.
        n_before = len(feats)
        feats = feats[np.isfinite(feats["z_err_frac"]) & (feats["z_err_frac"] > 0)]
        feats = feats.dropna().reset_index(drop=True)
        n_dropped = n_before - len(feats)
        if n_dropped:
            warnings.append(
                f"Dropped {n_dropped} row(s) with non-finite features/target "
                "before validation split."
            )
        if len(feats) < 5:
            return AgentResult.failed(
                "Fewer than 5 usable (finite) feature rows after dropping "
                f"{n_dropped} non-finite row(s).",
                elapsed=time.time() - t0,
            )

        # ── held-out validation split (row-level, seeded) ───────────────────
        rng = np.random.default_rng(0)
        n = len(feats)
        idx = rng.permutation(n)
        n_val = max(1, int(n * self.val_frac))
        val_idx, train_idx = idx[:n_val], idx[n_val:]
        train_df = feats.iloc[train_idx].reset_index(drop=True)
        val_df = feats.iloc[val_idx].reset_index(drop=True)

        try:
            cal_val = UncertaintyCalibrator(hidden=self.hidden)
            cal_val.fit(train_df, epochs=self.epochs, seed=0, verbose=False)
            y_true = val_df["z_err_frac"].to_numpy()
            y_pred = cal_val.transform(val_df)
            rmse = float(np.sqrt(np.mean((y_true - y_pred) ** 2)))
            ss_res = float(np.sum((y_true - y_pred) ** 2))
            ss_tot = float(np.sum((y_true - y_true.mean()) ** 2)) + 1e-24
            r2 = 1.0 - ss_res / ss_tot
        except Exception as exc:
            return AgentResult.failed(
                f"UncertaintyCalibrator validation fit failed: {exc}",
                elapsed=time.time() - t0,
            )

        # ── deployed calibrator: fit on the full table ──────────────────────
        try:
            cal = UncertaintyCalibrator(hidden=self.hidden)
            cal.fit(feats, epochs=self.epochs, seed=0, verbose=False)
            table = cal.predict_table(sites)
            calibrated_sites = cal.apply(sites, inplace=False)
        except Exception as exc:
            return AgentResult.failed(
                f"UncertaintyCalibrator fit/apply failed: {exc}",
                elapsed=time.time() - t0,
            )

        # ── figures ───────────────────────────────────────────────────────
        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}
        try:
            from ..ai.processing import plot_uncertainty_summary

            fig = plot_uncertainty_summary(
                table,
                y_true,
                y_pred,
                history=cal,
                suptitle="UncertaintyCalibrator -- error-floor recalibration",
            )
            figures["uncertainty_summary"] = fig
            p = self._save_figure(
                fig, output_dir, "uncertainty_summary", warnings_list=warnings
            )
            if p:
                fig_paths["uncertainty_summary"] = p
        except Exception as exc:
            warnings.append(f"Uncertainty summary figure: {exc}")

        # ── LLM interpretation ──────────────────────────────────────────────
        interp: str | None = None
        if self.llm_available:
            agg = table.groupby("station")[
                ["z_err_frac", "z_err_frac_calibrated"]
            ].median()
            top5 = agg.assign(
                delta=(agg["z_err_frac_calibrated"] - agg["z_err_frac"]).abs()
            ).nlargest(5, "delta")
            prompt = (
                f"Uncertainty calibration summary:\n"
                f"  Rows calibrated: {len(table)}\n"
                f"  Held-out RMSE: {rmse:.4f}, R^2: {r2:.3f}\n"
                f"  Median field z_err_frac: {feats['z_err_frac'].median():.4f}\n"
                f"  Largest per-station recalibration (median, abs delta):\n"
                f"{top5.to_string()}\n\n"
                "Assess calibration quality and recommend follow-up."
            )
            interp = self.query_llm(prompt, max_tokens=220)

        elapsed = time.time() - t0
        return AgentResult(
            status="success",
            summary=(
                f"Recalibrated error floor for {len(table)} row(s) "
                f"(held-out RMSE={rmse:.4f}, R^2={r2:.3f}). "
                f"{len(figures)} figure(s)."
            ),
            data={
                "calibrated_sites": calibrated_sites,
                "table": table,
                "rmse": rmse,
                "r2": r2,
                "figures": figures,
                "figure_paths": fig_paths,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


__all__ = ["UncertaintyCalibrationAgent"]
