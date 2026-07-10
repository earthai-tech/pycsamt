"""
pycsamt.agents.anomaly_agent
=============================

:class:`AnomalyDetectionAgent` — Unsupervised MT data anomaly detection.

Wraps :class:`~pycsamt.ai.processing.anomaly.AnomalyDetector`:

A convolutional autoencoder (CAE) is trained on the clean portions of the
dataset in an unsupervised manner.  Observations whose reconstruction error
exceeds the *threshold_percentile* are flagged as anomalies.

This agent detects anomalous (station, frequency) cells that classical
rule-based QC may miss — e.g. coherent noise from a nearby source, subtle
sensor drift, or 3-D scattering that deviates from the expected frequency
dependence.

Requires PyTorch **or** TensorFlow.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np

from ._base import AgentResult, BaseAgent

_SYSTEM_PROMPT = """\
You are an expert in unsupervised anomaly detection for MT/AMT data.
Given an anomaly detection result, write 3-4 sentences that:
1. State how many observations were flagged as anomalous and their distribution.
2. Identify which stations or frequency bands are most affected.
3. Diagnose the likely source (powerline harmonics, near-field, 3-D, instrument).
4. Recommend whether to mask flagged data or apply targeted filtering.
Reply in plain English.
"""


class AnomalyDetectionAgent(BaseAgent):
    """Detect anomalous (station, frequency) observations in MT data.

    Parameters
    ----------
    api_key, model, llm_provider : str
    threshold_percentile : float
        Percentile of reconstruction errors used as the flagging threshold
        (default 95 — top 5 % are anomalies).
    latent_dim : int
        CAE latent space dimension (default 32).
    epochs : int
        Training epochs (default 50).

    Input keys
    ----------
    ``sites`` / ``path`` : Sites or str
    ``output_dir`` : str, optional

    Output data keys
    ----------------
    ``anomaly_scores``    ndarray — per-(station, freq) reconstruction error
    ``flags``             ndarray bool — True = anomalous
    ``flag_table``        pandas DataFrame {station, freq, score, flagged}
    ``n_flagged``         int
    ``flagged_stations``  list[str]
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
        threshold_percentile: float = 95.0,
        latent_dim: int = 32,
        epochs: int = 50,
    ) -> None:
        super().__init__(
            "AnomalyDetectionAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="pseudosection",
        )
        self.threshold_percentile = threshold_percentile
        self.latent_dim           = latent_dim
        self.epochs               = epochs

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        # ── backend check ──────────────────────────────────────────────────────
        try:
            from ..ai.processing.anomaly import (
                AnomalyDetector,
            )
            from ..backends import get_backend_instance
            if get_backend_instance() is None:
                raise ImportError("No DL backend.")
        except ImportError as exc:
            return AgentResult.failed(
                f"AnomalyDetectionAgent requires PyTorch or TensorFlow: {exc}",
                hint="pip install torch  or  pip install tensorflow",
                elapsed=time.time() - t0,
            )

        from ..emtools._core import (
            _get_z_block,
            _iter_items,
            _name,
            ensure_sites,
        )
        from .ai_inversion import _z_to_features

        sites_raw = input_data.get("sites") or input_data.get("path")
        if sites_raw is None:
            return AgentResult.failed("No 'sites' or 'path'.", elapsed=time.time()-t0)
        try:
            sites = ensure_sites(sites_raw, verbose=0)
        except Exception as exc:
            return AgentResult.failed(str(exc), elapsed=time.time()-t0)

        output_dir = input_data.get("output_dir")
        freqs = np.logspace(-4, 3, 40)

        # ── collect features ───────────────────────────────────────────────────
        station_names: list[str] = []
        feat_list: list[np.ndarray] = []

        for i, ed in enumerate(_iter_items(sites)):
            nm = _name(ed, i)
            _, z, fr = _get_z_block(ed)
            if z is None:
                continue
            feat = _z_to_features(z, fr, freqs)
            if feat is None:
                warnings.append(f"{nm}: skipped.")
                continue
            station_names.append(nm)
            feat_list.append(feat)   # (n_freqs, n_comp)

        if len(feat_list) < 3:
            return AgentResult.failed(
                "Need ≥ 3 stations for anomaly detection.",
                elapsed=time.time() - t0,
            )

        n_sta = len(station_names)
        # X: (n_obs=n_sta, n_features=n_freqs*n_comp)
        n_freq, n_comp = feat_list[0].shape
        X = np.stack(feat_list, axis=0).reshape(n_sta, -1).astype(np.float32)
        n_feat = X.shape[1]

        # ── fit AnomalyDetector ────────────────────────────────────────────────
        try:
            detector = AnomalyDetector(
                n_features=n_feat,
                latent_dim=self.latent_dim,
                threshold_percentile=self.threshold_percentile,
            )
            detector.fit(X, epochs=self.epochs, batch_size=max(4, n_sta // 2),
                         verbose=False)
        except Exception as exc:
            return AgentResult.failed(
                f"AnomalyDetector.fit failed: {exc}", elapsed=time.time() - t0,
            )

        # ── score + flag ───────────────────────────────────────────────────────
        try:
            scores = detector.transform(X)    # (n_sta,) reconstruction errors
            flags  = detector.flag_anomalies(X)  # (n_sta,) bool
        except Exception as exc:
            return AgentResult.failed(
                f"AnomalyDetector.transform failed: {exc}",
                elapsed=time.time() - t0,
            )

        flagged_stations = [nm for nm, f in zip(station_names, flags) if f]
        n_flagged = int(np.sum(flags))

        # ── flag table ─────────────────────────────────────────────────────────
        try:
            import pandas as pd
            flag_table = pd.DataFrame({
                "station": station_names,
                "score":   scores,
                "flagged": flags.astype(bool),
            })
        except Exception:
            flag_table = None

        # ── figures ───────────────────────────────────────────────────────────
        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}

        try:
            import matplotlib.pyplot as plt
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4))

            # anomaly score bar chart
            colors = ["#e74c3c" if f else "#2ecc71" for f in flags]
            ax1.bar(range(n_sta), scores, color=colors, edgecolor="none", alpha=0.8)
            thresh = np.percentile(scores, self.threshold_percentile)
            ax1.axhline(thresh, color="#e74c3c", lw=1.5, ls="--",
                        label=f"threshold ({self.threshold_percentile:.0f}th pct)")
            ax1.set_xticks(range(n_sta))
            ax1.set_xticklabels(station_names, rotation=90, fontsize=6.5)
            ax1.set_ylabel("Reconstruction error", fontsize=9)
            ax1.set_title("Anomaly scores per station", fontsize=9, fontweight="bold")
            ax1.legend(fontsize=8)
            ax1.tick_params(labelsize=8)

            # flag map on pseudosection
            from ..emtools.inspect import pseudosection as _ps
            try:
                ax2 = _ps(sites, quantity="rho_xy", ax=ax2)
                # overlay flagged stations
                for si, (nm, flag) in enumerate(zip(station_names, flags)):
                    if flag:
                        ax2.axvline(si, color="red", lw=0.8, alpha=0.4)
                ax2.set_title("ρa section  (red = flagged)", fontsize=9, fontweight="bold")
            except Exception:
                ax2.text(0.5, 0.5, "pseudosection unavailable", ha="center")

            fig.suptitle(
                f"Anomaly detection — {n_flagged}/{n_sta} flagged "
                f"({100*n_flagged/max(n_sta,1):.0f}%)",
                fontsize=10, fontweight="bold",
            )
            fig.tight_layout()
            figures["anomaly_map"] = fig
            p = self._save_figure(fig, output_dir, "anomaly_detection",
                                  warnings_list=warnings)
            if p:
                fig_paths["anomaly_map"] = p
        except Exception as exc:
            warnings.append(f"Anomaly figure: {exc}")

        # ── LLM interpretation ────────────────────────────────────────────────
        interp: str | None = None
        if self.api_key:
            top5 = sorted(zip(station_names, scores),
                          key=lambda x: -x[1])[:5]
            prompt = (
                f"Anomaly detection summary:\n"
                f"  Stations: {n_sta}, flagged: {n_flagged} "
                f"({100*n_flagged/max(n_sta,1):.0f}%)\n"
                f"  Threshold percentile: {self.threshold_percentile:.0f}\n"
                f"  Highest anomaly scores: "
                f"{[(nm, f'{sc:.3f}') for nm, sc in top5]}\n"
                f"  Flagged: {flagged_stations}\n\n"
                "Diagnose the anomaly sources and recommend remediation."
            )
            interp = self.query_llm(prompt, max_tokens=200)

        elapsed = time.time() - t0
        return AgentResult(
            status="success",
            summary=(
                f"Anomaly detection: {n_flagged}/{n_sta} stations flagged "
                f"({100*n_flagged/max(n_sta,1):.0f}%). "
                f"{len(figures)} figures."
            ),
            data={
                "anomaly_scores":   scores,
                "flags":            flags,
                "flag_table":       flag_table,
                "n_flagged":        n_flagged,
                "flagged_stations": flagged_stations,
                "figures":          figures,
                "figure_paths":     fig_paths,
                "sites":            sites,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


__all__ = ["AnomalyDetectionAgent"]
