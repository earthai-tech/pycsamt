"""
pycsamt.agents.ensemble_agent
==============================

:class:`EnsembleAgent` — Ensemble 1-D inversion with uncertainty quantification.

Wraps :class:`~pycsamt.ai.inversion.ensemble.EnsembleInverter`:

* Trains *N* independent :class:`~pycsamt.ai.inversion.inv1d.EMInverter1D`
  models on a shared synthetic dataset using different random seeds.
* Predicts **mean** resistivity and **uncertainty** (std / confidence intervals)
  for every observed station.
* Optionally **calibrates** the intervals using a held-out conformal set.
* Reports empirical **coverage** (fraction of true values inside the interval)
  as a reliability metric.

The outputs feed directly into the :class:`ReportAgent` and are used as a
rigorous uncertainty-aware alternative to single-model AI inversion.

Requires PyTorch **or** TensorFlow.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np

from ._base import AgentResult, BaseAgent
from .ai_inversion import _default_thicknesses, _z_to_features

_SYSTEM_PROMPT = """\
You are an expert in Bayesian and ensemble methods for geophysical inversion.
Given an ensemble inversion result with uncertainty quantification, write 4-5
sentences that:
1. Describe the ensemble configuration (N models, architecture, training data).
2. State the prediction quality (mean RMS, uncertainty magnitude).
3. Assess the calibration: are the confidence intervals reliable?
4. Identify depth ranges or stations where uncertainty is largest.
5. Recommend whether the uncertainty is small enough for geological interpretation.
Reply in plain scientific English.
"""


class EnsembleAgent(BaseAgent):
    """Ensemble 1-D MT inversion with uncertainty bands.

    Parameters
    ----------
    api_key, model, llm_provider : str
    n_estimators : int
        Number of independent models in the ensemble (default 5).
    arch : str
        Network architecture for each estimator (default ``"resnet"``).
    n_layers : int
        Number of model layers (default 5).
    n_train_samples : int
        Synthetic training samples per estimator (default 2 000).
    epochs : int
        Training epochs per estimator (default 30).
    calibrate : bool
        Apply conformal calibration using 20 % of training data (default True).

    Input keys
    ----------
    ``sites`` / ``path`` : Sites or str
    ``output_dir`` : str, optional

    Output data keys
    ----------------
    ``ensemble``            EnsembleInverter
    ``pred_mean``           dict {station: ndarray}  — mean log₁₀ ρ
    ``pred_std``            dict {station: ndarray}  — std of log₁₀ ρ
    ``pred_lo``             dict {station: ndarray}  — 5th-percentile bound
    ``pred_hi``             dict {station: ndarray}  — 95th-percentile bound
    ``coverage``            float  — empirical 90 % interval coverage
    ``rms_global``          float
    ``figures``             dict
    ``figure_paths``        dict

    Examples
    --------
    >>> agent  = EnsembleAgent(n_estimators=3, epochs=20)
    >>> result = agent.execute({"path": "/data/L22PLT", "output_dir": "/out/ens"})
    >>> result["coverage"]   # should be ≈ 0.90 after calibration
    0.88
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        n_estimators: int = 5,
        arch: str = "resnet",
        n_layers: int = 5,
        n_train_samples: int = 2_000,
        epochs: int = 30,
        calibrate: bool = True,
    ) -> None:
        super().__init__(
            "EnsembleAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="inversion",
        )
        self.n_estimators = n_estimators
        self.arch = arch
        self.n_layers = n_layers
        self.n_train_samples = n_train_samples
        self.epochs = epochs
        self.calibrate = calibrate

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        # ── backend check ──────────────────────────────────────────────────────
        try:
            from ..ai.inversion.ensemble import (
                EnsembleInverter,
            )
            from ..ai.inversion.inv1d import EMInverter1D
            from ..backends import get_backend_instance
            from ..forward.batch import generate_dataset

            if get_backend_instance() is None:
                raise ImportError("No DL backend (torch / tensorflow).")
        except ImportError as exc:
            return AgentResult.failed(
                f"EnsembleAgent requires PyTorch or TensorFlow: {exc}",
                hint="pip install torch  or  pip install tensorflow",
                elapsed=time.time() - t0,
            )

        from ..emtools._core import (
            _get_z_block,
            _iter_items,
            _name,
            ensure_sites,
        )

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
        import os

        if output_dir:
            os.makedirs(output_dir, exist_ok=True)

        freqs = np.logspace(-4, 3, 40)

        # ── generate synthetic dataset ────────────────────────────────────────
        self._log.info(
            "Generating %d synthetic samples…", self.n_train_samples
        )
        try:
            dataset = generate_dataset(
                solver="mt1d",
                n_samples=self.n_train_samples,
                freqs=freqs,
                n_layers=self.n_layers,
                noise_level=0.03,
                seed=42,
                n_jobs=1,
                verbose=False,
            )
            X_all, y_all = dataset.X, dataset.y

            # split: 80 % train, 20 % calibration
            n_cal = max(50, int(len(X_all) * 0.2))
            X_train, y_train = X_all[n_cal:], y_all[n_cal:]
            X_cal, y_cal = X_all[:n_cal], y_all[:n_cal]
        except Exception as exc:
            return AgentResult.failed(
                f"Dataset generation failed: {exc}",
                elapsed=time.time() - t0,
            )

        # ── build base estimator + ensemble ───────────────────────────────────
        self._log.info(
            "Training ensemble of %d × %s models…",
            self.n_estimators,
            self.arch,
        )
        try:
            base = EMInverter1D(
                arch=self.arch,
                n_layers=self.n_layers,
                solver="mt1d",
                include_phase=True,
            )
            ens = EnsembleInverter(
                base_estimator=base,
                n_estimators=self.n_estimators,
            )
            ens.fit(
                X_train,
                y_train,
                epochs=self.epochs,
                batch_size=min(256, len(X_train) // 4),
                patience=max(5, self.epochs // 5),
                verbose=False,
            )
        except Exception as exc:
            return AgentResult.failed(
                f"Ensemble training failed: {exc}",
                elapsed=time.time() - t0,
            )

        # ── calibrate ─────────────────────────────────────────────────────────
        if self.calibrate:
            try:
                ens.calibrate(X_cal, y_cal, alpha=0.10)
                self._log.info("Conformal calibration done.")
            except Exception as exc:
                warnings.append(f"Calibration failed: {exc}")

        # ── coverage on held-out set ──────────────────────────────────────────
        coverage = np.nan
        try:
            coverage = float(ens.coverage(X_cal, y_cal))
        except Exception as exc:
            warnings.append(f"Coverage computation: {exc}")

        # ── predict on observed stations ──────────────────────────────────────
        pred_mean: dict[str, np.ndarray] = {}
        pred_std: dict[str, np.ndarray] = {}
        pred_lo: dict[str, np.ndarray] = {}
        pred_hi: dict[str, np.ndarray] = {}
        rms_list: list[float] = []

        for i, ed in enumerate(_iter_items(sites)):
            nm = _name(ed, i)
            _, z, fr = _get_z_block(ed)
            if z is None:
                continue
            X_obs = _z_to_features(ed, z, fr, freqs)
            if X_obs is None:
                warnings.append(f"{nm}: could not build feature vector.")
                continue
            try:
                X_in = X_obs[None, :]
                mu, sigma = ens.predict_with_uncertainty(X_in)
                # predict_intervals() returns (center, lower, upper); the
                # base estimator's raw vector is
                # [log10(rho) (n_layers) | log10(h) (n_layers - 1)], so
                # only the first n_layers entries are log10(rho).
                center, lo_arr, hi_arr = ens.predict_intervals(X_in)
                n = self.n_layers
                pred_mean[nm] = mu[0][:n]
                pred_std[nm] = sigma[0][:n]
                pred_lo[nm] = lo_arr[0][:n]
                pred_hi[nm] = hi_arr[0][:n]
                # compute forward RMS
                rms = _forward_rms(ed, mu[0][:n], freqs, self.n_layers)
                if rms is not None:
                    rms_list.append(rms)
            except Exception as exc:
                warnings.append(f"Prediction for {nm}: {exc}")

        rms_global = float(np.mean(rms_list)) if rms_list else np.nan
        n_pred = len(pred_mean)

        # ── save ensemble ─────────────────────────────────────────────────────
        if output_dir and n_pred:
            try:
                ens.save(os.path.join(output_dir, "ensemble_inverter.pkl"))
            except Exception as exc:
                warnings.append(f"Could not save ensemble: {exc}")

        # ── figures ───────────────────────────────────────────────────────────
        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}

        if pred_mean and pred_std:
            try:
                fig_unc = _plot_uncertainty_section(
                    pred_mean,
                    pred_std,
                    pred_lo,
                    pred_hi,
                    self.n_layers,
                    freqs,
                )
                if fig_unc is not None:
                    figures["uncertainty_section"] = fig_unc
                    p = self._save_figure(
                        fig_unc,
                        output_dir,
                        "ensemble_uncertainty_section",
                        warnings_list=warnings,
                    )
                    if p:
                        fig_paths["uncertainty_section"] = p
            except Exception as exc:
                warnings.append(f"Uncertainty section figure: {exc}")

        # ── uncertainty profile for first station ─────────────────────────────
        if pred_mean:
            first_station = next(iter(pred_mean))
            try:
                X_first = None
                for ii, ed in enumerate(_iter_items(sites)):
                    nm = _name(ed, ii)
                    if nm == first_station:
                        _, z, fr = _get_z_block(ed)
                        X_first = _z_to_features(ed, z, fr, freqs)
                        break
                if X_first is not None:
                    fig_prof = ens.plot_uncertainty_profile(
                        X_first[None, :],
                        sample_idx=0,
                    )
                    if fig_prof is not None:
                        f = (
                            fig_prof
                            if hasattr(fig_prof, "savefig")
                            else (
                                fig_prof.get_figure()
                                if hasattr(fig_prof, "get_figure")
                                else None
                            )
                        )
                        if f is not None:
                            figures["uncertainty_profile"] = f
                            p = self._save_figure(
                                f,
                                output_dir,
                                "ensemble_uncertainty_profile",
                                warnings_list=warnings,
                            )
                            if p:
                                fig_paths["uncertainty_profile"] = p
            except Exception as exc:
                warnings.append(f"Uncertainty profile figure: {exc}")

        # ── LLM interpretation ────────────────────────────────────────────────
        interp: str | None = None
        if self.api_key and n_pred:
            mean_std = (
                float(np.nanmean([v.mean() for v in pred_std.values()]))
                if pred_std
                else np.nan
            )
            prompt = (
                f"Ensemble inversion summary:\n"
                f"  N estimators: {self.n_estimators}, arch: {self.arch}\n"
                f"  Stations predicted: {n_pred}\n"
                f"  Global RMS: {rms_global:.3f}\n"
                f"  Mean uncertainty (std log₁₀ ρ): {mean_std:.3f}\n"
                f"  Empirical 90% coverage: {coverage:.1%}\n"
                f"  Calibrated: {self.calibrate}\n"
                f"  Warnings: {warnings[:3] if warnings else 'none'}\n\n"
                "Evaluate the ensemble inversion quality and uncertainty."
            )
            interp = self.query_llm(prompt, max_tokens=250)

        elapsed = time.time() - t0
        cov_str = f"{coverage:.1%}" if not np.isnan(coverage) else "N/A"
        rms_str = f"{rms_global:.3f}" if not np.isnan(rms_global) else "N/A"
        return AgentResult(
            status="success" if n_pred > 0 else "needs_review",
            summary=(
                f"Ensemble inversion ({self.n_estimators}× {self.arch}): "
                f"{n_pred} stations. RMS={rms_str}. "
                f"90% coverage={cov_str}. {len(figures)} figures."
            ),
            data={
                "ensemble": ens,
                "pred_mean": pred_mean,
                "pred_std": pred_std,
                "pred_lo": pred_lo,
                "pred_hi": pred_hi,
                "coverage": coverage,
                "rms_global": rms_global,
                "freqs": freqs,
                "figures": figures,
                "figure_paths": fig_paths,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


# ── private helpers ───────────────────────────────────────────────────────────


def _forward_rms(
    ed: Any, log_rho: np.ndarray, freqs: np.ndarray, n_layers: int
) -> float | None:
    """Re-run forward on predicted model and compute log-ρa RMS."""
    try:
        from ..emtools._core import _get_z_block
        from ..forward import LayeredModel, MT1DForward

        _, z, fr = _get_z_block(ed)
        if z is None:
            return None
        rhos = 10**log_rho
        ths = _default_thicknesses(n_layers, freqs)
        lm = LayeredModel(
            resistivity=rhos,
            thickness=ths[: n_layers - 1],
        )
        fwd = MT1DForward(freqs=freqs)
        resp = fwd.run(lm)
        rho_fwd = np.asarray(resp.rho_a)
        rho_xy = rho_fwd[:, 0, 1] if rho_fwd.ndim == 3 else rho_fwd
        rho_raw = getattr(ed, "rho", None)
        rho_obs = (
            rho_raw[:, 0, 1]
            if rho_raw is not None
            else (0.2 / np.where(fr == 0, np.nan, fr))
            * np.abs(z[:, 0, 1]) ** 2
        )
        per = 1.0 / np.where(fr == 0, np.nan, fr)
        per_fwd = 1.0 / np.where(freqs == 0, np.nan, freqs)
        mask = np.isfinite(per) & (rho_obs > 0)
        if mask.sum() < 2:
            return None
        interp = np.interp(
            np.log10(per[mask]),
            np.log10(per_fwd[np.isfinite(per_fwd)]),
            np.log10(np.clip(rho_xy[np.isfinite(per_fwd)], 1e-6, None)),
        )
        obs_log = np.log10(np.clip(rho_obs[mask], 1e-6, None))
        return float(np.sqrt(np.mean((obs_log - interp) ** 2)))
    except Exception:
        return None


def _plot_uncertainty_section(
    pred_mean: dict,
    pred_std: dict,
    pred_lo: dict,
    pred_hi: dict,
    n_layers: int,
    freqs: np.ndarray,
) -> Any:
    """Plot mean ± 2σ uncertainty section for all stations."""
    import matplotlib.pyplot as plt

    from ..api.section import PYCSAMT_SECTION
    from ..api.station import PYCSAMT_STATION_RENDERING

    station_names = list(pred_mean.keys())
    n_st = len(station_names)
    if n_st == 0:
        return None

    mat_mu = np.full((n_layers, n_st), np.nan)
    mat_std = np.full((n_layers, n_st), np.nan)
    for si, nm in enumerate(station_names):
        n = min(len(pred_mean[nm]), n_layers)
        mat_mu[:n, si] = pred_mean[nm][:n]
        if nm in pred_std:
            mat_std[:n, si] = pred_std[nm][:n]

    ths = _default_thicknesses(n_layers, freqs)
    depths = np.concatenate([[0], np.cumsum(ths)]) / 1000.0  # km

    section = PYCSAMT_SECTION.style_for("inversion")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ext = (-0.5, n_st - 0.5, depths[-1], depths[0])

    vv = mat_mu[np.isfinite(mat_mu)]
    vmin = float(np.percentile(vv, 5)) if vv.size else 0.0
    vmax = float(np.percentile(vv, 95)) if vv.size else 4.0

    im0 = axes[0].imshow(
        mat_mu,
        aspect="auto",
        origin="upper",
        extent=ext,
        cmap="jet_r",
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
    )
    axes[0].set_title(
        "Ensemble mean  $\\log_{10}\\rho$", fontsize=9, fontweight="bold"
    )

    sv = mat_std[np.isfinite(mat_std)]
    s_vmax = float(np.percentile(sv, 95)) if sv.size else 0.5
    im1 = axes[1].imshow(
        mat_std,
        aspect="auto",
        origin="upper",
        extent=ext,
        cmap="Oranges",
        vmin=0.0,
        vmax=s_vmax,
        interpolation="nearest",
    )
    axes[1].set_title(
        "Uncertainty  $\\sigma(\\log_{10}\\rho)$",
        fontsize=9,
        fontweight="bold",
    )

    for ax, im, lbl in [
        (axes[0], im0, "$\\log_{10}\\rho$ (Ω·m)"),
        (axes[1], im1, "$\\sigma$ (Ω·m)"),
    ]:
        PYCSAMT_STATION_RENDERING.apply(
            ax,
            np.arange(n_st, dtype=float),
            station_names,
            preset="inversion",
            xlim=(-0.5, n_st - 0.5),
        )
        ax.set_ylabel("Depth (km)", fontsize=9)
        ax.tick_params(axis="y", labelsize=8)
        section.add_colorbar(im, ax, label=lbl)

    fig.suptitle(
        "Ensemble inversion — mean and uncertainty",
        fontsize=11,
        fontweight="bold",
    )
    fig.tight_layout()
    return fig


__all__ = ["EnsembleAgent"]
