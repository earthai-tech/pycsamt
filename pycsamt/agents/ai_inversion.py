"""
pycsamt.agents.ai_inversion
============================

:class:`AIInversionAgent` — End-to-end AI 1-D MT inversion.

Workflow
--------
1. **Generate synthetic training data** — :func:`~pycsamt.forward.batch.generate_dataset`
   builds a ``ForwardDataset`` of (response, model) pairs covering a
   realistic resistivity range.
2. **Train** — :class:`~pycsamt.ai.inversion.inv1d.EMInverter1D` is fitted
   on the synthetic data using the selected network architecture.
3. **Predict** — the trained inverter predicts a layered resistivity model
   for every station in the observed dataset.
4. **Evaluate** — per-station RMS between observed and re-computed forward
   response validates the prediction quality.
5. **Visualise** — convergence curve + predicted model sections.

Requires PyTorch **or** TensorFlow (lazy import — no hard dependency at
import time).  When neither is available the agent returns a clear error
message with an installation hint.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any

import numpy as np

from ._base import AgentResult, BaseAgent

_SYSTEM_PROMPT = """\
You are an expert in AI-based MT inversion and deep learning for geophysics.
Given an AI inversion result, write 4–5 sentences that:
1. Describe the neural network architecture used and training convergence.
2. State the prediction quality (RMS, layer count, depth range).
3. Identify stations where the AI prediction is most / least reliable.
4. Compare AI results with classical Bostick depth estimates if available.
5. Recommend next steps (fine-tuning, ensemble, switch to 2-D inversion).
Reply in plain scientific English.
"""

# ── default training config ───────────────────────────────────────────────────
_DEFAULT_FREQS  = np.logspace(-4, 3, 40)  # 40 frequencies 10⁻⁴ – 10³ Hz
_DEFAULT_N_SAMP = 2_000                   # fast default; 10 000 for production
_DEFAULT_EPOCHS = 30                      # fast default; 100+ for production


class AIInversionAgent(BaseAgent):
    """Train an AI inverter on synthetic data then predict on observed sites.

    Parameters
    ----------
    api_key, model, llm_provider : str
    arch : {"resnet", "cnn1d", "fcn"}
        Neural network architecture.
    n_layers : int
        Number of model layers the inverter will predict (default 5).
    n_train_samples : int
        Number of synthetic training samples (default 2 000).
    epochs : int
        Training epochs (default 30).  Increase for better models.
    freqs : array-like or None
        Frequencies used for both training synthesis and observed data
        interpolation.  Default: 40 log-spaced 10⁻⁴–10³ Hz.
    pretrained : str or None
        Path to a pre-trained model checkpoint.  When set, skips training.

    Input keys
    ----------
    ``sites`` / ``path`` : Sites or str — observed data
    ``output_dir`` : str, optional
    ``arch``, ``epochs``, ``n_train_samples`` : optional overrides

    Output data keys
    ----------------
    ``inverter``          :class:`~pycsamt.ai.inversion.inv1d.EMInverter1D`
    ``predictions``       dict {station: ndarray of log₁₀ ρ values}
    ``best_model``        dict with "resistivity" and "thickness" for first station
    ``rms_per_station``   dict {station: float}
    ``rms_global``        float
    ``train_history``     dict (loss curves)
    ``figures``           dict
    ``figure_paths``      dict

    Examples
    --------
    >>> agent  = AIInversionAgent(arch="resnet", n_layers=5, epochs=30)
    >>> result = agent.execute({
    ...     "path":       "/data/L22PLT",
    ...     "output_dir": "/out/ai_inv",
    ... })
    >>> result["rms_global"]
    0.24
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    @classmethod
    def from_pretrained(
        cls,
        model_name: str,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        cache_dir: str | None = None,
        force_download: bool = False,
    ) -> "AIInversionAgent":
        """Return an :class:`AIInversionAgent` pre-loaded with a zoo checkpoint.

        Parameters
        ----------
        model_name : str
            Registry name — see :func:`~pycsamt.ai._zoo.list_pretrained`.
        cache_dir : str or None
            Override default cache ``~/.pycsamt/model_zoo/``.
        force_download : bool
            Re-download even if cached.

        Examples
        --------
        >>> agent = AIInversionAgent.from_pretrained("mt1d-resnet-5layer-v1")
        >>> result = agent.execute({"path": "/data/L22PLT"})
        """
        from ..ai._zoo import download_checkpoint, get_pretrained_info

        info      = get_pretrained_info(model_name)
        ckpt_path = download_checkpoint(
            model_name, cache_dir=cache_dir, force=force_download, verbose=True,
        )
        return cls(
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            arch=info.get("arch", "resnet"),
            n_layers=int(info.get("n_layers", 5)),
            pretrained=str(ckpt_path),
        )

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        arch: str = "resnet",
        n_layers: int = 5,
        n_train_samples: int = _DEFAULT_N_SAMP,
        epochs: int = _DEFAULT_EPOCHS,
        freqs: Any = None,
        pretrained: str | None = None,
    ) -> None:
        super().__init__(
            "AIInversionAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="inversion",
        )
        self.arch             = arch
        self.n_layers         = n_layers
        self.n_train_samples  = n_train_samples
        self.epochs           = epochs
        self._freqs_cfg       = freqs
        self.pretrained       = pretrained

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        # ── check AI backend ──────────────────────────────────────────────────
        try:
            from ..ai.inversion.inv1d import EMInverter1D
            from ..forward.batch import generate_dataset, ForwardDataset
            from ..backends import get_backend_instance
            if get_backend_instance() is None:
                raise ImportError("No DL backend available (torch or tensorflow).")
        except ImportError as exc:
            return AgentResult.failed(
                f"AI inversion requires PyTorch or TensorFlow: {exc}",
                hint=(
                    "Install PyTorch:     pip install torch\n"
                    "Install TensorFlow:  pip install tensorflow"
                ),
                elapsed=time.time() - t0,
            )

        from ..emtools._core import ensure_sites, _iter_items, _name, _get_z_block

        sites_raw = input_data.get("sites") or input_data.get("path")
        if sites_raw is None:
            return AgentResult.failed("No 'sites' or 'path'.", elapsed=time.time()-t0)
        try:
            sites = ensure_sites(sites_raw, verbose=0)
        except Exception as exc:
            return AgentResult.failed(str(exc), elapsed=time.time()-t0)

        arch             = str(input_data.get("arch",             self.arch))
        epochs           = int(input_data.get("epochs",           self.epochs))
        n_train          = int(input_data.get("n_train_samples",  self.n_train_samples))
        output_dir       = input_data.get("output_dir")
        freqs            = np.asarray(input_data.get("freqs") or self._freqs_cfg
                                      or _DEFAULT_FREQS, dtype=float)

        # ── build inverter ────────────────────────────────────────────────────
        inverter = EMInverter1D(
            arch=arch,
            n_layers=self.n_layers,
            solver="mt1d",
            include_phase=True,
        )

        # ── load pre-trained or train ─────────────────────────────────────────
        train_history: dict = {}
        if self.pretrained and Path(self.pretrained).exists():
            try:
                inverter = EMInverter1D.load(self.pretrained)
                self._log.info("Loaded pre-trained model from %s", self.pretrained)
            except Exception as exc:
                warnings.append(f"Could not load pre-trained model: {exc}. Training fresh.")
                self.pretrained = None

        if not self.pretrained:
            # generate synthetic training set
            self._log.info(
                "Generating %d synthetic training samples (freqs=%d) …",
                n_train, len(freqs),
            )
            try:
                dataset = generate_dataset(
                    solver="mt1d",
                    n_samples=n_train,
                    freqs=freqs,
                    n_layers=self.n_layers,
                    noise_level=0.03,
                    seed=42,
                    n_jobs=1,
                    verbose=False,
                )
                self._log.info("Training %s for %d epochs …", arch, epochs)
                inverter.fit(
                    dataset.X, dataset.y,
                    epochs=epochs,
                    batch_size=min(256, n_train // 4),
                    patience=max(5, epochs // 5),
                    verbose=False,
                )
                # collect training history
                if hasattr(inverter, "history_"):
                    train_history = dict(inverter.history_)
            except Exception as exc:
                return AgentResult.failed(
                    f"AI training failed: {exc}",
                    hint="Try reducing n_train_samples or epochs.",
                    elapsed=time.time() - t0,
                )

        # ── predict on observed stations ──────────────────────────────────────
        predictions: dict[str, np.ndarray] = {}
        rms_per: dict[str, float] = {}
        best_model: dict = {}

        for i, ed in enumerate(_iter_items(sites)):
            nm = _name(ed, i)
            _, z, fr = _get_z_block(ed)
            if z is None or fr is None:
                continue
            try:
                # interpolate Z onto training freqs
                X_obs = _z_to_features(z, fr, freqs)
                if X_obs is None:
                    warnings.append(f"{nm}: could not build feature vector.")
                    continue
                log_rho_pred = inverter.predict(X_obs[None, :])[0]  # (n_layers,)
                predictions[nm] = log_rho_pred

                # compute forward RMS for this prediction
                from ..forward import MT1DForward, LayeredModel
                rhos = 10 ** log_rho_pred
                ths  = _default_thicknesses(self.n_layers, freqs)
                try:
                    lm   = LayeredModel(
                        resistivity=rhos,
                        thickness=ths[:self.n_layers - 1],
                    )
                    fwd  = MT1DForward(freqs=freqs)
                    resp = fwd.run(lm)
                    rho_fwd = np.asarray(resp.rho_a)
                    rho_fwd_xy = rho_fwd[:, 0, 1] if rho_fwd.ndim == 3 else rho_fwd

                    per_obs = 1.0 / np.where(fr == 0, np.nan, fr)
                    per_fwd = 1.0 / np.where(freqs == 0, np.nan, freqs)
                    rho_obs_raw = getattr(ed, "rho", None)
                    rho_obs = (
                        rho_obs_raw[:, 0, 1] if rho_obs_raw is not None
                        else (0.2 / np.where(fr==0,np.nan,fr)) * np.abs(z[:,0,1])**2
                    )
                    mask = np.isfinite(per_obs) & (rho_obs > 0)
                    if mask.sum() >= 2:
                        interp = np.interp(
                            np.log10(per_obs[mask]),
                            np.log10(per_fwd[np.isfinite(per_fwd)]),
                            np.log10(np.clip(rho_fwd_xy[np.isfinite(per_fwd)], 1e-6, None)),
                        )
                        obs_log = np.log10(np.clip(rho_obs[mask], 1e-6, None))
                        rms_per[nm] = float(np.sqrt(np.mean((obs_log - interp)**2)))
                except Exception:
                    pass

                if not best_model:
                    best_model = {
                        "station":     nm,
                        "resistivity": list(rhos),
                        "thickness":   list(ths[:self.n_layers - 1]),
                        "log_rho":     list(log_rho_pred),
                    }
            except Exception as exc:
                warnings.append(f"Prediction failed for {nm}: {exc}")

        rms_global = float(np.nanmean(list(rms_per.values()))) if rms_per else np.nan
        n_pred = len(predictions)

        # ── save model ────────────────────────────────────────────────────────
        if output_dir:
            import os
            os.makedirs(output_dir, exist_ok=True)
            try:
                model_path = f"{output_dir}/ai_inverter_{arch}.pkl"
                inverter.save(model_path)
            except Exception as exc:
                warnings.append(f"Could not save inverter: {exc}")

        # ── figures ───────────────────────────────────────────────────────────
        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}

        # training convergence
        if train_history:
            try:
                from ..ai.plot.convergence import plot_convergence
                fig_cv = plot_convergence(train_history)
                if fig_cv is not None:
                    f = fig_cv if hasattr(fig_cv,"savefig") else (
                        fig_cv.get_figure() if hasattr(fig_cv,"get_figure") else None)
                    if f is not None:
                        figures["convergence"] = f
                        p = self._save_figure(f, output_dir, "ai_convergence",
                                              warnings_list=warnings)
                        if p:
                            fig_paths["convergence"] = p
            except Exception as exc:
                warnings.append(f"Convergence plot: {exc}")

        # predicted model section (station × layer colour-coded by log ρ)
        if predictions:
            try:
                fig_sec = _plot_ai_section(predictions, self.n_layers, freqs)
                if fig_sec is not None:
                    figures["ai_section"] = fig_sec
                    p = self._save_figure(fig_sec, output_dir, "ai_model_section",
                                          warnings_list=warnings)
                    if p:
                        fig_paths["ai_section"] = p
            except Exception as exc:
                warnings.append(f"AI section plot: {exc}")

        # ── LLM interpretation ────────────────────────────────────────────────
        interp: str | None = None
        if self.api_key and n_pred:
            n_hi = sum(1 for v in rms_per.values() if v > 0.5)
            prompt = (
                f"AI inversion summary:\n"
                f"  Architecture: {arch}, layers: {self.n_layers}\n"
                f"  Stations predicted: {n_pred}\n"
                f"  Global RMS: {rms_global:.3f} log₁₀(Ω·m)\n"
                f"  Stations with high RMS (> 0.5): {n_hi}\n"
                f"  Best model (station {best_model.get('station', '?')}): "
                f"  ρ = {[f'{r:.0f}' for r in best_model.get('resistivity', [])]} Ω·m\n\n"
                "Evaluate this AI inversion and recommend next steps."
            )
            interp = self.query_llm(prompt, max_tokens=250)

        elapsed = time.time() - t0
        rms_str = f"RMS {rms_global:.3f}" if not np.isnan(rms_global) else "RMS N/A"
        return AgentResult(
            status="success" if n_pred > 0 else "needs_review",
            summary=(
                f"AI inversion ({arch}, {self.n_layers} layers): "
                f"{n_pred} stations predicted. {rms_str}. "
                f"{len(figures)} figure(s)."
            ),
            data={
                "inverter":         inverter,
                "predictions":      predictions,
                "best_model":       best_model,
                "rms_per_station":  rms_per,
                "rms_global":       rms_global,
                "train_history":    train_history,
                "figures":          figures,
                "figure_paths":     fig_paths,
                "freqs":            freqs,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


# ── helpers ───────────────────────────────────────────────────────────────────

def _z_to_features(
    z: np.ndarray,
    fr: np.ndarray,
    freqs_target: np.ndarray,
) -> np.ndarray | None:
    """Interpolate Z to target freqs and return a flat feature vector."""
    try:
        per     = 1.0 / np.where(fr == 0, np.nan, fr)
        per_t   = 1.0 / np.where(freqs_target == 0, np.nan, freqs_target)
        n_t     = len(freqs_target)
        feats   = []
        for ri, ci in [(0,1), (1,0)]:
            zxy = np.abs(z[:, ri, ci])
            pha = np.degrees(np.angle(z[:, ri, ci]))
            mask = np.isfinite(per) & np.isfinite(zxy)
            if mask.sum() < 2:
                return None
            log_zxy = np.log10(np.clip(zxy[mask], 1e-30, None))
            interp_z = np.interp(
                np.log10(per_t[np.isfinite(per_t)]),
                np.log10(per[mask]),
                log_zxy,
            )
            interp_p = np.interp(
                np.log10(per_t[np.isfinite(per_t)]),
                np.log10(per[mask]),
                pha[mask],
            )
            feats.append(interp_z[:n_t])
            feats.append(interp_p[:n_t])
        arr = np.column_stack(feats)   # (n_t, 4)
        return arr.astype(np.float32)
    except Exception:
        return None


def _default_thicknesses(n_layers: int, freqs: np.ndarray) -> np.ndarray:
    """Return log-spaced thicknesses spanning Bostick depth range."""
    rho_ref = 100.0   # Ω·m reference
    per_max = float(1.0 / freqs.min())
    d_max   = np.sqrt(rho_ref * per_max / (4 * np.pi * 1e-7 * 2 * np.pi))
    ths = np.logspace(np.log10(max(d_max / 100, 50)), np.log10(d_max), n_layers - 1)
    return ths.astype(float)


def _plot_ai_section(
    predictions: dict[str, np.ndarray],
    n_layers: int,
    freqs: np.ndarray,
) -> Any:
    """Plot a station × layer colour-coded predicted resistivity section."""
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    from ..api.section import PYCSAMT_SECTION

    station_names = list(predictions.keys())
    n_st = len(station_names)
    if n_st == 0:
        return None

    mat = np.full((n_layers, n_st), np.nan)
    for si, nm in enumerate(station_names):
        log_rho = predictions[nm]
        n = min(len(log_rho), n_layers)
        mat[:n, si] = log_rho[:n]

    ths = _default_thicknesses(n_layers, freqs)
    depths = np.concatenate([[0], np.cumsum(ths)]) / 1000.0   # km

    section = PYCSAMT_SECTION.style_for("inversion")
    fig_w, fig_h = section.figsize_for(n_stations=n_st, n_y=n_layers)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    im = ax.imshow(
        mat,
        aspect="auto",
        origin="upper",
        extent=(-0.5, n_st - 0.5, depths[-1], depths[0]),
        cmap="jet_r",
        vmin=np.nanpercentile(mat, 5),
        vmax=np.nanpercentile(mat, 95),
        interpolation="nearest",
    )
    from ..api.station import PYCSAMT_STATION_RENDERING
    PYCSAMT_STATION_RENDERING.apply(
        ax, np.arange(n_st, dtype=float), station_names,
        preset="inversion", xlim=(-0.5, n_st - 0.5),
    )
    ax.set_ylabel("Depth (km)", fontsize=9)
    ax.tick_params(axis="y", labelsize=8)
    section.add_colorbar(im, ax, label="$\\log_{10}\\rho$ (Ω·m)")
    ax.set_title("AI-predicted resistivity section", fontsize=10, fontweight="bold")
    fig.tight_layout()
    return fig


__all__ = ["AIInversionAgent"]
