"""
pycsamt.agents.inv2d_agent
===========================

:class:`Inv2DAgent` — U-Net based 2-D MT profile inversion.

Wraps :class:`~pycsamt.ai.inversion.inv2d.EMInverter2D`:

* Assembles a **2-D pseudosection** image from the observed Sites (station ×
  frequency × impedance component) as the U-Net input.
* Generates a matching synthetic **2-D training dataset** by tiling 1-D
  forward models into profile-shaped arrays.
* Trains the U-Net and produces a **resistivity section** output:
  (n_depth × n_stations) in log₁₀ Ω·m.
* Visualises the result with
  :func:`~pycsamt.ai.plot.inversion.plot_inversion_result_2d`.

The U-Net treats the whole profile at once, so it naturally captures
lateral continuity — a key advantage over station-by-station 1-D inversion.

Requires PyTorch **or** TensorFlow.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np

from ._base import AgentResult, BaseAgent
from .ai_inversion import _z_to_features, _default_thicknesses

_SYSTEM_PROMPT = """\
You are an expert in 2-D MT inversion using deep learning (U-Net architecture).
Given a 2-D AI inversion result, write 4-5 sentences that:
1. Describe the input pseudosection geometry (stations × frequencies).
2. Interpret the dominant structural features in the resistivity section.
3. Assess lateral continuity and compare to classical smoothness-constrained results.
4. Identify artefacts or stations with poor convergence.
5. Recommend follow-up (regularisation, 3-D verification, drilling targets).
Reply in plain scientific English.
"""

_DEFAULT_FREQS_2D = np.logspace(-4, 3, 32)   # 32 frequencies — U-Net input size


class Inv2DAgent(BaseAgent):
    """2-D MT profile inversion using a U-Net convolutional architecture.

    Parameters
    ----------
    api_key, model, llm_provider : str
    n_depth : int
        Number of depth cells in the output section (default 40).
    n_freqs : int
        Number of input frequencies (default 32).
    n_components : int
        Number of impedance components in input (default 4: Re/Im × xy/yx).
    arch : str
        U-Net variant (default ``"unet"``).
    n_train_profiles : int
        Number of synthetic 2-D profiles for training (default 200).
    n_stations_per_profile : int
        Stations per synthetic profile (default 20).
    epochs : int
        Training epochs (default 30).

    Input keys
    ----------
    ``sites`` / ``path`` : Sites or str
    ``output_dir`` : str, optional
    ``period_range`` : [T_min, T_max], optional

    Output data keys
    ----------------
    ``pred_section``      ndarray (n_depth × n_stations) — log₁₀ ρ
    ``depths_km``         ndarray — depth axis (km)
    ``station_names``     list[str]
    ``rms_global``        float
    ``inverter``          EMInverter2D
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
        n_depth: int = 40,
        n_freqs: int = 32,
        n_components: int = 4,
        arch: str = "unet",
        n_train_profiles: int = 200,
        n_stations_per_profile: int = 20,
        epochs: int = 30,
    ) -> None:
        super().__init__(
            "Inv2DAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="inversion",
        )
        self.n_depth               = n_depth
        self.n_freqs               = n_freqs
        self.n_components          = n_components
        self.arch                  = arch
        self.n_train_profiles      = n_train_profiles
        self.n_stations_per_profile= n_stations_per_profile
        self.epochs                = epochs

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        # ── backend check ──────────────────────────────────────────────────────
        try:
            from ..ai.inversion.inv2d import EMInverter2D
            from ..forward.batch      import generate_dataset
            from ..backends           import get_backend_instance
            if get_backend_instance() is None:
                raise ImportError("No DL backend.")
        except ImportError as exc:
            return AgentResult.failed(
                f"Inv2DAgent requires PyTorch or TensorFlow: {exc}",
                hint="pip install torch  or  pip install tensorflow",
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

        output_dir = input_data.get("output_dir")
        import os
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)

        n_sta_obs = sum(1 for _ in _iter_items(sites))
        n_sta_use = min(n_sta_obs, self.n_stations_per_profile)
        freqs = _DEFAULT_FREQS_2D[:self.n_freqs]

        # ── build observed pseudosection ───────────────────────────────────────
        station_names: list[str] = []
        obs_feats: list[np.ndarray] = []  # each: (n_freqs, n_components)

        for i, ed in enumerate(_iter_items(sites)):
            if len(station_names) >= n_sta_use:
                break
            nm = _name(ed, i)
            _, z, fr = _get_z_block(ed)
            if z is None:
                continue
            feat = _z_to_features(z, fr, freqs)
            if feat is None:
                warnings.append(f"{nm}: skipped (bad data).")
                continue
            station_names.append(nm)
            obs_feats.append(feat)

        if len(station_names) < 3:
            return AgentResult.failed(
                "Fewer than 3 usable stations — cannot run 2-D inversion.",
                elapsed=time.time() - t0,
            )

        n_sta = len(station_names)
        # X_obs: (n_components, n_freqs, n_sta)
        X_obs = np.stack(obs_feats, axis=2)          # (n_freqs, n_components, n_sta)
        X_obs = X_obs.transpose(1, 0, 2)             # (n_components, n_freqs, n_sta)
        X_obs_4d = X_obs[None, ...]                  # (1, n_components, n_freqs, n_sta)

        # ── generate synthetic 2-D training data ───────────────────────────────
        self._log.info(
            "Generating %d synthetic 2-D profiles (%d stations × %d freqs)…",
            self.n_train_profiles, n_sta, self.n_freqs,
        )
        try:
            n_layers = max(3, self.n_depth // 8)
            n_1d = self.n_train_profiles * n_sta
            ds1d = generate_dataset(
                solver="mt1d", n_samples=n_1d, freqs=freqs,
                n_layers=n_layers, noise_level=0.03, seed=0,
                n_jobs=1, verbose=False,
            )
            # X1d: (n_1d, n_freqs, n_comp=4)  →  reshape to profiles
            X1d = ds1d.X          # (n_1d, n_freqs, 4)
            y1d = ds1d.y          # (n_1d, n_layers)

            n_samp = n_1d // n_sta
            X1d    = X1d[:n_samp * n_sta]
            y1d    = y1d[:n_samp * n_sta]

            # reshape: (n_samp, n_sta, n_freqs, 4)
            X2d_raw = X1d.reshape(n_samp, n_sta, self.n_freqs, 4)
            # → (n_samp, 4, n_freqs, n_sta)
            X2d = X2d_raw.transpose(0, 3, 2, 1)

            y2d_raw = y1d.reshape(n_samp, n_sta, n_layers)
            # → (n_samp, n_layers, n_sta)  then upsample to n_depth
            y2d = y2d_raw.transpose(0, 2, 1)  # (n_samp, n_layers, n_sta)
            if n_layers < self.n_depth:
                from scipy.ndimage import zoom
                y2d = zoom(y2d,
                           (1, self.n_depth / n_layers, 1),
                           order=1)  # (n_samp, n_depth, n_sta)

        except Exception as exc:
            return AgentResult.failed(
                f"2-D dataset assembly failed: {exc}", elapsed=time.time() - t0,
            )

        # ── train U-Net ────────────────────────────────────────────────────────
        self._log.info("Training EMInverter2D (%s) for %d epochs…",
                       self.arch, self.epochs)
        try:
            inv2d = EMInverter2D(
                n_components=self.n_components,
                n_depth=self.n_depth,
                n_stations=n_sta,
                n_freqs=self.n_freqs,
                arch=self.arch,
            )
            inv2d.fit(
                X2d, y2d,
                epochs=self.epochs,
                batch_size=max(4, min(16, n_samp // 10)),
                patience=max(5, self.epochs // 5),
                verbose=False,
            )
        except Exception as exc:
            return AgentResult.failed(
                f"U-Net training failed: {exc}", elapsed=time.time() - t0,
            )

        # ── predict ────────────────────────────────────────────────────────────
        try:
            pred_2d = inv2d.predict(X_obs_4d)[0]   # (n_depth, n_sta)
        except Exception as exc:
            return AgentResult.failed(
                f"2-D prediction failed: {exc}", elapsed=time.time() - t0,
            )

        ths    = _default_thicknesses(self.n_depth, freqs)
        depths = np.concatenate([[0.0], np.cumsum(ths)]) / 1000.0  # km

        # ── figures ───────────────────────────────────────────────────────────
        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}

        try:
            from ..ai.plot.inversion import plot_inversion_result_2d
            station_pos = np.arange(n_sta, dtype=float) * 0.5  # 0.5 km spacing
            fig_inv = plot_inversion_result_2d(
                pred_2d,
                depths=depths * 1000.0,         # expects metres
                stations=station_pos,
                station_labels=station_names,
                depth_max=float(depths[-1]) * 1000.0,
                show_misfit=False,
                show_convergence=False,
                show_rmse=False,
                suptitle="2-D AI inversion (U-Net)",
            )
            if fig_inv is not None:
                figures["inv2d_section"] = fig_inv
                p = self._save_figure(fig_inv, output_dir, "inv2d_section",
                                      warnings_list=warnings)
                if p:
                    fig_paths["inv2d_section"] = p
        except Exception as exc:
            warnings.append(f"plot_inversion_result_2d: {exc}")
            # fallback: simple imshow
            try:
                import matplotlib.pyplot as plt
                from ..api.station import PYCSAMT_STATION_RENDERING
                fig, ax = plt.subplots(figsize=(12, 5))
                vv = pred_2d[np.isfinite(pred_2d)]
                im = ax.imshow(
                    pred_2d, aspect="auto", origin="upper",
                    extent=(-0.5, n_sta - 0.5, depths[-1], depths[0]),
                    cmap="jet_r",
                    vmin=float(np.percentile(vv, 5)) if vv.size else 0,
                    vmax=float(np.percentile(vv, 95)) if vv.size else 4,
                    interpolation="bilinear",
                )
                PYCSAMT_STATION_RENDERING.apply(
                    ax, np.arange(n_sta, dtype=float), station_names,
                    preset="inversion", xlim=(-0.5, n_sta - 0.5),
                )
                ax.set_ylabel("Depth (km)", fontsize=9)
                ax.tick_params(axis="y", labelsize=8)
                self._section.add_colorbar(im, ax, label="$\\log_{10}\\rho$ (Ω·m)")
                ax.set_title("2-D AI inversion (U-Net)", fontsize=10, fontweight="bold")
                fig.tight_layout()
                figures["inv2d_section"] = fig
                p = self._save_figure(fig, output_dir, "inv2d_section",
                                      warnings_list=warnings)
                if p:
                    fig_paths["inv2d_section"] = p
            except Exception:
                pass

        # ── LLM interpretation ────────────────────────────────────────────────
        interp: str | None = None
        if self.api_key:
            rho_mean = float(np.nanmean(10**pred_2d))
            rho_std  = float(np.nanstd(10**pred_2d))
            prompt = (
                f"2-D AI inversion (U-Net) summary:\n"
                f"  Profile: {n_sta} stations × {self.n_freqs} frequencies\n"
                f"  Section: {self.n_depth} depth cells, "
                f"  max depth {depths[-1]:.1f} km\n"
                f"  Mean resistivity: {rho_mean:.0f} Ω·m ± {rho_std:.0f}\n"
                f"  Warnings: {warnings[:3] if warnings else 'none'}\n\n"
                "Interpret the 2-D resistivity section."
            )
            interp = self.query_llm(prompt, max_tokens=250)

        # ── data-space RMS ────────────────────────────────────────
        rms_global = _compute_rms_2d(
            X_obs, pred_2d, ths, freqs
        )

        elapsed = time.time() - t0
        return AgentResult(
            status="success",
            summary=(
                f"2-D AI inversion (U-Net): {n_sta} stations x "
                f"{self.n_depth} depth cells. "
                f"RMS={rms_global:.3f}. "
                f"{len(figures)} figures."
            ),
            data={
                "pred_section":  pred_2d,
                "depths_km":     depths,
                "station_names": station_names,
                "rms_global":    rms_global,
                "inverter":      inv2d,
                "figures":       figures,
                "figure_paths":  fig_paths,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


def _compute_rms_2d(
    X_obs: np.ndarray,
    pred_2d: np.ndarray,
    thicknesses: np.ndarray,
    freqs: np.ndarray,
) -> float:
    r"""
    Data-space RMS for the 2-D section.

    Uses the Bostick approximation to convert the
    predicted log10-resistivity section back to
    apparent resistivity and phase, then compares to
    the observed log10(rho_a) stored in X_obs.

    Parameters
    ----------
    X_obs : ndarray, shape (n_components, n_freqs, n_sta)
        Observed features; component 0 assumed to be
        log10(rho_a).
    pred_2d : ndarray, shape (n_depth, n_sta)
        Predicted log10(rho) section.
    thicknesses : ndarray, shape (n_depth - 1,) or (n_depth,)
        Layer thicknesses in metres.
    freqs : ndarray, shape (n_freqs,)
        Frequency array in Hz.

    Returns
    -------
    float
        Global normalised RMS in log-resistivity space.
    """
    try:
        n_comp, n_freqs, n_sta = X_obs.shape
        n_depth = pred_2d.shape[0]

        depths_m = np.concatenate(
            [[0.0], np.cumsum(thicknesses[:n_depth])]
        )
        periods = 1.0 / np.maximum(freqs, 1e-9)

        obs_log_rho = X_obs[0]  # (n_freqs, n_sta)

        # Bostick: rho_Bostick(T) ~ rho_a(T) * (phi/45 - 1)
        # Here we simply read off the predicted profile at the
        # Bostick depth d_B = 503 * sqrt(rho_a / f)
        rms_vals: list[float] = []
        for s in range(n_sta):
            profile = pred_2d[:, s]  # log10(rho), n_depth cells
            pred_log_rho_a = np.zeros(n_freqs)
            for fi, T in enumerate(periods):
                rho_a_obs = 10.0 ** obs_log_rho[fi, s]
                d_b = 503.0 * np.sqrt(
                    max(rho_a_obs, 1.0) * T
                )
                idx = int(
                    np.searchsorted(depths_m, d_b)
                )
                idx = min(idx, n_depth - 1)
                pred_log_rho_a[fi] = profile[idx]
            diff = pred_log_rho_a - obs_log_rho[:, s]
            finite = np.isfinite(diff)
            if finite.any():
                rms_vals.append(
                    float(np.sqrt(np.mean(diff[finite] ** 2)))
                )

        return float(np.mean(rms_vals)) if rms_vals else np.nan

    except Exception:
        return np.nan


__all__ = ["Inv2DAgent"]
