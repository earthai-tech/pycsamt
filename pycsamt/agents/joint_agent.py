"""
pycsamt.agents.joint_agent
===========================

:class:`JointInversionAgent` — Multi-modal joint MT inversion via DRCNN.

Wraps :class:`~pycsamt.ai.inversion.joint.JointInverter`:

* Fuses a **primary** MT dataset with a **secondary** modality (TEM, CSAMT,
  gravity proxy, or a second MT profile at a different frequency range).
* Both modalities are observed at the same stations.  When no secondary
  dataset is supplied the agent synthesises a complementary low-frequency
  response from the same :class:`~pycsamt.forward.synthetic.LayeredModel`
  to demonstrate the joint-inversion pipeline.
* Produces a **joint resistivity section** (station × depth) that benefits
  from the complementary depth sensitivities of the two modalities.

Architecture
------------
The :class:`~pycsamt.ai.nets.drcnn.DRCNNNet` (Dense-Residual CNN) is used
as the shared feature extractor; each modality has its own encoding branch
before the fused prediction head.

Requires PyTorch **or** TensorFlow.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np

from ._base import AgentResult, BaseAgent
from .ai_inversion import _z_to_features, _default_thicknesses

_SYSTEM_PROMPT = """\
You are an expert in multi-modal geophysical joint inversion using deep learning.
Given a joint MT inversion result, write 4-5 sentences that:
1. Describe the two modalities fused and their complementary depth sensitivities.
2. Assess the joint prediction quality (RMS, depth range, station count).
3. Compare the joint result to a single-modality approach where possible.
4. Identify where the secondary modality most improved the primary inversion.
5. Recommend validation (borehole, gravity, seismic) and next modelling steps.
Reply in plain scientific English.
"""

_DEFAULT_FREQS_MT  = np.logspace(-4, 3, 40)   # primary:  1e-4 – 1e3 Hz  (40 pts)
_DEFAULT_FREQS_SEC = np.logspace(-4, 1, 20)    # secondary: 1e-4 – 10 Hz (20 pts)
_N_COMP_MT  = 4    # Re/Im of Zxy + Zyx
_N_COMP_SEC = 2    # |Z| + phase for one component


class JointInversionAgent(BaseAgent):
    """Multi-modal MT joint inversion using DRCNN.

    Parameters
    ----------
    api_key, model, llm_provider : str
    modalities : list[str]
        Names of the two modalities, e.g. ``["mt", "tem"]``.
        The first entry is the primary modality (loaded from ``sites``/``path``);
        the second is loaded from ``secondary_path`` or synthesised when absent.
    n_layers : int
        Number of depth layers in the output model (default 5).
    n_freqs_primary : int
        Frequencies for the primary MT response features (default 40).
    n_freqs_secondary : int
        Frequencies for the secondary modality features (default 20).
    n_train_samples : int
        Synthetic training samples shared across both modalities (default 2000).
    epochs : int
        Training epochs (default 30).
    growth_rate : int
        DRCNN dense-block growth rate (default 32).

    Input keys
    ----------
    ``sites`` / ``path`` : Sites or str — primary MT dataset
    ``secondary_path`` : str, optional — secondary modality EDI/TEM path
    ``output_dir`` : str, optional
    ``period_range`` : [T_min, T_max], optional

    Output data keys
    ----------------
    ``inverter``          JointInverter
    ``predictions``       dict {station: ndarray}  — mean log₁₀ ρ per layer
    ``rms_per_station``   dict {station: float}
    ``rms_global``        float
    ``modalities``        list[str]
    ``figures``           dict
    ``figure_paths``      dict

    Examples
    --------
    >>> agent = JointInversionAgent(modalities=["mt", "tem"], n_layers=5, epochs=20)
    >>> result = agent.execute({"path": "/data/L22PLT"})
    >>> result["rms_global"]
    0.31
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        modalities: list[str] | None = None,
        n_layers: int = 5,
        n_freqs_primary: int = 40,
        n_freqs_secondary: int = 20,
        n_train_samples: int = 2_000,
        epochs: int = 30,
        growth_rate: int = 32,
    ) -> None:
        super().__init__(
            "JointInversionAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="inversion",
        )
        self.modalities        = modalities or ["mt", "tem"]
        self.n_layers          = n_layers
        self.n_freqs_primary   = n_freqs_primary
        self.n_freqs_secondary = n_freqs_secondary
        self.n_train_samples   = n_train_samples
        self.epochs            = epochs
        self.growth_rate       = growth_rate

    # ── public ────────────────────────────────────────────────────────────────

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        # ── backend check ─────────────────────────────────────────────────────
        try:
            from ..ai.inversion.joint import JointInverter
            from ..forward.batch      import generate_dataset
            from ..backends           import get_backend_instance
            if get_backend_instance() is None:
                raise ImportError("No DL backend.")
        except ImportError as exc:
            return AgentResult.failed(
                f"JointInversionAgent requires PyTorch or TensorFlow: {exc}",
                hint="pip install torch  or  pip install tensorflow",
                elapsed=time.time() - t0,
            )

        from ..emtools._core import ensure_sites, _iter_items, _name, _get_z_block

        # ── load primary MT sites ─────────────────────────────────────────────
        sites_raw = input_data.get("sites") or input_data.get("path")
        if sites_raw is None:
            return AgentResult.failed("No 'sites' or 'path'.", elapsed=time.time() - t0)
        try:
            sites = ensure_sites(sites_raw, verbose=0)
        except Exception as exc:
            return AgentResult.failed(str(exc), elapsed=time.time() - t0)

        output_dir = input_data.get("output_dir")
        import os
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)

        freqs_mt  = _DEFAULT_FREQS_MT[:self.n_freqs_primary]
        freqs_sec = _DEFAULT_FREQS_SEC[:self.n_freqs_secondary]

        n_feat_mt  = self.n_freqs_primary  * _N_COMP_MT
        n_feat_sec = self.n_freqs_secondary * _N_COMP_SEC

        # ── load optional secondary sites ─────────────────────────────────────
        sec_path = input_data.get("secondary_path")
        sec_sites = None
        if sec_path:
            try:
                sec_sites = ensure_sites(sec_path, verbose=0)
            except Exception as exc:
                warnings.append(f"Could not load secondary sites: {exc}. "
                                "Synthesising secondary responses.")

        # ── collect observed primary features ─────────────────────────────────
        station_names: list[str] = []
        X_mt_obs: list[np.ndarray] = []

        for i, ed in enumerate(_iter_items(sites)):
            nm = _name(ed, i)
            _, z, fr = _get_z_block(ed)
            if z is None:
                continue
            feat = _z_to_features(z, fr, freqs_mt)
            if feat is None:
                warnings.append(f"{nm}: skipped (bad MT data).")
                continue
            station_names.append(nm)
            X_mt_obs.append(feat.reshape(-1))   # (n_feat_mt,)

        if len(station_names) < 2:
            return AgentResult.failed(
                "Need ≥ 2 usable stations for joint inversion.",
                elapsed=time.time() - t0,
            )

        n_sta = len(station_names)
        X_mt_obs_arr = np.stack(X_mt_obs, axis=0).astype(np.float32)  # (n_sta, n_feat_mt)

        # ── collect or synthesise secondary features ──────────────────────────
        X_sec_obs_arr = _collect_secondary_features(
            station_names, sites, sec_sites, freqs_mt, freqs_sec,
            n_feat_sec, warnings,
        )

        # ── generate synthetic training data ──────────────────────────────────
        self._log.info(
            "Generating %d synthetic samples for joint training…",
            self.n_train_samples,
        )
        try:
            ds_mt = generate_dataset(
                solver="mt1d",
                n_samples=self.n_train_samples,
                freqs=freqs_mt,
                n_layers=self.n_layers,
                noise_level=0.03,
                seed=42,
                n_jobs=1,
                verbose=False,
            )
            X_mt_train  = ds_mt.X.reshape(len(ds_mt.X), -1).astype(np.float32)
            y_train     = ds_mt.y.astype(np.float32)

            # secondary: same models, low-frequency range, 2 components
            ds_sec = generate_dataset(
                solver="mt1d",
                n_samples=self.n_train_samples,
                freqs=freqs_sec,
                n_layers=self.n_layers,
                noise_level=0.05,
                seed=123,
                n_jobs=1,
                verbose=False,
            )
            # keep only first 2 components (log|Zxy|, phase Zxy)
            X_sec_raw   = ds_sec.X   # (n_samp, n_freqs_sec, 4)
            X_sec_train = X_sec_raw[:, :, :_N_COMP_SEC].reshape(
                len(X_sec_raw), -1
            ).astype(np.float32)

            # align feature sizes if mismatch
            X_mt_train  = _pad_or_trim(X_mt_train,  n_feat_mt)
            X_sec_train = _pad_or_trim(X_sec_train, n_feat_sec)

        except Exception as exc:
            return AgentResult.failed(
                f"Synthetic dataset generation failed: {exc}",
                elapsed=time.time() - t0,
            )

        # ── train JointInverter ───────────────────────────────────────────────
        self._log.info(
            "Training JointInverter (DRCNN) — modalities %s — %d epochs…",
            self.modalities, self.epochs,
        )
        try:
            inverter = JointInverter(
                n_features_list=(n_feat_mt, n_feat_sec),
                n_layers=self.n_layers,
                growth_rate=self.growth_rate,
            )
            inverter.fit(
                [X_mt_train, X_sec_train],
                y_train,
                epochs=self.epochs,
                batch_size=min(256, self.n_train_samples // 4),
                patience=max(5, self.epochs // 5),
                verbose=False,
            )
        except Exception as exc:
            return AgentResult.failed(
                f"JointInverter training failed: {exc}",
                elapsed=time.time() - t0,
            )

        # ── predict on observed stations ──────────────────────────────────────
        predictions: dict[str, np.ndarray] = {}
        rms_per: dict[str, float] = {}

        X_mt_obs_arr  = _pad_or_trim(X_mt_obs_arr,  n_feat_mt)
        X_sec_obs_arr = _pad_or_trim(X_sec_obs_arr, n_feat_sec)

        try:
            y_pred_all = inverter.predict([X_mt_obs_arr, X_sec_obs_arr])  # (n_sta, n_layers)
        except Exception as exc:
            return AgentResult.failed(
                f"Joint prediction failed: {exc}", elapsed=time.time() - t0,
            )

        for si, nm in enumerate(station_names):
            log_rho = y_pred_all[si]
            predictions[nm] = log_rho
            rms = _forward_rms_joint(sites, nm, si, log_rho, freqs_mt, self.n_layers)
            if rms is not None:
                rms_per[nm] = rms

        rms_global = float(np.nanmean(list(rms_per.values()))) if rms_per else np.nan
        n_pred = len(predictions)

        # ── figures ───────────────────────────────────────────────────────────
        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}

        try:
            fig_sec = _plot_joint_section(
                predictions, self.n_layers, freqs_mt, station_names,
                self.modalities,
            )
            if fig_sec is not None:
                figures["joint_section"] = fig_sec
                p = self._save_figure(fig_sec, output_dir, "joint_inversion_section",
                                      warnings_list=warnings)
                if p:
                    fig_paths["joint_section"] = p
        except Exception as exc:
            warnings.append(f"Joint section figure: {exc}")

        # ── LLM interpretation ────────────────────────────────────────────────
        interp: str | None = None
        if self.api_key and n_pred:
            rms_str = f"{rms_global:.3f}" if not np.isnan(rms_global) else "N/A"
            prompt = (
                f"Joint inversion summary:\n"
                f"  Modalities: {' + '.join(self.modalities)}\n"
                f"  Stations predicted: {n_pred}\n"
                f"  Layers: {self.n_layers}, max depth: "
                f"{_default_thicknesses(self.n_layers, freqs_mt).sum()/1000:.1f} km\n"
                f"  Global RMS: {rms_str} log₁₀(Ω·m)\n"
                f"  Secondary synthesised: {sec_sites is None}\n"
                f"  Warnings: {warnings[:3] if warnings else 'none'}\n\n"
                "Interpret the joint inversion result and compare to single-modality approaches."
            )
            interp = self.query_llm(prompt, max_tokens=250)

        elapsed = time.time() - t0
        rms_disp = f"RMS {rms_global:.3f}" if not np.isnan(rms_global) else "RMS N/A"
        return AgentResult(
            status="success" if n_pred > 0 else "needs_review",
            summary=(
                f"Joint inversion ({' + '.join(self.modalities)}): "
                f"{n_pred} stations, {self.n_layers} layers. "
                f"{rms_disp}. {len(figures)} figure(s)."
            ),
            data={
                "inverter":        inverter,
                "predictions":     predictions,
                "rms_per_station": rms_per,
                "rms_global":      rms_global,
                "modalities":      self.modalities,
                "freqs_mt":        freqs_mt,
                "freqs_secondary": freqs_sec,
                "figures":         figures,
                "figure_paths":    fig_paths,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


# ── private helpers ───────────────────────────────────────────────────────────

def _collect_secondary_features(
    station_names: list[str],
    primary_sites: Any,
    secondary_sites: Any | None,
    freqs_mt: np.ndarray,
    freqs_sec: np.ndarray,
    n_feat_sec: int,
    warnings: list[str],
) -> np.ndarray:
    """Build secondary feature matrix (n_sta, n_feat_sec).

    When *secondary_sites* are available they are matched by station index.
    Otherwise a low-frequency sub-set of the primary impedance is used to
    simulate a complementary modality (TEM-like amplitude + phase).
    """
    from ..emtools._core import _iter_items, _name, _get_z_block

    n_sta = len(station_names)
    X_sec = np.zeros((n_sta, n_feat_sec), dtype=np.float32)

    if secondary_sites is not None:
        sec_iter = list(_iter_items(secondary_sites))
        for si in range(min(n_sta, len(sec_iter))):
            ed = sec_iter[si]
            _, z, fr = _get_z_block(ed)
            if z is None:
                continue
            # log|Zxy| + phase Zxy at secondary frequencies
            feat = _extract_sec_features(z, fr, freqs_sec)
            if feat is not None:
                X_sec[si, :len(feat)] = feat[:n_feat_sec]
        return X_sec

    # fallback: derive from primary sites at low-frequency subset
    primary_iter = list(_iter_items(primary_sites))
    for si, ed in enumerate(primary_iter):
        if si >= n_sta:
            break
        _, z, fr = _get_z_block(ed)
        if z is None:
            continue
        feat = _extract_sec_features(z, fr, freqs_sec)
        if feat is not None:
            X_sec[si, :len(feat)] = feat[:n_feat_sec]

    warnings.append(
        "No secondary dataset provided — low-frequency MT sub-band used as "
        "proxy for secondary modality."
    )
    return X_sec


def _extract_sec_features(
    z: np.ndarray,
    fr: np.ndarray,
    freqs_target: np.ndarray,
) -> np.ndarray | None:
    """Return (n_freqs_target * 2,) — log|Zxy| and phase at target freqs."""
    try:
        per     = 1.0 / np.where(fr == 0, np.nan, fr)
        per_t   = 1.0 / np.where(freqs_target == 0, np.nan, freqs_target)
        n_t     = len(freqs_target)
        zxy     = np.abs(z[:, 0, 1])
        pha     = np.degrees(np.angle(z[:, 0, 1]))
        mask    = np.isfinite(per) & np.isfinite(zxy) & (zxy > 0)
        if mask.sum() < 2:
            return None
        log_z = np.log10(np.clip(zxy[mask], 1e-30, None))
        log_z_interp = np.interp(
            np.log10(per_t[np.isfinite(per_t)]),
            np.log10(per[mask]),
            log_z,
        )
        pha_interp = np.interp(
            np.log10(per_t[np.isfinite(per_t)]),
            np.log10(per[mask]),
            pha[mask],
        )
        arr = np.concatenate([log_z_interp[:n_t], pha_interp[:n_t]])
        return arr.astype(np.float32)
    except Exception:
        return None


def _pad_or_trim(X: np.ndarray, target_cols: int) -> np.ndarray:
    """Ensure X has exactly *target_cols* columns."""
    n, c = X.shape
    if c == target_cols:
        return X
    if c > target_cols:
        return X[:, :target_cols]
    # pad with zeros
    pad = np.zeros((n, target_cols - c), dtype=X.dtype)
    return np.concatenate([X, pad], axis=1)


def _forward_rms_joint(
    sites: Any,
    station_name: str,
    station_idx: int,
    log_rho: np.ndarray,
    freqs: np.ndarray,
    n_layers: int,
) -> float | None:
    """Compute forward-response RMS for one predicted station."""
    try:
        from ..forward import MT1DForward, LayeredModel
        from ..emtools._core import _iter_items, _name, _get_z_block

        for i, ed in enumerate(_iter_items(sites)):
            if _name(ed, i) != station_name and i != station_idx:
                continue
            _, z, fr = _get_z_block(ed)
            if z is None:
                return None
            rhos = 10 ** log_rho
            ths  = _default_thicknesses(n_layers, freqs)
            lm   = LayeredModel(resistivity=rhos, thickness=ths[:n_layers - 1])
            fwd  = MT1DForward(freqs=freqs)
            resp = fwd.run(lm)
            rho_fwd = np.asarray(resp.rho_a)
            rho_xy  = rho_fwd[:, 0, 1] if rho_fwd.ndim == 3 else rho_fwd
            rho_raw = getattr(ed, "rho", None)
            rho_obs = (
                rho_raw[:, 0, 1] if rho_raw is not None
                else (0.2 / np.where(fr == 0, np.nan, fr)) * np.abs(z[:, 0, 1]) ** 2
            )
            per     = 1.0 / np.where(fr == 0, np.nan, fr)
            per_fwd = 1.0 / np.where(freqs == 0, np.nan, freqs)
            mask    = np.isfinite(per) & (rho_obs > 0)
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


def _plot_joint_section(
    predictions: dict[str, np.ndarray],
    n_layers: int,
    freqs: np.ndarray,
    station_names: list[str],
    modalities: list[str],
) -> Any:
    """Plot the joint-inversion predicted log₁₀ρ section."""
    import matplotlib.pyplot as plt
    from ..api.section import PYCSAMT_SECTION
    from ..api.station import PYCSAMT_STATION_RENDERING

    n_st = len(station_names)
    if n_st == 0:
        return None

    mat = np.full((n_layers, n_st), np.nan)
    for si, nm in enumerate(station_names):
        v = predictions.get(nm)
        if v is None:
            continue
        n = min(len(v), n_layers)
        mat[:n, si] = v[:n]

    ths    = _default_thicknesses(n_layers, freqs)
    depths = np.concatenate([[0], np.cumsum(ths)]) / 1000.0  # km

    section   = PYCSAMT_SECTION.style_for("inversion")
    fig_w, fig_h = section.figsize_for(n_stations=n_st, n_y=n_layers)
    fig, ax   = plt.subplots(figsize=(fig_w, fig_h))

    vv   = mat[np.isfinite(mat)]
    vmin = float(np.percentile(vv, 5))  if vv.size else 0.0
    vmax = float(np.percentile(vv, 95)) if vv.size else 4.0

    im = ax.imshow(
        mat,
        aspect="auto",
        origin="upper",
        extent=(-0.5, n_st - 0.5, depths[-1], depths[0]),
        cmap="jet_r",
        vmin=vmin,
        vmax=vmax,
        interpolation="bilinear",
    )
    PYCSAMT_STATION_RENDERING.apply(
        ax, np.arange(n_st, dtype=float), station_names,
        preset="inversion", xlim=(-0.5, n_st - 0.5),
    )
    ax.set_ylabel("Depth (km)", fontsize=9)
    ax.tick_params(axis="y", labelsize=8)
    section.add_colorbar(im, ax, label="$\\log_{10}\\rho$ (Ω·m)")
    title = f"Joint inversion ({' + '.join(modalities).upper()}) — predicted section"
    ax.set_title(title, fontsize=10, fontweight="bold")
    fig.tight_layout()
    return fig


__all__ = ["JointInversionAgent"]
