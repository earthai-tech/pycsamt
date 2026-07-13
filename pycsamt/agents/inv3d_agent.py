"""
pycsamt.agents.inv3d_agent
===========================

:class:`Inv3DAgent` — Graph-convolutional 3-D MT spatial inversion.

Wraps :class:`~pycsamt.ai.inversion.inv3d.GCNInverter3D`:

* Represents the survey network as a **spatial graph** whose edges connect
  stations within a configurable radius.  Spectral GCN message-passing
  propagates information between neighbouring stations so the resulting
  3-D resistivity volume is spatially coherent — artefacts from
  station-by-station 1-D inversion are suppressed.

* Trains on **synthetic 3-D profiles** assembled by tiling independent
  1-D forward models across a virtual station grid, then predicts on the
  observed :class:`~pycsamt.site.Sites` dataset.

* Outputs log₁₀ρ per depth layer **and** log₁₀h per interface for every
  station, giving a full layered earth model that can be gridded into a
  3-D resistivity volume.

* Optionally runs **MC-dropout uncertainty** (``n_mc`` stochastic passes)
  to produce depth-resolved confidence maps alongside the main prediction.

Requires PyTorch **or** TensorFlow.

Architecture
------------
:class:`~pycsamt.ai.nets.gcn.GCNNet` — spectral graph convolutional network
(Kipf & Welling 2017).  No external graph library is required.

References
----------
.. [1] Kipf, T. N. & Welling, M. (2017). Semi-supervised classification
       with graph convolutional networks. *ICLR 2017*.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np

from ._base import AgentResult, BaseAgent
from .ai_inversion import _default_thicknesses, _z_to_features

_SYSTEM_PROMPT = """\
You are an expert in 3-D MT inversion using graph-convolutional deep learning.
Given a GCN-based 3-D inversion result, write 4-5 sentences that:
1. Describe the survey geometry (station count, spatial extent, adjacency radius).
2. Interpret the dominant 3-D resistivity structures and their spatial continuity.
3. Assess prediction quality (RMS, depth range) relative to station spacing.
4. Compare the GCN spatial result to independent 1-D predictions where possible.
5. Recommend geological follow-up and areas with highest uncertainty.
Reply in plain scientific English.
"""

_DEFAULT_FREQS = np.logspace(-4, 3, 32)  # 32 frequencies 10⁻⁴ – 10³ Hz
_N_COMP = 4  # Re/Im of Zxy and Zyx


class Inv3DAgent(BaseAgent):
    """3-D MT profile inversion using a graph-convolutional network (GCN).

    Parameters
    ----------
    api_key, model, llm_provider : str
    n_layers : int
        Number of depth layers per station (default 5).
    n_freqs : int
        Number of frequencies used for feature extraction (default 32).
    n_train_profiles : int
        Number of synthetic 3-D training profiles (default 150).
    epochs : int
        Training epochs (default 30).
    radius : float
        Maximum inter-station edge distance in metres for the adjacency graph
        (default 5 000 m).  Stations farther apart than *radius* are
        disconnected in the graph.
    hidden : tuple of int
        GCN hidden-layer sizes (default (256, 128, 64)).
    dropout : float
        Dropout probability (default 0.1); also used for MC uncertainty.
    n_mc : int
        Number of Monte-Carlo dropout passes for uncertainty estimation.
        Set to 0 to skip uncertainty (faster, default 20).

    Input keys
    ----------
    ``sites`` / ``path`` : Sites or str — observed MT dataset
    ``coords`` : ndarray (n_stations, 2), optional — station (x, y) in metres.
        Auto-extracted from EDI lat/lon when absent.
    ``adjacency`` : ndarray (n_stations, n_stations), optional — pre-computed
        normalised adjacency; overrides *radius* when supplied.
    ``output_dir`` : str, optional
    ``period_range`` : [T_min, T_max], optional

    Output data keys
    ----------------
    ``pred_rho``          ndarray (n_sta, n_layers)  — log₁₀ρ
    ``pred_thick``        ndarray (n_sta, n_layers-1) — log₁₀h (metres)
    ``pred_uncertainty``  ndarray (n_sta, n_layers) or None  — MC-dropout std
    ``depths_km``         ndarray — depth axis at station midpoints (km)
    ``station_names``     list[str]
    ``station_coords``    ndarray (n_sta, 2)  — metres
    ``adjacency``         ndarray (n_sta, n_sta)
    ``rms_global``        float
    ``inverter``          GCNInverter3D
    ``figures``           dict
    ``figure_paths``      dict

    Examples
    --------
    >>> agent = Inv3DAgent(n_layers=5, epochs=20, n_mc=10)
    >>> result = agent.execute({
    ...     "path":       "/data/WILLY_EDIs",
    ...     "output_dir": "/out/inv3d",
    ... })
    >>> result["rms_global"]
    0.28
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        n_layers: int = 5,
        n_freqs: int = 32,
        n_train_profiles: int = 150,
        epochs: int = 30,
        radius: float = 5_000.0,
        hidden: tuple[int, ...] = (256, 128, 64),
        dropout: float = 0.1,
        n_mc: int = 20,
    ) -> None:
        super().__init__(
            "Inv3DAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="inversion",
        )
        self.n_layers = n_layers
        self.n_freqs = n_freqs
        self.n_train_profiles = n_train_profiles
        self.epochs = epochs
        self.radius = radius
        self.hidden = tuple(hidden)
        self.dropout = dropout
        self.n_mc = n_mc

    # ── public ────────────────────────────────────────────────────────────────

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        # ── backend check ─────────────────────────────────────────────────────
        try:
            from ..ai.inversion.inv3d import GCNInverter3D
            from ..ai.nets.gcn import build_adjacency
            from ..backends import get_backend_instance
            from ..forward.batch import generate_dataset

            if get_backend_instance() is None:
                raise ImportError("No DL backend.")
        except ImportError as exc:
            return AgentResult.failed(
                f"Inv3DAgent requires PyTorch or TensorFlow: {exc}",
                hint="pip install torch  or  pip install tensorflow",
                elapsed=time.time() - t0,
            )

        from ..emtools._core import (
            _get_z_block,
            _iter_items,
            _name,
            ensure_sites,
        )

        # ── load sites ────────────────────────────────────────────────────────
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

        freqs = _DEFAULT_FREQS[: self.n_freqs]
        n_features = self.n_freqs * _N_COMP  # flat feature vector per station
        n_out = 2 * self.n_layers - 1  # log10(ρ) layers + log10(h) interfaces

        # ── collect observed features + station coordinates ────────────────────
        station_names: list[str] = []
        feat_list: list[np.ndarray] = []
        coord_list: list[np.ndarray] = []  # (x_m, y_m) per station

        ext_coords = input_data.get("coords")  # user-supplied (n_sta, 2)

        for i, ed in enumerate(_iter_items(sites)):
            nm = _name(ed, i)
            z_obj, z, fr = _get_z_block(ed)
            if z is None:
                warnings.append(f"{nm}: no Z data, skipped.")
                continue
            feat = _z_to_features(z_obj, z, fr, freqs)
            if feat is None:
                warnings.append(f"{nm}: feature extraction failed, skipped.")
                continue
            station_names.append(nm)
            feat_list.append(feat.reshape(-1))  # (n_features,)

            # coordinates: try EDI header, then sequential fallback
            if ext_coords is None:
                xy = _extract_station_xy(ed, i)
                coord_list.append(xy)

        n_sta = len(station_names)
        if n_sta < 2:
            return AgentResult.failed(
                f"Only {n_sta} usable station(s) — need ≥ 2 for GCN.",
                elapsed=time.time() - t0,
            )

        X_obs = np.stack(feat_list, axis=0).astype(
            np.float32
        )  # (n_sta, n_feat)
        X_obs = _pad_or_trim(X_obs, n_features)

        # station coordinates (metres)
        if ext_coords is not None:
            coords_m = np.asarray(ext_coords, dtype=np.float64)[:n_sta]
        else:
            coords_m = np.stack(coord_list, axis=0)  # (n_sta, 2) in metres

        # ── build adjacency ───────────────────────────────────────────────────
        user_adj = input_data.get("adjacency")
        if user_adj is not None:
            A = np.asarray(user_adj, dtype=np.float32)
            if A.shape != (n_sta, n_sta):
                warnings.append(
                    f"Supplied adjacency shape {A.shape} ≠ ({n_sta},{n_sta}); "
                    "rebuilding from coordinates."
                )
                A = build_adjacency(coords_m, radius=self.radius)
        else:
            A = build_adjacency(coords_m, radius=self.radius)

        # warn if graph is fully disconnected
        off_diag_edges = int((A > 0).sum()) - n_sta
        if off_diag_edges == 0:
            warnings.append(
                f"No inter-station edges within radius={self.radius:.0f} m.  "
                "The GCN will act as independent 1-D inversion.  "
                "Try increasing radius."
            )

        # ── generate synthetic 3-D training data ───────────────────────────────
        n_1d_total = self.n_train_profiles * n_sta
        self._log.info(
            "Generating %d synthetic profiles (%d×%d) for 3-D GCN training…",
            n_1d_total,
            self.n_train_profiles,
            n_sta,
        )
        try:
            ds = generate_dataset(
                solver="mt1d",
                n_samples=n_1d_total,
                freqs=freqs,
                n_layers=self.n_layers,
                noise_level=0.03,
                seed=42,
                n_jobs=1,
                verbose=False,
            )
            # X_1d: (n_1d, n_freqs, 4) → flatten → (n_1d, n_features)
            X_1d = ds.X.reshape(n_1d_total, -1).astype(np.float32)
            X_1d = _pad_or_trim(X_1d, n_features)

            # y_1d: (n_1d, n_layers) — log10(ρ) only
            # y target for GCN: (n_1d, 2*n_layers-1) = log10(ρ) + log10(h)
            y_log_rho = ds.y[:, : self.n_layers].astype(np.float32)
            ths = _default_thicknesses(self.n_layers, freqs)  # (n_layers-1,)
            log_h = np.log10(np.clip(ths, 1.0, None)).astype(np.float32)
            log_h_tile = np.tile(
                log_h[None, :], (n_1d_total, 1)
            )  # (n_1d, n_layers-1)
            y_1d = np.concatenate(
                [y_log_rho, log_h_tile], axis=1
            )  # (n_1d, n_out)

            # reshape into 3-D profile tensors
            n_samp = n_1d_total // n_sta
            X_1d = X_1d[: n_samp * n_sta]
            y_1d = y_1d[: n_samp * n_sta]
            X_3d = X_1d.reshape(
                n_samp, n_sta, n_features
            )  # (n_samp, n_sta, n_feat)
            y_3d = y_1d.reshape(
                n_samp, n_sta, n_out
            )  # (n_samp, n_sta, n_out)

            # synthetic adjacency: same A for all profiles
            # (GCNInverter3D uses one A across the mini-batch)

        except Exception as exc:
            return AgentResult.failed(
                f"Synthetic dataset generation failed: {exc}",
                elapsed=time.time() - t0,
            )

        # ── train GCNInverter3D ───────────────────────────────────────────────
        self._log.info(
            "Training GCNInverter3D (GCN hidden=%s) for %d epochs…",
            self.hidden,
            self.epochs,
        )
        try:
            inverter = GCNInverter3D(
                n_features=n_features,
                n_layers=self.n_layers,
                hidden=self.hidden,
                dropout=self.dropout,
            )
            inverter.fit(
                X_3d,
                y_3d,
                adjacency=A,
                epochs=self.epochs,
                batch_size=max(4, min(16, n_samp // 10)),
                patience=max(5, self.epochs // 5),
                verbose=False,
            )
        except Exception as exc:
            return AgentResult.failed(
                f"GCNInverter3D training failed: {exc}",
                elapsed=time.time() - t0,
            )

        # ── predict on observed stations ──────────────────────────────────────
        pred_uncertainty: np.ndarray | None = None
        try:
            y_pred = inverter.predict(X_obs, adjacency=A)  # (n_sta, n_out)
        except Exception as exc:
            return AgentResult.failed(
                f"3-D prediction failed: {exc}",
                elapsed=time.time() - t0,
            )

        pred_rho = y_pred[:, : self.n_layers]  # (n_sta, n_layers)
        pred_thick = y_pred[:, self.n_layers :]  # (n_sta, n_layers-1)

        # optional MC-dropout uncertainty
        if self.n_mc > 0:
            try:
                mu, sigma = inverter.predict_with_uncertainty(
                    X_obs,
                    adjacency=A,
                    n_mc=self.n_mc,
                )
                pred_uncertainty = sigma[
                    :, : self.n_layers
                ]  # (n_sta, n_layers)
            except Exception as exc:
                warnings.append(f"MC-dropout uncertainty failed: {exc}")

        # ── depth axis (in km) ────────────────────────────────────────────────
        ths = _default_thicknesses(self.n_layers, freqs)
        depths = np.concatenate([[0.0], np.cumsum(ths)]) / 1000.0  # km

        # ── per-station forward RMS ────────────────────────────────────────────
        rms_list: list[float] = []
        for si, (nm, ed) in enumerate(
            _iter_station_items(sites, station_names)
        ):
            rms = _forward_rms_3d(ed, pred_rho[si], ths, freqs)
            if rms is not None:
                rms_list.append(rms)
        rms_global = float(np.nanmean(rms_list)) if rms_list else np.nan

        # ── figures ───────────────────────────────────────────────────────────
        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}

        try:
            fig_slices = _plot_depth_slices(
                pred_rho,
                coords_m,
                station_names,
                depths,
                self.n_layers,
            )
            if fig_slices is not None:
                figures["depth_slices"] = fig_slices
                p = self._save_figure(
                    fig_slices,
                    output_dir,
                    "inv3d_depth_slices",
                    warnings_list=warnings,
                )
                if p:
                    fig_paths["depth_slices"] = p
        except Exception as exc:
            warnings.append(f"Depth slice figure: {exc}")

        try:
            fig_sec = _plot_resistivity_section(
                pred_rho,
                station_names,
                depths,
                coords_m,
            )
            if fig_sec is not None:
                figures["resistivity_section"] = fig_sec
                p = self._save_figure(
                    fig_sec,
                    output_dir,
                    "inv3d_resistivity_section",
                    warnings_list=warnings,
                )
                if p:
                    fig_paths["resistivity_section"] = p
        except Exception as exc:
            warnings.append(f"Resistivity section figure: {exc}")

        if pred_uncertainty is not None:
            try:
                fig_unc = _plot_uncertainty_depth_map(
                    pred_uncertainty,
                    coords_m,
                    station_names,
                    depths,
                )
                if fig_unc is not None:
                    figures["uncertainty_map"] = fig_unc
                    p = self._save_figure(
                        fig_unc,
                        output_dir,
                        "inv3d_uncertainty",
                        warnings_list=warnings,
                    )
                    if p:
                        fig_paths["uncertainty_map"] = p
            except Exception as exc:
                warnings.append(f"Uncertainty map figure: {exc}")

        # ── LLM interpretation ────────────────────────────────────────────────
        interp: str | None = None
        if self.api_key:
            rho_mean = float(np.nanmean(10**pred_rho))
            rho_std = float(np.nanstd(10**pred_rho))
            extent_km = (
                float(
                    np.max(
                        np.linalg.norm(
                            coords_m - coords_m.mean(axis=0), axis=1
                        )
                    )
                )
                / 1000.0
            )
            rms_str = (
                f"{rms_global:.3f}" if not np.isnan(rms_global) else "N/A"
            )
            prompt = (
                f"3-D GCN inversion summary:\n"
                f"  Stations: {n_sta}, extent: ~{extent_km:.1f} km\n"
                f"  Adjacency radius: {self.radius / 1000:.1f} km, "
                f"  edges per station: {off_diag_edges}/{n_sta}\n"
                f"  Layers: {self.n_layers}, max depth: {depths[-1]:.2f} km\n"
                f"  Mean resistivity: {rho_mean:.0f} Ω·m ± {rho_std:.0f}\n"
                f"  Global RMS: {rms_str} log₁₀(Ω·m)\n"
                f"  MC uncertainty: {'computed' if pred_uncertainty is not None else 'skipped'}\n"
                f"  Warnings: {warnings[:3] if warnings else 'none'}\n\n"
                "Interpret the 3-D resistivity volume geologically."
            )
            interp = self.query_llm(prompt, max_tokens=280)

        elapsed = time.time() - t0
        rms_disp = (
            f"RMS {rms_global:.3f}" if not np.isnan(rms_global) else "RMS N/A"
        )
        unc_disp = ", MC σ computed" if pred_uncertainty is not None else ""
        return AgentResult(
            status="success",
            summary=(
                f"3-D GCN inversion: {n_sta} stations × {self.n_layers} layers. "
                f"{rms_disp}{unc_disp}. {len(figures)} figures."
            ),
            data={
                "pred_rho": pred_rho,
                "pred_thick": pred_thick,
                "pred_uncertainty": pred_uncertainty,
                "depths_km": depths,
                "station_names": station_names,
                "station_coords": coords_m,
                "adjacency": A,
                "rms_global": rms_global,
                "inverter": inverter,
                "n_edges": off_diag_edges,
                "figures": figures,
                "figure_paths": fig_paths,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


# ── private helpers ───────────────────────────────────────────────────────────


def _extract_station_xy(ed: Any, idx: int) -> np.ndarray:
    """Return (x_m, y_m) for one station from EDI header lat/lon."""
    try:
        lat = float(getattr(ed, "lat", None) or getattr(ed, "latitude", 0.0))
        lon = float(getattr(ed, "lon", None) or getattr(ed, "longitude", 0.0))
    except Exception:
        lat, lon = 0.0, float(idx) * 0.001  # fallback: spaced 1 m apart

    # approximate UTM-like local projection (metres)
    x_m = lon * 111_320.0 * np.cos(np.radians(lat))
    y_m = lat * 110_574.0
    return np.array([x_m, y_m], dtype=np.float64)


def _pad_or_trim(X: np.ndarray, target_cols: int) -> np.ndarray:
    """Ensure X has exactly *target_cols* columns."""
    n = X.shape[0]
    c = X.shape[1]
    if c == target_cols:
        return X
    if c > target_cols:
        return X[:, :target_cols]
    pad = np.zeros((n, target_cols - c), dtype=X.dtype)
    return np.concatenate([X, pad], axis=1)


def _iter_station_items(sites: Any, names: list[str]):
    """Yield (name, EDI_object) only for stations in *names*."""
    from ..emtools._core import _iter_items, _name

    wanted = set(names)
    for i, ed in enumerate(_iter_items(sites)):
        nm = _name(ed, i)
        if nm in wanted:
            yield nm, ed


def _forward_rms_3d(
    ed: Any,
    log_rho: np.ndarray,  # (n_layers,)
    ths: np.ndarray,  # (n_layers-1,)
    freqs: np.ndarray,
) -> float | None:
    """Compute forward RMS for one predicted station."""
    try:
        from ..emtools._core import _get_z_block
        from ..forward import LayeredModel, MT1DForward

        _, z, fr = _get_z_block(ed)
        if z is None or fr is None:
            return None

        rhos = 10.0**log_rho
        lm = LayeredModel(resistivity=rhos, thickness=ths)
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


def _plot_depth_slices(
    pred_rho: np.ndarray,  # (n_sta, n_layers)
    coords_m: np.ndarray,  # (n_sta, 2)
    station_names: list[str],
    depths: np.ndarray,  # (n_layers+1,) in km
    n_layers: int,
) -> Any:
    """Horizontal depth-slice maps at 3 representative depths."""
    import matplotlib.pyplot as plt
    from matplotlib.tri import Triangulation

    from ..api.section import PYCSAMT_SECTION

    n_sta = len(station_names)
    PYCSAMT_SECTION.style_for("inversion")
    n_slices = min(3, n_layers)
    layer_idxs = np.linspace(0, n_layers - 1, n_slices, dtype=int)

    fig, axes = plt.subplots(1, n_slices, figsize=(5.0 * n_slices, 5.0))
    if n_slices == 1:
        axes = [axes]

    x_m = coords_m[:, 0]
    y_m = coords_m[:, 1]

    # normalise to km for display
    cx = x_m / 1000.0
    cy = y_m / 1000.0

    vv = pred_rho[np.isfinite(pred_rho)]
    vmin = float(np.percentile(vv, 5)) if vv.size else 0.0
    vmax = float(np.percentile(vv, 95)) if vv.size else 4.0

    for ax, li in zip(axes, layer_idxs):
        layer_vals = pred_rho[:, li]
        depth_km = float(depths[li])

        if n_sta >= 4:
            # Triangulated interpolation
            try:
                tri = Triangulation(cx, cy)
                tcf = ax.tricontourf(
                    tri,
                    layer_vals,
                    levels=20,
                    cmap="jet_r",
                    vmin=vmin,
                    vmax=vmax,
                )
                plt.colorbar(tcf, ax=ax, label="log₁₀ρ (Ω·m)", shrink=0.8)
            except Exception:
                sc = ax.scatter(
                    cx,
                    cy,
                    c=layer_vals,
                    cmap="jet_r",
                    vmin=vmin,
                    vmax=vmax,
                    s=80,
                    edgecolors="k",
                    lw=0.4,
                )
                plt.colorbar(sc, ax=ax, label="log₁₀ρ (Ω·m)", shrink=0.8)
        else:
            sc = ax.scatter(
                cx,
                cy,
                c=layer_vals,
                cmap="jet_r",
                vmin=vmin,
                vmax=vmax,
                s=120,
                edgecolors="k",
                lw=0.5,
            )
            plt.colorbar(sc, ax=ax, label="log₁₀ρ (Ω·m)", shrink=0.8)

        # station labels
        for nm, xi, yi in zip(station_names, cx, cy):
            ax.text(
                xi,
                yi,
                nm,
                fontsize=5.5,
                ha="center",
                va="bottom",
                color="white",
                fontweight="bold",
            )

        ax.set_xlabel("E–W (km)", fontsize=8)
        ax.set_ylabel("N–S (km)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.set_title(
            f"Depth slice  {depth_km:.2f} km", fontsize=8, fontweight="bold"
        )
        ax.set_aspect("equal")

    fig.suptitle(
        "3-D GCN inversion — horizontal depth slices",
        fontsize=10,
        fontweight="bold",
    )
    fig.tight_layout()
    return fig


def _plot_resistivity_section(
    pred_rho: np.ndarray,  # (n_sta, n_layers)
    station_names: list[str],
    depths: np.ndarray,  # km
    coords_m: np.ndarray,  # (n_sta, 2)
) -> Any:
    """Pseudo-section along the survey profile (distance vs depth)."""
    import matplotlib.pyplot as plt

    from ..api.section import PYCSAMT_SECTION
    from ..api.station import PYCSAMT_STATION_RENDERING

    n_sta, n_layers = pred_rho.shape
    if n_sta == 0:
        return None

    # project stations onto the dominant profile direction
    cx = coords_m[:, 0] / 1000.0
    cy = coords_m[:, 1] / 1000.0
    if n_sta > 1:
        vec = np.array([cx[-1] - cx[0], cy[-1] - cy[0]])
        vec_len = np.linalg.norm(vec) + 1e-9
        vec /= vec_len
        dist = np.array(
            [
                np.dot([cx[i] - cx[0], cy[i] - cy[0]], vec)
                for i in range(n_sta)
            ]
        )
    else:
        dist = np.zeros(1)

    section = PYCSAMT_SECTION.style_for("inversion")
    fig_w, fig_h = section.figsize_for(n_stations=n_sta, n_y=n_layers)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    mat = pred_rho.T  # (n_layers, n_sta)
    vv = mat[np.isfinite(mat)]
    vmin = float(np.percentile(vv, 5)) if vv.size else 0.0
    vmax = float(np.percentile(vv, 95)) if vv.size else 4.0

    im = ax.imshow(
        mat,
        aspect="auto",
        origin="upper",
        extent=(dist[0] - 0.25, dist[-1] + 0.25, depths[-1], depths[0]),
        cmap="jet_r",
        vmin=vmin,
        vmax=vmax,
        interpolation="bilinear",
    )
    PYCSAMT_STATION_RENDERING.apply(
        ax,
        dist,
        station_names,
        preset="inversion",
        xlim=(dist[0] - 0.25, dist[-1] + 0.25),
    )
    ax.set_ylabel("Depth (km)", fontsize=9)
    ax.set_xlabel("Profile distance (km)", fontsize=9)
    ax.tick_params(labelsize=8)
    section.add_colorbar(im, ax, label="$\\log_{10}\\rho$ (Ω·m)")
    ax.set_title(
        "3-D GCN inversion — resistivity section",
        fontsize=10,
        fontweight="bold",
    )
    fig.tight_layout()
    return fig


def _plot_uncertainty_depth_map(
    pred_unc: np.ndarray,  # (n_sta, n_layers)
    coords_m: np.ndarray,
    station_names: list[str],
    depths: np.ndarray,  # km
) -> Any:
    """Uncertainty map at the shallowest and deepest depth slices."""
    import matplotlib.pyplot as plt
    from matplotlib.tri import Triangulation

    from ..api.section import PYCSAMT_SECTION

    n_sta, n_layers = pred_unc.shape
    if n_sta == 0:
        return None

    PYCSAMT_SECTION.style_for("inversion")
    n_panels = min(2, n_layers)
    layer_idxs = [0, n_layers - 1] if n_panels == 2 else [0]
    titles = ["Shallow uncertainty", "Deep uncertainty"]

    fig, axes = plt.subplots(1, n_panels, figsize=(5.0 * n_panels, 4.5))
    if n_panels == 1:
        axes = [axes]

    cx = coords_m[:, 0] / 1000.0
    cy = coords_m[:, 1] / 1000.0
    sv = pred_unc[np.isfinite(pred_unc)]
    s_vmax = float(np.percentile(sv, 95)) if sv.size else 1.0

    for ax, li, title in zip(axes, layer_idxs, titles):
        unc_vals = pred_unc[:, li]
        depth_km = float(depths[li])

        if n_sta >= 4:
            try:
                tri = Triangulation(cx, cy)
                tcf = ax.tricontourf(
                    tri,
                    unc_vals,
                    levels=15,
                    cmap="Oranges",
                    vmin=0.0,
                    vmax=s_vmax,
                )
                plt.colorbar(tcf, ax=ax, label="σ (log₁₀ρ)", shrink=0.8)
            except Exception:
                sc = ax.scatter(
                    cx,
                    cy,
                    c=unc_vals,
                    cmap="Oranges",
                    vmin=0.0,
                    vmax=s_vmax,
                    s=80,
                    edgecolors="k",
                    lw=0.4,
                )
                plt.colorbar(sc, ax=ax, label="σ (log₁₀ρ)", shrink=0.8)
        else:
            sc = ax.scatter(
                cx,
                cy,
                c=unc_vals,
                cmap="Oranges",
                vmin=0.0,
                vmax=s_vmax,
                s=120,
                edgecolors="k",
                lw=0.5,
            )
            plt.colorbar(sc, ax=ax, label="σ (log₁₀ρ)", shrink=0.8)

        ax.set_xlabel("E–W (km)", fontsize=8)
        ax.set_ylabel("N–S (km)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.set_title(
            f"{title}  ({depth_km:.2f} km)", fontsize=8, fontweight="bold"
        )
        ax.set_aspect("equal")

    fig.suptitle(
        "3-D GCN — MC-dropout uncertainty (σ log₁₀ρ)",
        fontsize=10,
        fontweight="bold",
    )
    fig.tight_layout()
    return fig


__all__ = ["Inv3DAgent"]
