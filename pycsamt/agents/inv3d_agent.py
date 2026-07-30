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

* Trains on synthetic 3-D profiles, using one of two ``physics`` modes
  (see below), then predicts on the observed
  :class:`~pycsamt.site.Sites` dataset.

* Outputs log₁₀ρ per depth layer **and** log₁₀h per interface for every
  station, giving a full layered earth model that can be gridded into a
  3-D resistivity volume.

* Optionally runs **MC-dropout uncertainty** (``n_mc`` stochastic passes)
  to produce depth-resolved confidence maps alongside the main prediction.

Physics modes
-------------
``physics="mt1d"`` (default)
    Tiles independent 1-D forward models across the real station
    positions (:func:`~pycsamt.forward.batch.generate_dataset`). This
    is the original smoke/demo path from the AI-inversion
    implementation plan: nothing enforces genuine lateral/3-D coupling
    between the tiled 1-D columns beyond what the GCN's spatial
    smoothing adds after the fact, so the plan explicitly does not
    call this genuine 3-D inversion. Kept unconditionally as the
    default and as an explicit fallback, per the plan's requirement.
``physics="mt3d"``
    Generates genuinely 3-D correlated geological volumes at the
    survey's own real station ``(x, y)`` positions and solves them
    with the research-only, small-grid
    :class:`~pycsamt.forward.maxwell.mt3d.MT3DAdapter`
    (:func:`~pycsamt.ai.training.dataset3d.generate_3d_maxwell_dataset`).
    Unlike ``Inv2DAgent(physics="mt2d")``, this uses the survey's
    *actual* station geometry rather than a synthetic uniform
    spacing, since the real coordinates are already needed to build
    the GCN adjacency graph. Each training profile's per-station
    target is the true 3-D volume's own vertical resistivity column
    at that station, resampled onto the agent's display depth grid
    (nearest-cell; the geological grid's own resolution is coarser
    than most production sections, per :class:`MT3DAdapter`'s
    small-cell-budget research-only status). Only the TE-like
    ``zxy`` response is requested and used, matching the observed
    feature pipeline (:func:`~pycsamt.agents.ai_inversion._z_to_features`,
    which is ``[log10(rho_a_xy), phase_xy]`` only, zero-padded to this
    agent's 4-wide feature slot — using the unused ``zyx``/diagonal
    slots for synthetic data would create a train/observation
    distribution mismatch, not add information). When a held-out
    synthetic split is available, a genuine (non-fabricated) recovery
    check against known-truth resistivity is added to the result via
    :func:`~pycsamt.ai.validation.recovery_report`, since the field
    survey has no ground truth to check against. This mode is far
    more expensive per training profile than ``"mt1d"`` (a real 3-D
    Maxwell solve vs. an independent 1-D solve per station) — lower
    ``n_train_profiles`` accordingly. Mesh accuracy degrades at higher
    frequencies (~10% error by 20-50 Hz vs ~1-2% at 1-2 Hz, at this
    agent's default ``cells_per_skin_depth=None``); passing an
    explicit ``cells_per_skin_depth`` (e.g. 8.0) narrows this to
    ~5-9% at the same cell budget but does not fully close it — see
    :mod:`~pycsamt.ai.training.dataset3d`'s module docstring for
    measured numbers and the cost/accuracy trade-off raising it
    together with ``max_mesh_cells`` buys. Prefer a ``freqs`` grid
    concentrated in the low-frequency range this solver's own
    benchmark validates when accuracy matters more
    than matching a specific field acquisition band.

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
from ._topography import resolve_agent_topography
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
    freqs : array-like or None
        Explicit positive frequency grid in hertz. When supplied, it replaces
        the legacy ``10^-4``--``10^3`` Hz grid and determines ``n_freqs``.
    depth_max : float or None
        Maximum cumulative model depth in metres. When supplied, layer
        thicknesses are geometrically graded and normalized to this depth.
        ``None`` preserves the legacy Bostick-derived display parameterization.
    n_train_profiles : int
        Number of synthetic 3-D training profiles (default 150). With
        ``physics="mt3d"``, each profile costs a real 3-D Maxwell
        solve, not a cheap independent 1-D solve — lower this
        substantially (e.g. 20-40) for that mode.
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
    physics : {"mt1d", "mt3d"}, default "mt1d"
        Synthetic training-data physics; see the module docstring.
    correlation_length_x_m, correlation_length_y_m, correlation_length_z_m
        : (float, float)
        ``physics="mt3d"`` only: horizontal/horizontal/vertical
        correlation length ranges forwarded to
        :class:`~pycsamt.ai.training.dataset3d.Maxwell3DDatasetConfig`.
    log_resistivity_mean, log_resistivity_std : float
        ``physics="mt3d"`` only: affine map from the standardized
        correlated field to ``log10(resistivity_ohm_m)``.
    mesh_safety_factor, max_mesh_cells : float, int
        ``physics="mt3d"`` only: forwarded to
        ``Maxwell3DDatasetConfig``. ``max_mesh_cells`` also builds the
        matching ``MT3DAdapter(max_cells=...)`` used to solve.
    cells_per_skin_depth : float or None, default None
        ``physics="mt3d"`` only: opt-in frequency-aware solver core
        resolution, forwarded to ``Maxwell3DDatasetConfig`` — see its
        docstring for the accuracy-vs-cell-cost trade-off this exists
        to address. ``None`` (the default) preserves this agent's
        original behavior (core resolution = the geological grid's
        own spacing); this agent's default ``freqs`` grid spans up to
        ~1000 Hz (``_DEFAULT_FREQS``), so enabling this without also
        raising ``max_mesh_cells`` and/or narrowing ``freqs`` risks a
        "needs N solver cells" error at high frequencies.
    geology_grid_nx_ny, geology_grid_nz : int, int or None
        ``physics="mt3d"`` only: resolution of the *geological*
        training grid (how finely the true 3-D resistivity volumes
        are drawn), independent of the solver mesh's own
        frequency-aware resolution above. ``geology_grid_nz=None``
        (default) uses ``min(max(n_layers, 4), 8)``. Raising these
        increases realism at the cost of more solver cells per
        realization; see ``max_mesh_cells`` for the resulting budget.

    Input keys
    ----------
    ``sites`` / ``path`` : Sites or str — observed MT dataset
    ``coords`` : ndarray (n_stations, 2), optional — station (x, y) in metres.
        Auto-extracted from EDI lat/lon when absent.
    ``adjacency`` : ndarray (n_stations, n_stations), optional — pre-computed
        normalised adjacency; overrides *radius* when supplied.
    ``output_dir`` : str, optional
    ``freqs`` : array-like, optional — execution-time frequency-grid override.
    ``depth_max`` : float, optional — cumulative finite-layer depth in metres.
    ``topography`` : bool or dict, optional
        ``True`` extracts elevation and chainage from ``sites`` and adds a
        terrain-draped section.  A mapping may instead provide
        ``elevation_m`` and optionally ``chainage_km``, plus
        ``exaggeration`` and ``interp_method`` rendering options.  This is a
        geometry/visualisation contract; it does not alter the MT forward
        solver or GCN prediction.
    ``period_range`` : [T_min, T_max], optional

    Output data keys
    ----------------
    ``pred_rho``          ndarray (n_sta, n_layers)  — log₁₀ρ
    ``pred_thick``        ndarray (n_sta, n_layers-1) — log₁₀h (metres)
    ``pred_uncertainty``  ndarray (n_sta, n_layers) or None  — MC-dropout std
    ``depths_km``         ndarray — depth axis at station midpoints (km)
    ``frequency_grid_hz`` ndarray — effective feature/forward frequency grid
    ``configured_depth_max_m`` float or None — explicit depth contract
    ``topography``       dict or None — resolved terrain metadata
    ``station_names``     list[str]
    ``station_coords``    ndarray (n_sta, 2)  — metres
    ``adjacency``         ndarray (n_sta, n_sta)
    ``rms_global``        float
    ``inverter``          GCNInverter3D
    ``figures``           dict
    ``figure_paths``      dict
    ``physics``           str — ``"mt1d"`` or ``"mt3d"``, mode used
    ``mt3d_recovery``     dict or None — held-out known-truth
                           recovery metrics (``physics="mt3d"`` only)

    Examples
    --------
    >>> agent = Inv3DAgent(n_layers=5, epochs=20, n_mc=10)
    >>> result = agent.execute(
    ...     {
    ...         "path": "/data/WILLY_EDIs",
    ...         "output_dir": "/out/inv3d",
    ...     }
    ... )
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
        freqs: Any = None,
        depth_max: float | None = None,
        n_train_profiles: int = 150,
        epochs: int = 30,
        radius: float = 5_000.0,
        hidden: tuple[int, ...] = (256, 128, 64),
        dropout: float = 0.1,
        n_mc: int = 20,
        physics: str = "mt1d",
        correlation_length_x_m: tuple[float, float] = (500.0, 2000.0),
        correlation_length_y_m: tuple[float, float] = (500.0, 2000.0),
        correlation_length_z_m: tuple[float, float] = (100.0, 500.0),
        log_resistivity_mean: float = 2.0,
        log_resistivity_std: float = 0.5,
        mesh_safety_factor: float = 3.0,
        max_mesh_cells: int = 6_000,
        cells_per_skin_depth: float | None = None,
        geology_grid_nx_ny: int = 6,
        geology_grid_nz: int | None = None,
    ) -> None:
        if physics not in ("mt1d", "mt3d"):
            raise ValueError("physics must be 'mt1d' or 'mt3d'.")
        super().__init__(
            "Inv3DAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="inversion",
        )
        self.n_layers = n_layers
        self.n_freqs = n_freqs
        self.freqs = None if freqs is None else np.asarray(freqs, dtype=float)
        self.depth_max = None if depth_max is None else float(depth_max)
        self.n_train_profiles = n_train_profiles
        self.epochs = epochs
        self.radius = radius
        self.hidden = tuple(hidden)
        self.dropout = dropout
        self.n_mc = n_mc
        self.physics = physics
        self.correlation_length_x_m = correlation_length_x_m
        self.correlation_length_y_m = correlation_length_y_m
        self.correlation_length_z_m = correlation_length_z_m
        self.log_resistivity_mean = float(log_resistivity_mean)
        self.log_resistivity_std = float(log_resistivity_std)
        self.mesh_safety_factor = float(mesh_safety_factor)
        self.max_mesh_cells = int(max_mesh_cells)
        self.cells_per_skin_depth = (
            None
            if cells_per_skin_depth is None
            else float(cells_per_skin_depth)
        )
        self.geology_grid_nx_ny = int(geology_grid_nx_ny)
        self.geology_grid_nz = (
            None if geology_grid_nz is None else int(geology_grid_nz)
        )

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

        freqs_cfg = input_data.get("freqs", self.freqs)
        if freqs_cfg is None:
            freqs = _DEFAULT_FREQS[: self.n_freqs]
        else:
            freqs = np.asarray(freqs_cfg, dtype=float).reshape(-1)
            if (
                freqs.size < 2
                or not np.all(np.isfinite(freqs))
                or np.any(freqs <= 0)
            ):
                return AgentResult.failed(
                    "'freqs' must contain at least two finite positive values.",
                    elapsed=time.time() - t0,
                )
            freqs = np.unique(freqs)
        n_freqs = int(freqs.size)
        n_features = n_freqs * _N_COMP  # flat feature vector per station
        n_out = 2 * self.n_layers - 1  # log10(ρ) layers + log10(h) interfaces
        depth_max_cfg = input_data.get("depth_max", self.depth_max)
        if depth_max_cfg is not None and float(depth_max_cfg) <= 0:
            return AgentResult.failed(
                "'depth_max' must be positive when supplied.",
                elapsed=time.time() - t0,
            )
        ths = _agent_thicknesses(
            self.n_layers,
            freqs,
            None if depth_max_cfg is None else float(depth_max_cfg),
        )

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
        mt3d_dataset = None
        if self.physics == "mt3d":
            self._log.info(
                "Generating %d synthetic 3-D Maxwell realizations "
                "(%d stations × %d freqs)…",
                self.n_train_profiles,
                n_sta,
                n_freqs,
            )
            try:
                mt3d_dataset, X_3d, y_3d = _generate_mt3d_training_data(
                    station_names=station_names,
                    coords_m=coords_m,
                    freqs=freqs,
                    n_layers=self.n_layers,
                    ths=ths,
                    depth_max_m=depth_max_cfg,
                    correlation_length_x_m=self.correlation_length_x_m,
                    correlation_length_y_m=self.correlation_length_y_m,
                    correlation_length_z_m=self.correlation_length_z_m,
                    log_resistivity_mean=self.log_resistivity_mean,
                    log_resistivity_std=self.log_resistivity_std,
                    mesh_safety_factor=self.mesh_safety_factor,
                    max_mesh_cells=self.max_mesh_cells,
                    n_realizations=self.n_train_profiles,
                    n_features=n_features,
                    geology_grid_nx_ny=self.geology_grid_nx_ny,
                    geology_grid_nz=self.geology_grid_nz,
                    cells_per_skin_depth=self.cells_per_skin_depth,
                )
                n_samp = len(mt3d_dataset.select("train"))
            except Exception as exc:
                return AgentResult.failed(
                    f"3-D Maxwell dataset assembly failed: {exc}",
                    elapsed=time.time() - t0,
                )
        else:
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

        # ── known-truth recovery check (physics="mt3d" only) ────────────────
        # The field survey has no ground truth; this instead checks the
        # trained network against held-out synthetic realizations whose
        # true resistivity is known, per the AI-inversion plan's M0
        # baseline-metrics requirement.
        mt3d_recovery: dict[str, Any] | None = None
        if mt3d_dataset is not None:
            try:
                mt3d_recovery = _mt3d_recovery_check(
                    inverter, mt3d_dataset, A, self.n_layers, ths, n_features
                )
            except Exception as exc:
                warnings.append(f"mt3d recovery check failed: {exc}")

        # ── depth axis (in km) ────────────────────────────────────────────────
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

        # Topography is deliberately applied after inversion.  The current
        # mt1d synthetic forward solver has no terrain-aware mesh, so changing
        # the predicted resistivity here would imply physics that was never
        # solved.  Instead, resolve station elevations and return/render the
        # model in an absolute-elevation coordinate system.
        topography = resolve_agent_topography(
            input_data.get("topography", False),
            sites=sites,
            station_names=station_names,
            coords_m=coords_m,
            warnings_list=warnings,
        )
        if topography is not None and topography["applied"]:
            try:
                from ..topo import plot_topo_section

                topo_model = {
                    "pred_rho": pred_rho,
                    "depths_km": depths,
                    "station_names": station_names,
                    "station_coords": coords_m,
                    "rms_global": rms_global,
                }
                ax = plot_topo_section(
                    topo_model,
                    elevation=topography["elevation_m"],
                    chainage=topography["chainage_km"],
                    station_names=station_names,
                    station_x=topography["chainage_km"],
                    topo_source="array",
                    model_unit="km",
                    depth_max=float(depths[-1]),
                    exaggeration=topography["exaggeration"],
                    interp_method=topography["interp_method"],
                    title="GCN inversion draped below station topography",
                )
                fig_topo = ax.get_figure()
                figures["topography_section"] = fig_topo
                p = self._save_figure(
                    fig_topo,
                    output_dir,
                    "inv3d_topography_section",
                    warnings_list=warnings,
                )
                if p:
                    fig_paths["topography_section"] = p
                topography["rendered"] = True
            except Exception as exc:
                topography["rendered"] = False
                warnings.append(f"Topography section: {exc}")

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
        recovery_note = ""
        if mt3d_recovery is not None:
            recovery_note = (
                f" Held-out recovery RMSE={mt3d_recovery['rmse']:.3f} "
                f"(log10 Ω·m, n={mt3d_recovery['n_samples']})."
            )
        return AgentResult(
            status="success",
            summary=(
                f"3-D GCN inversion (physics={self.physics}): "
                f"{n_sta} stations × {self.n_layers} layers. "
                f"{rms_disp}{unc_disp}. {len(figures)} figures."
                f"{recovery_note}"
            ),
            data={
                "pred_rho": pred_rho,
                "pred_thick": pred_thick,
                "pred_uncertainty": pred_uncertainty,
                "depths_km": depths,
                "frequency_grid_hz": freqs,
                "configured_depth_max_m": depth_max_cfg,
                "topography": topography,
                "station_elevation_m": (
                    None if topography is None else topography["elevation_m"]
                ),
                "station_chainage_km": (
                    None if topography is None else topography["chainage_km"]
                ),
                "station_names": station_names,
                "station_coords": coords_m,
                "adjacency": A,
                "rms_global": rms_global,
                "inverter": inverter,
                "n_edges": off_diag_edges,
                "figures": figures,
                "figure_paths": fig_paths,
                "physics": self.physics,
                "mt3d_recovery": mt3d_recovery,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


# ── private helpers ───────────────────────────────────────────────────────────


def _agent_thicknesses(
    n_layers: int,
    freqs: np.ndarray,
    depth_max: float | None,
) -> np.ndarray:
    """Return legacy or explicitly depth-limited interface thicknesses."""
    if depth_max is None:
        return _default_thicknesses(n_layers, freqs)
    n_finite = max(int(n_layers) - 1, 1)
    weights = np.geomspace(1.0, 3.0, n_finite)
    return (float(depth_max) * weights / weights.sum()).astype(float)


def _resolve_agent_topography(
    config: Any,
    *,
    sites: Any,
    station_names: list[str],
    coords_m: np.ndarray,
    warnings_list: list[str],
) -> dict[str, Any] | None:
    """Resolve and align optional terrain data to usable inversion stations."""
    if config is None or config is False:
        return None
    if config is True:
        cfg: dict[str, Any] = {}
    elif isinstance(config, dict):
        cfg = dict(config)
        if not cfg.get("enabled", True):
            return None
    else:
        warnings_list.append(
            "'topography' must be a bool or mapping; ignored."
        )
        return None

    exaggeration = float(cfg.get("exaggeration", 1.0))
    if not np.isfinite(exaggeration) or exaggeration <= 0:
        warnings_list.append(
            "Topography exaggeration must be finite and positive; ignored."
        )
        return None
    interp_method = str(cfg.get("interp_method", "linear")).lower()
    if interp_method not in {"linear", "nearest", "cubic"}:
        warnings_list.append(
            "Topography interp_method must be linear, nearest, or cubic; ignored."
        )
        return None

    explicit_elev = cfg.get("elevation_m")
    explicit_chain = cfg.get("chainage_km")
    source = "array" if explicit_elev is not None else "sites"
    if explicit_elev is None:
        from ..topo.extract import (
            extract_chainage,
            extract_elevation,
            extract_station_names,
        )

        all_names = extract_station_names(sites)
        all_elev = extract_elevation(sites)
        all_chain = extract_chainage(sites)
        # Name alignment is essential because stations without impedance data
        # may have been skipped above.  Positional truncation would attach the
        # wrong hill/valley to every following station.
        lookup: dict[str, int] = {}
        for i, name in enumerate(all_names):
            lookup.setdefault(str(name).strip().casefold(), i)
        missing = [
            n for n in station_names if str(n).strip().casefold() not in lookup
        ]
        if missing:
            warnings_list.append(
                "Topography station-name alignment failed for: "
                + ", ".join(missing[:5])
                + (" ..." if len(missing) > 5 else "")
            )
            return {
                "requested": True,
                "applied": False,
                "rendered": False,
                "source": source,
                "affects_forward_physics": False,
                "vertical_datum": "metres above sea level",
                "elevation_m": None,
                "chainage_km": None,
                "exaggeration": exaggeration,
                "interp_method": interp_method,
            }
        idx = np.asarray(
            [lookup[str(n).strip().casefold()] for n in station_names]
        )
        elevation = np.asarray(all_elev, dtype=float)[idx]
        chainage = np.asarray(all_chain, dtype=float)[idx]
        chainage = chainage - chainage[0]
    else:
        elevation = np.asarray(explicit_elev, dtype=float).reshape(-1)
        if elevation.size != len(station_names):
            warnings_list.append(
                "Explicit elevation_m length must equal the usable station count; "
                "topography ignored."
            )
            return None
        if explicit_chain is None:
            xy = np.asarray(coords_m, dtype=float)
            seg = np.sqrt((np.diff(xy, axis=0) ** 2).sum(axis=1))
            chainage = np.concatenate([[0.0], np.cumsum(seg)]) / 1000.0
        else:
            chainage = np.asarray(explicit_chain, dtype=float).reshape(-1)

    if chainage.size != elevation.size or not np.all(np.isfinite(chainage)):
        warnings_list.append(
            "Topography chainage is invalid; topography ignored."
        )
        return None
    if not np.all(np.isfinite(elevation)) or not np.any(elevation != 0.0):
        warnings_list.append(
            "No finite, non-zero station elevations were available; "
            "topography was not applied."
        )
        return {
            "requested": True,
            "applied": False,
            "rendered": False,
            "source": source,
            "affects_forward_physics": False,
            "vertical_datum": "metres above sea level",
            "elevation_m": elevation,
            "chainage_km": chainage,
            "exaggeration": exaggeration,
            "interp_method": interp_method,
        }
    if np.any(np.diff(chainage) <= 0):
        warnings_list.append(
            "Topography chainage must increase in station order; topography ignored."
        )
        return None

    return {
        "requested": True,
        "applied": True,
        "rendered": False,
        "source": source,
        "affects_forward_physics": False,
        "vertical_datum": "metres above sea level",
        "station_names": list(station_names),
        "elevation_m": elevation,
        "chainage_km": chainage,
        "exaggeration": exaggeration,
        "interp_method": interp_method,
    }


def _extract_station_xy(ed: Any, idx: int) -> np.ndarray:
    """Return (x_m, y_m) for one station from EDI header lat/lon."""
    try:
        coords = getattr(ed, "coords", None)
        if coords is not None and len(coords) >= 2:
            lat, lon = float(coords[0]), float(coords[1])
        else:
            lat = float(
                getattr(ed, "lat", None) or getattr(ed, "latitude", 0.0)
            )
            lon = float(
                getattr(ed, "lon", None) or getattr(ed, "longitude", 0.0)
            )
    except Exception:
        lat, lon = 0.0, 0.0
    if lat == 0.0 and lon == 0.0:
        return np.array(
            [float(idx) * 1.0, 0.0], dtype=np.float64
        )  # fallback: spaced 1 m apart

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


def _generate_mt3d_training_data(
    *,
    station_names: list[str],
    coords_m: np.ndarray,
    freqs: np.ndarray,
    n_layers: int,
    ths: np.ndarray,
    depth_max_m: float | None,
    correlation_length_x_m: tuple[float, float],
    correlation_length_y_m: tuple[float, float],
    correlation_length_z_m: tuple[float, float],
    log_resistivity_mean: float,
    log_resistivity_std: float,
    mesh_safety_factor: float,
    max_mesh_cells: int,
    n_realizations: int,
    n_features: int,
    geology_grid_nx_ny: int,
    geology_grid_nz: int | None,
    cells_per_skin_depth: float | None,
):
    """Build a genuinely 3-D correlated training set at the survey's
    own real station positions and convert it to
    :class:`~pycsamt.ai.inversion.inv3d.GCNInverter3D`'s tensors.

    Unlike ``Inv2DAgent(physics="mt2d")``'s synthetic uniform station
    spacing, this reuses *coords_m* directly (already needed for the
    GCN adjacency graph), so training profiles share the real
    station geometry.

    Returns
    -------
    dataset : Maxwell3DDataset
        Full dataset (train/validation/test), kept for a later
        known-truth recovery check.
    X : ndarray (n_train, n_stations, n_features)
    y : ndarray (n_train, n_stations, 2*n_layers-1)
    """
    from ..ai.geology import GeologyGrid
    from ..ai.training.dataset3d import (
        Maxwell3DDatasetConfig,
        generate_3d_maxwell_dataset,
    )

    xs, ys = coords_m[:, 0], coords_m[:, 1]
    x_span = max(float(xs.max() - xs.min()), 1.0)
    y_span = max(float(ys.max() - ys.min()), 1.0)
    margin_x = max(0.25 * x_span, 300.0)
    margin_y = max(0.25 * y_span, 300.0)
    x_lo, x_hi = float(xs.min()) - margin_x, float(xs.max()) + margin_x
    y_lo, y_hi = float(ys.min()) - margin_y, float(ys.max()) + margin_y

    depth_total_m = (
        float(depth_max_m) if depth_max_m is not None else float(np.sum(ths))
    )
    n_core_xy = int(geology_grid_nx_ny)
    n_core_z = (
        min(max(n_layers, 4), 8)
        if geology_grid_nz is None
        else int(geology_grid_nz)
    )
    grid = GeologyGrid.regular_3d(
        nx=n_core_xy,
        ny=n_core_xy,
        nz=n_core_z,
        dx_m=(x_hi - x_lo) / n_core_xy,
        dy_m=(y_hi - y_lo) / n_core_xy,
        dz_m=depth_total_m / n_core_z,
        x_origin_m=x_lo,
        y_origin_m=y_lo,
    )
    config = Maxwell3DDatasetConfig(
        dataset_id="inv3dagent-mt3d",
        grid=grid,
        correlation_length_x_m=correlation_length_x_m,
        correlation_length_y_m=correlation_length_y_m,
        correlation_length_z_m=correlation_length_z_m,
        frequencies_hz=freqs,
        station_xy_m=coords_m,
        n_realizations=n_realizations,
        seed=0,
        log_resistivity_mean=log_resistivity_mean,
        log_resistivity_std=log_resistivity_std,
        components=("zxy",),
        mesh_safety_factor=mesh_safety_factor,
        max_mesh_cells=max_mesh_cells,
        cells_per_skin_depth=cells_per_skin_depth,
        validation_fraction=0.1,
        test_fraction=0.1,
    )
    dataset = generate_3d_maxwell_dataset(config)
    X, y = _maxwell3d_samples_to_gcn_arrays(
        dataset.select("train"),
        dataset.grid,
        coords_m,
        freqs,
        n_layers,
        ths,
        n_features,
    )
    return dataset, X, y


def _nearest_index(axis: np.ndarray, value: float) -> int:
    """Return the index of *axis*'s entry closest to *value*."""
    return int(np.argmin(np.abs(axis - value)))


def _maxwell3d_samples_to_gcn_arrays(
    samples,
    grid,
    coords_m: np.ndarray,
    freqs: np.ndarray,
    n_layers: int,
    ths: np.ndarray,
    n_features: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert ``Maxwell3DSample`` objects into
    ``GCNInverter3D``'s ``(X, y)`` convention.

    ``X`` uses the same ``[log10(rho_a_xy), phase_xy]`` layout as
    :func:`~pycsamt.agents.ai_inversion._z_to_features` (zero-padded
    to *n_features*), so synthetic training data and the observed
    features share one convention. ``y``'s resistivity half is each
    station's true vertical column, nearest-cell resampled from the
    3-D volume's own (coarser) depth grid onto the agent's display
    depth grid (*ths*); the thickness half is the fixed display grid
    itself, exactly as the ``physics="mt1d"`` path already uses.
    """
    if not samples:
        raise ValueError("samples must contain at least one entry.")
    n_sta = coords_m.shape[0]
    n_freq = len(freqs)
    mu0 = 4.0e-7 * np.pi
    log_h = np.log10(np.clip(ths, 1.0, None)).astype(np.float32)
    layer_depths_m = np.concatenate([[0.0], np.cumsum(ths)])[:n_layers]

    x_idx = np.asarray([_nearest_index(grid.x_m, x) for x in coords_m[:, 0]])
    y_idx = np.asarray([_nearest_index(grid.y_m, y) for y in coords_m[:, 1]])
    z_idx = np.asarray(
        [_nearest_index(grid.z_m, d) for d in layer_depths_m]
    )  # (n_layers,)

    n_out = 2 * n_layers - 1
    X = np.zeros((len(samples), n_sta, n_features), dtype=np.float32)
    y = np.empty((len(samples), n_sta, n_out), dtype=np.float32)
    for i, sample in enumerate(samples):
        zxy = sample.survey.impedance[:, :, 0]  # (n_sta, n_freq)
        f = sample.survey.frequencies_hz[None, :]
        rho_a = np.abs(zxy) ** 2 / (2.0 * np.pi * f * mu0)
        phase = np.degrees(np.angle(zxy))
        feat = np.concatenate(
            [
                np.log10(np.clip(rho_a, 1e-12, None)),
                phase,
            ],
            axis=1,
        ).astype(np.float32)  # (n_sta, 2 * n_freq)
        X[i, :, : min(2 * n_freq, n_features)] = feat[
            :, : min(2 * n_freq, n_features)
        ]

        columns = sample.resistivity_ohm_m[:, y_idx, x_idx]  # (nz, n_sta)
        log_rho_cols = np.log10(columns[z_idx, :]).T  # (n_sta, n_layers)
        y[i, :, :n_layers] = log_rho_cols
        y[i, :, n_layers:] = log_h[None, :]
    return X, y


def _mt3d_recovery_check(
    inverter,
    dataset,
    adjacency: np.ndarray,
    n_layers: int,
    ths: np.ndarray,
    n_features: int,
) -> dict[str, Any] | None:
    """Evaluate known-truth recovery on a held-out dataset partition.

    Prefers the ``test`` partition, falling back to ``validation``;
    returns ``None`` when neither has samples (e.g. too few
    realizations were requested to split three ways).
    """
    from ..ai.validation import recovery_report

    held_out = dataset.select("test") or dataset.select("validation")
    if not held_out:
        return None
    coords_m = np.asarray(held_out[0].survey.coordinates_m)[:, :2]
    freqs = held_out[0].survey.frequencies_hz
    rmse_vals, mae_vals, r2_vals = [], [], []
    for sample in held_out:
        X_i, _ = _maxwell3d_samples_to_gcn_arrays(
            [sample], dataset.grid, coords_m, freqs, n_layers, ths, n_features
        )
        pred = inverter.predict(X_i[0], adjacency=adjacency)  # (n_sta, n_out)
        pred_log_rho = pred[:, :n_layers]

        x_idx = np.asarray(
            [_nearest_index(dataset.grid.x_m, x) for x in coords_m[:, 0]]
        )
        y_idx = np.asarray(
            [_nearest_index(dataset.grid.y_m, y) for y in coords_m[:, 1]]
        )
        layer_depths_m = np.concatenate([[0.0], np.cumsum(ths)])[:n_layers]
        z_idx = np.asarray(
            [_nearest_index(dataset.grid.z_m, d) for d in layer_depths_m]
        )
        columns = sample.resistivity_ohm_m[:, y_idx, x_idx]
        true_log_rho = np.log10(columns[z_idx, :]).T  # (n_sta, n_layers)

        report = recovery_report(
            pred_log_rho, true_log_rho, compute_ssim=False
        )
        rmse_vals.append(report.rmse)
        mae_vals.append(report.mae)
        r2_vals.append(report.r2)
    return {
        "rmse": float(np.mean(rmse_vals)),
        "mae": float(np.mean(mae_vals)),
        "r2": float(np.nanmean(r2_vals)),
        "n_samples": len(held_out),
    }


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


def _profile_distance_km(coords_m: np.ndarray) -> np.ndarray:
    """Project station coordinates onto the dominant profile direction."""
    n_sta = coords_m.shape[0]
    cx = coords_m[:, 0] / 1000.0
    cy = coords_m[:, 1] / 1000.0
    if n_sta <= 1:
        return np.zeros(max(n_sta, 1))
    vec = np.array([cx[-1] - cx[0], cy[-1] - cy[0]])
    vec_len = np.linalg.norm(vec) + 1e-9
    vec /= vec_len
    return np.array(
        [np.dot([cx[i] - cx[0], cy[i] - cy[0]], vec) for i in range(n_sta)]
    )


def _is_profile_geometry(
    coords_m: np.ndarray, ratio_threshold: float = 0.15
) -> bool:
    """Return True when *coords_m* is essentially a 1-D line (a profile),
    rather than a genuine areal 2-D footprint.

    Many AMT/CSAMT surveys (e.g. WILLY L18) place every station on one
    line, so a plan-view E-W/N-S map (Delaunay-triangulated) degenerates
    to a zero-width sliver with nothing to fill outside it.  This is
    detected via the ratio of the two principal spatial extents (SVD of
    the mean-centred coordinates): a ratio below *ratio_threshold* means
    the cross-line spread is negligible next to the along-line spread.
    """
    if coords_m.shape[0] < 3:
        return True
    centred = coords_m - coords_m.mean(axis=0, keepdims=True)
    singular_values = np.linalg.svd(centred, compute_uv=False)
    if singular_values[0] <= 0:
        return True
    return bool((singular_values[1] / singular_values[0]) < ratio_threshold)


def _smooth_depth_axis(
    mat: np.ndarray, depths: np.ndarray, n_upsample: int = 200
) -> tuple[np.ndarray, np.ndarray]:
    """Resample a (n_layers, n_sta) section onto a finer depth axis.

    The agent's layer count (``n_layers``) is a network-capacity/training
    choice, typically kept small (5-10); rendering it directly with
    ``imshow`` therefore produces a coarse, blocky section.  This applies
    a monotone cubic (PCHIP) interpolation independently per station so
    the displayed section reads as a continuous resistivity profile
    without inventing new inversion information between layers.
    """
    from scipy.interpolate import PchipInterpolator

    n_layers, n_sta = mat.shape
    # 'depths' has one entry per layer (depths[0] == 0.0, the survey
    # surface); it is used directly as each row's sample location, same
    # as the pre-upsampling imshow's implicit even spacing between
    # depths[0] and depths[-1].
    layer_depths = np.asarray(depths, dtype=float)[:n_layers]
    if n_layers < 2 or layer_depths.size < 2:
        return mat, depths
    fine_depths = np.linspace(layer_depths[0], layer_depths[-1], n_upsample)
    fine_mat = np.empty((n_upsample, n_sta), dtype=mat.dtype)
    for si in range(n_sta):
        column = mat[:, si]
        finite = np.isfinite(column)
        if finite.sum() < 2:
            fine_mat[:, si] = np.nanmean(column) if finite.any() else np.nan
            continue
        fine_mat[:, si] = PchipInterpolator(
            layer_depths[finite], column[finite]
        )(fine_depths)
    return fine_mat, fine_depths


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

    dist = _profile_distance_km(coords_m)

    section = PYCSAMT_SECTION.style_for("inversion")
    fig_w, fig_h = section.figsize_for(n_stations=n_sta, n_y=n_layers)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    mat = pred_rho.T  # (n_layers, n_sta)
    mat, depth_axis = _smooth_depth_axis(mat, depths)
    vv = mat[np.isfinite(mat)]
    vmin = float(np.percentile(vv, 5)) if vv.size else 0.0
    vmax = float(np.percentile(vv, 95)) if vv.size else 4.0

    im = ax.imshow(
        mat,
        aspect="auto",
        origin="upper",
        extent=(
            dist[0] - 0.25,
            dist[-1] + 0.25,
            depth_axis[-1],
            depth_axis[0],
        ),
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


def _plot_uncertainty_section(
    pred_unc: np.ndarray,  # (n_sta, n_layers)
    coords_m: np.ndarray,
    station_names: list[str],
    depths: np.ndarray,  # km
) -> Any:
    """Uncertainty pseudo-section (profile distance vs depth).

    Used instead of :func:`_plot_uncertainty_depth_map`'s plan-view map
    when the station layout is a profile (see
    :func:`_is_profile_geometry`): a Delaunay-triangulated E-W/N-S map of
    near-collinear stations degenerates to a zero-width sliver, so the
    same distance-vs-depth projection the resistivity section already
    uses is far more informative here.
    """
    import matplotlib.pyplot as plt

    from ..api.section import PYCSAMT_SECTION
    from ..api.station import PYCSAMT_STATION_RENDERING

    n_sta, n_layers = pred_unc.shape
    if n_sta == 0:
        return None

    dist = _profile_distance_km(coords_m)
    section = PYCSAMT_SECTION.style_for("inversion")
    fig_w, fig_h = section.figsize_for(n_stations=n_sta, n_y=n_layers)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    mat = pred_unc.T  # (n_layers, n_sta)
    mat, depth_axis = _smooth_depth_axis(mat, depths)
    vv = mat[np.isfinite(mat)]
    vmax = float(np.percentile(vv, 95)) if vv.size else 1.0

    im = ax.imshow(
        mat,
        aspect="auto",
        origin="upper",
        extent=(
            dist[0] - 0.25,
            dist[-1] + 0.25,
            depth_axis[-1],
            depth_axis[0],
        ),
        cmap="Oranges",
        vmin=0.0,
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
    section.add_colorbar(im, ax, label="σ (log₁₀ρ)")
    ax.set_title(
        "3-D GCN — MC-dropout uncertainty (σ log₁₀ρ)",
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
    """Uncertainty map at the shallowest and deepest depth slices.

    Falls back to :func:`_plot_uncertainty_section` for profile-shaped
    surveys, where a plan-view map has no real cross-line information.
    """
    import matplotlib.pyplot as plt
    from matplotlib.tri import Triangulation

    from ..api.section import PYCSAMT_SECTION

    n_sta, n_layers = pred_unc.shape
    if n_sta == 0:
        return None

    if _is_profile_geometry(coords_m):
        return _plot_uncertainty_section(
            pred_unc, coords_m, station_names, depths
        )

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
