"""
pycsamt.agents.inv2d_agent
===========================

:class:`Inv2DAgent` — U-Net based 2-D MT profile inversion.

Wraps :class:`~pycsamt.ai.inversion.inv2d.EMInverter2D`:

* Assembles a **2-D pseudosection** image from the observed Sites (station ×
  frequency × impedance component) as the U-Net input.
* Generates a matching synthetic **2-D training dataset**, using one of two
  ``physics`` modes (see below).
* Trains the U-Net and produces a **resistivity section** output:
  (n_depth × n_stations) in log₁₀ Ω·m.
* Visualises the result with
  :func:`~pycsamt.ai.plot.inversion.plot_inversion_result_2d`.

The U-Net treats the whole profile at once, so it naturally captures
lateral continuity — a key advantage over station-by-station 1-D inversion.

Physics modes
-------------
``physics="mt1d"`` (default)
    Tiles independent 1-D forward models into profile-shaped arrays
    (:func:`~pycsamt.forward.batch.generate_dataset`). This is the
    original smoke/demo path from the AI-inversion implementation
    plan: nothing enforces lateral continuity between the tiled 1-D
    columns, so the plan explicitly does not call it genuine 2-D
    inversion. Kept unconditionally as the default and as an explicit
    fallback, per the plan's requirement.
``physics="mt2d"``
    Generates genuinely 2-D correlated training models and solves
    them with the verified
    :class:`~pycsamt.forward.maxwell.mt2d.MT2DAdapter`
    (:func:`~pycsamt.ai.training.dataset2d.generate_2d_maxwell_dataset`),
    and can add spatial regularization
    (:meth:`~pycsamt.ai.inversion.inv2d.EMInverter2D.fit`'s
    ``lambda_x``/``lambda_z``/``lambda_tv``). Only the TE-mode
    (``zxy``) response is requested here, even though
    :mod:`pycsamt.ai.training.dataset2d` now validates the TM-mode
    (``zyx``) response too (see its module docstring): this agent's
    feature layout is still TE-only, so requesting ``zyx`` here would
    only cost an extra, unused forward solve per realization, not
    change accuracy. When a held-out synthetic split is available,
    a genuine (non-fabricated) recovery check against known-truth
    resistivity is added to the result via
    :func:`~pycsamt.ai.validation.recovery_report`, since the field
    survey has no ground truth to check against.
``physics="mt2d_tri"``
    Generates triangular-mesh correlated training models
    (:func:`~pycsamt.ai.training.dataset2d_tri.generate_2d_tri_maxwell_dataset`)
    solved with :class:`~pycsamt.forward.maxwell.mare2dem.Mare2DEMAdapter`
    -- **unverified physics**, see that adapter's own module docstring;
    it needs a real compiled MARE2DEM binary this environment does not
    have. Uses the same synthetic uniform ``station_spacing_m`` geometry
    ``physics="mt2d"`` already documents (real station coordinates are
    not read here either). Output is a per-triangle
    ``log10(resistivity)`` field, not a ``(depth, station)`` grid, so it
    is trained and rendered differently from the other two modes: a
    graph-convolutional inverter
    (:class:`~pycsamt.ai.inversion.inv3d.GCNInverter3D` with
    ``n_layers=1``, reused rather than duplicated -- its architecture
    has nothing 3-D/station-specific about it, only its ``n_layers``
    naming, here repurposed to mean "one scalar per graph node" instead
    of "one depth column per station") over the mesh's
    triangle-centroid adjacency graph, rendered with
    :func:`~pycsamt.api.mesh.draw_tri_mesh`. Not plotted via
    :func:`~pycsamt.ai.plot.inversion.plot_inversion_result_2d` (that
    function is rectilinear-only) or draped with real topography (that
    needs real station coordinates, which this mode -- like
    ``physics="mt2d"`` -- does not read).

Requires PyTorch **or** TensorFlow (``physics="mt2d"``'s spatial
regularization is PyTorch-only, and ``physics="mt2d_tri"``'s GCN
inverter always needs one; see ``EMInverter2D.fit``/``GCNInverter3D``).
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np

from ..compat.sklearn import validate_params
from ._base import AgentResult, BaseAgent
from ._param_schemas import INV2D_PARAM_CONSTRAINTS
from .ai_inversion import _default_thicknesses, _z_to_features

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

_DEFAULT_FREQS_2D = np.logspace(-4, 3, 32)  # 32 frequencies — U-Net input size


class Inv2DAgent(BaseAgent):
    """2-D MT profile inversion using a U-Net convolutional architecture.

    Parameters
    ----------
    api_key, model, llm_provider : str
    n_depth : int
        Number of depth cells in the output section (default 40).
    n_freqs : int
        Number of input frequencies (default 32).
    freqs : array-like or None
        Explicit positive frequency grid in hertz.
    depth_max : float or None
        Maximum cumulative model depth in metres. ``None`` preserves the
        legacy frequency-derived parameterization.
    n_components : int
        Number of channels in the input pseudosection (default 2:
        log10 apparent resistivity and phase for the xy component).
    arch : str
        U-Net variant (default ``"unet"``).
    n_train_profiles : int
        Number of synthetic 2-D profiles for training (default 200).
    n_stations_per_profile : int
        Stations per synthetic profile (default 20).
    epochs : int
        Training epochs (default 30).
    patience : int or None
        Early-stopping patience. ``None`` uses one fifth of ``epochs`` with
        a minimum of five. Set greater than ``epochs`` only for a controlled
        fixed-epoch experiment; validation-based stopping is normally safer.
    physics : {"mt1d", "mt2d", "mt2d_tri"}, default "mt1d"
        Synthetic training-data physics; see the module docstring.
    station_spacing_m : float, default 500.0
        Uniform synthetic station spacing used by ``physics="mt2d"``
        and ``physics="mt2d_tri"``; the actual survey's real station
        geometry is not read by either.
    correlation_length_x_m, correlation_length_z_m : (float, float)
        ``physics="mt2d"``/``"mt2d_tri"`` only: horizontal/vertical
        correlation length ranges forwarded to
        :class:`~pycsamt.ai.training.dataset2d.Maxwell2DDatasetConfig`/
        :class:`~pycsamt.ai.training.dataset2d_tri.MaxwellTri2DDatasetConfig`.
    log_resistivity_mean, log_resistivity_std : float
        ``physics="mt2d"``/``"mt2d_tri"`` only: affine map from the
        standardized correlated field to ``log10(resistivity_ohm_m)``.
    mesh_safety_factor, max_mesh_cells : float, int
        ``physics="mt2d"`` only: forwarded to
        ``Maxwell2DDatasetConfig``.
    lambda_x, lambda_z, lambda_tv : float, default 0.0
        Spatial-regularization weights forwarded to
        :meth:`~pycsamt.ai.inversion.inv2d.EMInverter2D.fit`. Zero by
        default, so nothing changes unless explicitly requested.
        ``physics="mt2d"`` only.
    mesh_target_cell_m, field_grid_cell_m : float
        ``physics="mt2d_tri"`` only: forwarded to
        ``MaxwellTri2DDatasetConfig``.
    topo_x_m, topo_z_m : array-like, optional
        ``physics="mt2d_tri"`` only: real topography (z positive down)
        forwarded to ``MaxwellTri2DDatasetConfig``/
        :func:`~pycsamt.forward.maxwell.tri_mesh_gen.build_graded_tri_mesh`.
        Both default to ``None`` (flat surface at ``z=0``, unchanged from
        before this parameter existed) -- when given, training stations
        sit at their true interpolated elevation instead. This builds a
        real topography-following *training mesh*; it is unrelated to
        the ``topography`` execute()-time input key above, which only
        re-renders an already-flat ``mt1d``/``mt2d`` prediction in an
        absolute-elevation display frame.
    gcn_hidden : tuple of int, default (64, 32, 16)
        ``physics="mt2d_tri"`` only: hidden-layer widths forwarded to
        :class:`~pycsamt.ai.inversion.inv3d.GCNInverter3D`.
    gcn_adjacency_radius_m : float, default 300.0
        ``physics="mt2d_tri"`` only: triangle-centroid adjacency
        radius forwarded to
        :func:`~pycsamt.ai.nets.gcn.build_adjacency`. Must be set
        relative to the actual mesh scale, not left at the default: it
        is a hard cutoff in the same metres as ``mesh_target_cell_m``,
        and if it is smaller than the typical distance between
        neighbouring triangle centroids, ``build_adjacency`` returns
        the identity matrix (every triangle connected only to itself)
        and the GCN silently degenerates into a per-triangle lookup
        with no spatial message-passing at all -- no error is raised.
        A radius of roughly 1.5-2x ``mesh_target_cell_m`` is a
        reasonable starting point; verify with
        ``build_adjacency(mesh.triangle_centroids_m, radius).sum() >
        mesh.n_triangles`` (more than just the self-loops) before
        trusting a training run.
    mare2dem_adapter : object, optional
        ``physics="mt2d_tri"`` only: pre-built
        :class:`~pycsamt.forward.maxwell.mare2dem.Mare2DEMAdapter`
        (e.g. pointed at a specific compiled binary). Defaults to
        ``Mare2DEMAdapter()``, resolved from the environment.

    Input keys
    ----------
    ``sites`` / ``path`` : Sites or str
    ``output_dir`` : str, optional
    ``topography`` : bool or dict, optional
        Extract terrain from ``sites`` and render the predicted section in an
        absolute-elevation frame. A mapping can provide ``elevation_m`` and
        ``chainage_km`` plus ``exaggeration`` and ``interp_method``.
    ``period_range`` : [T_min, T_max], optional
    ``freqs`` : array-like, optional — execution-time frequency override
    ``depth_max`` : float, optional — cumulative model depth in metres

    Output data keys
    ----------------
    ``pred_section``      ndarray (n_depth × n_stations) — log₁₀ ρ
                           (``physics="mt1d"``/``"mt2d"`` only)
    ``pred_triangles``    dict with ``"mesh"``/``"log10_resistivity"``
                           (``physics="mt2d_tri"`` only)
    ``depths_km``         ndarray — depth axis (km)
                           (``physics="mt1d"``/``"mt2d"`` only)
    ``station_names``     list[str]
    ``topography``        dict or None — resolved terrain metadata
                           (``physics="mt1d"``/``"mt2d"`` only)
    ``rms_global``        float (``physics="mt1d"``/``"mt2d"`` only)
    ``inverter``          EMInverter2D or GCNInverter3D
    ``figures``           dict
    ``figure_paths``      dict
    ``physics``           str — ``"mt1d"``, ``"mt2d"``, or
                           ``"mt2d_tri"``, mode used
    ``mt2d_recovery``     dict or None — held-out known-truth
                           recovery metrics (``physics="mt2d"`` only)
    ``mt2d_tri_recovery`` dict or None — held-out known-truth
                           recovery metrics (``physics="mt2d_tri"`` only)
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    @validate_params(INV2D_PARAM_CONSTRAINTS, prefer_skip_nested_validation=True)
    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        n_depth: int = 40,
        n_freqs: int = 32,
        freqs: Any = None,
        depth_max: float | None = None,
        n_components: int = 2,
        arch: str = "unet",
        n_train_profiles: int = 200,
        n_stations_per_profile: int = 20,
        epochs: int = 30,
        patience: int | None = None,
        physics: str = "mt1d",
        station_spacing_m: float = 500.0,
        correlation_length_x_m: tuple[float, float] = (500.0, 2000.0),
        correlation_length_z_m: tuple[float, float] = (100.0, 500.0),
        log_resistivity_mean: float = 2.0,
        log_resistivity_std: float = 0.5,
        mesh_safety_factor: float = 8.0,
        max_mesh_cells: int = 200_000,
        lambda_x: float = 0.0,
        lambda_z: float = 0.0,
        lambda_tv: float = 0.0,
        mesh_target_cell_m: float = 100.0,
        field_grid_cell_m: float = 50.0,
        topo_x_m: Any = None,
        topo_z_m: Any = None,
        gcn_hidden: tuple[int, ...] = (64, 32, 16),
        gcn_adjacency_radius_m: float = 300.0,
        mare2dem_adapter: Any | None = None,
        verbose: bool | int | str = False,
    ) -> None:
        super().__init__(
            "Inv2DAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="inversion",
            verbose=verbose,
        )
        self.n_depth = n_depth
        self.n_freqs = n_freqs
        self.freqs = None if freqs is None else np.asarray(freqs, dtype=float)
        self.depth_max = None if depth_max is None else float(depth_max)
        self.n_components = n_components
        self.arch = arch
        self.n_train_profiles = n_train_profiles
        self.n_stations_per_profile = n_stations_per_profile
        self.epochs = epochs
        self.patience = (
            None if patience is None else max(1, int(patience))
        )
        self.physics = physics
        self.station_spacing_m = float(station_spacing_m)
        self.correlation_length_x_m = correlation_length_x_m
        self.correlation_length_z_m = correlation_length_z_m
        self.log_resistivity_mean = float(log_resistivity_mean)
        self.log_resistivity_std = float(log_resistivity_std)
        self.mesh_safety_factor = float(mesh_safety_factor)
        self.max_mesh_cells = int(max_mesh_cells)
        self.lambda_x = float(lambda_x)
        self.lambda_z = float(lambda_z)
        self.lambda_tv = float(lambda_tv)
        self.mesh_target_cell_m = float(mesh_target_cell_m)
        self.field_grid_cell_m = float(field_grid_cell_m)
        self.topo_x_m = topo_x_m
        self.topo_z_m = topo_z_m
        self.gcn_hidden = tuple(int(h) for h in gcn_hidden)
        self.gcn_adjacency_radius_m = float(gcn_adjacency_radius_m)
        self.mare2dem_adapter = mare2dem_adapter

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        # ── backend check ──────────────────────────────────────────────────────
        try:
            from ..ai.inversion.inv2d import EMInverter2D
            from ..backends import get_backend_instance
            from ..forward.batch import generate_dataset

            if get_backend_instance() is None:
                raise ImportError("No DL backend.")
        except ImportError as exc:
            return AgentResult.failed(
                f"Inv2DAgent requires PyTorch or TensorFlow: {exc}",
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

        n_sta_obs = sum(1 for _ in _iter_items(sites))
        n_sta_use = min(n_sta_obs, self.n_stations_per_profile)
        freqs_cfg = input_data.get("freqs", self.freqs)
        if freqs_cfg is None:
            freqs = _DEFAULT_FREQS_2D[: self.n_freqs]
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
        depth_max_cfg = input_data.get("depth_max", self.depth_max)
        if depth_max_cfg is not None and float(depth_max_cfg) <= 0:
            return AgentResult.failed(
                "'depth_max' must be positive when supplied.",
                elapsed=time.time() - t0,
            )

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
            feat = _z_to_features(ed, z, fr, freqs)
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
        # each observed feature is flat [log10(rho_xy) | phase_xy]
        # (length 2 * n_freqs): fold back into (n_components=2, n_freqs)
        comp_feats = [
            np.asarray(f, dtype=np.float32).reshape(2, n_freqs)
            for f in obs_feats
        ]
        X_obs = np.stack(comp_feats, axis=2)  # (n_components, n_freqs, n_sta)
        X_obs_4d = X_obs[None, ...]  # (1, n_components, n_freqs, n_sta)

        # ── triangular-mesh path: separate output shape, separate return ──────
        if self.physics == "mt2d_tri":
            return self._execute_mt2d_tri(
                station_names=station_names,
                comp_feats=comp_feats,
                freqs=freqs,
                output_dir=output_dir,
                warnings=warnings,
                t0=t0,
            )

        # ── generate synthetic 2-D training data ───────────────────────────────
        mt2d_dataset = None
        if self.physics == "mt2d":
            self._log.info(
                "Generating %d synthetic 2-D Maxwell realizations "
                "(%d stations × %d freqs)…",
                self.n_train_profiles,
                n_sta,
                n_freqs,
            )
            try:
                mt2d_dataset, X2d, y2d = _generate_mt2d_training_data(
                    n_sta=n_sta,
                    freqs=freqs,
                    n_depth=self.n_depth,
                    depth_max_m=depth_max_cfg,
                    station_spacing_m=self.station_spacing_m,
                    correlation_length_x_m=self.correlation_length_x_m,
                    correlation_length_z_m=self.correlation_length_z_m,
                    log_resistivity_mean=self.log_resistivity_mean,
                    log_resistivity_std=self.log_resistivity_std,
                    mesh_safety_factor=self.mesh_safety_factor,
                    max_mesh_cells=self.max_mesh_cells,
                    n_realizations=self.n_train_profiles,
                    verbose=self.verbose,
                )
                n_samp = len(mt2d_dataset.select("train"))
            except Exception as exc:
                return AgentResult.failed(
                    f"2-D Maxwell dataset assembly failed: {exc}",
                    elapsed=time.time() - t0,
                )
        else:
            self._log.info(
                "Generating %d synthetic 2-D profiles (%d stations × %d freqs)…",
                self.n_train_profiles,
                n_sta,
                n_freqs,
            )
            try:
                n_layers = max(3, self.n_depth // 8)
                n_1d = self.n_train_profiles * n_sta
                ds1d = generate_dataset(
                    solver="mt1d",
                    n_samples=n_1d,
                    freqs=freqs,
                    n_layers=n_layers,
                    noise_level=0.03,
                    seed=0,
                    n_jobs=1,
                    verbose=self.verbose,
                )
                # X1d: (n_1d, n_freqs * 2) flat [log10(rho_xy) | phase_xy]
                # → reshape to profiles
                X1d = ds1d.X  # (n_1d, n_freqs * 2)
                # ds1d.y is [log10(rho) (n_layers) | log10(h) (n_layers - 1)];
                # the 2-D section target is log10(rho) by depth only.
                y1d = ds1d.y[:, :n_layers]  # (n_1d, n_layers)

                n_samp = n_1d // n_sta
                X1d = X1d[: n_samp * n_sta]
                y1d = y1d[: n_samp * n_sta]

                # reshape: (n_samp, n_sta, n_freqs, 2)
                X2d_raw = X1d.reshape(n_samp, n_sta, n_freqs, 2)
                # → (n_samp, 2, n_freqs, n_sta)
                X2d = X2d_raw.transpose(0, 3, 2, 1)

                y2d_raw = y1d.reshape(n_samp, n_sta, n_layers)
                # → (n_samp, n_layers, n_sta)  then upsample to n_depth
                y2d = y2d_raw.transpose(0, 2, 1)  # (n_samp, n_layers, n_sta)
                if n_layers < self.n_depth:
                    from scipy.ndimage import zoom

                    y2d = zoom(
                        y2d, (1, self.n_depth / n_layers, 1), order=1
                    )  # (n_samp, n_depth, n_sta)

            except Exception as exc:
                return AgentResult.failed(
                    f"2-D dataset assembly failed: {exc}",
                    elapsed=time.time() - t0,
                )

        # ── train U-Net ────────────────────────────────────────────────────────
        self._log.info(
            "Training EMInverter2D (%s) for %d epochs…",
            self.arch,
            self.epochs,
        )
        try:
            inv2d = EMInverter2D(
                n_components=self.n_components,
                n_depth=self.n_depth,
                n_stations=n_sta,
                n_freqs=n_freqs,
                arch=self.arch,
            )
            inv2d.fit(
                X2d,
                y2d,
                epochs=self.epochs,
                batch_size=max(4, min(16, n_samp // 10)),
                patience=(
                    self.patience
                    if self.patience is not None
                    else max(5, self.epochs // 5)
                ),
                verbose=self.verbose,
                lambda_x=self.lambda_x,
                lambda_z=self.lambda_z,
                lambda_tv=self.lambda_tv,
            )
        except Exception as exc:
            return AgentResult.failed(
                f"U-Net training failed: {exc}",
                elapsed=time.time() - t0,
            )

        # ── predict ────────────────────────────────────────────────────────────
        try:
            pred_2d = inv2d.predict(X_obs_4d)[0]  # (n_depth, n_sta)
        except Exception as exc:
            return AgentResult.failed(
                f"2-D prediction failed: {exc}",
                elapsed=time.time() - t0,
            )

        # ── known-truth recovery check (physics="mt2d" only) ────────────────
        # The field survey has no ground truth; this instead checks the
        # trained network against held-out synthetic realizations whose
        # true resistivity is known, per the AI-inversion plan's M0
        # baseline-metrics requirement.
        mt2d_recovery: dict[str, Any] | None = None
        if mt2d_dataset is not None:
            try:
                mt2d_recovery = _mt2d_recovery_check(inv2d, mt2d_dataset)
            except Exception as exc:
                warnings.append(f"mt2d recovery check failed: {exc}")

        ths = _inv2d_thicknesses(
            self.n_depth,
            freqs,
            None if depth_max_cfg is None else float(depth_max_cfg),
        )
        depths = np.concatenate([[0.0], np.cumsum(ths)]) / 1000.0  # km

        # ── figures ───────────────────────────────────────────────────────────
        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}

        try:
            from ..ai.plot.inversion import (
                plot_inversion_result_2d,
            )

            station_pos = np.arange(n_sta, dtype=float) * 0.5  # 0.5 km spacing
            fig_inv = plot_inversion_result_2d(
                pred_2d,
                depths=depths * 1000.0,  # expects metres
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
                p = self._save_figure(
                    fig_inv,
                    output_dir,
                    "inv2d_section",
                    warnings_list=warnings,
                )
                if p:
                    fig_paths["inv2d_section"] = p
        except Exception as exc:
            warnings.append(f"plot_inversion_result_2d: {exc}")
            # fallback: simple imshow
            try:
                import matplotlib.pyplot as plt

                from ..api.station import (
                    PYCSAMT_STATION_RENDERING,
                )

                fig, ax = plt.subplots(figsize=(12, 5))
                vv = pred_2d[np.isfinite(pred_2d)]
                im = ax.imshow(
                    pred_2d,
                    aspect="auto",
                    origin="upper",
                    extent=(-0.5, n_sta - 0.5, depths[-1], depths[0]),
                    cmap="jet_r",
                    vmin=float(np.percentile(vv, 5)) if vv.size else 0,
                    vmax=float(np.percentile(vv, 95)) if vv.size else 4,
                    interpolation="bilinear",
                )
                PYCSAMT_STATION_RENDERING.apply(
                    ax,
                    np.arange(n_sta, dtype=float),
                    station_names,
                    preset="inversion",
                    xlim=(-0.5, n_sta - 0.5),
                )
                ax.set_ylabel("Depth (km)", fontsize=9)
                ax.tick_params(axis="y", labelsize=8)
                self._section.add_colorbar(
                    im, ax, label="$\\log_{10}\\rho$ (Ω·m)"
                )
                ax.set_title(
                    "2-D AI inversion (U-Net)", fontsize=10, fontweight="bold"
                )
                fig.tight_layout()
                figures["inv2d_section"] = fig
                p = self._save_figure(
                    fig, output_dir, "inv2d_section", warnings_list=warnings
                )
                if p:
                    fig_paths["inv2d_section"] = p
            except Exception:
                pass

        from ._topography import resolve_agent_topography

        topography = resolve_agent_topography(
            input_data.get("topography", False),
            sites=sites,
            station_names=station_names,
            coords_m=None,
            warnings_list=warnings,
        )
        if topography is not None and topography["applied"]:
            try:
                from ..topo import plot_topo_section

                # The shared terrain adapter expects station x depth.  The
                # U-Net prediction is depth x station, hence the transpose.
                topo_model = {
                    "pred_rho": pred_2d.T,
                    "depths_km": depths,
                    "station_names": station_names,
                    "rms_global": np.nan,
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
                    title="2-D U-Net inversion draped below station topography",
                )
                fig_topo = ax.get_figure()
                figures["topography_section"] = fig_topo
                p = self._save_figure(
                    fig_topo,
                    output_dir,
                    "inv2d_topography_section",
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
            rho_mean = float(np.nanmean(10**pred_2d))
            rho_std = float(np.nanstd(10**pred_2d))
            prompt = (
                f"2-D AI inversion (U-Net) summary:\n"
                f"  Profile: {n_sta} stations × {n_freqs} frequencies\n"
                f"  Section: {self.n_depth} depth cells, "
                f"  max depth {depths[-1]:.1f} km\n"
                f"  Mean resistivity: {rho_mean:.0f} Ω·m ± {rho_std:.0f}\n"
                f"  Warnings: {warnings[:3] if warnings else 'none'}\n\n"
                "Interpret the 2-D resistivity section."
            )
            interp = self.query_llm(prompt, max_tokens=250)

        # ── data-space RMS ────────────────────────────────────────
        rms_global = _compute_rms_2d(X_obs, pred_2d, ths, freqs)

        elapsed = time.time() - t0
        recovery_note = ""
        if mt2d_recovery is not None:
            recovery_note = (
                f" Held-out recovery RMSE={mt2d_recovery['rmse']:.3f} "
                f"(log10 Ω·m, n={mt2d_recovery['n_samples']})."
            )
        return AgentResult(
            status="success",
            summary=(
                f"2-D AI inversion (U-Net, physics={self.physics}): "
                f"{n_sta} stations x "
                f"{self.n_depth} depth cells. "
                f"RMS={rms_global:.3f}. "
                f"{len(figures)} figures."
                f"{recovery_note}"
            ),
            data={
                "pred_section": pred_2d,
                "depths_km": depths,
                "frequency_grid_hz": freqs,
                "configured_depth_max_m": depth_max_cfg,
                "station_names": station_names,
                "topography": topography,
                "station_elevation_m": (
                    None if topography is None else topography["elevation_m"]
                ),
                "station_chainage_km": (
                    None if topography is None else topography["chainage_km"]
                ),
                "rms_global": rms_global,
                "inverter": inv2d,
                "figures": figures,
                "figure_paths": fig_paths,
                "physics": self.physics,
                "mt2d_recovery": mt2d_recovery,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )

    def _execute_mt2d_tri(
        self,
        *,
        station_names: list[str],
        comp_feats: list[np.ndarray],
        freqs: np.ndarray,
        output_dir: str | None,
        warnings: list[str],
        t0: float,
    ) -> AgentResult:
        """``physics="mt2d_tri"`` branch: GCN inversion on a triangular mesh.

        Separate from the shared U-Net path above because its output is
        a per-triangle field, not a ``(depth, station)`` grid -- see the
        module docstring's ``physics="mt2d_tri"`` section for why this
        return is self-contained rather than merged into the shared
        ``AgentResult`` at the end of :meth:`execute`.
        """
        n_sta = len(station_names)
        n_freqs = len(freqs)
        station_x = np.arange(n_sta, dtype=float) * self.station_spacing_m

        self._log.info(
            "Generating %d synthetic triangular-mesh 2-D Maxwell "
            "realizations (%d stations x %d freqs)…",
            self.n_train_profiles,
            n_sta,
            n_freqs,
        )
        try:
            dataset, mesh, X_train, y_train = _generate_mt2d_tri_training_data(
                station_x_m=station_x,
                freqs=freqs,
                depth_max_m=self.depth_max,
                correlation_length_x_m=self.correlation_length_x_m,
                correlation_length_z_m=self.correlation_length_z_m,
                log_resistivity_mean=self.log_resistivity_mean,
                log_resistivity_std=self.log_resistivity_std,
                mesh_target_cell_m=self.mesh_target_cell_m,
                field_grid_cell_m=self.field_grid_cell_m,
                topo_x_m=self.topo_x_m,
                topo_z_m=self.topo_z_m,
                n_realizations=self.n_train_profiles,
                adapter=self.mare2dem_adapter,
                verbose=self.verbose,
            )
        except Exception as exc:
            return AgentResult.failed(
                f"Triangular-mesh Maxwell dataset assembly failed: {exc}",
                elapsed=time.time() - t0,
            )

        self._log.info(
            "Training GCNInverter3D (n_layers=1) for %d epochs…", self.epochs
        )
        try:
            from ..ai.inversion.inv3d import GCNInverter3D

            gcn = GCNInverter3D(
                n_features=2 * n_freqs + N_POSITION_FEATURES,
                n_layers=1,
                hidden=self.gcn_hidden,
            )
            gcn.fit(
                X_train,
                y_train,
                coords=mesh.triangle_centroids_m,
                radius=self.gcn_adjacency_radius_m,
                epochs=self.epochs,
                batch_size=max(2, min(8, len(X_train) // 5)),
                patience=(
                    self.patience
                    if self.patience is not None
                    else max(5, self.epochs // 5)
                ),
                verbose=self.verbose,
            )
        except Exception as exc:
            return AgentResult.failed(
                f"GCN training failed: {exc}", elapsed=time.time() - t0
            )

        try:
            X_obs = _maxwell_tri_samples_observed_features(
                comp_feats, mesh, station_x
            )
            # No coords=/adjacency= here: reuses the exact adjacency built
            # from this same shared mesh during fit() (self._A_stored),
            # rather than rebuilding one at predict()'s different default
            # radius (5_000.0 m) than gcn_adjacency_radius_m used to train.
            pred_log_rho = gcn.predict(X_obs[None, ...])[0, :, 0]
        except Exception as exc:
            return AgentResult.failed(
                f"Triangular-mesh prediction failed: {exc}",
                elapsed=time.time() - t0,
            )

        mt2d_tri_recovery: dict[str, Any] | None = None
        try:
            mt2d_tri_recovery = _mt2d_tri_recovery_check(
                gcn, dataset, mesh, station_x
            )
        except Exception as exc:
            warnings.append(f"mt2d_tri recovery check failed: {exc}")

        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}
        try:
            import matplotlib.pyplot as plt

            from ..api.mesh import draw_tri_mesh
            from ..api.station import StationAxisStyle, StationMarkerStyle

            fig, ax = plt.subplots(figsize=(10, 5))
            fill, _edges = draw_tri_mesh(
                ax, mesh, pred_log_rho, preset="review", cmap="jet_r"
            )
            station_z = (
                np.interp(station_x, self.topo_x_m, self.topo_z_m)
                if self.topo_x_m is not None
                else np.zeros_like(station_x)
            )
            ax.scatter(
                station_x,
                station_z,
                label="Stations",
                **StationMarkerStyle().kwargs(),
            )
            y0, y1 = ax.get_ylim()
            label_pad = 0.03 * abs(y1 - y0)
            # Thin labels the same way StationAxisStyle does for every other
            # station-axis figure: every one of the n markers stays visible,
            # but only a subset gets a text label, so this scales from a
            # ten-station line (all labelled) to a fifty-station one
            # (crowding and title collisions otherwise) without a separate
            # code path.
            visible = StationAxisStyle().label_indices(
                station_names, figwidth_in=fig.get_figwidth()
            )
            for i in visible:
                ax.annotate(
                    station_names[i],
                    (station_x[i], station_z[i] - label_pad),
                    ha="center",
                    va="bottom",
                    fontsize=7,
                    rotation=90,
                    clip_on=False,
                    zorder=6,
                )
            ax.invert_yaxis()

            # The triangular mesh's own topography-following top edge
            # already marks the surface; a plain box-frame line there
            # besides is redundant and, once rotated station labels are
            # added above it, is what they visually collide with -- the
            # same fix applied to pycsamt.topo.overlay.draw_topo_section
            # for every flat-grid topography-draped section in the
            # package. The y-axis is inverted here (depth increases
            # downward), so "the surface side" is the *shallow* end of
            # ylim, not simply ylim[1].
            ax.spines["top"].set_visible(False)
            if len(visible) > 0:
                try:
                    fig.canvas.draw()
                    renderer = fig.canvas.get_renderer()
                    inv_transform = ax.transData.inverted()
                    shallow_y = min(
                        inv_transform.transform(
                            (t.get_window_extent(renderer=renderer).x0,
                             t.get_window_extent(renderer=renderer).y1)
                        )[1]
                        for t in ax.texts
                    )
                    cur_lo, cur_hi = ax.get_ylim()
                    shallow_bound = min(cur_lo, cur_hi)
                    if np.isfinite(shallow_y) and shallow_y < shallow_bound:
                        pad = 0.02 * abs(cur_hi - cur_lo)
                        new_shallow = shallow_y - pad
                        if cur_lo <= cur_hi:
                            ax.set_ylim(new_shallow, cur_hi)
                        else:
                            ax.set_ylim(cur_lo, new_shallow)
                except Exception:
                    pass

            if ((self.depth_max is not None and self.depth_max >= 10_000.0)
                    or np.nanmax(station_x) >= 10_000.0):
                from matplotlib.ticker import FuncFormatter

                km = FuncFormatter(lambda value, _pos: f"{value / 1000:g}")
                ax.xaxis.set_major_formatter(km)
                ax.yaxis.set_major_formatter(km)
                ax.set_xlabel("Profile distance (km)")
                ax.set_ylabel("Depth below elevation datum (km)")
            else:
                ax.set_xlabel("x (m)")
                ax.set_ylabel("Depth (m)")
            ax.legend(loc="lower right", frameon=True)
            fig.colorbar(fill, ax=ax, label=r"$\log_{10}\rho$ (Ohm.m)")
            # Rotated station labels sit outside the axes box (clip_on=False)
            # above each marker, and PlotConfig's default bbox_inches="tight"
            # crops away any empty margin reserved via subplots_adjust, so
            # only a real points-based title pad (not axes-fraction margin)
            # reliably keeps the two apart across both a ten- and a
            # fifty-station line.
            ax.set_title(
                "2-D AI inversion (triangular mesh, GCN)",
                fontsize=10,
                fontweight="bold",
                pad=90,
            )
            figures["inv2d_tri_section"] = fig
            p = self._save_figure(
                fig, output_dir, "inv2d_tri_section", warnings_list=warnings
            )
            if p:
                fig_paths["inv2d_tri_section"] = p
        except Exception as exc:
            warnings.append(f"draw_tri_mesh: {exc}")

        interp: str | None = None
        if self.api_key:
            rho_mean = float(np.nanmean(10.0**pred_log_rho))
            rho_std = float(np.nanstd(10.0**pred_log_rho))
            prompt = (
                f"2-D AI inversion (triangular mesh, GCN) summary:\n"
                f"  Profile: {n_sta} stations x {n_freqs} frequencies\n"
                f"  Mesh: {mesh.n_triangles} triangles\n"
                f"  Mean resistivity: {rho_mean:.0f} Ohm.m +/- {rho_std:.0f}\n"
                f"  Warnings: {warnings[:3] if warnings else 'none'}\n\n"
                "Interpret the triangular-mesh resistivity section."
            )
            interp = self.query_llm(prompt, max_tokens=250)

        elapsed = time.time() - t0
        recovery_note = ""
        if mt2d_tri_recovery is not None:
            recovery_note = (
                f" Held-out recovery RMSE={mt2d_tri_recovery['rmse']:.3f} "
                f"(log10 Ohm.m, n={mt2d_tri_recovery['n_samples']})."
            )
        return AgentResult(
            status="success",
            summary=(
                f"2-D AI inversion (GCN, physics=mt2d_tri): {n_sta} "
                f"stations x {mesh.n_triangles} triangles. "
                f"{len(figures)} figures.{recovery_note}"
            ),
            data={
                "pred_triangles": {
                    "mesh": mesh,
                    "log10_resistivity": pred_log_rho,
                },
                "frequency_grid_hz": freqs,
                "station_names": station_names,
                "inverter": gcn,
                "figures": figures,
                "figure_paths": fig_paths,
                "physics": self.physics,
                "mt2d_tri_recovery": mt2d_tri_recovery,
                "training_history": dict(getattr(gcn, "_history", {})),
                "epochs_completed": len(
                    getattr(gcn, "_history", {}).get("train_loss", [])
                ),
                "best_validation_loss": float(
                    np.nanmin(getattr(gcn, "_history", {}).get("val_loss", [np.nan]))
                ),
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


def _generate_mt2d_training_data(
    *,
    n_sta: int,
    freqs: np.ndarray,
    n_depth: int,
    depth_max_m: float | None,
    station_spacing_m: float,
    correlation_length_x_m: tuple[float, float],
    correlation_length_z_m: tuple[float, float],
    log_resistivity_mean: float,
    log_resistivity_std: float,
    mesh_safety_factor: float,
    max_mesh_cells: int,
    n_realizations: int,
    verbose: bool | int | str = False,
):
    """Build a genuinely 2-D correlated training set and convert it
    to :class:`~pycsamt.ai.inversion.inv2d.EMInverter2D`'s tensors.

    Uses :func:`~pycsamt.ai.training.dataset2d.\
generate_2d_maxwell_dataset`
    on a synthetic grid whose ``nx`` stations sit exactly
    ``station_spacing_m`` apart (the real survey's actual station
    geometry is not read at this stage) and whose ``nz`` matches
    ``n_depth`` exactly, so no resize is needed to align the label
    with :meth:`EMInverter2D.predict`'s output shape.

    Returns
    -------
    dataset : Maxwell2DDataset
        Full dataset (train/validation/test), kept for a later
        known-truth recovery check.
    X : ndarray (n_train, 2, n_freqs, n_stations)
    y : ndarray (n_train, n_depth, n_stations)
    """
    from ..ai.geology import GeologyGrid
    from ..ai.training.dataset2d import (
        Maxwell2DDatasetConfig,
        generate_2d_maxwell_dataset,
    )

    depth_total_m = (
        float(depth_max_m)
        if depth_max_m is not None
        else float(np.sum(_default_thicknesses(n_depth, freqs)))
    )
    grid = GeologyGrid.regular_2d(
        nx=n_sta,
        nz=n_depth,
        dx_m=station_spacing_m,
        dz_m=depth_total_m / n_depth,
    )
    config = Maxwell2DDatasetConfig(
        dataset_id="inv2dagent-mt2d",
        grid=grid,
        correlation_length_x_m=correlation_length_x_m,
        correlation_length_z_m=correlation_length_z_m,
        frequencies_hz=freqs,
        station_x_m=grid.x_m,
        n_realizations=n_realizations,
        seed=0,
        log_resistivity_mean=log_resistivity_mean,
        log_resistivity_std=log_resistivity_std,
        components=("zxy",),
        mesh_safety_factor=mesh_safety_factor,
        max_mesh_cells=max_mesh_cells,
        validation_fraction=0.1,
        test_fraction=0.1,
        verbose=verbose,
    )
    dataset = generate_2d_maxwell_dataset(config)
    X, y = _maxwell_samples_to_unet_arrays(dataset.select("train"), freqs)
    return dataset, X, y


def _maxwell_samples_to_unet_arrays(
    samples,
    freqs: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert ``Maxwell2DSample`` objects into ``EMInverter2D``'s
    ``(X, y)`` convention.

    Channels 0/1 of ``X`` are log10(apparent resistivity)/phase
    (degrees) for ``zxy``, matching :func:`_z_to_features`'s feature
    layout so synthetic training data and the observed pseudosection
    share one convention.
    """
    if not samples:
        raise ValueError("samples must contain at least one entry.")
    n_freq = len(freqs)
    n_depth, n_sta = samples[0].resistivity_ohm_m.shape
    X = np.empty((len(samples), 2, n_freq, n_sta), dtype=np.float32)
    y = np.empty((len(samples), n_depth, n_sta), dtype=np.float32)
    for i, sample in enumerate(samples):
        zxy = sample.survey.impedance[:, :, 0]  # (n_sta, n_freq)
        f = sample.survey.frequencies_hz[None, :]
        mu0 = 4.0e-7 * np.pi
        rho_a = np.abs(zxy) ** 2 / (2.0 * np.pi * f * mu0)
        phase = np.degrees(np.angle(zxy))
        X[i, 0] = np.log10(np.clip(rho_a, 1e-12, None)).T
        X[i, 1] = phase.T
        y[i] = np.log10(sample.resistivity_ohm_m)
    return X, y


def _mt2d_recovery_check(inv2d, dataset) -> dict[str, Any] | None:
    """Evaluate known-truth recovery on a held-out dataset partition.

    Prefers the ``test`` partition, falling back to ``validation``;
    returns ``None`` when neither has samples (e.g. too few
    realizations were requested to split three ways).
    """
    from ..ai.validation import recovery_report

    held_out = dataset.select("test") or dataset.select("validation")
    if not held_out:
        return None
    freqs = held_out[0].survey.frequencies_hz
    rmse_vals, mae_vals, r2_vals = [], [], []
    for sample in held_out:
        X_i, _ = _maxwell_samples_to_unet_arrays([sample], freqs)
        pred_log_rho = inv2d.predict(X_i)[0]  # (n_depth, n_sta)
        true_log_rho = np.log10(sample.resistivity_ohm_m)
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


def _generate_mt2d_tri_training_data(
    *,
    station_x_m: np.ndarray,
    freqs: np.ndarray,
    depth_max_m: float | None,
    correlation_length_x_m: tuple[float, float],
    correlation_length_z_m: tuple[float, float],
    log_resistivity_mean: float,
    log_resistivity_std: float,
    mesh_target_cell_m: float,
    field_grid_cell_m: float,
    n_realizations: int,
    adapter: Any | None,
    topo_x_m: Any = None,
    topo_z_m: Any = None,
    verbose: bool | int | str = False,
):
    """Build a triangular-mesh training set and convert it to
    :class:`~pycsamt.ai.inversion.inv3d.GCNInverter3D`'s ``(X, y)``
    convention.

    Returns
    -------
    dataset : MaxwellTri2DDataset
        Full dataset (train/validation/test), kept for a later
        known-truth recovery check.
    mesh : TriMesh
        Shared triangular geometry (same object dataset.mesh).
    X : ndarray (n_train, n_triangles, 2 * n_freqs + N_POSITION_FEATURES)
    y : ndarray (n_train, n_triangles, 1)
    """
    from ..ai.training.dataset2d_tri import (
        MaxwellTri2DDatasetConfig,
        generate_2d_tri_maxwell_dataset,
    )

    n_sta = len(station_x_m)
    depth_total_m = (
        float(depth_max_m)
        if depth_max_m is not None
        else float(np.sum(_default_thicknesses(max(n_sta, 8), freqs)))
    )
    x_max = float(station_x_m[-1]) if n_sta > 1 else float(station_x_m[0]) + 1.0
    config = MaxwellTri2DDatasetConfig(
        dataset_id="inv2dagent-mt2dtri",
        x_range_m=(0.0, x_max),
        z_range_m=(0.0, depth_total_m),
        correlation_length_x_m=correlation_length_x_m,
        correlation_length_z_m=correlation_length_z_m,
        frequencies_hz=freqs,
        station_x_m=station_x_m,
        n_realizations=n_realizations,
        seed=0,
        log_resistivity_mean=log_resistivity_mean,
        log_resistivity_std=log_resistivity_std,
        components=("zxy",),
        mesh_target_cell_m=mesh_target_cell_m,
        field_grid_cell_m=field_grid_cell_m,
        topo_x_m=topo_x_m,
        topo_z_m=topo_z_m,
        validation_fraction=0.1,
        test_fraction=0.1,
        verbose=verbose,
    )
    dataset = generate_2d_tri_maxwell_dataset(config, adapter=adapter)
    X, y = _maxwell_tri_samples_to_gcn_arrays(
        dataset.select("train"), dataset.mesh, station_x_m
    )
    return dataset, dataset.mesh, X, y


def _nearest_station_features(
    comp_feats: list[np.ndarray] | np.ndarray,
    mesh: Any,
    station_x_m: np.ndarray,
) -> np.ndarray:
    """Broadcast each triangle's nearest station's ``(2, n_freq)`` feature
    block onto that triangle, flattened to ``(n_triangles, 2 * n_freq)``.

    Nearest is by x-distance only, since receivers sit on the mesh's own
    ``z=0`` surface and every triangle's centroid has a well-defined x.

    Every triangle assigned to the same nearest station receives an
    *identical* feature block here, regardless of its own depth or exact
    position -- this function alone gives a GCN nothing to key a
    depth-varying prediction on. See :func:`_triangle_position_features`
    for the per-triangle signal that fixes that.
    """
    station_x = np.asarray(station_x_m, dtype=float)
    centroid_x = mesh.triangle_centroids_m[:, 0]
    nearest = np.argmin(
        np.abs(centroid_x[:, None] - station_x[None, :]), axis=1
    )
    flat = np.stack(
        [np.asarray(f, dtype=np.float32).reshape(-1) for f in comp_feats],
        axis=0,
    )  # (n_sta, 2 * n_freq)
    return flat[nearest]


def _triangle_position_features(mesh: Any) -> np.ndarray:
    """Return each triangle's own ``[x_norm, z_norm]`` position, both
    normalized to the mesh's own bounding box (``[0, 1]``).

    :func:`_nearest_station_features` broadcasts the *same* sounding
    curve onto every triangle nearest a given station, so two triangles
    at very different depths under the same station otherwise look
    identical to the network. Appending each triangle's own normalized
    position gives it something to actually vary a depth-dependent
    prediction on, rather than relying solely on message-passing
    between differently-stationed neighbours to break the symmetry.
    """
    centroids = mesh.triangle_centroids_m
    lo = centroids.min(axis=0)
    span = np.maximum(centroids.max(axis=0) - lo, 1.0)
    return ((centroids - lo) / span).astype(np.float32)


N_POSITION_FEATURES = 2


def _maxwell_tri_samples_observed_features(
    comp_feats: list[np.ndarray],
    mesh: Any,
    station_x_m: np.ndarray,
) -> np.ndarray:
    """Return ``(n_triangles, 2 * n_freq + 2)`` observed features for
    prediction: the nearest station's sounding curve plus the triangle's
    own normalized ``[x, z]`` position.
    """
    nearest_feats = _nearest_station_features(comp_feats, mesh, station_x_m)
    return np.concatenate(
        [nearest_feats, _triangle_position_features(mesh)], axis=1
    )


def _maxwell_tri_samples_to_gcn_arrays(
    samples,
    mesh: Any,
    station_x_m: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert ``MaxwellTri2DSample`` objects into
    :class:`~pycsamt.ai.inversion.inv3d.GCNInverter3D`'s ``(X, y)``
    convention: per-triangle nearest-station ``[log10(rho_a), phase]``
    features for ``zxy`` concatenated with the triangle's own normalized
    ``[x, z]`` position, and per-triangle ``log10(resistivity)`` targets.
    """
    if not samples:
        raise ValueError("samples must contain at least one entry.")
    n_tri = mesh.n_triangles
    n_freq = len(samples[0].survey.frequencies_hz)
    n_feat = 2 * n_freq + N_POSITION_FEATURES
    pos_feats = _triangle_position_features(mesh)
    X = np.empty((len(samples), n_tri, n_feat), dtype=np.float32)
    y = np.empty((len(samples), n_tri, 1), dtype=np.float32)
    for i, sample in enumerate(samples):
        comp_idx = sample.survey.components.index("zxy")
        zxy = sample.survey.impedance[:, :, comp_idx]  # (n_sta, n_freq)
        f = sample.survey.frequencies_hz[None, :]
        mu0 = 4.0e-7 * np.pi
        rho_a = np.abs(zxy) ** 2 / (2.0 * np.pi * f * mu0)
        phase = np.degrees(np.angle(zxy))
        log_rho_a = np.log10(np.clip(rho_a, 1e-12, None))
        per_station = np.concatenate([log_rho_a, phase], axis=1)  # (n_sta, 2*n_freq)
        nearest_feats = _nearest_station_features(per_station, mesh, station_x_m)
        X[i] = np.concatenate([nearest_feats, pos_feats], axis=1)
        y[i, :, 0] = np.log10(sample.resistivity_ohm_m)
    return X, y


def _mt2d_tri_recovery_check(gcn, dataset, mesh, station_x_m) -> dict[str, Any] | None:
    """Evaluate known-truth recovery on a held-out dataset partition.

    Prefers the ``test`` partition, falling back to ``validation``;
    returns ``None`` when neither has samples.
    """
    from ..ai.validation import recovery_report

    held_out = dataset.select("test") or dataset.select("validation")
    if not held_out:
        return None
    rmse_vals, mae_vals, r2_vals = [], [], []
    for sample in held_out:
        X_i, _ = _maxwell_tri_samples_to_gcn_arrays([sample], mesh, station_x_m)
        pred = gcn.predict(X_i)[0, :, 0]  # reuses self._A_stored, see above
        true = np.log10(sample.resistivity_ohm_m)
        # recovery_report requires a 2-D/3-D grid; a triangle mesh has no
        # such shape, so reshape to a (n_triangles, 1) column -- degenerate
        # but valid, and the per-cell metrics (rmse/mae/r2) are unaffected.
        report = recovery_report(
            pred.reshape(-1, 1), true.reshape(-1, 1), compute_ssim=False
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


def _inv2d_thicknesses(
    n_depth: int,
    freqs: np.ndarray,
    depth_max: float | None,
) -> np.ndarray:
    """Return an n-depth display grid with an optional explicit bottom."""
    if depth_max is None:
        return _default_thicknesses(n_depth, freqs)
    n_intervals = max(int(n_depth) - 1, 1)
    weights = np.geomspace(1.0, 3.0, n_intervals)
    return (float(depth_max) * weights / weights.sum()).astype(float)


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

        depths_m = np.concatenate([[0.0], np.cumsum(thicknesses[:n_depth])])
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
                d_b = 503.0 * np.sqrt(max(rho_a_obs, 1.0) * T)
                idx = int(np.searchsorted(depths_m, d_b))
                idx = min(idx, n_depth - 1)
                pred_log_rho_a[fi] = profile[idx]
            diff = pred_log_rho_a - obs_log_rho[:, s]
            finite = np.isfinite(diff)
            if finite.any():
                rms_vals.append(float(np.sqrt(np.mean(diff[finite] ** 2))))

        return float(np.mean(rms_vals)) if rms_vals else np.nan

    except Exception:
        return np.nan


__all__ = ["Inv2DAgent"]
