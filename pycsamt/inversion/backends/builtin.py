# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Built-in dependency-light EM inversion backend.

The built-in backend provides a small SciPy-based inversion engine for
routine smoke tests, demos, and lightweight studies where optional
physics packages such as SimPEG or pyGIMLi are not available. It supports
natural-source 1-D soundings, stitched 2-D profiles, finite-difference
2-D MT profile inversion, and 1-D/stitched 2-D TDEM runs.

The backend parameterizes layered soundings in log10 resistivity and
log10 thickness. Profile mode either runs each station independently and
stitches the recovered columns, or, when ``profile_mode="fd2d"`` is
selected, solves a finite-difference 2-D natural-source problem against
the local :mod:`pycsamt.forward` solver.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from ..base import BaseInversionBackend
from ..mesh import (
    InversionMesh,
    build_fd2d_grid,
    core_rho_from_start,
)
from ..model import StartingModel
from ..objective import (
    component_errors,
    component_mask,
    weighted_rms,
)
from ..regularization import (
    regularization_from_config,
    regularization_residual,
    regularization_weight,
)
from ..results import (
    InversionHistory,
    InversionResult,
    InversionUncertainty,
)

__all__ = ["Builtin1DBackend"]


class Builtin1DBackend(BaseInversionBackend):
    name = "builtin"
    supports = (
        ("mt", "1d"),
        ("mt", "2d"),
        ("amt", "1d"),
        ("amt", "2d"),
        ("csamt", "1d"),
        ("csamt", "2d"),
        ("emap", "2d"),
        ("tdem", "1d"),
        ("tdem", "2d"),
    )

    def run(self, data: Any | None = None) -> InversionResult:
        """Run the configured built-in inversion.

        Parameters
        ----------
        data : mapping, object, sequence, or path-like, optional
            Optional data override for this run. When omitted, the backend
            uses ``self.config.data``. Values are coerced through
            :class:`pycsamt.inversion.data.EMData` by the base backend.

        Returns
        -------
        InversionResult
            Backend-neutral result containing the recovered model, predicted
            response, RMS, objective value, uncertainty proxy, convergence
            history, and backend metadata.

        Raises
        ------
        ValueError
            If the configured method does not receive the required data
            components. Natural-source runs require frequencies plus apparent
            resistivity and/or phase. TDEM runs require times plus values.
        ImportError
            If SciPy is unavailable when a runnable least-squares inversion is
            requested.
        NotImplementedError
            If the configured method/dimension pair is not listed in
            :attr:`supports`.

        Examples
        --------
        >>> from pycsamt.inversion.backends.builtin import Builtin1DBackend
        >>> from pycsamt.inversion.config import InversionConfig
        >>> cfg = InversionConfig(
        ...     method="mt",
        ...     dimension="1d",
        ...     backend="builtin",
        ...     data={"freqs": [1.0, 10.0],
        ...           "rho_a": [100.0, 120.0],
        ...           "phase": [45.0, 47.0]},
        ...     max_iter=4,
        ... )
        >>> result = Builtin1DBackend(cfg).run()  # doctest: +SKIP
        """
        self.check_supported()
        cfg = self.config
        em_data = self.prepare_data(data)
        if cfg.method == "tdem":
            if not em_data.has_tdem_response:
                raise ValueError("builtin TDEM backend requires times plus values.")
        elif not em_data.has_mt_response:
            raise ValueError(
                "builtin backend requires frequencies plus rho_a and/or phase."
            )
        if cfg.dimension == "2d" and cfg.method != "tdem":
            mode = str(cfg.backend_options.get("profile_mode", "")).lower()
            if mode in {"fd2d", "2d", "physics", "finite_difference"}:
                return self._run_mt_2d_physics(em_data)
        if cfg.dimension == "2d":
            return self._run_profile(em_data)
        return self._run_sounding(em_data, station_index=None)

    def _run_profile(self, em_data) -> InversionResult:
        cfg = self.config
        n_st = em_data.n_stations
        station_names = _station_names(em_data, n_st)
        station_x = _station_x(em_data, n_st)
        columns: list[np.ndarray] = []
        station_results: list[InversionResult] = []
        used_indices: list[int] = []
        warnings: list[str] = []
        z_centers = None

        for idx in range(n_st):
            sounding = _station_data(em_data, idx)
            try:
                result = self._run_sounding(sounding, station_index=idx)
            except Exception as exc:
                warnings.append(f"{station_names[idx]}: inversion failed: {exc}")
                continue
            station_results.append(result)
            used_indices.append(idx)
            columns.append(np.log10(result.model.resistivities))
            if z_centers is None and result.mesh is not None:
                z_centers = np.asarray(result.mesh.z_centers, dtype=float)

        if not columns:
            raise RuntimeError("all station inversions failed.")

        rho_2d = np.stack(columns, axis=1)
        if z_centers is None:
            z_centers = np.arange(rho_2d.shape[0], dtype=float)
        used_names = [station_names[i] for i in used_indices]
        used_x = station_x[used_indices]
        model = {
            "rho_2d": rho_2d,
            "x_centers": used_x,
            "z_centers": z_centers,
            "station_x": used_x,
            "station_names": used_names,
        }
        mesh = InversionMesh(
            dimension="2d",
            x_centers=used_x,
            z_centers=z_centers,
        )
        rms_values = np.asarray([r.rms for r in station_results], dtype=float)
        rms = float(np.nanmean(rms_values)) if rms_values.size else float("nan")
        uncertainty = _profile_uncertainty(station_results)
        history = _profile_history(station_results)
        status = "success" if not warnings else "needs_review"
        return InversionResult(
            method=cfg.method,
            dimension=cfg.dimension,
            backend=self.name,
            status=status,
            model=model,
            mesh=mesh,
            data=em_data,
            predicted=[r.predicted for r in station_results],
            rms=rms,
            objective=float(np.nansum([r.objective for r in station_results])),
            n_iter=int(np.sum([r.n_iter for r in station_results])),
            workdir=cfg.workdir,
            native=station_results,
            uncertainty=uncertainty,
            history=history,
            warnings=warnings,
            metadata={
                **cfg.metadata,
                "profile_mode": "stitched_station_1d",
                "station_rms": rms_values.tolist(),
            },
        )

    def _run_sounding(
        self,
        em_data,
        *,
        station_index: int | None,
    ) -> InversionResult:
        cfg = self.config

        if cfg.method == "tdem":
            return self._run_tdem_sounding(em_data, station_index=station_index)
        return self._run_mt_sounding(em_data, station_index=station_index)

    def _run_mt_2d_physics(self, em_data) -> InversionResult:
        """Run a finite-difference 2-D MT inversion using the forward solver."""
        cfg = self.config

        try:
            from scipy.optimize import least_squares
        except ImportError as exc:
            raise ImportError("builtin 2-D inversion requires scipy.optimize.") from exc

        from ...forward import (
            Grid2D,
            MT2DForward,
            make_padding,
        )

        opts = cfg.backend_options
        start = StartingModel.coerce(cfg.starting_model, n_layers=cfg.n_layers)
        station_x = _station_x(em_data, em_data.n_stations)
        station_names = _station_names(em_data, em_data.n_stations)
        base_grid, core_shape = build_fd2d_grid(
            start,
            station_x,
            opts,
            Grid2D,
            make_padding_func=make_padding,
        )
        nz_core, nx_core = core_shape
        x0 = np.log10(core_rho_from_start(start, nz_core, nx_core))
        lower, upper = _rho_grid_bounds(cfg.bounds, x0.size)
        observed, errors, components = _pack_2d_observed(em_data, cfg, opts)
        reg = regularization_from_config(cfg)
        reg_weight = regularization_weight(cfg, default=0.0)
        ref = np.log10(core_rho_from_start(start, nz_core, nx_core))
        history = _HistoryRecorder(reg_weight=reg_weight)

        def residual(params: np.ndarray) -> np.ndarray:
            grid = _grid_with_core_model(base_grid, params.reshape(nz_core, nx_core))
            response = MT2DForward(
                np.asarray(em_data.frequencies, dtype=float),
                grid,
                verbose=bool(opts.get("forward_verbose", False)),
            ).run()
            predicted = _pack_2d_predicted(response, components)
            data_res = (predicted - observed) / errors
            parts = [data_res]
            reg_res = np.array([], dtype=float)
            reg_res = regularization_residual(
                params.reshape(nz_core, nx_core),
                reference=ref,
                regularization=reg,
                blocky_eps=float(opts.get("blocky_eps", 1e-2)),
                axes=("x", "z"),
            )
            if reg_res.size and reg_weight > 0.0:
                parts.append(np.sqrt(reg_weight) * reg_res)
            residuals = np.concatenate(parts)
            history.add(data_res, reg_res, params)
            return residuals

        opt = least_squares(
            residual,
            x0.reshape(-1),
            bounds=(lower, upper),
            max_nfev=cfg.max_iter,
            xtol=cfg.tol,
            ftol=cfg.tol,
            gtol=cfg.tol,
        )

        recovered_log = opt.x.reshape(nz_core, nx_core)
        recovered_rho = 10.0 ** recovered_log
        recovered_grid = _grid_with_core_model(base_grid, recovered_log)
        response = MT2DForward(
            np.asarray(em_data.frequencies, dtype=float),
            recovered_grid,
            verbose=bool(opts.get("forward_verbose", False)),
        ).run()
        predicted = _pack_2d_predicted(response, components)
        rms = weighted_rms(observed, predicted, errors)
        uncertainty = _uncertainty_from_lsq(
            opt,
            shape=(nz_core, nx_core),
            x_centers=None,
            z_centers=None,
            metadata={"parameter": "log10_rho", "mode": "fd2d"},
        )
        core_z = recovered_grid.z_centers[:nz_core]
        pad = int(recovered_grid.n_pad)
        core_x = recovered_grid.x_centers[pad: pad + nx_core]
        x_offset = float(getattr(recovered_grid, "_pycsamt_x_offset", 0.0))
        core_x = core_x - x_offset
        mesh = InversionMesh(
            dimension="2d",
            x_centers=core_x,
            z_centers=core_z,
            native=recovered_grid,
            metadata={
                "engine": "builtin_fd2d",
                "components": components,
                "grid_shape": [int(nz_core), int(nx_core)],
                "profile_mode": "fd2d",
            },
        )
        status = "converged" if opt.success else "needs_review"
        warnings = [] if opt.success else [str(opt.message)]
        return InversionResult(
            method=cfg.method,
            dimension="2d",
            backend=self.name,
            status=status,
            model={
                "rho_2d": np.log10(recovered_rho),
                "x_centers": core_x,
                "z_centers": core_z,
                "station_x": station_x,
                "station_names": station_names,
                "rho_2d_linear": recovered_rho,
            },
            mesh=mesh,
            data=em_data,
            predicted=response,
            rms=rms,
            objective=float(np.sum(opt.fun ** 2)),
            n_iter=int(opt.nfev),
            workdir=cfg.workdir,
            native=opt,
            uncertainty=uncertainty,
            history=history.to_history(metadata={"mode": "fd2d"}),
            warnings=warnings,
            metadata={
                **cfg.metadata,
                "profile_mode": "fd2d",
                "optimizer_message": str(opt.message),
                "starting_model": start.to_dict(),
                "components": components,
                "regularization_weight": reg_weight,
                "regularization": reg.to_dict(),
            },
        )

    def _run_mt_sounding(
        self,
        em_data,
        *,
        station_index: int | None,
    ) -> InversionResult:
        cfg = self.config

        try:
            from scipy.optimize import least_squares
        except ImportError as exc:
            raise ImportError("builtin inversion requires scipy.optimize.") from exc

        from ...forward import MT1DForward

        start = StartingModel.coerce(cfg.starting_model, n_layers=cfg.n_layers)
        n_layers = start.n_layers
        x0 = _pack(start)
        lower, upper = _bounds(cfg.bounds, n_layers)
        reg = regularization_from_config(cfg)
        reg_weight = regularization_weight(cfg, default=0.0)
        ref_log_rho = _reference_log10_rho(cfg, start)
        history = _HistoryRecorder(reg_weight=reg_weight)

        freqs = np.asarray(em_data.frequencies, dtype=float)
        obs_parts: list[np.ndarray] = []
        err_parts: list[np.ndarray] = []

        if em_data.rho_a is not None:
            rho_obs = np.log10(np.maximum(em_data.rho_a, 1e-12))
            obs_parts.append(rho_obs)
            err_parts.append(component_errors(
                em_data.rho_a,
                cfg,
                component="rho",
                explicit=em_data.errors,
                relative=True,
            ))
        if cfg.include_phase and em_data.phase is not None:
            obs_parts.append(np.asarray(em_data.phase, dtype=float))
            err_parts.append(component_errors(em_data.phase, cfg, component="phase"))

        observed = np.concatenate(obs_parts)
        errors = np.concatenate(err_parts)
        fwd = MT1DForward(freqs=freqs)

        def residual(params: np.ndarray) -> np.ndarray:
            model = _unpack(params, n_layers)
            response = fwd.run(model.to_layered_model())
            pred_parts: list[np.ndarray] = []
            if em_data.rho_a is not None:
                pred_parts.append(np.log10(np.maximum(response.rho_a, 1e-12)))
            if cfg.include_phase and em_data.phase is not None:
                pred_parts.append(np.asarray(response.phase, dtype=float))
            predicted = np.concatenate(pred_parts)
            data_res = (predicted - observed) / errors
            parts = [data_res]
            reg_res = np.array([], dtype=float)
            reg_res = regularization_residual(
                params[:n_layers],
                reference=ref_log_rho,
                regularization=reg,
                blocky_eps=float(cfg.backend_options.get("blocky_eps", 1e-2)),
                axes=("z",),
            )
            if reg_res.size and reg_weight > 0.0:
                parts.append(np.sqrt(reg_weight) * reg_res)
            residuals = np.concatenate(parts)
            history.add(data_res, reg_res, params)
            return residuals

        opt = least_squares(
            residual,
            x0,
            bounds=(lower, upper),
            max_nfev=cfg.max_iter,
            xtol=cfg.tol,
            ftol=cfg.tol,
            gtol=cfg.tol,
        )

        recovered = _unpack(opt.x, n_layers)
        response = fwd.run(recovered.to_layered_model())
        pred_parts = []
        if em_data.rho_a is not None:
            pred_parts.append(np.log10(np.maximum(response.rho_a, 1e-12)))
        if cfg.include_phase and em_data.phase is not None:
            pred_parts.append(np.asarray(response.phase, dtype=float))
        predicted = np.concatenate(pred_parts)
        rms = weighted_rms(observed, predicted, errors)
        z_centers = _layer_centers(recovered.thicknesses)
        uncertainty = _uncertainty_from_lsq(
            opt,
            shape=(n_layers, 1),
            parameter_indices=np.arange(n_layers),
            x_centers=np.array([0.0]),
            z_centers=z_centers,
            metadata={"parameter": "log10_rho", "mode": "mt_1d"},
        )

        mesh = InversionMesh.for_1d(z_centers)
        status = "converged" if opt.success else "needs_review"
        warnings = [] if opt.success else [str(opt.message)]
        return InversionResult(
            method=cfg.method,
            dimension="1d",
            backend=self.name,
            status=status,
            model=recovered,
            mesh=mesh,
            data=em_data,
            predicted=response,
            rms=rms,
            objective=float(np.sum(opt.fun ** 2)),
            n_iter=int(opt.nfev),
            workdir=cfg.workdir,
            native=opt,
            uncertainty=uncertainty,
            history=history.to_history(metadata={"mode": "mt_1d"}),
            warnings=warnings,
            metadata={
                **cfg.metadata,
                "optimizer_message": str(opt.message),
                "starting_model": start.to_dict(),
                "station_index": station_index,
                "regularization": reg.to_dict(),
                "regularization_weight": reg_weight,
            },
        )

    def _run_tdem_sounding(
        self,
        em_data,
        *,
        station_index: int | None,
    ) -> InversionResult:
        cfg = self.config

        try:
            from scipy.optimize import least_squares
        except ImportError as exc:
            raise ImportError("builtin inversion requires scipy.optimize.") from exc

        from ...forward import TEM1DForward

        start = StartingModel.coerce(cfg.starting_model, n_layers=cfg.n_layers)
        n_layers = start.n_layers
        x0 = _pack(start)
        lower, upper = _bounds(cfg.bounds, n_layers)
        reg = regularization_from_config(cfg)
        reg_weight = regularization_weight(cfg, default=0.0)
        ref_log_rho = _reference_log10_rho(cfg, start)
        history = _HistoryRecorder(reg_weight=reg_weight)
        times = np.asarray(em_data.times, dtype=float)
        observed = _log_abs(np.asarray(em_data.values, dtype=float))
        errors = component_errors(
            em_data.values,
            cfg,
            component="tdem",
            explicit=em_data.errors,
            relative=True,
        )

        opts = cfg.backend_options
        fwd = TEM1DForward(
            times=times,
            loop_radius=float(opts.get("loop_radius", 50.0)),
            moment=float(opts.get("moment", 1.0)),
            n_freqs=int(opts.get("n_freqs", 32)),
            n_lam=int(opts.get("n_lam", 48)),
        )

        def residual(params: np.ndarray) -> np.ndarray:
            model = _unpack(params, n_layers)
            response = fwd.run(model.to_layered_model())
            predicted = _log_abs(response.dBz_dt)
            data_res = (predicted - observed) / errors
            parts = [data_res]
            reg_res = np.array([], dtype=float)
            reg_res = regularization_residual(
                params[:n_layers],
                reference=ref_log_rho,
                regularization=reg,
                blocky_eps=float(cfg.backend_options.get("blocky_eps", 1e-2)),
                axes=("z",),
            )
            if reg_res.size and reg_weight > 0.0:
                parts.append(np.sqrt(reg_weight) * reg_res)
            residuals = np.concatenate(parts)
            history.add(data_res, reg_res, params)
            return residuals

        opt = least_squares(
            residual,
            x0,
            bounds=(lower, upper),
            max_nfev=cfg.max_iter,
            xtol=cfg.tol,
            ftol=cfg.tol,
            gtol=cfg.tol,
        )

        recovered = _unpack(opt.x, n_layers)
        response = fwd.run(recovered.to_layered_model())
        predicted = _log_abs(response.dBz_dt)
        rms = weighted_rms(observed, predicted, errors)
        z_centers = _layer_centers(recovered.thicknesses)
        uncertainty = _uncertainty_from_lsq(
            opt,
            shape=(n_layers, 1),
            parameter_indices=np.arange(n_layers),
            x_centers=np.array([0.0]),
            z_centers=z_centers,
            metadata={"parameter": "log10_rho", "mode": "tdem_1d"},
        )

        mesh = InversionMesh.for_1d(z_centers)
        status = "converged" if opt.success else "needs_review"
        warnings = [] if opt.success else [str(opt.message)]
        return InversionResult(
            method=cfg.method,
            dimension="1d",
            backend=self.name,
            status=status,
            model=recovered,
            mesh=mesh,
            data=em_data,
            predicted=response,
            rms=rms,
            objective=float(np.sum(opt.fun ** 2)),
            n_iter=int(opt.nfev),
            workdir=cfg.workdir,
            native=opt,
            uncertainty=uncertainty,
            history=history.to_history(metadata={"mode": "tdem_1d"}),
            warnings=warnings,
            metadata={
                **cfg.metadata,
                "optimizer_message": str(opt.message),
                "starting_model": start.to_dict(),
                "station_index": station_index,
                "regularization": reg.to_dict(),
                "regularization_weight": reg_weight,
                "tdem_options": {
                    "loop_radius": fwd.loop_radius,
                    "moment": fwd.moment,
                    "n_freqs": fwd.n_freqs,
                    "n_lam": fwd.n_lam,
                },
            },
        )


def _pack(model: StartingModel) -> np.ndarray:
    return np.r_[
        np.log10(model.resistivities),
        np.log10(model.thicknesses),
    ]


def _unpack(params: np.ndarray, n_layers: int) -> StartingModel:
    params = np.asarray(params, dtype=float)
    resistivities = 10.0 ** params[:n_layers]
    thicknesses = 10.0 ** params[n_layers:]
    return StartingModel(resistivities, thicknesses, name="builtin_1d")


def _reference_log10_rho(cfg: Any, start: StartingModel) -> np.ndarray:
    if cfg.reference_model is None:
        return np.log10(start.resistivities)
    reference = StartingModel.coerce(cfg.reference_model, n_layers=start.n_layers)
    if reference.n_layers != start.n_layers:
        reference = start
    return np.log10(reference.resistivities)


def _uncertainty_from_lsq(
    opt: Any,
    *,
    shape: tuple[int, int],
    parameter_indices: np.ndarray | None = None,
    x_centers: np.ndarray | None = None,
    z_centers: np.ndarray | None = None,
    metadata: dict[str, Any] | None = None,
) -> InversionUncertainty | None:
    jac = getattr(opt, "jac", None)
    fun = getattr(opt, "fun", None)
    if jac is None or fun is None:
        return None
    jac = np.asarray(jac, dtype=float)
    fun = np.asarray(fun, dtype=float)
    if jac.ndim != 2 or jac.size == 0:
        return None
    if parameter_indices is None:
        parameter_indices = np.arange(jac.shape[1])
    parameter_indices = np.asarray(parameter_indices, dtype=int)
    j_sel = jac[:, parameter_indices]
    n_params = int(parameter_indices.size)
    dof = max(int(fun.size - n_params), 1)
    variance = float(np.sum(fun ** 2) / dof)
    try:
        cov = np.linalg.pinv(j_sel.T @ j_sel) * variance
        cov_diag = np.maximum(np.diag(cov), 0.0)
    except Exception:
        cov_diag = np.full(n_params, np.nan, dtype=float)
    sensitivity = np.sqrt(np.sum(j_sel ** 2, axis=0))
    model_std = np.sqrt(cov_diag)
    confidence = _confidence_from_std(model_std, sensitivity)
    model_std_map = model_std.reshape(shape)
    sensitivity_map = sensitivity.reshape(shape)
    confidence_map = confidence.reshape(shape)
    return InversionUncertainty(
        model_std=model_std_map,
        covariance_diag=cov_diag.reshape(shape),
        sensitivity=sensitivity_map,
        confidence=confidence_map,
        station_confidence=np.nanmean(confidence_map, axis=0),
        depth_confidence=np.nanmean(confidence_map, axis=1),
        metadata={
            "type": "least_squares_jacobian_proxy",
            "dof": dof,
            "variance": variance,
            "x_centers": None if x_centers is None else np.asarray(x_centers, dtype=float).tolist(),
            "z_centers": None if z_centers is None else np.asarray(z_centers, dtype=float).tolist(),
            **dict(metadata or {}),
        },
    )


def _confidence_from_std(std: np.ndarray, sensitivity: np.ndarray) -> np.ndarray:
    std = np.asarray(std, dtype=float)
    sensitivity = np.asarray(sensitivity, dtype=float)
    if not np.any(np.isfinite(std)):
        std_score = np.zeros_like(std, dtype=float)
    else:
        finite_std = std[np.isfinite(std)]
        scale = float(np.nanmedian(finite_std)) if finite_std.size else 1.0
        scale = max(scale, 1e-12)
        std_score = 1.0 / (1.0 + np.nan_to_num(std, nan=scale) / scale)
    if not np.any(np.isfinite(sensitivity)) or np.nanmax(sensitivity) <= 0:
        sens_score = np.ones_like(sensitivity, dtype=float)
    else:
        sens_score = sensitivity / np.nanmax(sensitivity)
    confidence = np.sqrt(np.clip(std_score, 0.0, 1.0) * np.clip(sens_score, 0.0, 1.0))
    return np.clip(confidence, 0.0, 1.0)


def _profile_uncertainty(results: list[InversionResult]) -> InversionUncertainty | None:
    available = [r.uncertainty for r in results if r.uncertainty is not None]
    if not available:
        return None

    def stack_attr(name: str) -> np.ndarray | None:
        cols = []
        for unc in available:
            value = getattr(unc, name, None)
            if value is None:
                return None
            arr = np.asarray(value, dtype=float)
            cols.append(arr.reshape(arr.shape[0], -1)[:, 0])
        return np.stack(cols, axis=1)

    confidence = stack_attr("confidence")
    return InversionUncertainty(
        model_std=stack_attr("model_std"),
        covariance_diag=stack_attr("covariance_diag"),
        sensitivity=stack_attr("sensitivity"),
        confidence=confidence,
        station_confidence=(
            None if confidence is None else np.nanmean(confidence, axis=0)
        ),
        depth_confidence=(
            None if confidence is None else np.nanmean(confidence, axis=1)
        ),
        metadata={"type": "stitched_station_uncertainty"},
    )


def _profile_history(results: list[InversionResult]) -> InversionHistory | None:
    histories = [r.history for r in results if r.history is not None]
    if not histories:
        return None
    records: list[dict[str, Any]] = []
    for station_index, history in enumerate(histories):
        for record in history.records:
            merged = dict(record)
            merged["station_index"] = float(station_index)
            records.append(merged)
    return InversionHistory(
        records=records,
        metadata={
            "type": "stitched_station_history",
            "n_station_histories": len(histories),
        },
    )


def _bounds(bounds: dict[str, Any], n_layers: int) -> tuple[np.ndarray, np.ndarray]:
    rho_bounds = bounds.get("log10_rho", (-1.0, 6.0))
    thick_bounds = bounds.get("log10_thickness", (0.0, 5.0))
    lower = np.r_[
        np.full(n_layers, float(rho_bounds[0])),
        np.full(n_layers - 1, float(thick_bounds[0])),
    ]
    upper = np.r_[
        np.full(n_layers, float(rho_bounds[1])),
        np.full(n_layers - 1, float(thick_bounds[1])),
    ]
    return lower, upper


class _HistoryRecorder:
    def __init__(self, *, reg_weight: float = 0.0) -> None:
        self.reg_weight = float(reg_weight)
        self.records: list[dict[str, float]] = []

    def add(
        self,
        data_residual: np.ndarray,
        reg_residual: np.ndarray,
        params: np.ndarray,
    ) -> None:
        data_residual = np.asarray(data_residual, dtype=float).reshape(-1)
        reg_residual = np.asarray(reg_residual, dtype=float).reshape(-1)
        params = np.asarray(params, dtype=float).reshape(-1)
        phi_d = float(np.sum(data_residual ** 2))
        phi_m = float(np.sum(reg_residual ** 2))
        self.records.append({
            "iteration": float(len(self.records)),
            "phi_d": phi_d,
            "phi_m": phi_m,
            "objective": phi_d + self.reg_weight * phi_m,
            "rms": float(np.sqrt(np.mean(data_residual ** 2)))
            if data_residual.size else float("nan"),
            "lambda": self.reg_weight,
            "model_norm": float(np.linalg.norm(params)),
        })

    def to_history(self, *, metadata: dict[str, Any] | None = None) -> InversionHistory:
        return InversionHistory(
            records=self.records,
            metadata={
                "type": "residual_evaluation",
                "regularization_weight": self.reg_weight,
                **dict(metadata or {}),
            },
        )


def _rho_grid_bounds(bounds: dict[str, Any], n_params: int) -> tuple[np.ndarray, np.ndarray]:
    rho_bounds = bounds.get("log10_rho", (-1.0, 6.0))
    lower = np.full(n_params, float(rho_bounds[0]), dtype=float)
    upper = np.full(n_params, float(rho_bounds[1]), dtype=float)
    return lower, upper


def _grid_with_core_model(base_grid: Any, log10_core_rho: np.ndarray) -> Any:
    rho_core = 10.0 ** np.asarray(log10_core_rho, dtype=float)
    pad = int(getattr(base_grid, "n_pad", 0))
    rho = np.asarray(base_grid.resistivity, dtype=float).copy()
    nz_core, nx_core = rho_core.shape
    rho[:nz_core, pad: pad + nx_core] = rho_core
    if pad > 0:
        rho[:nz_core, :pad] = rho_core[:, :1]
        rho[:nz_core, pad + nx_core:] = rho_core[:, -1:]
        rho[nz_core:, :] = rho[nz_core - 1:nz_core, :]
    grid = base_grid.__class__(
        dx=base_grid.dx.copy(),
        dz=base_grid.dz.copy(),
        resistivity=rho,
        x_stations=base_grid.x_stations.copy(),
        n_pad=pad,
        name="builtin_fd2d",
    )
    grid._pycsamt_x_offset = float(getattr(base_grid, "_pycsamt_x_offset", 0.0))
    return grid


def _pack_2d_observed(
    em_data,
    cfg: Any,
    opts: dict[str, Any],
) -> tuple[np.ndarray, np.ndarray, tuple[str, ...]]:
    raw_components = opts.get("components", opts.get("component", ("te",)))
    if isinstance(raw_components, str):
        raw_components = (raw_components,)
    modes = tuple(
        str(comp).lower()
        for comp in raw_components
        if str(comp).lower() in {"te", "tm"}
    ) or ("te",)
    obs_parts: list[np.ndarray] = []
    err_parts: list[np.ndarray] = []
    fields: list[str] = []
    rho = None if em_data.rho_a is None else _station_freq_matrix(em_data.rho_a)
    phase = None if em_data.phase is None else _station_freq_matrix(em_data.phase)
    explicit_err = None if em_data.errors is None else _station_freq_matrix(em_data.errors)
    for component in modes:
        if rho is not None:
            mask = component_mask(rho, cfg, component="rho")
            obs_parts.append(np.log10(np.maximum(rho, 1e-12)).reshape(-1))
            if explicit_err is not None:
                rel = component_errors(rho, cfg, component="rho", explicit=explicit_err, relative=True)
                err_parts.append(rel.reshape(-1))
            else:
                err_parts.append(component_errors(rho, cfg, component="rho", relative=True).reshape(-1))
            fields.append(f"{component}_rho")
            if not np.all(mask):
                masked = err_parts[-1].reshape(rho.shape)
                masked[~mask] = 1e30
                err_parts[-1] = masked.reshape(-1)
        if cfg.include_phase and phase is not None:
            mask = component_mask(phase, cfg, component="phase")
            obs_parts.append(phase.reshape(-1))
            err_parts.append(component_errors(phase, cfg, component="phase").reshape(-1))
            fields.append(f"{component}_phase")
            if not np.all(mask):
                masked = err_parts[-1].reshape(phase.shape)
                masked[~mask] = 1e30
                err_parts[-1] = masked.reshape(-1)
    if not obs_parts:
        raise ValueError("finite-difference 2-D inversion needs rho_a and/or phase data.")
    return (
        np.concatenate(obs_parts),
        np.concatenate(err_parts),
        tuple(fields),
    )


def _pack_2d_predicted(response: Any, components: tuple[str, ...]) -> np.ndarray:
    parts: list[np.ndarray] = []
    for component in components:
        mode, field = component.split("_", 1)
        if field == "rho":
            rho = getattr(response, f"rho_a_{mode}")
            parts.append(
                np.log10(np.maximum(np.asarray(rho, dtype=float).T, 1e-12)).reshape(-1)
            )
        elif field == "phase":
            phase = getattr(response, f"phase_{mode}")
            parts.append(np.asarray(phase, dtype=float).T.reshape(-1))
    return np.concatenate(parts)


def _station_freq_matrix(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.ndim == 1:
        return values.reshape(1, -1)
    return values


def _layer_centers(thicknesses: np.ndarray) -> np.ndarray:
    tops = np.r_[0.0, np.cumsum(thicknesses)]
    last = thicknesses[-1] if thicknesses.size else 1.0
    bottoms = np.r_[tops[1:], tops[-1] + last]
    return 0.5 * (tops + bottoms)


def _log_abs(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    return np.log10(np.maximum(np.abs(values), 1e-30))


def _station_data(em_data, idx: int):
    from ..data import EMData

    rho = None if em_data.rho_a is None else _row(em_data.rho_a, idx)
    phase = None if em_data.phase is None else _row(em_data.phase, idx)
    errors = None if em_data.errors is None else _row(em_data.errors, idx)
    names = _station_names(em_data, em_data.n_stations)
    xs = _station_x(em_data, em_data.n_stations)
    return EMData(
        method=em_data.method,
        frequencies=em_data.frequencies,
        times=em_data.times,
        rho_a=rho,
        phase=phase,
        values=None if em_data.values is None else _row(em_data.values, idx),
        errors=errors,
        station_names=[names[idx]],
        station_x=np.array([xs[idx]], dtype=float),
        source=em_data.source,
        metadata=em_data.metadata_dict(),
    )


def _row(arr: np.ndarray, idx: int) -> np.ndarray:
    arr = np.asarray(arr, dtype=float)
    if arr.ndim == 1:
        return arr.copy()
    return arr[idx, :].copy()


def _station_names(em_data, n_st: int) -> list[str]:
    if em_data.station_names:
        return list(em_data.station_names)
    return [f"S{i:03d}" for i in range(n_st)]


def _station_x(em_data, n_st: int) -> np.ndarray:
    if em_data.station_x is not None:
        return np.asarray(em_data.station_x, dtype=float)
    return np.arange(n_st, dtype=float)


Builtin1DBackend.__doc__ = r"""
Run PyCSAMT's built-in EM inversion engine.

``Builtin1DBackend`` is the dependency-light backend used when an inversion
should run without SimPEG, pyGIMLi, Occam2D, or ModEM. Despite the historical
class name, it supports 1-D soundings and selected 2-D profile workflows for
MT, AMT, CSAMT, EMAP, and TDEM data.

The backend minimizes a weighted residual of the form

.. math::

    \Phi(m) =
        \| W_d (d_\mathrm{pred}(m) - d_\mathrm{obs}) \|_2^2
        + \lambda \Phi_m(m),

where :math:`W_d` is built from component errors and
:math:`\Phi_m` is supplied by
:mod:`pycsamt.inversion.regularization`. Natural-source 1-D runs invert
log10 apparent resistivity and optional phase. TDEM runs invert log10 absolute
transient response. Stitched 2-D mode repeats the 1-D inversion station by
station, while finite-difference 2-D mode uses the local MT profile forward
solver.

Parameters
----------
config : pycsamt.inversion.config.InversionConfig
    Normalized inversion configuration. The backend reads ``method``,
    ``dimension``, ``data``, ``starting_model``, ``reference_model``,
    ``bounds``, ``regularization``, ``max_iter``, ``tol``, ``workdir``, and
    ``backend_options`` from this object.

Attributes
----------
name : str
    Backend registry name, always ``"builtin"``.
supports : tuple of tuple
    Supported ``(method, dimension)`` combinations.

Notes
-----
Important ``backend_options`` entries include:

``profile_mode``
    For natural-source 2-D runs, use ``"fd2d"`` or
    ``"finite_difference"`` to select true finite-difference profile
    inversion. Any other value uses stitched station-wise 1-D inversion.
``components``
    Natural-source 2-D finite-difference components. Supported values are
    ``"te"``, ``"tm"``, or a sequence containing both.
``blocky_eps``
    Small stabilizer used by blocky regularization.
``loop_radius``, ``moment``, ``n_freqs``, ``n_lam``
    TDEM forward-model controls passed to
    :class:`pycsamt.forward.TEM1DForward`.

Examples
--------
Run a small MT sounding through the backend directly::

    >>> from pycsamt.inversion.backends.builtin import Builtin1DBackend
    >>> from pycsamt.inversion.config import InversionConfig
    >>> cfg = InversionConfig(
    ...     method="mt",
    ...     dimension="1d",
    ...     backend="builtin",
    ...     data={"freqs": [1.0, 10.0],
    ...           "rho_a": [100.0, 120.0],
    ...           "phase": [45.0, 47.0]},
    ...     max_iter=4,
    ... )
    >>> result = Builtin1DBackend(cfg).run()  # doctest: +SKIP

Run the same engine through the high-level workflow::

    >>> from pycsamt.inversion.workflow import run_inversion
    >>> result = run_inversion(
    ...     method="tdem",
    ...     dimension="1d",
    ...     backend="builtin",
    ...     data={"times": [1e-5, 1e-4],
    ...           "values": [1e-8, 2e-9]},
    ...     max_iter=4,
    ... )  # doctest: +SKIP

Select finite-difference 2-D MT profile mode::

    >>> from pycsamt.inversion.workflow import run_inversion
    >>> result = run_inversion(
    ...     method="mt",
    ...     dimension="2d",
    ...     backend="builtin",
    ...     data={"freqs": [1.0, 10.0],
    ...           "rho_a": [[100.0, 120.0], [90.0, 110.0]],
    ...           "phase": [[45.0, 47.0], [44.0, 46.0]],
    ...           "station_x": [0.0, 500.0]},
    ...     backend_options={"profile_mode": "fd2d",
    ...                      "components": ("te",)},
    ...     max_iter=4,
    ... )  # doctest: +SKIP

See Also
--------
pycsamt.inversion.workflow.InversionWorkflow
    High-level orchestration API that usually instantiates this backend.
pycsamt.inversion.config.InversionConfig
    Configuration object consumed by the backend.
pycsamt.inversion.regularization.Regularization
    Shared regularization description wired into the residual.
pycsamt.inversion.results.InversionResult
    Backend-neutral result returned by :meth:`run`.

References
----------
.. [1] Aster, R. C., Borchers, B. and Thurber, C. H. (2018). *Parameter
       Estimation and Inverse Problems*, 3rd edition. Elsevier.
.. [2] Constable, S. C., Parker, R. L. and Constable, C. G. (1987). Occam's
       inversion: A practical algorithm for generating smooth models from
       electromagnetic sounding data. *Geophysics*, 52(3), 289-300.
.. [3] Chave, A. D. and Jones, A. G. (2012). *The Magnetotelluric Method:
       Theory and Practice*. Cambridge University Press.
"""
