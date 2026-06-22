# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""SimPEG backend for physics-based EM inversion.

SimPEG is optional and imported only when this backend is selected. The
backend wraps SimPEG's natural-source electromagnetic API for MT, AMT, and
CSAMT inversion. It supports 1-D soundings, stitched station-by-station 2-D
profiles, and a 3-D primary-secondary natural-source path when the installed
SimPEG version exposes the required classes.

The model parameter is log conductivity on SimPEG mesh cells. PyCSAMT converts
user-facing layered resistivity starting models into log-conductivity cells,
runs the inversion, and converts recovered models back to resistivity-oriented
``InversionResult`` payloads.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from ..base import BaseInversionBackend
from ..data import EMData
from ..mesh import InversionMesh, build_1d_tensor_mesh, build_3d_tensor_mesh
from ..model import StartingModel
from ..objective import component_errors, component_mask, weighted_rms
from ..regularization import regularization_from_config
from ..results import InversionResult

__all__ = ["SimPEGBackend"]


@dataclass
class _SimPEGModules:
    discretize: Any
    maps: Any
    data: Any
    data_misfit: Any
    directives: Any
    inverse_problem: Any
    inversion: Any
    optimization: Any
    regularization: Any
    nsem: Any


class SimPEGBackend(BaseInversionBackend):
    name = "simpeg"
    supports = (
        ("mt", "1d"),
        ("mt", "2d"),
        ("mt", "3d"),
        ("amt", "1d"),
        ("amt", "2d"),
        ("amt", "3d"),
        ("csamt", "1d"),
        ("csamt", "2d"),
        ("csamt", "3d"),
    )

    def run(self, data: Any | None = None) -> InversionResult:
        """Run a SimPEG natural-source EM inversion.

        Parameters
        ----------
        data : mapping, object, sequence, or path-like, optional
            Optional data override for this call. When omitted, the backend
            uses ``self.config.data``. Values are coerced through
            :class:`pycsamt.inversion.data.EMData`.

        Returns
        -------
        InversionResult
            Backend-neutral result. For 1-D runs, ``model`` is a recovered
            :class:`pycsamt.inversion.model.StartingModel`. For stitched 2-D
            runs, ``model`` is a dictionary containing ``rho_2d`` and profile
            coordinates. For 3-D runs, ``model`` contains ``rho_3d`` and
            ``x_centers``, ``y_centers``, and ``z_centers`` arrays.

        Raises
        ------
        ImportError
            If SimPEG or discretize is not installed.
        ValueError
            If data do not contain frequencies plus apparent resistivity
            and/or phase.
        NotImplementedError
            If the configured method/dimension pair is unsupported by this
            backend or by the installed SimPEG natural-source API.

        Examples
        --------
        >>> from pycsamt.inversion.backends.simpeg import SimPEGBackend
        >>> from pycsamt.inversion.config import InversionConfig
        >>> cfg = InversionConfig(
        ...     method="mt",
        ...     dimension="1d",
        ...     backend="simpeg",
        ...     data={"freqs": [1.0, 10.0],
        ...           "rho_a": [100.0, 120.0],
        ...           "phase": [45.0, 47.0]},
        ...     max_iter=8,
        ... )
        >>> result = SimPEGBackend(cfg).run()  # doctest: +SKIP
        """
        self.check_supported()
        modules = _load_simpeg()
        em_data = self.prepare_data(data)
        if not em_data.has_mt_response:
            raise ValueError(
                "SimPEG natural-source backend requires frequencies plus "
                "rho_a and/or phase."
            )
        if self.config.dimension == "3d":
            return self._run_3d(em_data, modules)
        if self.config.dimension == "2d":
            return self._run_profile(em_data, modules)
        return self._run_sounding(em_data, modules, station_index=None)

    def _run_3d(self, em_data: EMData, modules: _SimPEGModules) -> InversionResult:
        cfg = self.config
        mesh, centers = _build_3d_mesh(em_data, cfg.backend_options, modules)
        survey = _build_nsem_survey(em_data, modules, dimension="3d")
        observed, errors = _pack_nsem_observations(em_data, cfg)

        active_map = modules.maps.IdentityMap(nP=mesh.nC)
        sigma_map = modules.maps.ExpMap(mesh) * active_map
        sigma_primary = float(cfg.backend_options.get("sigma_primary", 1.0 / 100.0))
        simulation = modules.nsem.Simulation3DPrimarySecondary(
            mesh,
            survey=survey,
            sigmaMap=sigma_map,
            sigmaPrimary=np.full(mesh.nC, sigma_primary, dtype=float),
        )
        simpeg_data = _build_simpeg_data(survey, observed, errors, modules)
        data_misfit = modules.data_misfit.L2DataMisfit(
            data=simpeg_data,
            simulation=simulation,
        )
        reg = _build_regularization(mesh, active_map, modules, cfg)
        opt = modules.optimization.InexactGaussNewton(
            maxIter=cfg.max_iter,
            maxIterLS=int(cfg.backend_options.get("max_iter_ls", 20)),
            tolX=cfg.tol,
            tolF=cfg.tol,
        )
        inv_problem = modules.inverse_problem.BaseInvProblem(data_misfit, reg, opt)
        beta0 = float(cfg.backend_options.get("beta0", 1.0))
        inv_problem.beta = beta0
        inv = modules.inversion.BaseInversion(
            inv_problem,
            directiveList=_build_directives(modules, cfg),
        )

        start = StartingModel.coerce(cfg.starting_model, n_layers=cfg.n_layers)
        m0 = _starting_3d_log_sigma(start, centers)
        recovered_model = inv.run(m0)
        predicted = np.asarray(simulation.dpred(recovered_model), dtype=float)
        rms = weighted_rms(observed, predicted, errors)
        rho_3d = _rho_3d_from_log_sigma(recovered_model, mesh)
        mesh_out = InversionMesh(
            dimension="3d",
            x_centers=centers["x"],
            z_centers=centers["z_depth"],
            native=mesh,
            metadata={
                "engine": "simpeg",
                "y_centers": centers["y"].tolist(),
                "mesh_shape": tuple(int(v) for v in _mesh_shape(mesh)),
            },
        )
        return InversionResult(
            method=cfg.method,
            dimension="3d",
            backend=self.name,
            status="success",
            model={
                "rho_3d": rho_3d,
                "x_centers": centers["x"],
                "y_centers": centers["y"],
                "z_centers": centers["z_depth"],
                "station_x": _station_x(em_data, em_data.n_stations),
                "station_y": _station_y(em_data, em_data.n_stations),
                "station_names": _station_names(em_data, em_data.n_stations),
            },
            mesh=mesh_out,
            data=em_data,
            predicted=predicted,
            rms=rms,
            objective=float(data_misfit(recovered_model)),
            n_iter=int(getattr(opt, "iter", 0)),
            workdir=cfg.workdir,
            native={
                "mesh": mesh,
                "survey": survey,
                "data": simpeg_data,
                "simulation": simulation,
                "data_misfit": data_misfit,
                "regularization": reg,
                "optimization": opt,
                "inverse_problem": inv_problem,
                "inversion": inv,
                "recovered_model": recovered_model,
            },
            metadata={
                **cfg.metadata,
                "engine": "simpeg",
                "simulation": "Simulation3DPrimarySecondary",
                "source": "PlanewaveXYPrimary",
                "model_parameter": "log_sigma_cells",
                "beta0": beta0,
                "sigma_primary": sigma_primary,
            },
        )

    def _run_profile(self, em_data: EMData, modules: _SimPEGModules) -> InversionResult:
        cfg = self.config
        n_st = em_data.n_stations
        names = _station_names(em_data, n_st)
        xs = _station_x(em_data, n_st)
        columns: list[np.ndarray] = []
        station_results: list[InversionResult] = []
        used: list[int] = []
        warnings: list[str] = []
        z_centers = None

        for idx in range(n_st):
            sounding = _station_data(em_data, idx)
            try:
                result = self._run_sounding(sounding, modules, station_index=idx)
            except Exception as exc:
                warnings.append(f"{names[idx]}: SimPEG inversion failed: {exc}")
                continue
            station_results.append(result)
            used.append(idx)
            columns.append(np.log10(result.model.resistivities))
            if z_centers is None and result.mesh is not None:
                z_centers = result.mesh.z_centers

        if not columns:
            raise RuntimeError("all SimPEG station inversions failed.")

        rho_2d = np.stack(columns, axis=1)
        used_x = xs[used]
        used_names = [names[i] for i in used]
        if z_centers is None:
            z_centers = np.arange(rho_2d.shape[0], dtype=float)
        mesh = InversionMesh(
            dimension="2d",
            x_centers=used_x,
            z_centers=z_centers,
            metadata={"engine": "simpeg", "profile_mode": "stitched_1d"},
        )
        rms_values = np.asarray([r.rms for r in station_results], dtype=float)
        return InversionResult(
            method=cfg.method,
            dimension="2d",
            backend=self.name,
            status="success" if not warnings else "needs_review",
            model={
                "rho_2d": rho_2d,
                "x_centers": used_x,
                "z_centers": z_centers,
                "station_x": used_x,
                "station_names": used_names,
            },
            mesh=mesh,
            data=em_data,
            predicted=[r.predicted for r in station_results],
            rms=float(np.nanmean(rms_values)),
            objective=float(np.nansum([r.objective for r in station_results])),
            n_iter=int(np.sum([r.n_iter for r in station_results])),
            workdir=cfg.workdir,
            native=station_results,
            warnings=warnings,
            metadata={
                **cfg.metadata,
                "engine": "simpeg",
                "profile_mode": "stitched_station_1d",
                "station_rms": rms_values.tolist(),
            },
        )

    def _run_sounding(
        self,
        em_data: EMData,
        modules: _SimPEGModules,
        *,
        station_index: int | None,
    ) -> InversionResult:
        cfg = self.config
        start = StartingModel.coerce(cfg.starting_model, n_layers=cfg.n_layers)
        mesh, z_centers = _build_1d_mesh(start, cfg.backend_options, modules)
        survey = _build_nsem_survey(em_data, modules)
        observed, errors = _pack_nsem_observations(em_data, cfg)

        active_map = modules.maps.IdentityMap(nP=mesh.nC)
        sigma_map = modules.maps.ExpMap(mesh) * active_map
        simulation = modules.nsem.Simulation1DElectricField(
            mesh,
            survey=survey,
            sigmaMap=sigma_map,
        )

        simpeg_data = _build_simpeg_data(survey, observed, errors, modules)
        data_misfit = modules.data_misfit.L2DataMisfit(
            data=simpeg_data,
            simulation=simulation,
        )
        reg = _build_regularization(mesh, active_map, modules, cfg)
        opt = modules.optimization.InexactGaussNewton(
            maxIter=cfg.max_iter,
            maxIterLS=int(cfg.backend_options.get("max_iter_ls", 20)),
            tolX=cfg.tol,
            tolF=cfg.tol,
        )
        inv_problem = modules.inverse_problem.BaseInvProblem(data_misfit, reg, opt)
        beta0 = float(cfg.backend_options.get("beta0", 1.0))
        inv_problem.beta = beta0
        directive_list = _build_directives(modules, cfg)
        inv = modules.inversion.BaseInversion(inv_problem, directiveList=directive_list)

        m0 = _starting_sigma_model(start, z_centers)
        recovered_model = inv.run(m0)
        predicted = np.asarray(simulation.dpred(recovered_model), dtype=float)
        rms = weighted_rms(observed, predicted, errors)
        recovered = _model_from_sigma_cells(recovered_model, z_centers, start.n_layers)
        mesh_out = InversionMesh(
            dimension="1d",
            x_centers=np.array([0.0]),
            z_centers=_layer_centers(recovered.thicknesses),
            native=mesh,
            metadata={
                "engine": "simpeg",
                "n_cells": int(mesh.nC),
                "cell_z_centers": z_centers.tolist(),
            },
        )
        return InversionResult(
            method=cfg.method,
            dimension="1d",
            backend=self.name,
            status="success",
            model=recovered,
            mesh=mesh_out,
            data=em_data,
            predicted=predicted,
            rms=rms,
            objective=float(data_misfit(recovered_model)),
            n_iter=int(getattr(opt, "iter", 0)),
            workdir=cfg.workdir,
            native={
                "mesh": mesh,
                "survey": survey,
                "data": simpeg_data,
                "simulation": simulation,
                "data_misfit": data_misfit,
                "regularization": reg,
                "optimization": opt,
                "inverse_problem": inv_problem,
                "inversion": inv,
                "recovered_model": recovered_model,
            },
            metadata={
                **cfg.metadata,
                "engine": "simpeg",
                "model_parameter": "log_sigma_cells",
                "station_index": station_index,
                "beta0": beta0,
            },
        )


def _load_simpeg() -> _SimPEGModules:
    try:
        import discretize
        from simpeg import (
            data,
            data_misfit,
            directives,
            inverse_problem,
            inversion,
            maps,
            optimization,
            regularization,
        )
        from simpeg.electromagnetics import natural_source as nsem
    except ImportError as exc:
        raise ImportError(
            "SimPEG backend selected, but SimPEG/discretize is not installed. "
            "Install SimPEG, or choose backend='builtin'/'occam2d'."
        ) from exc
    return _SimPEGModules(
        discretize=discretize,
        maps=maps,
        data=data,
        data_misfit=data_misfit,
        directives=directives,
        inverse_problem=inverse_problem,
        inversion=inversion,
        optimization=optimization,
        regularization=regularization,
        nsem=nsem,
    )


def _build_1d_mesh(
    start: StartingModel,
    options: dict[str, Any],
    modules: _SimPEGModules,
):
    return build_1d_tensor_mesh(start, options, modules.discretize.TensorMesh)


def _build_3d_mesh(
    em_data: EMData,
    options: dict[str, Any],
    modules: _SimPEGModules,
):
    station_x = _station_x(em_data, em_data.n_stations)
    station_y = _station_y(em_data, em_data.n_stations)
    return build_3d_tensor_mesh(
        station_x,
        station_y,
        options,
        modules.discretize.TensorMesh,
    )


def _build_nsem_survey(
    em_data: EMData,
    modules: _SimPEGModules,
    *,
    dimension: str = "1d",
):
    receivers = modules.nsem.receivers
    sources = modules.nsem.sources
    survey_cls = getattr(modules.nsem, "Survey", None)
    if survey_cls is None:
        from simpeg import survey as survey_mod
        survey_cls = survey_mod.Survey
    freqs = np.asarray(em_data.frequencies, dtype=float)
    if dimension == "3d":
        location = _station_locations(em_data)
    else:
        location = np.array([[0.0, 0.0, 0.0]], dtype=float)
    source_list = []
    for freq in freqs:
        rx_list = []
        if em_data.rho_a is not None:
            rx_list.append(receivers.PointNaturalSource(
                location,
                orientation="xy",
                component="apparent_resistivity",
            ))
        if em_data.phase is not None:
            rx_list.append(receivers.PointNaturalSource(
                location,
                orientation="xy",
                component="phase",
            ))
        if dimension == "3d" and hasattr(sources, "PlanewaveXYPrimary"):
            source_list.append(sources.PlanewaveXYPrimary(rx_list, frequency=float(freq)))
        else:
            source_list.append(sources.Planewave(rx_list, frequency=float(freq)))
    return survey_cls(source_list)


def _build_simpeg_data(
    survey: Any,
    observed: np.ndarray,
    errors: np.ndarray,
    modules: _SimPEGModules,
):
    try:
        return modules.data.Data(
            survey,
            dobs=observed,
            standard_deviation=errors,
        )
    except TypeError:
        return modules.nsem.survey.Data(
            survey,
            dobs=observed,
            standard_deviation=errors,
        )


def _pack_nsem_observations(
    em_data: EMData,
    cfg: Any,
) -> tuple[np.ndarray, np.ndarray]:
    values: list[float] = []
    errors: list[float] = []
    rho = None if em_data.rho_a is None else np.asarray(em_data.rho_a, dtype=float)
    phase = None if em_data.phase is None else np.asarray(em_data.phase, dtype=float)
    raw_err = None if em_data.errors is None else np.asarray(em_data.errors, dtype=float)
    n = em_data.n_samples
    for i in range(n):
        if rho is not None:
            rho_i = rho[i] if rho.ndim == 1 else rho[:, i]
            mask_i = component_mask(rho_i, cfg, component="rho").reshape(-1)
            values.extend(np.asarray(rho_i, dtype=float).reshape(-1).tolist())
            err_i = raw_err[i] if raw_err is not None and raw_err.ndim == 1 else (
                raw_err[:, i] if raw_err is not None else None
            )
            err = component_errors(rho_i, cfg, component="rho", explicit=err_i).reshape(-1)
            err[~mask_i] = 1e30
            errors.extend(err.tolist())
        if phase is not None:
            phase_i = phase[i] if phase.ndim == 1 else phase[:, i]
            phase_i = np.asarray(phase_i, dtype=float).reshape(-1)
            mask_i = component_mask(phase_i, cfg, component="phase").reshape(-1)
            values.extend(phase_i.tolist())
            err = component_errors(phase_i, cfg, component="phase").reshape(-1)
            err[~mask_i] = 1e30
            errors.extend(err.tolist())
    return np.asarray(values, dtype=float), np.asarray(errors, dtype=float)


def _build_regularization(
    mesh: Any,
    mapping: Any,
    modules: _SimPEGModules,
    cfg: Any,
):
    reg_cfg = regularization_from_config(cfg)
    kind = reg_cfg.kind
    if hasattr(modules.regularization, "WeightedLeastSquares"):
        reg = modules.regularization.WeightedLeastSquares(mesh, mapping=mapping)
    else:
        reg = modules.regularization.Simple(mesh, mapping=mapping)
    if kind == "none":
        _set_if_present(reg, "alpha_s", 0.0)
        _set_if_present(reg, "alpha_x", 0.0)
        _set_if_present(reg, "alpha_y", 0.0)
        _set_if_present(reg, "alpha_z", 0.0)
    elif kind == "damped":
        _set_if_present(reg, "alpha_s", reg_cfg.alpha_s)
        _set_if_present(reg, "alpha_x", 0.0)
        _set_if_present(reg, "alpha_y", 0.0)
        _set_if_present(reg, "alpha_z", 0.0)
    else:
        _set_if_present(reg, "alpha_s", reg_cfg.alpha_s)
        _set_if_present(reg, "alpha_x", reg_cfg.alpha_x)
        _set_if_present(reg, "alpha_y", float(cfg.backend_options.get("alpha_y", reg_cfg.alpha_x)))
        _set_if_present(reg, "alpha_z", reg_cfg.alpha_z)
    reference = _simpeg_reference_model(cfg, mesh)
    if reference is not None:
        _set_if_present(reg, "reference_model", reference)
        _set_if_present(reg, "mref", reference)
    return reg


def _set_if_present(obj: Any, name: str, value: Any) -> None:
    if hasattr(obj, name):
        setattr(obj, name, value)


def _simpeg_reference_model(cfg: Any, mesh: Any) -> np.ndarray | None:
    value = cfg.backend_options.get("reference_model", cfg.reference_model)
    if value is None:
        return None
    raw = getattr(value, "resistivities", value)
    try:
        arr = np.asarray(raw, dtype=float).reshape(-1)
    except Exception:
        return None
    n_cells = int(getattr(mesh, "nC", arr.size))
    if arr.size != n_cells:
        return None
    if np.all(arr > 0.0):
        return np.log(1.0 / arr)
    return arr


def _build_directives(modules: _SimPEGModules, cfg: Any) -> list[Any]:
    directives = []
    if bool(cfg.backend_options.get("estimate_beta", True)):
        beta_cls = getattr(modules.directives, "BetaEstimate_ByEig", None)
        if beta_cls is not None:
            directives.append(beta_cls(beta0_ratio=float(cfg.backend_options.get("beta0_ratio", 1.0))))
    target_cls = getattr(modules.directives, "TargetMisfit", None)
    if target_cls is not None:
        directives.append(target_cls(chifact=float(cfg.backend_options.get("target_chifact", 1.0))))
    beta_schedule = getattr(modules.directives, "BetaSchedule", None)
    if beta_schedule is not None:
        directives.append(beta_schedule(
            coolingFactor=float(cfg.backend_options.get("cooling_factor", 2.0)),
            coolingRate=int(cfg.backend_options.get("cooling_rate", 1)),
        ))
    return directives


def _starting_sigma_model(start: StartingModel, z_centers: np.ndarray) -> np.ndarray:
    rho = np.asarray(start.resistivities, dtype=float)
    depths = np.r_[0.0, np.cumsum(start.thicknesses)]
    out = np.empty_like(z_centers, dtype=float)
    for i, z in enumerate(z_centers):
        layer = int(np.searchsorted(depths[1:], z, side="right"))
        out[i] = np.log(1.0 / rho[min(layer, rho.size - 1)])
    return out


def _starting_3d_log_sigma(
    start: StartingModel,
    centers: dict[str, np.ndarray],
) -> np.ndarray:
    z_depth = np.asarray(centers["z_depth"], dtype=float)
    one_d = _starting_sigma_model(start, z_depth)
    nx = len(centers["x"])
    ny = len(centers["y"])
    nz = len(centers["z"])
    model = np.empty((nx, ny, nz), dtype=float)
    for iz in range(nz):
        model[:, :, iz] = one_d[iz]
    return model.reshape(-1, order="F")


def _rho_3d_from_log_sigma(log_sigma: np.ndarray, mesh: Any) -> np.ndarray:
    shape = _mesh_shape(mesh)
    sigma = np.exp(np.asarray(log_sigma, dtype=float)).reshape(shape, order="F")
    return 1.0 / np.maximum(sigma, 1e-12)


def _mesh_shape(mesh: Any) -> tuple[int, int, int]:
    if hasattr(mesh, "shape_cells"):
        shape = tuple(int(v) for v in mesh.shape_cells)
    elif hasattr(mesh, "vnC"):
        shape = tuple(int(v) for v in mesh.vnC)
    else:
        raise AttributeError("Cannot determine SimPEG mesh cell shape.")
    if len(shape) != 3:
        raise ValueError("Expected a 3-D mesh.")
    return shape


def _model_from_sigma_cells(
    log_sigma: np.ndarray,
    z_centers: np.ndarray,
    n_layers: int,
) -> StartingModel:
    log_sigma = np.asarray(log_sigma, dtype=float)
    z_centers = np.asarray(z_centers, dtype=float)
    if n_layers < 2:
        n_layers = 2
    edges = np.linspace(float(z_centers[0]), float(z_centers[-1]), n_layers + 1)
    resistivities = []
    thicknesses = []
    for layer in range(n_layers):
        lo = edges[layer]
        hi = edges[layer + 1]
        mask = (z_centers >= lo) & (z_centers <= hi if layer == n_layers - 1 else z_centers < hi)
        if not np.any(mask):
            idx = int(np.argmin(np.abs(z_centers - 0.5 * (lo + hi))))
            mask[idx] = True
        sigma = np.exp(float(np.nanmean(log_sigma[mask])))
        resistivities.append(1.0 / max(sigma, 1e-12))
        if layer < n_layers - 1:
            thicknesses.append(max(float(hi - lo), 1.0))
    return StartingModel(resistivities, thicknesses, name="simpeg_1d")


def _layer_centers(thicknesses: np.ndarray) -> np.ndarray:
    tops = np.r_[0.0, np.cumsum(thicknesses)]
    last = thicknesses[-1] if thicknesses.size else 1.0
    bottoms = np.r_[tops[1:], tops[-1] + last]
    return 0.5 * (tops + bottoms)


def _station_data(em_data: EMData, idx: int) -> EMData:
    return EMData(
        method=em_data.method,
        frequencies=em_data.frequencies,
        rho_a=None if em_data.rho_a is None else _row(em_data.rho_a, idx),
        phase=None if em_data.phase is None else _row(em_data.phase, idx),
        errors=None if em_data.errors is None else _row(em_data.errors, idx),
        station_names=[_station_names(em_data, em_data.n_stations)[idx]],
        station_x=np.array([_station_x(em_data, em_data.n_stations)[idx]], dtype=float),
        source=em_data.source,
        metadata=em_data.metadata_dict(),
    )


def _row(arr: np.ndarray, idx: int) -> np.ndarray:
    arr = np.asarray(arr, dtype=float)
    if arr.ndim == 1:
        return arr.copy()
    return arr[idx, :].copy()


def _station_names(em_data: EMData, n_st: int) -> list[str]:
    if em_data.station_names:
        return list(em_data.station_names)
    return [f"S{i:03d}" for i in range(n_st)]


def _station_x(em_data: EMData, n_st: int) -> np.ndarray:
    if em_data.station_x is not None:
        return np.asarray(em_data.station_x, dtype=float)
    return np.arange(n_st, dtype=float)


def _station_y(em_data: EMData, n_st: int) -> np.ndarray:
    meta = em_data.metadata_dict()
    for key in ("station_y", "y", "northing"):
        if key in meta:
            arr = np.asarray(meta[key], dtype=float)
            if arr.size == n_st:
                return arr
    return np.zeros(n_st, dtype=float)


def _station_locations(em_data: EMData) -> np.ndarray:
    n_st = em_data.n_stations
    return np.c_[
        _station_x(em_data, n_st),
        _station_y(em_data, n_st),
        np.zeros(n_st, dtype=float),
    ]


SimPEGBackend.__doc__ = r"""
Run optional SimPEG natural-source EM inversions.

``SimPEGBackend`` adapts PyCSAMT's inversion configuration objects to SimPEG's
natural-source electromagnetic machinery. It is intended for users who want a
physics-backed inversion engine while keeping the same
:class:`pycsamt.inversion.workflow.InversionWorkflow` and
:class:`pycsamt.inversion.results.InversionResult` API used by the built-in and
external backends.

The backend supports MT, AMT, and CSAMT. One-dimensional runs use
``Simulation1DElectricField``. Two-dimensional runs are stitched
station-by-station 1-D inversions, not full 2-D SimPEG EM inversion. Three-
dimensional runs use ``Simulation3DPrimarySecondary`` with plane-wave natural-
source survey objects when available in the installed SimPEG version.

Parameters
----------
config : pycsamt.inversion.config.InversionConfig
    Inversion configuration. Important fields are ``method``, ``dimension``,
    ``data``, ``starting_model``, ``reference_model``, ``regularization``,
    ``max_iter``, ``tol``, ``error_floor``, ``phase_error``,
    ``backend_options``, ``workdir``, and ``metadata``.

Attributes
----------
name : str
    Backend registry name, always ``"simpeg"``.
supports : tuple of tuple
    Supported ``(method, dimension)`` pairs.

Notes
-----
The inversion model is log conductivity, :math:`m = \log(\sigma)`, on SimPEG
mesh cells. Starting and recovered models exposed through PyCSAMT are
resistivity oriented. For 1-D output, recovered cell conductivities are
compressed back into a layered
:class:`pycsamt.inversion.model.StartingModel`. For 3-D output, the result
stores ``rho_3d`` as linear resistivity.

``backend_options`` can contain:

``mesh`` controls
    Options consumed by :func:`pycsamt.inversion.mesh.build_1d_tensor_mesh`
    and :func:`pycsamt.inversion.mesh.build_3d_tensor_mesh`, such as core
    cell widths, padding, depth extent, and station padding.
``max_iter_ls``
    Maximum line-search iterations passed to
    ``optimization.InexactGaussNewton``.
``beta0``, ``estimate_beta``, ``beta0_ratio``
    Initial trade-off and beta-estimation controls.
``target_chifact``
    Target misfit directive chi factor.
``cooling_factor`` and ``cooling_rate``
    Beta schedule directive controls.
``sigma_primary``
    Background conductivity for 3-D primary-secondary simulation.
``reference_model``
    Cell-sized reference model. Positive values are treated as resistivity and
    converted to log conductivity; otherwise values are assumed already in the
    SimPEG model domain.
``alpha_y``
    Optional y-direction smoothness weight for 3-D regularization.

Examples
--------
Run a 1-D MT sounding directly through the backend::

    >>> from pycsamt.inversion.backends.simpeg import SimPEGBackend
    >>> from pycsamt.inversion.config import InversionConfig
    >>> cfg = InversionConfig(
    ...     method="mt",
    ...     dimension="1d",
    ...     backend="simpeg",
    ...     data={"freqs": [1.0, 10.0],
    ...           "rho_a": [100.0, 120.0],
    ...           "phase": [45.0, 47.0]},
    ...     max_iter=8,
    ... )
    >>> result = SimPEGBackend(cfg).run()  # doctest: +SKIP

Run a stitched 2-D AMT profile through the high-level workflow::

    >>> from pycsamt.inversion.workflow import run_inversion
    >>> result = run_inversion(
    ...     method="amt",
    ...     dimension="2d",
    ...     backend="simpeg",
    ...     data={"freqs": [10.0, 100.0],
    ...           "rho_a": [[80.0, 100.0], [90.0, 110.0]],
    ...           "phase": [[42.0, 45.0], [43.0, 46.0]],
    ...           "station_x": [0.0, 250.0],
    ...           "station_names": ["A01", "A02"]},
    ...     max_iter=6,
    ... )  # doctest: +SKIP
    >>> result.metadata["profile_mode"]  # doctest: +SKIP
    'stitched_station_1d'

Configure a small 3-D natural-source run::

    >>> from pycsamt.inversion.workflow import run_inversion
    >>> result = run_inversion(
    ...     method="mt",
    ...     dimension="3d",
    ...     backend="simpeg",
    ...     data={"freqs": [1.0],
    ...           "rho_a": [[100.0], [120.0]],
    ...           "phase": [[45.0], [47.0]],
    ...           "station_x": [0.0, 500.0],
    ...           "station_names": ["M01", "M02"],
    ...           "metadata": {"station_y": [0.0, 250.0]}},
    ...     backend_options={
    ...         "sigma_primary": 0.01,
    ...         "beta0": 1.0,
    ...         "target_chifact": 1.0,
    ...     },
    ...     max_iter=4,
    ... )  # doctest: +SKIP

Pass SimPEG directive and optimization controls::

    >>> from pycsamt.inversion.workflow import run_inversion
    >>> result = run_inversion(
    ...     method="csamt",
    ...     dimension="1d",
    ...     backend="simpeg",
    ...     data={"freqs": [1.0, 10.0], "rho_a": [100.0, 120.0]},
    ...     backend_options={
    ...         "estimate_beta": True,
    ...         "beta0_ratio": 2.0,
    ...         "cooling_factor": 2.0,
    ...         "cooling_rate": 1,
    ...         "max_iter_ls": 10,
    ...     },
    ... )  # doctest: +SKIP

See Also
--------
pycsamt.inversion.workflow.InversionWorkflow
    High-level entry point that instantiates this backend.
pycsamt.inversion.backends.builtin.Builtin1DBackend
    Dependency-light fallback backend.
pycsamt.inversion.backends.pygimli.PyGIMLiBackend
    Optional pyGIMLi backend for 1-D EM inversions.
pycsamt.inversion.mesh.build_1d_tensor_mesh
    Mesh helper used by the 1-D SimPEG path.
pycsamt.inversion.mesh.build_3d_tensor_mesh
    Mesh helper used by the 3-D SimPEG path.
pycsamt.inversion.results.InversionResult
    Backend-neutral result returned by :meth:`run`.

References
----------
.. [1] Cockett, R., Kang, S., Heagy, L. J., Pidlisecky, A. and Oldenburg,
       D. W. (2015). SimPEG: An open source framework for simulation and
       gradient based parameter estimation in geophysical applications.
       *Computers & Geosciences*, 85, 142-154.
.. [2] Heagy, L. J., Cockett, R., Kang, S., Rosenkjaer, G. K. and Oldenburg,
       D. W. (2017). A framework for simulation and inversion in
       electromagnetics. *Computers & Geosciences*, 107, 1-19.
.. [3] Haber, E. (2015). *Computational Methods in Geophysical
       Electromagnetics*. SIAM.
"""
