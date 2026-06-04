# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""SimPEG backend for physics-based EM inversion.

The implementation keeps SimPEG optional and imports it only when this
backend is selected. The supported first target is natural-source 1-D
MT/AMT/CSAMT inversion using SimPEG's ``natural_source`` module. Profile
data are handled as stitched station-by-station 1-D inversions so the
result still converts into :class:`pycsamt.interp.ResistivityModel`.

The code follows the documented SimPEG natural-source API: 1-D simulations
use ``Simulation1DElectricField``; MT observations are represented by
``PointNaturalSource`` receivers with ``apparent_resistivity`` and
``phase`` components; frequencies are represented by ``Planewave`` sources.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from ..base import BaseInversionBackend
from ..data import EMData
from ..mesh import InversionMesh
from ..model import StartingModel
from ..objective import weighted_rms
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
    """Run optional SimPEG natural-source EM inversions."""

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
    n_cells = int(options.get("n_cells", max(32, start.n_layers * 12)))
    depth_max = float(options.get(
        "depth_max",
        max(float(np.sum(start.thicknesses)) * 3.0, float(start.thicknesses[-1]) * 4.0),
    ))
    min_cell = float(options.get("min_cell_size", max(depth_max / n_cells / 4.0, 1.0)))
    growth = float(options.get("growth_factor", 1.08))
    widths = min_cell * growth ** np.arange(n_cells, dtype=float)
    widths *= depth_max / np.sum(widths)
    mesh = modules.discretize.TensorMesh([widths], origin="0")
    z_centers = np.cumsum(widths) - 0.5 * widths
    return mesh, z_centers


def _build_3d_mesh(
    em_data: EMData,
    options: dict[str, Any],
    modules: _SimPEGModules,
):
    station_x = _station_x(em_data, em_data.n_stations)
    station_y = _station_y(em_data, em_data.n_stations)
    x_pad = float(options.get("x_pad", 1000.0))
    y_pad = float(options.get("y_pad", 1000.0))
    depth_max = float(options.get("depth_max", 3000.0))
    nx = int(options.get("nx", max(8, min(32, em_data.n_stations * 4))))
    ny = int(options.get("ny", max(8, min(32, em_data.n_stations * 4))))
    nz = int(options.get("nz", 16))
    x_min, x_max = float(np.min(station_x) - x_pad), float(np.max(station_x) + x_pad)
    y_min, y_max = float(np.min(station_y) - y_pad), float(np.max(station_y) + y_pad)
    hx = np.full(nx, (x_max - x_min) / max(nx, 1), dtype=float)
    hy = np.full(ny, (y_max - y_min) / max(ny, 1), dtype=float)
    hz = _depth_widths(depth_max, nz, options)
    origin = np.array([x_min, y_min, -float(np.sum(hz))], dtype=float)
    mesh = modules.discretize.TensorMesh([hx, hy, hz], origin=origin)
    centers = {
        "x": x_min + np.cumsum(hx) - 0.5 * hx,
        "y": y_min + np.cumsum(hy) - 0.5 * hy,
        "z": origin[2] + np.cumsum(hz) - 0.5 * hz,
    }
    centers["z_depth"] = -centers["z"]
    return mesh, centers


def _depth_widths(depth_max: float, nz: int, options: dict[str, Any]) -> np.ndarray:
    min_cell = float(options.get("min_cell_size", max(depth_max / nz / 4.0, 1.0)))
    growth = float(options.get("growth_factor", 1.08))
    widths = min_cell * growth ** np.arange(nz, dtype=float)
    widths *= depth_max / np.sum(widths)
    return widths


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
            values.extend(np.asarray(rho_i, dtype=float).reshape(-1).tolist())
            if raw_err is not None:
                err_i = raw_err[i] if raw_err.ndim == 1 else raw_err[:, i]
                floor_i = np.abs(rho_i) * cfg.error_floor
                errors.extend(np.maximum(err_i, np.maximum(floor_i, 1e-12)).reshape(-1).tolist())
            else:
                errors.extend(np.maximum(np.abs(rho_i) * cfg.error_floor, 1e-12).reshape(-1).tolist())
        if phase is not None:
            phase_i = phase[i] if phase.ndim == 1 else phase[:, i]
            phase_i = np.asarray(phase_i, dtype=float).reshape(-1)
            values.extend(phase_i.tolist())
            errors.extend(np.full(phase_i.shape, cfg.phase_error, dtype=float).tolist())
    return np.asarray(values, dtype=float), np.asarray(errors, dtype=float)


def _build_regularization(
    mesh: Any,
    mapping: Any,
    modules: _SimPEGModules,
    cfg: Any,
):
    kind = str(cfg.regularization).lower()
    if hasattr(modules.regularization, "WeightedLeastSquares"):
        reg = modules.regularization.WeightedLeastSquares(mesh, mapping=mapping)
    else:
        reg = modules.regularization.Simple(mesh, mapping=mapping)
    if kind == "none":
        reg.alpha_s = 0.0
        if hasattr(reg, "alpha_x"):
            reg.alpha_x = 0.0
    elif kind == "damped":
        reg.alpha_s = float(cfg.backend_options.get("alpha_s", 1.0))
        if hasattr(reg, "alpha_x"):
            reg.alpha_x = 0.0
    else:
        reg.alpha_s = float(cfg.backend_options.get("alpha_s", 1.0))
        if hasattr(reg, "alpha_x"):
            reg.alpha_x = float(cfg.backend_options.get("alpha_x", 1.0))
    return reg


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
