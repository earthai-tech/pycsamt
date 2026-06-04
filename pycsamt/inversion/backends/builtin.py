# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Built-in lightweight 1-D EM inversion backend."""

from __future__ import annotations

from typing import Any

import numpy as np

from ..base import BaseInversionBackend
from ..mesh import InversionMesh
from ..model import StartingModel
from ..objective import weighted_rms
from ..results import InversionResult

__all__ = ["Builtin1DBackend"]


class Builtin1DBackend(BaseInversionBackend):
    """Small SciPy least-squares inversion for 1-D soundings and profiles."""

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

        freqs = np.asarray(em_data.frequencies, dtype=float)
        obs_parts: list[np.ndarray] = []
        err_parts: list[np.ndarray] = []

        if em_data.rho_a is not None:
            rho_obs = np.log10(np.maximum(em_data.rho_a, 1e-12))
            obs_parts.append(rho_obs)
            if em_data.errors is not None:
                err = np.asarray(em_data.errors, dtype=float)
                err_parts.append(np.maximum(err / np.maximum(em_data.rho_a, 1e-12), 1e-6))
            else:
                err_parts.append(np.full_like(rho_obs, max(cfg.error_floor, 1e-6)))
        if cfg.include_phase and em_data.phase is not None:
            obs_parts.append(np.asarray(em_data.phase, dtype=float))
            err_parts.append(np.full_like(em_data.phase, cfg.phase_error, dtype=float))

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
            return (predicted - observed) / errors

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
            warnings=warnings,
            metadata={
                **cfg.metadata,
                "optimizer_message": str(opt.message),
                "starting_model": start.to_dict(),
                "station_index": station_index,
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
        times = np.asarray(em_data.times, dtype=float)
        observed = _log_abs(np.asarray(em_data.values, dtype=float))
        if em_data.errors is not None:
            raw_err = np.asarray(em_data.errors, dtype=float)
            errors = np.maximum(
                raw_err / np.maximum(np.abs(em_data.values), 1e-30),
                1e-6,
            )
        else:
            errors = np.full_like(observed, max(cfg.error_floor, 1e-6))

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
            return (predicted - observed) / errors

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
            warnings=warnings,
            metadata={
                **cfg.metadata,
                "optimizer_message": str(opt.message),
                "starting_model": start.to_dict(),
                "station_index": station_index,
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
