# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pyGIMLi backend for 1-D EM inversion.

pyGIMLi is optional. This backend imports it only when selected. The first
supported targets are 1-D MT/AMT/CSAMT and TDEM inversions using
``pygimli.physics.em`` modelling operators and ``pg.Inversion``.
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

__all__ = ["PyGIMLiBackend"]


@dataclass
class _PyGIMLiModules:
    pg: Any
    em: Any


class PyGIMLiBackend(BaseInversionBackend):
    """Run optional pyGIMLi 1-D EM inversions."""

    name = "pygimli"
    supports = (
        ("mt", "1d"),
        ("mt", "2d"),
        ("amt", "1d"),
        ("amt", "2d"),
        ("csamt", "1d"),
        ("csamt", "2d"),
        ("tdem", "1d"),
        ("tdem", "2d"),
    )

    def run(self, data: Any | None = None) -> InversionResult:
        self.check_supported()
        modules = _load_pygimli()
        em_data = self.prepare_data(data)
        if self.config.method == "tdem":
            if not em_data.has_tdem_response:
                raise ValueError("pyGIMLi TDEM backend requires times plus values.")
        elif not em_data.has_mt_response:
            raise ValueError(
                "pyGIMLi MT backend requires frequencies plus rho_a and/or phase."
            )
        if self.config.dimension == "2d":
            return self._run_profile(em_data, modules)
        return self._run_sounding(em_data, modules, station_index=None)

    def _run_profile(self, em_data: EMData, modules: _PyGIMLiModules) -> InversionResult:
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
            try:
                result = self._run_sounding(_station_data(em_data, idx), modules, station_index=idx)
            except Exception as exc:
                warnings.append(f"{names[idx]}: pyGIMLi inversion failed: {exc}")
                continue
            station_results.append(result)
            used.append(idx)
            columns.append(np.log10(result.model.resistivities))
            if z_centers is None and result.mesh is not None:
                z_centers = result.mesh.z_centers

        if not columns:
            raise RuntimeError("all pyGIMLi station inversions failed.")

        rho_2d = np.stack(columns, axis=1)
        used_x = xs[used]
        used_names = [names[i] for i in used]
        if z_centers is None:
            z_centers = np.arange(rho_2d.shape[0], dtype=float)
        mesh = InversionMesh(
            dimension="2d",
            x_centers=used_x,
            z_centers=z_centers,
            metadata={"engine": "pygimli", "profile_mode": "stitched_1d"},
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
                "engine": "pygimli",
                "profile_mode": "stitched_station_1d",
                "station_rms": rms_values.tolist(),
            },
        )

    def _run_sounding(
        self,
        em_data: EMData,
        modules: _PyGIMLiModules,
        *,
        station_index: int | None,
    ) -> InversionResult:
        if self.config.method == "tdem":
            return self._run_tdem_sounding(em_data, modules, station_index=station_index)
        return self._run_mt_sounding(em_data, modules, station_index=station_index)

    def _run_mt_sounding(
        self,
        em_data: EMData,
        modules: _PyGIMLiModules,
        *,
        station_index: int | None,
    ) -> InversionResult:
        cfg = self.config
        pg = modules.pg
        em = modules.em
        start = StartingModel.coerce(cfg.starting_model, n_layers=cfg.n_layers)
        periods = 1.0 / np.asarray(em_data.frequencies, dtype=float)
        thk = np.asarray(start.thicknesses, dtype=float)
        observed = _pack_mt_observations(em_data)
        errors = _pack_mt_errors(em_data, cfg)

        fop_cls = getattr(em, "MT1dSmoothModelling", None)
        if fop_cls is None:
            fop_cls = getattr(em, "MT1dBlockModelling")
            fop = fop_cls(T=periods, nLayers=start.n_layers, verbose=False)
            start_model = np.r_[thk, start.resistivities]
        else:
            fop = fop_cls(T=periods, thk=thk, verbose=False)
            start_model = start.resistivities

        inv = pg.Inversion(fop=fop, verbose=bool(cfg.backend_options.get("verbose", False)))
        if hasattr(pg, "trans"):
            inv.modelTrans = pg.trans.TransLog()
        lam = float(cfg.backend_options.get("lam", cfg.backend_options.get("lambda", 20.0)))
        recovered_raw = inv.run(
            observed,
            startModel=start_model,
            errorVals=errors,
            lam=lam,
            maxIter=cfg.max_iter,
            verbose=bool(cfg.backend_options.get("verbose", False)),
        )
        response = np.asarray(fop.response(recovered_raw), dtype=float)
        recovered = _recover_mt_model(recovered_raw, start)
        rms = weighted_rms(observed, response, errors)
        return _result_from_pygimli(
            cfg,
            em_data,
            recovered,
            response,
            rms,
            inv,
            fop,
            station_index,
            extra={"lam": lam, "mode": "mt"},
        )

    def _run_tdem_sounding(
        self,
        em_data: EMData,
        modules: _PyGIMLiModules,
        *,
        station_index: int | None,
    ) -> InversionResult:
        cfg = self.config
        pg = modules.pg
        em = modules.em
        start = StartingModel.coerce(cfg.starting_model, n_layers=cfg.n_layers)
        thk = np.asarray(start.thicknesses, dtype=float)
        values = np.asarray(em_data.values, dtype=float)
        errors = _tdem_errors(em_data, cfg)
        opts = cfg.backend_options
        tx_area = float(opts.get("tx_area", opts.get("txArea", np.pi * 50.0 ** 2)))
        rx_area = opts.get("rx_area", opts.get("rxArea", None))

        fop = em.TDEMSmoothModelling(
            thk=thk,
            times=np.asarray(em_data.times, dtype=float),
            txArea=tx_area,
            rxArea=rx_area,
        )
        inv = pg.Inversion(fop=fop, verbose=bool(opts.get("verbose", False)))
        if hasattr(pg, "trans"):
            inv.dataTrans = pg.trans.TransLog()
            inv.modelTrans = pg.trans.TransLog()
        lam = float(opts.get("lam", opts.get("lambda", 20.0)))
        recovered_raw = inv.run(
            values,
            startModel=start.resistivities,
            errorVals=errors,
            lam=lam,
            maxIter=cfg.max_iter,
            verbose=bool(opts.get("verbose", False)),
        )
        response = np.asarray(fop.response(recovered_raw), dtype=float)
        recovered = StartingModel(recovered_raw, thk, name="pygimli_tdem_1d")
        rms = weighted_rms(np.log10(np.maximum(np.abs(values), 1e-30)),
                           np.log10(np.maximum(np.abs(response), 1e-30)),
                           errors)
        return _result_from_pygimli(
            cfg,
            em_data,
            recovered,
            response,
            rms,
            inv,
            fop,
            station_index,
            extra={"lam": lam, "mode": "tdem", "tx_area": tx_area, "rx_area": rx_area},
        )


def _load_pygimli() -> _PyGIMLiModules:
    try:
        import pygimli as pg
        from pygimli.physics import em
    except ImportError as exc:
        raise ImportError(
            "pyGIMLi backend selected, but pygimli is not installed. "
            "Install pyGIMLi, or choose backend='builtin'/'simpeg'."
        ) from exc
    return _PyGIMLiModules(pg=pg, em=em)


def _pack_mt_observations(em_data: EMData) -> np.ndarray:
    values = []
    if em_data.rho_a is not None:
        values.extend(np.asarray(em_data.rho_a, dtype=float).tolist())
    if em_data.phase is not None:
        values.extend(np.asarray(em_data.phase, dtype=float).tolist())
    return np.asarray(values, dtype=float)


def _pack_mt_errors(em_data: EMData, cfg: Any) -> np.ndarray:
    errors = []
    if em_data.rho_a is not None:
        rho = np.asarray(em_data.rho_a, dtype=float)
        if em_data.errors is not None:
            errors.extend(np.asarray(em_data.errors, dtype=float).tolist())
        else:
            errors.extend(np.full(rho.shape, cfg.error_floor, dtype=float).tolist())
    if em_data.phase is not None:
        phase = np.asarray(em_data.phase, dtype=float)
        # pyGIMLi inversion error values are relative, so convert degrees
        # to a conservative relative phase error.
        phase_rel = np.maximum(cfg.phase_error / np.maximum(np.abs(phase), 1.0), 1e-3)
        errors.extend(phase_rel.tolist())
    return np.asarray(errors, dtype=float)


def _tdem_errors(em_data: EMData, cfg: Any) -> np.ndarray:
    if em_data.errors is not None:
        raw = np.asarray(em_data.errors, dtype=float)
        values = np.asarray(em_data.values, dtype=float)
        return np.maximum(raw / np.maximum(np.abs(values), 1e-30), 1e-3)
    return np.full_like(np.asarray(em_data.values, dtype=float), cfg.error_floor, dtype=float)


def _recover_mt_model(raw: np.ndarray, start: StartingModel) -> StartingModel:
    raw = np.asarray(raw, dtype=float)
    if raw.size == start.n_layers:
        return StartingModel(raw, start.thicknesses, name="pygimli_mt_1d")
    if raw.size == 2 * start.n_layers - 1:
        thk = raw[: start.n_layers - 1]
        rho = raw[start.n_layers - 1:]
        return StartingModel(rho, thk, name="pygimli_mt_1d")
    return StartingModel(raw[-start.n_layers:], start.thicknesses, name="pygimli_mt_1d")


def _result_from_pygimli(
    cfg: Any,
    em_data: EMData,
    recovered: StartingModel,
    response: np.ndarray,
    rms: float,
    inv: Any,
    fop: Any,
    station_index: int | None,
    *,
    extra: dict[str, Any],
) -> InversionResult:
    mesh = InversionMesh.for_1d(_layer_centers(recovered.thicknesses))
    mesh.metadata.update({"engine": "pygimli"})
    return InversionResult(
        method=cfg.method,
        dimension="1d",
        backend="pygimli",
        status="success",
        model=recovered,
        mesh=mesh,
        data=em_data,
        predicted=response,
        rms=rms,
        objective=float(getattr(inv, "chi2", lambda: np.nan)()
                        if callable(getattr(inv, "chi2", None))
                        else getattr(inv, "chi2", np.nan)),
        n_iter=int(getattr(inv, "iter", 0)),
        workdir=cfg.workdir,
        native={"inversion": inv, "fop": fop},
        metadata={
            **cfg.metadata,
            "engine": "pygimli",
            "station_index": station_index,
            **extra,
        },
    )


def _layer_centers(thicknesses: np.ndarray) -> np.ndarray:
    tops = np.r_[0.0, np.cumsum(thicknesses)]
    last = thicknesses[-1] if thicknesses.size else 1.0
    bottoms = np.r_[tops[1:], tops[-1] + last]
    return 0.5 * (tops + bottoms)


def _station_data(em_data: EMData, idx: int) -> EMData:
    return EMData(
        method=em_data.method,
        frequencies=em_data.frequencies,
        times=em_data.times,
        rho_a=None if em_data.rho_a is None else _row(em_data.rho_a, idx),
        phase=None if em_data.phase is None else _row(em_data.phase, idx),
        values=None if em_data.values is None else _row(em_data.values, idx),
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
