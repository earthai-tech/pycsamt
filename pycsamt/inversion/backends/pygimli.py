# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pyGIMLi backend for EM inversion.

pyGIMLi is optional and imported only when this backend is selected. The
backend wraps ``pygimli.physics.em`` modelling operators and ``pg.Inversion``
for 1-D MT/AMT/CSAMT and TDEM soundings. Two-dimensional runs are represented
as stitched station-by-station 1-D inversions so their output can still be
used by PyCSAMT interpretation and plotting tools.

The implementation deliberately accepts several pyGIMLi operator names and
constructor signatures because the EM API has varied across pyGIMLi releases.
Backend options let users force a specific operator when automatic discovery
is not enough.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from ..base import BaseInversionBackend
from ..data import EMData
from ..mesh import InversionMesh
from ..model import StartingModel
from ..objective import component_errors, weighted_rms
from ..regularization import (
    pygimli_lambda,
    regularization_from_config,
)
from ..results import InversionResult

__all__ = ["PyGIMLiBackend"]


@dataclass
class _PyGIMLiModules:
    pg: Any
    em: Any


class PyGIMLiBackend(BaseInversionBackend):
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
        """Run a pyGIMLi-backed EM inversion.

        Parameters
        ----------
        data : mapping, object, sequence, or path-like, optional
            Optional data override for this call. When omitted, the backend
            uses ``self.config.data``. Values are coerced through
            :class:`pycsamt.inversion.data.EMData`.

        Returns
        -------
        InversionResult
            Backend-neutral result. For 1-D runs the result model is a
            :class:`pycsamt.inversion.model.StartingModel` containing recovered
            resistivities and thicknesses. For stitched 2-D profile runs the
            model is a dictionary with ``rho_2d``, station positions, and depth
            centres.

        Raises
        ------
        ImportError
            If ``pygimli`` is not installed.
        ValueError
            If the selected method does not receive the required data
            components. TDEM requires times and values. Natural-source methods
            require frequencies plus apparent resistivity and/or phase.
        NotImplementedError
            If the configured method/dimension pair is not supported, or if
            the installed pyGIMLi version does not expose a usable EM
            modelling operator.

        Examples
        --------
        >>> from pycsamt.inversion.backends.pygimli import PyGIMLiBackend
        >>> from pycsamt.inversion.config import InversionConfig
        >>> cfg = InversionConfig(
        ...     method="mt",
        ...     dimension="1d",
        ...     backend="pygimli",
        ...     data={"freqs": [1.0, 10.0],
        ...           "rho_a": [100.0, 120.0],
        ...           "phase": [45.0, 47.0]},
        ...     max_iter=8,
        ... )
        >>> result = PyGIMLiBackend(cfg).run()  # doctest: +SKIP
        """
        self.check_supported()
        modules = _load_pygimli()
        em_data = self.prepare_data(data)
        if self.config.method == "tdem":
            if not em_data.has_tdem_response:
                raise ValueError(
                    "pyGIMLi TDEM backend requires times plus values."
                )
        elif not em_data.has_mt_response:
            raise ValueError(
                "pyGIMLi MT backend requires frequencies plus rho_a and/or phase."
            )
        if self.config.dimension == "2d":
            return self._run_profile(em_data, modules)
        return self._run_sounding(em_data, modules, station_index=None)

    def _run_profile(
        self, em_data: EMData, modules: _PyGIMLiModules
    ) -> InversionResult:
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
                result = self._run_sounding(
                    _station_data(em_data, idx), modules, station_index=idx
                )
            except Exception as exc:
                warnings.append(
                    f"{names[idx]}: pyGIMLi inversion failed: {exc}"
                )
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
            objective=float(
                np.nansum([r.objective for r in station_results])
            ),
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
            return self._run_tdem_sounding(
                em_data, modules, station_index=station_index
            )
        return self._run_mt_sounding(
            em_data, modules, station_index=station_index
        )

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
        start = StartingModel.coerce(
            cfg.starting_model, n_layers=cfg.n_layers
        )
        periods = 1.0 / np.asarray(em_data.frequencies, dtype=float)
        np.asarray(start.thicknesses, dtype=float)
        observed = _pack_mt_observations(em_data)
        errors = _pack_mt_errors(em_data, cfg)

        opts = cfg.backend_options
        verbose = bool(opts.get("verbose", False))
        fop, start_model, operator_name, parameterization = (
            _make_mt_forward_operator(
                em,
                periods,
                start,
                opts,
            )
        )

        inv = _make_inversion(pg, fop, verbose=verbose)
        if hasattr(pg, "trans"):
            inv.modelTrans = pg.trans.TransLog()
        reg = regularization_from_config(cfg)
        lam = pygimli_lambda(cfg, default=20.0)
        recovered_raw = _run_inversion(
            inv,
            observed,
            startModel=start_model,
            errorVals=errors,
            lam=lam,
            maxIter=cfg.max_iter,
            verbose=verbose,
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
            extra={
                "lam": lam,
                "mode": "mt",
                "operator": operator_name,
                "parameterization": parameterization,
                "regularization": reg.to_dict(),
            },
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
        start = StartingModel.coerce(
            cfg.starting_model, n_layers=cfg.n_layers
        )
        thk = np.asarray(start.thicknesses, dtype=float)
        values = np.asarray(em_data.values, dtype=float)
        errors = _tdem_errors(em_data, cfg)
        opts = cfg.backend_options
        tx_area = float(
            opts.get("tx_area", opts.get("txArea", np.pi * 50.0**2))
        )
        rx_area = opts.get("rx_area", opts.get("rxArea", None))

        verbose = bool(opts.get("verbose", False))
        fop, operator_name = _make_tdem_forward_operator(
            em,
            np.asarray(em_data.times, dtype=float),
            thk,
            opts,
            tx_area=tx_area,
            rx_area=rx_area,
        )
        inv = _make_inversion(pg, fop, verbose=verbose)
        if hasattr(pg, "trans"):
            inv.dataTrans = pg.trans.TransLog()
            inv.modelTrans = pg.trans.TransLog()
        reg = regularization_from_config(cfg)
        lam = pygimli_lambda(cfg, default=20.0)
        recovered_raw = _run_inversion(
            inv,
            values,
            startModel=start.resistivities,
            errorVals=errors,
            lam=lam,
            maxIter=cfg.max_iter,
            verbose=verbose,
        )
        response = np.asarray(fop.response(recovered_raw), dtype=float)
        recovered = StartingModel(recovered_raw, thk, name="pygimli_tdem_1d")
        rms = weighted_rms(
            np.log10(np.maximum(np.abs(values), 1e-30)),
            np.log10(np.maximum(np.abs(response), 1e-30)),
            errors,
        )
        return _result_from_pygimli(
            cfg,
            em_data,
            recovered,
            response,
            rms,
            inv,
            fop,
            station_index,
            extra={
                "lam": lam,
                "mode": "tdem",
                "operator": operator_name,
                "tx_area": tx_area,
                "rx_area": rx_area,
                "regularization": reg.to_dict(),
            },
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


def _make_mt_forward_operator(
    em: Any,
    periods: np.ndarray,
    start: StartingModel,
    opts: dict[str, Any],
) -> tuple[Any, np.ndarray, str, str]:
    """Create a pyGIMLi MT forward operator across common API variants."""
    thk = np.asarray(start.thicknesses, dtype=float)
    verbose = bool(opts.get("verbose", False))
    candidates = _candidate_names(
        opts,
        "mt_operator",
        (
            "MT1dSmoothModelling",
            "MT1DSmoothModelling",
            "MT1dBlockModelling",
            "MT1DBlockModelling",
            "MT1dModelling",
            "MT1DModelling",
        ),
    )
    cls, operator_name = _resolve_operator(em, candidates, "MT/AMT/CSAMT 1-D")
    block = "block" in operator_name.lower()
    attempts = []
    if block:
        attempts.extend(
            [
                {"T": periods, "nLayers": start.n_layers, "verbose": verbose},
                {
                    "periods": periods,
                    "nLayers": start.n_layers,
                    "verbose": verbose,
                },
                {"T": periods, "nlay": start.n_layers, "verbose": verbose},
                {"T": periods, "thk": thk, "verbose": verbose},
            ]
        )
    else:
        attempts.extend(
            [
                {"T": periods, "thk": thk, "verbose": verbose},
                {"periods": periods, "thk": thk, "verbose": verbose},
                {"T": periods, "nLayers": start.n_layers, "verbose": verbose},
                {"T": periods, "verbose": verbose},
            ]
        )
    fop = _construct_operator(cls, attempts, operator_name)
    if block:
        start_model = np.r_[thk, start.resistivities]
        parameterization = "thickness_resistivity"
    else:
        start_model = start.resistivities
        parameterization = "resistivity"
    return (
        fop,
        np.asarray(start_model, dtype=float),
        operator_name,
        parameterization,
    )


def _make_tdem_forward_operator(
    em: Any,
    times: np.ndarray,
    thk: np.ndarray,
    opts: dict[str, Any],
    *,
    tx_area: float,
    rx_area: Any,
) -> tuple[Any, str]:
    """Create a pyGIMLi TDEM forward operator across common API variants."""
    verbose = bool(opts.get("verbose", False))
    candidates = _candidate_names(
        opts,
        "tdem_operator",
        (
            "TDEMSmoothModelling",
            "TDEMBlockModelling",
            "TDEM1dModelling",
            "TDEM1DModelling",
        ),
    )
    cls, operator_name = _resolve_operator(em, candidates, "TDEM 1-D")
    attempts = [
        {
            "thk": thk,
            "times": times,
            "txArea": tx_area,
            "rxArea": rx_area,
            "verbose": verbose,
        },
        {
            "thk": thk,
            "t": times,
            "txArea": tx_area,
            "rxArea": rx_area,
            "verbose": verbose,
        },
        {
            "thk": thk,
            "times": times,
            "tx_area": tx_area,
            "rx_area": rx_area,
            "verbose": verbose,
        },
        {"thk": thk, "times": times, "verbose": verbose},
    ]
    return _construct_operator(cls, attempts, operator_name), operator_name


def _candidate_names(
    opts: dict[str, Any],
    key: str,
    defaults: tuple[str, ...],
) -> tuple[str, ...]:
    configured = opts.get(key)
    if configured is None:
        return defaults
    if isinstance(configured, str):
        return (configured,)
    return tuple(str(name) for name in configured)


def _resolve_operator(
    em: Any,
    candidates: tuple[str, ...],
    label: str,
) -> tuple[Any, str]:
    for name in candidates:
        operator = getattr(em, name, None)
        if operator is not None:
            return operator, name
    available = sorted(
        name
        for name in dir(em)
        if "modelling" in name.lower()
        or name.lower().startswith(("mt", "tdem"))
    )
    raise NotImplementedError(
        f"pyGIMLi does not expose a supported {label} modelling operator. "
        f"Tried {', '.join(candidates)}. Available EM operators: "
        f"{', '.join(available) if available else 'none'}."
    )


def _construct_operator(
    cls: Any, attempts: list[dict[str, Any]], name: str
) -> Any:
    errors = []
    for kwargs in attempts:
        clean = {
            key: value for key, value in kwargs.items() if value is not None
        }
        try:
            return cls(**clean)
        except TypeError as exc:
            errors.append(f"{sorted(clean)}: {exc}")
    raise TypeError(
        f"could not construct pyGIMLi operator {name!r}; attempted signatures "
        + "; ".join(errors)
    )


def _make_inversion(pg: Any, fop: Any, *, verbose: bool) -> Any:
    try:
        return pg.Inversion(fop=fop, verbose=verbose)
    except TypeError:
        return pg.Inversion(fop, verbose=verbose)


def _run_inversion(
    inv: Any, observed: np.ndarray, **kwargs: Any
) -> np.ndarray:
    try:
        return np.asarray(inv.run(observed, **kwargs), dtype=float)
    except TypeError:
        fallback = dict(kwargs)
        if "errorVals" in fallback:
            fallback["relativeError"] = fallback.pop("errorVals")
        try:
            return np.asarray(inv.run(observed, **fallback), dtype=float)
        except TypeError:
            minimal = {
                key: value
                for key, value in fallback.items()
                if key in {"startModel", "lam", "maxIter", "verbose"}
            }
            return np.asarray(inv.run(observed, **minimal), dtype=float)


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
        errors.extend(
            component_errors(
                rho,
                cfg,
                component="rho",
                explicit=em_data.errors,
                relative=True,
            ).tolist()
        )
    if em_data.phase is not None:
        phase = np.asarray(em_data.phase, dtype=float)
        # pyGIMLi inversion error values are relative, so convert degrees
        # to a conservative relative phase error.
        errors.extend(
            component_errors(
                phase, cfg, component="phase", relative=True
            ).tolist()
        )
    return np.asarray(errors, dtype=float)


def _tdem_errors(em_data: EMData, cfg: Any) -> np.ndarray:
    return component_errors(
        em_data.values,
        cfg,
        component="tdem",
        explicit=em_data.errors,
        relative=True,
    )


def _recover_mt_model(raw: np.ndarray, start: StartingModel) -> StartingModel:
    raw = np.asarray(raw, dtype=float)
    if raw.size == start.n_layers:
        return StartingModel(raw, start.thicknesses, name="pygimli_mt_1d")
    if raw.size == 2 * start.n_layers - 1:
        thk = raw[: start.n_layers - 1]
        rho = raw[start.n_layers - 1 :]
        return StartingModel(rho, thk, name="pygimli_mt_1d")
    return StartingModel(
        raw[-start.n_layers :], start.thicknesses, name="pygimli_mt_1d"
    )


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
        objective=float(
            getattr(inv, "chi2", lambda: np.nan)()
            if callable(getattr(inv, "chi2", None))
            else getattr(inv, "chi2", np.nan)
        ),
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
        station_x=np.array(
            [_station_x(em_data, em_data.n_stations)[idx]], dtype=float
        ),
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


PyGIMLiBackend.__doc__ = r"""
Run optional pyGIMLi EM inversions.

``PyGIMLiBackend`` adapts
:mod:`pycsamt.inversion` configurations to pyGIMLi's electromagnetic
modelling and inversion APIs. It is an optional backend: importing
:mod:`pycsamt.inversion` does not require pyGIMLi, but selecting
``backend="pygimli"`` does.

The backend supports two execution patterns:

* 1-D MT/AMT/CSAMT soundings with apparent resistivity and optional phase;
* 1-D TDEM soundings with transient time gates and decay values.

For ``dimension="2d"``, the backend performs stitched station-by-station 1-D
inversions: each station is inverted independently with the same 1-D machinery
and the recovered log10 resistivity columns are assembled into a profile
model. This is useful for quick profile interpretation, but it is not a full
2-D finite-element pyGIMLi EM inversion.

Parameters
----------
config : pycsamt.inversion.config.InversionConfig
    Inversion configuration. Important fields are ``method``, ``dimension``,
    ``data``, ``starting_model``, ``regularization``, ``max_iter``, ``tol``,
    ``error_floor``, ``phase_error``, ``backend_options``, ``workdir``, and
    ``metadata``.

Attributes
----------
name : str
    Backend registry name, always ``"pygimli"``.
supports : tuple of tuple
    Supported ``(method, dimension)`` pairs.

Notes
-----
Natural-source observations are passed to pyGIMLi in the order apparent
resistivity followed by phase when both are present. pyGIMLi expects relative
error values, so PyCSAMT converts the shared component error model into the
form required by ``pg.Inversion``.

``backend_options`` can contain:

``mt_operator``
    Name or ordered sequence of names to try for MT/AMT/CSAMT modelling.
    Defaults include ``MT1dSmoothModelling``, ``MT1DSmoothModelling``,
    ``MT1dBlockModelling``, ``MT1DBlockModelling``, ``MT1dModelling``, and
    ``MT1DModelling``.
``tdem_operator``
    Name or ordered sequence of names to try for TDEM modelling. Defaults
    include ``TDEMSmoothModelling``, ``TDEMBlockModelling``,
    ``TDEM1dModelling``, and ``TDEM1DModelling``.
``lam``
    pyGIMLi regularization strength. If omitted, PyCSAMT maps the shared
    :class:`pycsamt.inversion.regularization.Regularization` settings through
    :func:`pycsamt.inversion.regularization.pygimli_lambda`.
``verbose``
    Forwarded to pyGIMLi operator and inversion construction when supported.
``tx_area`` / ``txArea`` and ``rx_area`` / ``rxArea``
    TDEM transmitter and receiver loop areas passed to the TDEM operator when
    supported by the installed pyGIMLi version.

Examples
--------
Run a 1-D MT sounding through the backend directly::

    >>> from pycsamt.inversion.backends.pygimli import PyGIMLiBackend
    >>> from pycsamt.inversion.config import InversionConfig
    >>> cfg = InversionConfig(
    ...     method="mt",
    ...     dimension="1d",
    ...     backend="pygimli",
    ...     data={"freqs": [1.0, 10.0],
    ...           "rho_a": [100.0, 120.0],
    ...           "phase": [45.0, 47.0]},
    ...     max_iter=8,
    ... )
    >>> result = PyGIMLiBackend(cfg).run()  # doctest: +SKIP

Run a TDEM sounding and pass loop geometry to pyGIMLi::

    >>> from pycsamt.inversion.workflow import run_inversion
    >>> result = run_inversion(
    ...     method="tdem",
    ...     dimension="1d",
    ...     backend="pygimli",
    ...     data={"times": [1e-5, 3e-5, 1e-4],
    ...           "values": [1e-8, 5e-9, 1e-9]},
    ...     backend_options={"tx_area": 7850.0, "rx_area": 100.0},
    ...     max_iter=8,
    ... )  # doctest: +SKIP

Build a stitched 2-D profile from station-wise MT inversions::

    >>> from pycsamt.inversion.workflow import run_inversion
    >>> result = run_inversion(
    ...     method="amt",
    ...     dimension="2d",
    ...     backend="pygimli",
    ...     data={"freqs": [10.0, 100.0],
    ...           "rho_a": [[80.0, 100.0], [90.0, 110.0]],
    ...           "phase": [[42.0, 45.0], [43.0, 46.0]],
    ...           "station_x": [0.0, 250.0],
    ...           "station_names": ["A01", "A02"]},
    ...     max_iter=6,
    ... )  # doctest: +SKIP
    >>> result.metadata["profile_mode"]  # doctest: +SKIP
    'stitched_station_1d'

Force a specific pyGIMLi operator name when needed::

    >>> from pycsamt.inversion.workflow import run_inversion
    >>> result = run_inversion(
    ...     method="mt",
    ...     dimension="1d",
    ...     backend="pygimli",
    ...     data={"freqs": [1.0], "rho_a": [100.0]},
    ...     backend_options={
    ...         "mt_operator": "MT1DModelling",
    ...         "lam": 15.0,
    ...         "verbose": False,
    ...     },
    ... )  # doctest: +SKIP

See Also
--------
pycsamt.inversion.workflow.InversionWorkflow
    High-level entry point that instantiates this backend.
pycsamt.inversion.backends.builtin.Builtin1DBackend
    Dependency-light fallback backend for local smoke inversions.
pycsamt.inversion.backends.simpeg.SimPEGBackend
    Optional SimPEG backend for natural-source physics inversion.
pycsamt.inversion.regularization.pygimli_lambda
    Helper that maps shared regularization settings to pyGIMLi ``lam``.
pycsamt.inversion.results.InversionResult
    Backend-neutral result returned by :meth:`run`.

References
----------
.. [1] Rücker, C., Günther, T. and Wagner, F. M. (2017). pyGIMLi: An
       open-source library for modelling and inversion in geophysics.
       *Computers & Geosciences*, 109, 106-123.
.. [2] Aster, R. C., Borchers, B. and Thurber, C. H. (2018). *Parameter
       Estimation and Inverse Problems*, 3rd edition. Elsevier.
.. [3] Chave, A. D. and Jones, A. G. (2012). *The Magnetotelluric Method:
       Theory and Practice*. Cambridge University Press.
"""
