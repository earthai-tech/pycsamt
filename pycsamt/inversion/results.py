# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Unified inversion result container."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ..api.property import MetadataMixin, PyCSAMTObject

from .doc import _inversion_param_docs
from .mesh import InversionMesh

__all__ = ["InversionHistory", "InversionResult", "InversionUncertainty"]


@dataclass
class InversionHistory(PyCSAMTObject, MetadataMixin):
    records: list[dict[str, Any]] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.records = [dict(record) for record in self.records]
        self.metadata = dict(self.metadata or {})

    @property
    def n_records(self) -> int:
        """Number of convergence records.

        Examples
        --------
        >>> from pycsamt.inversion.results import InversionHistory
        >>> InversionHistory(records=[{"iteration": 0}]).n_records
        1
        """
        return len(self.records)

    def arrays(self) -> dict[str, np.ndarray]:
        """Return numeric history fields as arrays.

        Non-numeric fields are skipped. Missing numeric fields are represented
        by ``nan`` so arrays have one value per record.

        Returns
        -------
        dict of ndarray
            Mapping from numeric history field name to a float array.

        Examples
        --------
        >>> from pycsamt.inversion.results import InversionHistory
        >>> history = InversionHistory(records=[
        ...     {"iteration": 0, "objective": 5.0},
        ...     {"iteration": 1, "objective": 2.5},
        ... ])
        >>> history.arrays()["objective"].tolist()
        [5.0, 2.5]
        """
        keys: set[str] = set()
        for record in self.records:
            keys.update(record)
        out: dict[str, np.ndarray] = {}
        for key in sorted(keys):
            values = [record.get(key, np.nan) for record in self.records]
            try:
                out[key] = np.asarray(values, dtype=float)
            except (TypeError, ValueError):
                continue
        return out


@dataclass
class InversionUncertainty(PyCSAMTObject, MetadataMixin):
    model_std: Any = None
    covariance_diag: Any = None
    sensitivity: Any = None
    confidence: Any = None
    station_confidence: Any = None
    depth_confidence: Any = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.model_std = _array_or_none(self.model_std)
        self.covariance_diag = _array_or_none(self.covariance_diag)
        self.sensitivity = _array_or_none(self.sensitivity)
        self.confidence = _array_or_none(self.confidence)
        self.station_confidence = _array_or_none(self.station_confidence)
        self.depth_confidence = _array_or_none(self.depth_confidence)
        self.metadata = dict(self.metadata or {})


@dataclass
class InversionResult(PyCSAMTObject, MetadataMixin):
    method: str
    dimension: str
    backend: str
    status: str = "success"
    model: Any = None
    mesh: InversionMesh | None = None
    data: Any = None
    predicted: Any = None
    rms: float = float("nan")
    objective: float = float("nan")
    n_iter: int = 0
    workdir: str | None = None
    files: dict[str, str] = field(default_factory=dict)
    native: Any = field(default=None, repr=False)
    uncertainty: InversionUncertainty | None = None
    history: InversionHistory | None = None
    warnings: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.method = str(self.method).lower()
        self.dimension = str(self.dimension).lower()
        self.backend = str(self.backend).lower()
        self.status = str(self.status)
        self.rms = float(self.rms)
        self.objective = float(self.objective)
        self.n_iter = int(self.n_iter)
        self.files = dict(self.files or {})
        if isinstance(self.uncertainty, dict):
            self.uncertainty = InversionUncertainty(**self.uncertainty)
        if isinstance(self.history, dict):
            self.history = InversionHistory(**self.history)
        elif isinstance(self.history, list):
            self.history = InversionHistory(records=self.history)
        self.warnings = [str(w) for w in self.warnings]
        self.metadata = dict(self.metadata or {})

    @property
    def converged(self) -> bool:
        """Whether the backend reported a usable result.

        Returns
        -------
        bool
            ``True`` for statuses ``"success"``, ``"converged"``,
            ``"prepared"``, and ``"loaded"``.

        Examples
        --------
        >>> from pycsamt.inversion.results import InversionResult
        >>> InversionResult("mt", "1d", "builtin", status="converged").converged
        True
        """
        return self.status in {"success", "converged", "prepared", "loaded"}

    def to_resistivity_model(self):
        """Convert the result to :class:`pycsamt.interp.ResistivityModel`.

        2-D result arrays are passed through directly. A recovered 1-D
        layered model is expanded to one column so the interpretation API can
        consume it uniformly.

        Returns
        -------
        pycsamt.interp.ResistivityModel
            Unified 2-D log10 resistivity model with profile and depth
            coordinates.

        Raises
        ------
        ValueError
            If no model is attached or if the attached model cannot be adapted
            to a resistivity section.

        Notes
        -----
        The conversion accepts three result forms:

        * backend-native Occam-like objects exposing ``rho_2d``;
        * dictionaries containing ``rho_2d`` plus coordinate arrays;
        * layered models exposing resistivities and thicknesses.

        Examples
        --------
        >>> import numpy as np
        >>> from pycsamt.inversion.results import InversionResult
        >>> result = InversionResult(
        ...     method="mt",
        ...     dimension="2d",
        ...     backend="builtin",
        ...     model={
        ...         "rho_2d": np.array([[2.0, 2.2]]),
        ...         "x_centers": np.array([0.0, 500.0]),
        ...         "z_centers": np.array([100.0]),
        ...     },
        ... )
        >>> result.to_resistivity_model().n_x
        2
        """
        from ..interp import ResistivityModel

        if hasattr(self.native, "rho_2d"):
            try:
                return ResistivityModel.from_occam2d(self.native)
            except Exception:
                pass
        if isinstance(self.model, dict) and "rho_2d" in self.model:
            rho = np.asarray(self.model["rho_2d"], dtype=float)
            x = np.asarray(self.model.get("x_centers", np.arange(rho.shape[1])), dtype=float)
            z = np.asarray(self.model.get("z_centers", np.arange(rho.shape[0])), dtype=float)
            return ResistivityModel.from_array(
                rho,
                x,
                z,
                station_x=np.asarray(self.model.get("station_x", x), dtype=float),
                station_names=list(self.model.get("station_names", [])) or None,
                method=f"{self.backend}:{self.method}",
                rms=self.rms,
            )
        if self.model is None:
            raise ValueError("result has no model to convert.")

        resistivity = np.asarray(
            getattr(self.model, "resistivities", getattr(self.model, "resistivity", [])),
            dtype=float,
        )
        thickness = np.asarray(
            getattr(self.model, "thicknesses", getattr(self.model, "thickness", [])),
            dtype=float,
        )
        if resistivity.size == 0 or thickness.size != resistivity.size - 1:
            raise ValueError("cannot convert model to ResistivityModel.")

        tops = np.r_[0.0, np.cumsum(thickness)]
        bottoms = np.r_[tops[1:], tops[-1] + thickness[-1]]
        z = 0.5 * (tops + bottoms)
        rho_2d = np.log10(resistivity).reshape(-1, 1)
        return ResistivityModel.from_array(
            rho_2d,
            np.array([0.0]),
            z,
            station_x=np.array([0.0]),
            station_names=["S000"],
            method=f"{self.backend}:{self.method}:1d",
            rms=self.rms,
        )

    def summary(self, *, max_fields: int | None = None) -> str:
        """Return a compact one-line result summary.

        Parameters
        ----------
        max_fields : int, optional
            Reserved for future expanded summaries. Currently ignored.

        Returns
        -------
        str
            Human-readable summary containing method, dimension, backend,
            status, and RMS.

        Examples
        --------
        >>> from pycsamt.inversion.results import InversionResult
        >>> InversionResult("mt", "1d", "builtin", rms=1.2).summary()
        "InversionResult(method='mt', dimension='1d', backend='builtin', status='success', rms=1.2)"
        """
        status = self.status
        rms = "nan" if not np.isfinite(self.rms) else f"{self.rms:.3g}"
        return (
            f"InversionResult(method={self.method!r}, dimension={self.dimension!r}, "
            f"backend={self.backend!r}, status={status!r}, rms={rms})"
        )


InversionHistory.__doc__ = rf"""
Common convergence-history container.

``InversionHistory`` stores backend-neutral iteration diagnostics. Built-in
solvers record fields such as ``iteration``, ``phi_d``, ``phi_m``,
``objective``, ``rms``, regularization weight, and model norm. Optional or
external backends can map native logs into the same record structure.

Parameters
----------
records : list of dict, optional
    One dictionary per iteration or residual evaluation. Numeric fields can be
    converted to arrays with :meth:`arrays`.
metadata : dict, optional
    Free-form provenance metadata such as backend mode or station index.

{_inversion_param_docs.result.history_examples}

See Also
--------
InversionResult.history
    Result field that carries convergence diagnostics.

{_inversion_param_docs.result.references}
"""

InversionUncertainty.__doc__ = rf"""
Backend-neutral uncertainty and sensitivity diagnostics.

``InversionUncertainty`` stores uncertainty products in a common shape for
exports, plots, and interpretation. Built-in solvers currently provide proxy
diagnostics from least-squares Jacobians; optional engines may attach richer
posterior or sensitivity products over time.

Parameters
----------
model_std : array-like, optional
    Per-model-cell standard-deviation proxy.
covariance_diag : array-like, optional
    Diagonal covariance or covariance-like proxy.
sensitivity : array-like, optional
    Relative model sensitivity map.
confidence : array-like, optional
    Normalized confidence map, typically scaled to ``[0, 1]``.
station_confidence : array-like, optional
    One confidence value per station.
depth_confidence : array-like, optional
    One confidence value per depth cell.
metadata : dict, optional
    Free-form uncertainty provenance.

{_inversion_param_docs.result.uncertainty_examples}

References
----------
.. [1] Aster, R. C., Borchers, B. and Thurber, C. H. (2018). *Parameter
   Estimation and Inverse Problems*, 3rd edition. Elsevier.
"""

InversionResult.__doc__ = rf"""
Backend-neutral post-inversion result.

``InversionResult`` is the common result object returned by every inversion
backend. It stores recovered models, predicted responses, misfit diagnostics,
uncertainty, convergence history, generated files, backend-native objects, and
free-form metadata while exposing conversion hooks for plotting, exporting, and
hydrogeophysical interpretation.

Parameters
----------
method : str
    EM method label such as ``"mt"``, ``"amt"``, ``"csamt"``, ``"emap"``, or
    ``"tdem"``.
dimension : {"1d", "2d", "3d"}
    Inversion dimensionality.
backend : str
    Backend name such as ``"builtin"``, ``"simpeg"``, ``"pygimli"``,
    ``"occam2d"``, or ``"modem"``.
status : str, default "success"
    Backend status. Common values include ``"success"``, ``"converged"``,
    ``"needs_review"``, ``"prepared"``, ``"ready"``, and ``"loaded"``.
model : object, optional
    Recovered model. Accepted conversion forms include a layered model or a
    dictionary with ``rho_2d`` and coordinate arrays.
mesh : InversionMesh, optional
    Common mesh descriptor.
data, predicted : object, optional
    Observed and predicted response containers.
rms : float, default nan
    Weighted RMS misfit.
objective : float, default nan
    Final objective-function value.
n_iter : int, default 0
    Number of backend iterations or residual evaluations.
workdir : str, optional
    Backend working directory.
files : dict, optional
    Generated or backend-native file paths.
native : object, optional
    Backend-native result object retained for advanced users.
uncertainty : InversionUncertainty or dict, optional
    Uncertainty diagnostics. Dictionaries are coerced automatically.
history : InversionHistory, dict, or list, optional
    Convergence history. Dictionaries and record lists are coerced
    automatically.
warnings : list of str, optional
    Backend warnings and lifecycle messages.
metadata : dict, optional
    Free-form result provenance.

{_inversion_param_docs.result.result_examples}

See Also
--------
pycsamt.inversion.export
    Export helpers that consume ``InversionResult``.
pycsamt.inversion.plot
    Plotting helpers that consume ``InversionResult``.
pycsamt.interp.ResistivityModel
    Interpretation model returned by :meth:`to_resistivity_model`.

{_inversion_param_docs.result.references}
"""


def _array_or_none(value: Any) -> np.ndarray | None:
    if value is None:
        return None
    return np.asarray(value, dtype=float)
