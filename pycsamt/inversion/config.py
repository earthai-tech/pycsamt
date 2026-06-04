# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Configuration for physics-based EM inversion workflows."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ..api.property import MetadataMixin, PyCSAMTObject
from ..models.config_io import (
    ConfigParameter,
    read_config_file,
    write_config_template,
)

__all__ = ["InversionConfig"]


_SCHEMA = [
    ConfigParameter(
        "method",
        "EM method to invert. Accepted values are 'mt', 'amt', 'csamt', "
        "'emap', and 'tdem'. The built-in runnable inversion targets "
        "MT/AMT/CSAMT impedance rho/phase data and TDEM decay data.",
        "Survey",
    ),
    ConfigParameter(
        "dimension",
        "Inversion dimensionality: '1d', '2d', or '3d'. Built-in paths "
        "cover layered 1-D and stitched station-by-station 2-D sections; "
        "external adapters cover Occam2D and ModEM workflows.",
        "Survey",
    ),
    ConfigParameter(
        "backend",
        "Backend name: 'builtin', 'occam2d', 'modem', 'simpeg', or "
        "'pygimli'. Optional backends are imported only when selected.",
        "Backend",
    ),
    ConfigParameter(
        "data",
        "Input data object, mapping, or path. Mappings can contain freqs, "
        "rho_a, phase, times, values, errors, station_names, and station_x.",
        "Survey",
    ),
    ConfigParameter(
        "workdir",
        "Working directory for generated files, logs, and backend outputs.",
        "Output",
    ),
    ConfigParameter(
        "output_dir",
        "Optional output directory for exported products. Defaults to workdir.",
        "Output",
    ),
    ConfigParameter(
        "n_layers",
        "Number of 1-D layers, including the bottom halfspace.",
        "Model",
    ),
    ConfigParameter(
        "starting_model",
        "Starting model mapping with resistivities and thicknesses. If null, "
        "a conservative default layered model is generated.",
        "Model",
    ),
    ConfigParameter(
        "reference_model",
        "Optional reference model used by backends that support reference "
        "regularization.",
        "Model",
    ),
    ConfigParameter(
        "error_floor",
        "Relative error floor used when explicit errors are unavailable.",
        "Data Errors",
    ),
    ConfigParameter(
        "include_phase",
        "Include phase observations in the built-in MT/AMT objective.",
        "Data Errors",
    ),
    ConfigParameter(
        "phase_error",
        "Absolute phase error in degrees used by the built-in objective.",
        "Data Errors",
    ),
    ConfigParameter(
        "regularization",
        "Regularization kind: 'none', 'smooth', 'damped', or 'blocky'.",
        "Regularization",
    ),
    ConfigParameter(
        "max_iter",
        "Maximum optimizer iterations for runnable local backends.",
        "Solver",
    ),
    ConfigParameter(
        "tol",
        "Optimizer termination tolerance.",
        "Solver",
    ),
    ConfigParameter(
        "bounds",
        "Optional parameter bounds. For the built-in 1-D backend use keys "
        "'log10_rho' and 'log10_thickness'.",
        "Solver",
    ),
    ConfigParameter(
        "run_external",
        "Whether external-code backends may launch their executable. Phase 1 "
        "defaults to preparation/loading only.",
        "Backend",
    ),
    ConfigParameter(
        "backend_options",
        "Backend-specific keyword options passed through unchanged.",
        "Backend",
    ),
    ConfigParameter(
        "metadata",
        "Free-form provenance metadata attached to results.",
        "Output",
    ),
]


@dataclass
class InversionConfig(PyCSAMTObject, MetadataMixin):
    """Source-of-truth configuration for :mod:`pycsamt.inversion`."""

    method: str = "mt"
    dimension: str = "1d"
    backend: str = "builtin"
    data: Any = None
    workdir: str = "pycsamt_inversion"
    output_dir: str | None = None
    n_layers: int = 4
    starting_model: Any = None
    reference_model: Any = None
    error_floor: float = 0.05
    include_phase: bool = True
    phase_error: float = 3.0
    regularization: str = "smooth"
    max_iter: int = 80
    tol: float = 1e-5
    bounds: dict[str, Any] = field(default_factory=dict)
    run_external: bool = False
    backend_options: dict[str, Any] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.method = str(self.method).lower()
        self.dimension = str(self.dimension).lower()
        self.backend = str(self.backend).lower()
        self.workdir = str(self.workdir)
        if self.output_dir is not None:
            self.output_dir = str(self.output_dir)
        self.n_layers = int(self.n_layers)
        self.error_floor = float(self.error_floor)
        self.include_phase = bool(self.include_phase)
        self.phase_error = float(self.phase_error)
        self.max_iter = int(self.max_iter)
        self.tol = float(self.tol)
        self.run_external = bool(self.run_external)
        self.bounds = dict(self.bounds or {})
        self.backend_options = dict(self.backend_options or {})
        self.metadata = dict(self.metadata or {})
        self.validate()

    @property
    def output_path(self) -> Path:
        """Resolved output directory."""
        return Path(self.output_dir or self.workdir)

    def validate(self) -> None:
        if self.method not in {"mt", "amt", "csamt", "emap", "tdem"}:
            raise ValueError("method must be mt/amt/csamt/emap/tdem.")
        if self.dimension not in {"1d", "2d", "3d"}:
            raise ValueError("dimension must be '1d', '2d', or '3d'.")
        if self.backend not in {"builtin", "occam2d", "modem", "simpeg", "pygimli"}:
            raise ValueError("unsupported inversion backend.")
        if self.n_layers < 2:
            raise ValueError("n_layers must be >= 2.")
        if self.error_floor < 0:
            raise ValueError("error_floor must be non-negative.")
        if self.phase_error <= 0:
            raise ValueError("phase_error must be positive.")
        if self.max_iter <= 0:
            raise ValueError("max_iter must be positive.")
        if self.tol <= 0:
            raise ValueError("tol must be positive.")

    def to_backend(self):
        """Instantiate the selected backend."""
        from .backends import get_backend

        return get_backend(self.backend)(self)

    def write_template(
        self,
        path: str | Path,
        *,
        fmt: str | None = None,
    ) -> Path:
        """Write an editable configuration template."""
        return write_config_template(
            path,
            self,
            _SCHEMA,
            fmt=fmt,
            title="PyCSAMT physics-based EM inversion configuration",
        )

    @classmethod
    def from_file(cls, path: str | Path, *, strict: bool = True) -> "InversionConfig":
        """Load a Python/JSON/YAML configuration file."""
        return cls(**read_config_file(path, cls, strict=strict))

    read = from_file
