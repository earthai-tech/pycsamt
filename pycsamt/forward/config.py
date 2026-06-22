# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Configuration for pyCSAMT 1-D forward modelling and dataset generation.

The module exposes :class:`ForwardConfig`, a dataclass that collects every
tuneable parameter used by :func:`~pycsamt.forward.batch.generate_dataset`
and the three 1-D forward solvers.  Users can:

1. Call :meth:`ForwardConfig.write_template` to generate a fully annotated
   source-of-truth file (Python, JSON, or YAML).
2. Edit the file to match their survey setup and earth-model priors.
3. Load the edited file with :meth:`ForwardConfig.from_file` and pass the
   resulting object to :meth:`ForwardConfig.to_dataset_kwargs`.

Quick start
-----------
Generate a default template, edit it, run dataset generation::

    from pycsamt.forward.config import ForwardConfig

    # 1 — write an annotated source-of-truth file
    ForwardConfig.write_template("my_forward.py")

    # 2 — edit values in my_forward.py …

    # 3 — load and run
    cfg = ForwardConfig.from_file("my_forward.py")
    cfg.validate()

    from pycsamt.forward.batch import generate_dataset
    ds = generate_dataset(**cfg.to_dataset_kwargs())

Alternatively, construct the config entirely in Python::

    cfg = ForwardConfig(
        solver="mt1d",
        n_samples=20_000,
        n_layers_min=3,
        n_layers_max=7,
        noise_level=0.03,
        seed=42,
    )
    ds = generate_dataset(**cfg.to_dataset_kwargs())
"""
from __future__ import annotations

import os
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import Any, Dict, Optional, Tuple, Union

import numpy as np

from ..models.config_io import (
    ConfigParameter,
    read_config_file,
    write_config_template,
)

__all__ = ["ForwardConfig"]


# ── parameter schema (documentation shown in generated config files) ────────

_FORWARD_CONFIG_SCHEMA: list[ConfigParameter] = [
    # ── Solver ──────────────────────────────────────────────────────────────
    ConfigParameter(
        "solver",
        "Forward solver to use.  Accepted values are 'mt1d' (plane-wave "
        "magnetotelluric), 'csamt1d' (controlled-source AMT, far-field "
        "approximation), and 'tem1d' (central-loop time-domain EM with "
        "step-off waveform).  This value also governs which frequency/time "
        "grid parameters are active.",
        "Solver",
    ),
    ConfigParameter(
        "freq_min",
        "Minimum frequency in Hz for the MT/CSAMT frequency grid. "
        "The grid is constructed as numpy.logspace(log10(freq_min), "
        "log10(freq_max), n_freqs).  Typical AMT lower bound: 1e-4 Hz. "
        "Ignored when solver is 'tem1d'.",
        "Solver",
    ),
    ConfigParameter(
        "freq_max",
        "Maximum frequency in Hz for the MT/CSAMT frequency grid. "
        "Typical AMT upper bound: 1e5 Hz (CSAMT surveys reach 8192 Hz). "
        "Must be strictly greater than freq_min.  Ignored for TEM.",
        "Solver",
    ),
    ConfigParameter(
        "n_freqs",
        "Number of frequency points in the logarithmic frequency grid. "
        "30 points per decade gives a smooth apparent-resistivity curve "
        "for most network architectures.  Minimum meaningful value is 4.",
        "Solver",
    ),
    ConfigParameter(
        "time_min",
        "Minimum time in seconds for the TEM step-off time-gate grid. "
        "Constructed as numpy.logspace(log10(time_min), log10(time_max), "
        "n_times).  Typical TEM early gate: 1e-6 s.  Ignored for MT/CSAMT.",
        "Solver",
    ),
    ConfigParameter(
        "time_max",
        "Maximum time in seconds for the TEM time-gate grid.  Typical "
        "TEM late gate: 1e-2 s.  Must be strictly greater than time_min.",
        "Solver",
    ),
    ConfigParameter(
        "n_times",
        "Number of time-gate points in the logarithmic TEM grid. "
        "25 points is a good compromise between model accuracy and feature "
        "vector length for neural-network inversion.",
        "Solver",
    ),
    ConfigParameter(
        "loop_radius",
        "Transmitter loop radius in metres for TEM1D.  The Hankel-transform "
        "kernel uses this as the horizontal source extent.  Ignored for "
        "MT and CSAMT solvers.",
        "Solver",
    ),
    # ── Earth model ─────────────────────────────────────────────────────────
    ConfigParameter(
        "n_layers_min",
        "Minimum number of earth layers (including the halfspace).  When "
        "n_layers_min equals n_layers_max the layer count is fixed; "
        "otherwise it is drawn uniformly at random over [min, max] for "
        "each synthetic sample.  Minimum useful value is 2.",
        "Earth Model",
    ),
    ConfigParameter(
        "n_layers_max",
        "Maximum number of earth layers (including the halfspace).  "
        "Must be greater than or equal to n_layers_min.  Samples with "
        "fewer layers than the maximum are NaN-padded in the target vector.",
        "Earth Model",
    ),
    ConfigParameter(
        "rho_min",
        "Lower resistivity bound in Ω·m for the log-uniform prior used "
        "by LayeredModel.random.  Must be strictly positive.  Typical "
        "geophysical surveys cover 1–10 000 Ω·m.",
        "Earth Model",
    ),
    ConfigParameter(
        "rho_max",
        "Upper resistivity bound in Ω·m.  Must be strictly greater than "
        "rho_min.  Very large values (> 100 000) extend the dynamic range "
        "but reduce the density of conductive training samples.",
        "Earth Model",
    ),
    ConfigParameter(
        "depth_max",
        "Maximum depth in metres spanned by all layers except the "
        "halfspace.  Layer thicknesses are drawn from a Dirichlet-like "
        "distribution that sums to this value, so it controls the total "
        "investigation depth of each random model.",
        "Earth Model",
    ),
    ConfigParameter(
        "geology",
        "Optional geological prior scenario name passed to "
        "LayeredModel.from_geology.  When set, layer counts, resistivity "
        "bounds, and depth ranges are drawn from the named prior instead "
        "of the generic rho_min/rho_max/depth_max values.  "
        "Available scenarios: 'sedimentary', 'crystalline', 'geothermal', "
        "'marine', 'permafrost'.  Set to null/None to use the generic "
        "random prior.",
        "Earth Model",
    ),
    # ── Dataset ─────────────────────────────────────────────────────────────
    ConfigParameter(
        "n_samples",
        "Total number of (model, response) pairs to generate.  Larger "
        "datasets improve generalisation but increase generation time "
        "linearly.  10 000–50 000 samples is typical for 1-D inversion "
        "with five to seven layers.",
        "Dataset",
    ),
    ConfigParameter(
        "noise_level",
        "Relative noise standard deviation added to synthetic data after "
        "the forward solve.  A value of 0.05 means five percent Gaussian "
        "noise on each data point.  Set to 0.0 for noise-free data (not "
        "recommended for training: the network will not generalise to "
        "real field data).",
        "Dataset",
    ),
    ConfigParameter(
        "noise_type",
        "Noise model applied to synthetic responses.  Options: "
        "'gaussian' adds independent normal noise scaled by noise_level; "
        "'multiplicative' scales each point by (1 + epsilon) where "
        "epsilon ~ N(0, noise_level); 'field' mimics field data with "
        "frequency-dependent noise and occasional spikes.",
        "Dataset",
    ),
    ConfigParameter(
        "include_phase",
        "Include impedance phase in the MT/CSAMT feature vector.  "
        "When True the feature vector has length 2*n_freqs "
        "(log10(rho_a) concatenated with phase in degrees).  "
        "When False the feature vector has length n_freqs (rho_a only). "
        "Setting this to False approximately halves network input size "
        "at the cost of lost phase information.  Ignored for TEM.",
        "Dataset",
    ),
    ConfigParameter(
        "seed",
        "Global random seed for reproducible dataset generation.  "
        "Individual worker seeds are derived from this value, so the same "
        "seed always produces the same dataset regardless of n_jobs.  "
        "Set to null/None for a time-based random seed.",
        "Dataset",
    ),
    ConfigParameter(
        "n_jobs",
        "Number of parallel worker processes used by "
        "generate_dataset.  Use -1 to consume all available CPU cores. "
        "Values greater than 1 require the script to be protected by "
        "an 'if __name__ == \"__main__\":' guard on Windows.",
        "Dataset",
    ),
    # ── Output ──────────────────────────────────────────────────────────────
    ConfigParameter(
        "output_dir",
        "Directory where the generated dataset (.npz file) is written. "
        "The directory is created automatically if it does not exist. "
        "Set to '.' to save in the current working directory.  "
        "Set to null/None to skip saving to disk (dataset is returned "
        "from generate_dataset but not persisted).",
        "Output",
    ),
    ConfigParameter(
        "output_name",
        "Base file name for the saved dataset, without extension.  "
        "The extension '.npz' is added automatically.  For example, "
        "'mt1d_train' produces 'mt1d_train.npz' in output_dir.  "
        "Ignored when output_dir is null.",
        "Output",
    ),
    ConfigParameter(
        "verbose",
        "Print a progress bar and summary statistics while generating "
        "the dataset.  Set to False for batch scripts where only "
        "the final dataset matters.",
        "Output",
    ),
]


# ── dataclass ───────────────────────────────────────────────────────────────

@dataclass
class ForwardConfig:
    """Collect settings that define a 1-D forward modelling / dataset run.

    ``ForwardConfig`` is the central configuration object for
    :func:`~pycsamt.forward.batch.generate_dataset` and the three 1-D
    solvers.  All fields are plain Python scalars so that the dataclass can
    be serialised to Python, JSON, and YAML without any numpy dependencies.

    The recommended workflow is:

    1. Generate a template with :meth:`write_template`.
    2. Edit the values in the generated file.
    3. Load the edited file with :meth:`from_file`.
    4. Optionally call :meth:`validate` to catch range errors early.
    5. Pass :meth:`to_dataset_kwargs` to :func:`generate_dataset`.

    Parameters
    ----------
    solver : {'mt1d', 'csamt1d', 'tem1d'}
        Forward solver.
    freq_min, freq_max : float
        Frequency grid bounds [Hz] for MT/CSAMT.
    n_freqs : int
        Number of logarithmically spaced frequency points.
    time_min, time_max : float
        Time-gate grid bounds [s] for TEM.
    n_times : int
        Number of logarithmically spaced time-gate points.
    loop_radius : float
        TEM transmitter loop radius [m].
    n_layers_min, n_layers_max : int
        Layer count range (inclusive).  Equal values fix the count.
    rho_min, rho_max : float
        Resistivity bounds [Ω·m] for the random prior.
    depth_max : float
        Maximum total depth [m] across all layers except the halfspace.
    geology : str or None
        Geological prior scenario name; overrides rho/depth bounds when set.
    n_samples : int
        Number of (model, response) pairs to generate.
    noise_level : float
        Relative noise standard deviation (0 = noise-free).
    noise_type : str
        Noise model: ``'gaussian'``, ``'multiplicative'``, ``'field'``.
    include_phase : bool
        Include impedance phase in the MT/CSAMT feature vector.
    seed : int or None
        Global random seed for reproducibility.
    n_jobs : int
        Parallel worker count (-1 = all cores).
    output_dir : str or None
        Directory for the saved ``.npz`` dataset.  ``None`` skips saving.
    output_name : str
        Base file name (without extension) for the saved dataset.
    verbose : bool
        Print progress and summary.

    Examples
    --------
    Default MT1D configuration::

        >>> cfg = ForwardConfig()
        >>> cfg.solver
        'mt1d'

    Fixed 5-layer configuration with geology prior::

        >>> cfg = ForwardConfig(
        ...     solver="mt1d",
        ...     n_layers_min=5,
        ...     n_layers_max=5,
        ...     geology="sedimentary",
        ...     n_samples=10_000,
        ...     seed=0,
        ... )

    Generate a template file and round-trip::

        >>> path = ForwardConfig.write_template("forward_config.py")
        >>> cfg = ForwardConfig.from_file(path)
        >>> cfg.solver
        'mt1d'
    """

    # ── Solver ───────────────────────────────────────────────────────────────
    solver:      str   = "mt1d"
    freq_min:    float = 1e-4
    freq_max:    float = 1e4
    n_freqs:     int   = 30
    time_min:    float = 1e-6
    time_max:    float = 1e-2
    n_times:     int   = 25
    loop_radius: float = 50.0

    # ── Earth model ──────────────────────────────────────────────────────────
    n_layers_min: int            = 3
    n_layers_max: int            = 7
    rho_min:      float          = 1.0
    rho_max:      float          = 10_000.0
    depth_max:    float          = 2_000.0
    geology:      Optional[str]  = None

    # ── Dataset ──────────────────────────────────────────────────────────────
    n_samples:     int   = 10_000
    noise_level:   float = 0.05
    noise_type:    str   = "gaussian"
    include_phase: bool  = True
    seed:          Optional[int] = None
    n_jobs:        int   = 1

    # ── Output ───────────────────────────────────────────────────────────────
    output_dir:  Optional[str] = "."
    output_name: str           = "forward_dataset"
    verbose:     bool          = True

    # ─────────────────────────────────────────────────────────────────────────
    # Validation
    # ─────────────────────────────────────────────────────────────────────────

    def validate(self) -> None:
        """Check parameter ranges and raise :class:`ValueError` on errors.

        Call before a long dataset generation run to catch obvious
        mistakes (e.g. ``freq_min >= freq_max``) before any compute is
        spent.

        Raises
        ------
        ValueError
            Descriptive message pointing to the offending parameter.
        """
        _VALID_SOLVERS = {"mt1d", "csamt1d", "tem1d"}
        if self.solver not in _VALID_SOLVERS:
            raise ValueError(
                f"solver must be one of {_VALID_SOLVERS!r}, got {self.solver!r}."
            )

        if self.solver in ("mt1d", "csamt1d"):
            if self.freq_min <= 0.0:
                raise ValueError("freq_min must be strictly positive.")
            if self.freq_max <= self.freq_min:
                raise ValueError("freq_max must be strictly greater than freq_min.")
            if self.n_freqs < 4:
                raise ValueError("n_freqs must be at least 4.")

        if self.solver == "tem1d":
            if self.time_min <= 0.0:
                raise ValueError("time_min must be strictly positive.")
            if self.time_max <= self.time_min:
                raise ValueError("time_max must be strictly greater than time_min.")
            if self.n_times < 4:
                raise ValueError("n_times must be at least 4.")
            if self.loop_radius <= 0.0:
                raise ValueError("loop_radius must be strictly positive.")

        if self.n_layers_min < 2:
            raise ValueError("n_layers_min must be at least 2 (min: 1 layer + halfspace).")
        if self.n_layers_max < self.n_layers_min:
            raise ValueError("n_layers_max must be >= n_layers_min.")

        if self.geology is None:
            if self.rho_min <= 0.0:
                raise ValueError("rho_min must be strictly positive.")
            if self.rho_max <= self.rho_min:
                raise ValueError("rho_max must be strictly greater than rho_min.")
            if self.depth_max <= 0.0:
                raise ValueError("depth_max must be strictly positive.")

        if self.n_samples < 1:
            raise ValueError("n_samples must be at least 1.")
        if not 0.0 <= self.noise_level <= 1.0:
            raise ValueError("noise_level must be in [0, 1].")

        _VALID_NOISE = {"gaussian", "multiplicative", "field"}
        if self.noise_type not in _VALID_NOISE:
            raise ValueError(
                f"noise_type must be one of {_VALID_NOISE!r}, "
                f"got {self.noise_type!r}."
            )

    # ─────────────────────────────────────────────────────────────────────────
    # Assemblers
    # ─────────────────────────────────────────────────────────────────────────

    def freq_grid(self) -> np.ndarray:
        """Return the logarithmically spaced frequency array [Hz].

        This is the array that would be passed as *freqs* to the MT/CSAMT
        solvers and to :func:`~pycsamt.forward.batch.generate_dataset`.

        Returns
        -------
        numpy.ndarray, shape (n_freqs,)
        """
        return np.logspace(
            np.log10(self.freq_min), np.log10(self.freq_max), self.n_freqs
        )

    def time_grid(self) -> np.ndarray:
        """Return the logarithmically spaced time-gate array [s].

        This is the array that would be passed as *times* to the TEM solver
        and to :func:`~pycsamt.forward.batch.generate_dataset`.

        Returns
        -------
        numpy.ndarray, shape (n_times,)
        """
        return np.logspace(
            np.log10(self.time_min), np.log10(self.time_max), self.n_times
        )

    def to_dataset_kwargs(self) -> Dict[str, Any]:
        """Assemble keyword arguments for :func:`~pycsamt.forward.batch.generate_dataset`.

        The returned dict is ready to be unpacked directly::

            from pycsamt.forward.batch import generate_dataset
            ds = generate_dataset(**cfg.to_dataset_kwargs())

        Returns
        -------
        dict
            All keyword arguments for ``generate_dataset``, with numpy
            arrays built from the stored scalar parameters.
        """
        kw: Dict[str, Any] = dict(
            solver=self.solver,
            n_samples=self.n_samples,
            n_layers=(self.n_layers_min
                      if self.n_layers_min == self.n_layers_max
                      else (self.n_layers_min, self.n_layers_max)),
            rho_range=(self.rho_min, self.rho_max),
            depth_max=self.depth_max,
            loop_radius=self.loop_radius,
            noise_level=self.noise_level,
            noise_type=self.noise_type,
            geology=self.geology,
            include_phase=self.include_phase,
            seed=self.seed,
            n_jobs=self.n_jobs,
            verbose=self.verbose,
        )

        if self.solver in ("mt1d", "csamt1d"):
            kw["freqs"] = self.freq_grid()
            kw["times"] = None
        else:
            kw["times"] = self.time_grid()
            kw["freqs"] = None

        if self.output_dir is not None:
            out_dir = Path(self.output_dir).expanduser()
            out_dir.mkdir(parents=True, exist_ok=True)
            kw["output"] = str(out_dir / f"{self.output_name}.npz")
        else:
            kw["output"] = None

        return kw

    # ─────────────────────────────────────────────────────────────────────────
    # Config file I/O
    # ─────────────────────────────────────────────────────────────────────────

    def to_template(
        self,
        path: Union[str, Path] = "forward_config.py",
        *,
        fmt: Optional[str] = None,
    ) -> Path:
        """Write this configuration to an annotated source-of-truth file.

        Parameters
        ----------
        path : path-like, default "forward_config.py"
            Destination file.  Suffixes ``.py``, ``.json``, ``.yml``,
            and ``.yaml`` select the output format automatically.
        fmt : {"py", "json", "yml", "yaml"}, optional
            Explicit format override.  When omitted, the suffix of *path*
            is used.  Paths without a suffix default to Python.

        Returns
        -------
        pathlib.Path
            Path of the generated file.
        """
        return write_config_template(
            path,
            self,
            _FORWARD_CONFIG_SCHEMA,
            fmt=fmt,
            title="PyCSAMT forward modelling configuration",
        )

    @classmethod
    def write_template(
        cls,
        path: Union[str, Path] = "forward_config.py",
        *,
        fmt: Optional[str] = None,
    ) -> Path:
        """Generate a documented source-of-truth configuration file.

        Creates a file with default parameter values and an inline comment
        for every parameter explaining its role.  Edit the generated file,
        then load it with :meth:`from_file`.

        Parameters
        ----------
        path : path-like, default "forward_config.py"
            Destination file.  The suffix determines the output format
            (``.py``, ``.json``, ``.yml``).
        fmt : {"py", "json", "yml", "yaml"}, optional
            Explicit format override.

        Returns
        -------
        pathlib.Path
            Path of the generated file.

        Examples
        --------
        >>> from pycsamt.forward.config import ForwardConfig
        >>> path = ForwardConfig.write_template("my_forward.yml")
        >>> path.suffix
        '.yml'
        """
        return cls().to_template(path, fmt=fmt)

    @classmethod
    def from_file(
        cls,
        path: Union[str, Path],
        *,
        strict: bool = True,
    ) -> "ForwardConfig":
        """Load a configuration from a source-of-truth file.

        Parameters
        ----------
        path : path-like
            Python, JSON, YML, or YAML file generated by
            :meth:`write_template` or following the same structure.
        strict : bool, default True
            If ``True``, unknown keys raise :class:`ValueError`.
            If ``False``, unknown keys are silently ignored.

        Returns
        -------
        ForwardConfig
            Configuration populated from the file values.

        Examples
        --------
        >>> ForwardConfig.write_template("forward_config.json")
        PosixPath('forward_config.json')
        >>> cfg = ForwardConfig.from_file("forward_config.json")
        >>> cfg.solver
        'mt1d'
        """
        values = read_config_file(path, cls, strict=strict)
        return cls(**values)

    #: Alias — matches the convention used by ModEmConfig and OccamConfig.
    read = from_file

    # ─────────────────────────────────────────────────────────────────────────
    # repr / summary
    # ─────────────────────────────────────────────────────────────────────────

    def summary(self) -> str:
        """Return a human-readable multi-line summary of the configuration."""
        lines = [
            "ForwardConfig",
            f"  {'solver':<20s} = {self.solver!r}",
        ]
        if self.solver in ("mt1d", "csamt1d"):
            lines += [
                f"  {'freq_min':<20s} = {self.freq_min:.3g} Hz",
                f"  {'freq_max':<20s} = {self.freq_max:.3g} Hz",
                f"  {'n_freqs':<20s} = {self.n_freqs}",
            ]
        else:
            lines += [
                f"  {'time_min':<20s} = {self.time_min:.3g} s",
                f"  {'time_max':<20s} = {self.time_max:.3g} s",
                f"  {'n_times':<20s} = {self.n_times}",
                f"  {'loop_radius':<20s} = {self.loop_radius} m",
            ]
        layer_s = (
            str(self.n_layers_min)
            if self.n_layers_min == self.n_layers_max
            else f"{self.n_layers_min}–{self.n_layers_max}"
        )
        lines += [
            f"  {'n_layers':<20s} = {layer_s}",
        ]
        if self.geology:
            lines.append(f"  {'geology':<20s} = {self.geology!r}")
        else:
            lines += [
                f"  {'rho_min':<20s} = {self.rho_min:.3g} Ω·m",
                f"  {'rho_max':<20s} = {self.rho_max:.3g} Ω·m",
                f"  {'depth_max':<20s} = {self.depth_max:.3g} m",
            ]
        lines += [
            f"  {'n_samples':<20s} = {self.n_samples:,}",
            f"  {'noise_level':<20s} = {self.noise_level}  ({self.noise_type})",
            f"  {'seed':<20s} = {self.seed!r}",
            f"  {'n_jobs':<20s} = {self.n_jobs}",
        ]
        if self.output_dir:
            lines.append(
                f"  {'output':<20s} = {self.output_dir}/{self.output_name}.npz"
            )
        return "\n".join(lines)

    def __repr__(self) -> str:
        return self.summary()
