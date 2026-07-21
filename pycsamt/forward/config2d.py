# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Configuration for pyCSAMT 2-D MT forward modelling.

:class:`ForwardConfig2D` is the configuration object for
:class:`~pycsamt.forward.em2d.MT2DForward`.  It records every parameter
needed to build a :class:`~pycsamt.forward.grid2d.Grid2D` and run the
finite-difference solver, and can be exported to Python, JSON, or YAML.

Quick start
-----------
::

    from pycsamt.forward.config2d import ForwardConfig2D

    # Write an annotated source-of-truth file
    ForwardConfig2D.write_template("my_2d.yml")

    # Edit my_2d.yml …

    # Load, validate, build grid, run solver
    cfg = ForwardConfig2D.from_file("my_2d.yml")
    cfg.validate()

    grid = cfg.to_grid()
    from pycsamt.forward.em2d import MT2DForward
    resp = MT2DForward(cfg.freq_grid(), grid).run()
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from ..models.config_io import (
    ConfigParameter,
    read_config_file,
    write_config_template,
)

__all__ = ["ForwardConfig2D"]

_FWD2D_SCHEMA: list[ConfigParameter] = [
    # ── Solver ───────────────────────────────────────────────────────────────
    ConfigParameter(
        "freq_min",
        "Minimum frequency in Hz. The frequency grid is "
        "np.logspace(log10(freq_min), log10(freq_max), n_freqs). "
        "Typical AMT lower bound: 1e-4 Hz.",
        "Solver",
    ),
    ConfigParameter(
        "freq_max",
        "Maximum frequency in Hz. Must be strictly greater than freq_min. "
        "Typical AMT upper bound: 1e4 Hz.",
        "Solver",
    ),
    ConfigParameter(
        "n_freqs",
        "Number of logarithmically spaced frequency points. "
        "15-30 points covers the AMT band adequately.",
        "Solver",
    ),
    # ── Grid ─────────────────────────────────────────────────────────────────
    ConfigParameter(
        "nx",
        "Number of core cells in the horizontal (x) direction, "
        "excluding padding. Minimum 20; 40-80 is typical.",
        "Grid",
    ),
    ConfigParameter(
        "nz",
        "Number of core cells in the vertical (z, depth) direction, "
        "excluding padding. Minimum 15; 30-60 is typical.",
        "Grid",
    ),
    ConfigParameter(
        "x_max",
        "Total horizontal extent of the core model region in metres. "
        "Should span all stations plus at least one skin-depth on each side.",
        "Grid",
    ),
    ConfigParameter(
        "z_max",
        "Total depth of the core model region in metres. "
        "Should be at least 3× the maximum investigation depth "
        "(≈ skin depth at the lowest frequency).",
        "Grid",
    ),
    ConfigParameter(
        "n_pad",
        "Number of exponentially growing padding cells added to both "
        "horizontal sides and the bottom.  8-12 cells is standard.",
        "Grid",
    ),
    ConfigParameter(
        "pad_factor",
        "Growth factor for padding cell sizes. Each successive padding "
        "cell is pad_factor times wider/deeper than the previous. "
        "Values 1.2-1.5 are appropriate; 1.3 is a safe default.",
        "Grid",
    ),
    # ── Earth model ──────────────────────────────────────────────────────────
    ConfigParameter(
        "bg_rho",
        "Background (starting) resistivity in Ω·m for the halfspace or "
        "anomaly model constructors. Ignored when model_type is 'random'.",
        "Earth Model",
    ),
    ConfigParameter(
        "model_type",
        "Model constructor to use. Options: 'halfspace' — uniform "
        "background; 'anomaly' — one rectangular conductive/resistive "
        "block on the background; 'random' — random layered model with "
        "optional lateral variation.",
        "Earth Model",
    ),
    ConfigParameter(
        "anomaly_rho",
        "Anomaly resistivity in Ω·m. Used only when model_type='anomaly'.",
        "Earth Model",
    ),
    ConfigParameter(
        "anomaly_x_lo",
        "Left edge of the anomaly in core x-coordinates [m]. "
        "Used only when model_type='anomaly'.",
        "Earth Model",
    ),
    ConfigParameter(
        "anomaly_x_hi",
        "Right edge of the anomaly [m].",
        "Earth Model",
    ),
    ConfigParameter(
        "anomaly_z_lo",
        "Top of the anomaly in depth [m].",
        "Earth Model",
    ),
    ConfigParameter(
        "anomaly_z_hi",
        "Bottom of the anomaly in depth [m].",
        "Earth Model",
    ),
    # ── Stations ─────────────────────────────────────────────────────────────
    ConfigParameter(
        "n_stations",
        "Number of MT stations distributed uniformly across the core region.",
        "Stations",
    ),
    ConfigParameter(
        "station_x_max",
        "x-coordinate of the rightmost station [m]. "
        "Defaults to x_max when null/None.",
        "Stations",
    ),
    # ── Output ───────────────────────────────────────────────────────────────
    ConfigParameter(
        "verbose",
        "Print per-frequency progress during the forward solve.",
        "Output",
    ),
]


@dataclass
class ForwardConfig2D:
    """Collect settings for a 2-D MT finite-difference forward run.

    The configuration covers four areas: frequency grid, FD grid geometry,
    starting earth model, and station layout.

    Parameters
    ----------
    freq_min, freq_max : float
        Frequency range [Hz].
    n_freqs : int
        Number of logarithmically spaced frequencies.
    nx, nz : int
        Core cell counts (x, z).
    x_max, z_max : float
        Core region extents [m].
    n_pad : int
        Padding cells per side.
    pad_factor : float
        Padding growth factor.
    bg_rho : float
        Background resistivity [Ω·m].
    model_type : {'halfspace', 'anomaly', 'random'}
        Grid constructor to use.
    anomaly_rho : float
        Anomaly resistivity [Ω·m] (model_type='anomaly' only).
    anomaly_x_lo, anomaly_x_hi : float
        Horizontal anomaly bounds [m].
    anomaly_z_lo, anomaly_z_hi : float
        Depth anomaly bounds [m].
    n_stations : int
        Number of equally spaced surface stations.
    station_x_max : float or None
        Right-most station position [m].
    verbose : bool

    Examples
    --------
    Default halfspace::

        >>> cfg = ForwardConfig2D()
        >>> cfg.model_type
        'halfspace'
        >>> cfg.validate()

    Template round-trip::

        >>> import tempfile, os
        >>> with tempfile.NamedTemporaryFile(suffix='.yml', delete=False) as f:
        ...     path = f.name
        >>> _ = ForwardConfig2D.write_template(path)
        >>> cfg2 = ForwardConfig2D.from_file(path)
        >>> os.unlink(path)
    """

    # ── Solver ───────────────────────────────────────────────────────────────
    freq_min: float = 1e-3
    freq_max: float = 1e3
    n_freqs: int = 20

    # ── Grid ─────────────────────────────────────────────────────────────────
    nx: int = 40
    nz: int = 30
    x_max: float = 10_000.0
    z_max: float = 6_000.0
    n_pad: int = 8
    pad_factor: float = 1.3

    # ── Earth model ──────────────────────────────────────────────────────────
    bg_rho: float = 100.0
    model_type: str = "halfspace"
    anomaly_rho: float = 5.0
    anomaly_x_lo: float = 2_000.0
    anomaly_x_hi: float = 6_000.0
    anomaly_z_lo: float = 300.0
    anomaly_z_hi: float = 1_500.0

    # ── Stations ─────────────────────────────────────────────────────────────
    n_stations: int = 10
    station_x_max: float | None = None

    # ── Output ───────────────────────────────────────────────────────────────
    verbose: bool = True

    # ─────────────────────────────────────────────────────────────────────────
    # Validation
    # ─────────────────────────────────────────────────────────────────────────

    def validate(self) -> None:
        """Raise :class:`ValueError` on invalid parameter values."""
        if self.freq_min <= 0 or self.freq_max <= self.freq_min:
            raise ValueError("freq_min must be > 0 and freq_max > freq_min.")
        if self.n_freqs < 2:
            raise ValueError("n_freqs must be at least 2.")
        if self.nx < 4 or self.nz < 4:
            raise ValueError("nx and nz must each be at least 4.")
        if self.x_max <= 0 or self.z_max <= 0:
            raise ValueError("x_max and z_max must be strictly positive.")
        if self.n_pad < 0:
            raise ValueError("n_pad must be non-negative.")
        if self.pad_factor <= 1.0:
            raise ValueError("pad_factor must be greater than 1.0.")
        if self.bg_rho <= 0:
            raise ValueError("bg_rho must be strictly positive.")
        _VALID = {"halfspace", "anomaly", "random"}
        if self.model_type not in _VALID:
            raise ValueError(f"model_type must be one of {_VALID!r}.")
        if self.model_type == "anomaly":
            if self.anomaly_rho <= 0:
                raise ValueError("anomaly_rho must be strictly positive.")
            if self.anomaly_x_lo >= self.anomaly_x_hi:
                raise ValueError("anomaly_x_lo must be < anomaly_x_hi.")
            if self.anomaly_z_lo >= self.anomaly_z_hi:
                raise ValueError("anomaly_z_lo must be < anomaly_z_hi.")
        if self.n_stations < 1:
            raise ValueError("n_stations must be at least 1.")

    # ─────────────────────────────────────────────────────────────────────────
    # Assemblers
    # ─────────────────────────────────────────────────────────────────────────

    def freq_grid(self) -> np.ndarray:
        """Return the logarithmically spaced frequency array [Hz]."""
        return np.logspace(
            np.log10(self.freq_min), np.log10(self.freq_max), self.n_freqs
        )

    def to_grid(self, seed=None) -> Grid2D:
        """Construct a :class:`~pycsamt.forward.grid2d.Grid2D` from this config.

        Parameters
        ----------
        seed : int or None
            Random seed (used only for ``model_type='random'``).

        Returns
        -------
        Grid2D
        """
        from .grid2d import Grid2D

        kw = dict(
            nx=self.nx,
            nz=self.nz,
            x_max=self.x_max,
            z_max=self.z_max,
            n_pad=self.n_pad,
            pad_factor=self.pad_factor,
            n_stations=self.n_stations,
            station_x_max=self.station_x_max,
        )
        if self.model_type == "halfspace":
            return Grid2D.halfspace(rho=self.bg_rho, **kw)
        if self.model_type == "anomaly":
            return Grid2D.with_anomaly(
                bg_rho=self.bg_rho,
                anomaly_rho=self.anomaly_rho,
                anomaly_bounds=(
                    self.anomaly_x_lo,
                    self.anomaly_x_hi,
                    self.anomaly_z_lo,
                    self.anomaly_z_hi,
                ),
                **kw,
            )
        # random: Grid2D.random() has no station_x_max parameter.
        kw.pop("station_x_max")
        return Grid2D.random(seed=seed, **kw)

    def to_solver_kwargs(self) -> dict[str, Any]:
        """Return kwargs for :class:`~pycsamt.forward.em2d.MT2DForward`."""
        return dict(freqs=self.freq_grid(), verbose=self.verbose)

    # ─────────────────────────────────────────────────────────────────────────
    # Config file I/O
    # ─────────────────────────────────────────────────────────────────────────

    def to_template(
        self,
        path: str | Path = "forward_config_2d.py",
        *,
        fmt: str | None = None,
    ) -> Path:
        """Write to an annotated source-of-truth file."""
        return write_config_template(
            path,
            self,
            _FWD2D_SCHEMA,
            fmt=fmt,
            title="PyCSAMT 2-D MT forward configuration",
        )

    @classmethod
    def write_template(
        cls,
        path: str | Path = "forward_config_2d.py",
        *,
        fmt: str | None = None,
    ) -> Path:
        """Generate a documented source-of-truth file with default values.

        Parameters
        ----------
        path : path-like
        fmt : {"py", "json", "yml", "yaml"}, optional

        Returns
        -------
        pathlib.Path
        """
        return cls().to_template(path, fmt=fmt)

    @classmethod
    def from_file(
        cls,
        path: str | Path,
        *,
        strict: bool = True,
    ) -> ForwardConfig2D:
        """Load from a generated source-of-truth file.

        Parameters
        ----------
        path : path-like
        strict : bool

        Returns
        -------
        ForwardConfig2D
        """
        values = read_config_file(path, cls, strict=strict)
        return cls(**values)

    read = from_file

    # ─────────────────────────────────────────────────────────────────────────
    # Summary
    # ─────────────────────────────────────────────────────────────────────────

    def summary(self) -> str:
        lines = [
            "ForwardConfig2D",
            f"  {'model_type':<22s} = {self.model_type!r}",
            f"  {'bg_rho':<22s} = {self.bg_rho} Ω·m",
            f"  {'freq range':<22s} = {self.freq_min:.3g}–{self.freq_max:.3g} Hz  ({self.n_freqs} pts)",
            f"  {'grid (nx × nz)':<22s} = {self.nx} × {self.nz}  (core, +{self.n_pad} pad)",
            f"  {'x_max':<22s} = {self.x_max:.0f} m",
            f"  {'z_max':<22s} = {self.z_max:.0f} m",
            f"  {'n_stations':<22s} = {self.n_stations}",
        ]
        if self.model_type == "anomaly":
            lines.append(
                f"  {'anomaly':<22s} = {self.anomaly_rho} Ω·m  "
                f"x=[{self.anomaly_x_lo},{self.anomaly_x_hi}]  "
                f"z=[{self.anomaly_z_lo},{self.anomaly_z_hi}]"
            )
        return "\n".join(lines)

    def __repr__(self) -> str:
        return self.summary()
