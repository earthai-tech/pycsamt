# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Configuration for pyCSAMT quasi-3D MT forward modelling.

:class:`ForwardConfig3D` follows the same source-of-truth file pattern as
:class:`~pycsamt.forward.config.ForwardConfig` and
:class:`~pycsamt.forward.config2d.ForwardConfig2D`.

Quick start
-----------
::

    from pycsamt.forward.config3d import ForwardConfig3D

    ForwardConfig3D.write_template("my_3d.yml")   # annotated defaults

    # edit my_3d.yml …

    cfg = ForwardConfig3D.from_file("my_3d.yml")
    cfg.validate()

    grid = cfg.to_grid()
    from pycsamt.forward.em3d import MT3DForward
    resp = MT3DForward(cfg.freq_grid(), grid, verbose=cfg.verbose).run()
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

__all__ = ["ForwardConfig3D"]


_FWD3D_SCHEMA: list[ConfigParameter] = [
    # ── Solver ───────────────────────────────────────────────────────────────
    ConfigParameter(
        "freq_min",
        "Minimum frequency [Hz] for the logarithmic frequency grid.  "
        "Grid is np.logspace(log10(freq_min), log10(freq_max), n_freqs).",
        "Solver",
    ),
    ConfigParameter(
        "freq_max",
        "Maximum frequency [Hz].  Must exceed freq_min.",
        "Solver",
    ),
    ConfigParameter(
        "n_freqs",
        "Number of frequency points.  15–25 covers the AMT band well.",
        "Solver",
    ),
    ConfigParameter(
        "method",
        "Forward solver method.  Currently only 'quasi3d' is available: "
        "runs MT2DForward on orthogonal XZ and YZ profile slices and "
        "averages the impedance tensor contributions.",
        "Solver",
    ),
    # ── Grid ─────────────────────────────────────────────────────────────────
    ConfigParameter(
        "nx",
        "Number of core cells in x (easting), excluding padding.",
        "Grid",
    ),
    ConfigParameter(
        "ny",
        "Number of core cells in y (northing), excluding padding.",
        "Grid",
    ),
    ConfigParameter(
        "nz",
        "Number of core cells in z (depth), excluding padding.",
        "Grid",
    ),
    ConfigParameter(
        "x_max",
        "Core horizontal extent in x [m].  Should comfortably span all "
        "stations plus one skin-depth on each side.",
        "Grid",
    ),
    ConfigParameter(
        "y_max",
        "Core horizontal extent in y [m].",
        "Grid",
    ),
    ConfigParameter(
        "z_max",
        "Total depth of the core model [m].  Should be ≥ 3× the maximum "
        "investigation depth (skin depth at the lowest frequency).",
        "Grid",
    ),
    ConfigParameter(
        "n_pad",
        "Padding cells on each side in x and y and at the bottom in z.  "
        "8–12 cells is standard for MT FD modelling.",
        "Grid",
    ),
    ConfigParameter(
        "pad_factor",
        "Exponential growth factor for padding cell sizes.  "
        "Values in [1.2, 1.5]; 1.3 is a safe default.",
        "Grid",
    ),
    # ── Earth model ──────────────────────────────────────────────────────────
    ConfigParameter(
        "bg_rho",
        "Background resistivity [Ω·m] for 'halfspace' and 'block_anomaly' "
        "model types.",
        "Earth Model",
    ),
    ConfigParameter(
        "model_type",
        "Model constructor: 'halfspace' — uniform background; "
        "'block_anomaly' — one 3-D rectangular anomaly on the background; "
        "'random_layered' — random horizontal layers with optional lateral "
        "variation driven by Gaussian random fields.",
        "Earth Model",
    ),
    ConfigParameter(
        "anomaly_rho",
        "Anomaly resistivity [Ω·m].  Used only when model_type='block_anomaly'.",
        "Earth Model",
    ),
    ConfigParameter(
        "anomaly_x_lo",
        "Left edge of the 3-D anomaly in core x-coordinates [m].",
        "Earth Model",
    ),
    ConfigParameter(
        "anomaly_x_hi",
        "Right edge of the anomaly [m].",
        "Earth Model",
    ),
    ConfigParameter(
        "anomaly_y_lo",
        "Near edge of the anomaly in core y-coordinates [m].",
        "Earth Model",
    ),
    ConfigParameter(
        "anomaly_y_hi",
        "Far edge of the anomaly [m].",
        "Earth Model",
    ),
    ConfigParameter(
        "anomaly_z_lo",
        "Top of the anomaly in depth [m].",
        "Earth Model",
    ),
    ConfigParameter(
        "anomaly_z_hi",
        "Bottom of the anomaly [m].",
        "Earth Model",
    ),
    ConfigParameter(
        "n_layers",
        "Number of horizontal layers for model_type='random_layered'.",
        "Earth Model",
    ),
    ConfigParameter(
        "lateral_variation",
        "Enable Gaussian-random-field lateral resistivity perturbation for "
        "model_type='random_layered'.  Produces realistic 3-D structure.",
        "Earth Model",
    ),
    ConfigParameter(
        "corr_length",
        "GRF spatial correlation length [m] for model_type='random_layered'.  "
        "Larger values give smoother lateral variation.",
        "Earth Model",
    ),
    # ── Stations ─────────────────────────────────────────────────────────────
    ConfigParameter(
        "nx_stations",
        "Number of stations in the x direction (regular grid).",
        "Stations",
    ),
    ConfigParameter(
        "ny_stations",
        "Number of stations in the y direction (regular grid).  "
        "Total stations = nx_stations × ny_stations.",
        "Stations",
    ),
    # ── Output ───────────────────────────────────────────────────────────────
    ConfigParameter(
        "verbose",
        "Print per-profile progress during the forward solve.",
        "Output",
    ),
]


@dataclass
class ForwardConfig3D:
    """Collect settings for a quasi-3D MT forward modelling run.

    Parameters
    ----------
    freq_min, freq_max : float
        Frequency range [Hz].
    n_freqs : int
    method : str
        Forward method ('quasi3d').
    nx, ny, nz : int
        Core cell counts.
    x_max, y_max, z_max : float
        Core extents [m].
    n_pad : int
    pad_factor : float
    bg_rho : float
    model_type : {'halfspace', 'block_anomaly', 'random_layered'}
    anomaly_rho : float
    anomaly_x_lo, anomaly_x_hi : float
    anomaly_y_lo, anomaly_y_hi : float
    anomaly_z_lo, anomaly_z_hi : float
    n_layers : int
    lateral_variation : bool
    corr_length : float
    nx_stations, ny_stations : int
    verbose : bool

    Examples
    --------
    Default halfspace::

        >>> cfg = ForwardConfig3D()
        >>> cfg.validate()
        >>> grid = cfg.to_grid()

    Block anomaly::

        >>> cfg = ForwardConfig3D(
        ...     model_type='block_anomaly',
        ...     anomaly_rho=5.0,
        ...     anomaly_x_lo=2000., anomaly_x_hi=6000.,
        ...     anomaly_y_lo=2000., anomaly_y_hi=6000.,
        ...     anomaly_z_lo=300.,  anomaly_z_hi=1500.,
        ... )
    """

    # ── Solver ───────────────────────────────────────────────────────────────
    freq_min: float = 1e-2
    freq_max: float = 1e3
    n_freqs:  int   = 15
    method:   str   = "quasi3d"

    # ── Grid ─────────────────────────────────────────────────────────────────
    nx:         int   = 20
    ny:         int   = 20
    nz:         int   = 15
    x_max:      float = 8_000.0
    y_max:      float = 8_000.0
    z_max:      float = 4_000.0
    n_pad:      int   = 8
    pad_factor: float = 1.3

    # ── Earth model ──────────────────────────────────────────────────────────
    bg_rho:         float          = 100.0
    model_type:     str            = "halfspace"
    anomaly_rho:    float          = 5.0
    anomaly_x_lo:   float          = 2_000.0
    anomaly_x_hi:   float          = 6_000.0
    anomaly_y_lo:   float          = 2_000.0
    anomaly_y_hi:   float          = 6_000.0
    anomaly_z_lo:   float          =   300.0
    anomaly_z_hi:   float          = 1_500.0
    n_layers:       int            = 4
    lateral_variation: bool        = True
    corr_length:    float          = 2_000.0

    # ── Stations ─────────────────────────────────────────────────────────────
    nx_stations: int = 5
    ny_stations: int = 5

    # ── Output ───────────────────────────────────────────────────────────────
    verbose: bool = True

    # ─────────────────────────────────────────────────────────────────────────
    # Validation
    # ─────────────────────────────────────────────────────────────────────────

    def validate(self) -> None:
        """Raise :class:`ValueError` on invalid parameters."""
        if self.freq_min <= 0 or self.freq_max <= self.freq_min:
            raise ValueError("freq_min must be > 0 and freq_max > freq_min.")
        if self.n_freqs < 2:
            raise ValueError("n_freqs must be at least 2.")
        if self.method not in {"quasi3d"}:
            raise ValueError(f"method must be 'quasi3d', got {self.method!r}.")
        if any(n < 4 for n in (self.nx, self.ny, self.nz)):
            raise ValueError("nx, ny, nz must each be at least 4.")
        if any(v <= 0 for v in (self.x_max, self.y_max, self.z_max)):
            raise ValueError("x_max, y_max, z_max must be strictly positive.")
        if self.n_pad < 0:
            raise ValueError("n_pad must be non-negative.")
        if self.pad_factor <= 1.0:
            raise ValueError("pad_factor must be > 1.0.")
        if self.bg_rho <= 0:
            raise ValueError("bg_rho must be strictly positive.")
        valid = {"halfspace", "block_anomaly", "random_layered"}
        if self.model_type not in valid:
            raise ValueError(f"model_type must be one of {valid!r}.")
        if self.model_type == "block_anomaly":
            if self.anomaly_rho <= 0:
                raise ValueError("anomaly_rho must be strictly positive.")
            for lo, hi, ax in (
                (self.anomaly_x_lo, self.anomaly_x_hi, "x"),
                (self.anomaly_y_lo, self.anomaly_y_hi, "y"),
                (self.anomaly_z_lo, self.anomaly_z_hi, "z"),
            ):
                if lo >= hi:
                    raise ValueError(f"anomaly_{ax}_lo must be < anomaly_{ax}_hi.")
        if self.nx_stations < 1 or self.ny_stations < 1:
            raise ValueError("nx_stations and ny_stations must be ≥ 1.")

    # ─────────────────────────────────────────────────────────────────────────
    # Assemblers
    # ─────────────────────────────────────────────────────────────────────────

    def freq_grid(self) -> np.ndarray:
        """Return the logarithmically spaced frequency array [Hz]."""
        return np.logspace(
            np.log10(self.freq_min), np.log10(self.freq_max), self.n_freqs
        )

    def to_grid(self, seed=None) -> Grid3D:
        """Construct a :class:`~pycsamt.forward.grid3d.Grid3D` from this config.

        Parameters
        ----------
        seed : int or None
            Random seed (used only for ``model_type='random_layered'``).

        Returns
        -------
        Grid3D
        """
        from .grid3d import Grid3D

        kw = dict(
            nx=self.nx, ny=self.ny, nz=self.nz,
            x_max=self.x_max, y_max=self.y_max, z_max=self.z_max,
            n_pad=self.n_pad, pad_factor=self.pad_factor,
            nx_stations=self.nx_stations,
            ny_stations=self.ny_stations,
        )
        if self.model_type == "halfspace":
            return Grid3D.halfspace(rho=self.bg_rho, **kw)

        if self.model_type == "block_anomaly":
            return Grid3D.block_anomaly(
                bg_rho=self.bg_rho,
                anomaly_rho=self.anomaly_rho,
                bounds=(
                    self.anomaly_x_lo, self.anomaly_x_hi,
                    self.anomaly_y_lo, self.anomaly_y_hi,
                    self.anomaly_z_lo, self.anomaly_z_hi,
                ),
                **kw,
            )

        # random_layered
        return Grid3D.random_layered(
            n_layers=self.n_layers,
            lateral_variation=self.lateral_variation,
            corr_length=self.corr_length,
            seed=seed,
            **kw,
        )

    def to_solver_kwargs(self) -> dict[str, Any]:
        """Return kwargs for :class:`~pycsamt.forward.em3d.MT3DForward`."""
        return dict(
            freqs=self.freq_grid(),
            method=self.method,
            verbose=self.verbose,
        )

    # ─────────────────────────────────────────────────────────────────────────
    # Config file I/O
    # ─────────────────────────────────────────────────────────────────────────

    def to_template(
        self,
        path: str | Path = "forward_config_3d.py",
        *,
        fmt: str | None = None,
    ) -> Path:
        """Write an annotated source-of-truth file."""
        return write_config_template(
            path, self, _FWD3D_SCHEMA, fmt=fmt,
            title="PyCSAMT quasi-3D MT forward configuration",
        )

    @classmethod
    def write_template(
        cls,
        path: str | Path = "forward_config_3d.py",
        *,
        fmt: str | None = None,
    ) -> Path:
        """Generate a documented source-of-truth configuration file.

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
    ) -> ForwardConfig3D:
        """Load from a source-of-truth file.

        Parameters
        ----------
        path : path-like
        strict : bool

        Returns
        -------
        ForwardConfig3D
        """
        return cls(**read_config_file(path, cls, strict=strict))

    read = from_file

    # ─────────────────────────────────────────────────────────────────────────

    def summary(self) -> str:
        n_tot = self.nx_stations * self.ny_stations
        lines = [
            "ForwardConfig3D",
            f"  {'method':<24s} = {self.method!r}",
            f"  {'model_type':<24s} = {self.model_type!r}  (bg={self.bg_rho} Ω·m)",
            f"  {'freq range':<24s} = {self.freq_min:.3g}–{self.freq_max:.3g} Hz"
            f"  ({self.n_freqs} pts)",
            f"  {'grid (nx×ny×nz)':<24s} = {self.nx}×{self.ny}×{self.nz}"
            f"  (core, +{self.n_pad} pad)",
            f"  {'extents (x,y,z)':<24s} = {self.x_max:.0f} m × {self.y_max:.0f} m"
            f" × {self.z_max:.0f} m",
            f"  {'stations':<24s} = {self.nx_stations}×{self.ny_stations}"
            f" = {n_tot} total",
        ]
        if self.model_type == "block_anomaly":
            lines.append(
                f"  {'anomaly':<24s} = {self.anomaly_rho} Ω·m  "
                f"x=[{self.anomaly_x_lo},{self.anomaly_x_hi}]  "
                f"y=[{self.anomaly_y_lo},{self.anomaly_y_hi}]  "
                f"z=[{self.anomaly_z_lo},{self.anomaly_z_hi}]"
            )
        return "\n".join(lines)

    def __repr__(self) -> str:
        return self.summary()
