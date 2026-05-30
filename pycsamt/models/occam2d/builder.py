# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""InputBuilder — produce all Occam2D input files from EDI data.

``InputBuilder`` is the main v2 entry point for Occam2D.  It consumes an
``EDICollection`` or ``Sites`` object (the single source of truth) and
writes the four files the Fortran binary needs:

    OccamDataFile.dat   (data)
    Occam2DMesh         (mesh)
    Occam2DModel        (model parameterisation)
    Startup             (inversion control)

Usage
-----
>>> from pycsamt.models.occam2d import InputBuilder
>>> from pycsamt.site import Sites
>>>
>>> sites   = Sites.from_dir("edi/")
>>> builder = InputBuilder(sites, workdir="occam_run/")
>>> builder.build(modes=["TE", "TM"], n_layers=30)
>>> # → writes four files to occam_run/
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Union

from .base     import OccamBase
from .config   import OccamConfig
from .data     import OccamData
from .mesh     import OccamMesh
from .model    import OccamModel
from .startup  import OccamStartup

PathLike = Union[str, Path]

__all__ = ["InputBuilder"]


class InputBuilder(OccamBase):
    """Build all Occam2D input files from an EDI source.

    Parameters
    ----------
    source : EDICollection | Sites
        Loaded EDI data.
    workdir : path-like
        Directory where output files are written (created if absent).
    config : OccamConfig, optional
        Full run configuration.  Keyword arguments passed to ``build()``
        override individual config fields.

    Attributes
    ----------
    data : OccamData | None
        Populated after calling ``build()``.
    mesh : OccamMesh | None
    model : OccamModel | None
    startup : OccamStartup | None
    """

    def __init__(
        self,
        source,
        workdir: PathLike = ".",
        config: Optional[OccamConfig] = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.source  = source
        self.workdir = Path(workdir)
        self.config  = config or OccamConfig()

        self.data:    Optional[OccamData]    = None
        self.mesh:    Optional[OccamMesh]    = None
        self.model:   Optional[OccamModel]   = None
        self.startup: Optional[OccamStartup] = None

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------
    def build(
        self,
        modes: Optional[List[str]] = None,
        n_layers: Optional[int]    = None,
        cell_size: Optional[float] = None,
        error_floor_rho: Optional[float]   = None,
        error_floor_phase: Optional[float] = None,
        freq_min: Optional[float] = None,
        freq_max: Optional[float] = None,
        title: str = "pycsamt Occam2D run",
        **kwargs,
    ) -> "InputBuilder":
        """Build and write all four Occam2D input files.

        Parameters
        ----------
        modes : list[str], optional
            Override ``config.modes``.
        n_layers : int, optional
            Override ``config.n_layers``.
        cell_size : float, optional
            Override ``config.cell_size_horizontal`` (metres).
        error_floor_rho, error_floor_phase : float, optional
            Override error-floor settings in config.
        freq_min, freq_max : float, optional
            Restrict to a frequency sub-band (Hz).
        title : str
            Title written into the data file header.

        Returns
        -------
        self
            Enables method chaining.
        """
        # Apply one-shot overrides
        cfg = self.config
        if modes            is not None: cfg.modes                  = modes
        if n_layers         is not None: cfg.n_layers               = n_layers
        if cell_size        is not None: cfg.cell_size_horizontal   = cell_size
        if error_floor_rho  is not None: cfg.error_floor_rho        = error_floor_rho
        if error_floor_phase is not None: cfg.error_floor_phase     = error_floor_phase
        if freq_min         is not None: cfg.freq_min               = freq_min
        if freq_max         is not None: cfg.freq_max               = freq_max

        self.workdir.mkdir(parents=True, exist_ok=True)

        # TODO: implement each step and remove NotImplementedError
        # Step 1 — data file
        self.data = OccamData.from_edi(self.source, config=cfg, title=title)
        self.data.write(self.workdir / cfg.data_file)

        # Step 2 — mesh
        self.mesh = OccamMesh.from_data(self.data, config=cfg)
        self.mesh.write(self.workdir / cfg.mesh_file)

        # Step 3 — model
        self.model = OccamModel.from_mesh(self.mesh, config=cfg)
        self.model.write(self.workdir / cfg.model_file)

        # Step 4 — startup
        self.startup = OccamStartup.from_model(self.model, config=cfg)
        self.startup.write(self.workdir / cfg.startup_file)

        if self.verbose:
            self.logger.info(
                "InputBuilder: wrote 4 files to %s", self.workdir
            )

        return self

    # ------------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------------
    @property
    def is_ready(self) -> bool:
        """True if all four files have been built."""
        return all(
            obj is not None
            for obj in (self.data, self.mesh, self.model, self.startup)
        )

    def summary(self) -> str:
        """Return a short human-readable build summary."""
        if not self.is_ready:
            return "InputBuilder: not yet built"
        return (
            f"InputBuilder summary\n"
            f"  workdir   : {self.workdir}\n"
            f"  sites     : {self.data.n_sites}\n"
            f"  freqs     : {self.data.n_frequencies}\n"
            f"  data pts  : {self.data.n_data}\n"
            f"  mesh      : {self.mesh.n_xcells} x {self.mesh.n_zcells} cells\n"
            f"  params    : {self.model.n_params}\n"
            f"  modes     : {self.config.modes}\n"
        )
