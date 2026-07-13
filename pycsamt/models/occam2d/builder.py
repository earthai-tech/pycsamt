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
from typing import Union

from .base import OccamBase
from .config import OccamConfig
from .data import OccamData
from .doc import _occam_param_docs as _params
from .mesh import OccamMesh
from .model import OccamModel
from .startup import OccamStartup

PathLike = Union[str, Path]

__all__ = ["InputBuilder"]


class InputBuilder(OccamBase):
    def __init__(
        self,
        source,
        workdir: PathLike = ".",
        config: OccamConfig | None = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.source = source
        self.workdir = Path(workdir)
        self.config = config or OccamConfig()

        self.data: OccamData | None = None
        self.mesh: OccamMesh | None = None
        self.model: OccamModel | None = None
        self.startup: OccamStartup | None = None

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------
    def build(
        self,
        modes: list[str] | None = None,
        n_layers: int | None = None,
        cell_size: float | None = None,
        error_floor_rho: float | None = None,
        error_floor_phase: float | None = None,
        freq_min: float | None = None,
        freq_max: float | None = None,
        title: str = "pycsamt Occam2D run",
        **kwargs,
    ) -> InputBuilder:
        # Apply one-shot overrides
        cfg = self.config
        if modes is not None:
            cfg.modes = modes
        if n_layers is not None:
            cfg.n_layers = n_layers
        if cell_size is not None:
            cfg.cell_size_horizontal = cell_size
        if error_floor_rho is not None:
            cfg.error_floor_rho = error_floor_rho
        if error_floor_phase is not None:
            cfg.error_floor_phase = error_floor_phase
        if freq_min is not None:
            cfg.freq_min = freq_min
        if freq_max is not None:
            cfg.freq_max = freq_max

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
        """Return ``True`` when all build objects are populated."""
        return all(
            obj is not None
            for obj in (self.data, self.mesh, self.model, self.startup)
        )

    def summary(self) -> str:
        """Return a compact, human-readable build summary."""
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


InputBuilder.__doc__ = rf"""
Build the complete input set for an Occam2D inversion.

``InputBuilder`` is the main construction entry point for
the v2 Occam2D workflow. It takes one survey source and writes
the four
files consumed by the Occam2DMT Fortran program:

* ``OccamDataFile.dat``: observed data and uncertainty rows.
* ``Occam2DMesh``: finite-element mesh geometry.
* ``Occam2DModel``: mapping from mesh cells to parameters.
* ``Startup``: inversion controls and initial model vector.

The build chain follows Occam smooth inversion [1]_, [2]_.
The inversion later seeks the smoothest model that reaches a
target
normalized RMS misfit:

.. math::

    \mathrm{{RMS}} =
    \sqrt{{\frac{{1}}{{N}}\sum_{{i=1}}^N r_i^2}}.

``InputBuilder`` does not run the inversion. It prepares a
self-contained working directory that can be passed to
:class:`~pycsamt.models.occam2d.OccamRunner`.

Parameters
----------
{_params.data.source}
{_params.common.workdir}
{_params.common.config}
{_params.common.verbose}
{_params.common.logger}

Attributes
----------
source : object
    Original survey source passed to the builder. It is kept
    so repeated :meth:`build` calls can rebuild outputs with
    different one-shot overrides.
workdir : pathlib.Path
    Output directory where Occam files are written.
config : OccamConfig
    Mutable run configuration used by the build steps.
data : OccamData or None
    Populated after :meth:`build`. It stores station names,
    offsets, frequencies, data-type codes, datum values, and
    errors.
mesh : OccamMesh or None
    Finite-element mesh generated from ``data.offsets`` and
    mesh-related configuration values.
model : OccamModel or None
    Model-parameter definition generated from ``mesh``.
startup : OccamStartup or None
    Iteration-zero startup object generated from ``model``.

Notes
-----
Keyword overrides supplied to :meth:`build` update the stored
``config`` object before files are written. This makes later
calls convenient, but overrides persist unless the
caller restores the configuration.

See Also
--------
OccamData.from_edi
    Convert the survey source into Occam data rows.
OccamMesh.from_data
    Build the finite-element mesh from station offsets.
OccamModel.from_mesh
    Convert the mesh into an inversion parameter mapping.
OccamStartup.from_model
    Build the initial log10-resistivity parameter vector.
OccamRunner
    Execute the compiled Occam2DMT binary in ``workdir``.
InversionResult
    Load iteration, response, model, mesh, and log outputs.

Examples
--------
Build default TE/TM run files from a directory of EDI files:

>>> from pycsamt.models.occam2d import InputBuilder
>>> from pycsamt.site import Sites
>>> sites = Sites.from_dir("edi")
>>> builder = InputBuilder(sites, workdir="occam_run")
>>> builder.build()
>>> builder.is_ready
True

Build only TM data over a restricted frequency band:

>>> builder = InputBuilder(sites, workdir="occam_tm")
>>> builder.build(
...     modes=["TM"],
...     freq_min=0.1,
...     freq_max=1000.0,
...     title="TM-only Occam2D test",
... )

Use a configuration object for a finer mesh:

>>> from pycsamt.models.occam2d import OccamConfig
>>> cfg = OccamConfig(n_layers=36)
>>> cfg.cell_size_horizontal = 50.0
>>> builder = InputBuilder(sites, workdir="fine_run")
>>> builder.config = cfg
>>> builder.build(error_floor_rho=0.07)

Chain build and summary calls in a script:

>>> summary = InputBuilder(sites, "run").build().summary()
>>> "data pts" in summary
True

References
----------
.. [1] Constable, S. C., Parker, R. L., and Constable,
   C. G., "Occam's inversion: A practical algorithm for
   generating smooth models from electromagnetic sounding
   data", Geophysics, 52(3), 289-300, 1987.
.. [2] deGroot-Hedlin, C., and Constable, S.,
   "Occam's inversion to generate smooth, two-dimensional
   models from magnetotelluric data", Geophysics, 55(12),
   1613-1624, 1990.
"""

InputBuilder.build.__doc__ = rf"""
Build and write all required Occam2D input files.

The method executes the input-file pipeline in fixed order:
data, mesh, model, then startup. Each step depends on the
previous object, so failures usually identify the first bad
input in the chain.

The generated files are written using names stored in
``self.config``:

* ``config.data_file``
* ``config.mesh_file``
* ``config.model_file``
* ``config.startup_file``

One-shot arguments update the builder configuration before
writing
files. For example, ``n_layers=40`` updates
``self.config.n_layers`` and affects mesh, model, and startup
objects created by this call.

Parameters
----------
{_params.data.modes}
{_params.mesh.n_layers}
{_params.mesh.cell_size}
{_params.data.error_floor_rho}
{_params.data.error_floor_phase}
{_params.data.freq_min}
{_params.data.freq_max}
{_params.data.title}

Returns
-------
InputBuilder
    Same builder instance, populated with ``data``, ``mesh``,
    ``model``, and ``startup`` objects.

Raises
------
ValueError
    If the survey source cannot produce sites, frequencies,
    modes, station offsets, or model parameters.
OSError
    If the output directory or files cannot be written.

See Also
--------
OccamData.write
    Write the generated data file.
OccamMesh.write
    Write the generated mesh file.
OccamModel.write
    Write the generated model file.
OccamStartup.write
    Write the generated startup file.

Examples
--------
Build a standard run:

>>> from pycsamt.models.occam2d import InputBuilder
>>> from pycsamt.site import Sites
>>> sites = Sites.from_dir("edi")
>>> builder = InputBuilder(sites, workdir="run")
>>> builder.build(modes=["TE", "TM"], n_layers=30)

Build a shallow, high-frequency test:

>>> builder.build(
...     modes=["TE"],
...     freq_min=10.0,
...     freq_max=10000.0,
...     n_layers=18,
...     cell_size=75.0,
... )

Inspect the generated objects without rereading files:

>>> builder.data.n_sites > 0
True
>>> builder.mesh.n_xcells > 0
True
>>> builder.model.n_params == builder.startup.n_params
True
"""
