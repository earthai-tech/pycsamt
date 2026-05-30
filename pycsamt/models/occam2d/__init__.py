# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.models.occam2d
======================

Python interface to the Occam2DMT v3.0 inversion code.

Typical workflow
----------------
1. Load EDI files via ``Sites`` or ``EDICollection`` (single source of truth).
2. Build all Occam input files with ``InputBuilder``.
3. Run the inversion with ``OccamRunner`` (calls the compiled Fortran binary).
4. Load and visualise results with ``InversionResult``.

Quick start
-----------
>>> from pycsamt.models.occam2d import InputBuilder, OccamRunner, InversionResult
>>> from pycsamt.site import Sites
>>>
>>> sites   = Sites.from_dir("path/to/edi/")
>>> builder = InputBuilder(sites, workdir="occam_run/")
>>> builder.build(modes=["TE", "TM"], n_layers=30, cell_size=100.)
>>>
>>> runner  = OccamRunner(workdir="occam_run/")
>>> runner.run(max_iter=100, target_misfit=1.0)
>>>
>>> result  = InversionResult(workdir="occam_run/")
>>> result.plot_model()

Fortran source
--------------
The Occam2DMT v3.0 Fortran source lives in ``_source/``.
Compile once with::

    cd pycsamt/models/occam2d/_source
    make FC90=gfortran

or let ``OccamRunner`` compile it automatically on first use.

References
----------
- Constable, Parker & Constable (1987) Geophysics 52, 289-300.
- deGroot-Hedlin & Constable (1990) Geophysics 55, 1613-1624.
- http://marineemlab.ucsd.edu/Projects/Occam/2DMT/
"""

from .builder  import InputBuilder
from .runner   import OccamRunner
from .results  import InversionResult
from .data     import OccamData
from .mesh     import OccamMesh
from .model    import OccamModel
from .startup  import OccamStartup, OccamIter
from .response import OccamResponse
from .log      import OccamLog
from .plot     import PlotModel, PlotResponse, PlotPseudo, PlotMisfit
from .config   import OccamConfig

__all__ = [
    "InputBuilder",
    "OccamRunner",
    "InversionResult",
    "OccamData",
    "OccamMesh",
    "OccamModel",
    "OccamStartup",
    "OccamIter",
    "OccamResponse",
    "OccamLog",
    "PlotModel",
    "PlotResponse",
    "PlotPseudo",
    "PlotMisfit",
    "OccamConfig",
]
