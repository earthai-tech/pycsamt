# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
Python interface to the Occam2DMT inversion workflow.

The :mod:`pycsamt.models.occam2d` subpackage builds, runs, reads, and
plots two-dimensional magnetotelluric and CSAMT inversions using the
Occam smooth-model approach [module-1]_, [module-2]_. The inversion seeks a model
that fits the data to a target normalized RMS while minimizing
roughness.

.. math::

    \phi_d =
    \sqrt{\\frac{1}{N}\sum_{i=1}^{N} r_i^2},
    \qquad
    \rho_i = 10^{m_i}.

Typical Workflow
----------------
1. Load EDI files with :class:`pycsamt.site.Sites`.
2. Build input files with :class:`InputBuilder`.
3. Run the Fortran executable with :class:`OccamRunner`.
4. Load results with :class:`InversionResult`.
5. Plot models, responses, pseudosections, and misfit curves.

Examples
--------
>>> from pycsamt.models.occam2d import (
...     InputBuilder,
...     OccamRunner,
...     InversionResult,
... )
>>> from pycsamt.site import Sites
>>> sites = Sites.from_any("edi")
>>> InputBuilder(sites, workdir="occam_run").build()
>>> OccamRunner("occam_run").run(target_misfit=1.0)
>>> result = InversionResult("occam_run")
>>> result.plot_model()

References
----------
.. [module-1] Constable, S. C., Parker, R. L., and Constable, C. G.,
   "Occam's inversion: A practical algorithm for generating
   smooth models from electromagnetic sounding data", Geophysics,
   52(3), 289-300, 1987.
.. [module-2] deGroot-Hedlin, C., and Constable, S.,
   "Occam's inversion to generate smooth, two-dimensional models
   from magnetotelluric data", Geophysics, 55(12), 1613-1624,
   1990.
"""

from .builder import InputBuilder
from .config import OccamConfig
from .data import OccamData
from .log import OccamLog
from .mesh import OccamMesh
from .model import OccamModel
from .plot import (
    PlotMisfit,
    PlotModel,
    PlotPseudo,
    PlotResponse,
    PlotResponseGrid,
    PlotSiteMisfit,
    PlotSounding1D,
    PlotStation1DFit,
    plot_station_1d_fit,
)
from .response import OccamResponse
from .results import InversionResult
from .runner import OccamRunner
from .startup import OccamIter, OccamStartup

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
    "PlotSounding1D",
    "PlotSiteMisfit",
    "PlotResponseGrid",
    "OccamConfig",
]
