"""Per-backend converters into :class:`pycsamt.format.schema.PCSFModel`.

Occam2D (Phase 2), ModEM 3-D (Phase 3), and MARE2DEM (Phase 4) are
implemented, per ``PYCSAMT-PCSF-INVERSION-FORMAT-PLAN.md``. DUHI has no
separate adapter — its own output only becomes a final resistivity model
once folded back into an Occam2D run, so it rides
:func:`occam2d_to_pcsf`.

:mod:`pycsamt.format.adapters.generic` is a fourth, solver-agnostic path
for any *other* AI/DL inversion result (a UNet, a GCN, a ResNet, or a
third party's own model) whose output is itself the final resistivity
model — no dependency on any pycsamt result class, only plain arrays.
"""

from .generic import grid2d_to_pcsf, grid3d_to_pcsf, mesh_to_pcsf
from .mare2dem import mare2dem_to_pcsf
from .modem3d import modem3d_to_pcsf
from .occam2d import occam2d_to_pcsf

__all__ = [
    "occam2d_to_pcsf",
    "modem3d_to_pcsf",
    "mare2dem_to_pcsf",
    "grid2d_to_pcsf",
    "grid3d_to_pcsf",
    "mesh_to_pcsf",
]
