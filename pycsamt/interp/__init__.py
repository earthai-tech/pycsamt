# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt.interp — geophysical-geological interpretation for EM methods.

This package converts 2-D EM resistivity models (from any inversion code
— Occam2D, ModEM, AI-based, or generic arrays) into geological
deliverables: pseudo-stratigraphic logs, borehole-calibrated models, and
industry-standard exports.

It supersedes the ``geodrill`` module from pyCSAMT v1, extending its
scope beyond OCCAM2D and CSAMT to cover MT, AMT, CSAMT, and CSEM.

Quickstart
----------
>>> from pycsamt.interp import ResistivityModel, ModelCalibrator
>>> from pycsamt.interp import export, plot as iplot
>>>
>>> # Adapt any inversion result
>>> model = ResistivityModel.from_occam2d(result)          # doctest: +SKIP
>>>
>>> # Optional: calibrate with borehole data
>>> from pycsamt.interp import Borehole
>>> bh = Borehole.from_csv("boreholes/Bo.csv", x=1050.0)  # doctest: +SKIP
>>> cal = ModelCalibrator(ptol=0.10).fit(model, [bh])      # doctest: +SKIP
>>>
>>> # Build stratigraphic logs
>>> logs = cal.stratigraphic_logs()                        # doctest: +SKIP
>>>
>>> # Export
>>> export.to_oasis_montaj_xyz(logs, "profile.xyz")        # doctest: +SKIP
>>> export.to_las(logs[0], "S17.las")                      # doctest: +SKIP
>>> export.to_csv(logs, "all_stations.csv")                # doctest: +SKIP
>>>
>>> # Plot
>>> fig = iplot.PlotStratigraphicLog(logs[0]).plot()       # doctest: +SKIP
>>> fig = iplot.PlotFenceDiagram(logs).plot()              # doctest: +SKIP
>>> fig = iplot.PlotCalibratedModel(model, cal.calibrated_model(),
...                                 cal.misfit_map()).plot()  # doctest: +SKIP

Package layout
--------------
``_base``
    :class:`ResistivityModel` — unified model container; adapters for
    Occam2D, ModEM, AI outputs, and raw numpy arrays.
``borehole``
    :class:`Borehole`, :class:`Interval` — ground-truth depth-interval
    data; readers for CSV and LAS 2.0.
``calibrate``
    :class:`ModelCalibrator` — generalised NM construction (Kouadio
    et al. 2022), borehole-constrained or database-only.
``lithology``
    :class:`RockDatabase`, :class:`StratigraphicLog`, :class:`Layer` —
    resistivity-to-lithology classification; built-in 25-rock EM
    database.
``export``
    :func:`~export.to_oasis_montaj_xyz`,
    :func:`~export.to_las`,
    :func:`~export.to_csv`,
    :func:`~export.to_vtk` — industry-format writers.
``plot``
    :class:`~plot.PlotStratigraphicLog`,
    :class:`~plot.PlotFenceDiagram`,
    :class:`~plot.PlotCalibratedModel` — matplotlib-based visuals.

References
----------
.. [1] Kouadio, K.L. et al. (2022). pyCSAMT: An alternative Python
   toolbox for groundwater exploration using controlled source
   audio-frequency magnetotelluric. *Journal of Applied Geophysics*,
   201, 104647. https://doi.org/10.1016/j.jappgeo.2022.104647
"""

from ._base      import ResistivityModel
from .borehole   import Borehole, Interval
from .calibrate  import ModelCalibrator
from .hydro      import (
    AquiferZone,
    HydroGeophysicalModel,
    HydroInterpreter,
    HydroUnit,
)
from .lithology  import RockDatabase, RockEntry, Layer, StratigraphicLog
from .           import export
from .           import plot

__all__ = [
    # core model
    "ResistivityModel",
    # borehole
    "Borehole",
    "Interval",
    # calibration
    "ModelCalibrator",
    # hydrogeology
    "AquiferZone",
    "HydroGeophysicalModel",
    "HydroInterpreter",
    "HydroUnit",
    # lithology
    "RockDatabase",
    "RockEntry",
    "Layer",
    "StratigraphicLog",
    # sub-modules
    "export",
    "plot",
]
