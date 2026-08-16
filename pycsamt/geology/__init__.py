# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt.geology — general-purpose geological domain knowledge.

This package holds earth-science concepts that are not themselves
electromagnetic: a resistivity-to-lithology classification engine, a
literature-compiled rock/fluid property table, ground-truth borehole logs,
and (as the package grows) structural-geology primitives. Nothing here
imports :class:`pycsamt.interp.ResistivityModel` or any other EM concept;
:mod:`pycsamt.interp` depends on this package, not the other way round.

For the interpretation workflow that turns an EM resistivity model into a
calibrated, geologically classified section — the part that *does* know
about resistivity models, boreholes-as-constraints, and misfit review —
see :mod:`pycsamt.interp` (:class:`~pycsamt.interp.ModelCalibrator` in
particular). ``pycsamt.interp.RockDatabase``, ``pycsamt.interp.Borehole``,
and the other classes below remain importable from :mod:`pycsamt.interp`
for backward compatibility; this package is their canonical home.

Quickstart
----------
>>> from pycsamt.geology import RockDatabase, Borehole, FaultTrace
>>>
>>> db = RockDatabase.default()
>>> db.classify(250.0).name
'Granite (weathered)'
>>>
>>> bh = Borehole.from_csv("boreholes/Bo.csv", x=1050.0)  # doctest: +SKIP
>>>
>>> fault = FaultTrace(x=500.0, dip_deg=70.0, downthrown_side="right")

Package layout
--------------
``lithology``
    :class:`RockDatabase`, :class:`RockEntry`, :class:`StratigraphicLog`,
    :class:`Layer` — resistivity-to-lithology classification engine only;
    the built-in table itself lives in ``rock_library``.
``rock_library``
    :data:`~rock_library.BUILTIN_ROCKS` — the literature-compiled rock/
    fluid resistivity table behind :meth:`RockDatabase.default`, kept
    separate so it can grow without touching classification logic.
``rock_providers``
    :class:`~rock_providers.RockPropertyProvider`,
    :class:`~rock_providers.LocalRockPropertyProvider`,
    :class:`~rock_providers.RemoteRockPropertyProvider` — pluggable
    rock-property sources behind :meth:`RockDatabase.from_url` and
    :meth:`RockDatabase.from_provider`, with local caching and
    fallback-to-default on failure.
``borehole``
    :class:`Borehole`, :class:`Interval` — ground-truth depth-interval
    data; readers for CSV and LAS 2.0.
``structural``
    :class:`StructuralMeasurement` (planar: strike/dip/dip-direction),
    :class:`LinearMeasurement` (linear: trend/plunge), and
    :class:`FaultTrace` (where a fault crosses the profile) — field
    structural evidence; :class:`StructuralModel` collects all three per
    profile with CSV I/O and nearest/within queries.

References
----------
.. [1] Palacky, G.J. (1988). Resistivity characteristics of geologic
   targets. In Nabighian, M.N. (ed.), *Electromagnetic Methods in Applied
   Geophysics*, Vol. 1, Society of Exploration Geophysicists.
.. [2] Telford, W.M., Geldart, L.P., Sheriff, R.E. (1990). *Applied
   Geophysics*, 2nd ed., Cambridge University Press.
"""

from .borehole import Borehole, Interval
from .lithology import Layer, RockDatabase, RockEntry, StratigraphicLog
from .structural import (
    FaultTrace,
    LinearMeasurement,
    StructuralMeasurement,
    StructuralModel,
)

__all__ = [
    # lithology
    "RockDatabase",
    "RockEntry",
    "Layer",
    "StratigraphicLog",
    # borehole
    "Borehole",
    "Interval",
    # structural
    "StructuralMeasurement",
    "LinearMeasurement",
    "FaultTrace",
    "StructuralModel",
]
