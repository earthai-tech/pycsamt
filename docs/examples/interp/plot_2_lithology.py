r"""
Lithology: rock classification and logs
=======================================

The first interpretation step is turning resistivity numbers into named
rock units. :mod:`pycsamt.interp.lithology` provides a
:class:`~pycsamt.interp.RockDatabase` that classifies a resistivity into a
lithology, and a :class:`~pycsamt.interp.lithology.StratigraphicLog` that
stacks those classifications into a borehole-style column. This example
classifies the synthetic section and draws single-station and multi-station
(fence) logs.
"""

# %%
# The rock database
# -----------------
# :meth:`RockDatabase.default <pycsamt.interp.RockDatabase.default>` ships a
# ready table of rock types with resistivity ranges. ``classify`` maps a
# value to its most likely rock.

from pycsamt.interp import RockDatabase

# Use the fence diagram (2nd figure) as the card thumbnail.
# sphinx_gallery_thumbnail_number = 2

db = RockDatabase.default()
for rho in (8.0, 40.0, 200.0, 1500.0):
    print(f"{rho:7.0f} ohm-m  ->  {db.classify(rho).name}")

# %%
# From model to stratigraphic logs
# --------------------------------
# The quickest way to get per-station logs is the hydro interpreter, which
# classifies every sounding and exposes the resulting
# :class:`~pycsamt.interp.lithology.StratigraphicLog` objects on ``.logs``.
# (The next example covers the hydrogeology it adds on top.)

from _interp_data import demo_model

from pycsamt.interp import HydroInterpreter

rm = demo_model()
hydro = HydroInterpreter(water_table_depth=20.0, aquifer_range=(30.0, 300.0),
                         clay_max=20.0, min_zone_thickness=8.0).fit(rm)
logs = hydro.logs
_mid = logs[len(logs) // 2]
print(f"{len(logs)} logs, e.g. {_mid.station_name!r} with "
      f"{len(_mid.layers)} layers")

# %%
# A single stratigraphic log
# --------------------------
# :class:`~pycsamt.interp.plot.PlotStratigraphicLog` draws the classic
# two-track column: hatched lithology on the left, the resistivity curve on
# the right. This is the figure you would hand a hydrogeologist.

from pycsamt.interp.plot import (
    PlotFenceDiagram,
    PlotStratigraphicLog,
)

PlotStratigraphicLog(logs[len(logs) // 2]).plot()

# %%
# A fence diagram across the line
# -------------------------------
# :class:`~pycsamt.interp.plot.PlotFenceDiagram` places every station's log
# side by side at its true distance, so laterally-continuous units (here the
# aquifer and clay) line up into correlatable layers — the standard
# cross-section deliverable.

PlotFenceDiagram(logs).plot()

# %%
# **Reading it.** The fence shows the four-unit sequence holding laterally
# with a gently undulating aquifer base, exactly the structure built into
# the model. Where a log's colours break the trend, that station's sounding
# is worth revisiting before trusting the correlation.
