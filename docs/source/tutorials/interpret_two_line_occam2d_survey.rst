.. _tutorial_interpret_two_line_occam2d_survey:

Interpret A Two-Line Occam2D Survey With interp And geology
================================================================

:doc:`build_two_line_occam2d_survey` produced two real, independently
inverted resistivity sections -- Line A and Line B, each crossing the
same kind of dipping fault at a different angle. Both sections show only
a broad, smooth resistivity tilt; neither shows the fault itself as a
sharp feature, because smoothness-regularized Occam2D is not built to
reproduce a discrete offset that way. This page is what turns those two
smooth images into an actual geological interpretation: calibrate each
one against a pair of boreholes, classify the calibrated model into
lithology with :mod:`pycsamt.geology`, record the fault as independent
structural evidence rather than expecting to read it off the resistivity
image, and then -- honestly -- check the result against two more
boreholes that were *not* used for calibration.

Every real number below comes from the same real Occam2D results Part 1
produced (loaded fresh, not re-inverted) and from the same true earth
model, queried directly with
:meth:`~pycsamt.forward.grid2d.Grid2D.value_at` for borehole ground
truth -- the same function, the same coordinate convention, used
consistently from the forward model all the way through validation.

A project-specific rock table
---------------------------------

The built-in :meth:`~pycsamt.geology.RockDatabase.default` table is a
literature compilation, not tuned to any one site --
:doc:`../user_guide/geology/rock_database` already shows the same query
value classified differently by two different tables. Classifying this
survey's own true resistivities (80, 300, 3000 Ω·m) against the built-in
table returns "Fractured zone", "Granite (weathered)", and "Dolomite" --
geologically plausible names, but not the actual overburden / weathered
basement / fresh basement sequence this survey's own boreholes describe.
A small regional table, built the same way
:doc:`../user_guide/geology/rock_database` demonstrates, fixes that:

.. code-block:: pycon

   >>> from pycsamt.geology import RockDatabase, RockEntry
   >>> regional_db = RockDatabase([
   ...     RockEntry(name="Overburden", rho_min=30, rho_max=180,
   ...               color="#D4AC0D", source="Site-specific, this survey"),
   ...     RockEntry(name="Weathered basement", rho_min=180, rho_max=900,
   ...               color="#A9780C", source="Site-specific, this survey"),
   ...     RockEntry(name="Fresh basement", rho_min=900, rho_max=8000,
   ...               color="#4A4A4A", source="Site-specific, this survey"),
   ... ])
   >>> [regional_db.classify(rho).name for rho in (80.0, 300.0, 3000.0)]
   ['Overburden', 'Weathered basement', 'Fresh basement']

Boreholes, sampled honestly from the true model
------------------------------------------------------

Every borehole below is built the same way: read the *true* earth model
Part 1 built -- not the inversion, not a guess -- at a chosen profile
position, with :meth:`~pycsamt.forward.grid2d.Grid2D.value_at`, exactly
as a real drilling program would measure it. Four boreholes per line,
split two ways: two for calibration, two held out entirely for
validation later on this page.

.. code-block:: pycon

   >>> from pycsamt.geology import Borehole, Interval
   >>> from pycsamt.forward.synthetic import LayeredModel
   >>> from pycsamt.forward.grid2d import Grid2D

   >>> true_layers = LayeredModel(resistivity=[80.0, 300.0, 3000.0], thickness=[60.0, 260.0])
   >>> grid_a = Grid2D.layered_with_fault(
   ...     true_layers, fault_x_m=1200.0, apparent_dip_deg=65.0,
   ...     throw_m=700.0, downthrown_side="right",
   ...     nx=60, nz=50, x_max=2400.0, z_max=1600.0, n_stations=25, name="Line A",
   ... )

Scanning depth at a fixed ``x`` and recording every resistivity change
gives the real interval boundaries -- close to, but not exactly, the
nominal 60/320 m layer depths, since a real log is quantised by the
grid's own cell size just as a real logging tool has finite sample
spacing:

.. code-block:: pycon

   >>> import numpy as np
   >>> def true_intervals(grid, x_m, max_depth_m, dz=1.0):
   ...     zs = np.arange(0.0, max_depth_m + dz, dz)
   ...     vals = np.array([grid.value_at(x_m, z, chainage=True) for z in zs])
   ...     out, start = [], 0
   ...     for i in range(1, len(zs)):
   ...         if vals[i] != vals[i - 1]:
   ...             out.append((float(zs[start]), float(zs[i]), float(vals[i - 1])))
   ...             start = i
   ...     out.append((float(zs[start]), float(zs[-1]), float(vals[-1])))
   ...     return out

   >>> true_intervals(grid_a, 300.0, 450.0)
   [(0.0, 65.0, 80.0), (65.0, 321.0, 300.0), (321.0, 450.0, 3000.0)]

Wrapping each interval in this way, with the true resistivity carried
alongside the lithology name classified against ``regional_db``,
produces the actual :class:`~pycsamt.geology.Borehole` objects the
calibrator -- and, later, the validation step -- will use:

.. code-block:: pycon

   >>> def borehole_from_true(grid, name, x_m, max_depth_m):
   ...     ivs = [
   ...         Interval(top=top, bottom=bottom, resistivity=rho,
   ...                  lithology=regional_db.classify(rho).name)
   ...         for top, bottom, rho in true_intervals(grid, x_m, max_depth_m)
   ...     ]
   ...     return Borehole(name, x=x_m, intervals=ivs)

   >>> bh_a1 = borehole_from_true(grid_a, "A1", 300.0, 450.0)
   >>> bh_a2 = borehole_from_true(grid_a, "A2", 2100.0, 1100.0)
   >>> [iv.lithology for iv in bh_a1.intervals]
   ['Overburden', 'Weathered basement', 'Fresh basement']

``A1`` (footwall, x=300 m) only needs to reach 450 m to confirm fresh
basement; ``A2`` (downthrown, x=2100 m) needs 1100 m for the same
confirmation, because the same basement contact sits roughly 700 m
deeper there -- a direct, physical consequence of the fault throw, not
an arbitrary choice. Both calibration boreholes sit comfortably inside
the footwall or downthrown block at every depth drilled -- the fault
plane itself migrates with depth (``x = fault_x + z / tan(apparent_dip)``),
and a borehole placed too close to the surface trace can cross from one
side to the other partway down, which would make "the true log at this
x" ambiguous. A1, A2, and the two validation boreholes below were all
checked against that migration and kept clear of it; :doc:`../user_guide/geology/structural`
covers the general apparent-dip geometry this follows.

The two validation boreholes -- ``A3`` (footwall, x=900 m) and ``A4``
(downthrown, x=1800 m) -- are built exactly the same way, but are never
passed to the calibrator. Line B's four boreholes mirror this with its
own fault geometry (``downthrown_side="left"``): ``B1``/``B2`` for
calibration at x=300/1900 m, ``B3``/``B4`` held out at x=100/2100 m.

See every borehole at a glance
--------------------------------------

Before any of these boreholes touch the inversion, it helps to look at
all four on a line side by side -- exactly the raw evidence the
calibrator and the validation step will each draw on, with nothing from
Occam2D involved yet. :class:`~pycsamt.interp.plot.PlotBoreholeFence`
draws this the same way :class:`~pycsamt.interp.plot.PlotFenceDiagram`
draws classified model logs, but straight from
:class:`~pycsamt.geology.Interval` data:

.. code-dropdown:: ../../scripts/generate_tutorial_interpret_two_line_occam2d.py
   :language: python
   :pyobject: make_borehole_fence_figures
   :linenos:
   :title: View the borehole-fence figure source code

.. figure:: /images/tutorials/interpret_two_line_occam2d_survey/borehole_fence_a.png
   :alt: Side-by-side comparison of boreholes A1, A3, A4, and A2 on Line A, ordered by profile position.
   :width: 100%

   Line A's four boreholes, ordered by ``x``. ``A1``/``A3`` (footwall)
   only needed 450 m to confirm fresh basement; ``A4``/``A2``
   (downthrown) needed 1100 m for the same confirmation -- the fault
   throw, read directly off real drilling depths rather than off a
   model. The blank space below ``A1``/``A3`` is left genuinely blank,
   with no panel border drawn into it: nothing was drilled there, so
   there is nothing to frame.

.. figure:: /images/tutorials/interpret_two_line_occam2d_survey/borehole_fence_b.png
   :alt: Side-by-side comparison of boreholes B3, B1, B2, and B4 on Line B, ordered by profile position.
   :width: 100%

   Line B's four boreholes, ordered by ``x``. The pattern mirrors Line A
   with the sides swapped -- ``downthrown_side="left"`` here, so the
   deep contact sits under ``B3``/``B1`` instead.

Calibrate against the two boreholes
------------------------------------------

:class:`~pycsamt.interp.ModelCalibrator` blends borehole evidence into
the raw calculated model (CRM): cells within ``ptol`` (10% here) of a
nearby borehole's true resistivity are replaced outright; everything
else is classified from the regional table alone. Stations farther than
``max_borehole_distance`` (500 m, the default) from any borehole get no
direct TRES influence at all:

.. code-block:: pycon

   >>> from pathlib import Path
   >>> from pycsamt.interp import ResistivityModel, ModelCalibrator
   >>> from pycsamt.models.occam2d import InversionResult

   >>> result_a = InversionResult(Path("runs/line_a_occam2d"))
   >>> model_a = ResistivityModel.from_occam2d(result_a).clip_to_stations()
   >>> cal_a = ModelCalibrator(ptol=0.10, db=regional_db, verbose=False)
   >>> _ = cal_a.fit(model_a, [bh_a1, bh_a2])
   >>> nm_a = cal_a.calibrated_model()
   >>> nm_a.method
   'occam2d+calibrated'

``A1`` is 600 m from station ``L09`` (x=900 m) and ``A2`` is 300 m from
``L18`` (x=1800 m) -- both within reach in principle, though whether TRES
replacement actually *fires* at a given station depends on whether the
raw CRM resistivity there ever comes close enough to the borehole's true
value, which the validation section below shows is not guaranteed just
because a borehole is nearby.

See where calibration actually changed the model
----------------------------------------------------------

``cal_a.misfit_map()`` gives a per-column G (%) figure quantifying how
much a station's calibrated column differs from the raw CRM:

.. code-block:: pycon

   >>> mm_a = cal_a.misfit_map()
   >>> mm_a.shape
   (32, 50)
   >>> round(float(mm_a.min()), 2), round(float(mm_a.max()), 2)
   (8.46, 12.24)

Every one of Line A's 50 columns sits at 8-12% misfit -- not just the
two calibrated ones. That is real, and it is not a sign that something
went wrong: :meth:`~pycsamt.interp.ModelCalibrator.fit` runs
:meth:`~pycsamt.geology.RockDatabase.classify`-based autolayer
reclassification on *every* column outside ``max_borehole_distance`` of
a borehole, not just a straight TRES swap near the two that are close
enough. Autolayering alone -- collapsing a smooth, continuous CRM
column into a handful of discrete rock-database bins -- already moves
most cells' resistivity by a similar margin everywhere it runs, so the
misfit floor sits well above zero across the whole line, not only away
from the boreholes.

Occam2D's own mesh reaches almost 10 km depth for numerical-boundary
reasons that have nothing to do with the shallow geology this survey
cares about; plotted at full depth, the interesting top 1600 m would be
squeezed into a sliver of the figure. The figures below crop that
directly with ``depth_max`` (shown next); the same crop is available
without any topography draping through the new
:meth:`~pycsamt.interp.ResistivityModel.clip_to_depth`, useful whenever
a plain flat-depth model is all that is needed:

.. code-block:: pycon

   >>> model_a.n_z, model_a.clip_to_depth(1600.0).n_z
   (32, 19)

Rather than the flat, no-terrain :class:`~pycsamt.interp.plot.PlotCalibratedModel`
panel, the figure below drapes CRM, NM, and the misfit map over this
line's real terrain with :func:`~pycsamt.topo.plot_topo_section` -- the
same function :doc:`build_two_line_occam2d_survey` used for the
topography figure in Part 1 -- plus the shared ``inversion`` station
marker (a black-edged, white-filled downward triangle) and every
station's real elevation, resolvable directly from a
:class:`~pycsamt.interp.ResistivityModel` since CRM and NM both are
one. The misfit array is not a resistivity grid, so it goes through the
generic :func:`~pycsamt.topo.plot_topo_array` instead, which drapes and
colour-maps any 2-D field without assuming a log10(rho) transform.
Both functions accept ``depth_max`` directly, cropping out Occam2D's
own numerical-boundary mesh (it reaches almost 10 km depth for reasons
that have nothing to do with this survey's shallow geology) down to the
top 1600 m that is actually worth looking at. Because the misfit floor
sits at 8-12% almost everywhere rather than the more typical case of a
mostly-untouched line, the colour scale is widened to ``vmax=15``
rather than the tighter 0-10% range that would suit that typical case:

.. code-dropdown:: ../../scripts/generate_tutorial_interpret_two_line_occam2d.py
   :language: python
   :pyobject: make_calibration_effect_figures
   :linenos:
   :title: View the calibration-effect figure source code

.. figure:: /images/tutorials/interpret_two_line_occam2d_survey/calibration_effect_a.png
   :alt: CRM vs calibrated NM vs misfit percentage for Line A, draped over real topography, top 1600 metres.
   :width: 100%

   Line A, draped over its real terrain (the ridge peaking near
   ``L14``-``L18``). The misfit panel is *lowest* (pale yellow, near
   8%) right at ``L00``/``L24`` where the calibration boreholes sit,
   and highest (orange, near 12%) at the columns farthest from both --
   exactly the gradient autolayer-dominated misfit should produce.

.. figure:: /images/tutorials/interpret_two_line_occam2d_survey/calibration_effect_b.png
   :alt: CRM vs calibrated NM vs misfit percentage for Line B, draped over real topography, top 1600 metres.
   :width: 100%

   Line B shows the same pattern with its own boreholes and its own
   (independent) terrain: misfit peaks (dark red, above 13%) roughly
   midway between ``B1`` (x=300 m) and ``B2`` (x=1900 m) -- the point
   on the line farthest from either calibration borehole -- and falls
   toward both ends.

Classify into lithology and draw the fence diagram
---------------------------------------------------------

``cal_a.stratigraphic_logs()`` classifies every station's calibrated
column into merged lithology layers with the regional table, the same
mechanism :doc:`../user_guide/interpretation/lithology` covers in depth:

.. code-block:: pycon

   >>> logs_a = cal_a.stratigraphic_logs()
   >>> len(logs_a)
   25
   >>> log_l00 = [log for log in logs_a if log.station_name == "L00"][0]
   >>> [(round(l.top, 1), round(l.bottom, 1), l.lithology) for l in log_l00.layers]
   [(0.0, 76.9, 'Overburden'), (74.9, 322.1, 'Weathered basement'), (317.2, 465.8, 'Fresh basement'), (459.2, 659.2, 'Weathered basement'), (650.3, 10684.7, 'Fresh basement')]

Station ``L00`` sits directly at ``A1``'s own position, and the shallow
boundary (overburden ending around 77 m, against a true 65 m) is close
-- but the layer sequence oscillates once, back to "Weathered basement"
around 460-660 m before returning to "Fresh basement", rather than
settling cleanly. That is genuine classification noise, not a display
bug: the soft-replace step only pulls individual *cells* toward the true
value within ``ptol``; it does not smooth the resulting sequence
afterward, so a column that dips in and out of the 10% tolerance band
near the true transition depth can flicker between two lithologies for
a few cells before committing to one.

The fence diagram below renders every one of Line A's 25 classified
logs side by side, so that flicker -- and the overall footwall-to-
downthrown transition -- is visible across the whole profile at once.
:class:`~pycsamt.interp.plot.PlotFenceDiagram` takes an ``elevation_m``
array (one value per log, same order) and draws a real terrain strip
above the panels, with the same shared ``inversion`` station marker
used everywhere else on this page:

.. code-dropdown:: ../../scripts/generate_tutorial_interpret_two_line_occam2d.py
   :language: python
   :pyobject: make_fence_diagrams
   :linenos:
   :title: View the fence-diagram figure source code

.. figure:: /images/tutorials/interpret_two_line_occam2d_survey/fence_a.png
   :alt: Fence diagram of classified pseudo-stratigraphic logs along Line A, with a real terrain strip and station markers above.
   :width: 100%

   :class:`~pycsamt.interp.plot.PlotFenceDiagram` for Line A, all 25
   stations, with its real terrain strip on top. Grey (fresh basement)
   sits shallow through roughly ``L00``-``L06`` -- the footwall, close
   to ``A1`` -- then the boundary steps down through the fault-crossing
   zone before the profile settles into the deeper, downthrown pattern
   from about ``L17`` onward, close to ``A2``. Note that the strip uses
   station *index* along the x-axis to line up with the equal-width
   panels below, not real chainage.

Record the fault as structural evidence
----------------------------------------------

Nothing in the classified section above *places* the fault -- the
lithology boundary just described is a smoothed reflection of it, not a
located pick. :class:`~pycsamt.geology.FaultTrace` records the fault
directly, from what is actually known about it here (the same true
geometry :doc:`build_two_line_occam2d_survey` used to build the earth
model in the first place, standing in for what a real project would
learn from surface mapping, a resistivity offset, or a nearby drill
result):

.. code-block:: pycon

   >>> from pycsamt.geology import FaultTrace
   >>> fault_a = FaultTrace(
   ...     x=1200.0, dip_deg=65.0, downthrown_side="right", throw_m=700.0,
   ...     sense="normal", evidence="true model (synthetic control)",
   ... )
   >>> fault_a
   FaultTrace(x=1200.0 m, dip=65 deg, down=right, throw=700.0 m)

Line B's own fault trace is identical in kind, built from its own
geometry (``x=1400.0``, ``dip_deg=47.0``, ``downthrown_side="left"``) --
the same true fault system Part 1 introduced, crossed at a shallower
apparent dip. Drawing both lines' calibrated sections with the true
fault trace and all four boreholes overlaid makes the relationship
between them concrete -- and doing it over real terrain, the way the
calibration-effect figures above already do, means the fault trace and
every borehole need converting from flat ``(x, depth)`` into
terrain-following ``(x, elevation)`` first. :func:`~pycsamt.topo.interp_elev`
(the same function this page already relies on for every terrain
lookup) gives the local elevation at any along-profile ``x``; the fault
trace's own ``x = fault_x + z / tan(apparent_dip)`` relationship then
becomes ``elevation = interp_elev(x) - z`` at each depth sample, and
each borehole's line runs from ``interp_elev(bh.x)`` down to
``interp_elev(bh.x) - bh.max_depth``:

.. code-dropdown:: ../../scripts/generate_tutorial_interpret_two_line_occam2d.py
   :language: python
   :pyobject: make_structural_section_figure
   :linenos:
   :title: View the structural-overlay figure source code

.. figure:: /images/tutorials/interpret_two_line_occam2d_survey/structural_overlay.png
   :alt: Calibrated resistivity sections for Line A and Line B, draped over real topography, with the true fault trace, station markers, and all four boreholes overlaid.
   :width: 100%

   Both lines draped over their own real terrain, with the ``inversion``
   station marker and ``L00``-``L24`` labels from
   :func:`~pycsamt.topo.plot_topo_section` along the top. Solid lines
   and bold labels are calibration boreholes (``A1``/``A2``,
   ``B1``/``B2``); dotted lines and italic labels are the held-out
   validation boreholes -- both now hanging from the real terrain
   surface rather than a flat datum. The sharp, blocky colouring
   immediately around each calibration borehole is the soft-replace
   step visibly overriding the smooth Occam2D image with real TRES
   values; everywhere else keeps the raw inversion's own smooth
   gradient.

Validate against the held-out boreholes
------------------------------------------------

This is the step that actually tests whether calibration and
classification recovered something real, rather than just producing a
plausible-looking picture. ``A3`` and ``A4`` were never shown to the
calibrator:

.. code-block:: pycon

   >>> bh_a3 = borehole_from_true(grid_a, "A3", 900.0, 450.0)
   >>> bh_a4 = borehole_from_true(grid_a, "A4", 1800.0, 1100.0)
   >>> logs_by_station = {log.station_name: log for log in logs_a}

   >>> log_l09 = logs_by_station["L09"]  # nearest station to A3 (exact: same x)
   >>> bh_a3.intervals[0].bottom, log_l09.layers[0].bottom   # overburden base, true vs classified
   (65.0, 76.93395000000001)
   >>> bh_a3.intervals[-1].top    # true basement top
   321.0
   >>> next(l.top for l in log_l09.layers if l.lithology == "Fresh basement")
   2345.512975

The shallow overburden boundary is close either way (77 m classified
against 65 m true) -- resolvable because Occam2D's own shallow
resolution is good almost everywhere on this line, calibration borehole
nearby or not. The basement top is not close at all: 2346 m classified
against 321 m true, off by a factor of more than seven. ``A4``
(downthrown, 300 m from calibration borehole ``A2``, well inside
``max_borehole_distance``) does not fare better:

.. code-block:: pycon

   >>> log_l18 = logs_by_station["L18"]  # nearest station to A4 (exact: same x)
   >>> bh_a4.intervals[0].bottom, log_l18.layers[0].bottom
   (769.0, 779.6850000000001)
   >>> bh_a4.intervals[-1].top
   1025.0
   >>> next(l.top for l in log_l18.layers if l.lithology == "Fresh basement")
   5849.26135

Proximity to a calibration borehole did not save this one either. The
reason is visible directly in the section figure above: the raw
Occam2D image on the downthrown side never actually recovers a
resistivity anywhere close to 3000 Ω·m at any depth this model
resolves, so the soft-replace step -- which only fires within ``ptol``
of the true value -- never has a matching cell to grab onto near the
true basement depth, calibration borehole nearby or not. Depth of
investigation, not distance to the nearest borehole, is the binding
constraint here.

Line B's footwall validation borehole, ``B4`` (x=2100 m, 300 m from
calibration borehole ``B2``), tells the opposite story -- because on
Line B's *footwall* side the true basement sits shallow (321 m, the
same depth ``A1``/``A3`` see on Line A), well within what this model
actually resolves:

.. code-block:: pycon

   >>> grid_b = Grid2D.layered_with_fault(
   ...     true_layers, fault_x_m=1400.0, apparent_dip_deg=47.0,
   ...     throw_m=700.0, downthrown_side="left",
   ...     nx=60, nz=50, x_max=2400.0, z_max=1600.0, n_stations=25, name="Line B",
   ... )
   >>> bh_b4 = borehole_from_true(grid_b, "B4", 2100.0, 450.0)
   >>> bh_b4.intervals[-1].top
   321.0

   >>> result_b = InversionResult(Path("runs/line_b_occam2d"))
   >>> model_b = ResistivityModel.from_occam2d(result_b).clip_to_stations()
   >>> bh_b1 = borehole_from_true(grid_b, "B1", 300.0, 1100.0)
   >>> bh_b2 = borehole_from_true(grid_b, "B2", 1900.0, 450.0)
   >>> cal_b = ModelCalibrator(ptol=0.10, db=regional_db, verbose=False)
   >>> _ = cal_b.fit(model_b, [bh_b1, bh_b2])
   >>> logs_b = {log.station_name: log for log in cal_b.stratigraphic_logs()}
   >>> log_l21 = logs_b["L21"]  # nearest station to B4
   >>> next(l.top for l in log_l21.layers if l.lithology == "Fresh basement")
   317.175025

317 m classified against 321 m true -- within 4 m, the best match on
either line. Four validation boreholes, two clear outcomes: the shallow
overburden contact classifies well everywhere on both lines regardless
of calibration proximity; the basement contact classifies well only
where it is shallow enough for Occam2D to have actually resolved it, and
badly wherever it is not, independent of how close a calibration
borehole happens to sit. Reporting both, rather than only the flattering
``B4`` result, is the entire point of holding boreholes out in the first
place.

What this adds up to
------------------------

The raw inversion result from :doc:`build_two_line_occam2d_survey`
showed a smooth resistivity gradient and nothing else -- no fault, no
named lithology, no way to tell a well-resolved boundary from a poorly
resolved one. Everything on this page came from combining that result
with independent evidence: :class:`~pycsamt.geology.RockDatabase` turned
resistivity into lithology names a project actually recognises;
:class:`~pycsamt.interp.ModelCalibrator` pulled the model toward real
measurements near the two calibration boreholes, and its ``misfit_map``
showed exactly how far that influence actually reaches;
:class:`~pycsamt.geology.FaultTrace` placed the fault using what is
actually known about it, rather than reading a position off a smoothed
image that was never going to show one; and the two held-out boreholes
turned "this looks reasonable" into an honest, checkable claim about
where the result can and cannot be trusted. The borehole fence,
calibration-effect, classified fence, and structural-overlay figures
along the way are all built the same way: real
:mod:`pycsamt.interp.plot` classes and, wherever real terrain and
station markers belong in the picture, real :mod:`pycsamt.topo`
functions -- called directly on this page's own real objects, never a
hand-drawn illustration.

See Also
--------

:doc:`build_two_line_occam2d_survey`
    Builds the two real Occam2D lines this page interprets.

:doc:`../user_guide/interpretation/workflow`
    The general calibrate-classify-review sequence this page follows,
    on a smaller synthetic fixture.

:doc:`../user_guide/geology/rock_database`
    Why a project-specific rock table changes the classified answer even
    when the resistivity does not.

:doc:`../user_guide/geology/borehole`
    ``Borehole``/``Interval`` mechanics used throughout this page.

:doc:`../user_guide/geology/structural`
    ``FaultTrace``'s apparent-dip geometry, covered in full.

:doc:`../user_guide/interpretation/uncertainty`
    Where to go next for a quantitative, sampled treatment of the same
    kind of uncertainty this page's validation section found by hand.
