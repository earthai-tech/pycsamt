.. _geology_structural:

Structural measurements
========================

:class:`~pycsamt.geology.StructuralMeasurement`,
:class:`~pycsamt.geology.LinearMeasurement`, and
:class:`~pycsamt.geology.FaultTrace` record field structural evidence
against a profile position -- the same role
:class:`~pycsamt.geology.Borehole` plays for lithology, but for strike,
dip, and fault geometry instead of depth intervals.
:class:`~pycsamt.geology.StructuralModel` collects all three per
profile. None of the four are electromagnetic, and none of them place
themselves anywhere except along the profile -- see :doc:`concepts`
for why ``x`` is always a profile-relative metre value here, never a
latitude/longitude pair, and for how ``clone()``/``update()``
re-validate an edited ``StructuralMeasurement`` (they do, because it
overrides ``validate()``; not every class in this package does).

Recording a planar measurement
----------------------------------

A planar feature -- bedding, foliation, a joint, a fault plane -- is
recorded as strike, dip, and dip direction, as a compass and
clinometer would actually report it:

.. code-block:: pycon

   >>> from pycsamt.geology import StructuralMeasurement
   >>> bedding = StructuralMeasurement(
   ...     x=500.0, kind="bedding",
   ...     strike_deg=45.0, dip_deg=30.0, dip_direction_deg=135.0,
   ... )
   >>> bedding
   StructuralMeasurement(x=500.0 m, 'bedding', 45/30->135)

``kind`` is free text -- ``'bedding'``, ``'foliation'``, ``'joint'``,
``'cleavage'``, ``'contact'``, ``'fault_plane'``, ``'unconformity'``,
or a project-specific label -- matching how
:attr:`~pycsamt.geology.borehole.Interval.lithology` is free text
rather than an enforced enumeration. ``strike_deg`` and
``dip_direction_deg`` are compass bearings in ``[0, 360)`` degrees, not
reduced modulo 180: this module intentionally does not reuse MT
geoelectric strike's axial ``(-90, 90]`` convention or the axial
mod-180 convention :mod:`pycsamt.ai.geology.lenses` uses for synthetic
ellipse geometry, because a raw compass reading is directed, not
axial, and reducing it early would throw that information away before
it can be cross-checked.

That cross-check is why ``dip_direction_deg`` is stored at all rather
than leaving strike to imply it through the right-hand rule:
:attr:`~pycsamt.geology.StructuralMeasurement.dip_azimuth_ok` requires
``dip_direction_deg`` to sit within ``dip_direction_tolerance_deg``
(20 degrees by default) of ``strike_deg + 90`` or ``strike_deg - 90``,
which catches a transposed field-notebook entry that recording only
one of the two numbers never would:

.. code-block:: pycon

   >>> bedding.dip_azimuth_ok
   True
   >>> StructuralMeasurement(
   ...     x=500.0, kind="bedding",
   ...     strike_deg=45.0, dip_deg=30.0, dip_direction_deg=200.0,
   ... )
   Traceback (most recent call last):
       ...
   ValueError: dip_direction_deg (200.0) is not within 20.0 deg of strike_deg (45.0) +/- 90 -- check for a transposed strike/dip-direction reading.

Drawing the strike line, the two acceptance wedges, and both the
accepted and rejected ``dip_direction_deg`` values on a compass rose
makes the geometry concrete:

.. code-dropdown:: ../../../scripts/generate_user_guide_geology_structural_figures.py
   :language: python
   :pyobject: make_structural_measurement_geometry
   :linenos:
   :title: View the strike/dip-direction geometry figure source code

.. figure:: /images/user_guide/geology/structural_measurement_geometry.png
   :alt: Compass diagram showing the strike line, the two dip-direction acceptance wedges, and an accepted versus a rejected dip direction.
   :width: 65%

   The black line is the strike (drawn both directions, since a strike
   line has no single sense). The green wedges are the only
   ``dip_direction_deg`` values ``dip_azimuth_ok`` accepts -- 20 degrees
   either side of ``strike +/- 90``. ``135`` lands inside a wedge;
   ``200`` does not, which is exactly why the second construction above
   raises.

For anyone who prefers recording only dip direction and dip -- the
common two-number field style --
:meth:`~pycsamt.geology.StructuralMeasurement.from_right_hand_rule`
derives strike as ``dip_direction_deg - 90`` under the right-hand-rule
convention (facing along strike, dip direction to your right), and
produces an identical, already-consistent measurement:

.. code-block:: pycon

   >>> rhr = StructuralMeasurement.from_right_hand_rule(
   ...     x=500.0, kind="bedding",
   ...     dip_direction_deg=135.0, dip_deg=30.0,
   ... )
   >>> rhr.strike_deg == bedding.strike_deg
   True

Recording a linear measurement
----------------------------------

A linear feature -- a fold axis, a lineation, a slickenline -- has no
dip direction to cross-check against, so
:class:`~pycsamt.geology.LinearMeasurement` is simpler: just a
compass trend and a plunge below horizontal.

.. code-block:: pycon

   >>> from pycsamt.geology import LinearMeasurement
   >>> fold_axis = LinearMeasurement(
   ...     x=500.0, kind="fold_axis", trend_deg=210.0, plunge_deg=15.0,
   ... )
   >>> fold_axis
   LinearMeasurement(x=500.0 m, 'fold_axis', 210/15)

Fault traces
---------------

:class:`~pycsamt.geology.FaultTrace` is coarser than the two leaf
measurements above -- it is not a stereonet reading but a statement
about where a fault crosses *this* profile, which side dropped, and by
how much:

.. code-block:: pycon

   >>> from pycsamt.geology import FaultTrace
   >>> fault = FaultTrace(
   ...     x=500.0, dip_deg=70.0, downthrown_side="right",
   ...     sense="normal", throw_m=12.0, evidence="resistivity offset",
   ... )
   >>> fault
   FaultTrace(x=500.0 m, dip=70 deg, down=right, throw=12.0 m)

``dip_deg`` here is the *apparent* dip in the 2-D section, not
necessarily the fault's true 3-D dip -- the angle a single EM profile
can actually constrain, unless the line happens to run perpendicular
to strike. Pass the independently known ``strike_deg`` (from surface
mapping or a borehole) alongside it when the true attitude matters;
``FaultTrace`` keeps the two separate rather than conflating them.
``throw_m`` is always non-negative -- direction is carried entirely by
``downthrown_side``, so a negative throw would be redundant with, and
could contradict, that field:

.. code-block:: pycon

   >>> FaultTrace(x=500.0, dip_deg=70.0, downthrown_side="right", throw_m=-5.0)
   Traceback (most recent call last):
       ...
   ValueError: throw_m (-5.0) must be >= 0; direction is carried by downthrown_side, not the sign of throw_m.

This is the class meant to plug into the structural-continuity
question :doc:`../interpretation/workflow` already asks during misfit
review -- "do apparent boundary offsets align with known structures?"
-- with an actual record rather than an unbacked checklist item:
``evidence`` names what supports the pick (``'resistivity offset'``,
``'borehole'``, ``'surface mapping'``), and ``sense``
(``'normal'``, ``'reverse'``, ``'strike_slip'``, or the default
``'unknown'``) records the kinematics where known.

Collecting evidence along a profile
---------------------------------------

:class:`~pycsamt.geology.StructuralModel` holds every planar
measurement, linear measurement, and fault trace for one profile
together:

.. code-block:: pycon

   >>> from pycsamt.geology import StructuralModel
   >>> model = StructuralModel()
   >>> model.add_planar(StructuralMeasurement(
   ...     x=200.0, kind="bedding",
   ...     strike_deg=40.0, dip_deg=25.0, dip_direction_deg=130.0,
   ... ))
   >>> model.add_planar(StructuralMeasurement(
   ...     x=650.0, kind="bedding",
   ...     strike_deg=50.0, dip_deg=35.0, dip_direction_deg=140.0,
   ... ))
   >>> model.add_fault(fault)
   >>> model.add_fault(FaultTrace(
   ...     x=900.0, dip_deg=60.0, downthrown_side="left",
   ...     evidence="surface mapping",
   ... ))
   >>> model
   StructuralModel(2 planar, 0 linear, 2 faults)

Drawn as an actual cross section -- fault traces tilted toward their
``downthrown_side`` at their apparent ``dip_deg``, bedding measurements
as strike/dip markers at the surface -- this is the picture
``StructuralModel`` is meant to back the structural-continuity review
mentioned above with:

.. code-dropdown:: ../../../scripts/generate_user_guide_geology_structural_figures.py
   :language: python
   :pyobject: make_structural_evidence_section
   :linenos:
   :title: View the structural evidence section figure source code

.. figure:: /images/user_guide/geology/structural_evidence_section.png
   :alt: Cross-section view of two fault traces and two bedding measurements along the profile.
   :width: 100%

   The normal fault at ``x=500`` dips toward its downthrown (right)
   side; the fault of unknown sense at ``x=900`` dips toward its
   downthrown (left) side. Both bedding markers sit at the surface
   (``z=None`` was not overridden), labelled with their recorded
   strike/dip/dip-direction.

:meth:`~pycsamt.geology.StructuralModel.within` restricts a model to a
profile span -- useful for reviewing only the segment around one
inversion panel or station cluster -- and returns a new
``StructuralModel`` rather than mutating in place:

.. code-block:: pycon

   >>> model.within(0.0, 600.0)
   StructuralModel(1 planar, 0 linear, 1 faults)

:meth:`~pycsamt.geology.StructuralModel.nearest` answers "what is the
closest piece of evidence of this kind to this position," with an
optional ``max_distance`` so a distant match is not returned silently
as if it were locally relevant:

.. code-block:: pycon

   >>> model.nearest(520.0, kind="faults")
   FaultTrace(x=500.0 m, dip=70 deg, down=right, throw=12.0 m)
   >>> model.nearest(520.0, kind="faults", max_distance=5.0) is None
   True
   >>> model.nearest(300.0, kind="planar")
   StructuralMeasurement(x=200.0 m, 'bedding', 40/25->130)

Like :class:`~pycsamt.geology.Borehole` and
:class:`~pycsamt.geology.RockDatabase`,
:meth:`~pycsamt.geology.StructuralModel.from_csv` loads field data
from disk -- here, up to three independent CSV files, one per evidence
type, so a project can maintain separate spreadsheets for planar
readings, linear readings, and fault picks and combine only the ones
that exist:

.. code-block:: text

   x,dip_deg,downthrown_side,sense,throw_m,evidence
   500.0,70.0,right,normal,12.0,resistivity offset
   900.0,60.0,left,unknown,,surface mapping

.. code-block:: pycon

   >>> model_csv = StructuralModel.from_csv(faults_path="structure/faults.csv")
   >>> model_csv
   StructuralModel(0 planar, 0 linear, 2 faults)
   >>> for f in model_csv.faults:
   ...     print(f)
   FaultTrace(x=500.0 m, dip=70 deg, down=right, throw=12.0 m)
   FaultTrace(x=900.0 m, dip=60 deg, down=left, throw=?)

The second row's blank ``throw_m`` becomes ``None`` (unknown), not
zero -- ``from_csv`` treats an empty field the same way
:meth:`~pycsamt.geology.Borehole.from_csv` and
:meth:`~pycsamt.geology.RockDatabase.from_csv` treat their own
optional columns elsewhere in this package. ``planar_path`` and
``linear_path`` work the same way, with the column sets from
`Recording a planar measurement`_ and `Recording a linear measurement`_
respectively; any path left as ``None`` simply yields an empty list
for that evidence type rather than an error.

Where to go next
-------------------

:doc:`../interpretation/workflow` is where fault and structural
evidence actually gets weighed against inversion results during misfit
review. :doc:`rock_database` and :doc:`borehole` cover the other two
data families in this package.
