.. _geology_borehole:

Borehole logs
=============

:doc:`../interpretation/workflow` already shows the everyday path for
a :class:`~pycsamt.geology.Borehole`: load one from CSV, give it a
profile position, and hand it to
:class:`~pycsamt.interp.ModelCalibrator` as a ground-truth constraint.
This page does not repeat that walkthrough. It covers what a
``Borehole`` and its :class:`~pycsamt.geology.borehole.Interval`
records actually offer beyond that one call -- querying a log
directly, building one without a CSV file, and reading LAS 2.0 well
logs, which :doc:`../interpretation/workflow` mentions but does not
demonstrate.

Querying a log directly
--------------------------

A :class:`~pycsamt.geology.borehole.Interval` is a leaf value object
in the sense introduced in :doc:`concepts`: ``bottom`` must exceed
``top``, checked once in ``__post_init__`` (see :doc:`concepts` for
why editing one after construction via ``clone()``/``update()`` does
not re-check this). ``thickness`` and ``contains()`` are the two
things every other method on this page is built from:

.. code-block:: pycon

   >>> from pycsamt.geology import Borehole, Interval
   >>> intervals = [
   ...     Interval(top=0.0, bottom=8.0, lithology="lateritic soil", resistivity=450.0),
   ...     Interval(top=8.0, bottom=31.0, lithology="clayey sand", resistivity=42.0),
   ...     Interval(top=31.0, bottom=67.0, lithology="weathered granite", resistivity=185.0),
   ...     Interval(top=67.0, bottom=120.0, lithology="fresh granite", resistivity=3200.0),
   ... ]
   >>> intervals[0].thickness
   8.0
   >>> intervals[0].contains(4.0), intervals[0].contains(8.0)
   (True, False)

A :class:`~pycsamt.geology.Borehole` collects them and is just as
usable built directly from a list as it is loaded from a file --
useful for a log entered by hand, generated from a database query, or
assembled in a test:

.. code-block:: pycon

   >>> bh = Borehole("BH01", x=500.0, intervals=intervals, collar_elevation=238.4)
   >>> bh
   Borehole('BH01', x=500.0 m, 4 intervals, depth=120.0 m)
   >>> bh.min_depth, bh.max_depth
   (0.0, 120.0)

``contains()`` is a half-open interval (``top <= z < bottom``), so
every depth belongs to exactly one interval except the log's own
lower bound. ``interval_at_depth()``, ``tres_at_depth()``, and
``lithology_at_depth()`` look a single depth up against it, returning
``None`` -- not raising -- outside the logged range:

.. code-block:: pycon

   >>> bh.interval_at_depth(50.0)
   Interval(top=31.0, bottom=67.0, lithology='weathered granite', resistivity=185.0)
   >>> bh.tres_at_depth(50.0)
   185.0
   >>> bh.lithology_at_depth(50.0)
   'weathered granite'
   >>> bh.tres_at_depth(500.0) is None
   True

``tres_column()`` vectorizes the same lookup over an array of depths
-- exactly what :class:`~pycsamt.interp.ModelCalibrator` uses
internally to compare a model column against this log, cell by cell,
before blending:

.. code-block:: pycon

   >>> import numpy as np
   >>> z = np.array([4.0, 20.0, 50.0, 100.0, 150.0])
   >>> bh.tres_column(z)
   array([ 450.,   42.,  185., 3200.,   nan])

The last depth (150 m) is below ``max_depth`` (120 m) and comes back
``nan`` rather than extrapolated -- a calibrator or any other consumer
of ``tres_column()`` should treat ``nan`` as "no ground truth here",
not as zero resistivity. Drawing the log as a resistivity-versus-depth
step curve and marking exactly where those five ``z`` values land
makes the boundary behaviour easy to see at a glance:

.. code-dropdown:: ../../../scripts/generate_user_guide_geology_borehole_figures.py
   :language: python
   :pyobject: make_borehole_tres_sampling
   :linenos:
   :title: View the tres_column sampling figure source code

.. figure:: /images/user_guide/geology/borehole_tres_sampling.png
   :alt: BH01 resistivity-versus-depth log with tres_column() samples marked, including one nan below max_depth.
   :width: 80%

   Each dot is one ``tres_column()`` result at the requested depth; the
   ``x`` marker at ``z=150`` sits below ``max_depth`` and comes back
   ``nan`` rather than a value read off the deepest interval.

Boreholes across a profile
------------------------------

A single log is only half the picture -- what makes a
:class:`~pycsamt.geology.Borehole` useful in pycsamt is having several
of them positioned along the same profile as independent ground-truth
control points. Adding a shallow log and a deeper one alongside BH01,
each with its own ``x``:

.. code-block:: pycon

   >>> bh03 = Borehole("BH03", x=150.0, intervals=[
   ...     Interval(top=0.0, bottom=6.0, lithology="topsoil", resistivity=300.0),
   ...     Interval(top=6.0, bottom=35.0, lithology="saprolite", resistivity=120.0),
   ...     Interval(top=35.0, bottom=60.0, lithology="fresh basement", resistivity=8000.0),
   ... ])
   >>> bh04 = Borehole("BH04", x=900.0, intervals=[
   ...     Interval(top=0.0, bottom=15.0, lithology="alluvium", resistivity=30.0),
   ...     Interval(top=15.0, bottom=45.0, lithology="aquifer sand", resistivity=90.0),
   ...     Interval(top=45.0, bottom=60.0, lithology="clay aquitard", resistivity=8.0),
   ...     Interval(top=60.0, bottom=90.0, lithology="bedrock", resistivity=4000.0),
   ... ])
   >>> [bh03.x, bh.x, bh04.x]
   [150.0, 500.0, 900.0]

There is nothing linking these three ``Borehole`` objects to each
other beyond sharing the same profile coordinate system -- no shared
container class is required to plot or reason about them together,
since ``x`` alone is enough to place each one correctly:

.. code-dropdown:: ../../../scripts/generate_user_guide_geology_borehole_figures.py
   :language: python
   :pyobject: make_borehole_profile_positions
   :linenos:
   :title: View the profile-position figure source code

.. figure:: /images/user_guide/geology/borehole_profile_positions.png
   :alt: Three boreholes drawn as vertical log strips at their real profile positions.
   :width: 100%

   BH03, BH01, and BH04 drawn at their real ``x`` positions along a
   shared profile axis -- the same spatial layout
   :class:`~pycsamt.interp.ModelCalibrator` sees when it blends
   several boreholes into one resistivity model in
   :doc:`../interpretation/workflow`. Depth increases downward in every
   strip; only the profile position changes between them.

Reading a LAS 2.0 log
------------------------

:meth:`~pycsamt.geology.Borehole.from_las` reads the industry-standard
LAS 2.0 well-log ASCII format directly, converting a continuous depth
curve into discrete intervals by grouping consecutive samples that
share the same lithology code. A minimal LAS 2.0 file has a well
section, a curve section naming each column by its mnemonic, and an
ASCII data section:

.. code-block:: text

   ~VERSION INFORMATION
   VERS.   2.0 : CWLS LOG ASCII STANDARD - VERSION 2.0
   WRAP.    NO : ONE LINE PER DEPTH STEP
   ~WELL INFORMATION
   WELL.        BH02 : WELL NAME
   NULL.    -9999.25 : NULL VALUE
   ~CURVE INFORMATION
   DEPT.M       : DEPTH
   RESD.OHMM    : DEEP RESISTIVITY
   LITH.        : LITHOLOGY CODE
   ~ASCII
   0.0     450.0   1
   2.0     430.0   1
   4.0     410.0   1
   6.0      60.0   2
   8.0      45.0   2
   10.0     40.0   2
   12.0    190.0   3
   14.0    180.0   3

Reading it groups the eight 2 m samples into three intervals, one per
run of a constant lithology code, with the interval's resistivity set
to the mean of the samples inside it:

.. code-block:: pycon

   >>> bh_las = Borehole.from_las("BH02.las", x=750.0)
   >>> bh_las
   Borehole('BH02', x=750.0 m, 3 intervals, depth=16.0 m)
   >>> for iv in bh_las.intervals:
   ...     print(iv)
   Interval(top=0.0, bottom=6.0, lithology='1', resistivity=430.0)
   Interval(top=6.0, bottom=12.0, lithology='2', resistivity=48.333333333333336)
   Interval(top=12.0, bottom=16.0, lithology='3', resistivity=185.0)

The default curve mnemonics are ``DEPT`` for depth and ``RESD`` for
resistivity; override them with ``depth_curve``/``resistivity_curve``
for a file that names them differently. ``lithology`` above is the
literal LAS code (``'1'``, ``'2'``, ``'3'``) rather than a rock name --
LAS 2.0 has no standard mapping from a numeric lithology code to a
name, so ``from_las`` does not invent one. Remap the codes to names
your project recognises before comparing this log's ``lithology``
field against a :class:`~pycsamt.geology.RockDatabase` classification;
comparing the ``resistivity`` values directly, as
:class:`~pycsamt.interp.ModelCalibrator` does, needs no such mapping.
Pass ``lithology_curve=None`` for a file with no lithology curve at
all -- every interval then gets the same generic label instead.

Plotting both logs side by side as a classic resistivity-versus-depth
track makes the difference between the two construction paths
concrete: BH01's four hand-specified intervals on the left, BH02's
eight raw 2 m samples collapsing into three ``from_las``-grouped
intervals on the right.

.. code-dropdown:: ../../../scripts/generate_user_guide_geology_borehole_figures.py
   :language: python
   :pyobject: make_borehole_log_comparison
   :linenos:
   :title: View the log-comparison figure source code

.. figure:: /images/user_guide/geology/borehole_log_comparison.png
   :alt: Resistivity-versus-depth log tracks for BH01 (direct construction) and BH02 (from_las).
   :width: 100%

   BH01 (left), colored by lithology from the intervals built directly
   above. BH02 (right), with the eight raw LAS samples as dots and the
   three grouped intervals as shaded bands -- the boundary between the
   blue and pink bands falls exactly where the lithology code changes
   from 1 to 2, between the samples at 4 m and 6 m.

Serializing a log
--------------------

``to_dataframe()`` (requires :mod:`pandas`) and ``to_dict()`` both
flatten a ``Borehole`` for export or storage; unlike most classes in
this package, ``Borehole`` writes both by hand rather than relying on
:meth:`~pycsamt.api.property.PyCSAMTObject.to_dict`, so that the
``thickness`` column below is always present without being stored
redundantly on every ``Interval``:

.. code-block:: pycon

   >>> bh_small = Borehole("BH01", x=500.0, intervals=intervals[:2])
   >>> bh_small.to_dataframe()
      top  bottom  thickness       lithology  resistivity
   0  0.0     8.0        8.0  lateritic soil        450.0
   1  8.0    31.0       23.0     clayey sand         42.0

Where to go next
-------------------

:doc:`../interpretation/workflow` covers loading a log from CSV and
using it as a :class:`~pycsamt.interp.ModelCalibrator` constraint end
to end. :doc:`rock_database` and :doc:`structural` cover the other two
data families in this package.
