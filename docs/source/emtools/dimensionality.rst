.. _emtools_dimensionality:

pycsamt.emtools.dimensionality — Rule-Based and Dictionary-Learned Dimensionality
=======================================================================================

:mod:`pycsamt.emtools.dimensionality` answers "how far from 1-D/2-D is
this sounding?" two independent ways. The classical route uses the
phase tensor's Bibby (2005) skew angle :math:`\beta` and ellipticity
against fixed thresholds. The second route standardizes four
impedance/phase-tensor features and learns a small sparse dictionary
(MOD + ISTA) over them with no thresholds at all, then labels each
station-frequency by its dominant learned atom. On top of both sit
masking helpers (drop frequencies judged too 3-D) and a strike-rotation
+ antisymmetrization step meant to prepare data for 2-D inversion.

.. warning::

   Building this example surfaced **four** independent, previously
   untested bugs, all now fixed:

   1. :func:`~pycsamt.emtools.dimensionality.plot_dim_map` looked up
      coordinates via ``ed.lat``/``.latitude``/``.lon``/``.longitude``
      — attributes the real ``Site`` class does not have (it exposes
      one ``.coords`` property returning ``(lat, lon, elev)``) — so it
      always reported "no coords", even for lines with valid ones.
   2. & 3. :func:`~pycsamt.emtools.dimensionality.project_to_2d` was
      broken in **both** of its modes: the auto-strike path called a
      ``pycsamt.site.edit.rotate_to_strike`` that does not exist
      (it lives in :mod:`pycsamt.emtools.tensor`), and the explicit-
      strike path passed a keyword (``angle=``) that
      ``pycsamt.site.edit.rotate`` does not accept (it is
      ``angle_deg``). Both now reuse ``tensor``'s already-correct
      ``rotate``/``rotate_to_strike`` wrappers.
   4. :func:`~pycsamt.emtools.dimensionality.learn_dim_dictionary`'s
      training loop transposed its dictionary update a second time
      after an internal helper had already returned it correctly
      oriented, flipping its shape and crashing on the second
      iteration with the default ``n_iter=40`` — every real training
      run.

   No test in the repository exercised any of these four code paths.

Functions
---------

- :func:`~pycsamt.emtools.dimensionality.phase_features_table`
- :func:`~pycsamt.emtools.dimensionality.classify_dimensionality`
- :func:`~pycsamt.emtools.dimensionality.mask_by_dimensionality`
- :func:`~pycsamt.emtools.dimensionality.project_to_2d`
- :func:`~pycsamt.emtools.dimensionality.learn_dim_dictionary`
- :func:`~pycsamt.emtools.dimensionality.encode_dimensionality`
- :func:`~pycsamt.emtools.dimensionality.mask_by_dictionary`
- :func:`~pycsamt.emtools.dimensionality.plot_dim_confidence_grid`
- :func:`~pycsamt.emtools.dimensionality.plot_dim_occupancy_area`
- :func:`~pycsamt.emtools.dimensionality.plot_dim_map`
- :func:`~pycsamt.emtools.dimensionality.plot_atom_psection`

Full signatures and parameter documentation live in the
:doc:`API reference <../api/emtools>`.

Worked Example
--------------

The example uses **L18PLT** (``data/AMT/WILLY_DATA/``), the same
CSAMT-band line as the ``anisotropy``/``csumt``/``diag`` examples. It
starts with one station's raw skew/ellipticity curves against the
classifier's own default thresholds, then shows the rule's geometry
directly by scattering every (station, frequency) pair in feature
space — at defaults, only 31 of 1484 pairs qualify as 1-D/2-D. Rather
than stop there, a threshold-sensitivity sweep asks whether that is a
threshold problem or a data problem: the 3-D fraction falls smoothly
but stays above 70% even at a skew threshold nearly 10x the default,
which is a property of this survey, not a knob to tune away. After the
module's own pseudo-section, occupancy, and map views (the last only
working thanks to the ``.coords`` fix), two advanced checks close the
page: whether ``project_to_2d``'s strike rotation actually reduces the
classified 3-D fraction (mean skew/ellipticity both improve slightly,
but the 3-D fraction does not — it barely moves), and whether the
rule-based and dictionary-learned labels agree (68% of the time —
meaningful, but far from total, and the dictionary approach masks a
much smaller fraction of the survey when used for the same
keep-1D/2D filtering).

.. include:: auto_examples/plot_dimensionality.rst
