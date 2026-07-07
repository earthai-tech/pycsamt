.. _emtools_tensor:

pycsamt.emtools.tensor — Phase-Tensor Analysis and Impedance-Tensor Editing
==============================================================================

:mod:`pycsamt.emtools.tensor` is the largest ``emtools`` module: the
phase-tensor invariants (Caldwell et al. 2004) computed from the
impedance tensor — principal axes, strike angle θ, skew β, and
ellipticity — a full family of pseudo-section, rose-diagram,
geographic-map, and per-station "ellipse strip" views built on them,
and a set of impedance-tensor editing operations (rotation,
antisymmetrisation, inversion, sensor-orientation correction,
outlier-clipping, off-diagonal balancing).

Functions
---------

Tensor editing:

- :func:`~pycsamt.emtools.tensor.rotate_to_strike`
- :func:`~pycsamt.emtools.tensor.rotate`
- :func:`~pycsamt.emtools.tensor.rotate_by_map`
- :func:`~pycsamt.emtools.tensor.antisymmetrize`
- :func:`~pycsamt.emtools.tensor.invert`
- :func:`~pycsamt.emtools.tensor.orient_from_sensors`
- :func:`~pycsamt.emtools.tensor.sigma_clip_z`
- :func:`~pycsamt.emtools.tensor.balance_offdiag`

Phase-tensor computation and visualisation:

- :func:`~pycsamt.emtools.tensor.build_phase_tensor_table`
- :func:`~pycsamt.emtools.tensor.plot_phase_tensor_psection`
- :func:`~pycsamt.emtools.tensor.plot_phase_tensor_skewmap`
- :func:`~pycsamt.emtools.tensor.plot_theta_vs_period`
- :func:`~pycsamt.emtools.tensor.plot_ellipticity_psection`
- :func:`~pycsamt.emtools.tensor.plot_dimensionality_psection`
- :func:`~pycsamt.emtools.tensor.plot_phase_tensor_rose`
- :func:`~pycsamt.emtools.tensor.plot_phase_tensor_map`
- :func:`~pycsamt.emtools.tensor.plot_phase_tensor_summary`
- :func:`~pycsamt.emtools.tensor.phase_tensor_legend`
- :func:`~pycsamt.emtools.tensor.plot_dimensionality_grid`
- :func:`~pycsamt.emtools.tensor.plot_theta_stability_stripe`
- :func:`~pycsamt.emtools.tensor.plot_skew_ellipt_density`
- :func:`~pycsamt.emtools.tensor.plot_theta_rose_grid`
- :func:`~pycsamt.emtools.tensor.plot_phase_tensor_strip`
- :func:`~pycsamt.emtools.tensor.plot_phase_tensor_strip_grid`

Full signatures and parameter documentation live in the
:doc:`API reference <../api/emtools>`.

Worked Example
--------------

Built on **L18PLT** (``data/AMT/WILLY_DATA/``), with **KAP03**
(``data/MT/kap03lmt_edis/``) brought in for the geographic map's
tipper overlay. Five real bugs turned up and were fixed along the way:
three impedance-rotation functions
(:func:`~pycsamt.emtools.tensor.rotate`,
:func:`~pycsamt.emtools.tensor.rotate_to_strike`,
:func:`~pycsamt.emtools.tensor.rotate_by_map`) were complete silent
no-ops — one handed a whole ``Sites`` collection to a helper built for
a single EDI item, the other two called that helper on a
freshly-wrapped ``Sites`` object instead of the underlying EDI item,
and ``rotate_to_strike`` compounded this with a strike-estimation call
that always returned an empty table;
:func:`~pycsamt.emtools.tensor.orient_from_sensors` raised a
``TypeError`` on every single call (wrong keyword names passed to its
own math helper); and :func:`~pycsamt.emtools.tensor.plot_phase_tensor_map`
crashed on any survey whose EDI files carry no per-station coordinates
in the standard header field — which is exactly KAP03's situation —
instead of showing the graceful message already written for that case.
None of the five had test coverage.

The example moves from the phase-tensor invariants that everything
else in the module is built on, through simple per-station and
pseudo-section views, the flagship ellipse pseudo-section (with a
dimensionality finding — 97.9% of cells classify as 3-D under the
module's own default thresholds, a third independent confirmation of
the strong 3-D distortion already documented via skew in :doc:`qc` and
via static shift in :doc:`ss`), rose diagrams that reproduce a number
already seen in :doc:`strike`, a stability stripe, joint
skew-ellipticity density, a combined summary figure, the geographic
map (exposing and working around the bug above), per-station ellipse
strips, and finishes with the impedance-editing operations —
antisymmetrisation, inversion, sensor-orientation correction (exposing
another bug above), sigma-clipping, off-diagonal balancing, and
rotation (exposing the three rotation bugs above) — each verified with
concrete before/after numbers on real data.

.. include:: ../examples/emtools/plot_tensor.rst
   :start-after: .. _sphx_glr_examples_emtools_plot_tensor.py:
   :end-before: .. _sphx_glr_download_examples_emtools_plot_tensor.py:
