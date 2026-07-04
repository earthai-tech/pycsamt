.. _emtools_strike:

pycsamt.emtools.strike — Geoelectric Strike Estimation and Visualization
===========================================================================

:mod:`pycsamt.emtools.strike` estimates the geoelectric strike
direction three independent ways — an impedance-tensor rotation
sweep, phase-tensor azimuth, and a weighted consensus blend — applies
that estimate to rotate stations onto strike, and renders it through
five plot styles: a per-frequency ribbon, an along-line profile,
geographic map-sticks, single and multi-line rose diagrams, and a
three-panel diagnostic comparable to MTPy's ``StrikeAnalysis`` plot.

Functions
---------

- :func:`~pycsamt.emtools.strike.estimate_strike_sweep`
- :func:`~pycsamt.emtools.strike.estimate_strike_phase_tensor`
- :func:`~pycsamt.emtools.strike.estimate_strike_consensus`
- :func:`~pycsamt.emtools.strike.rotate_to_strike`
- :func:`~pycsamt.emtools.strike.strike_curve_sweep`
- :func:`~pycsamt.emtools.strike.plot_strike_rose_by_line`
- :func:`~pycsamt.emtools.strike.plot_strike_rose`
- :func:`~pycsamt.emtools.strike.plot_strike_ribbon`
- :func:`~pycsamt.emtools.strike.plot_strike_profile`
- :func:`~pycsamt.emtools.strike.plot_strike_mapsticks`
- :func:`~pycsamt.emtools.strike.plot_strike_analysis`

Full signatures and parameter documentation live in the
:doc:`API reference <../api/emtools>`.

Worked Example
--------------

Built on **L18PLT** (``data/AMT/WILLY_DATA/``), with its sibling line
**L22PLT** brought in for the multi-line rose comparison. Two real
bugs turned up and were fixed along the way: :func:`~pycsamt.emtools.strike.rotate_to_strike`
located each station's angle correctly but then rotated the wrong
object — a freshly-wrapped ``Sites`` collection instead of the
underlying EDI item — so it could never find the section it needed to
rotate and silently returned data identical to its input; and
:func:`~pycsamt.emtools.strike.plot_strike_profile`'s
``sort_by="lon"``/``"lat"``/``"auto"`` checked only flat ``.lon``/``.lat``
attributes that real ``Site`` objects do not carry (the same recurring
pattern already found in :doc:`ss` and :doc:`dimensionality`), so
station order silently collapsed to alphabetical-by-name. The example
moves from the three strike estimators and an investigation into why
naively correlating their raw angles is misleading for axial data,
through rotating data onto strike (with the before/after numbers that
expose the bug above), a per-frequency ribbon view, single- and
multi-line rose diagrams — including an honest look at why this
survey's station-naming convention defeats the module's automatic
line-grouping heuristic — a frequency-band decomposition, a
geographic map-sticks view, an along-line profile (exposing the second
bug above), and finishes with the combined three-panel Strike/PT/Tipper
diagnostic.

.. include:: auto_examples/plot_strike.rst
