.. _emtools_anisotropy:

pycsamt.emtools.anisotropy — CSAMT Axial-Anisotropy Diagnostics
====================================================================

:mod:`pycsamt.emtools.anisotropy` implements the axial-anisotropy metrics
from Wang & Tan (2017, *J. Appl. Geophys.* 146, 27-36): the two
independent Cagniard apparent resistivities :math:`\rho_{xy}` and
:math:`\rho_{yx}` computed from the off-diagonal impedance tensor, their
log-ratio :math:`\Lambda = \log_{10}(\rho_{xy}/\rho_{yx})` (zero for a
perfectly isotropic 1-D earth), and the classic Swift (1967) skew and
strike angle from the full tensor. All functions accept the same
flexible ``sites`` input as the rest of ``emtools`` — a path, an
``EDIFile``/``EDICollection``, an ``APISurvey``, or an iterable of
site-like objects.

Because the method relies on the off-axis Cagniard resistivities being
well-resolved across a broad, closely sampled frequency sweep, it is
best suited to CSAMT/AMT-style data (dense coverage from a few Hz up to
several kHz) rather than sparsely sampled long-period natural-source MT.

Functions
---------

- :func:`~pycsamt.emtools.anisotropy.analyze_anisotropy`
- :func:`~pycsamt.emtools.anisotropy.anisotropy_table`
- :func:`~pycsamt.emtools.anisotropy.plot_anisotropy`

Full signatures and parameter documentation live in the
:doc:`API reference <../api/emtools>`.

Worked Example
--------------

The example below uses **L18PLT** and **L22PLT**, two real AMT/CSAMT-band
lines (1 Hz-10.4 kHz, 53 frequencies) bundled at
``data/AMT/WILLY_DATA/`` (see its ``README.md``) — a full, energetic
impedance tensor is exactly what this method needs, unlike KAP03's
long-period natural-source band used in the ``tf`` example.

It moves from a single station's raw ratio/skew-vs-period curves, to a
per-station ranking, to the module's core station x period
pseudo-section shown for three of its four selectable metrics, and
finishes with two genuinely data-driven findings rather than just
repeating the same plot: the ratio and skew rankings disagree at the
extremes (``18-016A`` vs ``18-007U``), and pooling both lines' 53
stations turns that into a quantified, moderate *negative* correlation
(Pearson r ≈ -0.5) between the two indicators — evidence that a real
anisotropy assessment needs both metrics, since neither alone tells the
full story. Along the way, one station's Swift-skew pseudo-section shows
a textbook numerical instability (a single-frequency spike above 40, from
the metric's denominator passing near zero) worth knowing to recognise
rather than mistake for a real anomaly.

.. include:: auto_examples/plot_anisotropy.rst
