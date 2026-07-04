.. _emtools_csumt:

pycsamt.emtools.csumt — Bostick Depth and CSUMT Survey Design
==================================================================

:mod:`pycsamt.emtools.csumt` is built around the Bostick depth transform
(Zhang et al. 2025, *Measurement*):

.. math::

    D(f) = 356 \sqrt{\rho_a(f) / f} \quad [\text{m}]

The module has two distinct halves. The **pure survey-planning
functions** — :func:`~pycsamt.emtools.csumt.bostick_depth_from_rho`,
:func:`~pycsamt.emtools.csumt.frequency_for_depth`,
:func:`~pycsamt.emtools.csumt.frequency_schedule`,
:func:`~pycsamt.emtools.csumt.vertical_resolution_pair` — need no EDI
data, only a resistivity estimate, and are meant for designing a CSUMT
acquisition before going to the field. The **sites-based analysis**
functions — :func:`~pycsamt.emtools.csumt.bostick_depth`,
:func:`~pycsamt.emtools.csumt.vertical_resolution`,
:func:`~pycsamt.emtools.csumt.depth_coverage_table`,
:func:`~pycsamt.emtools.csumt.plot_depth_section` — apply the same
transform to measured apparent resistivity from real EDI data.

Functions
---------

- :func:`~pycsamt.emtools.csumt.bostick_depth_from_rho`
- :func:`~pycsamt.emtools.csumt.frequency_for_depth`
- :func:`~pycsamt.emtools.csumt.frequency_schedule`
- :func:`~pycsamt.emtools.csumt.vertical_resolution_pair`
- :func:`~pycsamt.emtools.csumt.bostick_depth`
- :func:`~pycsamt.emtools.csumt.vertical_resolution`
- :func:`~pycsamt.emtools.csumt.depth_coverage_table`
- :func:`~pycsamt.emtools.csumt.plot_depth_section`

Full signatures and parameter documentation live in the
:doc:`API reference <../api/emtools>`.

Worked Example
--------------

The Bostick transform itself does not care what frequency band it is
given, so this example applies it to **L18PLT**/**L22PLT**
(``data/AMT/WILLY_DATA/``, 1 Hz-10.4 kHz) even though CSUMT proper
targets a much higher band (9.6 kHz-614.4 kHz) for very shallow,
high-resolution work — that contrast is itself part of what the example
shows.

It starts with the pure, data-free planning API: the Bostick
relationship across resistivities, the depth range the CSUMT band can
actually reach (only ~14-115 m even at 1000 Ω·m), and a frequency
schedule built from target depths — including the real trap that
:func:`~pycsamt.emtools.csumt.frequency_schedule` silently drops any
target outside the instrument's band rather than warning about it. It
then switches to the real-data half — one station's depth-vs-period
curve, a per-station depth-coverage ranking, the module's pseudo-section
— and finishes with two advanced checks: whether measured vertical
resolution really coarsens with depth the way the formula predicts
(binned across all 28 stations, compared against the purely analytic
curve for the same line's median resistivity), and a same-survey,
two-line comparison against L22PLT.

.. include:: auto_examples/plot_csumt.rst
