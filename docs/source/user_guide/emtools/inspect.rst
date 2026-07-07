.. _emtools_inspect:

pycsamt.emtools.inspect — First-Look Survey Inspection
===========================================================

:mod:`pycsamt.emtools.inspect` is the module to reach for right after
loading a survey: per-site summaries, missing-section checks, frequency
coverage, resistivity/phase and tipper curves, pseudo-sections, and a
single-station "everything at once" response view with an optional
model overlay.

Functions
---------

- :func:`~pycsamt.emtools.inspect.sites_summary`
- :func:`~pycsamt.emtools.inspect.list_missing_sections`
- :func:`~pycsamt.emtools.inspect.frequency_coverage`
- :func:`~pycsamt.emtools.inspect.plot_coverage`
- :func:`~pycsamt.emtools.inspect.plot_rhoa_phi`
- :func:`~pycsamt.emtools.inspect.plot_tipper_components`
- :func:`~pycsamt.emtools.inspect.pseudosection`
- :func:`~pycsamt.emtools.inspect.plot_station_response`

Full signatures and parameter documentation live in the
:doc:`API reference <../api/emtools>`.

Worked Example
--------------

Uses all three bundled datasets — **L18PLT** and **L22PLT**
(``data/AMT/WILLY_DATA/``) and **KAP03** (``data/MT/kap03lmt_edis``) —
since a good first look is precisely where their differences (tipper
presence, frequency band, per-station grid) matter most. It moves from
a per-site summary and missing-section check, through frequency
coverage, resistivity/phase and tipper curves (using small,
representative station subsets rather than the full survey, since an
all-station legend quickly becomes unreadable), a pseudo-section, and
finishes with the module's richest view: the full single-station
response for ``kap151`` — the same station singled out in the ``tf``
example for its unusually strong tipper — shown first alone and then
with a smoothed "model" overlay to demonstrate the built-in per-component
RMS misfit reporting.

.. include:: ../examples/emtools/plot_inspect.rst
   :start-after: .. _sphx_glr_examples_emtools_plot_inspect.py:
   :end-before: .. _sphx_glr_download_examples_emtools_plot_inspect.py:
