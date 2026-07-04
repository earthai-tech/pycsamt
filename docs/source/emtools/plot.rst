.. _emtools_plot:

pycsamt.emtools.plot — Multi-Station Diagnostic Panels
=============================================================

:mod:`pycsamt.emtools.plot` is the "many stations, one figure" layer of
``emtools``: dense apparent-resistivity/phase panel grids for a whole
survey, full-tensor raw-data diagnostics, combined response/tipper
quality-control figures, side-by-side before/after comparisons, and
measured-vs-predicted fit grids with per-component RMS annotations.
Every function reads its axis labels, phase-wrapping convention, and
resistivity display scale from the shared
:data:`pycsamt.api.control.PYCSAMT_CONTROL`, so one context override
changes every figure's display consistently.

Functions
---------

- :func:`~pycsamt.emtools.plot.plot_sites_panels`
- :func:`~pycsamt.emtools.plot.plot_raw_sites_1d`
- :func:`~pycsamt.emtools.plot.plot_response_tipper`
- :func:`~pycsamt.emtools.plot.plot_sites_compare`
- :func:`~pycsamt.emtools.plot.plot_sites_fit_grid`

Full signatures and parameter documentation live in the
:doc:`API reference <../api/emtools>`.

Worked Example
--------------

Built on **L18PLT** (``data/AMT/WILLY_DATA/``) for the tensor-only
functions — reusing the same four cross-module stations
(``18-001A``, ``18-007U``, ``18-016A``, ``18-018A``) flagged in the
``anisotropy``/``impedance``/``lcurve`` examples — and on **KAP03**
(``data/MT/kap03lmt_edis``) for :func:`~pycsamt.emtools.plot.plot_response_tipper`,
the one function here that needs real tipper. It moves from the
simplest call (a station panel grid) through a raw, deliberately
styleless full-tensor diagnostic, a combined response/tipper
quality-control figure, a genuine before/after comparison using
:func:`~pycsamt.emtools.frequency.smooth_mavg` as an honest stand-in
for "processed" data, a measured-vs-predicted fit grid with
per-component RMS (explained carefully: the RMS is error-normalized,
not a raw log10(rho) difference, which matters once the "prediction"
is an un-fit smoothing rather than a genuine inversion response), and
finishes with the same display re-rendered under a different
:data:`~pycsamt.api.control.PYCSAMT_CONTROL` context to show every
label and limit tracking the active control.

.. include:: auto_examples/plot_plot.rst
