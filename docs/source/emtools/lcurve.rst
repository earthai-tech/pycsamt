.. _emtools_lcurve:

pycsamt.emtools.lcurve — L-Curve Regularization Selection
=============================================================

:mod:`pycsamt.emtools.lcurve` is the one ``emtools`` module that is not
about EDI data at all: it is a generic Tikhonov-regularization
diagnostic. Given any (misfit, roughness, :math:`\lambda`) triple —
from a 1-D sounding inversion, a 2-D/3-D model, or any other
regularized least-squares problem — it locates the "corner" of the
L-shaped trade-off curve between fitting the data and keeping the
model simple, and plots it.

Functions
---------

- :func:`~pycsamt.emtools.lcurve.lcurve_table`
- :func:`~pycsamt.emtools.lcurve.plot_lcurve`

Full signatures and parameter documentation live in the
:doc:`API reference <../api/emtools>`.

Worked Example
--------------

Because this module needs a real regularization *sweep* rather than
EDI values directly, the example builds one honestly instead of
fabricating a textbook curve: a small, genuine 1-D Tikhonov smoothing
problem — recover a smooth model from noisy real apparent resistivity
— solved in closed form at 60 :math:`\lambda` values for stations from
**L18PLT** (``data/AMT/WILLY_DATA/``), the same line used by the
``anisotropy``/``impedance`` examples. It starts with the mechanics on
a hand-built synthetic curve, then moves to the real sweep: what the
corner protects against (under- vs. over-regularized models plotted
against the data), a three-station comparison that reuses the stations
flagged elsewhere for unusual anisotropy/skew, a robustness check of
the two corner-picking methods (``curvature`` vs. ``maxdist``) across
smoothing levels, and finally a direction-of-travel view along the
curve using ``arrow_every``.

.. include:: auto_examples/plot_lcurve.rst
