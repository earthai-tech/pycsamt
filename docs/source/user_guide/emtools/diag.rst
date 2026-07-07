.. _emtools_diag:

pycsamt.emtools.diag — Polar Uncertainty Diagnostics
========================================================

:mod:`pycsamt.emtools.diag` adapts the "k-diagram" polar-uncertainty
framework (Kouadio 2025, *JOSS* 10(116)) to CSAMT apparent resistivity.
Given an observed sounding and a **predicted quantile interval**
[L, U] at each frequency, it checks empirical coverage, how interval
width drifts across the frequency band (a proxy for how uncertainty
grows with depth), and the rose of relative residuals against a point
prediction. Unlike the other ``emtools`` modules covered so far, this
one needs a *prediction* to evaluate against, not just the observed
data.

.. warning::

   While building this example, ``rho_obs`` in this module was found to
   be off by a ~10^5-10^6 factor from the physically realistic values
   :mod:`pycsamt.emtools.csumt` computes for the same data — a missing
   practical-units conversion factor in ``_rho_a_from_z``. It has been
   fixed to match ``csumt``'s convention (``0.2|Z|²/f``). A codebase-wide
   check found the same pattern in two more places, now fixed the same
   way:

   - :mod:`pycsamt.emtools.anisotropy`'s ``_rho_and_phase`` — its own
     ``ratio_log10`` results were unaffected since the missing factor
     cancels in a ratio, but the absolute ``rho_xy_ohmm``/``rho_yx_ohmm``
     columns were wrong before this fix.
   - :mod:`pycsamt.emtools.remove_noise`'s ``correct_static_shift`` —
     likewise unaffected in practice (its correction factor is a
     same-frequency ``sqrt(rho_smooth/rho_obs)`` ratio), fixed for
     consistency.

   Other ``_MU0``/skin-depth-style formulas found in ``advanced.py``,
   ``fieldzone.py``, ``source_array.py``, ``source_effects.py``, and the
   AI-inversion PINN operators (``pycsamt.ai.inversion._pinn_ops_*``)
   were checked and are unaffected: they either consume an
   already-correct precomputed ``rho`` (e.g. ``Site.rho``) or, in the
   PINN case, forward-model Z from resistivity in a self-contained SI
   system where the matching SI apparent-resistivity formula is correct.

Functions
---------

- :func:`~pycsamt.emtools.diag.coverage_score`
- :func:`~pycsamt.emtools.diag.rho_coverage`
- :func:`~pycsamt.emtools.diag.rho_error_stats`
- :func:`~pycsamt.emtools.diag.coverage_table`
- :func:`~pycsamt.emtools.diag.plot_polar_coverage`
- :func:`~pycsamt.emtools.diag.plot_polar_errors`
- :func:`~pycsamt.emtools.diag.plot_width_drift`

Full signatures and parameter documentation live in the
:doc:`API reference <../api/emtools>`.

Worked Example
--------------

There is no real forecasting model bundled with this documentation, so
the example builds one honestly: a smoothed (rolling-median) version of
each station's own real observed resistivity from **L18PLT**
(``data/AMT/WILLY_DATA/``) stands in for a "model", with synthetic
quantile bounds around it whose width is designed to grow with period.
It starts with the pure ``coverage_score`` concept on a tiny toy array,
moves to one real station's observed-vs-predicted curve, a per-station
coverage ranking, the module's polar coverage view, a width-drift check
that confirms the intentionally-designed growing-uncertainty trend, and
an error rose against the smoothed model. It finishes with the point of
the whole module: three variants of the same bounds — sensible,
artificially narrowed ("overconfident"), and artificially widened
("underconfident") — show the diagnostics correctly catching a bad
calibration (0/28 stations pass when narrowed) and the coverage/width
trade-off of an overly cautious one (nearly 3x wider intervals for
~8 points of extra coverage).

.. include:: ../examples/emtools/plot_diag.rst
   :start-after: .. _sphx_glr_examples_emtools_plot_diag.py:
   :end-before: .. _sphx_glr_download_examples_emtools_plot_diag.py:
