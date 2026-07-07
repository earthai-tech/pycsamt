.. _emtools_skew:

pycsamt.emtools.skew — Phase-Tensor and Bahr Skew Diagnostics
=================================================================

:mod:`pycsamt.emtools.skew` measures how far a station departs from a
simple 1-D/2-D structure using two independent formulas — the classic
Bahr (1988) invariant :math:`\eta` computed directly from :math:`Z`,
and the Caldwell-Bibby-Bahr phase-tensor skew angle :math:`\beta` (via
:mod:`pycsamt.emtools.tensor`) — then offers several ways to mask,
trim, or vote on frequencies where that skew is acceptably low.

Functions
---------

- :func:`~pycsamt.emtools.skew.bahr_skewness`
- :func:`~pycsamt.emtools.skew.skew_table`
- :func:`~pycsamt.emtools.skew.mask_by_skew`
- :func:`~pycsamt.emtools.skew.keep_longest_low_skew`
- :func:`~pycsamt.emtools.skew.close_skew_gaps`
- :func:`~pycsamt.emtools.skew.select_low_skew_band`
- :func:`~pycsamt.emtools.skew.plot_skew_traffic_psection`
- :func:`~pycsamt.emtools.skew.plot_skew_percentile_ribbon`
- :func:`~pycsamt.emtools.skew.plot_skew_vote_band`
- :func:`~pycsamt.emtools.skew.plot_skewness`

Full signatures and parameter documentation live in the
:doc:`API reference <../api/emtools>`.

Worked Example
--------------

Built on **L18PLT** (``data/AMT/WILLY_DATA/``), the same line the
``qc`` example already flagged ``high_skew`` at every one of its 28
stations under a strict default threshold. That finding is traced to
its source here (median :math:`|\beta|` 41 degrees, mean 44.5, survey
-wide) and through its consequences in every skew-based masking
function in the module: near-total masking at textbook-strict
defaults, a "kept everything" result from
:func:`~pycsamt.emtools.skew.keep_longest_low_skew` that turns out to
mean its fallback silently triggered rather than that a station
passed, and a survey-wide vote that keeps *zero* rows at a 60%
threshold because each station's acceptable band sits at a different
frequency. Two stations recur throughout — ``18-016A`` and
``18-018A``, the calmest and most phase-tensor-skewed of the cast from
earlier examples — including a section showing that Bahr :math:`\eta`
and phase-tensor :math:`\beta` do not always rank stations the same
way, and a closing section reconciling why a vote-band plot and a
band-selection function can report very different numbers from what
looks like the same threshold.

.. include:: ../examples/emtools/plot_skew.rst
   :start-after: .. _sphx_glr_examples_emtools_plot_skew.py:
   :end-before: .. _sphx_glr_download_examples_emtools_plot_skew.py:
