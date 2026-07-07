.. _emtools_source_effects:

pycsamt.emtools.source_effects — CSAMT Source Overprint and Near-Field Effects
=================================================================================

:mod:`pycsamt.emtools.source_effects` quantifies how the artificial
CSAMT transmitter contaminates a sounding, using two independent
formulas from two different papers: the Yan & Fu (2004) / Da et al.
(2016) ground-wave/surface-wave amplitude ratio :math:`\beta_{Ey}` (a
source-overprint index, threshold 3%), and the Wang & Lin (2023)
skin-depth field-zone classification with a near-field correction
built on the equatorial HED transfer function :math:`F(p)`.

Functions
---------

- :func:`~pycsamt.emtools.source_effects.overprint_beta`
- :func:`~pycsamt.emtools.source_effects.detect_source_overprint`
- :func:`~pycsamt.emtools.source_effects.source_overprint_table`
- :func:`~pycsamt.emtools.source_effects.plot_overprint_section`
- :func:`~pycsamt.emtools.source_effects.normalize_response`
- :func:`~pycsamt.emtools.source_effects.correct_near_field`
- :func:`~pycsamt.emtools.source_effects.plot_normalized_response`

Full signatures and parameter documentation live in the
:doc:`API reference <../api/emtools>`.

Worked Example
--------------

Both formulas need the source-receiver offset :math:`r`, which no
real EDI file carries — the same gap already handled honestly in the
``fieldzone`` example, whose representative 2 km offset (and 500 m /
8 km sensitivity pair) is reused here for **L18PLT**
(``data/AMT/WILLY_DATA/``) so the two pages stay consistent. It moves
from the pure-math :math:`\beta_{Ey}` formula's offset dependence,
through per-frequency and per-station overprint detection (including
the da2016 low-/high-frequency slope criterion, and why the module's
default frequency split leaves every low-frequency slope undefined on
this particular line), the overprint pseudo-section, the Wang & Lin
normalized-response/field-zone view and its two-panel pseudo-section,
a near-field correction whose multi-decade swing turns out to be
genuine physics rather than a bug, and finishes by merging the two
independent formulas point-by-point — finding they agree completely
on which frequencies are contaminated, despite coming from unrelated
papers and physical arguments.

.. include:: ../examples/emtools/plot_source_effects.rst
   :start-after: .. _sphx_glr_examples_emtools_plot_source_effects.py:
   :end-before: .. _sphx_glr_download_examples_emtools_plot_source_effects.py:
