.. _emtools_impedance:

pycsamt.emtools.impedance — Impedance-Tensor Diagnostics
=============================================================

:mod:`pycsamt.emtools.impedance` works directly with the complex
impedance tensor rather than derived apparent resistivity: a polar
"phasor wheel" view of individual tensor components, a pseudo-section
of how far the off-diagonal terms depart from perfect antisymmetry
(:math:`Z_{xy} \approx -Z_{yx}`, the 1-D/2-D-friendly case), and a
determinant track with an error-propagated confidence band.

Functions
---------

- :func:`~pycsamt.emtools.impedance.plot_phasor_wheel`
- :func:`~pycsamt.emtools.impedance.plot_offdiag_antisym_residual`
- :func:`~pycsamt.emtools.impedance.plot_determinant_track`

Full signatures and parameter documentation live in the
:doc:`API reference <../api/emtools>`.

Worked Example
--------------

Built on **L18PLT** (``data/AMT/WILLY_DATA/``), the same CSAMT-band
line as the ``anisotropy``/``csumt``/``diag``/``fieldzone``/``gradient_imaging``
examples — which lets this page connect directly back to those pages'
own station rankings rather than starting from scratch. It moves from
the phasor wheel for one station (all off-diagonal, then split by
period band, then all four tensor components including the diagonal),
to the antisymmetry-residual pseudo-section and a full per-station
ranking, to a genuine cross-module check: does this residual actually
relate to the ``anisotropy`` example's Swift skew and :math:`\Lambda`
ratio rankings, or is it measuring something independent? (Answer:
strongly correlated with :math:`|\Lambda|` (r ≈ 0.72) and moderately
anti-correlated with skew (r ≈ -0.59) — the same two stations,
``18-016A`` and ``18-018A``, top both this ranking and the
``anisotropy`` ratio ranking, while the skew-ranking leader,
``18-007U``, sits at the *bottom* here.) It finishes with the
determinant track for that same contrasting pair, a look at how the
reported confidence level changes the band width, and a same-survey
comparison against neighbouring line L22PLT.

.. include:: auto_examples/plot_impedance.rst
