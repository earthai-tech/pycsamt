.. _emtools_fieldzone:

pycsamt.emtools.fieldzone — CSAMT Field-Zone Classification
================================================================

:mod:`pycsamt.emtools.fieldzone` answers a question specific to
controlled-source AMT: at a given frequency and source-receiver offset
*r*, is a measurement in the plane-wave (far) zone the standard
apparent-resistivity formula assumes, or has it drifted into the
near/transition zone where that assumption breaks down? Both
:func:`~pycsamt.emtools.fieldzone.classify_field_zones` and
:func:`~pycsamt.emtools.fieldzone.near_field_factor` reduce to the same
dimensionless parameter, :math:`|k \cdot r| = r / \delta_B` — the
source-receiver distance measured in Bostick skin depths (Chen & Yan,
2005).

Functions
---------

- :func:`~pycsamt.emtools.fieldzone.classify_field_zones`
- :func:`~pycsamt.emtools.fieldzone.near_field_factor`
- :func:`~pycsamt.emtools.fieldzone.plot_field_zones`

Full signatures and parameter documentation live in the
:doc:`API reference <../api/emtools>`.

Worked Example
--------------

This is a genuinely CSAMT-specific concept — it needs a real
source-receiver offset, which natural-source **L18PLT**
(``data/AMT/WILLY_DATA/``) does not record. As in the ``csumt``
example, this page applies the same physics to L18PLT's real,
CSAMT-band (1 Hz-10.4 kHz) apparent resistivity using a few
representative assumed offsets (500 m, 2 km, 8 km) — a legitimate way
to see how the method behaves, while being upfront that the exact
offset numbers are chosen for illustration, not read from survey
metadata.

It opens with the pure |k.r| relationship (no data), then one real
station's measured |k.r| curve, a direct cross-check between the
threshold-based zoning and the continuous near-field correction factor
(they agree closely — |F| averages 0.99 in the far zone and explodes to
a mean above 900, peaking near 17,000, in the near zone), the module's
pseudo-section, and finishes with two advanced points: how much the
*assumed* offset alone changes the zone mix for the exact same data
(the near-field fraction runs from 35% down to 0.3% across the three
offsets tested), and an illustrative near-field-corrected sounding
curve showing where an uncorrected plane-wave inversion would be
misled by this station's own longest-period data.

.. include:: auto_examples/plot_fieldzone.rst
