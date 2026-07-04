.. _emtools_source_array:

pycsamt.emtools.source_array — Phased-Array CSAMT Transmitter Design
=========================================================================

:mod:`pycsamt.emtools.source_array` is the one ``emtools`` module with
no site or EDI data at all: it implements the phased-array
transmitting source (PAS) antenna theory of Fan, Zhang & Wang (2022) —
element pattern, array factor, beam steering, directivity, and SNR
gain for an N-element linear array of CSAMT dipole sources, replacing
the traditional single-dipole antenna source (SDAS).

Functions
---------

- :func:`~pycsamt.emtools.source_array.wavenumber`
- :func:`~pycsamt.emtools.source_array.sdas_element_pattern`
- :func:`~pycsamt.emtools.source_array.array_factor`
- :func:`~pycsamt.emtools.source_array.pas_pattern`
- :func:`~pycsamt.emtools.source_array.beam_steer`
- :func:`~pycsamt.emtools.source_array.steering_angles`
- :func:`~pycsamt.emtools.source_array.sdas_directivity`
- :func:`~pycsamt.emtools.source_array.snr_gain_db`
- :func:`~pycsamt.emtools.source_array.plot_radiation_pattern`

Full signatures and parameter documentation live in the
:doc:`API reference <../api/emtools>`.

Worked Example
--------------

Every number in this example comes directly from the formulas
themselves at representative CSAMT frequencies and resistivities — no
bundled survey is needed or used. It moves from the single-dipole
element pattern and the earth-vs-free-space wavenumber contrast,
through the array factor at a low CSAMT frequency (where a realistic
2 km element spacing turns out to be far too small, in wavelengths,
for classical nulls to form) and again at a higher one (where the same
physical array develops real nulls, side lobes, and a genuine grating
lobe), beam steering and its exact numerical confirmation, the
combined element-times-array pattern, directivity's non-monotonic
dependence on dipole length versus SNR gain's clean :math:`20\log_{10}N`
law, and finishes with a concrete 8-element design whose 18 dB SNR
gain turns out to be split almost evenly between the intended beam and
an equally strong grating lobe — a real caveat caught only by checking
the numbers rather than trusting the first, more flattering read of
the figure.

.. include:: auto_examples/plot_source_array.rst
