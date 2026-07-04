.. _emtools_spectra:

pycsamt.emtools.spectra — Cross-Spectra Analysis and Visualization
=====================================================================

:mod:`pycsamt.emtools.spectra` works one level below the impedance
tensor used everywhere else in this gallery: it operates directly on
the raw cross-power spectral matrix stored in a
:class:`~pycsamt.seg.spectra.Spectra` object — one ``(n_chan, n_chan)``
complex matrix per frequency, read from an EDI file's
``>=SPECTRASECT`` block — rather than a pre-computed ``Z`` tensor.
That block is a richer, less common EDI structure than the
impedance-only files used elsewhere (``data/AMT/WILLY_DATA/``,
``data/MT/kap03lmt_edis/``), which only carry the final processed
estimates.

Functions
---------

Analysis:

- :func:`~pycsamt.emtools.spectra.coherence_matrix`
- :func:`~pycsamt.emtools.spectra.psd_table`
- :func:`~pycsamt.emtools.spectra.coherence_table`
- :func:`~pycsamt.emtools.spectra.snr_table`
- :func:`~pycsamt.emtools.spectra.band_select`
- :func:`~pycsamt.emtools.spectra.mask_low_coherence`
- :func:`~pycsamt.emtools.spectra.spectra_summary`

Visualization:

- :func:`~pycsamt.emtools.spectra.plot_psd`
- :func:`~pycsamt.emtools.spectra.plot_coherence`
- :func:`~pycsamt.emtools.spectra.plot_spectra_matrix`
- :func:`~pycsamt.emtools.spectra.plot_z_from_spectra`
- :func:`~pycsamt.emtools.spectra.plot_tipper_from_spectra`
- :func:`~pycsamt.emtools.spectra.plot_psd_section`
- :func:`~pycsamt.emtools.spectra.plot_coherence_section`

Full signatures and parameter documentation live in the
:doc:`API reference <../api/emtools>`.

Worked Example
--------------

Built on two de-identified real-field spectra EDI files bundled under
``data/MT/SPECTRA/`` (see ``data/MT/SPECTRA/README.md`` for the
anonymization notes) — one short-period/AMT-band station and one
broadband/long-period station, chosen precisely because they carry the
raw ``>=SPECTRASECT`` block this module needs rather than only a
processed ``Z`` tensor. It moves from a single channel's power
spectrum, through pairwise coherence between the electric and magnetic
channels, tabular SNR and per-frequency summaries, a practical
band-selection recipe that turns "9 of 51 frequencies fail a coherence
threshold" into "restrict the band and all 27 remaining frequencies
pass", the full 7x7 cross-power matrix at a single frequency spanning
nearly ten orders of magnitude, apparent resistivity/phase and
induction tipper recovered directly from the spectra (confirming the
same structural complexity seen via the impedance tensor elsewhere in
this gallery), and finishes with the two functions built for multiple
stations at once — PSD and coherence pseudo-sections that
automatically restrict themselves to the overlapping frequency band
between stations with very different coverage.

.. include:: auto_examples/plot_spectra.rst
