.. _emtools_tf:

pycsamt.emtools.tf — Induction Vectors and Tipper Diagnostics
================================================================

:mod:`pycsamt.emtools.tf` turns the vertical-field (tipper) component of a
transfer function into induction-vector diagnostics: arrow maps and
profiles, polar and rose plots of arrow azimuth, and side-by-side panels
comparing the Parkinson and Wiese sign conventions. All functions accept
the same flexible ``sites`` input as the rest of ``emtools`` — a path, an
``EDIFile``/``EDICollection``, or an iterable of site-like objects — via
:func:`~pycsamt.emtools._core.ensure_sites`.

Induction arrows point toward laterally conductive structure (in the
Parkinson convention) and their length grows with how strongly the tipper
responds at a given period. Plotting them across a profile or a rose
diagram is one of the fastest qualitative checks for regional strike and
localized conductors before running a full 2-D/3-D inversion.

Functions
---------

- :func:`~pycsamt.emtools.tf.plot_tipper_hodograms`
- :func:`~pycsamt.emtools.tf.plot_induction_arrows`
- :func:`~pycsamt.emtools.tf.plot_induction_map`
- :func:`~pycsamt.emtools.tf.plot_induction_section`
- :func:`~pycsamt.emtools.tf.plot_induction_convention`
- :func:`~pycsamt.emtools.tf.plot_tipper_polar`
- :func:`~pycsamt.emtools.tf.plot_induction_rose`
- :func:`~pycsamt.emtools.tf.plot_induction_multiperiod_map`
- :func:`~pycsamt.emtools.tf.plot_induction_map_from_spectra`
- :func:`~pycsamt.emtools.tf.plot_tipper_polar_from_spectra`
- :func:`~pycsamt.emtools.tf.plot_induction_rose_from_spectra`

Full signatures and parameter documentation live in the
:doc:`API reference <../api/emtools>`.

Worked Example
--------------

The example below uses **KAP03**, a real long-period MT profile from
SAMTEX bundled at ``data/MT/kap03lmt_edis`` (see ``data/MT/README.md`` for
origin and required attribution) — unlike the AMT lines in
``data/AMT/WILLY_DATA/``, it was recorded with a real vertical-field
sensor, so every arrow below comes from genuine tipper data. The code is
executed when these docs are built, so the figures always match the
current ``pycsamt.emtools.tf`` code.

It works through the module from simple to composite — one station's raw
tipper (hodograms, polar plot), then single- and multi-period arrow
maps, sign conventions, azimuth roses, a period pseudo-section, and
finally the publication-style multi-panel map — and each figure is
followed by a short analysis of what it actually shows for this survey,
not just a description of the plot type. The two headline findings that
recur across several of the views: station ``kap151`` has an
exceptionally strong but band-limited response (roughly 50-1600 s, see
the hodogram, polar plot, and pseudo-section), while the profile's SW
end (``kap103``-``kap121``) carries a weaker but far more *directionally
consistent* response that sharpens at the longest recorded periods (see
the rose comparison and the multi-period map) — a good candidate for a
real, laterally continuous deep conductor.

For a future module page that needs site-like data but has no suitable
bundled real dataset, ``docs/examples/emtools/_synthetic.py``
provides a reusable synthetic-tipper generator (real or made-up station
geometry, a single-conductor response) — see its docstring.

.. include:: ../examples/emtools/plot_tf.rst
   :start-after: .. _sphx_glr_examples_emtools_plot_tf.py:
   :end-before: .. _sphx_glr_download_examples_emtools_plot_tf.py:

Advanced: working directly from Spectra objects
-------------------------------------------------

:func:`~pycsamt.emtools.tf.plot_induction_map_from_spectra`,
:func:`~pycsamt.emtools.tf.plot_tipper_polar_from_spectra`, and
:func:`~pycsamt.emtools.tf.plot_induction_rose_from_spectra` are not
demonstrated above because they serve a different, earlier stage of
processing: they read tipper straight out of a
:class:`~pycsamt.seg.spectra.Spectra` object (cross-power/spectral
estimates), before ``Z``/``Tipper`` have been assembled into an EDI. Use
these when you are working from raw processed spectra rather than EDI
files — station coordinates are supplied explicitly via a ``coords``
dict since a bare ``Spectra`` carries no position information::

    from pycsamt.emtools.tf import plot_induction_map_from_spectra

    plot_induction_map_from_spectra(
        {"site_a": spectra_a, "site_b": spectra_b},
        coords={"site_a": (0.0, 0.0), "site_b": (500.0, 0.0)},
        period=100.0,
    )
