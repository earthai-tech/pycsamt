.. _user_guide_emtf:

EMTF
====

:mod:`pycsamt.emtf` is the format-neutral scientific core for
electromagnetic transfer functions: matrix-oriented
:class:`~pycsamt.emtf.transfer.TransferFunction` and
:class:`~pycsamt.emtf.estimates.StatisticalEstimate` objects held
inside one :class:`~pycsamt.emtf.document.EMTF` document, together
with the adapters that read and write it as historical SEG EDI, EDI
SPECTRA, or modern EMTF XML. The central design rule, following
Kelbert (2020)'s EMTF XML paper: **EDI and EMTF XML are formats; the
transfer function is the scientific object.** Impedance, tipper,
variance, and full covariance are represented once, independently of
which file produced or will store them, and a conversion between
formats never silently discards or fabricates content -- a missing
element stays missing rather than becoming a zero or an invented
default.

This section assumes the metadata objects covered in
:doc:`../metadata/index` -- ``SiteMeta``, ``ProvenanceMeta``,
``SiteLayout``, ``OrientationMeta``, ``ProcessingMeta``,
``TransferFunctionQuality`` -- since ``EMTF`` holds exactly those
objects as its document metadata rather than reinventing a
format-specific metadata model. Use this section when the question is
about the transfer function itself: its matrix shape and channel
axes, how it survives a round trip between EDI and XML, how EDI
SPECTRA's cross-power matrices are recovered into full covariance, or
how rotation is applied without disturbing the physical site layout
that produced it. For the complete callable reference, see
:doc:`../../api/emtf`.

.. toctree::
   :maxdepth: 3
   :class: pycsamt-guide-toc

   document
   transfer_functions
   edi_interop
   spectra
   xml
   rotation
