.. _user_guide_topo:

Topography
==========

The :mod:`pycsamt.topo` package embeds real station elevation into 2-D
resistivity sections so ``pcolormesh``/``imshow`` plots drape correctly
over terrain instead of assuming a flat ``z = 0`` datum. It underlies
terrain-aware section rendering used from :doc:`../map/index`,
:doc:`../inversion/index`, and :doc:`../ai_inversion/index`, and also
exposes a one-call section function that accepts raw arrays, a
:class:`pycsamt.interp.ResistivityModel`, a
:class:`pycsamt.inversion.results.InversionResult`, native Occam2D or
ModEM results, or AI inversion results directly.

Use this section when a section plot -- pseudosection or inversion
depth section -- needs to reflect the real ground surface instead of a
flat datum, or when you need fine control over how terrain is
extracted, draped, and rendered.

.. toctree::
   :maxdepth: 3
   :class: pycsamt-guide-toc

   concepts
   extract
   drape
   overlay
   section
   api
