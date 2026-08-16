.. _user_guide_geology:

Geology
=======

:mod:`pycsamt.geology` holds the earth-science domain knowledge that
:mod:`pycsamt.interp` builds on but does not itself define: rock
resistivity classification, ground-truth borehole logs, and field
structural measurements. Nothing here depends on a resistivity model,
an inversion result, or any other electromagnetic concept -- these
pages describe geological data and reasoning that would be exactly as
useful on a project with no EM survey at all.

Use this section when you need to build or extend the rock-property
table behind a classification, load ground-truth borehole logs, or
record field structural evidence against a profile. For the workflow
that turns an EM resistivity model into a calibrated, geologically
classified section -- the part that *does* know about resistivity
models and misfit review -- see :doc:`../interpretation/index`. For
the complete callable reference, see :doc:`../../api/geology`.

.. toctree::
   :maxdepth: 3
   :class: pycsamt-guide-toc

   concepts
   rock_database
   borehole
   structural
