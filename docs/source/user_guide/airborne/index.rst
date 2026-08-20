.. _user_guide_airborne:

Airborne EM Guide
=================

:mod:`pycsamt.airborne` is the technology-neutral airborne data model:
flight lines and datasets built from :class:`~pycsamt.emtf.EMTF`
documents, with no EDI bridge and no impedance requirement, since a
genuine airborne EM measurement rarely has one. It sits alongside
:mod:`pycsamt.site` (ground, EDI-shaped) rather than inside it, and
this section covers that data model on its own terms: building and
inspecting flight lines and datasets, reading them through the
``Sites``-shaped :class:`~pycsamt.airborne.site.AirborneSite` view,
registering or detecting a technology and a native format, and
running structural QC. Start with :doc:`overview` for the full
picture of how AFMAG, ZTEM, and MobileMT map onto this shared model.

.. toctree::
   :maxdepth: 3
   :class: pycsamt-guide-toc

   overview
   data_model
   site
   registry_and_io
   quality_control
