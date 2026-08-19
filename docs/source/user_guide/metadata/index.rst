.. _user_guide_metadata:

Metadata
========

:mod:`pycsamt.metadata` holds the format-neutral metadata objects that
describe *how a transfer function came to be*, independently of
whether it is currently stored as EDI, EMTF XML, or an in-memory
:class:`~pycsamt.core.base.TFBundle`. A station's coordinates, the
instrument that recorded it, the channels and orientation that give a
:math:`2\times2` impedance matrix physical meaning, who processed it
and with what remote reference, and how good the result actually is
are all scientific facts about the measurement -- not serialization
details of one file format. Keeping them in one reusable package
means EDI workflows, the EMTF XML reader/writer, and forward-modelling
priors all share the same classes instead of each format inventing its
own metadata hierarchy.

Use this section when you need to build, inspect, or round-trip any of
these objects directly -- constructing a :class:`~pycsamt.metadata.survey.SurveyMeta`
for a new campaign, comparing two :class:`~pycsamt.metadata.instrument.InstrumentMeta`
presets, or reading what a :class:`~pycsamt.metadata.quality.DataQuality`
flag actually means. When you only need the *effect* these objects have
on a concrete workflow, see :doc:`../site/index` for station-level EDI
and EMTF-XML metadata editing, or :doc:`../emtf/index` for how
:class:`~pycsamt.metadata.provenance.ProvenanceMeta`,
:class:`~pycsamt.metadata.channels.SiteLayout`, and
:class:`~pycsamt.metadata.orientation.OrientationMeta` map onto EMTF
XML's ``<Provenance>``, ``<SiteLayout>``, and ``<Orientation>``
elements. For the complete callable reference, see
:doc:`../../api/metadata`.

.. toctree::
   :maxdepth: 3
   :class: pycsamt-guide-toc

   site_and_survey
   provenance_and_bibliography
   channels_and_orientation
   instrument
   processing_and_quality
   frequency_bands
