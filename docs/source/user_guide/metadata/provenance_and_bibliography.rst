.. _user_guide_metadata_provenance_and_bibliography:

Provenance and Bibliography
===========================

.. currentmodule:: pycsamt.metadata

:class:`ProvenanceMeta` records how a transfer-function document came
to be -- creation time, creating application, creator, and submitter
-- so that history survives a format conversion instead of being reset
every time a file is re-saved. :class:`Reference`, :class:`CopyrightInfo`,
:class:`Person`, and :class:`Software` supply the bibliographic detail
that a provenance or copyright record is built from: who to cite,
under what license, and with which software. Both groups are
:term:`format-neutral metadata`: :attr:`~pycsamt.emtf.document.EMTF.provenance`
and :attr:`~pycsamt.emtf.document.EMTF.copyright` are read from EDI or
EMTF XML through the same objects, and the same objects are written
back out regardless of target format.

Two Real Provenance Records
---------------------------

Two real, unrelated surveys in the repository exercise ``ProvenanceMeta``
from opposite directions. ``data/gv_data/xml/gv100.xml`` is an EMTF
XML document already stored in the repository for a 2020 USGS
magnetotelluric survey in Gabbs Valley, Nevada; it was produced by
converting the genuine field EDI (``data/gv_data/gv_final_edi/gv100.edi``,
itself written by `MTpy <https://github.com/MTgeophysics/mtpy>`_
during USGS processing) through pycsamt's own EDI adapter, so reading
it back exercises the XML reader on real content rather than a
hand-written fixture. ``data/AMT/WILLY_DATA/L18PLT`` is the historical
EDI line used throughout :doc:`site_and_survey`, produced by a
different, unrelated processing chain (``MTPROC``), and read here
straight from its native EDI form.

.. note::

   The Gabbs Valley data is a public-domain USGS release; attribution
   is required for any use. Cite:

       Peacock, J. R., Siler, D. L., Dean, B., Zielinski, L., 2021,
       Magnetotelluric Data from Gabbs Valley Region, Nevada, 2020:
       U.S. Geological Survey data release,
       `<https://doi.org/10.5066/P9GZ9Z56>`__

.. code-block:: pycon

   >>> from pathlib import Path
   >>> from pycsamt.emtf import EMTF
   >>> from pycsamt.seg.edi import EDIFile

   >>> tf_xml = EMTF.from_xml("data/gv_data/xml/gv100.xml")
   >>> print(tf_xml.provenance)
   ProvenanceMeta(create_time='2021-06-24T20:24:04.077672+00:00', creating_application='1.1.5', creator=<Person name='Jared Peacock' email=None>, submitter=None, extra=dict(len=0, keys=[]))
   >>> print(tf_xml.copyright)
   None

   >>> edi = EDIFile(Path("data/AMT/WILLY_DATA/L18PLT/18-001A.edi"))
   >>> tf_edi = EMTF.from_edi(edi)
   >>> print(tf_edi.provenance)
   ProvenanceMeta(create_time='11/9/23', creating_application='MTPROC V1.0.7', creator=None, submitter=None, extra=dict(len=0, keys=[]))
   >>> print(tf_edi.copyright)
   None

Both ``create_time`` values are kept exactly as the source file wrote
them -- an ISO-8601 timestamp with microseconds from the XML file, a
bare ``"11/9/23"`` from the EDI ``HEAD.FILEDATE`` field -- rather than
being parsed into a stricter type. :class:`ProvenanceMeta` documents
this choice explicitly: inventing a timezone or a day/month order for
an ambiguous historical date would add false precision, so the string
is preserved verbatim and left for the caller to interpret. The EDI
record also has no ``creator``: ``HEAD.FILEBY`` is genuinely absent
from this file, and the EDI adapter leaves the field ``None`` rather
than guessing an author from context such as ``HEAD.ACQBY``.

Missing Copyright Is Not Fabricated
-----------------------------------

Both records above show ``copyright`` as ``None``, and this is
information, not a parsing failure. Neither source file contains a
``<Copyright>``-equivalent block, and pycsamt does not synthesize one
to hold a plausible-looking citation. The Gabbs Valley file makes the
same point from the opposite direction: it *does* carry rotation
history at the document root,

.. code-block:: pycon

   >>> print(tf_xml.orientation.rotation_info)
   EDI ZROT is constant at 347.5 deg relative to HEAD.COORDSYS='Geomagnetic North'; it was not promoted to an angle_to_geographic_north because historical EDI rotation metadata can be ambiguous.

and the EMTF XML writer preserves exactly that shape on a round trip
-- a bare ``<RotationInfo>`` element, with no ``<Copyright>`` wrapper
invented to hold it:

.. code-block:: pycon

   >>> from tempfile import TemporaryDirectory

   >>> with TemporaryDirectory() as tmp:
   ...     out = Path(tmp) / "roundtrip.xml"
   ...     tf_xml.write_xml(out)
   ...     text = out.read_text(encoding="utf-8")
   ...     print("<Copyright>" in text, "<RotationInfo>" in text)
   ...     reloaded = EMTF.from_xml(out)
   ...     print(reloaded.copyright)
   ...     print(reloaded.orientation.rotation_info == tf_xml.orientation.rotation_info)
   False True
   None
   True

This is deliberate, not incidental: the serializer only opens a
``<Copyright>`` element when ``document.copyright`` is set, and writes
``RotationInfo`` directly under the document root otherwise, precisely
so that attaching real rotation history never forces a fabricated
citation into existence.

Building a Citation by Hand
---------------------------

pycsamt cannot invent a citation for the Gabbs Valley file, but the
public USGS data release
(``data/gv_data/gv_survey_v2.xml``, an FGDC metadata record distributed
with this dataset) has everything needed to build one by hand and
attach it explicitly:

.. code-block:: pycon

   >>> from pycsamt.metadata import Reference, CopyrightInfo, Person, Software
   >>> from pycsamt.metadata import MetadataError

   >>> reference = Reference(
   ...     author="Peacock, J. R.; Siler, D. L.; Dean, B. J.; Zielinski, L. A.",
   ...     title="Magnetotelluric Data from the Gabbs Valley Region, Nevada, 2020",
   ...     year=2021,
   ...     doi="10.5066/P9GZ9Z56",
   ... )
   >>> print(reference)
   <Reference title='Magnetotelluric Data from the Gabbs Valley Region, Nevada, 2020' author='Peacock, J. R.; Siler, D. L.; Dean, B. J.; Zielinski, L. A.' year=2021>
   >>> print(reference.to_bibtex("Peacock2021"))
   @article{Peacock2021,
     author = {Peacock, J. R.; Siler, D. L.; Dean, B. J.; Zielinski, L. A.},
     title = {Magnetotelluric Data from the Gabbs Valley Region, Nevada, 2020},
     year = {2021},
     doi = {10.5066/P9GZ9Z56},
   }

``author``, ``title``, ``year``, and ``doi`` above are copied directly
from the FGDC ``<citeinfo>`` block (``origin``, ``title``, ``pubdate``,
``onlink``); nothing is invented. ``doi`` is validated against the
real DOI syntax at construction, so a malformed identifier is rejected
immediately instead of silently stored:

.. code-block:: pycon

   >>> try:
   ...     Reference(author="x", title="y", doi="not-a-doi")
   ... except MetadataError as exc:
   ...     print(type(exc).__name__, exc)
   MetadataError Invalid DOI format: not-a-doi

``MetadataError`` is a plain ``Exception`` subclass, not a
``ValueError`` -- unlike the ``ValueError`` raised by
:class:`LocationMeta`/:class:`SiteMeta` in :doc:`site_and_survey`, so
catch it explicitly rather than assuming every metadata class in this
package signals invalid input the same way.

:class:`Person` and :class:`Software` fill out the picture: a contact
for the release, and the software that produced it (recovered here
from the plain-text breadcrumbs the same file left in ``FieldNotes``
-- ``provenance.software.name = MTpy``, ``provenance.software.author =
Jared Peacock`` -- since the structured ``<Provenance>`` block only
captured the creator's name, not the software that ran):

.. code-block:: pycon

   >>> contact = Person(
   ...     name="Jared R. Peacock",
   ...     email="jpeacock@usgs.gov",
   ...     organization="U.S. Geological Survey",
   ... )
   >>> print(contact)
   <Person name='Jared R. Peacock' email='jpeacock@usgs.gov'>

   >>> copyright_info = CopyrightInfo(
   ...     release_status="Unrestricted Release",
   ...     conditions_of_use=(
   ...         "Any use of trade, firm, or product names is for descriptive "
   ...         "purposes only and does not imply endorsement by the U.S. Government."
   ...     ),
   ...     reference=reference,
   ... )
   >>> print(copyright_info)
   <CopyrightInfo status='Unrestricted Release'>

   >>> software = Software(
   ...     name="MTpy",
   ...     version="metadata branch",
   ...     author=Person(name="Jared Peacock"),
   ... )
   >>> print(software)
   <Software name='MTpy' version='metadata branch'>

``conditions_of_use`` is condensed from the FGDC ``<useconst>`` field,
which in full also disclaims USGS liability for downstream use --
shortened here for the example, but a real delivery should carry the
complete text. Unlike :class:`ProvenanceMeta`, which inherits
:class:`~pycsamt.api.property.PyCSAMTObject`, these four bibliographic
classes are plain dataclasses with their own hand-written ``__repr__``
and ``to_dict``; that is why their compact display does not follow the
same "cap at eight fields, append ``...``" convention documented in
:doc:`site_and_survey`.

Enriching and Re-Writing a Document
-----------------------------------

Attaching the hand-built citation to the document and writing it back
out round-trips cleanly, and now a real ``<Copyright>`` element is
present because a real ``CopyrightInfo`` was actually supplied:

.. code-block:: pycon

   >>> tf_xml.copyright = copyright_info
   >>> with TemporaryDirectory() as tmp:
   ...     out = Path(tmp) / "gv100_enriched.xml"
   ...     tf_xml.write_xml(out)
   ...     text = out.read_text(encoding="utf-8")
   ...     print("<Copyright>" in text)
   ...     reloaded = EMTF.from_xml(out)
   ...     print(reloaded.copyright.reference.doi)
   ...     print(reloaded.copyright.reference.author)
   ...     print(reloaded.copyright.release_status)
   ...     print(reloaded.orientation.rotation_info == tf_xml.orientation.rotation_info)
   True
   10.5066/P9GZ9Z56
   Peacock, J. R.; Siler, D. L.; Dean, B. J.; Zielinski, L. A.
   Unrestricted Release
   True

The pre-existing ``RotationInfo`` at the document root survives
unchanged alongside the new ``<Copyright>`` element -- attaching a
citation later does not disturb metadata that was already there.
``tf_xml.copyright = copyright_info`` mutates the in-memory document
directly; nothing is written to disk until :meth:`~pycsamt.emtf.document.EMTF.write_xml`
is called, so the enrichment step and the export step remain two
separate, inspectable decisions.

A Second, Unrelated Software/Copyright Model
--------------------------------------------

.. note::

   :class:`pycsamt.seg.property.Software` and
   :class:`pycsamt.seg.property.Copyright` are **not** the classes
   documented on this page, despite the identical names. They are the
   historical EDI ``>HEAD``/``>INFO`` section holders that
   :class:`~pycsamt.seg.edi.EDIFile` uses internally (for example, the
   ``PYCSAMT``-named :class:`~pycsamt.seg.property.Software` instance
   :class:`~pycsamt.seg.heads.Head` stamps into ``ProcessingSoftware``
   when writing a file) and predate the format-neutral
   :mod:`pycsamt.metadata` package. The EMTF-XML implementation plan
   already flags this as an open reconciliation item for
   :class:`~pycsamt.metadata.processing.ProcessingMeta`, deferred
   rather than resolved in passing. Import from :mod:`pycsamt.metadata`
   for anything read from or written to an :class:`EMTF` document; the
   ``pycsamt.seg.property`` classes stay internal to the EDI writer.

Choosing the Right Object
-------------------------

.. list-table::
   :header-rows: 1
   :widths: 32 32 36

   * - Need
     - Object
     - Notes
   * - Who/when/what created a document
     - :class:`ProvenanceMeta`
     - ``create_time`` kept verbatim; maps to EMTF XML ``<Provenance>``.
   * - A citation for a dataset
     - :class:`Reference`
     - Validates DOI syntax; :meth:`~Reference.to_bibtex` for export.
   * - License/usage terms plus a citation
     - :class:`CopyrightInfo`
     - Maps to EMTF XML ``<Copyright>``; never fabricated when absent.
   * - Contact information for an author, submitter, or point of contact
     - :class:`Person`
     - Reused inside both ``ProvenanceMeta`` and ``CopyrightInfo``.
   * - The software that produced a document
     - :class:`Software`
     - Standalone record; not the same class as ``pycsamt.seg.property.Software``.

Next Steps
----------

* :doc:`site_and_survey` covers station identity and campaign-level
  aggregation;
* :doc:`channels_and_orientation` covers the channel geometry and
  orientation metadata that give a transfer-function matrix physical
  meaning;
* :doc:`processing_and_quality` covers how a transfer function was
  estimated and how good the result is.
