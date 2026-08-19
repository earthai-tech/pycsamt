.. _user_guide_emtf_xml:

Reading and Writing EMTF XML
============================

.. currentmodule:: pycsamt.emtf

:class:`EMTFXMLReader` and :class:`EMTFXMLWriter` implement the modern
EMTF XML format on top of the Python standard library's
``xml.etree.ElementTree`` -- no mandatory XML dependency beyond what
already ships with Python. They are the mechanism behind
:meth:`EMTF.from_xml <pycsamt.emtf.document.EMTF.from_xml>` and
:meth:`EMTF.write_xml <pycsamt.emtf.document.EMTF.write_xml>` /
:meth:`EMTF.to_xml <pycsamt.emtf.document.EMTF.to_xml>`, and this page
exercises them directly rather than through those convenience wrappers,
against the same real, externally-produced-then-pycsamt-converted file
already used across :doc:`../metadata/index` --
``data/gv_data/xml/gv100.xml``, public-domain USGS data
(https://doi.org/10.5066/P9GZ9Z56). ``EMTFXMLReader``,
``EMTFXMLWriter``, ``EMTFXMLParseError``, and ``EMTFXMLSerializationError``
all import from the top level: ``from pycsamt.emtf import (EMTF,
EMTFXMLReader, EMTFXMLWriter, EMTFXMLParseError,
EMTFXMLSerializationError, read_emtf_xml, write_emtf_xml)``.

Reading a Real Document Directly
--------------------------------

:meth:`EMTFXMLReader.read` accepts a path, a string, bytes, or an open
stream, and returns a fully populated :class:`~pycsamt.emtf.document.EMTF`
-- the same object :doc:`../metadata/provenance_and_bibliography` and
:doc:`../metadata/channels_and_orientation` already built through the
``EMTF.from_xml`` shortcut:

.. code-block:: pycon

   >>> from pycsamt.emtf import EMTFXMLReader

   >>> doc = EMTFXMLReader(strict=True).read("data/gv_data/xml/gv100.xml")
   >>> print(doc.site.site_id, doc.n_periods)
   gv100 48

Strict Rejects, Permissive Recovers
-----------------------------------

``strict=True`` (the default) rejects malformed content outright. A
document whose root element is not ``<EM_TF>`` fails immediately:

.. code-block:: pycon

   >>> from pycsamt.emtf import EMTFXMLParseError

   >>> bad_root = "<NOT_EM_TF><Description>x</Description></NOT_EM_TF>"
   >>> try:
   ...     EMTFXMLReader(strict=True).read(bad_root)
   ... except EMTFXMLParseError as exc:
   ...     print(type(exc).__name__, exc)
   EMTFXMLParseError <string>: root element must be <EM_TF>, got <NOT_EM_TF>

``strict=False`` turns the same problem into a warning and keeps
going, reading whatever the (technically wrong) document still
contains:

.. code-block:: pycon

   >>> import warnings

   >>> with warnings.catch_warnings(record=True) as caught:
   ...     warnings.simplefilter("always")
   ...     doc_bad_root = EMTFXMLReader(strict=False).read(bad_root)
   ...     print(len(caught), caught[0].category.__name__)
   1 EMTFXMLWarning
   >>> print(doc_bad_root.description)
   x

.. note::

   Building this exact permissive-mode example surfaced a real gap:
   an invalid ``<Period>`` value used to warn and then still crash
   with an unrelated, uncaught ``ValueError`` even under
   ``strict=False``, because the reader substituted ``NaN`` for the
   bad period but :class:`~pycsamt.emtf.document.EMTF`'s own
   validation rejects any ``NaN`` period unconditionally. This has
   been fixed: permissive mode now drops the invalid period -- and
   only that period's data row -- instead, so it returns a genuinely
   usable, shorter document rather than crashing. Strict mode is
   unaffected; it still rejects the file immediately, as shown next.

A two-period document with one negative period demonstrates both
sides of that fix directly. Strict mode still refuses the whole file:

.. code-block:: pycon

   >>> minimal_two = '''<?xml version="1.0"?>
   ... <EM_TF>
   ...   <Site><Id>TEST001</Id></Site>
   ...   <Data count="2">
   ...     <Period value="-1.0" units="secs">
   ...       <Z type="complex" size="2 2" units="[mV/km]/[nT]">
   ...         <value name="Zxx" output="Ex" input="Hx">1.0 2.0</value>
   ...         <value name="Zxy" output="Ex" input="Hy">3.0 4.0</value>
   ...         <value name="Zyx" output="Ey" input="Hx">5.0 6.0</value>
   ...         <value name="Zyy" output="Ey" input="Hy">7.0 8.0</value>
   ...       </Z>
   ...     </Period>
   ...     <Period value="10.0" units="secs">
   ...       <Z type="complex" size="2 2" units="[mV/km]/[nT]">
   ...         <value name="Zxx" output="Ex" input="Hx">9.0 10.0</value>
   ...         <value name="Zxy" output="Ex" input="Hy">11.0 12.0</value>
   ...         <value name="Zyx" output="Ey" input="Hx">13.0 14.0</value>
   ...         <value name="Zyy" output="Ey" input="Hy">15.0 16.0</value>
   ...       </Z>
   ...     </Period>
   ...   </Data>
   ... </EM_TF>'''
   >>> try:
   ...     EMTFXMLReader(strict=True).read(minimal_two)
   ... except EMTFXMLParseError as exc:
   ...     print(type(exc).__name__, exc)
   EMTFXMLParseError <string>: Period[0] must be positive

Permissive mode drops the bad period and its ``Z`` row, and returns
the remaining one cleanly -- the second period's real values, not a
placeholder:

.. code-block:: pycon

   >>> with warnings.catch_warnings(record=True) as caught:
   ...     warnings.simplefilter("always")
   ...     doc_dropped = EMTFXMLReader(strict=False).read(minimal_two)
   ...     print(len(caught), caught[0].category.__name__)
   ...     print(caught[0].message)
   1 EMTFXMLWarning
   <string>: Period[0] must be positive
   >>> print(doc_dropped.n_periods, doc_dropped.periods)
   1 [10.]
   >>> print(doc_dropped.z[0])
   [[ 9.+10.j 11.+12.j]
    [13.+14.j 15.+16.j]]

Not every problem is treated the same way even within one mode. An
unrecognized period unit is *always* just a warning, even under
``strict=True``, because pycsamt has an unambiguous fallback
(assume seconds) rather than a genuine ambiguity to reject:

.. code-block:: pycon

   >>> minimal_units = minimal_two.replace(
   ...     'value="-1.0" units="secs">',
   ...     'value="1.0" units="days">',
   ... )
   >>> with warnings.catch_warnings(record=True) as caught:
   ...     warnings.simplefilter("always")
   ...     doc_units = EMTFXMLReader(strict=True).read(minimal_units)
   ...     print(len(caught), caught[0].category.__name__)
   1 EMTFXMLWarning
   >>> print(doc_units.n_periods)
   2

Missing Stays Missing on Round Trip
-----------------------------------

:doc:`transfer_functions` already showed that gv100's real tipper has
14 genuinely missing (``NaN``, not zero) periods out of 48, split
between the shortest and longest ends of the band. Writing that same
document to XML does not paper over the gap -- it omits the ``<T>``
element for those periods entirely, not just its ``<value>`` children:

.. code-block:: pycon

   >>> import tempfile
   >>> from pathlib import Path
   >>> import numpy as np
   >>> import xml.etree.ElementTree as ET
   >>> from pycsamt.emtf import EMTF

   >>> tf_gv = EMTF.from_edi("data/gv_data/gv_final_edi/gv100.edi")
   >>> tip = tf_gv.get_transfer_function("tipper")
   >>> finite = np.isfinite(tip.data.real) & np.isfinite(tip.data.imag)
   >>> nan_idx = np.where(~finite.all(axis=(1, 2)))[0]
   >>> print(nan_idx.size)
   14

   >>> with tempfile.TemporaryDirectory() as tmp:
   ...     out = Path(tmp) / "gv100.xml"
   ...     tf_gv.write_xml(out)
   ...     root = ET.fromstring(out.read_text(encoding="utf-8"))
   ...     period_nodes = root.find("Data").findall("Period")
   ...     print(period_nodes[nan_idx[0]].find("T") is None)
   ...     print(period_nodes[20].find("T") is not None)
   ...     reloaded = EMTF.from_xml(out)
   ...     print(np.allclose(tip.data, reloaded.tipper, equal_nan=True))
   True
   True
   True

The first ``print`` confirms no ``<T>`` element at all at a period
where tipper is missing; the second confirms one *is* present at an
ordinary period; the third confirms the full array -- gaps included
-- survives the round trip exactly. This is the same discipline
:doc:`edi_interop` and :doc:`transfer_functions` describe for EDI and
for the in-memory model: a missing value is information, and writing
it as anything else -- zero, an empty tag, a placeholder -- would
destroy that information permanently.

Deterministic Element Order and Precision
-----------------------------------------

:class:`~pycsamt.emtf.xml.serializer.EMTFXMLSerializer` (used
internally by ``EMTFXMLWriter``) always emits top-level elements in
the same fixed order, so two runs of the same document produce a
comparable diff instead of an arbitrarily shuffled one:

.. code-block:: pycon

   >>> order = [child.tag for child in root]
   >>> print(order[:8], "...", order[-4:])
   ['SubType', 'Tags', 'Provenance', 'RotationInfo', 'Site', 'FieldNotes', 'FieldNotes', 'FieldNotes'] ... ['DataTypes', 'SiteLayout', 'Data', 'PeriodRange']

(gv100's real ``<FieldNotes>`` repeats once per raw EDI ``>INFO`` line
-- dozens of them -- which is why the middle of that list is elided
above.) ``precision`` (default 17, the digits needed to round-trip an
IEEE float64 exactly) controls response *data* values only -- period
values always use a fixed 15 significant digits regardless of
``precision``, since a period is a coordinate, not a noisy
measurement:

.. code-block:: pycon

   >>> from pycsamt.emtf import EMTFXMLWriter

   >>> xml_17 = EMTFXMLWriter(strict=True, precision=17).dumps(tf_gv)
   >>> xml_6 = EMTFXMLWriter(strict=True, precision=6).dumps(tf_gv)
   >>> import re
   >>> print(re.search(r'<value name="Zxx"[^>]*>([^<]+)</value>', xml_17).group(1))
   -86.575890000000001 -714.197
   >>> print(re.search(r'<value name="Zxx"[^>]*>([^<]+)</value>', xml_6).group(1))
   -86.5759 -714.197
   >>> print(re.search(r'<Period value="([^"]+)"', xml_17).group(1))
   0.00130209994867122
   >>> print(re.search(r'<Period value="([^"]+)"', xml_6).group(1))
   0.00130209994867122

When Writing Refuses
--------------------

Serialization has its own unconditional checks, independent of
``strict``. A document that holds a transfer function but no period
vector at all cannot be written as response data, because there is
nothing to attach each row's period to:

.. code-block:: pycon

   >>> from pycsamt.emtf import TransferFunction, EMTFXMLSerializationError

   >>> tf_noperiods = TransferFunction(
   ...     name="impedance", data=np.ones((1, 2, 2), dtype=complex),
   ...     input_channels=("Hx", "Hy"), output_channels=("Ex", "Ey"),
   ...     periods=None,
   ... )
   >>> doc_empty = EMTF()
   >>> doc_empty.add_transfer_function(tf_noperiods)
   >>> try:
   ...     EMTFXMLWriter(strict=True).dumps(doc_empty)
   ... except EMTFXMLSerializationError as exc:
   ...     print(type(exc).__name__, exc)
   EMTFXMLSerializationError cannot serialize transfer-function data without periods

``strict=False`` on the writer downgrades this to a warning and
simply omits the ``<Data>``/``<PeriodRange>`` elements rather than
inventing a period vector -- consistent with every other "missing is
not fabricated" example on this page.

Full Round-Trip Verification
----------------------------

``XML -> EMTF -> XML`` on an already-canonical file is not just
numerically close -- it is a fixed point. Reading gv100.xml, writing
it back out, reading *that*, and writing again produces byte-identical
text the second time:

.. code-block:: pycon

   >>> from pycsamt.emtf import read_emtf_xml, write_emtf_xml

   >>> with tempfile.TemporaryDirectory() as tmp:
   ...     doc1 = EMTF.from_xml("data/gv_data/xml/gv100.xml")
   ...     out1 = Path(tmp) / "a.xml"
   ...     write_emtf_xml(doc1, out1)
   ...     doc2 = read_emtf_xml(out1)
   ...     out2 = Path(tmp) / "b.xml"
   ...     write_emtf_xml(doc2, out2)
   ...     print(out1.read_text() == out2.read_text())
   True

That guarantee only holds once a document has already passed through
the writer once (a freshly EDI-loaded document may reorder or
reformat fields the first time, exactly as :doc:`edi_interop` shows
for EDI's own ``HEAD.LOC`` casing). From the second write onward,
there is nothing left to normalize.

Choosing the Right Function
---------------------------

.. list-table::
   :header-rows: 1
   :widths: 32 30 38

   * - Need
     - Function/class
     - Notes
   * - Read from a document (usual case)
     - :meth:`EMTF.from_xml <pycsamt.emtf.document.EMTF.from_xml>` / :func:`read_emtf_xml`
     - Both wrap :class:`EMTFXMLReader`.
   * - Read with explicit control over strictness
     - :class:`EMTFXMLReader`
     - Reuse one instance to read several files under the same policy.
   * - Write to a path or stream
     - :meth:`EMTF.write_xml <pycsamt.emtf.document.EMTF.write_xml>` / :func:`write_emtf_xml`
     - Both wrap :class:`EMTFXMLWriter`.
   * - Write to an in-memory string or element tree
     - :meth:`EMTF.to_xml <pycsamt.emtf.document.EMTF.to_xml>` / ``EMTFXMLWriter.dumps``/``to_element``
     - No filesystem access required.

Next Steps
----------

This closes the current pass through :doc:`index`. From here:

* :doc:`edi_interop` covers the EDI adapter that produced
  ``gv100.xml`` in the first place, and the same missing-value
  discipline applied to EDI instead of XML;
* :doc:`transfer_functions` covers the in-memory
  ``TransferFunction``/``StatisticalEstimate`` model this reader
  populates;
* :doc:`../metadata/index` covers every metadata object
  (``SiteMeta``, ``ProvenanceMeta``, ``SiteLayout``, ``OrientationMeta``,
  ``ProcessingMeta``) this reader and writer carry to and from XML.
