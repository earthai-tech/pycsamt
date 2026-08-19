.. _user_guide_emtf_edi_interop:

EDI Interoperability
====================

.. currentmodule:: pycsamt.emtf

:meth:`EMTF.from_edi <pycsamt.emtf.document.EMTF.from_edi>` and
:meth:`EMTF.to_edi <pycsamt.emtf.document.EMTF.to_edi>` bridge the
historical :class:`~pycsamt.seg.edi.EDIFile` model and the
format-neutral :class:`EMTF` document. The direction into ``EMTF`` is
the one already exercised piecemeal throughout :doc:`../metadata/index`
and :doc:`document`; this page assembles the full picture, then covers
the return trip -- which pycsamt treats as potentially lossy on
purpose, since standard EDI genuinely cannot represent everything an
``EMTF`` document can hold. Everything on this page imports directly
from the top level: ``from pycsamt.emtf import EMTF, EMTFEDIConversionError,
DataLossWarning, write_edi, bundle_to_emtf, emtf_to_bundle``.

What EDI Extraction Preserves
-----------------------------

Loading the real Gabbs Valley station (public-domain USGS data,
https://doi.org/10.5066/P9GZ9Z56, already used throughout
:doc:`../metadata/index`) populates four of ``EMTF``'s seven metadata
slots from one EDI file, in a single call --
:attr:`~pycsamt.emtf.document.EMTF.provenance`,
:attr:`~pycsamt.emtf.document.EMTF.site`,
:attr:`~pycsamt.emtf.document.EMTF.site_layout`, and
:attr:`~pycsamt.emtf.document.EMTF.orientation`; ``copyright``,
``processing``, and ``quality`` all stay ``None`` because this
particular file has nothing explicit for any of them:

.. code-block:: pycon

   >>> from pathlib import Path
   >>> from pycsamt.emtf import EMTF

   >>> tf_gv = EMTF.from_edi(Path("data/gv_data/gv_final_edi/gv100.edi"))
   >>> print(tf_gv.site.project, tf_gv.site.survey, tf_gv.site.name)
   Energy Resources Program GV2020 Gabbs Valley
   >>> print(tf_gv.processing)
   None
   >>> print(tf_gv.site_layout.input_names, tf_gv.site_layout.output_names)
   ('Hx', 'Hy') ('Hz', 'Ex', 'Ey')
   >>> print(tf_gv.orientation.rotation_info)
   EDI ZROT is constant at 347.5 deg relative to HEAD.COORDSYS='Geomagnetic North'; it was not promoted to an angle_to_geographic_north because historical EDI rotation metadata can be ambiguous.

``tf_gv.processing`` is ``None`` for a reason already established in
:doc:`../metadata/processing_and_quality`: this particular file's raw
``>INFO`` text has no structured processing keys at all, so there is
nothing explicit to extract. ``site``, ``site_layout``, and
``orientation`` all come from structured EDI sections
(``HEAD``/``MTSECT``/``DEFINEMEAS``/``ZROT``) that this file *does*
populate. Every one of these extractions follows the same rule
demonstrated repeatedly across :doc:`../metadata/index`: a field is
set only when the corresponding EDI key is explicitly present, never
inferred from context or defaulted from pycsamt's own conventions.

Variance and the Legacy ``z_err`` Convention
--------------------------------------------

:doc:`document` showed that ``tf_gv.z_err`` is ``None`` even though a
``VAR`` estimate is attached, because pycsamt's historical ``z_err``
field and the EMTF ``VAR`` estimate are related but not
interchangeable. The adapter's own docstring states the exact
relationship: pycsamt has historically stored
``Z.z_err = sqrt(EDI complex variance)``, so ``edi_to_emtf``
reconstructs ``VAR`` as ``z_err ** 2`` exactly -- meaning
``np.sqrt(variance)`` (as done in :doc:`document`) recovers that
legacy value precisely, not approximately. The *per-component* real or
imaginary standard error, when that specific statistic is needed, is
``sqrt(VAR / 2)`` instead -- a different, more refined quantity than
the legacy ``z_err`` convention computes, and worth knowing before
reaching for one when the other was actually wanted.

Converting Back to EDI
----------------------

A document that only ever held what the source EDI already expressed
round-trips through :meth:`~pycsamt.emtf.document.EMTF.to_edi` without
any warning at all:

.. code-block:: pycon

   >>> import warnings

   >>> with warnings.catch_warnings(record=True) as caught:
   ...     warnings.simplefilter("always")
   ...     edi_back = tf_gv.to_edi()
   ...     print(len(caught))
   0

Attaching something EDI cannot hold -- a citation, in this case, built
the same way as in :doc:`../metadata/provenance_and_bibliography` --
makes the same call emit a :class:`DataLossWarning` naming exactly
what will not survive:

.. code-block:: pycon

   >>> from pycsamt.metadata import Reference, CopyrightInfo

   >>> reference = Reference(
   ...     author="Peacock, J. R.; Siler, D. L.; Dean, B. J.; Zielinski, L. A.",
   ...     title="Magnetotelluric Data from the Gabbs Valley Region, Nevada, 2020",
   ...     year=2021, doi="10.5066/P9GZ9Z56",
   ... )
   >>> tf_gv.copyright = CopyrightInfo(
   ...     release_status="Unrestricted Release",
   ...     conditions_of_use="See USGS terms.",
   ...     reference=reference,
   ... )
   >>> with warnings.catch_warnings(record=True) as caught:
   ...     warnings.simplefilter("always")
   ...     edi_back2 = tf_gv.to_edi()
   ...     print(len(caught), caught[0].category.__name__)
   ...     print(caught[0].message)
   1 DataLossWarning
   current SEG EDI writer has no lossless mapping for EMTF Copyright/Citation metadata

``on_loss`` controls what happens at that same point instead of
warning: ``"raise"`` turns it into an :exc:`EMTFEDIConversionError`
before anything is written, and ``"ignore"`` proceeds silently:

.. code-block:: pycon

   >>> from pycsamt.emtf import EMTFEDIConversionError

   >>> try:
   ...     tf_gv.to_edi(on_loss="raise")
   ... except EMTFEDIConversionError as exc:
   ...     print(type(exc).__name__, exc)
   EMTFEDIConversionError current SEG EDI writer has no lossless mapping for EMTF Copyright/Citation metadata

   >>> with warnings.catch_warnings(record=True) as caught:
   ...     warnings.simplefilter("always")
   ...     tf_gv.to_edi(on_loss="ignore")
   ...     print(len(caught))
   0

``"ignore"`` is for a caller that has already reviewed what will be
dropped and does not want the warning repeated on every call -- not a
default to reach for casually, since the dropped content is otherwise
unrecoverable from the written EDI file.

When EDI Writing Refuses Outright
---------------------------------

Some conditions are not a loss-policy question at all: standard EDI
has no impedance-shaped slot to omit gracefully, so writing a document
with no impedance transfer function fails unconditionally, regardless
of ``on_loss``:

.. code-block:: pycon

   >>> empty_doc = EMTF()
   >>> try:
   ...     empty_doc.to_edi()
   ... except EMTFEDIConversionError as exc:
   ...     print(type(exc).__name__, exc)
   EMTFEDIConversionError standard EDI writing requires an impedance transfer function

The same unconditional check applies to shape and channel identity:
the impedance matrix must be ``(n_period, 2, 2)`` with input channels
exactly ``Hx``/``Hy`` and output channels exactly ``Ex``/``Ey`` (and
tipper, when present, ``(n_period, 1, 2)`` with output ``Hz``) --
standard EDI's ``Z``/``TIP`` blocks have no other layout to write into,
so a document built around a different channel convention (the custom
``"admittance"`` type registered in :doc:`transfer_functions`, for
example) cannot be written as EDI at all, loss policy or not.

Round-Trip Verification
-----------------------

Writing and reloading a clean document (no attached ``copyright``,
which was cleared before this step) preserves every number exactly:

.. code-block:: pycon

   >>> import tempfile
   >>> import numpy as np
   >>> from pycsamt.emtf import write_edi

   >>> tf_gv.copyright = None
   >>> with tempfile.TemporaryDirectory() as tmp:
   ...     out_path = Path(tmp) / "gv100_roundtrip.edi"
   ...     written = write_edi(tf_gv, out_path)
   ...     tf_reloaded = EMTF.from_edi(written)
   ...     print(np.allclose(tf_gv.frequency, tf_reloaded.frequency))
   ...     print(np.allclose(tf_gv.z, tf_reloaded.z))
   ...     print(np.allclose(tf_gv.tipper, tf_reloaded.tipper, equal_nan=True))
   ...     print(tf_gv.site.name, "->", tf_reloaded.site.name)
   True
   True
   True
   Gabbs Valley -> GABBS VALLEY

The frequency grid, impedance, and tipper (``NaN`` gaps included,
via ``equal_nan=True``) all compare exactly equal. ``site.name`` does
not: the EDI writer upper-cases ``HEAD.LOC`` on write, so
``"Gabbs Valley"`` comes back as ``"GABBS VALLEY"``. This is a real,
minor round-trip artifact -- a text-casing convention of the EDI
writer, not a scientific change -- and a reminder that "numerically
exact" and "byte-identical" are different, narrower claims than
"unchanged."

Converting through XML instead of back to EDI keeps everything,
including the ``NaN`` gaps, with no analogous casing surprise (XML
element text is not case-folded the way EDI ``HEAD`` fields are):

.. code-block:: pycon

   >>> with tempfile.TemporaryDirectory() as tmp:
   ...     xml_path = Path(tmp) / "gv100_roundtrip.xml"
   ...     tf_gv.write_xml(xml_path)
   ...     tf_from_xml = EMTF.from_xml(xml_path)
   ...     print(np.allclose(tf_gv.z, tf_from_xml.z))
   ...     print(np.allclose(tf_gv.tipper, tf_from_xml.tipper, equal_nan=True))
   True
   True

:doc:`xml` covers the ``EMTFXMLReader``/``EMTFXMLWriter`` pair used
above in full.

The TFBundle Bridge
-------------------

:func:`bundle_to_emtf` and :func:`emtf_to_bundle` are thin, explicitly
typed wrappers around :meth:`EMTF.from_bundle
<pycsamt.emtf.document.EMTF.from_bundle>` and :meth:`EMTF.to_bundle
<pycsamt.emtf.document.EMTF.to_bundle>`, already covered in full --
including exactly what survives and what does not -- in
:doc:`document`:

.. code-block:: pycon

   >>> from pycsamt.emtf import bundle_to_emtf, emtf_to_bundle

   >>> bundle = emtf_to_bundle(tf_gv)
   >>> back = bundle_to_emtf(bundle)
   >>> print(type(bundle).__name__, type(back).__name__)
   TFBundle EMTF
   >>> print(np.allclose(back.z, tf_gv.z))
   True

   >>> try:
   ...     emtf_to_bundle("not-an-emtf")
   ... except TypeError as exc:
   ...     print(type(exc).__name__, exc)
   TypeError document must be a pycsamt.emtf.EMTF

The only reason these wrappers exist alongside the methods is explicit
type-checking at a stable function boundary -- useful when the caller
only has an arbitrary object and wants a clear error rather than an
``AttributeError`` several calls later.

Choosing the Right Function
---------------------------

.. list-table::
   :header-rows: 1
   :widths: 32 30 38

   * - Need
     - Function
     - Notes
   * - Load a document from EDI
     - :meth:`EMTF.from_edi <pycsamt.emtf.document.EMTF.from_edi>`
     - Prefers EDI SPECTRA recovery when present -- see :doc:`spectra`.
   * - Write a document as EDI
     - :meth:`EMTF.to_edi <pycsamt.emtf.document.EMTF.to_edi>`
     - Loss policy via ``on_loss``; some failures are unconditional.
   * - Write directly to a file path
     - :func:`write_edi`
     - Wraps ``to_edi`` plus the historical ``EDIFile.write``.
   * - Explicit TFBundle interop with a type check
     - :func:`bundle_to_emtf` / :func:`emtf_to_bundle`
     - See :doc:`document` for what the bundle bridge preserves.

Next Steps
----------

* :doc:`document` covers ``EMTF.to_bundle``/``from_bundle`` and the
  ``z_err`` accessor in full;
* :doc:`spectra` covers the richer covariance recovery ``from_edi``
  prefers whenever a file has a ``SPECTRA`` section;
* :doc:`xml` covers ``EMTFXMLReader``/``EMTFXMLWriter`` used in the
  cross-format round trip above.
