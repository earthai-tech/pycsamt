.. _user_guide_metadata_processing_and_quality:

Processing and Quality
======================

.. currentmodule:: pycsamt.metadata

:class:`ProcessingMeta` and :class:`RemoteReferenceMeta` record *how*
a transfer function was estimated -- sign convention, processing
software, and remote-reference bookkeeping -- while
:func:`normalize_sign_convention` reconciles the spelling differences
that otherwise silently flip impedance phases between processing
codes. :class:`DataQuality` and :class:`TransferFunctionQuality` then
answer a different question, *how good* the result is, from two
different sources: ``DataQuality`` is a coverage/SNR assessment pycsamt
computes directly from the response arrays, while
``TransferFunctionQuality`` preserves an *external* quality judgment
-- an archive's EMTF 0-5 rating, an accepted period range, or expert
comments -- that pycsamt did not compute and must not overwrite with
its own opinion. All five classes are :term:`format-neutral metadata`,
read from EDI or EMTF XML through :attr:`~pycsamt.emtf.document.EMTF.processing`
and :attr:`~pycsamt.emtf.document.EMTF.quality`.

Processing Metadata From a Real Station
---------------------------------------

The WILLY L18 EDI line used throughout this section never states a
processing method in its free-text ``>INFO`` block, but its ``MTSECT``
section carries explicit reference-channel coordinates -- and the EDI
adapter treats that as scientifically meaningful even though no INFO
line ever says the word "remote":

.. code-block:: pycon

   >>> from pathlib import Path
   >>> from pycsamt.seg.edi import EDIFile
   >>> from pycsamt.emtf import EMTF

   >>> edi_willy = EDIFile(Path("data/AMT/WILLY_DATA/L18PLT/18-001A.edi"))
   >>> tf_willy = EMTF.from_edi(edi_willy)
   >>> print(tf_willy.processing)
   ProcessingMeta(sign_convention=None, processed_by=None, software=None, remote_reference=RemoteReferenceMeta(reference_type='Remote Reference', site=None, extra=dict(len=2, keys=['edi_rx', 'edi_ry'])), processing_tag=None, run_list=None, extra=dict(len=0, keys=[]))

``remote_reference.reference_type`` reads ``"Remote Reference"`` as a
generated label, not a value copied from the file: ``MTSECT.RX``/``RY``
(``10016.001``/``10017.001`` here) are explicit reference-channel
identifiers, and the adapter infers that their presence *means*
remote-reference processing was used even though the historical
``INFO`` block never named the method outright. The Gabbs Valley EDI
used in :doc:`provenance_and_bibliography` and :doc:`channels_and_orientation`
has none of this -- no sign convention, no software, no remote
reference -- and correctly reports nothing rather than a guess:

.. code-block:: pycon

   >>> edi_gv = EDIFile(Path("data/gv_data/gv_final_edi/gv100.edi"))
   >>> tf_gv = EMTF.from_edi(edi_gv)
   >>> print(tf_gv.processing)
   None

.. note::

   ``gv100.edi`` above is public-domain USGS data; see
   :doc:`provenance_and_bibliography` for the required citation.

A Naive Read Would Mislead
--------------------------

WILLY's raw ``>INFO`` text is free-form field-log prose with no
structured processing keys at all:

.. code-block:: pycon

   >>> info = edi_willy.get_section("info")
   >>> for line in info.info_text[:4]:
   ...     print(line)
                RUN INFORMATION                     STATION 1
       PROCESSED FROM DFT TIME SERIES     STN Number: 23-18-001A
       COMPANY:                           Lat  32:07:131N Long 119:07:438E
       START-UP: 2023/11/09 01:03:00      Site Layout by:

Yet asking the historical EDI section object directly for its
processing metadata returns a plausible-looking, fully populated
answer:

.. code-block:: pycon

   >>> print(info.Processing.signconvention)
   exp(+i ω t)
   >>> print(info.Processing.ProcessingSoftware.name)
   PYCSAMT

Neither value came from the file. ``info.Processing`` is
:class:`pycsamt.seg.property.Processing` -- the historical, EDI-internal
holder described in :doc:`provenance_and_bibliography`'s note on the
two separate metadata models -- and it hands back *its own* default
placeholders (the running package's name and a stock sign convention)
when nothing has ever been assigned, rather than raising or returning
``None``. This is exactly why the format-neutral
:func:`~pycsamt.emtf.converters.edi.edi_to_emtf` adapter never reads
``info.Processing`` directly: it separately confirms each field's key
is *explicitly present* in the raw ``INFO`` text before trusting it,
which is why ``tf_willy.processing.sign_convention`` and ``.software``
both came back ``None`` above despite ``info.Processing`` claiming
otherwise. Prefer ``EMTF.processing`` over the raw EDI section object
whenever the question is "what does this file actually say," not
"what would pycsamt do by default."

Sign Convention Normalization
-----------------------------

:func:`normalize_sign_convention` reconciles the handful of ways
people spell the two Fourier sign conventions without inventing a
convention for text it does not recognize:

.. code-block:: pycon

   >>> from pycsamt.metadata.processing import normalize_sign_convention

   >>> print(normalize_sign_convention("exp(+iwt)"))
   exp(+i ω t)
   >>> print(normalize_sign_convention("EXP(-I OMEGA T)"))
   exp(-i ω t)
   >>> print(normalize_sign_convention("unknown-convention"))
   unknown-convention
   >>> print(normalize_sign_convention(""))
   None

An unrecognized non-empty string is returned unchanged rather than
coerced to a default -- deliberately stricter than the historical SEG
helper it replaces, which silently defaulted anything unrecognized to
the positive convention. Getting the sign wrong flips every impedance
phase by 180°, so a convention pycsamt cannot confidently parse should
surface as-is for a human to check, not disappear into a guess.

Coverage-Based Data Quality
---------------------------

:meth:`~quality.DataQuality.from_site` computes coverage directly from
a station's response arrays -- the fraction of frequencies with a
finite value, per component. Every WILLY L18 station is complete:

.. code-block:: pycon

   >>> from pycsamt.site import Sites
   >>> from pycsamt.metadata.quality import DataQuality

   >>> data_root = Path("data/AMT/WILLY_DATA")
   >>> l18_paths = sorted((data_root / "L18PLT").glob("*.edi"))
   >>> sites = Sites([EDIFile(p) for p in l18_paths])
   >>> dq_willy = DataQuality.from_site(sites[0])
   >>> print(dq_willy.summary())
   Station : 18-001A
     Frequencies : 53  [1.01 – 1.04e+04 Hz]
     Overall     : GOOD
     Zxx      ████████  100%  [good]
     Zxy      ████████  100%  [good]
     Zyx      ████████  100%  [good]
     Zyy      ████████  100%  [good]

Gabbs Valley is not, and its gap lands exactly where the survey's own
release notes say to expect it -- longer periods, tipper specifically:

.. code-block:: pycon

   >>> from pycsamt.site.base import Site

   >>> site_gv = Site(edi_gv)
   >>> dq_gv = DataQuality.from_site(site_gv)
   >>> print(dq_gv.summary())
   Station : gv100
     Frequencies : 48  [0.000488 – 768 Hz]
     Overall     : PARTIAL
     Zxx      ████████  100%  [good]
     Zxy      ████████  100%  [good]
     Zyx      ████████  100%  [good]
     Zyy      ████████  100%  [good]
     Tipper   ██████░░   71%  [partial]

:attr:`~quality.DataQuality.overall` is the *worst* component flag,
so one degraded ``Tipper`` column is enough to pull the whole
station's rating down from ``GOOD`` to ``PARTIAL`` even though all
four impedance components are complete --
:meth:`~quality.QualityFlag.worst` deliberately does not average
components together, since a station usable for impedance-only
inversion but not for induction-vector interpretation should say so
plainly rather than reporting a blended, ambiguous score.

:func:`~quality.assess_collection` and :func:`~quality.quality_dataframe`
scale the same computation to a whole collection:

.. code-block:: pycon

   >>> from pycsamt.metadata.quality import quality_dataframe

   >>> df = quality_dataframe(sites)
   >>> print(df["overall"].value_counts().to_string())
   overall
   good    28

All 28 WILLY stations come back ``GOOD`` -- a real, clean, single-site
survey rather than a contrived example, and worth knowing about this
particular dataset before assuming every station needs individual
review.

Archival Quality Versus Computed Quality
----------------------------------------

:class:`~quality.TransferFunctionQuality` records a judgment pycsamt
did not compute. The Gabbs Valley release notes describe good data
quality over periods of 0.007-2048 s, with some estimates less robust
at the longest periods -- exactly the kind of expert assessment this
class exists to carry, separately from the coverage numbers above:

.. code-block:: pycon

   >>> from pycsamt.metadata.quality import TransferFunctionQuality, QualityComment

   >>> tfq = TransferFunctionQuality(
   ...     rating=4,
   ...     good_from_period=0.007,
   ...     good_to_period=2048.0,
   ...     comments=[
   ...         QualityComment(
   ...             text="Tipper degraded at the longest periods per the USGS release notes.",
   ...             author="pycsamt-docs",
   ...         )
   ...     ],
   ... )
   >>> print(tfq.rating_label)
   good
   >>> print(tfq.has_warning)
   False

A separate warning-flagged record demonstrates the other branch of
:attr:`~quality.TransferFunctionQuality.has_warning`:

.. code-block:: pycon

   >>> flagged = TransferFunctionQuality(
   ...     rating=2,
   ...     warning_flag=1,
   ...     warnings=[QualityComment(text="Static shift suspected at this site.")],
   ... )
   >>> print(flagged.has_warning, flagged.rating_label)
   True serious_issues

Each numeric field is validated against its real-world domain rather
than accepted as an arbitrary number:

.. code-block:: pycon

   >>> try:
   ...     TransferFunctionQuality(rating=7)
   ... except ValueError as exc:
   ...     print(type(exc).__name__, exc)
   ValueError transfer-function quality rating must be 0..5

   >>> try:
   ...     TransferFunctionQuality(good_from_period=100.0, good_to_period=1.0)
   ... except ValueError as exc:
   ...     print(type(exc).__name__, exc)
   ValueError good_from_period must not exceed good_to_period

   >>> try:
   ...     QualityComment(text="   ")
   ... except ValueError as exc:
   ...     print(type(exc).__name__, exc)
   ValueError quality comment text must be non-empty

``rating`` follows the EMTF archive convention: 5 is "great," 0 is
"not_assessed," and everything between is available through
:attr:`~quality.TransferFunctionQuality.rating_label` without
memorizing the numeric scale.

Choosing the Right Object
-------------------------

.. list-table::
   :header-rows: 1
   :widths: 32 30 38

   * - Need
     - Object
     - Notes
   * - Sign convention, software, remote reference
     - :class:`ProcessingMeta` / :class:`RemoteReferenceMeta`
     - Trust ``EMTF.processing``, not a raw EDI section's defaults.
   * - Reconcile sign-convention spelling
     - :func:`normalize_sign_convention`
     - Unrecognized text is preserved, never guessed.
   * - Coverage/SNR computed from response arrays
     - :class:`DataQuality` / :class:`ComponentQuality` / :class:`QualityFlag`
     - Worst-component rule; build with :meth:`~quality.DataQuality.from_site`.
   * - Coverage across a whole collection
     - :func:`~quality.assess_collection` / :func:`~quality.quality_dataframe`
     - Same computation, one row per station.
   * - An external archive's own quality judgment
     - :class:`TransferFunctionQuality` / :class:`QualityComment`
     - Never inferred by pycsamt; attach what a human or archive decided.

Next Steps
----------

* :doc:`site_and_survey` covers station identity and campaign-level
  aggregation;
* :doc:`provenance_and_bibliography` covers who created a document and
  the separate ``pycsamt.seg.property`` metadata model referenced
  above;
* :doc:`instrument` covers the acquisition system behind these
  measurements.
