.. _user_guide_airborne_registry_and_io:

Technologies, Formats, and Native I/O
=====================================

:mod:`pycsamt.airborne.registry` separates two concepts that a single
"format" tag would blur together: a **technology** describes
scientific semantics -- what physical quantity a response represents,
whether a fixed ground reference is required, which transfer-function
names are unique enough to infer it from -- while a **format**
describes one concrete way that technology's response is delivered on
disk. :mod:`pycsamt.airborne.io` then dispatches
:func:`~pycsamt.airborne.read_airborne`/
:func:`~pycsamt.airborne.write_airborne` through the format half of
that registry. Every name used on this page is imported straight from
the top-level :mod:`pycsamt.airborne` package -- ``from pycsamt.airborne
import ...`` -- since the registry and I/O dispatcher are both part of
its public surface; only the tiny illustrative ZTEM line built partway
through needs a submodule import, because ``build_ztem_line`` is a
technology-specific constructor that :mod:`pycsamt.airborne` itself
deliberately does not re-export.

This registry is not the one that reads the ``data/ZTEM/``,
``data/AFMAG/``, and ``data/mobileMT/`` sample surveys used elsewhere
in this guide -- that path is
:func:`~pycsamt.airborne.site.ensure_asites`, covered in :doc:`site`,
reading EMTF-XML directly with no format lookup at all. The registry
on this page exists for a different, still-open problem: a genuine
*native* vendor delivery (whatever binary or proprietary text format a
MobileMT, ZTEM, or AFMAG system actually writes in the field), which
no vendor has yet supplied a verified sample of.

Built-In Technologies
---------------------

Four technology contracts are registered the moment
:mod:`pycsamt.airborne` is imported --
:func:`~pycsamt.airborne.list_airborne_technologies` lists them in
registration order:

.. code-block:: pycon

   >>> from pycsamt.airborne import list_airborne_technologies, get_airborne_technology
   >>> for t in list_airborne_technologies():
   ...     print(
   ...         t.name, "|", t.label, "| family=", t.family,
   ...         "| aliases=", t.aliases,
   ...         "| reference_required=", t.reference_required,
   ...         "| infer_from_tf=", t.infer_from_tf,
   ...     )
   mobilemt | MobileMT | family= natural_field_airborne_em | aliases= ('mobile_mt',) | reference_required= True | infer_from_tf= True
   ztem | ZTEM | family= natural_field_airborne_em | aliases= ('z_tem',) | reference_required= True | infer_from_tf= False
   afmag | AFMAG (original comparator) | family= afmag | aliases= ('original_afmag', 'comparator_afmag') | reference_required= False | infer_from_tf= True
   airmt | AirMt / tensor AFMAG | family= afmag | aliases= ('tensor_afmag', 'afmag_tensor') | reference_required= True | infer_from_tf= False

Two fields decide how much a technology can be inferred rather than
declared. ``reference_required`` records a scientific fact --
:mod:`~pycsamt.airborne.qc` treats a missing fixed ground reference as
a QC issue only for technologies where one is genuinely needed, which
is every technology here except the original comparator. ``afmag`` and
``mobilemt`` allow the stronger ``infer_from_tf=True``, meaning
:func:`~pycsamt.airborne.identify_airborne_technologies` may recognize
them purely from a matching transfer-function name (``afmag_tilt``,
``mobilemt_admittance``); ``ztem`` and ``airmt`` deliberately do not,
because their signature transfer-function names -- a standard tipper,
an interstation transfer function -- are not unique enough by
themselves to be trusted as proof.

Lookup is by canonical name or by any registered alias, and an
unregistered name returns ``None`` rather than raising, since callers
use this for best-effort inference over untrusted metadata:

.. code-block:: pycon

   >>> get_airborne_technology("mobile_mt") is get_airborne_technology("mobilemt")
   True
   >>> get_airborne_technology("nope") is None
   True

Identifying A Technology
------------------------

:func:`~pycsamt.airborne.identify_airborne_technologies` inspects an
already-built object -- an ``EMTF`` document, record, line, or dataset
-- and returns every technology it can safely recognize.
:func:`~pycsamt.airborne.detect_airborne_technology` collapses that to
one name, or raises when more than one is found. Running it on one
real ``EMTF`` document from each of the four sample surveys makes a
naming distinction worth being deliberate about:

.. code-block:: pycon

   >>> from pycsamt.airborne import ensure_asites, identify_airborne_technologies
   >>> docs = {
   ...     "abitibi_on": ensure_asites("data/AFMAG/abitibi_on")[0].emtf,
   ...     "yulong_belt_cn": ensure_asites("data/AFMAG/yulong_belt_cn")[0].emtf,
   ...     "gold_springs_nv": ensure_asites("data/ZTEM/gold_springs_nv")[0].emtf,
   ...     "flammefjeld_greenland": ensure_asites("data/mobileMT/flammefjeld_greenland")[0].emtf,
   ... }
   >>> for survey, doc in docs.items():
   ...     print(survey, "subtype=", doc.subtype, "identify=", identify_airborne_technologies(doc))
   abitibi_on subtype= afmag_original identify= ('afmag',)
   yulong_belt_cn subtype= afmag_airmt identify= ('airmt',)
   gold_springs_nv subtype= ztem identify= ('ztem',)
   flammefjeld_greenland subtype= mobilemt identify= ('mobilemt',)

``EMTF.subtype`` -- the same string :doc:`site`'s
:attr:`~pycsamt.airborne.site.AirborneSite.technology` reads -- and
the registry's canonical technology name agree for ``ztem`` and
``mobilemt``, but not for either AFMAG generation: ``subtype``
carries ``"afmag_original"``/``"afmag_airmt"``, while the registry
canonicalizes those to ``"afmag"``/``"airmt"`` (an internal mapping
inside :mod:`pycsamt.airborne.registry`). Code that compares
``AirborneSite.technology`` against a name returned by
``identify_airborne_technologies``/``detect_airborne_technology``
needs to know both spellings exist for AFMAG; ``ztem`` and
``mobilemt`` happen to make that mistake invisible, which is exactly
why it is easy to miss until an AFMAG survey hits it.

Identification does not just trust one tag blindly -- it walks
``attrs["technology"]``, ``EMTF.subtype``, and per-transfer-function
tags, and *accumulates* every match rather than stopping at the
first. Deliberately mislabelling a real ZTEM document's owning record
demonstrates what a genuine conflict looks like:

.. code-block:: pycon

   >>> from pycsamt.airborne import (
   ...     AirborneEMRecord, AirborneTechnologyAmbiguityError, detect_airborne_technology,
   ... )
   >>> mislabeled = AirborneEMRecord(
   ...     sample_id="S00", emtf=docs["gold_springs_nv"], attrs={"technology": "mobilemt"},
   ... )
   >>> identify_airborne_technologies(mislabeled)
   ('mobilemt', 'ztem')
   >>> try:
   ...     detect_airborne_technology(mislabeled)
   ... except AirborneTechnologyAmbiguityError as exc:
   ...     print("AirborneTechnologyAmbiguityError:", exc)
   AirborneTechnologyAmbiguityError: multiple airborne technologies are present: mobilemt, ztem
   >>> detect_airborne_technology(mislabeled, strict=False) is None
   True

The explicit ``attrs`` tag (``"mobilemt"``) and the document's own
``subtype``-derived signal (``"ztem"``) are both genuine, so
``identify_airborne_technologies`` reports both instead of silently
preferring one; ``detect_airborne_technology`` turns that into a hard
error by default, and only swallows it to ``None`` when the caller has
explicitly opted into ``strict=False``. This is precisely the kind of
inconsistency :doc:`quality_control` checks for structurally, on a
whole dataset rather than one hand-built record.

No Native Reader Yet
--------------------

The *format* half of the registry starts empty, because no vendor has
supplied pyCSAMT with a verified native delivery to build a reader
against -- registering a decoder from a published system description
alone risks silently guessing a schema wrong:

.. code-block:: pycon

   >>> from pycsamt.airborne import (
   ...     list_airborne_formats, available_airborne_readers, available_airborne_writers,
   ...     detect_airborne_format, read_airborne, AirborneIOError,
   ... )
   >>> list_airborne_formats()
   ()
   >>> available_airborne_readers(), available_airborne_writers()
   ((), ())
   >>> detect_airborne_format("data/ZTEM/gold_springs_nv/GO_L1_001.xml") is None
   True
   >>> try:
   ...     read_airborne("data/ZTEM/gold_springs_nv/GO_L1_001.xml")
   ... except AirborneIOError as exc:
   ...     print("AirborneIOError:", exc)
   AirborneIOError: no registered native airborne format recognized the source; register a verified reader after obtaining a representative delivery file

Note that this is *not* the same failure a corrupt or unsupported file
would raise elsewhere in pyCSAMT -- ``AirborneIOError`` is a
``RuntimeError``, not a ``ValueError``, precisely because every
failure reachable here today means "no reader has been registered
yet," a capability gap rather than bad input. Passing that same
``.xml`` path to :func:`~pycsamt.airborne.site.ensure_asites` instead
reads it immediately, since EMTF-XML flows through an entirely
separate path that was never routed through this registry in the
first place.

Registering A Format
--------------------

Once a verified sample exists, a format is registered by pairing a
reader/writer with a way to recognize it -- a ``detector`` callable, a
tuple of file extensions, or both. The example below stands in for
that real future reader with a small in-memory ZTEM line (three
stations, two frequencies, built the same way :doc:`data_model`
constructs one from scratch) so the mechanics can be shown without a
real vendor file to point at:

.. code-block:: pycon

   >>> import numpy as np
   >>> from pycsamt.airborne import NavigationTrack, AirborneEMDataset
   >>> from pycsamt.airborne.ztem import build_ztem_line, ZTEMSystemSpec
   >>> nav = NavigationTrack(
   ...     sample_ids=("S00", "S01", "S02"),
   ...     easting=np.array([0.0, 50.0, 100.0]),
   ...     northing=np.zeros(3),
   ... )
   >>> tip = np.zeros((3, 2, 2), dtype=complex)
   >>> tip[:, :, 0] = 0.01 + 0.002j
   >>> tip[:, :, 1] = 0.003 - 0.001j
   >>> demo_line = build_ztem_line(
   ...     "DEMO01", nav, tip, frequency=np.array([90.0, 180.0]),
   ...     system_spec=ZTEMSystemSpec(),
   ... )
   >>> demo_dataset = AirborneEMDataset(name="demo_survey", lines={"DEMO01": demo_line})

A detector gets first refusal at recognizing a source; an extension is
only consulted as a fallback hint, and only once no detector claims
it. Registering two formats under the same technology -- one that
recognizes a header string, one that only knows the extension --
makes that ordering concrete:

.. code-block:: pycon

   >>> from pathlib import Path
   >>> from pycsamt.airborne import AirborneFormatDefinition, register_airborne_format
   >>> def _has_pycsamt_header(source):
   ...     try:
   ...         with open(source) as f:
   ...             return f.readline().strip() == "PYCSAMT-DEMO-ZTEM-V1"
   ...     except OSError:
   ...         return False
   >>> def _demo_reader(source, **kwargs):
   ...     return demo_dataset
   >>> def _demo_writer(dataset, target, **kwargs):
   ...     with open(target, "w") as f:
   ...         f.write("PYCSAMT-DEMO-ZTEM-V1\n")
   ...         f.write(f"n_lines={dataset.n_lines} n_records={dataset.n_records}\n")
   ...     return Path(target)
   >>> _ = register_airborne_format(AirborneFormatDefinition(
   ...     name="demo_ztem_detected", technology="ztem",
   ...     reader=_demo_reader, writer=_demo_writer, detector=_has_pycsamt_header,
   ...     description="Detector-recognized illustrative format (has the header).",
   ... ))
   >>> _ = register_airborne_format(AirborneFormatDefinition(
   ...     name="demo_ztem_by_ext", technology="ztem",
   ...     reader=_demo_reader, extensions=(".demo",),
   ...     description="Extension-only fallback illustrative format.",
   ... ))
   >>> [(f.name, f.readable, f.writable, f.extensions) for f in list_airborne_formats(technology="ztem")]
   [('demo_ztem_detected', True, True, ()), ('demo_ztem_by_ext', True, False, ('.demo',))]
   >>> available_airborne_readers(technology="ztem")
   ('demo_ztem_detected', 'demo_ztem_by_ext')
   >>> available_airborne_writers(technology="ztem")
   ('demo_ztem_detected',)

Two files sharing the exact same ``.demo`` extension resolve to two
different formats depending only on their content, confirming the
detector really is checked first rather than the extension winning by
registration order:

.. code-block:: pycon

   >>> import tempfile
   >>> from pycsamt.airborne import write_airborne
   >>> with tempfile.TemporaryDirectory() as tmp:
   ...     tmp = Path(tmp)
   ...     with_header = tmp / "survey_a.demo"
   ...     without_header = tmp / "survey_b.demo"
   ...     _ = with_header.write_text("PYCSAMT-DEMO-ZTEM-V1\nfake bytes\n")
   ...     _ = without_header.write_text("not a pycsamt demo file\n")
   ...     detected_with = detect_airborne_format(with_header)
   ...     detected_without = detect_airborne_format(without_header)
   ...     result = read_airborne(with_header)
   ...     out_path = tmp / "roundtrip.demo"
   ...     written = write_airborne(demo_dataset, out_path, format="demo_ztem_detected")
   ...     written_text = written.read_text()
   >>> detected_with, detected_without
   ('demo_ztem_detected', 'demo_ztem_by_ext')
   >>> result is demo_dataset
   True
   >>> print(written_text)
   PYCSAMT-DEMO-ZTEM-V1
   n_lines=1 n_records=3
   <BLANKLINE>

``with_header`` and ``without_header`` differ only in their first
line, not their extension, yet resolve to different formats -- exactly
the detector-before-extension priority the module docstring promises.
Registering a second extension-only format that collides with the
first, with no detector to break the tie, turns that same lookup into
an error instead of a silent guess:

.. code-block:: pycon

   >>> from pycsamt.airborne import AirborneFormatDetectionError
   >>> _ = register_airborne_format(AirborneFormatDefinition(
   ...     name="demo_ztem_alt_ext", technology="ztem", extensions=(".demo",),
   ...     description="A second, deliberately colliding extension-only format.",
   ... ))
   >>> try:
   ...     detect_airborne_format("another_survey.demo")
   ... except AirborneFormatDetectionError as exc:
   ...     print("AirborneFormatDetectionError:", exc)
   AirborneFormatDetectionError: file extension is ambiguous across airborne formats: demo_ztem_by_ext, demo_ztem_alt_ext

An ambiguous extension is a real failure mode, not a corner case
invented for this page: it is exactly what would happen if two
technologies' native deliveries genuinely shared a generic extension
like ``.csv``, which is one more reason a verified ``detector`` is
worth writing alongside a reader rather than relying on extensions
alone once real formats are registered here.
