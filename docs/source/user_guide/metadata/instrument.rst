.. _user_guide_metadata_instrument:

Instrument Metadata
===================

.. currentmodule:: pycsamt.metadata.instrument

:class:`InstrumentMeta` describes the acquisition system that recorded
a station -- the data logger, its sensors, and the acquisition
software -- in a form that round-trips to and from EDI ``>HEAD``
fields as well as JSON/YAML. :class:`SensorSpec` describes one sensor
channel (electrode pair or magnetic coil) within that system. Because
the same handful of commercial systems -- Phoenix, Metronix, LEMI,
Zonge -- recur across most MT/AMT/CSAMT surveys, :data:`KNOWN_SYSTEMS`
ships ready-made presets accessible through :func:`known_system` and
:func:`list_presets`, so a station's instrument metadata does not need
to be typed out by hand for a standard configuration.

Unlike :class:`~pycsamt.metadata.site.SiteMeta` or
:class:`~pycsamt.metadata.provenance.ProvenanceMeta`,
``InstrumentMeta`` and ``SensorSpec`` are plain dataclasses -- they do
not inherit :class:`~pycsamt.api.property.PyCSAMTObject`, so their
compact display comes from a hand-written ``__repr__``/``__str__``
rather than the shared field-capping convention documented in
:doc:`site_and_survey`.

Sensor Coverage
---------------

A :class:`SensorSpec` records the transducer type, model, and rated
frequency band. :meth:`~SensorSpec.covers` is a direct, practical use
of that band: does this sensor's rating actually span the frequencies
a survey needs?

.. code-block:: pycon

   >>> from pycsamt.metadata.instrument import SensorSpec

   >>> coil = SensorSpec("induction_coil", "MTC-150H", (1e-4, 2e3))
   >>> print(coil)
   <SensorSpec 'MTC-150H' [0.0001–2e+03 Hz]>

The Gabbs Valley survey documented in :doc:`provenance_and_bibliography`
reports good data quality over roughly 0.0005-770 Hz. This coil's
rated band comfortably covers that range, but not an arbitrarily low
frequency:

.. code-block:: pycon

   >>> print(coil.covers(0.0005), coil.covers(770.0))
   True True
   >>> print(coil.covers(0.00001))
   False

A ``frequency_range`` of ``None`` means "unrated/unspecified", which
:meth:`~SensorSpec.covers` treats as universal coverage (always
``True``) rather than as "unknown, therefore fails" -- appropriate for
an electrode, which has no meaningful frequency cutoff of its own, but
worth remembering if a magnetic sensor's range is left unset by
mistake: an omitted rating silently stops acting as a check at all.

Built-In Presets
----------------

:data:`KNOWN_SYSTEMS` currently covers eight common field systems.
:func:`list_presets` returns their keys, and :func:`known_system`
builds a ready-made :class:`InstrumentMeta` from one:

.. code-block:: pycon

   >>> from pycsamt.metadata.instrument import known_system, list_presets

   >>> print(list_presets())
   ['generic_fluxgate', 'geometrics_stratagem', 'lemi_424', 'metronix_adu07', 'metronix_adu08', 'phoenix_mtx', 'phoenix_v8', 'zonge_gdp32']

   >>> phoenix = known_system("phoenix_v8")
   >>> print(phoenix.summary())
   System       : Phoenix V8
   Serial       : —
   Software     : SSMT2000
   Mag. sensor  : <SensorSpec 'MTC-150H' [0.0001–2e+03 Hz]>
   Elec. sensor : <SensorSpec 'Pb-PbCl2 porous pot' [full band]>
   Notes        : Standard Phoenix V8 MT/AMT system

   >>> try:
   ...     known_system("not_a_real_system")
   ... except KeyError as exc:
   ...     print(type(exc).__name__, exc)
   KeyError "Unknown preset 'not_a_real_system'. Available: generic_fluxgate, geometrics_stratagem, lemi_424, metronix_adu07, metronix_adu08, phoenix_mtx, phoenix_v8, zonge_gdp32"

``known_system`` is a thin shortcut for
:meth:`InstrumentMeta.from_preset`; the preset key is normalized
(spaces and hyphens both become underscores, matched case-insensitively)
before the lookup, so ``"Phoenix V8"`` and ``"phoenix-v8"`` resolve the
same way as ``"phoenix_v8"``.

Presets are a starting point, not a claim about any specific station.
Neither of the real files used elsewhere in this section names its
acquisition hardware in a way pycsamt can recognize automatically --
see the ``ACQBY`` heuristic below -- so treat a preset as "the closest
documented system," to be confirmed against the field log, and set
``serial`` explicitly once the actual unit is known:

.. code-block:: pycon

   >>> from pycsamt.metadata.instrument import InstrumentMeta

   >>> inst = known_system("phoenix_v8")
   >>> inst.serial = "V8-20473"
   >>> print(inst.label)
   Phoenix V8 / V8-20473

EDI ``>HEAD`` Round-Trip
------------------------

:meth:`~InstrumentMeta.to_head_fields` maps ``system``/``serial`` to
``acqby`` and ``software_version`` to ``progvers`` -- falling back to
the running pycsamt version when no software version was recorded, so
the field is never left blank by omission:

.. code-block:: pycon

   >>> print(inst.to_head_fields())
   {'acqby': 'Phoenix V8 / V8-20473', 'progvers': 'pyCSAMT 2.3.1'}

The reverse direction, :meth:`~InstrumentMeta.from_head`, is a
best-effort heuristic, not a guaranteed inverse: it splits ``ACQBY`` on
the first ``"/"`` and calls the left half ``system``. That works
exactly when a file's ``ACQBY`` was written as ``"System / Serial"``
in the first place. Neither real EDI used throughout this section was:

.. note::

   ``gv100.edi`` below is public-domain USGS data; see
   :doc:`provenance_and_bibliography` for the required citation.

.. code-block:: pycon

   >>> from pathlib import Path
   >>> from pycsamt.seg.edi import EDIFile

   >>> edi_willy = EDIFile(Path("data/AMT/WILLY_DATA/L18PLT/18-001A.edi"))
   >>> head_willy = edi_willy.get_section("head")
   >>> print(repr(head_willy.acqby), repr(head_willy.progvers))
   None 'MTPROC V1.0.7'
   >>> inst_willy = InstrumentMeta.from_head(head_willy)
   >>> print(inst_willy)
   <InstrumentMeta system='' serial=None>

   >>> edi_gv = EDIFile(Path("data/gv_data/gv_final_edi/gv100.edi"))
   >>> head_gv = edi_gv.get_section("head")
   >>> print(repr(head_gv.acqby), repr(head_gv.progvers))
   'Jared Peacock' '1.1.5'
   >>> inst_gv = InstrumentMeta.from_head(head_gv)
   >>> print(inst_gv)
   <InstrumentMeta system='Jared Peacock' serial=None>

WILLY's ``ACQBY`` was empty, so ``system`` comes back empty too --
harmless. Gabbs Valley's ``ACQBY`` holds the *person* who acquired the
data, not an instrument string, and ``from_head`` has no way to tell
the difference: it confidently reports the acquiring scientist's name
as the ``system``. Neither result is a parsing failure; ``from_head``
faithfully reflects what the heuristic can extract from ``ACQBY``
alone, and that heuristic simply does not hold for every real file.
Treat its output as a suggestion to verify, not as ground truth, and
prefer :func:`known_system` or a manual :class:`InstrumentMeta` when
the real hardware is documented elsewhere (a field log, the survey's
own metadata release, or -- as with the presets above -- domain
knowledge of what that project typically ran).

Serialization
-------------

:class:`InstrumentMeta` round-trips through JSON (and YAML, when
PyYAML is installed) the same way
:class:`~pycsamt.metadata.survey.SurveyMeta` does in
:doc:`site_and_survey`, nesting each :class:`SensorSpec` as a plain
dict:

.. code-block:: pycon

   >>> manual = InstrumentMeta(
   ...     system="Phoenix V8",
   ...     serial="V8-20473",
   ...     magnetic_sensor=SensorSpec(
   ...         sensor_type="induction_coil",
   ...         model="MTC-150H",
   ...         frequency_range=(1e-4, 2e3),
   ...     ),
   ...     electric_sensor=SensorSpec(
   ...         sensor_type="electrode",
   ...         model="Pb-PbCl2 porous pot",
   ...     ),
   ... )
   >>> text = manual.to_json()
   >>> reloaded = InstrumentMeta.from_json(text)
   >>> print(reloaded.magnetic_sensor)
   <SensorSpec 'MTC-150H' [0.0001–2e+03 Hz]>
   >>> print(reloaded == manual)
   True

``save``/``load`` wrap the same JSON/YAML methods and pick the format
from the file extension, so a preset built in memory can be written
once (``inst.save("phoenix_v8.json")``) and reloaded later without
re-specifying every sensor field.

Next Steps
----------

* :doc:`site_and_survey` covers station identity and campaign-level
  aggregation;
* :doc:`channels_and_orientation` covers the physical channel geometry
  these sensors are wired into;
* :doc:`processing_and_quality` covers how the resulting transfer
  function was estimated and how good the result is.
