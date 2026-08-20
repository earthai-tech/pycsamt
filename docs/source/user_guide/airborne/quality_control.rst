.. _user_guide_airborne_quality_control:

Structural Quality Control
==========================

:mod:`pycsamt.airborne.qc` assesses representation quality,
completeness, and metadata consistency for an
:class:`~pycsamt.airborne.AirborneEMDataset` -- how much navigation
has an attached response, whether a fixed ground reference is present
where the technology genuinely needs one, whether a flight line's
frequency axis is even physically valid -- without inventing a
universal geophysical threshold for signal-to-noise ratio, frequency
coverage, or anomaly quality. Those judgments are technology-specific
and belong in :mod:`pycsamt.emtools`, built on top of this common
structural layer rather than duplicated inside it. Every name below
imports from the top-level :mod:`pycsamt.airborne` package.

That restraint shows up directly in how findings are scored.
:func:`~pycsamt.airborne.assess_airborne_qc` classifies a finding as
``"error"`` only when the data are internally inconsistent in a way
nothing downstream could sensibly interpret -- today, exactly one
thing: a non-positive or non-finite frequency axis. A missing EM
record, absent coordinates, or missing reference-station metadata can
matter a great deal, but they describe *incomplete* data, not
*invalid* data, so they are reported at ``"info"``/``"warning"``
severity instead. Knowing that distinction in advance makes the
report's :attr:`~pycsamt.airborne.qc.AirborneQCReport.status` legible
without having to read every finding: ``"error"`` means something is
scientifically broken, ``"warning"``/``"pass"`` means the data are
merely as complete (or not) as they are.

A Compact Inventory
-------------------

:func:`~pycsamt.airborne.inspect_airborne` accepts a dataset, line,
record, or bare ``EMTF`` document and returns one
:class:`~pycsamt.airborne.qc.AirborneInspection`, filling in only the
fields meaningful at that level:

.. code-block:: pycon

   >>> from pycsamt.airborne import ensure_asites, inspect_airborne
   >>> site = ensure_asites("data/ZTEM/gold_springs_nv")[0]
   >>> inspect_airborne(site.record)
   AirborneInspection(object_type='record', technologies=tuple(['ztem']), n_lines=0, n_samples=0, n_records=1, transfer_function_names=tuple(['tipper']), bbox=None, attrs=dict(len=1, keys=[sample_id]))
   >>> inspect_airborne(site.emtf).attrs
   {'product_id': 'gold_springs_nv_L1.GO_L1_001', 'subtype': 'ztem'}

A record's inventory has no line or sample counts to speak of --
``n_lines`` and ``n_samples`` stay at their defaults -- while
``n_records`` and ``transfer_function_names`` are meaningful at every
level, up to and including a whole dataset. Building one from a full
survey means first grouping the flat
:class:`~pycsamt.airborne.site.AirborneSites` view :doc:`site` reads
back into the line/dataset structure
:func:`~pycsamt.airborne.assess_airborne_qc` actually consumes -- a
short, reusable step:

.. code-block:: pycon

   >>> import numpy as np
   >>> from pycsamt.airborne import NavigationTrack, AirborneEMLine, AirborneEMDataset
   >>> def make_line(line_label, members):
   ...     nav = NavigationTrack(
   ...         sample_ids=tuple(s.sample_id for s in members),
   ...         latitude=np.array([s.coords[0] for s in members]),
   ...         longitude=np.array([s.coords[1] for s in members]),
   ...         terrain_elevation=np.array([s.coords[2] for s in members]),
   ...     )
   ...     line = AirborneEMLine(line_id=f"gold_springs_nv_{line_label}", navigation=nav)
   ...     for s in members:
   ...         line.add_emtf(s.sample_id, s.emtf)
   ...     return line
   >>> ztem_sites = ensure_asites("data/ZTEM/gold_springs_nv")
   >>> line_labels = sorted(set(s.name.split("_")[1] for s in ztem_sites))
   >>> lines = {
   ...     f"gold_springs_nv_{lbl}": make_line(
   ...         lbl, [s for s in ztem_sites if s.name.split("_")[1] == lbl],
   ...     )
   ...     for lbl in line_labels
   ... }
   >>> clean_dataset = AirborneEMDataset(name="gold_springs_nv", lines=lines)
   >>> insp = inspect_airborne(clean_dataset)
   >>> insp.object_type, insp.technologies, insp.n_lines, insp.n_samples, insp.n_records
   ('dataset', ('ztem',), 7, 105, 105)

A Clean Survey
--------------

Running :func:`~pycsamt.airborne.assess_airborne_qc` on that same,
untouched real survey establishes what a genuinely complete submission
looks like -- the baseline every later finding below is a deliberate
departure from:

.. code-block:: pycon

   >>> from pycsamt.airborne import assess_airborne_qc
   >>> clean_report = assess_airborne_qc(clean_dataset)
   >>> clean_report.status
   'pass'
   >>> len(clean_report.issues)
   0
   >>> clean_report.metrics["record_coverage_fraction"]
   1.0
   >>> clean_report.metrics["reference_metadata_fraction"]
   1.0
   >>> clean_report.metrics["variance_tf_fraction"]
   1.0
   >>> clean_report.metrics["inverse_signal_covariance_tf_fraction"]
   0.0

Every navigation sample has an attached record, every ZTEM record
carries its fixed-ground reference-station tag, and every primary
transfer function has a variance estimate attached -- this synthetic
survey was built to be complete. The last line is the honest
exception: ``inverse_signal_covariance_tf_fraction`` is ``0.0`` because
none of these records carry a full inverse-signal-covariance matrix,
which is normal for a variance-only delivery and is reported as a
plain metric rather than manufactured into a warning nobody asked for.

Five Real Findings
------------------

Nothing above exercises the interesting code paths, because nothing
is actually wrong with ``gold_springs_nv``. Building a second, smaller
dataset out of the *same* real records -- one navigation sample with
no attached record, one record with its reference-station metadata
stripped, one record with a corrupted frequency axis, and one real
AFMAG record spliced onto an otherwise all-ZTEM line -- exercises five
independent finding codes at once, each from a genuine structural
defect rather than a fabricated one:

.. code-block:: pycon

   >>> import copy
   >>> l1_members = [s for s in ztem_sites if s.name.split("_")[1] == "L1"]
   >>> l2_members = [s for s in ztem_sites if s.name.split("_")[1] == "L2"]
   >>> line1 = make_line("L1", l1_members)
   >>> _ = line1.records.pop("GO_L1_008")
   >>> line2 = make_line("L2", l2_members)
   >>> no_ref_doc = copy.deepcopy(line2.get_record("GO_L2_003").emtf)
   >>> no_ref_doc.processing.remote_reference = None
   >>> _ = line2.add_emtf("GO_L2_003", no_ref_doc, replace=True)
   >>> bad_freq_doc = copy.deepcopy(line2.get_record("GO_L2_010").emtf)
   >>> bad_freq_doc.periods = bad_freq_doc.periods.copy()
   >>> bad_freq_doc.periods[0] = -1.0
   >>> _ = line2.add_emtf("GO_L2_010", bad_freq_doc, replace=True)
   >>> afmag_site = ensure_asites("data/AFMAG/abitibi_on")[0]
   >>> line2.navigation = NavigationTrack(
   ...     sample_ids=tuple(line2.navigation.sample_ids) + ("AB_001",),
   ...     latitude=np.append(line2.navigation.latitude, afmag_site.coords[0]),
   ...     longitude=np.append(line2.navigation.longitude, afmag_site.coords[1]),
   ...     terrain_elevation=np.append(
   ...         line2.navigation.terrain_elevation, afmag_site.coords[2],
   ...     ),
   ... )
   >>> _ = line2.add_emtf("AB_001", afmag_site.emtf)
   >>> messy_dataset = AirborneEMDataset(
   ...     name="gold_springs_nv_demo",
   ...     lines={"gold_springs_nv_L1": line1, "gold_springs_nv_L2": line2},
   ... )
   >>> messy_report = assess_airborne_qc(messy_dataset)
   >>> messy_report.status
   'error'
   >>> for issue in messy_report.issues:
   ...     print(f"[{issue.severity:<7}] {issue.code:<25} line={issue.line_id!r} sample={issue.sample_id!r}")
   [info   ] missing_em_records        line='gold_springs_nv_L1' sample=None
   [warning] mixed_line_technology     line='gold_springs_nv_L2' sample=None
   [warning] missing_reference_station line='gold_springs_nv_L2' sample='GO_L2_003'
   [error  ] invalid_frequency_axis    line='gold_springs_nv_L2' sample='GO_L2_010'
   [info   ] mixed_dataset_technology  line=None sample=None

Each finding traces back to exactly one deliberate change. Removing
``GO_L1_008``'s record but leaving its navigation sample in place --
``line.records.pop(...)`` rather than rebuilding the line without it
-- produces ``missing_em_records`` scoped to line ``L1`` alone. Setting
``processing.remote_reference`` to ``None`` on a copy of
``GO_L2_003``'s document reproduces exactly the check
:doc:`registry_and_io` describes: ZTEM is one of the technologies
whose ``reference_required`` contract is ``True``, so a ZTEM record
with no reference tag is flagged, sample-scoped. Corrupting
``GO_L2_010``'s frequency axis with one negative period is the only
``"error"``-severity finding here, exactly as the module's severity
philosophy above promises -- everything else is *incomplete*, this
alone is *invalid*. Splicing a real AFMAG record onto a ZTEM line's
navigation produces two distinct findings from one change: line ``L2``
itself is flagged as ``mixed_line_technology``, and because that
mixture is now visible at dataset scope too,
:func:`~pycsamt.airborne.identify_airborne_technologies` -- the same
function :doc:`registry_and_io` uses directly -- also drives
``mixed_dataset_technology`` on the whole report.

:attr:`~pycsamt.airborne.qc.AirborneQCReport.errors` and
:attr:`~pycsamt.airborne.qc.AirborneQCReport.warnings` filter
:attr:`~pycsamt.airborne.qc.AirborneQCReport.issues` by severity, and
:attr:`~pycsamt.airborne.qc.AirborneQCReport.line_metrics` keeps the
same coverage-fraction metrics :attr:`~pycsamt.airborne.qc.AirborneQCReport.metrics`
reports dataset-wide, but broken out per line:

.. code-block:: pycon

   >>> len(messy_report.errors), len(messy_report.warnings)
   (1, 2)
   >>> messy_report.line_metrics["gold_springs_nv_L1"]["record_coverage_fraction"]
   0.9333333333333333
   >>> messy_report.line_metrics["gold_springs_nv_L2"]["transfer_function_names"]
   ('afmag_tilt', 'tipper')

``L1``'s coverage fraction, ``14/15``, isolates exactly which line
lost a record without having to search the dataset-wide fraction for
it, and ``L2``'s transfer-function names show the AFMAG splice from a
completely different angle than the issues list did: this line now
carries two distinct transfer-function types where a genuine ZTEM
line would carry exactly one. Neither of these two examples is
exhaustive -- a line with no finite navigation coordinates at all
raises a sixth code, ``missing_navigation_coordinates``, not
reproduced here because every sample in both surveys used above
already has a real position.
