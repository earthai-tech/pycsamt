.. _user_guide_iot_provenance:

Provenance and Reproducibility
==============================

Provenance records explain what happened during acquisition and make a
:term:`field session` reproducible after the instruments have left the
site. In an IoT workflow, that :term:`provenance manifest` should
include station occupation metadata, devices, :term:`quality control`
decisions, raw-file integrity hashes, processing steps, and the
environment used to create the manifest itself.

This page uses one real file from the repository demo data,
``data/AMT/WILLY_DATA/L18PLT/18-001A.edi``, for the raw-file hash. The
IoT :term:`telemetry packet`\ s are synthetic because provenance is
documenting acquisition events around the data, not reprocessing the
:term:`EDI` file itself.

Hash A Raw File
---------------

Use :func:`~pycsamt.iot.provenance.hash_raw_file` to create an integrity
record. It streams the file in fixed-size chunks through SHA-256 rather
than loading it whole, so :math:`\mathrm{digest} =
\mathrm{SHA256}(\mathrm{bytes})` comes out identically however the file
is read — in one pass or in chunks, on this machine or another. The full
digest is stored in manifests; the example prints only a prefix to keep
the documentation readable.

.. code-block:: pycon

   >>> import json
   >>> from pathlib import Path
   >>> from pycsamt.iot import hash_raw_file
   >>> raw_file = Path("data/AMT/WILLY_DATA/L18PLT/18-001A.edi")
   >>> raw_hash = hash_raw_file(str(raw_file))
   >>> print(json.dumps({
   ...     "name": raw_hash["name"],
   ...     "bytes": raw_hash["bytes"],
   ...     "algo": raw_hash["algo"],
   ...     "digest_prefix": raw_hash["digest"][:16],
   ... }, indent=2))
   {
     "name": "18-001A.edi",
     "bytes": 22758,
     "algo": "sha256",
     "digest_prefix": "a7bc0c7b290c1847"
   }

Log QC Decisions
----------------

Use :func:`~pycsamt.iot.provenance.log_qc_decision` for normalized audit
records. It standardizes station IDs, channel names, decisions, reasons,
timestamps, and optional acquisition windows into the same shape every
time, which matters later: these records are exactly what feeds the
:term:`hash chain` built further down this page.

.. code-block:: pycon

   >>> from pycsamt.iot import log_qc_decision
   >>> qc_reject = log_qc_decision(
   ...     "001A", "reject", channel="ey",
   ...     reasons=["finite_coverage_below_threshold"],
   ...     operator="pyCSAMT documentation", timestamp=1_700_000_120.0,
   ...     window=(1_700_000_060.0, 1_700_000_120.0),
   ... )
   >>> print(json.dumps(qc_reject, indent=2))
   {
     "station": "001A",
     "channel": "ey",
     "decision": "reject",
     "reasons": [
       "finite_coverage_below_threshold"
     ],
     "operator": "pyCSAMT documentation",
     "logged_utc": "2023-11-14T22:15:20+00:00",
     "window": [
       1700000060.0,
       1700000120.0
     ]
   }

Building A Session Manifest
---------------------------

The easiest route is to build a :term:`field session`, add devices,
stations, and packets, then call
:meth:`~pycsamt.iot.session.FieldSession.to_manifest`. The manifest can
still be enriched afterward with raw-file hashes, extra QC decisions,
and processing steps. ``created_utc`` and each decision's ``logged_utc``
are overridden below to fixed timestamps purely so this page's output is
reproducible; a real deployment would leave them at their real wall-clock
values.

.. code-block:: pycon

   >>> from pycsamt.iot import (
   ...     DeviceConfig, FieldSession, MonitoringConfig, StationConfig, TelemetryPacket,
   ... )
   >>> device = DeviceConfig(
   ...     "l18-node-01", station="001A", protocol="file", sample_rate_hz=512.0,
   ...     channels=["ex", "ey", "hx", "hy"], role="recorder",
   ...     metadata={
   ...         "instrument_serial": "DOC-REC-001", "firmware": "demo-2.0",
   ...         "source_edi": raw_file.as_posix(),
   ...     },
   ... )
   >>> station = StationConfig(
   ...     "001A", profile="L18", position_m=0.0, channels=["ex", "ey", "hx", "hy"],
   ...     dipole_length_m=50.0, ex_azimuth_deg=90.0, ey_azimuth_deg=0.0,
   ...     device_ids=[device.device_id], operator="pyCSAMT documentation",
   ...     notes="Station id and raw EDI file are from the repository demo data.",
   ... )
   >>> session = FieldSession(
   ...     "WILLY-L18-PROVENANCE-DEMO", devices=[device], stations=[station],
   ...     method="amt", operator="pyCSAMT documentation",
   ...     monitoring_config=MonitoringConfig(
   ...         method="amt", required_channels=["ex", "ey", "hx", "hy"],
   ...         frequency_band_hz=(1.0, 1000.0),
   ...     ),
   ...     metadata={"real_data_path": "data/AMT/WILLY_DATA/L18PLT"},
   ... )
   >>> for index, (accepted, decision, reasons) in enumerate(
   ...     [(True, "accept", []), (False, "reject", ["finite_coverage_below_threshold"])]
   ... ):
   ...     _ = session.add_packet(
   ...         TelemetryPacket.from_device(
   ...             device, timestamp=1_700_000_000.0 + 60.0 * index,
   ...             survey_id=session.survey_id, kind="qc",
   ...             payload={
   ...                 "method": "amt", "station": "001A",
   ...                 "channels": ["ex", "ey", "hx", "hy"],
   ...                 "frequency_band_hz": [1.0, 1000.0], "accepted": accepted,
   ...                 "decision": decision, "reasons": reasons,
   ...                 "battery_v": 12.4 - 0.1 * index,
   ...             },
   ...         )
   ...     )
   >>> manifest = session.to_manifest()
   >>> manifest.created_utc = "2023-11-14T22:13:20+00:00"
   >>> manifest.environment = {
   ...     "python": "3.10", "platform": "documentation", "machine": "example",
   ... }
   >>> for record in manifest.records:
   ...     for item in record.qc_decisions:
   ...         item["logged_utc"] = "2023-11-14T22:13:20+00:00"
   >>> _ = manifest.add_file(str(raw_file))
   >>> _ = manifest.add_qc_decision(
   ...     station="001A", decision="reject", channel="ey",
   ...     reasons=["finite_coverage_below_threshold"],
   ...     operator="pyCSAMT documentation", timestamp=1_700_000_120.0,
   ...     window=(1_700_000_060.0, 1_700_000_120.0),
   ... )
   >>> manifest.add_processing_step("loaded station metadata from L18 demo line")
   >>> manifest.add_processing_step("recorded synthetic IoT QC packet decisions")
   >>> manifest_dict = manifest.as_dict()
   >>> print(json.dumps({
   ...     "survey_id": manifest_dict["survey_id"],
   ...     "method": manifest_dict["method"],
   ...     "operator": manifest_dict["operator"],
   ...     "tool_version": manifest_dict["tool_version"],
   ...     "n_records": len(manifest_dict["records"]),
   ...     "n_devices": len(manifest_dict["devices"]),
   ...     "n_files": len(manifest_dict["files"]),
   ...     "n_qc_decisions": len(manifest_dict["qc_decisions"]),
   ...     "content_hash_prefix": manifest_dict["content_hash"][:16],
   ... }, indent=2))
   {
     "survey_id": "WILLY-L18-PROVENANCE-DEMO",
     "method": "amt",
     "operator": "pyCSAMT documentation",
     "tool_version": "2.2.1",
     "n_records": 1,
     "n_devices": 1,
     "n_files": 1,
     "n_qc_decisions": 3,
     "content_hash_prefix": "0d109730c14707c4"
   }

``n_qc_decisions`` is 3, not 2, even though only one explicit
``add_qc_decision`` call happened at this level: ``to_manifest()``
already folded the two packet-derived QC decisions into the station
record, and ``add_qc_decision`` above appended a third, manifest-level
entry alongside them — the two lists are related but not the same one.

That ``content_hash`` is what makes the manifest self-checking. Every
field of the payload — including ``tool_version`` above — is serialised
to a *canonical* JSON encoding :math:`C(\cdot)`: keys sorted, no
surrounding or separating whitespace, so the same content always
produces the same bytes regardless of field-insertion order. The
:term:`content hash` is then

.. math::

   H(M) = \mathrm{SHA256}\bigl(C(M)\bigr), \qquad
   M_{\mathrm{payload}}[\texttt{"content\_hash"}] = H(M_{\mathrm{payload}}
   \setminus \{\texttt{"content\_hash"}\}),

i.e. the hash covers every field of the manifest *except* itself, which
is what lets a reviewer strip ``content_hash`` back out, recompute
:math:`H`, and compare.

Verifying The Content Hash
--------------------------

Doing exactly that independently of pyCSAMT's own bookkeeping confirms
the manifest above was not altered after export:

.. code-block:: pycon

   >>> import hashlib
   >>> payload_no_hash = {k: v for k, v in manifest_dict.items() if k != "content_hash"}
   >>> canonical = json.dumps(
   ...     payload_no_hash, sort_keys=True, default=str, separators=(",", ":")
   ... )
   >>> recomputed = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
   >>> print(f"recomputed == content_hash: {recomputed == manifest_dict['content_hash']}")
   recomputed == content_hash: True
   >>> print(f"canonical JSON length: {len(canonical)} bytes")
   canonical JSON length: 2362 bytes

Because ``tool_version`` is part of what gets hashed, the same
acquisition inputs will still produce a *different* ``content_hash``
after a pyCSAMT upgrade — the digest is a snapshot of the manifest
bytes, not a content-only fingerprint of the survey. Pin or record
``tool_version`` alongside any digest you archive for long-term
comparison, rather than expecting it to stay constant across releases.
This dependency is sharper than it first looks: ``tool_version`` comes
from :func:`importlib.metadata.version`, which resolves package metadata
by scanning installed distributions — a stale duplicate ``.dist-info``
directory left behind by a previous install can make that lookup
ambiguous and version-dependent output like this genuinely
non-deterministic on an otherwise unchanged machine. Keep an
environment's package metadata clean for exactly this reason before
archiving a digest for long-term comparison.

Inspect A Station Record
------------------------

Each station record collects occupation metadata and QC decisions for
that station. Session-derived records are convenient when packets
already carry station, method, channel, and battery metadata.

.. code-block:: pycon

   >>> record = manifest.as_dict()["records"][0]
   >>> print(json.dumps({
   ...     "station_id": record["station_id"],
   ...     "operator": record["operator"],
   ...     "sample_rate_hz": record["sample_rate_hz"],
   ...     "battery_status": record["battery_status"],
   ...     "accepted_band_hz": record["accepted_band_hz"],
   ...     "qc_decisions": record["qc_decisions"],
   ... }, indent=2))
   {
     "station_id": "001A",
     "operator": "pyCSAMT documentation",
     "sample_rate_hz": null,
     "battery_status": null,
     "accepted_band_hz": [
       1.0,
       1000.0
     ],
     "qc_decisions": [
       {
         "station": "001A",
         "channel": null,
         "decision": "accept",
         "reasons": [],
         "operator": "pyCSAMT documentation",
         "logged_utc": "2023-11-14T22:13:20+00:00",
         "window": null
       },
       {
         "station": "001A",
         "channel": null,
         "decision": "reject",
         "reasons": [
           "finite_coverage_below_threshold"
         ],
         "operator": "pyCSAMT documentation",
         "logged_utc": "2023-11-14T22:13:20+00:00",
         "window": null
       }
     ]
   }

``sample_rate_hz`` and ``battery_status`` are ``null`` here for a more
specific reason than "the data is missing": both packets above were sent
with ``kind="qc"``, and neither ``sample_rate_hz`` nor ``battery_v`` is a
field the QC payload schema recognises — ``sample_rate_hz`` belongs to
the acquisition (``data``) payload and ``battery_v`` to the health
payload, so those keys are silently ignored when folded into the station
aggregate even though the payload dictionaries above do carry a
``battery_v`` entry. Sending at least one ``health`` packet, the way :doc:`telemetry` does
alongside its ``qc`` packets, is what actually populates
``battery_status``.

Building A Manual Audit
-----------------------

Use :class:`~pycsamt.iot.provenance.ProvenanceRecord` and
:func:`~pycsamt.iot.provenance.build_acquisition_manifest` when you are
assembling a manifest from a station sheet, field notes, and file hashes
instead of a live :term:`field session`. The resulting manifest gets the
same :math:`H(\cdot)` treatment as above:

.. code-block:: pycon

   >>> from pycsamt.iot import ProvenanceRecord, build_acquisition_manifest
   >>> qc_accept = log_qc_decision(
   ...     "001A", "accept", channel="ex",
   ...     reasons=["finite_coverage_ok", "spike_fraction_ok"],
   ...     operator="pyCSAMT documentation", timestamp=1_700_000_060.0,
   ...     window=(1_700_000_000.0, 1_700_000_060.0),
   ... )
   >>> manual_record = ProvenanceRecord(
   ...     station_id="001A", instrument_serial="DOC-REC-001", firmware="demo-2.0",
   ...     operator="pyCSAMT documentation", ex_azimuth_deg=90.0, ey_azimuth_deg=0.0,
   ...     occupation_start=1_700_000_000.0, occupation_end=1_700_000_120.0,
   ...     sample_rate_hz=512.0, gps_quality="locked", battery_status="12.3 V",
   ...     accepted_band_hz=(1.0, 1000.0),
   ...     field_notes="Manual audit record built from station sheet and QC logs.",
   ...     qc_decisions=[qc_accept, qc_reject],
   ...     processing_steps=["edge finite-coverage screening", "operator review"],
   ... )
   >>> _ = manual_record.add_raw_file(str(raw_file))
   >>> manual_manifest = build_acquisition_manifest(
   ...     "WILLY-L18-MANUAL-AUDIT", records=[manual_record], devices=[device.as_dict()],
   ...     qc_decisions=[qc_accept, qc_reject],
   ...     processing_steps=["manual audit assembled for documentation"],
   ...     method="amt", operator="pyCSAMT documentation",
   ...     metadata={"source": "documentation example"},
   ... )
   >>> manual_manifest.created_utc = "2023-11-14T22:13:20+00:00"
   >>> manual_manifest.environment = {
   ...     "python": "3.10", "platform": "documentation", "machine": "example",
   ... }
   >>> print(manual_manifest.as_dict()["content_hash"][:16])
   8df723d18df316c7

The manual and session-derived manifests hash to different prefixes even
though they describe the same station, same instrument, and overlapping
QC decisions — ``survey_id`` alone (``WILLY-L18-MANUAL-AUDIT`` versus
``WILLY-L18-PROVENANCE-DEMO``) is enough to guarantee that, since it is
one of the hashed fields.

Hash-Chained QC Logs
--------------------

A single ``content_hash`` tells a reviewer whether *the whole manifest*
still matches what was exported. A sharper question needs its own tool:
which QC entry, if any, was altered inside a running log?
:func:`~pycsamt.iot.provenance.hash_chain` answers it by folding each
entry's hash into the next, the same way a chain of *given* blocks
would, using the same canonical-JSON :math:`H(\cdot)` from above. The
first entry chains from an empty genesis hash:

.. math::

   h_0 = H\bigl(e_0 \cup \{\texttt{seq}: 0,\ \texttt{prev\_hash}: {\texttt{""}}\}\bigr),
   \qquad
   h_i = H\bigl(e_i \cup \{\texttt{seq}: i,\ \texttt{prev\_hash}: h_{i-1}\}\bigr),
   \quad i \ge 1.

Because :math:`h_i` depends on :math:`h_{i-1}`, editing, deleting, or
reordering any entry :math:`e_j` changes :math:`h_j` and therefore every
hash after it — :func:`~pycsamt.iot.provenance.verify_hash_chain` just
recomputes the recurrence and checks it still matches what was stored:

.. code-block:: pycon

   >>> from pycsamt.iot.provenance import hash_chain, verify_hash_chain
   >>> chain = hash_chain([qc_accept, qc_reject])
   >>> for entry in chain:
   ...     print(
   ...         f"seq={entry['seq']} "
   ...         f"prev_hash={entry['prev_hash'][:8] or '(genesis)'} "
   ...         f"entry_hash={entry['entry_hash'][:8]}"
   ...     )
   seq=0 prev_hash=(genesis) entry_hash=f6d7c1d0
   seq=1 prev_hash=f6d7c1d0 entry_hash=4e91a4f8
   >>> print(f"verify_hash_chain (untouched): {verify_hash_chain(chain)}")
   verify_hash_chain (untouched): True
   >>> tampered = [dict(e) for e in chain]
   >>> tampered[1]["reasons"] = ["tampered_after_the_fact"]
   >>> print(f"verify_hash_chain (tampered):  {verify_hash_chain(tampered)}")
   verify_hash_chain (tampered):  False

Editing only ``chain[1]["reasons"]`` was enough to break verification —
notice that entry 0 was never touched, but the chain still fails as a
whole because entry 1's stored ``entry_hash`` no longer matches what the
recurrence recomputes for its (now-different) contents.

Signed Manifests
----------------

A hash chain proves internal consistency, but not authorship: anyone can
recompute :math:`H(\cdot)` and produce a chain that verifies. The
remaining question is whether someone other than the exporter can prove
*they* produced this manifest.
:func:`~pycsamt.iot.provenance.sign_mapping` answers that by wrapping the
same canonical JSON in an HMAC rather than a plain hash,

.. math::

   \mathrm{signature} = \mathrm{HMAC\text{-}SHA256}\bigl(\mathrm{key}, C(M)\bigr),

so only someone holding ``key`` could have produced that signature, and
:func:`~pycsamt.iot.provenance.verify_manifest` checks both the
:term:`manifest signature` and the wrapped :term:`content hash` in one
call:

.. code-block:: pycon

   >>> from pycsamt.iot.provenance import verify_manifest
   >>> signed = manual_manifest.sign("demo-shared-key")
   >>> print(f"signature_algo: {signed['signature_algo']}")
   signature_algo: hmac-sha256
   >>> print(f"signature prefix: {signed['signature'][:16]}")
   signature prefix: 17c7bf2d87a88db9
   >>> print(f"verify_manifest (untouched): {verify_manifest(signed, 'demo-shared-key')}")
   verify_manifest (untouched): True
   >>> tampered_signed = json.loads(json.dumps(signed))
   >>> tampered_signed["manifest"]["operator"] = "someone else"
   >>> print(f"verify_manifest (tampered):  {verify_manifest(tampered_signed, 'demo-shared-key')}")
   verify_manifest (tampered):  False
   >>> print(f"verify_manifest (wrong key): {verify_manifest(signed, 'not-the-key')}")
   verify_manifest (wrong key): False

Both failure modes matter for different reasons: a tampered payload
fails because its recomputed :term:`content hash` no longer matches,
while a correct payload checked against the wrong key fails the HMAC
comparison even though the JSON itself is untouched. A :term:`hash
chain` is enough for a single-party audit trail that only needs tamper
*evidence*; add a :term:`manifest signature` when a second party needs
to verify *who* produced the manifest, not just that it is internally
consistent.

Exporting A Reproducibility Bundle
----------------------------------

Use :func:`~pycsamt.iot.provenance.export_reproducibility_bundle` to
write the manifest and one JSON audit per station. Set
``zip_bundle=True`` when you also want a zip archive for hand-off. Set
``include_raw=True`` only when you intentionally want raw files copied
into the bundle.

.. code-block:: pycon

   >>> from pycsamt.iot import export_reproducibility_bundle, export_station_audit
   >>> out_dir = "docs/source/user_guide/iot/_provenance_demo_bundle"
   >>> bundle = export_reproducibility_bundle(manual_manifest, out_dir, zip_bundle=True)
   >>> audit_path = export_station_audit(
   ...     manual_record, f"{out_dir}/single_station_audit.json",
   ... )
   >>> print(json.dumps({
   ...     "manifest_name": Path(bundle["manifest"]).name,
   ...     "audit_count": len(bundle["audits"]),
   ...     "zip_name": Path(bundle["zip"]).name,
   ...     "single_station_audit": Path(audit_path).name,
   ... }, indent=2))
   {
     "manifest_name": "acquisition_manifest.json",
     "audit_count": 1,
     "zip_name": "_provenance_demo_bundle.zip",
     "single_station_audit": "single_station_audit.json"
   }

The manifest gives reviewers enough information to understand and
verify the acquisition trail: which station was occupied, which device
was used, which QC decisions were made, which raw file was referenced,
and whether the manifest content changed after export. Raw-file hashes
protect the integrity of data files; a :term:`content hash` protects the
manifest that describes them; a :term:`hash chain` protects the order
and content of a running QC log; and a :term:`manifest signature`
protects who is allowed to be believed as the source of all of the
above. Station audits preserve the field context around those files
independently of any of the four.

For production surveys, write the manifest beside the processed
products and keep the raw-file hash records even when raw files are too
large to bundle, and record the exporting ``tool_version`` alongside any
archived digest. That makes later reporting, reprocessing, and external
review much less fragile.
