.. _user_guide_iot_provenance:

Provenance and Reproducibility
==============================

Provenance records explain what happened during acquisition and make a
field session reproducible after the instruments have left the site. In an
IoT workflow, that audit trail should include station occupation metadata,
devices, QC decisions, raw-file integrity hashes, processing steps, and
the environment used to create the manifest.

This page uses one real file from the repository demo data,
``data/AMT/WILLY_DATA/L18PLT/18-001A.edi``, for the raw-file hash. The IoT
QC packets are synthetic because provenance is documenting acquisition
events around the data, not reprocessing the EDI file itself.

Hash A Raw File
---------------

Use :func:`pycsamt.iot.hash_raw_file` to create an integrity record. The
full digest is stored in manifests; the example prints only a prefix to
keep the documentation readable.

.. code-block:: python
   :linenos:

   from pathlib import Path

   from pycsamt.iot import hash_raw_file

   raw_file = Path("data/AMT/WILLY_DATA/L18PLT/18-001A.edi")
   raw_hash = hash_raw_file(str(raw_file))
   print(
       {
           "name": raw_hash["name"],
           "bytes": raw_hash["bytes"],
           "algo": raw_hash["algo"],
           "digest_prefix": raw_hash["digest"][:16],
       }
   )

Output:

.. code-block:: text

   {
     "name": "18-001A.edi",
     "bytes": 22758,
     "algo": "sha256",
     "digest_prefix": "a7bc0c7b290c1847"
   }

Log QC Decisions
----------------

Use :func:`pycsamt.iot.log_qc_decision` for normalized audit records. It
standardizes station IDs, channel names, decisions, reasons, timestamps,
and optional acquisition windows.

.. code-block:: python
   :linenos:

   from pycsamt.iot import log_qc_decision

   qc_reject = log_qc_decision(
       "001A",
       "reject",
       channel="ey",
       reasons=["finite_coverage_below_threshold"],
       operator="pyCSAMT documentation",
       timestamp=1_700_000_120.0,
       window=(1_700_000_060.0, 1_700_000_120.0),
   )
   print(qc_reject)

Output:

.. code-block:: text

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

Build A Session Manifest
------------------------

The easiest route is to build a :class:`pycsamt.iot.FieldSession`, add
devices, stations, and packets, then call
:meth:`pycsamt.iot.FieldSession.to_manifest`. The manifest can still be
enriched afterward with raw-file hashes, extra QC decisions, and processing
steps.

.. code-block:: python
   :linenos:

   from pycsamt.iot import (
       DeviceConfig,
       FieldSession,
       MonitoringConfig,
       StationConfig,
       TelemetryPacket,
   )

   device = DeviceConfig(
       "l18-node-01",
       station="001A",
       protocol="file",
       sample_rate_hz=512.0,
       channels=["ex", "ey", "hx", "hy"],
       role="recorder",
       metadata={
           "instrument_serial": "DOC-REC-001",
           "firmware": "demo-2.0",
           "source_edi": raw_file.as_posix(),
       },
   )
   station = StationConfig(
       "001A",
       profile="L18",
       position_m=0.0,
       channels=["ex", "ey", "hx", "hy"],
       dipole_length_m=50.0,
       ex_azimuth_deg=90.0,
       ey_azimuth_deg=0.0,
       device_ids=[device.device_id],
       operator="pyCSAMT documentation",
       notes="Station id and raw EDI file are from the repository demo data.",
   )
   session = FieldSession(
       "WILLY-L18-PROVENANCE-DEMO",
       devices=[device],
       stations=[station],
       method="amt",
       operator="pyCSAMT documentation",
       monitoring_config=MonitoringConfig(
           method="amt",
           required_channels=["ex", "ey", "hx", "hy"],
           frequency_band_hz=(1.0, 1000.0),
       ),
       metadata={"real_data_path": "data/AMT/WILLY_DATA/L18PLT"},
   )

   for index, (accepted, decision, reasons) in enumerate(
       [
           (True, "accept", []),
           (False, "reject", ["finite_coverage_below_threshold"]),
       ]
   ):
       session.add_packet(
           TelemetryPacket.from_device(
               device,
               timestamp=1_700_000_000.0 + 60.0 * index,
               survey_id=session.survey_id,
               kind="qc",
               payload={
                   "method": "amt",
                   "station": "001A",
                   "channels": ["ex", "ey", "hx", "hy"],
                   "frequency_band_hz": [1.0, 1000.0],
                   "accepted": accepted,
                   "decision": decision,
                   "reasons": reasons,
                   "battery_v": 12.4 - 0.1 * index,
               },
           )
       )

   manifest = session.to_manifest()
   manifest.created_utc = "2023-11-14T22:13:20+00:00"
   manifest.environment = {
       "python": "3.10",
       "platform": "documentation",
       "machine": "example",
   }
   for record in manifest.records:
       for item in record.qc_decisions:
           item["logged_utc"] = "2023-11-14T22:13:20+00:00"
   manifest.add_file(str(raw_file))
   manifest.add_qc_decision(
       station="001A",
       decision="reject",
       channel="ey",
       reasons=["finite_coverage_below_threshold"],
       operator="pyCSAMT documentation",
       timestamp=1_700_000_120.0,
       window=(1_700_000_060.0, 1_700_000_120.0),
   )
   manifest.add_processing_step("loaded station metadata from L18 demo line")
   manifest.add_processing_step("recorded synthetic IoT QC packet decisions")

   manifest_dict = manifest.as_dict()
   print(
       {
           "survey_id": manifest_dict["survey_id"],
           "method": manifest_dict["method"],
           "operator": manifest_dict["operator"],
           "n_records": len(manifest_dict["records"]),
           "n_devices": len(manifest_dict["devices"]),
           "n_files": len(manifest_dict["files"]),
           "n_qc_decisions": len(manifest_dict["qc_decisions"]),
           "content_hash_prefix": manifest_dict["content_hash"][:16],
       }
   )

Output:

.. code-block:: text

   {
     "survey_id": "WILLY-L18-PROVENANCE-DEMO",
     "method": "amt",
     "operator": "pyCSAMT documentation",
     "n_records": 1,
     "n_devices": 1,
     "n_files": 1,
     "n_qc_decisions": 3,
     "content_hash_prefix": "f609c0a3c301389a"
   }

Inspect A Station Record
------------------------

Each station record collects occupation metadata and QC decisions for that
station. Session-derived records are convenient when packets already carry
station, method, channel, and battery metadata.

.. code-block:: python
   :linenos:

   record = manifest.as_dict()["records"][0]
   print(
       {
           "station_id": record["station_id"],
           "operator": record["operator"],
           "sample_rate_hz": record["sample_rate_hz"],
           "battery_status": record["battery_status"],
           "accepted_band_hz": record["accepted_band_hz"],
           "qc_decisions": record["qc_decisions"],
       }
   )

Output:

.. code-block:: text

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

Build A Manual Audit
--------------------

Use :class:`pycsamt.iot.ProvenanceRecord` and
:func:`pycsamt.iot.build_acquisition_manifest` when you are assembling a
manifest from a station sheet, field notes, and file hashes instead of a
live :class:`~pycsamt.iot.FieldSession`.

.. code-block:: python
   :linenos:

   from pycsamt.iot import ProvenanceRecord, build_acquisition_manifest

   qc_accept = log_qc_decision(
       "001A",
       "accept",
       channel="ex",
       reasons=["finite_coverage_ok", "spike_fraction_ok"],
       operator="pyCSAMT documentation",
       timestamp=1_700_000_060.0,
       window=(1_700_000_000.0, 1_700_000_060.0),
   )
   manual_record = ProvenanceRecord(
       station_id="001A",
       instrument_serial="DOC-REC-001",
       firmware="demo-2.0",
       operator="pyCSAMT documentation",
       ex_azimuth_deg=90.0,
       ey_azimuth_deg=0.0,
       occupation_start=1_700_000_000.0,
       occupation_end=1_700_000_120.0,
       sample_rate_hz=512.0,
       gps_quality="locked",
       battery_status="12.3 V",
       accepted_band_hz=(1.0, 1000.0),
       field_notes="Manual audit record built from station sheet and QC logs.",
       qc_decisions=[qc_accept, qc_reject],
       processing_steps=["edge finite-coverage screening", "operator review"],
   )
   manual_record.add_raw_file(str(raw_file))

   manual_manifest = build_acquisition_manifest(
       "WILLY-L18-MANUAL-AUDIT",
       records=[manual_record],
       devices=[device.as_dict()],
       qc_decisions=[qc_accept, qc_reject],
       processing_steps=["manual audit assembled for documentation"],
       method="amt",
       operator="pyCSAMT documentation",
       metadata={"source": "documentation example"},
   )
   manual_manifest.created_utc = "2023-11-14T22:13:20+00:00"
   manual_manifest.environment = {
       "python": "3.10",
       "platform": "documentation",
       "machine": "example",
   }
   print(manual_manifest.as_dict()["content_hash"][:16])

Output:

.. code-block:: text

   12d8f383d727517a

Export A Reproducibility Bundle
-------------------------------

Use :func:`pycsamt.iot.export_reproducibility_bundle` to write the
manifest and one JSON audit per station. Set ``zip_bundle=True`` when you
also want a zip archive for hand-off. Set ``include_raw=True`` only when
you intentionally want raw files copied into the bundle.

.. code-block:: python
   :linenos:

   from pycsamt.iot import (
       export_reproducibility_bundle,
       export_station_audit,
   )

   out_dir = "docs/source/user_guide/iot/_provenance_demo_bundle"
   bundle = export_reproducibility_bundle(
       manual_manifest,
       out_dir,
       zip_bundle=True,
   )
   audit_path = export_station_audit(
       manual_record,
       f"{out_dir}/single_station_audit.json",
   )
   print(
       {
           "manifest_name": Path(bundle["manifest"]).name,
           "audit_count": len(bundle["audits"]),
           "zip_name": Path(bundle["zip"]).name,
           "single_station_audit": Path(audit_path).name,
       }
   )

Output:

.. code-block:: text

   {
     "manifest_name": "acquisition_manifest.json",
     "audit_count": 1,
     "zip_name": "_provenance_demo_bundle.zip",
     "single_station_audit": "single_station_audit.json"
   }

Field Interpretation
--------------------

The manifest gives reviewers enough information to understand and verify
the acquisition trail: which station was occupied, which device was used,
which QC decisions were made, which raw file was referenced, and whether
the manifest content changed after export. The ``content_hash`` is a
stable digest of the manifest payload. Raw-file hashes protect the
integrity of data files, while station audits preserve the field context
around those files.

For production surveys, write the manifest beside the processed products
and keep the raw-file hash records even when raw files are too large to
bundle. That makes later reporting, reprocessing, and external review much
less fragile.
