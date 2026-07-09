pycsamt.iot
===========

IoT field telemetry, edge processing, synchronisation, power budgeting,
security, provenance, and simulation helpers for AMT, MT, and CSAMT
deployments.

The :mod:`pycsamt.iot` package is the public API for IoT-enabled field
acquisition. It is intentionally separate from the EDI/impedance science
API: IoT objects describe what happened in the field, what packets were
sent, what edge decisions were made, how clocks and power behaved, and how
the acquisition can be audited. Processed transfer functions still flow
through the normal :mod:`pycsamt.emtools`, :mod:`pycsamt.site`, and
:mod:`pycsamt.z` layers.

Implemented Capabilities
------------------------

The package implements concrete capabilities that can be reported in a
methods section:

* device and station metadata, including station coordinates, line
  position, electric dipole geometry, channels, and recorder IDs;
* canonical telemetry packets and payload schemas for data, QC, health,
  power, synchronisation, event, and acquisition messages;
* dry-run and real telemetry clients for file, HTTP(S), MQTT, serial, and
  WebSocket transports, with a generic recorder fallback for protocols
  such as LoRa;
* generic edge QC: decimation, finite coverage, RMS/statistics, robust
  spike fraction, per-channel accept/reject decisions, and QC telemetry
  packet export;
* AMT/CSAMT edge QC: powerline harmonic detection, channel SNR, saturation,
  contact-resistance proxy checks, frequency coverage, live spectra,
  impedance stability, and sensor dropout checks;
* telemetry stream monitoring for AMT, MT, CSAMT, TEM, and TDEM, including
  packet success, edge acceptance, latency, packet gaps, clock offset,
  battery voltage, required channels, and frequency-band checks;
* power-budget estimation with duty cycle, battery reserve, regulator and
  charge efficiency, solar harvesting, telemetry load, edge-processing
  load, auxiliary load, status classification, and power packet export;
* clock synchronisation audit: offset, drift, jitter, GPS lock, quality
  grade, batch status tables, and GPS-dropout detection;
* acquisition provenance: station audit records, QC decision logs, raw-file
  hashes, manifests, and reproducibility bundles;
* transport security configuration with TLS options, credential schemes,
  and secret redaction; and
* reproducible synthetic field-network simulation for documentation,
  tests, and demonstrations without hardware.
* acquisition visualisation with a compact field dashboard for station
  health, edge-QC acceptance, power/synchronisation status, and telemetry
  packet timing, plus a dedicated edge-QC summary for channel coverage,
  spike fraction, decisions, rejection reasons, and a power-budget summary
  for daily load, harvest, runtime, and energy-management state, and a
  synchronisation-quality summary for offset, drift, jitter, GPS lock, and
  quality grades.

Minimal Example
---------------

.. code-block:: python

   import numpy as np
   from pycsamt.iot import (
       DeviceConfig,
       EdgeProcessingConfig,
       EdgeProcessor,
       FieldSession,
       MonitoringConfig,
   )

   device = DeviceConfig(
       "node-1",
       station="S01",
       protocol="mqtt",
       sample_rate_hz=256.0,
       channels=["ex", "hy"],
   )

   edge = EdgeProcessor(
       EdgeProcessingConfig(
           finite_threshold=0.95,
           channel_names=["ex", "hy"],
           spike_threshold=6.0,
       )
   )
   result = edge.process(np.random.randn(1024, 2))
   packet = result.to_packet(device, timestamp=1_700_000_000.0,
                             survey_id="survey-a")

   session = FieldSession(
       "survey-a",
       devices=[device],
       monitoring_config=MonitoringConfig(
           method="amt",
           required_channels=["ex", "hy"],
           min_edge_acceptance_rate=0.85,
       ),
   )
   session.add_packet(packet)
   status = session.assess()

Packet and Table API
--------------------

Most functions return plain Python objects or pandas-like tables wrapped by
the pyCSAMT API view layer when requested. Common reporting helpers include:

* :func:`pycsamt.iot.deployment_report`
* :func:`pycsamt.iot.station_table`
* :func:`pycsamt.iot.edge_summary_table`
* :func:`pycsamt.iot.amt_edge_table`
* :func:`pycsamt.iot.packet_table`
* :func:`pycsamt.iot.telemetry_summary`
* :func:`pycsamt.iot.monitoring_status_table`
* :func:`pycsamt.iot.power_summary_table`
* :func:`pycsamt.iot.estimate_deployment_energy`
* :func:`pycsamt.iot.sync_status_table`
* :func:`pycsamt.iot.plot_field_dashboard`
* :func:`pycsamt.iot.plot_edge_qc_summary`
* :func:`pycsamt.iot.plot_power_budget`
* :func:`pycsamt.iot.plot_sync_quality`

.. automodule:: pycsamt.iot
   :members:
   :show-inheritance:

Modules
-------

.. autosummary::
   :toctree: generated

   pycsamt.iot.core
   pycsamt.iot.station
   pycsamt.iot.schemas
   pycsamt.iot.session
   pycsamt.iot.protocols
   pycsamt.iot.protocols.base
   pycsamt.iot.protocols.file
   pycsamt.iot.protocols.http
   pycsamt.iot.protocols.mqtt
   pycsamt.iot.protocols.serial
   pycsamt.iot.protocols.websocket
   pycsamt.iot.edge
   pycsamt.iot.edge_amt
   pycsamt.iot.sync
   pycsamt.iot.power
   pycsamt.iot.monitoring
   pycsamt.iot.plot
   pycsamt.iot.provenance
   pycsamt.iot.security
   pycsamt.iot.sim
