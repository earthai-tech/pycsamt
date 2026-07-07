.. _user_guide_iot:

IoT-Enabled Field Acquisition
=============================

The :mod:`pycsamt.iot` subpackage documents and audits IoT-enabled AMT,
MT, and CSAMT acquisition. It does not replace the normal impedance,
quality-control, dimensionality, or inversion workflow. Instead, it records
the operational layer around those data: devices, stations, packets, edge
QC decisions, power, clock synchronisation, security configuration, and
provenance.

Use this page to configure and audit the operational layer of an
IoT-enabled field survey: telemetry, edge QC, monitoring, power,
synchronisation, security, and provenance.

Core Concepts
-------------

``DeviceConfig``
   A recorder, gateway, remote-reference node, or sensor node. It stores
   protocol, sample rate, station assignment, channels, and role.

``StationConfig``
   A field occupation location. It stores coordinates, profile chainage,
   channel list, dipole geometry, orientation, operator, and attached
   devices.

``TelemetryPacket``
   A timestamped message with a canonical topic and payload. Packet kinds
   include ``data``, ``qc``, ``health``, ``sync``, ``power``, and ``event``.

``FieldSession``
   A stateful survey container that links devices, stations, packets,
   monitoring, station tables, pipeline inputs, and provenance manifests.

Basic Session
-------------

.. code-block:: python

   from pycsamt.iot import (
       DeviceConfig,
       FieldSession,
       MonitoringConfig,
       StationConfig,
   )

   device = DeviceConfig(
       "node-1",
       station="S01",
       protocol="mqtt",
       sample_rate_hz=256.0,
       channels=["ex", "ey", "hx", "hy"],
       role="recorder",
   )

   station = StationConfig(
       "S01",
       lat=6.82,
       lon=-5.28,
       profile="L1",
       position_m=120.0,
       dipole_length_m=50.0,
       ex_azimuth_deg=90.0,
       ey_azimuth_deg=0.0,
   )

   session = FieldSession(
       "survey-a",
       devices=[device],
       stations=[station],
       monitoring_config=MonitoringConfig(
           method="amt",
           required_channels=["ex", "ey", "hx", "hy"],
           min_edge_acceptance_rate=0.85,
       ),
   )

Generic Edge QC
---------------

Use :class:`pycsamt.iot.EdgeProcessor` for lightweight node-side checks on
time-series windows before telemetry transmission or storage.

.. code-block:: python

   import numpy as np
   from pycsamt.iot import EdgeProcessingConfig, EdgeProcessor

   data = np.random.randn(2048, 4)
   processor = EdgeProcessor(
       EdgeProcessingConfig(
           decimation=2,
           finite_threshold=0.95,
           channel_names=["ex", "ey", "hx", "hy"],
           spike_threshold=6.0,
           max_spike_fraction=0.05,
       )
   )

   result = processor.process(data)
   qc_packet = result.to_packet(device, timestamp=1_700_000_000.0,
                                survey_id="survey-a")
   session.add_packet(qc_packet)

The result contains global metrics, per-channel summaries, rejection
reasons, and an accept/reject decision. This is the correct layer for
recording field-side decimation, finite-data coverage, spike screening, and
QC telemetry metadata.

AMT/CSAMT Edge Diagnostics
--------------------------

The :mod:`pycsamt.iot.edge_amt` module adds diagnostics specific to field
electromagnetics:

* ``detect_powerline_harmonics`` checks 50/60 Hz harmonics and reports the
  contaminated spectral power ratio.
* ``estimate_channel_snr`` estimates channel SNR from live windows.
* ``check_channel_saturation`` flags clipping or ADC saturation.
* ``check_contact_resistance`` provides a proxy check for electric dipole
  contact quality.
* ``estimate_frequency_coverage`` reports which frequency band is
  resolvable from the window length and sample rate.
* ``assess_impedance_stability`` compares apparent impedance stability
  across windows.
* ``detect_sensor_dropout`` identifies missing or flat sensor intervals.

.. code-block:: python

   from pycsamt.iot import (
       detect_powerline_harmonics,
       estimate_frequency_coverage,
   )

   harmonics = detect_powerline_harmonics(
       data[:, 0],
       sample_rate=256.0,
       mains_hz=50.0,
       n_harmonics=5,
       threshold_ratio=0.05,
   )

   coverage = estimate_frequency_coverage(
       data[:, 0],
       sample_rate=256.0,
       target_bands=[(1.0, 1000.0)],
   )

Telemetry Transports
--------------------

Use :func:`pycsamt.iot.build_telemetry_client` to construct a client. All
clients support ``dry_run=True`` so field workflows can be tested without
opening a network connection.

.. code-block:: python

   from pycsamt.iot import build_telemetry_client

   client = build_telemetry_client(
       "file",
       endpoint="survey-a-telemetry.jsonl",
       dry_run=False,
   )
   ack = client.send(qc_packet)

Supported concrete transports are file, HTTP(S), MQTT, serial, and
WebSocket. MQTT, serial, and WebSocket require their optional third-party
client libraries when ``dry_run=False``.

Monitoring
----------

Telemetry monitoring turns packets into an operational status. It checks
packet success, edge acceptance, latency, packet gaps, battery voltage,
clock offset, required channels, and frequency band.

.. code-block:: python

   from pycsamt.iot import TelemetryMonitor, monitoring_status_table

   status = TelemetryMonitor(session.monitoring_config).assess(
       session.packets,
       now=1_700_000_030.0,
   )
   table = monitoring_status_table(status)

The status level is one of ``ok``, ``warning``, ``critical``, or
``no_data``. The ``issues`` field gives machine-readable operational
reasons such as ``edge_acceptance_rate_below_threshold``,
``battery_below_threshold``, ``packet_gap_above_threshold``, or
``required_channels_missing``.

Visualization
-------------

Use :func:`pycsamt.iot.plot_field_dashboard` to summarise an acquisition
session in one figure. The dashboard combines station health/profile,
edge-QC acceptance, power or runtime telemetry, synchronisation labels, and
packet timing. Use :func:`pycsamt.iot.plot_edge_qc_summary` when the focus
is specifically on edge-QC channel metrics and rejection reasons.
Use :func:`pycsamt.iot.plot_power_budget` to inspect energy-management
metrics across devices or power telemetry packets, and
:func:`pycsamt.iot.plot_sync_quality` to inspect clock offset, drift,
jitter, GPS lock, and synchronisation grades.

.. code-block:: python

   from pycsamt.iot import (
       plot_edge_qc_summary,
       plot_field_dashboard,
       plot_power_budget,
       plot_sync_quality,
   )

   fig = plot_field_dashboard(session, station_axis="auto")
   fig.savefig("survey-a-iot-dashboard.png", dpi=150)

   qc_fig = plot_edge_qc_summary(session)
   qc_fig.savefig("survey-a-edge-qc.png", dpi=150)

   power_fig = plot_power_budget(session)
   power_fig.savefig("survey-a-power-budget.png", dpi=150)

   sync_fig = plot_sync_quality(session)
   sync_fig.savefig("survey-a-sync-quality.png", dpi=150)

The returned figure has a ``pycsamt_iot_dashboard`` attribute containing
the station, packet, and monitoring data used to draw the panels.
The QC summary figure similarly stores normalised rows in
``pycsamt_iot_edge_qc``.
The power-budget figure stores normalised rows in
``pycsamt_iot_power_budget``.
The synchronisation figure stores normalised rows in
``pycsamt_iot_sync_quality``.

Power Management
----------------

Energy management is represented by :class:`pycsamt.iot.EnergyConfig` and
:func:`pycsamt.iot.estimate_energy_budget`. The estimate accounts for
battery reserve, active/sleep duty cycle, regulator efficiency, solar
harvesting, charge efficiency, telemetry radio windows, edge-processing
load, auxiliary load, and minimum runtime targets.

.. code-block:: python

   from pycsamt.iot import EnergyConfig, estimate_energy_budget

   estimate = estimate_energy_budget(
       EnergyConfig(
           battery_wh=120.0,
           active_power_w=2.0,
           sleep_power_w=0.2,
           duty_cycle=0.25,
           reserve_fraction=0.2,
           regulator_efficiency=0.8,
           solar_wh_per_day=10.0,
           charge_efficiency=0.9,
           telemetry_power_w=3.0,
           telemetry_seconds_per_day=600.0,
           edge_power_w=0.5,
           edge_duty_cycle=0.5,
           min_runtime_days=8.0,
       )
   )

   power_packet = estimate.to_packet(device, timestamp=1_700_000_000.0,
                                     survey_id="survey-a")

Use ``power_summary_table`` or ``estimate_deployment_energy`` for
deployment-scale summaries.

Clock Synchronisation
---------------------

The synchronisation module estimates clock offset, drift, and jitter
against a reference clock such as GPS. It grades each device with
``excellent``, ``good``, ``fair``, ``poor``, or ``unknown``.

.. code-block:: python

   from pycsamt.iot import ClockSynchronizer, SyncConfig

   sync = ClockSynchronizer(
       SyncConfig(tolerance_ms=1.0, max_drift_ppm=10.0, max_jitter_ms=1.0)
   )
   status = sync.assess(
       "node-1",
       local_timestamps=[1.0005, 2.0005, 3.0006],
       reference_timestamps=[1.0, 2.0, 3.0],
       gps_lock=True,
   )

Provenance and Reproducibility
------------------------------

Use the provenance helpers to create a reproducible audit trail for a
field session:

* station occupation records;
* operator, firmware, instrument serial, orientation, and GPS metadata;
* accepted/rejected QC decisions and reasons;
* processing steps;
* raw-file hashes; and
* JSON manifests or zipped reproducibility bundles.

.. code-block:: python

   manifest = session.to_manifest()

Security
--------

Security configuration is serialisable and redacts secrets by default.
Transport encryption and authentication are delegated to the concrete
transport layer.

.. code-block:: python

   from pycsamt.iot import AuthScheme, Credential, SecurityConfig, TLSConfig

   security = SecurityConfig(
       tls=TLSConfig(enabled=True, ca_cert="ca.pem"),
       credential=Credential(
           scheme=AuthScheme.BEARER,
           token="secret-token",
       ),
   )
   options = security.client_options()

Simulation
----------

The simulator generates repeatable synthetic IoT/AMT field data for tests
and demonstrations:

* AMT channel windows with configurable SNR, powerline noise, and dropouts;
* station and device inventories;
* packet loss;
* GPS drift;
* battery decay.

.. code-block:: python

   from pycsamt.iot import simulate_amt_channel, simulate_iot_network

   ex = simulate_amt_channel(
       4096,
       sample_rate=256.0,
       snr_db=20.0,
       powerline_amplitude=0.1,
       seed=42,
   )

Scope Note
----------

The IoT layer records operational evidence around acquisition. It does not
change the electromagnetic inversion itself.

Edge packets can record finite-data coverage, spike and harmonic
contamination indicators, frequency coverage, channel health, and
accept/reject decisions. These metrics help identify poor acquisition
windows and operational faults before downstream AMT/CSAMT processing.
