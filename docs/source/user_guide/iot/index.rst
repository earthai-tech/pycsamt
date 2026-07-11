.. _user_guide_iot:

IoT-Enabled Field Acquisition
=============================

The :mod:`pycsamt.iot` subpackage documents and audits IoT-enabled AMT,
MT, CSAMT, and CSEM acquisition. It does not replace the normal impedance,
quality-control, dimensionality, or inversion workflow. Instead, it records
the operational layer around those data: devices, stations, packets, edge
QC decisions, power, clock synchronisation, security configuration, and
provenance. It also bridges back into the science API, turning edge
impedance estimates into a :class:`pycsamt.z.z.Z` and a preliminary EDI,
and seeding a re-occupation session from an existing EDI survey.

Use this guide to configure and audit the operational layer of an
IoT-enabled field survey: telemetry, edge QC, monitoring, power,
synchronisation, security, and provenance.

Guide Sections
--------------

.. grid:: 1 1 2 3
   :gutter: 3
   :class-container: cta-tiles

   .. grid-item-card:: Basic Session
      :link: basic_session
      :link-type: doc
      :img-top: ../../_static/icons/session.svg
      :class-card: pycsamt-card sd-text-center

      Build a ``FieldSession`` from devices, stations, telemetry packets,
      station tables, and pipeline hand-off metadata.

   .. grid-item-card:: Edge QC
      :link: edge_qc
      :link-type: doc
      :img-top: ../../_static/icons/edge.svg
      :class-card: pycsamt-card sd-text-center

      Reduce live windows, score finite coverage and spikes, and encode edge
      QC decisions as telemetry.

   .. grid-item-card:: AMT Diagnostics
      :link: amt_diagnostics
      :link-type: doc
      :img-top: ../../_static/icons/diagnostic.svg
      :class-card: pycsamt-card sd-text-center

      Check SNR, powerline harmonics, frequency coverage, sensor dropout,
      contact resistance, and impedance stability.

   .. grid-item-card:: Telemetry
      :link: telemetry
      :link-type: doc
      :img-top: ../../_static/icons/telemetry.svg
      :class-card: pycsamt-card sd-text-center

      Create canonical packets, use dry-run transports, write JSONL logs,
      replay messages, and inspect packet tables.

   .. grid-item-card:: Monitoring
      :link: monitoring
      :link-type: doc
      :img-top: ../../_static/icons/monitoring.svg
      :class-card: pycsamt-card sd-text-center

      Turn telemetry packets into survey status, packet success rates,
      latency checks, and issue lists.

   .. grid-item-card:: Visualization
      :link: visualization
      :link-type: doc
      :img-top: ../../_static/icons/map-and-profile.svg
      :class-card: pycsamt-card sd-text-center

      Plot dashboards, edge-QC summaries, power budgets, and clock quality
      panels from IoT sessions and packets.

   .. grid-item-card:: Power Management
      :link: power_management
      :link-type: doc
      :img-top: ../../_static/icons/power.svg
      :class-card: pycsamt-card sd-text-center

      Estimate load, harvest, runtime, autonomy, power state, and power
      telemetry for field devices.

   .. grid-item-card:: Clock Sync
      :link: clock_sync
      :link-type: doc
      :img-top: ../../_static/icons/clock.svg
      :class-card: pycsamt-card sd-text-center

      Assess clock offset, drift, jitter, reference support, GPS lock, and
      synchronisation quality.

   .. grid-item-card:: Provenance
      :link: provenance
      :link-type: doc
      :img-top: ../../_static/icons/provenance.svg
      :class-card: pycsamt-card sd-text-center

      Record acquisition manifests, file hashes, QC decisions, station audits,
      and reproducibility bundles.

   .. grid-item-card:: Security
      :link: security
      :link-type: doc
      :img-top: ../../_static/icons/security.svg
      :class-card: pycsamt-card sd-text-center

      Configure credentials, TLS settings, authentication modes, secret
      redaction, and deployment security summaries.

   .. grid-item-card:: Simulation
      :link: simulation
      :link-type: doc
      :img-top: ../../_static/icons/simulation.svg
      :class-card: pycsamt-card sd-text-center

      Generate deterministic AMT-like channels, stations, networks, packet
      loss, GPS drift, and battery decay for demos and tests.

.. toctree::
   :maxdepth: 2
   :caption: IoT Topics
   :hidden:

   basic_session
   edge_qc
   amt_diagnostics
   controlled_source
   method_profiles
   data_bridge
   telemetry
   monitoring
   visualization
   power_management
   clock_sync
   provenance
   security
   simulation

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
   include ``data``, ``qc``, ``health``, ``sync``, ``power``, ``event``,
   and ``source`` (controlled-source transmitter telemetry).

``FieldSession``
   A stateful survey container that links devices, stations, packets,
   monitoring, station tables, pipeline inputs, and provenance manifests.

``MethodProfile``
   Canonical acquisition characteristics for an :class:`~pycsamt.iot.monitoring.EMMethod`
   (AMT, MT, CSAMT, CSEM, TDEM/TEM): frequency band, required channels,
   nominal sample rate, and controlled-source / powerline-sensitivity
   flags. Profiles drive method-aware QC.

Scope Note
----------

The IoT layer records operational evidence around acquisition. It does not
change the electromagnetic inversion itself.

Edge packets can record finite-data coverage, spike and harmonic
contamination indicators, frequency coverage, channel health, and
accept/reject decisions. These metrics help identify poor acquisition
windows and operational faults before downstream AMT/CSAMT processing.
