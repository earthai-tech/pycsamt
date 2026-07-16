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

.. toctree::
   :maxdepth: 3
   :class: pycsamt-guide-toc

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
