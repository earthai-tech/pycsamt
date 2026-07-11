.. _user_guide_iot_telemetry:

Telemetry Transports
====================

Telemetry is the message layer between IoT field hardware and the rest of
the pyCSAMT acquisition workflow. A telemetry packet carries a device id,
timestamp, canonical topic, packet kind, QoS flag, retained flag, and a
JSON-serialisable payload. The same packet can be sent to MQTT, HTTP,
serial, WebSocket, a JSON-lines file, or a dry-run recorder.

Use dry-run clients while developing notebooks, dashboards, and examples:
they validate packets and record what would have been sent without opening
network or hardware connections. Use the file transport when you want an
offline replay log that can be read back later.

Create Canonical Packets
------------------------

Start by declaring the field devices. A
:class:`pycsamt.iot.DeviceConfig` knows the station, acquisition channels,
preferred protocol, and canonical topic layout. Use
:meth:`pycsamt.iot.TelemetryPacket.from_device` so the topic is built
consistently.

.. code-block:: python
   :linenos:

   from pycsamt.iot import DeviceConfig, PacketKind, TelemetryPacket

   devices = [
       DeviceConfig(
           "l18-node-01",
           station="L18-001",
           protocol="mqtt",
           sample_rate_hz=256.0,
           channels=["ex", "ey", "hx", "hy"],
       ),
       DeviceConfig(
           "l18-node-02",
           station="L18-002",
           protocol="mqtt",
           sample_rate_hz=256.0,
           channels=["ex", "ey", "hx", "hy"],
       ),
   ]

   packets = [
       TelemetryPacket.from_device(
           devices[0],
           timestamp=1_700_000_000.0,
           survey_id="WILLY-L18-TELEMETRY",
           kind=PacketKind.HEALTH,
           payload={
               "station": "L18-001",
               "battery_v": 12.62,
               "temperature_c": 31.4,
               "firmware": "iot-demo-0.4.0",
           },
           qos=1,
       ),
       TelemetryPacket.from_device(
           devices[0],
           timestamp=1_700_000_006.0,
           survey_id="WILLY-L18-TELEMETRY",
           kind=PacketKind.QC,
           payload={
               "station": "L18-001",
               "method": "amt",
               "channels": ["ex", "ey", "hx", "hy"],
               "accepted": True,
               "decision": "accept",
               "finite_coverage": 0.994,
               "frequency_band_hz": [2.0, 512.0],
           },
           qos=1,
       ),
       TelemetryPacket.from_device(
           devices[1],
           timestamp=1_700_000_012.0,
           survey_id="WILLY-L18-TELEMETRY",
           kind=PacketKind.HEALTH,
           payload={
               "station": "L18-002",
               "battery_v": 11.38,
               "temperature_c": 33.1,
               "firmware": "iot-demo-0.4.0",
           },
           qos=1,
       ),
       TelemetryPacket.from_device(
           devices[1],
           timestamp=1_700_000_018.0,
           survey_id="WILLY-L18-TELEMETRY",
           kind=PacketKind.QC,
           payload={
               "station": "L18-002",
               "method": "amt",
               "channels": ["ex", "ey", "hx", "hy"],
               "accepted": False,
               "decision": "reject",
               "finite_coverage": 0.861,
               "frequency_band_hz": [2.0, 512.0],
               "reasons": ["finite_coverage_below_threshold"],
           },
           qos=1,
       ),
       TelemetryPacket.from_device(
           devices[1],
           timestamp=1_700_000_024.0,
           survey_id="WILLY-L18-TELEMETRY",
           kind=PacketKind.SYNC,
           payload={
               "station": "L18-002",
               "quality": "fair",
               "clock_offset_ms": 4.8,
               "drift_ppm": 2.1,
           },
           qos=0,
       ),
   ]

   first = packets[0]
   print(f"device_id: {first.device_id}")
   print(f"kind: {first.kind.value}")
   print(f"topic: {first.topic}")
   print(f"qos: {first.qos}")
   print(f"retained: {first.retained}")
   print(f"payload_keys: {', '.join(sorted(first.payload))}")

Output:

.. code-block:: text

   device_id: l18-node-01
   kind: health
   topic: pycsamt/WILLY-L18-TELEMETRY/L18-001/l18-node-01/health
   qos: 1
   retained: False
   payload_keys: battery_v, firmware, station, temperature_c

Record A Dry-Run Transport
--------------------------

Use :func:`pycsamt.iot.build_telemetry_client` for every transport. The
default ``dry_run=True`` is deliberate: examples and tests can use MQTT,
HTTP, serial, WebSocket, or LoRa-style names without needing a broker,
port, modem, or optional client library.

.. code-block:: python
   :linenos:

   from pycsamt.iot import build_telemetry_client

   client = build_telemetry_client(
       "mqtt",
       endpoint="mqtt://field-gateway.local:1883",
       dry_run=True,
   )
   client.subscribe("pycsamt/WILLY-L18-TELEMETRY/+/+/qc")

   for packet in packets[:2]:
       ack = client.send(packet)

   print(f"protocol: {ack.protocol}")
   print(f"ack_ok: {ack.ok}")
   print(f"ack_detail: {ack.detail}")
   print(f"recorded_packets: {len(client.sent)}")
   print(f"subscriptions: {', '.join(client.subscriptions)}")
   print(f"healthcheck: {client.healthcheck()}")

Output:

.. code-block:: text

   protocol: mqtt
   ack_ok: True
   ack_detail: dry-run packet recorded
   recorded_packets: 2
   subscriptions: pycsamt/WILLY-L18-TELEMETRY/+/+/qc
   healthcheck: True

Write And Replay JSONL Telemetry
--------------------------------

The file transport is the most useful offline transport. It writes one
packet dictionary per line, can be read back by another file client, and
does not require any optional dependencies.

.. code-block:: python
   :linenos:

   from pathlib import Path

   from pycsamt.iot import build_telemetry_client

   result_dir = Path("results")
   result_dir.mkdir(exist_ok=True)
   jsonl_path = result_dir / "iot-telemetry-demo.jsonl"

   if jsonl_path.exists():
       jsonl_path.unlink()

   with build_telemetry_client(
       "file",
       endpoint=str(jsonl_path),
       dry_run=False,
       append=False,
   ) as file_client:
       for packet in packets:
           ack = file_client.send(packet)
       print(f"protocol: {ack.protocol}")
       print(f"ack_ok: {ack.ok}")
       print(f"detail: {ack.detail}")
       print(f"sent_packets: {len(file_client.sent)}")

   reader = build_telemetry_client(
       "file",
       endpoint=str(jsonl_path),
       dry_run=False,
   )
   rows = reader.read_all()

   print(f"jsonl_rows: {len(rows)}")
   print(f"first_row_kind: {rows[0]['kind']}")
   print(f"first_row_station: {rows[0]['payload']['station']}")

Output:

.. code-block:: text

   protocol: file
   ack_ok: True
   detail: appended to results\iot-telemetry-demo.jsonl
   sent_packets: 5
   jsonl_rows: 5
   first_row_kind: health
   first_row_station: L18-001

Inspect Packet Tables
---------------------

Use :func:`pycsamt.iot.packet_table` for a flat packet inventory and
:func:`pycsamt.iot.telemetry_summary` for counts by device/topic. These
tables are useful before monitoring because they reveal missing packet
kinds, wrong topics, or devices that stopped reporting.

.. code-block:: python
   :linenos:

   from pycsamt.iot import packet_table, telemetry_summary

   packet_df = packet_table(packets)
   print(
       packet_df[
           ["device_id", "kind", "timestamp", "payload_keys"]
       ].to_string(index=False)
   )

   print()

   summary = telemetry_summary(packets)
   print(summary[["device_id", "topic", "n_packet"]].to_string(index=False))

Output:

.. code-block:: text

     device_id   kind    timestamp                                                                        payload_keys
   l18-node-01 health 1700000000.0                                            battery_v;firmware;station;temperature_c
   l18-node-01     qc 1700000006.0         accepted;channels;decision;finite_coverage;frequency_band_hz;method;station
   l18-node-02 health 1700000012.0                                            battery_v;firmware;station;temperature_c
   l18-node-02     qc 1700000018.0 accepted;channels;decision;finite_coverage;frequency_band_hz;method;reasons;station
   l18-node-02   sync 1700000024.0                                           clock_offset_ms;drift_ppm;quality;station

     device_id                                                  topic  n_packet
   l18-node-01 pycsamt/WILLY-L18-TELEMETRY/L18-001/l18-node-01/health         1
   l18-node-01     pycsamt/WILLY-L18-TELEMETRY/L18-001/l18-node-01/qc         1
   l18-node-02 pycsamt/WILLY-L18-TELEMETRY/L18-002/l18-node-02/health         1
   l18-node-02     pycsamt/WILLY-L18-TELEMETRY/L18-002/l18-node-02/qc         1
   l18-node-02   pycsamt/WILLY-L18-TELEMETRY/L18-002/l18-node-02/sync         1

Assess Telemetry In A Session
-----------------------------

Add packets to a :class:`pycsamt.iot.FieldSession` when you want telemetry
to interact with station metadata, monitoring thresholds, and dashboards.
The example below is intentionally synthetic and deterministic: station
``L18-002`` rejects one QC packet, so the session is marked critical.

.. code-block:: python
   :linenos:

   from pycsamt.iot import (
       FieldSession,
       MonitoringConfig,
       StationConfig,
   )

   stations = [
       StationConfig(
           "L18-001",
           profile="L18",
           position_m=0.0,
           lat=7.501,
           lon=-5.201,
           channels=["ex", "ey", "hx", "hy"],
       ),
       StationConfig(
           "L18-002",
           profile="L18",
           position_m=50.0,
           lat=7.502,
           lon=-5.198,
           channels=["ex", "ey", "hx", "hy"],
       ),
   ]

   session = FieldSession(
       "WILLY-L18-TELEMETRY",
       devices=devices,
       stations=stations,
       monitoring_config=MonitoringConfig(
           method="amt",
           expected_interval_s=6.0,
           max_gap_s=18.0,
           min_edge_acceptance_rate=0.75,
           min_battery_v=11.2,
           required_channels=["ex", "ey", "hx", "hy"],
       ),
   )
   session.add_packets(packets)
   status = session.assess(now=1_700_000_030.0)

   print(f"level: {status.level.value}")
   print(f"n_packet: {status.n_packet}")
   print(f"packet_success_rate: {status.packet_success_rate:.3f}")
   print(f"edge_acceptance_rate: {status.edge_acceptance_rate:.3f}")
   print(f"max_gap_s: {status.max_gap_s:.1f}")
   print(f"issues: {', '.join(status.issues) or '-'}")

Output:

.. code-block:: text

   level: critical
   n_packet: 5
   packet_success_rate: 1.000
   edge_acceptance_rate: 0.500
   max_gap_s: 6.0
   issues: edge_acceptance_rate_below_threshold

Plot A Telemetry Dashboard
--------------------------

Use :func:`pycsamt.iot.plot_field_dashboard` to visualise the same
session. This figure combines station status, edge-QC acceptance, power
or synchronisation status, and packet timing.

.. code-block:: python
   :linenos:

   from pathlib import Path

   from pycsamt.iot import plot_field_dashboard

   out_dir = Path("docs/source/images/user_guide/iot")
   out_dir.mkdir(parents=True, exist_ok=True)

   fig = plot_field_dashboard(
       session,
       now=1_700_000_030.0,
       station_axis="profile",
       title="Telemetry transport demo: WILLY-L18",
   )
   fig.savefig(
       out_dir / "user-guide-iot-telemetry-01.png",
       dpi=180,
   )

.. image:: ../../images/user_guide/iot/user-guide-iot-telemetry-01.png
   :width: 100%

Transport Choices
-----------------

The transport name controls the client returned by
:func:`pycsamt.iot.build_telemetry_client`.

.. list-table::
   :header-rows: 1
   :widths: 18 34 48

   * - Protocol
     - Client
     - Typical use
   * - ``file``
     - :class:`pycsamt.iot.FileTelemetryClient`
     - Offline logging, replay, tests, documentation.
   * - ``mqtt``
     - :class:`pycsamt.iot.MQTTTelemetryClient`
     - Field broker publishing and topic subscriptions.
   * - ``http``
     - :class:`pycsamt.iot.HTTPTelemetryClient`
     - Gateway POST requests to an API endpoint.
   * - ``serial``
     - :class:`pycsamt.iot.SerialTelemetryClient`
     - UART links to local acquisition hardware.
   * - ``websocket``
     - :class:`pycsamt.iot.WebSocketTelemetryClient`
     - Live dashboard streams.
   * - ``lora``
     - :class:`pycsamt.iot.TelemetryClient`
     - Dry-run recorder unless a project-specific transport is added.

MQTT, serial, and WebSocket use optional third-party libraries only when a
real connection is opened. Keep ``dry_run=True`` in notebooks, CI jobs,
and documentation examples unless you intentionally want to contact a
field endpoint.

Surviving an intermittent uplink
--------------------------------

Remote nodes lose their link. Wrap any client in a
:class:`~pycsamt.iot.StoreAndForwardClient` so a failed send queues the
packet instead of dropping it; :meth:`flush` later drains the backlog in
order, and an optional spool file carries it across a restart:

.. code-block:: python

   from pycsamt.iot import StoreAndForwardClient, build_telemetry_client

   client = build_telemetry_client("http", endpoint="https://hub/ingest",
                                   dry_run=False)
   buffered = StoreAndForwardClient(client, spool_path="spool.jsonl",
                                    max_queue=5000)

   buffered.send(packet)        # delivered, or queued if the link is down
   # ... later, when connectivity returns ...
   buffered.flush()             # drains the queue in order (at-least-once)

Delivery is at-least-once and order-preserving, ``max_queue`` bounds the
buffer with a drop-oldest policy, and ``next_retry_delay_s`` provides an
exponential-backoff hint for a scheduling loop.

Field Interpretation
--------------------

Telemetry tables answer operational questions before geophysical
processing starts: which node reported, what kind of packet was sent,
which station it belongs to, and whether the payload contains health,
QC, sync, or power information. Monitoring then turns those packets into
a survey status. In this example the packet stream is complete, but the
edge acceptance rate is only ``0.500``, so the session is flagged
``critical`` even though packet delivery itself is healthy.

For real surveys, keep topic names stable across the acquisition campaign
and persist file logs or broker exports with the processing report. That
makes later QC, provenance, and troubleshooting much easier to audit.
