.. _user_guide_iot_monitoring:

Monitoring
==========

Telemetry monitoring turns IoT packets into an operational survey status.
It checks whether packets are arriving, whether packet acknowledgements
are healthy, whether edge QC is accepting enough windows, and whether
field thresholds for latency, packet gaps, battery, clock offset, required
channels, and frequency coverage are being respected.

The example below uses synthetic telemetry packets for three stations on
the L18 demo line. Synthetic packets are appropriate here because the
monitor operates on live IoT messages, not EDI impedance files. The packet
stream is deliberately mixed: one station is healthy, one has a rejected
edge-QC window and an acknowledgement failure, and one has clock and
frequency-band issues.

Create A Monitored Session
--------------------------

Start with devices and a :class:`pycsamt.iot.MonitoringConfig`. The config
records the operational contract for the stream.

.. code-block:: python
   :linenos:

   from pycsamt.iot import DeviceConfig, FieldSession, MonitoringConfig

   base_time = 1_700_000_000.0
   devices = [
       DeviceConfig(
           "l18-node-01",
           station="001A",
           protocol="file",
           sample_rate_hz=512.0,
           channels=["ex", "ey", "hx", "hy"],
           role="recorder",
       ),
       DeviceConfig(
           "l18-node-02",
           station="002U",
           protocol="file",
           sample_rate_hz=512.0,
           channels=["ex", "ey", "hx", "hy"],
           role="recorder",
       ),
       DeviceConfig(
           "l18-node-03",
           station="003A",
           protocol="file",
           sample_rate_hz=512.0,
           channels=["ex", "ey", "hx", "hy"],
           role="recorder",
       ),
   ]

   config = MonitoringConfig(
       method="amt",
       expected_interval_s=60.0,
       max_gap_s=120.0,
       max_latency_s=10.0,
       min_packet_success_rate=0.95,
       min_edge_acceptance_rate=0.80,
       min_battery_v=11.2,
       max_clock_offset_ms=5.0,
       required_channels=["ex", "ey", "hx", "hy"],
       frequency_band_hz=(1.0, 1000.0),
   )

   session = FieldSession(
       "WILLY-L18-MONITORING-DEMO",
       devices=devices,
       method="amt",
       monitoring_config=config,
   )

Add Synthetic Telemetry Packets
-------------------------------

Each packet is a ``qc`` message with fields that the monitor knows how to
enrich: station, method, channel list, accepted/rejected decision,
acknowledgement status, latency, battery voltage, clock offset, and
frequency band.

.. code-block:: python
   :linenos:

   from pycsamt.iot import TelemetryPacket

   specs = [
       (devices[0], 0.0, True, True, 2.1, 12.4, 0.8,
        ["ex", "ey", "hx", "hy"], [1.0, 1000.0]),
       (devices[0], 60.0, True, True, 2.5, 12.3, 0.7,
        ["ex", "ey", "hx", "hy"], [1.0, 1000.0]),
       (devices[0], 120.0, True, True, 3.0, 12.2, 0.9,
        ["ex", "ey", "hx", "hy"], [1.0, 1000.0]),
       (devices[1], 180.0, True, True, 4.2, 11.8, 1.6,
        ["ex", "ey", "hx", "hy"], [1.0, 1000.0]),
       (devices[1], 240.0, False, True, 13.5, 11.7, 2.0,
        ["ex", "ey", "hx", "hy"], [1.0, 1000.0]),
       (devices[1], 420.0, True, False, 9.8, 10.9, 3.2,
        ["ex", "ey", "hx", "hy"], [1.0, 1000.0]),
       (devices[2], 480.0, True, False, 8.5, 11.4, 7.5,
        ["ex", "hx", "hy"], [1.0, 1000.0]),
       (devices[2], 540.0, True, False, 6.0, 11.3, 8.1,
        ["ex", "hx", "hy"], [0.2, 1200.0]),
   ]

   for device, dt, accepted, ack_ok, latency, battery, clock, channels, band in specs:
       session.add_packet(
           TelemetryPacket.from_device(
               device,
               timestamp=base_time + dt,
               survey_id=session.survey_id,
               kind="qc",
               payload={
                   "method": "amt",
                   "station": device.station,
                   "channels": channels,
                   "frequency_band_hz": band,
                   "accepted": accepted,
                   "decision": "accept" if accepted else "reject",
                   "ack_ok": ack_ok,
                   "latency_s": latency,
                   "battery_v": battery,
                   "clock_offset_ms": clock,
               },
           )
       )

Packet Tables And Counts
------------------------

Use :func:`pycsamt.iot.packet_table` for raw packet inventory and
:func:`pycsamt.iot.telemetry_summary` for packet counts by device/topic.

.. code-block:: python
   :linenos:

   from pycsamt.iot import packet_table, telemetry_summary

   packets = packet_table(session.packets)
   print(
       packets[
           ["device_id", "kind", "timestamp", "payload_keys"]
       ].head(4).to_string(index=False)
   )
   print()

   summary = telemetry_summary(session.packets)
   print(summary[["device_id", "topic", "n_packet"]].to_string(index=False))

Output:

.. code-block:: text

     device_id kind    timestamp                                                                                           payload_keys
   l18-node-01   qc 1700000000.0 accepted;ack_ok;battery_v;channels;clock_offset_ms;decision;frequency_band_hz;latency_s;method;station
   l18-node-01   qc 1700000060.0 accepted;ack_ok;battery_v;channels;clock_offset_ms;decision;frequency_band_hz;latency_s;method;station
   l18-node-01   qc 1700000120.0 accepted;ack_ok;battery_v;channels;clock_offset_ms;decision;frequency_band_hz;latency_s;method;station
   l18-node-02   qc 1700000180.0 accepted;ack_ok;battery_v;channels;clock_offset_ms;decision;frequency_band_hz;latency_s;method;station

     device_id                                                 topic  n_packet
   l18-node-01 pycsamt/WILLY-L18-MONITORING-DEMO/001A/l18-node-01/qc         3
   l18-node-02 pycsamt/WILLY-L18-MONITORING-DEMO/002U/l18-node-02/qc         3
   l18-node-03 pycsamt/WILLY-L18-MONITORING-DEMO/003A/l18-node-03/qc         2

Inspect Enriched Rows
---------------------

The monitor normalises payload fields into analysis-ready columns. This is
the best view when you need to debug why a status became warning or
critical.

.. code-block:: python
   :linenos:

   from pycsamt.iot import TelemetryMonitor

   monitor = TelemetryMonitor(config)
   enriched = monitor.table(session.packets, now=base_time + 600.0)
   print(
       enriched[
           [
               "device_id",
               "station",
               "edge_accepted",
               "ack_ok",
               "latency_s",
               "battery_v",
               "clock_offset_ms",
               "frequency_min_hz",
               "frequency_max_hz",
           ]
       ].head(6).to_string(index=False)
   )

Output:

.. code-block:: text

     device_id station  edge_accepted  ack_ok  latency_s  battery_v  clock_offset_ms  frequency_min_hz  frequency_max_hz
   l18-node-01    001A           True    True        2.1       12.4              0.8               1.0            1000.0
   l18-node-01    001A           True    True        2.5       12.3              0.7               1.0            1000.0
   l18-node-01    001A           True    True        3.0       12.2              0.9               1.0            1000.0
   l18-node-02    002U           True    True        4.2       11.8              1.6               1.0            1000.0
   l18-node-02    002U          False    True       13.5       11.7              2.0               1.0            1000.0
   l18-node-02    002U           True   False        9.8       10.9              3.2               1.0            1000.0

Assess Stream Status
--------------------

The status level is one of ``ok``, ``warning``, ``critical``, or
``no_data``. Critical issues include packet success below threshold,
edge-acceptance failure, low battery, high clock offset, method mismatch,
or missing required channels.

.. code-block:: python
   :linenos:

   from pycsamt.iot import assess_telemetry, monitoring_status_table

   status = monitor.assess(session.packets, now=base_time + 600.0)
   status_df = monitoring_status_table(status)
   print(
       status_df[
           [
               "level",
               "n_packet",
               "packet_success_rate",
               "edge_acceptance_rate",
               "mean_latency_s",
               "max_gap_s",
               "battery_min_v",
               "clock_offset_max_ms",
               "issues",
           ]
       ].to_string(index=False)
   )

   status2 = assess_telemetry(
       session.packets,
       config=config,
       now=base_time + 600.0,
   )
   print(f"Convenience wrapper level: {status2.level.value}")

Output:

.. code-block:: text

      level  n_packet  packet_success_rate  edge_acceptance_rate  mean_latency_s  max_gap_s  battery_min_v  clock_offset_max_ms                                                                                                                                                                                     issues
   critical         8                0.625                 0.875             6.2      180.0           10.9                  8.1 battery_below_threshold;clock_offset_above_threshold;frequency_outside_configured_band;packet_gap_above_threshold;packet_gap_exceeds_expected_interval;packet_success_rate_below_threshold
   Convenience wrapper level: critical

Plot A Monitoring Audit
-----------------------

The figure below is built from the enriched monitor table and status
metrics. It shows which thresholds were violated without requiring the
reader to inspect every packet by hand.

.. code-block:: python
   :linenos:

   from pathlib import Path

   import matplotlib.pyplot as plt
   import numpy as np

   out_dir = Path("docs/source/images/user_guide/iot")
   out_dir.mkdir(parents=True, exist_ok=True)

   metrics = status.as_dict()
   fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.4),
                            constrained_layout=True)
   fig.suptitle(
       "Telemetry monitoring audit: WILLY-L18-MONITORING-DEMO",
       fontsize=14,
   )

   station_acceptance = (
       enriched.groupby("station")["edge_accepted"].mean().sort_index()
   )
   colors = [
       "#2ca25f" if value >= config.min_edge_acceptance_rate else "#de2d26"
       for value in station_acceptance
   ]
   axes[0, 0].bar(station_acceptance.index, station_acceptance.values,
                  color=colors)
   axes[0, 0].axhline(config.min_edge_acceptance_rate,
                      ls="--", color="#de2d26")
   axes[0, 0].set_ylim(0, 1.05)
   axes[0, 0].set_ylabel("Acceptance rate")
   axes[0, 0].set_title("Edge acceptance by station")

   summary_metrics = {
       "packet\nsuccess": metrics["packet_success_rate"],
       "edge\nacceptance": metrics["edge_acceptance_rate"],
   }
   axes[0, 1].bar(list(summary_metrics), list(summary_metrics.values()),
                  color=["#de2d26", "#2ca25f"])
   axes[0, 1].axhline(config.min_packet_success_rate, ls="--",
                      color="#756bb1", label="min packet success")
   axes[0, 1].axhline(config.min_edge_acceptance_rate, ls=":",
                      color="#de2d26", label="min edge acceptance")
   axes[0, 1].set_ylim(0, 1.05)
   axes[0, 1].set_title("Stream-level rates")
   axes[0, 1].legend(loc="lower left", fontsize=8)

   packet_minutes = (
       enriched["timestamp"] - enriched["timestamp"].min()
   ) / 60.0
   packet_colors = [
       "#2ca25f" if ok else "#de2d26"
       for ok in enriched["ack_ok"]
   ]
   axes[1, 0].scatter(
       packet_minutes,
       enriched["station"],
       c=packet_colors,
       s=70,
       edgecolor="black",
       linewidth=0.4,
   )
   axes[1, 0].set_xlabel("Minutes since first packet")
   axes[1, 0].set_ylabel("Station")
   axes[1, 0].set_title("Packet acknowledgements")

   ops = enriched.groupby("station").agg(
       battery_min_v=("battery_v", "min"),
       clock_offset_max_ms=(
           "clock_offset_ms",
           lambda s: float(np.max(np.abs(s))),
       ),
   )
   x = np.arange(len(ops))
   axes[1, 1].bar(x - 0.18, ops["battery_min_v"], width=0.36,
                  label="battery V")
   axes[1, 1].bar(x + 0.18, ops["clock_offset_max_ms"], width=0.36,
                  label="clock offset ms")
   axes[1, 1].axhline(config.min_battery_v, ls="--",
                      color="#de2d26", label="min battery")
   axes[1, 1].axhline(config.max_clock_offset_ms, ls=":",
                      color="#756bb1", label="max clock offset")
   axes[1, 1].set_xticks(x, ops.index)
   axes[1, 1].set_title("Power and clock thresholds")
   axes[1, 1].legend(loc="upper right", fontsize=8)

   fig.savefig(out_dir / "user-guide-iot-monitoring-01.png", dpi=180)

.. image:: ../../images/user_guide/iot/user-guide-iot-monitoring-01.png
   :width: 100%

Field Interpretation
--------------------

The stream is ``critical`` because packet acknowledgement success is below
the configured threshold, one gap is too long, the minimum battery voltage
falls below 11.2 V, the clock offset exceeds 5 ms, and one packet reports a
frequency band outside the configured AMT band. The edge-acceptance rate
is still above the stream-level threshold, but station ``002U`` is visibly
weaker than the others and should be reviewed.

In a live deployment, these status rows should be logged with the
acquisition manifest. They explain why a station was accepted, repeated,
or excluded before downstream impedance, dimensionality, or inversion
workflows begin.
