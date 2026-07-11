.. _user_guide_iot_simulation:

Simulation
==========

The IoT simulator creates deterministic field-like data for documentation,
tests, demos, and pipeline development when hardware is not available. It
can generate AMT channel windows, station/device inventories, telemetry
packets, packet loss, GPS clock drift, and battery decay. Every simulator
accepts a ``seed`` argument so examples can be reproduced exactly.

These examples are synthetic by design. They are not substitutes for real
EDI or logger data; they are controlled inputs for testing the IoT layer
around acquisition.

Simulate AMT Channels
---------------------

Use :func:`pycsamt.iot.simulate_amt_channel` to create one live channel
window. The signal includes band-limited AMT-like energy, Gaussian noise,
optional powerline harmonics, and optional dropout gaps.

.. code-block:: python
   :linenos:

   import numpy as np

   from pycsamt.iot import (
       EdgeProcessingConfig,
       EdgeProcessor,
       detect_powerline_harmonics,
       estimate_channel_snr,
       simulate_amt_channel,
   )

   sample_rate = 256.0
   n_samples = 4096

   ex = simulate_amt_channel(
       n_samples,
       sample_rate,
       snr_db=18.0,
       mains_hz=50.0,
       powerline_amplitude=0.18,
       dropout_rate=0.015,
       seed=42,
   )
   ey = simulate_amt_channel(
       n_samples,
       sample_rate,
       snr_db=12.0,
       mains_hz=50.0,
       powerline_amplitude=0.06,
       dropout_rate=0.04,
       seed=43,
   )

   window = np.column_stack([ex, ey])
   harmonics = detect_powerline_harmonics(
       ex,
       sample_rate,
       mains_hz=50.0,
       threshold_ratio=0.02,
   )
   snr_ex = estimate_channel_snr(
       ex,
       sample_rate,
       signal_band_hz=(2.0, 40.0),
   )
   processor = EdgeProcessor(
       EdgeProcessingConfig(
           decimation=4,
           finite_threshold=0.90,
           warn_finite_threshold=0.98,
           channel_names=["ex", "ey"],
           spike_threshold=5.0,
       )
   )
   qc = processor.process(window)

   print(f"n_samples: {n_samples}")
   print(f"sample_rate_hz: {sample_rate:.1f}")
   print(f"finite_coverage: {qc.metrics['finite_coverage']:.3f}")
   print(f"qc_decision: {qc.decision.value}")
   print(f"ex_snr_2_40_hz_db: {snr_ex:.2f}")
   print(f"powerline_contaminated: {harmonics.contaminated}")
   print(f"powerline_total_ratio: {harmonics.total_ratio:.3f}")

Output:

.. code-block:: text

   n_samples: 4096
   sample_rate_hz: 256.0
   finite_coverage: 0.975
   qc_decision: warning
   ex_snr_2_40_hz_db: -6.76
   powerline_contaminated: True
   powerline_total_ratio: 0.371

Simulate One Station
--------------------

Use :func:`pycsamt.iot.simulate_amt_station` when you need station
metadata, a device config, channel arrays, and basic health/QC packets in
one object.

.. code-block:: python
   :linenos:

   from pycsamt.iot import simulate_amt_station

   station = simulate_amt_station(
       "L18-S001",
       sample_rate=256.0,
       n_samples=1024,
       snr_db=16.0,
       powerline_amplitude=0.1,
       dropout_rate=0.02,
       survey_id="SIM-L18",
       profile="L18",
       position_m=0.0,
       seed=8,
   )

   print(f"station_id: {station['station'].station_id}")
   print(f"device_id: {station['device'].device_id}")
   print(f"channels: {', '.join(station['station'].channels)}")
   print(f"packets: {len(station['packets'])}")
   print(f"qc_decision: {station['packets'][1].payload['decision']}")

Output:

.. code-block:: text

   station_id: L18-S001
   device_id: node-L18-S001
   channels: ex, ey, hx, hy
   packets: 2
   qc_decision: accept

Simulate A Network
------------------

Use :func:`pycsamt.iot.simulate_iot_network` to create many stations
across one or more profiles. With ``detail=True``, the return value
contains station dictionaries and a flat packet list. Use
:func:`pycsamt.iot.simulate_packet_loss` to drop a reproducible fraction of
packets.

.. code-block:: python
   :linenos:

   from pycsamt.iot import (
       packet_table,
       simulate_iot_network,
       simulate_packet_loss,
       telemetry_summary,
   )

   network = simulate_iot_network(
       n_stations=6,
       profiles=["L18", "L22"],
       sample_rate=128.0,
       n_samples=512,
       snr_db=14.0,
       dropout_rate=0.03,
       survey_id="SIM-WILLY",
       station_spacing_m=50.0,
       seed=9,
       detail=True,
   )
   packets = network["packets"]
   kept = simulate_packet_loss(packets, dropout_rate=0.25, seed=10)

   print(f"original_packets: {len(packets)}")
   print(f"after_packet_loss: {len(kept)}")
   print(
       packet_table(kept)[
           ["device_id", "kind", "timestamp", "payload_keys"]
       ].head(6).to_string(index=False)
   )
   print()
   print(
       telemetry_summary(kept)[
           ["device_id", "topic", "n_packet"]
       ].head(6).to_string(index=False)
   )

Output:

.. code-block:: text

   original_packets: 12
   after_packet_loss: 9
       device_id   kind    timestamp                                                                payload_keys
   node-L18-S001 health 1700000000.0                                    battery_v;firmware;station;temperature_c
   node-L22-S001 health 1700000005.0                                    battery_v;firmware;station;temperature_c
   node-L18-S002 health 1700000010.0                                    battery_v;firmware;station;temperature_c
   node-L22-S002 health 1700000015.0                                    battery_v;firmware;station;temperature_c
   node-L22-S002     qc 1700000019.0 accepted;channels;decision;finite_coverage;frequency_band_hz;method;station
   node-L18-S003 health 1700000020.0                                    battery_v;firmware;station;temperature_c

       device_id                                           topic  n_packet
   node-L18-S001 pycsamt/SIM-WILLY/L18-S001/node-L18-S001/health         1
   node-L18-S002 pycsamt/SIM-WILLY/L18-S002/node-L18-S002/health         1
   node-L18-S003 pycsamt/SIM-WILLY/L18-S003/node-L18-S003/health         1
   node-L18-S003     pycsamt/SIM-WILLY/L18-S003/node-L18-S003/qc         1
   node-L22-S001 pycsamt/SIM-WILLY/L22-S001/node-L22-S001/health         1
   node-L22-S002 pycsamt/SIM-WILLY/L22-S002/node-L22-S002/health         1

Monitor Simulated Packets
-------------------------

Simulated packets can be fed directly into the monitoring layer. This is
useful when testing dashboard behavior or packet-loss handling before
hardware exists.

.. code-block:: python
   :linenos:

   from pycsamt.iot import MonitoringConfig, TelemetryMonitor

   monitor = TelemetryMonitor(
       MonitoringConfig(
           method="amt",
           expected_interval_s=5.0,
           max_gap_s=20.0,
           min_packet_success_rate=0.90,
           min_edge_acceptance_rate=0.80,
           required_channels=["ex", "ey", "hx", "hy"],
       )
   )
   status = monitor.assess(kept, now=1_700_000_080.0)

   print(f"level: {status.level.value}")
   print(f"n_packet: {status.n_packet}")
   print(f"edge_acceptance_rate: {status.edge_acceptance_rate:.3f}")
   print(f"max_gap_s: {status.max_gap_s:.1f}")
   print(f"issues: {', '.join(status.issues) or '-'}")

Output:

.. code-block:: text

   level: warning
   n_packet: 9
   edge_acceptance_rate: 1.000
   max_gap_s: 5.0
   issues: latency_above_threshold

Simulate GPS Drift And Battery Decay
------------------------------------

Use :func:`pycsamt.iot.simulate_gps_drift` to generate paired local and
reference timestamps for synchronisation tests. Use
:func:`pycsamt.iot.simulate_battery_decay` for power and monitoring demos.

.. code-block:: python
   :linenos:

   from pycsamt.iot import (
       ClockSynchronizer,
       simulate_battery_decay,
       simulate_gps_drift,
   )
   import numpy as np

   gps = simulate_gps_drift(
       120,
       sample_interval_s=1.0,
       drift_ppm=8.0,
       jitter_ms=0.25,
       offset_ms=0.7,
       dropout_rate=0.1,
       seed=11,
   )
   sync = ClockSynchronizer().assess(
       "l18-node-01",
       gps["local"],
       gps["reference"],
       gps_lock=bool(np.mean(gps["gps_lock"]) > 0.9),
   )
   battery = simulate_battery_decay(
       120,
       initial_v=13.1,
       final_v=10.8,
       seed=12,
   )

   print(f"gps_lock_fraction: {np.mean(gps['gps_lock']):.3f}")
   print(f"offset_ms: {sync.offset_ms:.3f}")
   print(f"drift_ppm: {sync.drift_ppm:.3f}")
   print(f"jitter_ms: {sync.jitter_ms:.3f}")
   print(f"sync_quality: {sync.quality.value}")
   print(f"battery_start_v: {battery[0]:.2f}")
   print(f"battery_end_v: {battery[-1]:.2f}")

Output:

.. code-block:: text

   gps_lock_fraction: 0.892
   offset_ms: 1.214
   drift_ppm: 9.354
   jitter_ms: 0.254
   sync_quality: fair
   battery_start_v: 13.10
   battery_end_v: 11.06

Plot Simulation Outputs
-----------------------

The first figure shows channel time series and a spectrum with simulated
mains contamination. The second shows the simulated network, packet counts
after loss, clock drift, and battery decay.

.. code-block:: python
   :linenos:

   from pathlib import Path

   import matplotlib.pyplot as plt
   import numpy as np

   out_dir = Path("docs/source/images/user_guide/iot")
   out_dir.mkdir(parents=True, exist_ok=True)

   t = np.arange(n_samples) / sample_rate
   freq = np.fft.rfftfreq(n_samples, d=1.0 / sample_rate)
   spec = np.abs(np.fft.rfft(np.nan_to_num(ex))) ** 2

   fig, axes = plt.subplots(
       2,
       1,
       figsize=(9.0, 6.2),
       constrained_layout=True,
   )
   axes[0].plot(t[:1200], ex[:1200], lw=0.9, label="Ex")
   axes[0].plot(t[:1200], ey[:1200], lw=0.9, alpha=0.75, label="Ey")
   axes[0].set_title("Simulated AMT channels")
   axes[0].set_xlabel("Time (s)")
   axes[0].set_ylabel("Amplitude")
   axes[0].legend(loc="upper right")
   axes[0].grid(alpha=0.25)

   axes[1].semilogy(freq[1:], spec[1:], color="tab:blue")
   for peak in harmonics.peaks:
       if peak.flagged:
           axes[1].axvline(
               peak.frequency_hz,
               color="tab:red",
               ls="--",
               lw=1.0,
           )
   axes[1].set_xlim(1.0, 128.0)
   axes[1].set_title("Ex spectrum with simulated mains harmonics")
   axes[1].set_xlabel("Frequency (Hz)")
   axes[1].set_ylabel("Power")
   axes[1].grid(alpha=0.25, which="both")
   fig.savefig(out_dir / "user-guide-iot-simulation-01.png", dpi=180)
   plt.close(fig)

   fig, axes = plt.subplots(
       2,
       2,
       figsize=(10.5, 7.0),
       constrained_layout=True,
   )
   profiles = [item["station"].profile for item in network["stations"]]
   positions = [item["station"].position_m for item in network["stations"]]
   accepted = [
       item["packets"][1].payload["accepted"]
       for item in network["stations"]
   ]
   colors = ["#2ca25f" if ok else "#de2d26" for ok in accepted]
   axes[0, 0].scatter(
       positions,
       profiles,
       c=colors,
       s=90,
       edgecolor="black",
   )
   axes[0, 0].set_title("Simulated station network")
   axes[0, 0].set_xlabel("Profile position (m)")
   axes[0, 0].set_ylabel("Profile")
   axes[0, 0].grid(alpha=0.25)

   kinds = {}
   for packet in kept:
       kinds[packet.kind.value] = kinds.get(packet.kind.value, 0) + 1
   axes[0, 1].bar(kinds.keys(), kinds.values(), color="#756bb1")
   axes[0, 1].set_title("Packets after simulated loss")
   axes[0, 1].set_ylabel("Count")
   axes[0, 1].grid(axis="y", alpha=0.25)

   elapsed = gps["reference"] - gps["reference"][0]
   offset_ms = (gps["local"] - gps["reference"]) * 1000.0
   axes[1, 0].plot(elapsed, offset_ms)
   axes[1, 0].scatter(
       elapsed[~gps["gps_lock"]],
       offset_ms[~gps["gps_lock"]],
       color="tab:red",
       s=18,
       label="GPS lost",
   )
   axes[1, 0].set_title("Simulated clock offset")
   axes[1, 0].set_xlabel("Elapsed time (s)")
   axes[1, 0].set_ylabel("Offset (ms)")
   axes[1, 0].legend(loc="upper left")
   axes[1, 0].grid(alpha=0.25)

   axes[1, 1].plot(np.arange(battery.size), battery, color="#2ca25f")
   axes[1, 1].axhline(
       11.0,
       ls="--",
       color="#de2d26",
       label="low battery",
   )
   axes[1, 1].set_title("Simulated battery decay")
   axes[1, 1].set_xlabel("Sample")
   axes[1, 1].set_ylabel("Voltage")
   axes[1, 1].legend(loc="upper right")
   axes[1, 1].grid(alpha=0.25)
   fig.savefig(out_dir / "user-guide-iot-simulation-02.png", dpi=180)
   plt.close(fig)

.. grid:: 1 1 2 2

   .. grid-item::

      .. image:: ../../images/user_guide/iot/user-guide-iot-simulation-01.png
         :width: 100%

   .. grid-item::

      .. image:: ../../images/user_guide/iot/user-guide-iot-simulation-02.png
         :width: 100%

Field Interpretation
--------------------

Simulation is useful because it makes edge cases reproducible. Here the
channel window has enough missing samples to trigger a warning, a strong
powerline signature, and a controlled SNR. The network example shows how
packet loss changes packet inventories while still leaving enough QC
packets for monitoring. The GPS and battery curves provide deterministic
inputs for synchronisation and power-management documentation.

For production work, clearly label simulated data in reports and examples.
Use simulation to test the IoT workflow, then replace simulated packets
with real telemetry when field hardware is available.
