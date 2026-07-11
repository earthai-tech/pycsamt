.. _user_guide_iot_visualization:

Visualization
=============

The IoT plotting layer lives in :mod:`pycsamt.iot.plot` and is re-exported
from :mod:`pycsamt.iot`. These figures are operational acquisition plots:
they show telemetry health, edge-QC decisions, power budget, and clock
synchronisation before the data are converted into impedance products or
inversion inputs.

The plotting functions accept high-level objects such as
:class:`pycsamt.iot.FieldSession`, telemetry packets, energy estimates,
sync status objects, or serialised mappings. Each returned Matplotlib
figure also carries the normalised data used to draw the panels in a
``fig.pycsamt_iot_*`` attribute. This makes report generation auditable:
the image and the plotted rows stay together.

Build A Visualisation Session
-----------------------------

This example is synthetic and deterministic. It creates three stations on
profile ``L18``: one healthy station, one station with poor finite
coverage, and one station with spike and timing problems. The goal is to
exercise all visual panels, not to represent a real field deployment.

.. code-block:: python
   :linenos:

   import numpy as np

   from pycsamt.iot import (
       DeviceConfig,
       EdgeProcessingConfig,
       EdgeProcessor,
       EnergyConfig,
       FieldSession,
       MonitoringConfig,
       PacketKind,
       StationConfig,
       TelemetryPacket,
       estimate_energy_budget,
   )

   survey_id = "WILLY-L18-VISUAL"
   devices = [
       DeviceConfig(
           "l18-node-01",
           station="L18-001",
           channels=["ex", "ey", "hx", "hy"],
           sample_rate_hz=256.0,
       ),
       DeviceConfig(
           "l18-node-02",
           station="L18-002",
           channels=["ex", "ey", "hx", "hy"],
           sample_rate_hz=256.0,
       ),
       DeviceConfig(
           "l18-node-03",
           station="L18-003",
           channels=["ex", "ey", "hx", "hy"],
           sample_rate_hz=256.0,
       ),
   ]
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
       StationConfig(
           "L18-003",
           profile="L18",
           position_m=100.0,
           lat=7.503,
           lon=-5.195,
           channels=["ex", "ey", "hx", "hy"],
       ),
   ]

   session = FieldSession(
       survey_id,
       devices=devices,
       stations=stations,
       monitoring_config=MonitoringConfig(
           method="amt",
           expected_interval_s=6.0,
           max_gap_s=18.0,
           min_edge_acceptance_rate=0.70,
           min_battery_v=11.2,
           required_channels=["ex", "ey", "hx", "hy"],
       ),
   )

   processor = EdgeProcessor(
       EdgeProcessingConfig(
           decimation=1,
           finite_threshold=0.85,
           warn_finite_threshold=0.95,
           channel_names=["ex", "ey", "hx", "hy"],
           spike_threshold=4.0,
           max_spike_fraction=0.12,
       )
   )

   rng = np.random.default_rng(42)
   for idx, device in enumerate(devices):
       data = rng.normal(size=(256, 4))
       if idx == 1:
           data[20:230, 1] = np.nan
       if idx == 2:
           data[80:125, 0] = 15.0

       edge = processor.process(data)
       qc_packet = edge.to_packet(
           device,
           timestamp=1_700_000_000.0 + idx * 6.0,
           survey_id=survey_id,
           qos=1,
       )
       qc_packet.payload["station"] = device.station
       qc_packet.payload["method"] = "amt"
       session.add_packet(qc_packet)

       health = TelemetryPacket.from_device(
           device,
           timestamp=1_700_000_002.0 + idx * 6.0,
           survey_id=survey_id,
           kind=PacketKind.HEALTH,
           payload={
               "station": device.station,
               "battery_v": [12.55, 11.84, 10.92][idx],
               "temperature_c": [30.1, 31.4, 33.2][idx],
           },
       )
       session.add_packet(health)

       power_config = EnergyConfig(
           battery_wh=[180.0, 120.0, 70.0][idx],
           active_power_w=[1.2, 1.6, 2.1][idx],
           sleep_power_w=0.18,
           duty_cycle=[0.25, 0.35, 0.55][idx],
           solar_wh_per_day=[8.0, 3.0, 0.0][idx],
           telemetry_power_w=2.0,
           telemetry_seconds_per_day=600.0,
           edge_power_w=0.35,
           edge_duty_cycle=0.20,
           min_runtime_days=5.0,
           device_id=device.device_id,
       )
       power_packet = estimate_energy_budget(power_config).to_packet(
           device,
           timestamp=1_700_000_004.0 + idx * 6.0,
           survey_id=survey_id,
       )
       power_packet.payload["station"] = device.station
       session.add_packet(power_packet)

   sync_payloads = [
       {
           "station": "L18-001",
           "offset_ms": 0.32,
           "drift_ppm": 1.4,
           "jitter_ms": 0.16,
           "gps_lock": True,
           "n_reference_points": 120,
           "quality": "excellent",
       },
       {
           "station": "L18-002",
           "offset_ms": 1.84,
           "drift_ppm": 11.2,
           "jitter_ms": 0.72,
           "gps_lock": True,
           "n_reference_points": 118,
           "quality": "fair",
       },
       {
           "station": "L18-003",
           "offset_ms": 7.62,
           "drift_ppm": 69.5,
           "jitter_ms": 2.15,
           "gps_lock": False,
           "n_reference_points": 61,
           "quality": "poor",
       },
   ]
   for idx, payload in enumerate(sync_payloads):
       session.add_packet(
           TelemetryPacket.from_device(
               devices[idx],
               timestamp=1_700_000_006.0 + idx * 6.0,
               survey_id=survey_id,
               kind=PacketKind.SYNC,
               payload=payload,
           )
       )

   print(f"survey_id: {session.survey_id}")
   print(f"n_devices: {session.n_devices}")
   print(f"n_stations: {session.n_stations}")
   print(f"n_packets: {session.n_packets}")

Output:

.. code-block:: text

   survey_id: WILLY-L18-VISUAL
   n_devices: 3
   n_stations: 3
   n_packets: 12

Plot The Field Dashboard
------------------------

Use :func:`pycsamt.iot.plot_field_dashboard` for an at-a-glance field
status. The four panels show station health, edge-QC acceptance, power or
synchronisation state, and packet timing. Set ``station_axis="profile"``
for line work and ``station_axis="map"`` when all stations have valid
coordinates.

.. code-block:: python
   :linenos:

   from pathlib import Path

   from pycsamt.iot import plot_field_dashboard

   out_dir = Path("docs/source/images/user_guide/iot")
   out_dir.mkdir(parents=True, exist_ok=True)

   dashboard = plot_field_dashboard(
       session,
       now=1_700_000_030.0,
       station_axis="profile",
       title="IoT field dashboard: WILLY L18",
   )
   dashboard.savefig(
       out_dir / "user-guide-iot-visualization-01.png",
       dpi=180,
   )

   dashboard_data = dashboard.pycsamt_iot_dashboard
   print(f"dashboard_stations: {len(dashboard_data['stations'])}")
   print(f"dashboard_packets: {len(dashboard_data['packets'])}")
   print(f"monitoring_level: {dashboard_data['monitoring']['level']}")
   print(f"issues: {', '.join(dashboard_data['issues']) or '-'}")

Output:

.. code-block:: text

   dashboard_stations: 3
   dashboard_packets: 12
   monitoring_level: critical
   issues: battery_below_threshold, edge_acceptance_rate_below_threshold, required_channels_missing

Plot Edge-QC Detail
-------------------

Use :func:`pycsamt.iot.plot_edge_qc_summary` when the dashboard indicates
that edge quality is the problem. The function accepts a
:class:`pycsamt.iot.FieldSession`, one or more QC telemetry packets, or
raw :class:`pycsamt.iot.EdgeProcessingResult` objects.

.. code-block:: python
   :linenos:

   from pycsamt.iot import plot_edge_qc_summary

   qc_fig = plot_edge_qc_summary(
       session,
       title="Edge QC visual summary: WILLY L18",
   )
   qc_fig.savefig(
       out_dir / "user-guide-iot-visualization-02.png",
       dpi=180,
   )

   qc_rows = qc_fig.pycsamt_iot_edge_qc
   print(f"qc_rows: {len(qc_rows)}")
   print(f"qc_decisions: {sorted({row['decision'] for row in qc_rows})}")
   print(f"qc_channels: {sorted({row['channel'] for row in qc_rows})}")

Output:

.. code-block:: text

   qc_rows: 12
   qc_decisions: ['accept', 'reject']
   qc_channels: ['ex', 'ey', 'hx', 'hy']

Plot Power Budgets
------------------

Use :func:`pycsamt.iot.plot_power_budget` to compare daily load, harvest,
runtime, and power states. The input can be a session containing power
packets, a list of :class:`pycsamt.iot.EnergyConfig` objects, power
telemetry packets, or energy estimates.

.. code-block:: python
   :linenos:

   import numpy as np

   from pycsamt.iot import plot_power_budget

   power_fig = plot_power_budget(
       session,
       title="Power budget visual summary: WILLY L18",
   )
   power_fig.savefig(
       out_dir / "user-guide-iot-visualization-03.png",
       dpi=180,
   )

   power_rows = power_fig.pycsamt_iot_power_budget
   print(f"power_rows: {len(power_rows)}")
   print(f"power_states: {sorted({row['state'] for row in power_rows})}")
   print(
       "runtime_days: "
       + ", ".join(
           "inf" if np.isinf(row["runtime_days"])
           else f"{row['runtime_days']:.2f}"
           for row in power_rows
       )
   )

Output:

.. code-block:: text

   power_rows: 3
   power_states: ['critical', 'ok']
   runtime_days: 40.42, 7.86, 2.21

Plot Clock Synchronisation
--------------------------

Use :func:`pycsamt.iot.plot_sync_quality` to inspect offset, drift, jitter,
reference support, GPS lock, and quality grades. Threshold arguments draw
reference lines but do not mutate the data.

.. code-block:: python
   :linenos:

   from pycsamt.iot import plot_sync_quality

   sync_fig = plot_sync_quality(
       session,
       title="Clock sync visual summary: WILLY L18",
       tolerance_ms=1.0,
       max_drift_ppm=10.0,
       max_jitter_ms=1.0,
   )
   sync_fig.savefig(
       out_dir / "user-guide-iot-visualization-04.png",
       dpi=180,
   )

   sync_rows = sync_fig.pycsamt_iot_sync_quality
   print(f"sync_rows: {len(sync_rows)}")
   print(f"sync_quality: {sorted({row['quality'] for row in sync_rows})}")
   print(f"gps_lock_values: {[row['gps_lock'] for row in sync_rows]}")

Output:

.. code-block:: text

   sync_rows: 3
   sync_quality: ['excellent', 'fair', 'poor']
   gps_lock_values: [True, True, False]

Generated Figures
-----------------

The figures are displayed in a two-column grid so the page remains compact
while still showing each diagnostic family clearly.

.. grid:: 1 1 2 2

   .. grid-item::

      .. image:: ../../images/user_guide/iot/user-guide-iot-visualization-01.png
         :width: 100%

   .. grid-item::

      .. image:: ../../images/user_guide/iot/user-guide-iot-visualization-02.png
         :width: 100%

   .. grid-item::

      .. image:: ../../images/user_guide/iot/user-guide-iot-visualization-03.png
         :width: 100%

   .. grid-item::

      .. image:: ../../images/user_guide/iot/user-guide-iot-visualization-04.png
         :width: 100%

Choosing The Right Plot
-----------------------

.. list-table::
   :header-rows: 1
   :widths: 28 32 40

   * - Function
     - Best input
     - Use it for
   * - :func:`pycsamt.iot.plot_field_dashboard`
     - :class:`pycsamt.iot.FieldSession`
     - Survey status, station health, packet timing, acceptance rates.
   * - :func:`pycsamt.iot.plot_edge_qc_summary`
     - Session, QC packets, or edge results
     - Rejection reasons, finite coverage, spike fractions.
   * - :func:`pycsamt.iot.plot_power_budget`
     - Session, power packets, energy configs, or estimates
     - Runtime, harvest/load balance, power states.
   * - :func:`pycsamt.iot.plot_sync_quality`
     - Session, sync packets, or sync status rows
     - Offset, drift, jitter, GPS lock, synchronisation grades.

Field Interpretation
--------------------

Use the dashboard first when reviewing a live or replayed acquisition
session. If the dashboard reports a telemetry issue, inspect packet tables
and monitoring summaries. If it reports edge acceptance problems, move to
the QC summary. If a station is losing voltage or runtime margin, use the
power plot. If packet timing or GPS lock is suspicious, use the sync plot.

These figures do not replace MT/AMT/CSAMT response plots. They explain the
condition of the acquisition system that produced the data, which is
exactly the context needed before deciding whether a station should enter
the geophysical processing workflow.
