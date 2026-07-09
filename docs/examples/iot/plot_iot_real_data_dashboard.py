r"""
IoT dashboard from bundled AMT station files
============================================

This example wraps the real WILLY AMT EDI files distributed in
``data/AMT/WILLY_DATA/`` with an IoT field-session layer.  The station
inventory comes from the real EDI filenames and survey-line directories;
the telemetry packets below are a reproducible operational overlay used
to demonstrate what pyCSAMT records during an IoT-enabled deployment:
edge QC decisions, recorder health, clock synchronisation, and power
budget state.

That separation is deliberate.  Historical EDI files usually preserve the
processed transfer functions, not live MQTT/serial telemetry.  The IoT
layer adds a machine-readable acquisition audit that can travel beside
the EDI data and explain which stations were accepted, which devices were
stable, and whether timing and energy margins were field-ready.
"""

# %%
# 1. Build a field session from real survey files
# -----------------------------------------------
# The helper below selects a few stations from the real WILLY ``L18PLT``
# line.  The EDI files are not parsed here; they are used as the station
# inventory that an IoT recorder fleet would report against.

from __future__ import annotations

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from pycsamt.iot import (
    DeviceConfig,
    DeviceRole,
    EnergyConfig,
    FieldSession,
    MonitoringConfig,
    PacketKind,
    StationConfig,
    TelemetryPacket,
    estimate_energy_budget,
    plot_edge_qc_summary,
    plot_field_dashboard,
    plot_power_budget,
    plot_sync_quality,
)


def repo_root() -> Path:
    """Return the repository root during docs builds and local runs."""
    env_root = os.environ.get("PYCSAMT_DOCS_REPO_ROOT")
    if env_root:
        return Path(env_root)
    # docs/examples/iot/plot_iot_real_data_dashboard.py -> repo root
    return Path(__file__).resolve().parents[3]


def real_willy_station_files(max_stations: int = 8) -> list[Path]:
    """Return real bundled WILLY AMT EDI files for one profile."""
    line_dir = repo_root() / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
    files = sorted(line_dir.glob("*.edi"))
    if not files:
        raise FileNotFoundError(f"No bundled EDI files found in {line_dir}")
    return files[:max_stations]


edi_files = real_willy_station_files()
sizes = np.asarray([path.stat().st_size for path in edi_files], dtype=float)
size_scale = max(float(np.nanmax(sizes)), 1.0)

devices = []
stations = []
for index, path in enumerate(edi_files, start=1):
    station_id = path.stem
    device_id = f"willy-node-{index:02d}"
    rel_path = path.relative_to(repo_root()).as_posix()
    channels = ["ex", "ey", "hx", "hy"]

    devices.append(
        DeviceConfig(
            device_id=device_id,
            station=station_id,
            protocol="mqtt",
            sample_rate_hz=256.0,
            channels=channels,
            role=DeviceRole.RECORDER,
            metadata={"source_file": rel_path},
        )
    )
    stations.append(
        StationConfig(
            station_id=station_id,
            profile="L18PLT",
            position_m=200.0 * (index - 1),
            channels=channels,
            dipole_length_m=50.0,
            ex_azimuth_deg=90.0,
            ey_azimuth_deg=0.0,
            device_ids=[device_id],
            metadata={
                "source_file": rel_path,
                "file_size_bytes": int(path.stat().st_size),
            },
        )
    )

session = FieldSession(
    "willy-amt-iot-demo",
    devices=devices,
    stations=stations,
    method="amt",
    monitoring_config=MonitoringConfig(
        method="amt",
        required_channels=["ex", "ey", "hx", "hy"],
        min_edge_acceptance_rate=0.85,
        min_battery_v=11.4,
        max_clock_offset_ms=5.0,
        frequency_band_hz=(1.0, 1000.0),
    ),
    metadata={
        "data_root": "data/AMT/WILLY_DATA/L18PLT",
        "telemetry_note": (
            "Operational telemetry is a deterministic demo overlay; "
            "station inventory is read from bundled real EDI files."
        ),
    },
)

print(session)
print(session.station_table().head())

# %%
# 2. Attach edge-QC, health, sync, and power telemetry
# ----------------------------------------------------
# File size is used only as a deterministic proxy for station-specific
# operational variation, so the example is reproducible on every machine.
# In a live deployment these packets would arrive through MQTT, serial,
# HTTP, or WebSocket clients in :mod:`pycsamt.iot.protocols`.

base_timestamp = 1_720_000_000.0

for index, (device, path, size) in enumerate(
    zip(devices, edi_files, sizes),
    start=1,
):
    rel_path = path.relative_to(repo_root()).as_posix()
    coverage = 0.86 + 0.12 * (size / size_scale)
    spike_fraction = 0.008 + 0.003 * ((index - 1) % 4)
    accepted = bool(coverage >= 0.90 and spike_fraction <= 0.02)
    rms_noise = 0.8 + 0.12 * index
    battery_v = 12.7 - 0.12 * index
    offset_ms = 0.35 + 0.42 * index
    drift_ppm = 0.7 + 0.9 * index
    jitter_ms = 0.08 + 0.03 * index

    qc_payload = {
        "station": device.station,
        "method": "amt",
        "channels": list(device.channels),
        "frequency_band_hz": [1.0, 1000.0],
        "accepted": accepted,
        "decision": "accept" if accepted else "review",
        "finite_coverage": round(float(min(coverage, 0.99)), 3),
        "spike_fraction": round(float(spike_fraction), 3),
        "rms": round(float(rms_noise), 3),
        "reasons": [] if accepted else ["coverage_review"],
        "source_file": rel_path,
    }
    session.add_packet(
        TelemetryPacket.from_device(
            device,
            timestamp=base_timestamp + 60.0 * index,
            payload=qc_payload,
            kind=PacketKind.QC,
            survey_id=session.survey_id,
        )
    )

    session.add_packet(
        TelemetryPacket.from_device(
            device,
            timestamp=base_timestamp + 60.0 * index + 5.0,
            payload={
                "station": device.station,
                "battery_v": round(float(battery_v), 2),
                "temperature_c": 25.0 + 0.4 * index,
                "latency_s": 1.2 + 0.1 * index,
                "packet_ok": True,
            },
            kind=PacketKind.HEALTH,
            survey_id=session.survey_id,
        )
    )

    session.add_packet(
        TelemetryPacket.from_device(
            device,
            timestamp=base_timestamp + 60.0 * index + 10.0,
            payload={
                "station": device.station,
                "offset_ms": round(float(offset_ms), 3),
                "within_tolerance": offset_ms <= 5.0,
                "reference": "gps",
                "drift_ppm": round(float(drift_ppm), 3),
                "jitter_ms": round(float(jitter_ms), 3),
                "gps_lock": index != len(devices),
                "n_reference_points": 24 - index,
                "quality": "good" if offset_ms <= 2.5 else "fair",
            },
            kind=PacketKind.SYNC,
            survey_id=session.survey_id,
        )
    )

    energy = estimate_energy_budget(
        EnergyConfig(
            battery_wh=120.0,
            active_power_w=1.6 + 0.05 * index,
            sleep_power_w=0.08,
            duty_cycle=0.68,
            solar_wh_per_day=18.0 - 0.6 * index,
            reserve_fraction=0.15,
            regulator_efficiency=0.9,
            telemetry_power_w=2.0,
            telemetry_seconds_per_day=900.0,
            edge_power_w=0.45,
            edge_duty_cycle=0.22,
            auxiliary_wh_per_day=0.8,
            min_runtime_days=3.0,
            device_id=device.device_id,
            metadata={"station": device.station, "source_file": rel_path},
        )
    )
    packet = energy.to_packet(
        device,
        timestamp=base_timestamp + 60.0 * index + 15.0,
        survey_id=session.survey_id,
    )
    packet.payload["station"] = device.station
    session.add_packet(packet)

status = session.assess(now=base_timestamp + 60.0 * len(devices) + 30.0)
print(status)

# %%
# 3. Visualise the IoT acquisition layer
# --------------------------------------
# ``plot_field_dashboard`` is the one-page field view.  The three
# companion plots expose the same packets as focused QC, power, and
# synchronisation summaries that can be saved into field reports.

plot_field_dashboard(
    session,
    now=base_timestamp + 60.0 * len(devices) + 30.0,
    station_axis="profile",
    title="WILLY L18 AMT real-data IoT field dashboard",
)

plot_edge_qc_summary(
    session,
    title="WILLY L18 AMT edge-QC telemetry",
)

plot_power_budget(
    session,
    title="WILLY L18 AMT recorder power budget",
)

plot_sync_quality(
    session,
    tolerance_ms=5.0,
    title="WILLY L18 AMT clock synchronisation",
)

plt.show()

