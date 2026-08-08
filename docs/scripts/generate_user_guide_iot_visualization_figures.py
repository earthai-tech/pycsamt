"""Generate the executed figures for user_guide/iot/visualization.rst."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.iot import (  # noqa: E402
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
    plot_edge_qc_summary,
    plot_field_dashboard,
    plot_power_budget,
    plot_sync_quality,
)

IMAGE_DIR = ROOT / "docs/source/images/user_guide/iot"

STATION_IDS = ["L18-001", "L18-002", "L18-003"]
POSITIONS_M = [0.0, 50.0, 100.0]
COORDS = [(7.501, -5.201), (7.502, -5.198), (7.503, -5.195)]
BATTERY_V = [12.55, 11.84, 10.92]
TEMPERATURE_C = [30.1, 31.4, 33.2]
BATTERY_WH = [180.0, 120.0, 70.0]
ACTIVE_POWER_W = [1.2, 1.6, 2.1]
DUTY_CYCLE = [0.25, 0.35, 0.55]
SOLAR_WH_PER_DAY = [8.0, 3.0, 0.0]
SYNC_SPECS = [
    {"offset_ms": 0.32, "drift_ppm": 1.4, "jitter_ms": 0.16,
     "gps_lock": True, "n_reference_points": 120, "quality": "excellent"},
    {"offset_ms": 1.84, "drift_ppm": 11.2, "jitter_ms": 0.72,
     "gps_lock": True, "n_reference_points": 118, "quality": "fair"},
    {"offset_ms": 7.62, "drift_ppm": 69.5, "jitter_ms": 2.15,
     "gps_lock": False, "n_reference_points": 61, "quality": "poor"},
]


def _build_session() -> FieldSession:
    """Rebuild the deterministic 3-station WILLY-L18-VISUAL session."""
    devices = [
        DeviceConfig(
            f"l18-node-{i + 1:02d}", station=sid,
            channels=["ex", "ey", "hx", "hy"], sample_rate_hz=256.0,
        )
        for i, sid in enumerate(STATION_IDS)
    ]
    stations = [
        StationConfig(
            sid, profile="L18", position_m=pos, lat=lat, lon=lon,
            channels=["ex", "ey", "hx", "hy"],
        )
        for sid, pos, (lat, lon) in zip(STATION_IDS, POSITIONS_M, COORDS)
    ]
    session = FieldSession(
        "WILLY-L18-VISUAL", devices=devices, stations=stations,
        monitoring_config=MonitoringConfig(
            method="amt", expected_interval_s=6.0, max_gap_s=18.0,
            min_edge_acceptance_rate=0.70, min_battery_v=11.2,
            required_channels=["ex", "ey", "hx", "hy"],
        ),
    )
    processor = EdgeProcessor(
        EdgeProcessingConfig(
            decimation=1, finite_threshold=0.85, warn_finite_threshold=0.95,
            channel_names=["ex", "ey", "hx", "hy"], spike_threshold=4.0,
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
            device, timestamp=1_700_000_000.0 + idx * 6.0,
            survey_id=session.survey_id, qos=1,
        )
        qc_packet.payload["station"] = device.station
        qc_packet.payload["method"] = "amt"
        session.add_packet(qc_packet)

        health = TelemetryPacket.from_device(
            device, timestamp=1_700_000_002.0 + idx * 6.0,
            survey_id=session.survey_id, kind=PacketKind.HEALTH,
            payload={
                "station": device.station, "battery_v": BATTERY_V[idx],
                "temperature_c": TEMPERATURE_C[idx],
            },
        )
        session.add_packet(health)

        power_config = EnergyConfig(
            battery_wh=BATTERY_WH[idx], active_power_w=ACTIVE_POWER_W[idx],
            sleep_power_w=0.18, duty_cycle=DUTY_CYCLE[idx],
            solar_wh_per_day=SOLAR_WH_PER_DAY[idx], telemetry_power_w=2.0,
            telemetry_seconds_per_day=600.0, edge_power_w=0.35,
            edge_duty_cycle=0.20, min_runtime_days=5.0,
            device_id=device.device_id,
        )
        power_packet = estimate_energy_budget(power_config).to_packet(
            device, timestamp=1_700_000_004.0 + idx * 6.0,
            survey_id=session.survey_id,
        )
        power_packet.payload["station"] = device.station
        session.add_packet(power_packet)

        sync_payload = dict(SYNC_SPECS[idx], station=device.station)
        session.add_packet(
            TelemetryPacket.from_device(
                device, timestamp=1_700_000_006.0 + idx * 6.0,
                survey_id=session.survey_id, kind=PacketKind.SYNC,
                payload=sync_payload,
            )
        )
    return session


def make_field_dashboard_profile() -> None:
    session = _build_session()
    fig = plot_field_dashboard(
        session, now=1_700_000_030.0, station_axis="profile",
        title="IoT field dashboard: WILLY L18",
    )
    fig.savefig(IMAGE_DIR / "user-guide-iot-visualization-01.png", dpi=180)


def make_field_dashboard_map() -> None:
    session = _build_session()
    fig = plot_field_dashboard(
        session, now=1_700_000_030.0, station_axis="map",
        title="IoT field dashboard (map view): WILLY L18",
    )
    fig.savefig(IMAGE_DIR / "user-guide-iot-visualization-05.png", dpi=180)


def make_edge_qc_summary() -> None:
    session = _build_session()
    fig = plot_edge_qc_summary(session, title="Edge QC visual summary: WILLY L18")
    fig.savefig(IMAGE_DIR / "user-guide-iot-visualization-02.png", dpi=180)


def make_power_budget() -> None:
    session = _build_session()
    fig = plot_power_budget(session, title="Power budget visual summary: WILLY L18")
    fig.savefig(IMAGE_DIR / "user-guide-iot-visualization-03.png", dpi=180)


def make_sync_quality() -> None:
    session = _build_session()
    fig = plot_sync_quality(
        session, title="Clock sync visual summary: WILLY L18",
        tolerance_ms=1.0, max_drift_ppm=10.0, max_jitter_ms=1.0,
    )
    fig.savefig(IMAGE_DIR / "user-guide-iot-visualization-04.png", dpi=180)


if __name__ == "__main__":
    IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    make_field_dashboard_profile()
    make_field_dashboard_map()
    make_edge_qc_summary()
    make_power_budget()
    make_sync_quality()
