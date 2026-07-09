from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pycsamt.iot import (
    DeviceConfig,
    EdgeProcessingConfig,
    EdgeProcessor,
    EnergyConfig,
    FieldSession,
    MonitoringConfig,
    PacketKind,
    StationConfig,
    SyncStatus,
    TelemetryPacket,
    plot_edge_qc_summary,
    plot_field_dashboard,
    plot_power_budget,
    plot_sync_quality,
)


def _dashboard_session() -> FieldSession:
    devices = [
        DeviceConfig("node-1", station="S01", channels=["ex", "hy"]),
        DeviceConfig("node-2", station="S02", channels=["ex", "hy"]),
    ]
    stations = [
        StationConfig(
            "S01",
            lat=6.5,
            lon=3.4,
            profile="L1",
            position_m=0.0,
            channels=["ex", "hy"],
        ),
        StationConfig(
            "S02",
            lat=6.6,
            lon=3.5,
            profile="L1",
            position_m=100.0,
            channels=["ex", "hy"],
        ),
    ]
    session = FieldSession(
        "survey-a",
        devices=devices,
        stations=stations,
        monitoring_config=MonitoringConfig(
            method="amt",
            required_channels=["ex", "hy"],
            min_edge_acceptance_rate=0.75,
            min_battery_v=11.0,
        ),
    )
    session.add_packet(
        TelemetryPacket(
            "node-1",
            1.0,
            "qc",
            {
                "station": "S01",
                "method": "amt",
                "channels": ["ex", "hy"],
                "accepted": True,
                "decision": "accept",
                "frequency_band_hz": [1.0, 1000.0],
            },
            kind=PacketKind.QC,
        )
    )
    session.add_packet(
        TelemetryPacket(
            "node-2",
            2.0,
            "qc",
            {
                "station": "S02",
                "method": "amt",
                "channels": ["ex", "hy"],
                "accepted": False,
                "decision": "reject",
                "reasons": ["spike_fraction_above_threshold"],
            },
            kind=PacketKind.QC,
        )
    )
    session.add_packet(
        TelemetryPacket(
            "node-1",
            3.0,
            "health",
            {"station": "S01", "battery_v": 12.4},
            kind=PacketKind.HEALTH,
        )
    )
    session.add_packet(
        TelemetryPacket(
            "node-2",
            4.0,
            "health",
            {"station": "S02", "battery_v": 10.8},
            kind=PacketKind.HEALTH,
        )
    )
    session.add_packet(
        TelemetryPacket(
            "node-2",
            5.0,
            "sync",
            {"station": "S02", "quality": "poor", "clock_offset_ms": 9.0},
            kind=PacketKind.SYNC,
        )
    )
    return session


def test_plot_field_dashboard_returns_figure_and_attaches_data():
    session = _dashboard_session()

    fig = plot_field_dashboard(session)

    assert isinstance(fig, plt.Figure)
    assert len(fig.axes) >= 4
    data = fig.pycsamt_iot_dashboard
    assert data["survey_id"] == "survey-a"
    assert data["n_packets"] == 5
    assert len(data["stations"]) == 2
    levels = {row["station_id"]: row["health_level"] for row in data["stations"]}
    assert levels["S01"] == "ok"
    assert levels["S02"] == "critical"
    plt.close(fig)


def test_plot_field_dashboard_accepts_mapping_and_map_axis():
    session = _dashboard_session()

    fig = plot_field_dashboard(session.to_dict(), station_axis="map")

    assert isinstance(fig, plt.Figure)
    assert fig.pycsamt_iot_dashboard["n_stations"] == 2
    plt.close(fig)


def test_plot_field_dashboard_handles_empty_session():
    session = FieldSession("empty")

    fig = plot_field_dashboard(session)

    assert isinstance(fig, plt.Figure)
    assert fig.pycsamt_iot_dashboard["n_packets"] == 0
    assert fig.pycsamt_iot_dashboard["issues"] == ["no_packets"]
    plt.close(fig)


def test_plot_edge_qc_summary_accepts_edge_results():
    proc = EdgeProcessor(
        EdgeProcessingConfig(
            channel_names=["ex", "hy"],
            finite_threshold=0.8,
            warn_finite_threshold=0.95,
            spike_threshold=2.0,
            max_spike_fraction=0.4,
        )
    )
    result = proc.process(
        [[1.0, 2.0], [1.1, 2.1], [float("nan"), 50.0], [1.0, 2.0]]
    )

    fig = plot_edge_qc_summary(result)

    assert isinstance(fig, plt.Figure)
    rows = fig.pycsamt_iot_edge_qc
    assert len(rows) == 2
    assert {row["channel"] for row in rows} == {"ex", "hy"}
    assert all("decision" in row for row in rows)
    plt.close(fig)


def test_plot_edge_qc_summary_accepts_session_packets():
    session = _dashboard_session()

    fig = plot_edge_qc_summary(session)

    assert isinstance(fig, plt.Figure)
    rows = fig.pycsamt_iot_edge_qc
    assert len(rows) == 4
    assert {row["station"] for row in rows} == {"S01", "S02"}
    assert {row["decision"] for row in rows} == {"accept", "reject"}
    assert {row["channel"] for row in rows} == {"ex", "hy"}
    plt.close(fig)


def test_plot_power_budget_accepts_configs_and_preserves_rows():
    configs = [
        EnergyConfig(
            battery_wh=100.0,
            active_power_w=1.0,
            sleep_power_w=0.1,
            duty_cycle=0.2,
            solar_wh_per_day=5.0,
            telemetry_power_w=2.0,
            telemetry_seconds_per_day=300.0,
            edge_power_w=0.2,
            edge_duty_cycle=0.25,
            device_id="node-1",
        ),
        EnergyConfig(
            battery_wh=20.0,
            active_power_w=2.0,
            sleep_power_w=0.2,
            duty_cycle=0.8,
            min_runtime_days=10.0,
            device_id="node-2",
        ),
    ]

    fig = plot_power_budget(configs)

    assert isinstance(fig, plt.Figure)
    rows = fig.pycsamt_iot_power_budget
    assert len(rows) == 2
    assert {row["device_id"] for row in rows} == {"node-1", "node-2"}
    assert all(row["load_wh_per_day"] > 0 for row in rows)
    assert rows[1]["state"] in {"warning", "critical"}
    plt.close(fig)


def test_plot_power_budget_accepts_session_power_packets():
    device = DeviceConfig("node-1", station="S01")
    estimate = EnergyConfig(
        battery_wh=100.0,
        active_power_w=1.0,
        duty_cycle=0.1,
        solar_wh_per_day=50.0,
        device_id="node-1",
    )
    from pycsamt.iot import estimate_energy_budget

    packet = estimate_energy_budget(estimate).to_packet(
        device,
        timestamp=10.0,
        survey_id="survey-a",
    )
    session = FieldSession("survey-a", devices=[device])
    session.add_packet(packet)

    fig = plot_power_budget(session)

    rows = fig.pycsamt_iot_power_budget
    assert len(rows) == 1
    assert rows[0]["device_id"] == "node-1"
    assert rows[0]["state"] == "sustaining"
    assert rows[0]["runtime_days"] == float("inf")
    plt.close(fig)


def test_plot_sync_quality_accepts_status_objects():
    statuses = [
        SyncStatus(
            "node-1",
            offset_ms=0.2,
            within_tolerance=True,
            drift_ppm=1.0,
            jitter_ms=0.1,
            gps_lock=True,
            n_reference_points=10,
            quality="excellent",
        ),
        SyncStatus(
            "node-2",
            offset_ms=8.0,
            within_tolerance=False,
            drift_ppm=80.0,
            jitter_ms=3.0,
            gps_lock=False,
            n_reference_points=8,
            quality="poor",
        ),
    ]

    fig = plot_sync_quality(statuses, tolerance_ms=1.0, max_drift_ppm=10.0)

    assert isinstance(fig, plt.Figure)
    rows = fig.pycsamt_iot_sync_quality
    assert len(rows) == 2
    assert {row["device_id"] for row in rows} == {"node-1", "node-2"}
    assert {row["quality"] for row in rows} == {"excellent", "poor"}
    assert rows[0]["gps_lock"] is True
    plt.close(fig)


def test_plot_sync_quality_accepts_session_packets():
    session = FieldSession("survey-a")
    session.add_packet(
        TelemetryPacket(
            "node-1",
            1.0,
            "sync",
            {
                "offset_ms": 0.5,
                "drift_ppm": 2.0,
                "jitter_ms": 0.2,
                "gps_lock": True,
                "n_reference_points": 12,
                "quality": "good",
            },
            kind=PacketKind.SYNC,
        )
    )

    fig = plot_sync_quality(session)

    rows = fig.pycsamt_iot_sync_quality
    assert len(rows) == 1
    assert rows[0]["device_id"] == "node-1"
    assert rows[0]["quality"] == "good"
    assert rows[0]["offset_ms"] == 0.5
    plt.close(fig)
