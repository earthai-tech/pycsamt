from __future__ import annotations

import math

import numpy as np
import pytest

from pycsamt.iot import (
    ClockSynchronizer,
    DeploymentConfig,
    DeviceConfig,
    DeviceRole,
    EdgeDecision,
    EdgeProcessingConfig,
    EdgeProcessor,
    EnergyConfig,
    DevicePowerProfile,
    EMMethod,
    IoTCapability,
    MonitoringConfig,
    MonitoringLevel,
    PacketKind,
    PowerState,
    TelemetryPacket,
    TelemetryMonitor,
    assess_telemetry,
    build_telemetry_client,
    deployment_report,
    edge_summary_table,
    estimate_clock_offset_ms,
    estimate_deployment_energy,
    estimate_energy_budget,
    monitoring_status_table,
    packet_table,
    power_summary_table,
    telemetry_summary,
)


def test_build_telemetry_client_dry_run_records_packet():
    client = build_telemetry_client("mqtt", endpoint="mqtt://broker")
    packet = TelemetryPacket(
        device_id="node-1",
        timestamp=10.0,
        topic="qc",
        payload={"snr": 12.5},
    )

    ack = client.send(packet)

    assert ack.ok is True
    assert ack.protocol == "mqtt"
    assert len(client.sent) == 1
    assert client.sent[0].kind == PacketKind.DATA


def test_http_client_requires_endpoint():
    from pycsamt.iot import TelemetryError

    client = build_telemetry_client("http", dry_run=False)
    packet = TelemetryPacket(
        device_id="node-1",
        timestamp=10.0,
        topic="qc",
        payload={},
    )

    with pytest.raises(TelemetryError):
        client.send(packet)


def test_unmapped_protocol_falls_back_to_recorder():
    # LoRa has no concrete transport; a real send must fail loudly.
    client = build_telemetry_client("lora", dry_run=False)
    packet = TelemetryPacket(
        device_id="node-1",
        timestamp=10.0,
        topic="qc",
        payload={},
    )

    with pytest.raises(NotImplementedError):
        client.send(packet)


def test_edge_processor_decimates_and_scores_coverage():
    proc = EdgeProcessor(
        EdgeProcessingConfig(decimation=2, finite_threshold=0.7)
    )
    result = proc.process(np.array([1.0, 2.0, np.nan, 4.0, 5.0]))

    assert result.data.shape[0] == 3
    assert result.metrics["original_samples"] == 5
    assert result.metrics["emitted_samples"] == 3
    assert result.metrics["finite_coverage"] == pytest.approx(2 / 3)
    assert result.accepted is False
    assert result.decision == EdgeDecision.REJECT


def test_edge_processor_channel_summaries_and_qc_packet():
    data = np.array(
        [
            [1.0, 10.0],
            [1.1, 10.1],
            [0.9, 9.9],
            [1.0, 60.0],
            [1.2, 10.2],
            [1.1, 10.1],
        ]
    )
    proc = EdgeProcessor(
        EdgeProcessingConfig(
            finite_threshold=0.9,
            channel_names=["EX", "HY"],
            spike_threshold=2.0,
            max_spike_fraction=0.1,
        )
    )
    device = DeviceConfig("node-1", station="S01", channels=["ex", "hy"])

    result = proc.process(data)
    packet = result.to_packet(
        device,
        timestamp=42.0,
        survey_id="survey-a",
        qos=1,
    )
    table = edge_summary_table(result)

    assert result.metrics["n_channel"] == 2
    assert result.channels[0].channel == "ex"
    assert result.channels[1].spike_fraction > 0
    assert result.accepted is False
    assert "spike_fraction_above_threshold" in result.reasons
    assert packet.kind == PacketKind.QC
    assert "survey-a/S01/node-1/qc" in packet.topic
    assert packet.payload["decision"] == "reject"
    assert table.loc[0, "channel"] == "ex"


def test_edge_processor_sample_axis_and_channel_override():
    data = np.array(
        [
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
        ]
    )
    proc = EdgeProcessor(
        EdgeProcessingConfig(decimation=2, sample_axis=1)
    )

    result = proc.process(data, channel_names=["ex", "hy"])

    assert result.data.shape == (2, 2)
    assert result.metrics["original_samples"] == 4
    assert result.metrics["emitted_samples"] == 2
    assert [ch.channel for ch in result.channels] == ["ex", "hy"]
    assert result.accepted is True


def test_edge_processor_rejects_invalid_config_and_channels():
    with pytest.raises(ValueError):
        EdgeProcessingConfig(decimation=0)
    with pytest.raises(ValueError):
        EdgeProcessingConfig(finite_threshold=1.5)
    with pytest.raises(ValueError):
        EdgeProcessingConfig(spike_threshold=0)

    proc = EdgeProcessor(EdgeProcessingConfig(channel_names=["ex", "hy"]))
    with pytest.raises(ValueError):
        proc.process(np.array([1.0, 2.0, 3.0]))


def test_clock_synchronizer_reports_offset_status():
    offset = estimate_clock_offset_ms([1.001, 2.001], [1.0, 2.0])
    status = ClockSynchronizer().assess(
        "node-1",
        [1.0005, 2.0005],
        [1.0, 2.0],
    )

    assert offset == pytest.approx(1.0)
    assert status.offset_ms == pytest.approx(0.5)
    assert status.within_tolerance is True


def test_energy_budget_runtime_with_solar_surplus_is_infinite():
    estimate = estimate_energy_budget(
        EnergyConfig(
            battery_wh=100.0,
            active_power_w=1.0,
            sleep_power_w=0.1,
            duty_cycle=0.1,
            solar_wh_per_day=20.0,
        )
    )

    assert estimate.average_power_w == pytest.approx(0.19)
    assert math.isinf(estimate.runtime_hours)
    assert estimate.state == PowerState.SUSTAINING


def test_energy_budget_accounts_for_reserve_and_iot_loads():
    estimate = estimate_energy_budget(
        EnergyConfig(
            battery_wh=120.0,
            active_power_w=2.0,
            sleep_power_w=0.2,
            duty_cycle=0.25,
            solar_wh_per_day=10.0,
            reserve_fraction=0.2,
            regulator_efficiency=0.8,
            charge_efficiency=0.9,
            telemetry_power_w=3.0,
            telemetry_seconds_per_day=600.0,
            edge_power_w=0.5,
            edge_duty_cycle=0.5,
            auxiliary_wh_per_day=1.0,
            min_runtime_days=8.0,
        )
    )
    table = power_summary_table(estimate, device_ids=["node-1"])
    packet = estimate.to_packet(
        DeviceConfig("node-1", station="S01"),
        timestamp=10.0,
        survey_id="survey-a",
    )

    assert estimate.usable_battery_wh == pytest.approx(96.0)
    assert estimate.telemetry_wh_per_day == pytest.approx(0.5)
    assert estimate.edge_wh_per_day == pytest.approx(6.0)
    assert estimate.net_wh_per_day > 0
    assert estimate.runtime_days < 8.0
    assert estimate.state in {PowerState.WARNING, PowerState.CRITICAL}
    assert "runtime_below_minimum" in estimate.issues
    assert table.loc[0, "device_id"] == "node-1"
    assert packet.kind == PacketKind.POWER
    assert packet.payload["state"] == estimate.state.value


def test_device_power_profile_and_deployment_energy_table():
    profile = DevicePowerProfile(
        "amt-node",
        active_power_w=1.2,
        sleep_power_w=0.1,
        telemetry_power_w=2.0,
        edge_power_w=0.2,
    )
    config = profile.apply(
        battery_wh=80.0,
        duty_cycle=0.5,
        telemetry_seconds_per_day=120.0,
        edge_duty_cycle=0.25,
        device_id="node-2",
    )

    table = estimate_deployment_energy([config])

    assert config.active_power_w == pytest.approx(1.2)
    assert table.loc[0, "device_id"] == "node-2"
    assert table.loc[0, "load_wh_per_day"] > 0


def test_energy_budget_rejects_invalid_inputs():
    with pytest.raises(ValueError):
        EnergyConfig(battery_wh=0.0, active_power_w=1.0)
    with pytest.raises(ValueError):
        EnergyConfig(
            battery_wh=1.0,
            active_power_w=1.0,
            duty_cycle=1.2,
        )
    with pytest.raises(ValueError):
        EnergyConfig(
            battery_wh=1.0,
            active_power_w=1.0,
            regulator_efficiency=0.0,
        )
    with pytest.raises(ValueError):
        DevicePowerProfile("", active_power_w=1.0)


def test_deployment_and_packet_reports_are_tables():
    deployment = DeploymentConfig(
        survey_id="survey-a",
        devices=[DeviceConfig("node-1", station="S01", channels=["ex"])],
        capabilities=[
            IoTCapability.TELEMETRY,
            IoTCapability.EDGE_PROCESSING,
        ],
    )
    packet = TelemetryPacket(
        device_id="node-1",
        timestamp=1.0,
        topic="health",
        payload={"battery_v": 12.2},
    )

    dep = deployment_report(deployment)
    packets = packet_table([packet])
    summary = telemetry_summary([packet])

    assert dep.loc[0, "survey_id"] == "survey-a"
    assert "telemetry" in dep.loc[0, "capabilities"]
    assert packets.loc[0, "payload_keys"] == "battery_v"
    assert summary.loc[0, "n_packet"] == 1


def test_telemetry_monitor_assesses_amt_stream_quality():
    packets = [
        TelemetryPacket(
            device_id="node-1",
            timestamp=100.0,
            topic="pycsamt/survey-a/S01/node-1/qc",
            kind=PacketKind.QC,
            payload={
                "method": "AMT",
                "station": "S01",
                "channels": ["EX", "HY"],
                "frequency_band_hz": [1.0, 1000.0],
                "accepted": True,
                "battery_v": 12.3,
                "clock_offset_ms": 2.0,
            },
        ),
        TelemetryPacket(
            device_id="node-1",
            timestamp=110.0,
            topic="pycsamt/survey-a/S01/node-1/qc",
            kind=PacketKind.QC,
            payload={
                "method": "amt",
                "station": "S01",
                "channels": ["ex", "hy"],
                "frequency_band_hz": [1.0, 1000.0],
                "accepted": False,
                "battery_v": 10.5,
                "clock_offset_ms": 8.0,
            },
        ),
    ]
    config = MonitoringConfig(
        method=EMMethod.AMT,
        expected_interval_s=10.0,
        min_edge_acceptance_rate=0.9,
        min_battery_v=11.0,
        max_clock_offset_ms=5.0,
        required_channels=["ex", "hy"],
        frequency_band_hz=(1.0, 1000.0),
    )

    status = TelemetryMonitor(config).assess(packets, now=115.0)
    table = monitoring_status_table(status)

    assert status.level == MonitoringLevel.CRITICAL
    assert status.edge_acceptance_rate == pytest.approx(0.5)
    assert "edge_acceptance_rate_below_threshold" in status.issues
    assert "battery_below_threshold" in status.issues
    assert "clock_offset_above_threshold" in status.issues
    assert status.methods == ["amt"]
    assert status.channels == ["ex", "hy"]
    assert table.loc[0, "level"] == "critical"


def test_telemetry_monitor_accepts_mapping_packets_and_enriched_table():
    packets = [
        {
            "device_id": "node-2",
            "timestamp": 1.0,
            "topic": "qc",
            "payload": {
                "survey_method": "csamt",
                "site": "S02",
                "channel": "hz",
                "frequency_hz": 128.0,
                "decision": "accept",
                "ack_ok": True,
                "latency_s": 2.0,
            },
        }
    ]
    monitor = TelemetryMonitor(
        MonitoringConfig(
            method="csamt",
            required_channels=["hz"],
            max_latency_s=5.0,
        )
    )

    status = assess_telemetry(packets, config=monitor.config)
    table = monitor.table(packets)

    assert status.ok is True
    assert status.frequency_min_hz == pytest.approx(128.0)
    assert table.loc[0, "method"] == "csamt"
    assert table.loc[0, "station"] == "S02"


def test_telemetry_monitor_rejects_invalid_config_and_no_data():
    with pytest.raises(ValueError):
        MonitoringConfig(min_packet_success_rate=1.2)
    with pytest.raises(ValueError):
        MonitoringConfig(frequency_band_hz=(100.0, 1.0))

    status = TelemetryMonitor().assess([])

    assert status.level == MonitoringLevel.NO_DATA
    assert status.n_packet == 0
    assert "no_packets" in status.issues


def test_core_objects_validate_clone_and_topics():
    device = DeviceConfig(
        "NODE-1",
        station="S01",
        protocol="MQTT",
        sample_rate_hz="128",
        channels=["EX", "HY"],
        role="gateway",
    )
    packet = TelemetryPacket.from_device(
        device,
        timestamp=12.0,
        payload={"battery_v": 12.5},
        kind="health",
        survey_id="survey-a",
        qos=1,
    )
    deployment = DeploymentConfig(
        survey_id="survey-a",
        devices=[device.as_dict()],
        capabilities=["telemetry", IoTCapability.HEALTH_MONITORING],
    )

    assert device.protocol == "mqtt"
    assert device.role == DeviceRole.GATEWAY
    assert device.channels == ["ex", "hy"]
    assert device.n_channels == 2
    assert "survey-a/S01/NODE-1/health" in packet.topic
    assert packet.kind == PacketKind.HEALTH
    assert deployment.n_devices == 1
    assert deployment.has_capability("telemetry")
    assert deployment.device("NODE-1").station == "S01"
    assert deployment.clone(notes="field test").notes == "field test"


def test_core_objects_reject_invalid_states():
    with pytest.raises(ValueError):
        DeviceConfig("", sample_rate_hz=1.0)
    with pytest.raises(ValueError):
        DeviceConfig("node", sample_rate_hz=0.0)
    with pytest.raises(ValueError):
        TelemetryPacket("node", 0.0, "topic", {}, qos=3)
    with pytest.raises(ValueError):
        DeploymentConfig(
            survey_id="survey",
            devices=[DeviceConfig("node"), DeviceConfig("node")],
        )
