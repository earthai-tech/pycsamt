"""Tests for the extended IoT subpackage: station, session, schemas,
AMT edge QC, sync audit, real transports, security, provenance, and the
field-network simulator.
"""

from __future__ import annotations

import json
import math

import numpy as np
import pytest

from pycsamt.iot import (
    AcquisitionManifest,
    AuthScheme,
    Credential,
    EdgeDecision,
    EdgeProcessingConfig,
    EdgeProcessor,
    FieldSession,
    FileTelemetryClient,
    HealthPayload,
    QCPayload,
    SecurityConfig,
    StationConfig,
    SyncQuality,
    TelemetryPacket,
    TLSConfig,
    assess_impedance_stability,
    assess_sync_quality,
    batch_assess_sync,
    build_acquisition_manifest,
    build_telemetry_client,
    check_channel_saturation,
    check_contact_resistance,
    compute_live_spectra,
    detect_gps_dropout,
    detect_powerline_harmonics,
    detect_sensor_dropout,
    estimate_channel_snr,
    estimate_clock_drift_ppm,
    estimate_frequency_coverage,
    export_reproducibility_bundle,
    hash_bytes,
    hash_mapping,
    hash_raw_file,
    log_qc_decision,
    parse_payload,
    redact_secret,
    simulate_amt_station,
    simulate_battery_decay,
    simulate_gps_drift,
    simulate_iot_network,
    simulate_packet_loss,
    station_table,
    validate_payload,
)


# ---------------------------------------------------------------------------
# core hardening
# ---------------------------------------------------------------------------
def test_packet_rejects_non_finite_and_negative_timestamp():
    with pytest.raises(ValueError):
        TelemetryPacket("node", float("nan"), "t", {})
    with pytest.raises(ValueError):
        TelemetryPacket("node", -1.0, "t", {})
    with pytest.raises(ValueError):
        TelemetryPacket("node", float("inf"), "t", {})
    # Zero (relative session clock) is allowed.
    assert TelemetryPacket("node", 0.0, "t", {}).timestamp == 0.0


def test_string_booleans_are_parsed_explicitly():
    from pycsamt.iot import DeviceConfig

    device = DeviceConfig("node", enabled="false")
    assert device.enabled is False
    packet = TelemetryPacket("node", 1.0, "t", {}, retained="0")
    assert packet.retained is False


# ---------------------------------------------------------------------------
# station
# ---------------------------------------------------------------------------
def test_station_config_validates_and_tabulates():
    station = StationConfig(
        "S01",
        lat=6.5,
        lon=3.4,
        elevation=120.0,
        profile="L1",
        position_m=50.0,
        channels=["EX", "HY"],
        device_ids=["node-1"],
    )
    assert station.has_location
    assert station.coords == (6.5, 3.4, 120.0)
    assert station.channels == ["ex", "hy"]

    with pytest.raises(ValueError):
        StationConfig("S02", lat=200.0)

    table = station_table([station])
    assert table.loc[0, "station_id"] == "S01"
    assert table.loc[0, "n_channels"] == 2


# ---------------------------------------------------------------------------
# schemas
# ---------------------------------------------------------------------------
def test_health_payload_canonicalises_aliases():
    payload = parse_payload(
        "health",
        {"battery_voltage": 12.3, "temp": 21.0, "site": "S01", "foo": "bar"},
    )
    assert isinstance(payload, HealthPayload)
    assert payload.battery_v == 12.3
    assert payload.temperature_c == 21.0
    assert payload.station == "S01"
    assert payload.extra == {"foo": "bar"}
    flat = payload.as_dict(drop_none=True)
    assert flat["battery_v"] == 12.3 and flat["foo"] == "bar"
    assert "uptime_s" not in flat  # dropped because None


def test_qc_payload_parses_string_bool_and_band():
    payload = validate_payload(
        "qc", {"ok": "false", "channel": "EX", "band_hz": [1.0, 1000.0]}
    )
    assert payload["accepted"] is False
    assert payload["channels"] == ["ex"]
    assert payload["frequency_band_hz"] == [1.0, 1000.0]

    assert isinstance(parse_payload("qc", {"snr": 15.0}), QCPayload)
    with pytest.raises(ValueError):
        validate_payload("qc", {"band_hz": [1000.0, 1.0]})  # unordered


# ---------------------------------------------------------------------------
# edge (warning state + robust bool)
# ---------------------------------------------------------------------------
def test_edge_processor_emits_warning_state():
    proc = EdgeProcessor(
        EdgeProcessingConfig(
            finite_threshold=0.5,
            warn_finite_threshold=0.95,
            compute_spikes=False,
        )
    )
    data = np.array([1.0, 2.0, np.nan, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, np.nan])
    result = proc.process(data)
    assert result.accepted is True
    assert result.decision == EdgeDecision.WARNING
    assert "finite_coverage_marginal" in result.metrics["warnings"]


# ---------------------------------------------------------------------------
# edge_amt
# ---------------------------------------------------------------------------
def test_detect_powerline_harmonics_flags_mains_tone():
    fs = 500.0
    n = 4096
    t = np.arange(n) / fs
    rng = np.random.default_rng(0)
    signal = np.sin(2 * np.pi * 50.0 * t) + 0.01 * rng.standard_normal(n)
    result = detect_powerline_harmonics(signal, fs, mains_hz=50.0)
    assert result.contaminated is True
    assert result.dominant is not None
    assert result.dominant.frequency_hz == pytest.approx(50.0)

    clean = rng.standard_normal(n)
    assert (
        detect_powerline_harmonics(clean, fs, mains_hz=50.0).contaminated
        is False
    )


def test_channel_snr_orders_clean_above_noisy():
    fs = 100.0
    t = np.arange(1000) / fs
    clean = np.sin(2 * np.pi * 1.0 * t)
    rng = np.random.default_rng(1)
    noisy = clean + 2.0 * rng.standard_normal(clean.size)
    snr_clean = estimate_channel_snr(clean)
    snr_noisy = estimate_channel_snr(noisy)
    assert snr_clean > snr_noisy
    assert snr_clean > 10.0


def test_saturation_and_contact_resistance():
    t = np.linspace(0, 1, 1000)
    clean = np.sin(2 * np.pi * 3 * t)
    clipped = np.clip(3.0 * clean, -1.0, 1.0)
    assert check_channel_saturation(clean)["saturated"] is False
    assert check_channel_saturation(clipped, limit=1.0)["saturated"] is True

    rng = np.random.default_rng(2)
    ex = 10.0 * rng.standard_normal(500)
    ey = 0.1 * rng.standard_normal(500)
    check = check_contact_resistance(ex, ey, noise_rms_threshold=1.0)
    assert check["ok"] is False
    assert "ex_ey_noise_imbalance" in check["flags"]


def test_frequency_coverage_and_live_spectra():
    fs = 500.0
    n = 4096
    t = np.arange(n) / fs
    tone = np.sin(2 * np.pi * 50.0 * t)
    coverage = estimate_frequency_coverage(
        tone, fs, target_bands=[(2.0, 100.0), (150.0, 200.0)]
    )
    assert coverage.nyquist_hz == pytest.approx(250.0)
    assert 0.0 < coverage.f_low_hz <= coverage.f_high_hz <= 250.0
    assert math.isfinite(coverage.n_decades)
    # Only the low band falls inside the resolved range of the tone.
    assert coverage.coverage_fraction == pytest.approx(0.5)
    assert (150.0, 200.0) in coverage.missing_bands

    spectra = compute_live_spectra(tone, fs)
    assert spectra["frequency_hz"].shape == spectra["psd"].shape
    assert spectra["frequency_hz"][0] == pytest.approx(0.0)


def test_impedance_stability_and_sensor_dropout():
    stable = np.full(12, 1.0 + 0.5j) + 1e-3 * (
        np.random.default_rng(3).standard_normal(12)
    )
    assert assess_impedance_stability(stable).stable is True
    unstable = np.random.default_rng(4).standard_normal(12) + 1j * (
        np.random.default_rng(5).standard_normal(12)
    )
    assert assess_impedance_stability(unstable).stable is False

    series = np.concatenate(
        [
            np.random.default_rng(6).standard_normal(50),
            np.full(30, 2.0),  # stuck sensor
        ]
    )
    series[5] = np.nan
    drop = detect_sensor_dropout(series)
    assert drop["dropout"] is True
    assert drop["longest_flat_run"] >= 30
    assert drop["n_nan"] == 1


# ---------------------------------------------------------------------------
# sync audit
# ---------------------------------------------------------------------------
def test_sync_drift_quality_and_batch():
    sim = simulate_gps_drift(
        200, sample_interval_s=1.0, drift_ppm=5.0, jitter_ms=0.02, seed=7
    )
    drift = estimate_clock_drift_ppm(sim["local"], sim["reference"])
    assert drift == pytest.approx(5.0, abs=1.0)

    quality = assess_sync_quality(offset_ms=0.2, drift_ppm=1.0, jitter_ms=0.1)
    assert quality == SyncQuality.EXCELLENT
    assert assess_sync_quality(offset_ms=float("nan")) == SyncQuality.UNKNOWN

    table = batch_assess_sync(
        {
            "node-1": (sim["local"], sim["reference"]),
            "node-2": {
                "local": sim["local"],
                "reference": sim["reference"],
                "gps_lock": True,
            },
        }
    )
    assert set(table["device_id"]) == {"node-1", "node-2"}


def test_detect_gps_dropout_intervals():
    flags = [True] * 10 + [False] * 5 + [True] * 5
    result = detect_gps_dropout(flags, timestamps=np.arange(20.0))
    assert result["n_dropout_events"] == 1
    assert result["longest_dropout_samples"] == 5
    assert result["lock_fraction"] == pytest.approx(0.75)
    assert result["ok"] is False


# ---------------------------------------------------------------------------
# protocols
# ---------------------------------------------------------------------------
def test_file_client_round_trip(tmp_path):
    path = tmp_path / "telemetry.jsonl"
    packets = [
        TelemetryPacket("node-1", 1.0, "t", {"a": 1}),
        TelemetryPacket("node-1", 2.0, "t", {"b": 2}),
    ]
    with FileTelemetryClient(str(path), dry_run=False) as client:
        for pkt in packets:
            ack = client.send(pkt)
            assert ack.ok is True
        assert client.healthcheck() is True

    reader = FileTelemetryClient(str(path), dry_run=False)
    rows = reader.read_all()
    assert len(rows) == 2
    assert rows[0]["payload"]["a"] == 1
    reader.connect()
    assert reader.receive()["timestamp"] == 1.0
    assert reader.receive()["timestamp"] == 2.0
    assert reader.receive() is None


def test_dry_run_client_records_without_io(tmp_path):
    path = tmp_path / "should_not_exist.jsonl"
    client = build_telemetry_client(
        "file", endpoint=str(path)
    )  # dry_run=True
    client.send(TelemetryPacket("node-1", 1.0, "t", {}))
    assert len(client.sent) == 1
    assert not path.exists()


# ---------------------------------------------------------------------------
# security
# ---------------------------------------------------------------------------
def test_credentials_redacted_and_headers():
    cred = Credential(scheme="bearer", token="super-secret-token")
    assert cred.headers()["Authorization"] == "Bearer super-secret-token"
    assert cred.as_dict()["token"] == "***REDACTED***"
    assert "super-secret-token" not in repr(cred)
    assert redact_secret("x") == "***REDACTED***"
    assert redact_secret(None) is None

    with pytest.raises(ValueError):
        Credential(scheme="bearer")  # missing token


def test_security_config_client_options_and_env(monkeypatch):
    config = SecurityConfig(
        tls=TLSConfig(enabled=True),
        credential=Credential(scheme=AuthScheme.BEARER, token="tok"),
        require_tls=True,
    )
    options = config.client_options()
    assert options["tls"] is True
    assert options["token"] == "tok"
    assert options["headers"]["Authorization"] == "Bearer tok"
    assert "***REDACTED***" in json.dumps(config.as_dict())

    monkeypatch.setenv("PYCSAMT_IOT_TOKEN", "envtok")
    from_env = SecurityConfig.from_env()
    assert from_env.credential.scheme == AuthScheme.BEARER
    assert from_env.credential.token == "envtok"


# ---------------------------------------------------------------------------
# provenance
# ---------------------------------------------------------------------------
def test_hashing_helpers(tmp_path):
    import hashlib

    assert hash_bytes(b"abc") == hashlib.sha256(b"abc").hexdigest()
    assert hash_mapping({"a": 1, "b": 2}) == hash_mapping({"b": 2, "a": 1})

    raw = tmp_path / "raw.bin"
    raw.write_bytes(b"0123456789")
    record = hash_raw_file(str(raw))
    assert record["bytes"] == 10
    assert record["digest"] == hashlib.sha256(b"0123456789").hexdigest()


def test_manifest_build_write_and_bundle(tmp_path):
    decision = log_qc_decision(
        "S01", "Accept", channel="EX", reasons=["snr_ok"]
    )
    assert decision["decision"] == "accept"
    manifest = build_acquisition_manifest(
        "SURVEY-1",
        records=[{"station_id": "S01", "lat": 6.5, "lon": 3.4}],
        qc_decisions=[decision],
        method="amt",
    )
    assert isinstance(manifest, AcquisitionManifest)
    path = tmp_path / "manifest.json"
    manifest.write(str(path))
    loaded = json.loads(path.read_text())
    assert loaded["survey_id"] == "SURVEY-1"
    assert "content_hash" in loaded

    bundle = export_reproducibility_bundle(manifest, str(tmp_path / "bundle"))
    assert bundle["manifest"].endswith("acquisition_manifest.json")
    assert len(bundle["audits"]) == 1


# ---------------------------------------------------------------------------
# simulator
# ---------------------------------------------------------------------------
def test_simulator_network_and_station():
    packets = simulate_iot_network(
        n_stations=6,
        profiles=["L1", "L3"],
        channels=["ex", "ey", "hx", "hy"],
        dropout_rate=0.05,
        seed=0,
    )
    assert len(packets) == 12  # one health + one qc per station
    assert all(isinstance(p, TelemetryPacket) for p in packets)

    station = simulate_amt_station("S01", n_samples=256, seed=1)
    assert set(station["data"].keys()) == {"ex", "ey", "hx", "hy"}
    assert len(station["packets"]) == 2

    fewer = simulate_packet_loss(packets, dropout_rate=0.5, seed=2)
    assert len(fewer) <= len(packets)

    battery = simulate_battery_decay(100, seed=3)
    assert battery[0] > battery[-1]


# ---------------------------------------------------------------------------
# field session (end to end)
# ---------------------------------------------------------------------------
def _qc_packet(station="S01", accepted=True, ts=100.0):
    return TelemetryPacket(
        device_id="node-1",
        timestamp=ts,
        topic=f"pycsamt/survey/{station}/node-1/qc",
        kind="qc",
        payload={
            "method": "amt",
            "station": station,
            "channels": ["ex", "hy"],
            "frequency_band_hz": [1.0, 1000.0],
            "accepted": accepted,
        },
    )


def test_field_session_end_to_end(tmp_path):
    from pycsamt.iot import DeviceConfig

    session = FieldSession(
        "SSL2026",
        devices=[
            DeviceConfig("node-1", station="S01", channels=["ex", "hy"])
        ],
        method="amt",
    )
    session.add_station(StationConfig("S01", lat=6.5, lon=3.4, profile="L1"))
    session.add_packet(_qc_packet(accepted=True, ts=100.0))
    session.add_packet(_qc_packet(accepted=False, ts=110.0))

    status = session.assess()
    assert status.n_packet == 2

    sites = session.to_sites()
    assert [s.station_id for s in sites] == ["S01"]
    assert "ex" in sites[0].channels

    pipeline = session.to_pipeline_input()
    station_out = pipeline["stations"][0]
    assert station_out["station_id"] == "S01"
    assert station_out["coords"] == (6.5, 3.4, 0.0)
    assert station_out["accepted_band_hz"] == [1.0, 1000.0]
    assert station_out["acceptance_rate"] == pytest.approx(0.5)

    # JSON round-trip preserves packets and stations.
    restored = FieldSession.from_json(session.to_json())
    assert restored.n_packets == 2
    assert restored.n_stations == 1

    manifest_path = tmp_path / "field_manifest.json"
    session.export_manifest(str(manifest_path))
    manifest = json.loads(manifest_path.read_text())
    assert manifest["survey_id"] == "SSL2026"
    assert manifest["records"][0]["station_id"] == "S01"
    assert len(manifest["records"][0]["qc_decisions"]) == 2
