"""Synthetic IoT/AMT field-network simulator.

These generators make the IoT subpackage easy to demo, document, and
test without hardware. They produce realistic-looking channel time
series (with configurable SNR, powerline contamination, and dropouts),
station/device configs, telemetry packets, GPS-drift clock pairs, and
battery-decay curves.

Randomness is controlled through a ``seed`` argument for reproducibility.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import (
    Any,
)

import numpy as np

from .core import DeviceConfig, PacketKind, TelemetryPacket
from .station import StationConfig

__all__ = [
    "simulate_powerline_noise",
    "simulate_amt_channel",
    "simulate_amt_station",
    "simulate_iot_network",
    "simulate_packet_loss",
    "simulate_gps_drift",
    "simulate_battery_decay",
]

_EPOCH_BASE = 1_700_000_000.0  # a fixed, plausible POSIX base time (2023-11)


def _rng(seed: int | None) -> np.random.Generator:
    return np.random.default_rng(seed)


def simulate_powerline_noise(
    n_samples: int,
    sample_rate: float,
    *,
    mains_hz: float = 50.0,
    amplitude: float = 1.0,
    n_harmonics: int = 3,
    seed: int | None = None,
) -> np.ndarray:
    """Return a powerline-noise time series (fundamental + harmonics)."""
    n_samples = int(n_samples)
    if n_samples <= 0:
        return np.empty(0)
    rng = _rng(seed)
    t = np.arange(n_samples) / float(sample_rate)
    noise = np.zeros(n_samples)
    for k in range(1, int(n_harmonics) + 1):
        amp = amplitude / k
        phase = rng.uniform(0.0, 2.0 * np.pi)
        noise += amp * np.sin(2.0 * np.pi * k * mains_hz * t + phase)
    return noise


def simulate_amt_channel(
    n_samples: int,
    sample_rate: float,
    *,
    snr_db: float = 20.0,
    mains_hz: float | None = 50.0,
    powerline_amplitude: float = 0.0,
    dropout_rate: float = 0.0,
    seed: int | None = None,
) -> np.ndarray:
    """Simulate one AMT channel: band-limited signal + noise + artefacts."""
    n_samples = int(n_samples)
    if n_samples <= 0:
        return np.empty(0)
    rng = _rng(seed)
    fs = float(sample_rate)
    t = np.arange(n_samples) / fs

    signal = np.zeros(n_samples)
    for _ in range(5):
        freq = rng.uniform(0.05, fs / 4.0)
        signal += rng.uniform(0.5, 1.5) * np.sin(
            2.0 * np.pi * freq * t + rng.uniform(0.0, 2.0 * np.pi)
        )
    sig_std = float(np.std(signal)) or 1.0
    signal /= sig_std

    snr_lin = 10.0 ** (float(snr_db) / 10.0)
    noise = rng.standard_normal(n_samples) / np.sqrt(max(snr_lin, 1e-9))
    x = signal + noise

    if mains_hz and powerline_amplitude > 0:
        x = x + simulate_powerline_noise(
            n_samples, fs, mains_hz=mains_hz,
            amplitude=powerline_amplitude, seed=int(rng.integers(0, 2**31)),
        )

    if dropout_rate and dropout_rate > 0:
        x = _inject_dropouts(x, dropout_rate, rng)
    return x


def _inject_dropouts(
    x: np.ndarray,
    dropout_rate: float,
    rng: np.random.Generator,
) -> np.ndarray:
    n = x.size
    n_drop = int(round(float(dropout_rate) * n))
    if n_drop <= 0:
        return x
    x = x.copy()
    # A few contiguous gaps rather than scattered single NaNs.
    n_gaps = max(1, n_drop // 16)
    for _ in range(n_gaps):
        length = max(1, n_drop // n_gaps)
        start = int(rng.integers(0, max(1, n - length)))
        x[start:start + length] = np.nan
    return x


def simulate_amt_station(
    station_id: str,
    *,
    channels: Sequence[str] = ("ex", "ey", "hx", "hy"),
    sample_rate: float = 128.0,
    n_samples: int = 1024,
    mains_hz: float = 50.0,
    snr_db: float = 20.0,
    powerline_amplitude: float = 0.2,
    dropout_rate: float = 0.0,
    survey_id: str | None = None,
    profile: str | None = None,
    position_m: float | None = None,
    lat: float | None = None,
    lon: float | None = None,
    timestamp: float | None = None,
    seed: int | None = None,
) -> dict[str, Any]:
    """Simulate a full AMT station: config, channel data, and packets.

    Returns
    -------
    dict
        Keys: ``station`` (:class:`StationConfig`), ``device``
        (:class:`DeviceConfig`), ``data`` (``{channel: ndarray}``), and
        ``packets`` (list of :class:`TelemetryPacket`: one health, one QC).
    """
    rng = _rng(seed)
    channels = [str(c).lower() for c in channels]
    device_id = f"node-{station_id}"
    ts = float(_EPOCH_BASE if timestamp is None else timestamp)

    device = DeviceConfig(
        device_id=device_id,
        station=station_id,
        protocol="mqtt",
        sample_rate_hz=sample_rate,
        channels=channels,
    )
    station = StationConfig(
        station_id=station_id,
        lat=lat,
        lon=lon,
        profile=profile,
        position_m=position_m,
        channels=channels,
        device_ids=[device_id],
    )

    data: dict[str, np.ndarray] = {}
    for channel in channels:
        data[channel] = simulate_amt_channel(
            n_samples, sample_rate,
            snr_db=snr_db, mains_hz=mains_hz,
            powerline_amplitude=powerline_amplitude,
            dropout_rate=dropout_rate,
            seed=int(rng.integers(0, 2**31)),
        )

    # Aggregate QC proxies from the synthetic data.
    coverage = float(np.mean([
        np.mean(np.isfinite(v)) for v in data.values()
    ])) if data else 1.0
    accepted = bool(coverage >= 0.95 and snr_db >= 6.0)
    battery_v = float(11.0 + rng.uniform(0.0, 2.0))

    health = TelemetryPacket.from_device(
        device,
        timestamp=ts,
        payload=dict(
            station=station_id,
            battery_v=round(battery_v, 3),
            temperature_c=round(float(20.0 + rng.uniform(-5, 10)), 2),
            firmware="sim-1.0",
        ),
        kind=PacketKind.HEALTH,
        survey_id=survey_id,
    )
    qc = TelemetryPacket.from_device(
        device,
        timestamp=ts + n_samples / float(sample_rate),
        payload=dict(
            method="amt",
            station=station_id,
            channels=channels,
            frequency_band_hz=[1.0, sample_rate / 2.0],
            finite_coverage=round(coverage, 4),
            accepted=accepted,
            decision="accept" if accepted else "reject",
        ),
        kind=PacketKind.QC,
        survey_id=survey_id,
    )

    return dict(
        station=station,
        device=device,
        data=data,
        packets=[health, qc],
    )


def simulate_iot_network(
    n_stations: int = 10,
    *,
    profiles: Sequence[str] = ("L1",),
    channels: Sequence[str] = ("ex", "ey", "hx", "hy"),
    sample_rate: float = 128.0,
    n_samples: int = 512,
    mains_hz: float = 50.0,
    snr_db: float = 20.0,
    dropout_rate: float = 0.05,
    survey_id: str = "SIM",
    station_spacing_m: float = 50.0,
    seed: int | None = None,
    detail: bool = False,
) -> Any:
    """Simulate a network of AMT stations across one or more profiles.

    Parameters
    ----------
    n_stations : int
        Total number of stations to generate.
    profiles : sequence of str
        Profile/line labels; stations are round-robin assigned to them.
    detail : bool
        When ``True`` return a dict with ``stations`` and ``packets``.
        When ``False`` (default) return just the flat list of packets, as
        in the documented example.

    Returns
    -------
    list of TelemetryPacket, or dict when ``detail=True``.
    """
    rng = _rng(seed)
    profiles = list(profiles) or ["L1"]
    n_stations = max(0, int(n_stations))
    per_profile_counts: dict[str, int] = {p: 0 for p in profiles}

    stations: list[dict[str, Any]] = []
    packets: list[TelemetryPacket] = []
    for i in range(n_stations):
        profile = profiles[i % len(profiles)]
        idx = per_profile_counts[profile]
        per_profile_counts[profile] = idx + 1
        station_id = f"{profile}-S{idx + 1:03d}"
        result = simulate_amt_station(
            station_id,
            channels=channels,
            sample_rate=sample_rate,
            n_samples=n_samples,
            mains_hz=mains_hz,
            snr_db=float(snr_db + rng.uniform(-4.0, 4.0)),
            dropout_rate=dropout_rate,
            survey_id=survey_id,
            profile=profile,
            position_m=idx * float(station_spacing_m),
            timestamp=_EPOCH_BASE + i * (n_samples / float(sample_rate) + 1.0),
            seed=int(rng.integers(0, 2**31)),
        )
        stations.append(result)
        packets.extend(result["packets"])

    if detail:
        return dict(stations=stations, packets=packets)
    return packets


def simulate_packet_loss(
    packets: Iterable[TelemetryPacket],
    dropout_rate: float = 0.05,
    *,
    seed: int | None = None,
) -> list[TelemetryPacket]:
    """Return *packets* with a fraction randomly dropped."""
    rng = _rng(seed)
    items = list(packets)
    if not 0.0 <= float(dropout_rate) <= 1.0:
        raise ValueError("dropout_rate must be between 0 and 1.")
    keep_mask = rng.random(len(items)) >= float(dropout_rate)
    return [pkt for pkt, keep in zip(items, keep_mask) if keep]


def simulate_gps_drift(
    n_samples: int,
    *,
    sample_interval_s: float = 1.0,
    drift_ppm: float = 5.0,
    jitter_ms: float = 0.2,
    offset_ms: float = 0.0,
    dropout_rate: float = 0.0,
    start_time: float = _EPOCH_BASE,
    seed: int | None = None,
) -> dict[str, np.ndarray]:
    """Simulate paired reference/local clocks with drift, jitter, dropout.

    Returns
    -------
    dict
        Keys ``reference`` and ``local`` (POSIX-second arrays) and
        ``gps_lock`` (bool array). Feed ``local``/``reference`` straight
        into :func:`~pycsamt.iot.sync.estimate_clock_drift_ppm` etc.
    """
    n_samples = int(n_samples)
    rng = _rng(seed)
    reference = start_time + np.arange(max(n_samples, 0)) * float(
        sample_interval_s
    )
    elapsed = reference - start_time
    drift_s = (float(drift_ppm) * 1e-6) * elapsed
    jitter_s = rng.normal(0.0, float(jitter_ms) / 1000.0, size=reference.size)
    local = reference + float(offset_ms) / 1000.0 + drift_s + jitter_s

    gps_lock = np.ones(reference.size, dtype=bool)
    if dropout_rate and dropout_rate > 0:
        gps_lock = rng.random(reference.size) >= float(dropout_rate)
        # Where GPS is lost, the local clock free-runs (no correction):
        # exaggerate drift on those samples for realism.
        local = np.where(gps_lock, local, local + drift_s * 0.5)
    return {"reference": reference, "local": local, "gps_lock": gps_lock}


def simulate_battery_decay(
    n_samples: int,
    *,
    initial_v: float = 13.2,
    final_v: float = 10.5,
    noise_v: float = 0.05,
    seed: int | None = None,
) -> np.ndarray:
    """Simulate a monotonic-ish battery discharge curve with noise."""
    n_samples = int(n_samples)
    if n_samples <= 0:
        return np.empty(0)
    rng = _rng(seed)
    frac = np.linspace(0.0, 1.0, n_samples)
    # Gentle exponential-ish sag from initial to final voltage.
    curve = final_v + (initial_v - final_v) * np.exp(-2.0 * frac)
    curve = curve + rng.normal(0.0, float(noise_v), size=n_samples)
    return curve
