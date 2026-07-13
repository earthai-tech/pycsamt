r"""
Telemetry transport, schemas, and security
==========================================

A complete tour of the IoT *infrastructure* layer — the part that moves QC
and health data off the recorder and into a processing host. It shows, on
a simulated but realistic telemetry stream:

1. **schemas** — folding messy, vendor-specific payload keys into a
   canonical form (:func:`~pycsamt.iot.validate_payload`);
2. **a file transport round-trip** — writing newline-delimited JSON with
   :class:`~pycsamt.iot.FileTelemetryClient` and reading it back;
3. **the unified client interface** — the same ``connect/send/receive/
   healthcheck`` API across MQTT, HTTP, serial, and WebSocket, all usable
   offline in ``dry_run`` mode;
4. **security** — TLS plus a bearer credential, with secrets redacted from
   logs but still delivered to the live client; and
5. **replay** — reading the log back into a :class:`~pycsamt.iot.FieldSession`
   and re-assessing it, proving the round-trip is loss-less.

Every transport works here without a broker, server, or serial device.
"""

# sphinx_gallery_thumbnail_number = 2

# %%
# 1. Canonical telemetry schemas
# ------------------------------
# Field fleets are heterogeneous: one node reports ``battery_voltage``,
# another ``voltage``. The schema layer folds those aliases into one
# canonical key and range-checks values, while preserving unknown fields.

from __future__ import annotations

import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from pycsamt.iot import (
    AuthScheme,
    Credential,
    FieldSession,
    FileTelemetryClient,
    PacketKind,
    SecurityConfig,
    TelemetryPacket,
    TLSConfig,
    build_telemetry_client,
    simulate_iot_network,
    validate_payload,
)

messy = {
    "battery_voltage": 12.42,
    "temp": 21.7,
    "site": "L1-S007",
    "rssi": -78,
    "fw": "2.3.1",
    "vendor": "acme-mk2",
}
clean = validate_payload("health", messy)
print("raw keys:   ", sorted(messy))
print("clean keys: ", sorted(clean))
print(
    f"battery_voltage -> battery_v = {clean['battery_v']}  "
    f"| temp -> temperature_c = {clean['temperature_c']}"
)
print(f"unknown field kept verbatim: vendor = {clean.get('vendor')!r}")

# %%
# 2. Build a multi-kind telemetry stream
# --------------------------------------
# A simulated 12-station network emits health and QC packets; we augment it
# with a few sync, power, and event packets so the stream exercises every
# packet kind.

result = simulate_iot_network(
    n_stations=12,
    profiles=["L1", "L3"],
    dropout_rate=0.04,
    survey_id="SIM",
    seed=5,
    detail=True,
)
packets = list(result["packets"])
t0 = min(p.timestamp for p in packets)
for offset, station in enumerate(result["stations"][:6]):
    device = station["device"]
    packets.append(
        TelemetryPacket.from_device(
            device,
            timestamp=t0 + 40 + 6 * offset,
            payload={
                "offset_ms": 0.4 + 0.2 * offset,
                "drift_ppm": 3 + offset,
                "gps_lock": offset != 5,
            },
            kind=PacketKind.SYNC,
            survey_id="SIM",
        )
    )
    packets.append(
        TelemetryPacket.from_device(
            device,
            timestamp=t0 + 42 + 6 * offset,
            payload={"battery_v": 12.6 - 0.1 * offset, "state": "ok"},
            kind=PacketKind.POWER,
            survey_id="SIM",
        )
    )
    if offset % 3 == 0:
        packets.append(
            TelemetryPacket.from_device(
                device,
                timestamp=t0 + 44 + 6 * offset,
                payload={
                    "event": "operator_note",
                    "severity": "info",
                    "message": "electrode re-watered",
                },
                kind=PacketKind.EVENT,
                survey_id="SIM",
            )
        )

kinds, counts = np.unique([p.kind.value for p in packets], return_counts=True)
print(
    f"{len(packets)} packets: "
    + ", ".join(f"{k}={c}" for k, c in zip(kinds, counts))
)

# %%
# 3. File transport round-trip
# ----------------------------
# ``FileTelemetryClient`` is a real, dependency-free transport: it appends
# each packet as one JSON line and can replay them. It doubles as the
# reference for offline logging and tests.

log_path = Path(tempfile.gettempdir()) / "iot_telemetry.jsonl"
if log_path.exists():
    log_path.unlink()
with FileTelemetryClient(str(log_path), dry_run=False) as client:
    for packet in packets:
        client.send(packet)
    healthy = client.healthcheck()
rows = FileTelemetryClient(str(log_path), dry_run=False).read_all()
print(
    f"wrote {len(packets)} packets -> {log_path.name}  "
    f"({log_path.stat().st_size / 1024:.1f} KiB), "
    f"read back {len(rows)}, healthcheck={healthy}"
)

# %%
# 4. The unified client interface (offline)
# -----------------------------------------
# Every protocol shares one interface and every client supports ``dry_run``,
# so the same code path works with or without a broker. Optional transport
# dependencies (paho-mqtt, pyserial, websocket-client) are imported only
# when a *real* connection is opened.

endpoints = {
    "file": str(log_path),
    "http": "https://collector.example/telemetry",
    "mqtt": "mqtt://broker.example:1883",
    "serial": "COM3",
    "websocket": "wss://collector.example/ws",
    "lora": None,
}
for proto, endpoint in endpoints.items():
    client = build_telemetry_client(proto, endpoint=endpoint, dry_run=True)
    ack = client.send(packets[0])
    print(
        f"  {proto:10s} send.ok={ack.ok!s:5s} protocol={ack.protocol:9s} "
        f"healthcheck={client.healthcheck()}"
    )

# %%
# The capability matrix documents what each transport implements. All
# support offline ``dry_run`` recording; the operations below are the live
# ones. Two transports need an optional dependency (noted on the axis).

CAP = {
    "File": dict(ops=[1, 1, 0, 0], dep=None),
    "HTTP": dict(ops=[1, 0, 0, 0], dep=None),
    "MQTT": dict(ops=[1, 1, 1, 1], dep="paho-mqtt"),
    "Serial": dict(ops=[1, 1, 0, 0], dep="pyserial"),
    "WebSocket": dict(ops=[1, 1, 0, 0], dep="websocket-client"),
}
op_labels = ["send", "receive", "subscribe", "listen"]
grid = np.array([CAP[t]["ops"] for t in CAP], dtype=float)
row_labels = [
    t if CAP[t]["dep"] is None else f"{t}\n(needs {CAP[t]['dep']})"
    for t in CAP
]

fig, ax = plt.subplots(figsize=(7.6, 4.6), constrained_layout=True)
ax.imshow(
    grid,
    cmap=plt.matplotlib.colors.ListedColormap(["#eceff1", "#009E73"]),
    vmin=0,
    vmax=1,
    aspect="auto",
)
for i in range(grid.shape[0]):
    for j in range(grid.shape[1]):
        ax.text(
            j,
            i,
            "✓" if grid[i, j] else "–",
            ha="center",
            va="center",
            color="white" if grid[i, j] else "#90a4ae",
            fontsize=14,
        )
ax.set_xticks(range(len(op_labels)), op_labels)
ax.set_yticks(range(len(row_labels)), row_labels)
ax.set_xticks(np.arange(-0.5, len(op_labels)), minor=True)
ax.set_yticks(np.arange(-0.5, len(row_labels)), minor=True)
ax.grid(which="minor", color="white", lw=2)
ax.tick_params(which="minor", length=0)
ax.set_title("Transport capability matrix (all support offline dry-run)")

# %%
# 5. Security: TLS and redacted credentials
# -----------------------------------------
# A :class:`~pycsamt.iot.SecurityConfig` carries TLS material and a
# credential. Secrets are **redacted** from ``repr`` and ``as_dict`` (so
# they never leak into a manifest or log) yet are still handed to the live
# client through ``client_options``.

security = SecurityConfig(
    tls=TLSConfig(enabled=True),
    credential=Credential(
        scheme=AuthScheme.BEARER, token="tok-SECRET-abc123"
    ),
    require_tls=True,
)
print("repr (safe):     ", repr(security.credential))
print("as_dict (safe):  ", security.as_dict()["credential"])
options = security.client_options()
secure_client = build_telemetry_client(
    "http",
    endpoint="https://collector.example/telemetry",
    dry_run=True,
    **options,
)
auth = secure_client.options.get("headers", {}).get("Authorization", "")
print(
    f"live client Authorization header present: {auth.startswith('Bearer ')}"
    f"  (token still hidden in logs)"
)

# %%
# 6. Replay the log into a session
# --------------------------------
# Reading the JSON-lines log back and feeding it to a fresh session
# reconstructs the deployment and re-runs QC — a loss-less round trip.

session = FieldSession("REPLAY-SIM")
for row in rows:
    session.add_packet(row)
status = session.assess()
print(
    f"replayed {session.n_packets} packets -> "
    f"{session.n_stations} stations, level={status.level.value}, "
    f"edge_acceptance_rate={status.edge_acceptance_rate:.2f}"
)

# %%
# Finally, the telemetry stream itself: each packet placed on a timeline by
# kind. Identity is carried by a legend, and the axis order groups the
# housekeeping kinds (health/sync/power/event) apart from QC.

KIND_ORDER = ["qc", "health", "sync", "power", "event"]
KIND_COLORS = {
    "qc": "#009E73",
    "health": "#56B4E9",
    "sync": "#E69F00",
    "power": "#CC79A7",
    "event": "#D55E00",
}
fig, ax = plt.subplots(figsize=(9.5, 4.2), constrained_layout=True)
present = [k for k in KIND_ORDER if k in set(p.kind.value for p in packets)]
for i, kind in enumerate(present):
    times = [(p.timestamp - t0) for p in packets if p.kind.value == kind]
    ax.scatter(
        times,
        np.full(len(times), i),
        s=42,
        color=KIND_COLORS[kind],
        edgecolor="white",
        linewidth=0.4,
        label=f"{kind} (n={len(times)})",
    )
ax.set_yticks(range(len(present)), present)
ax.set_xlabel("seconds since first packet")
ax.set_title("Telemetry stream by packet kind")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.legend(
    frameon=False,
    ncol=len(present),
    loc="upper center",
    bbox_to_anchor=(0.5, -0.18),
)

plt.show()
