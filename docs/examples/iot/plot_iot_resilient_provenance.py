r"""
Resilient telemetry and a tamper-evident provenance trail
=========================================================

Two things separate a field survey you can publish from one you cannot:
the telemetry has to survive a flaky uplink without losing packets, and the
audit trail has to be provably untampered. This example exercises both on
the **real** WILLY AMT line (``data/AMT/WILLY_DATA/L18PLT``): it streams a
QC packet per station through a link that drops out mid-line, buffers the
backlog with :class:`~pycsamt.iot.StoreAndForwardClient`, then signs the
acquisition manifest and hash-chains the QC decisions so any later edit is
detectable.

The workflow is:

1. **seed** a field session from the archived survey (real stations);
2. **stream** a QC packet per station over an intermittent link, buffering
   through the outage and flushing when it recovers;
3. **sign** the acquisition manifest (HMAC) and verify it, including a
   tamper check; and
4. **hash-chain** the QC decisions so reordering or edits are caught.
"""

# sphinx_gallery_thumbnail_number = 1

from __future__ import annotations

import copy
import os
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from pycsamt.iot import (
    StoreAndForwardClient,
    TelemetryPacket,
    build_acquisition_manifest,
    field_session_from_edis,
    hash_chain,
    impedance_to_z,
    log_qc_decision,
    verify_hash_chain,
    verify_manifest,
    z_to_edi,
)
from pycsamt.iot.protocols.base import (
    BaseTelemetryClient,
    IoTProtocol,
    TelemetryAck,
    TelemetryError,
    _coerce_packet,
)

warnings.filterwarnings("ignore")

DELIVERED, BUFFERED, FLUSHED = "#009E73", "#D55E00", "#0072B2"


def style_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, color="#000000", alpha=0.07, lw=0.7)
    ax.set_axisbelow(True)


def repo_root() -> Path:
    env_root = os.environ.get("PYCSAMT_DOCS_REPO_ROOT")
    if env_root:
        return Path(env_root)
    return Path(__file__).resolve().parents[3]


class IntermittentLink(BaseTelemetryClient):
    """A stand-in cellular/LoRa uplink that can be toggled offline.

    Wrapping any client like this in a :class:`StoreAndForwardClient` is the
    whole point: the transport can fail, and no packet is lost.
    """

    protocol = IoTProtocol.HTTP

    def __init__(self):
        super().__init__("https://hub/ingest", dry_run=False)
        self.online = True
        self.delivered: list[TelemetryPacket] = []

    def send(self, packet):
        pkt = _coerce_packet(packet)
        if not self.online:
            raise TelemetryError("uplink down")
        self.delivered.append(pkt)
        return TelemetryAck(ok=True, protocol="http", packet_id="x")


# %%
# 1. Seed a session from the archived survey
# ------------------------------------------
# The reverse bridge rebuilds the real station list; a synthetic stand-in
# keeps the example runnable on a checkout without the WILLY data.

survey_dir = repo_root() / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
if not (survey_dir.is_dir() and any(survey_dir.glob("*.edi"))):  # pragma: no cover
    import tempfile

    print("NOTE: bundled WILLY survey absent - using a synthetic stand-in.")
    tmp = Path(tempfile.mkdtemp())
    freq = np.logspace(4, 0, 24)
    for i in range(12):
        z = impedance_to_z((1 + 1j) * np.sqrt(freq), freq, station=f"S{i:02d}")
        z_to_edi(z, station=f"S{i:02d}", lat=32.12 + 0.002 * i,
                 lon=119.13, elevation=20.0, savepath=str(tmp))
    survey_dir = tmp

session = field_session_from_edis(
    str(survey_dir), survey_id="WILLY-L18", method="amt", operator="crew",
)
stations = [s.station_id for s in session.to_sites()]
print(f"session: {len(stations)} stations from {survey_dir.name}")

# %%
# 2. Stream QC packets through an intermittent link
# -------------------------------------------------
# The link is offline for a block of stations in the middle of the line.
# Every send during the outage is queued (never dropped); the moment the
# link returns, ``flush`` drains the backlog in order. The queue is spooled
# to disk, so even a logger reboot mid-outage would keep the backlog.

spool = Path(os.environ.get("TMPDIR", "/tmp")) / "willy_l18_spool.jsonl"
if spool.exists():
    spool.unlink()
link = IntermittentLink()
buffered_client = StoreAndForwardClient(
    link, spool_path=str(spool), max_queue=10_000
)

outage = range(len(stations) // 3, 2 * len(stations) // 3)  # middle third
timeline = []  # (station_index, state, pending)
for i, station in enumerate(stations):
    link.online = i not in outage
    packet = TelemetryPacket(
        device_id=f"{station}-node", timestamp=1_700_000_000.0 + i,
        topic=f"pycsamt/WILLY-L18/{station}/qc", kind="qc",
        payload={"station": station, "accepted": True,
                 "channels": ["ex", "ey", "hx", "hy"]},
    )
    before = buffered_client.pending
    ack = buffered_client.send(packet)
    state = "delivered" if ack.ok else "buffered"
    # if we just came back online, drain whatever had piled up
    flushed_now = 0
    if link.online and buffered_client.pending:
        flushed_now = buffered_client.flush()
    timeline.append((i, state, buffered_client.pending, flushed_now))

n_delivered_live = sum(1 for _, s, _, _ in timeline if s == "delivered")
n_buffered = sum(1 for _, s, _, _ in timeline if s == "buffered")
print(f"live-delivered: {n_delivered_live}, buffered during outage: "
      f"{n_buffered}")
print(f"finally delivered: {len(link.delivered)}/{len(stations)} "
      f"(pending: {buffered_client.pending})")
assert len(link.delivered) == len(stations)  # nothing lost

# %%
# 3. Sign the acquisition manifest
# --------------------------------
# The manifest records every station's provenance and QC decision. An HMAC
# signature makes it tamper-evident: a reviewer holding the survey key can
# confirm nothing changed between the field and the archive.

decisions = [
    log_qc_decision(station=s, decision="accept", operator="crew")
    for s in stations
]
manifest = build_acquisition_manifest(
    "WILLY-L18", records=list(session.to_manifest().records),
    qc_decisions=decisions, method="amt", operator="crew",
)
key = "willy-2026-survey-key"
signed = manifest.sign(key)
print(f"manifest signed ({signed['signature_algo']})")
print(f"verify (correct key): {verify_manifest(signed, key)}")
print(f"verify (wrong key):   {verify_manifest(signed, 'guessed')}")

tampered = copy.deepcopy(signed)
tampered["manifest"]["operator"] = "impostor"
print(f"verify (tampered):    {verify_manifest(tampered, key)}")

# %%
# 4. Hash-chain the QC decisions
# ------------------------------
# A hash chain over the QC log makes the *sequence* tamper-evident too:
# altering, inserting, or reordering any decision breaks the chain from
# that point on.

chain = hash_chain(decisions)
print(f"QC decision chain: {len(chain)} links, intact: "
      f"{verify_hash_chain(chain)}")
reordered = list(reversed(copy.deepcopy(chain)))
print(f"reordered chain verifies: {verify_hash_chain(reordered)}")

# %%
# 5. The resilience and integrity picture
# ---------------------------------------
# Top: the telemetry timeline -- green packets went out live, orange ones
# were buffered through the outage (grey band), and the blue line is the
# backlog draining to zero the instant the link returns. Bottom: the
# provenance integrity checks, all behaving as they must.

idx = np.array([t[0] for t in timeline])
pending = np.array([t[2] for t in timeline])
states = [t[1] for t in timeline]

fig, (ax_t, ax_i) = plt.subplots(
    2, 1, figsize=(10.5, 7.0), constrained_layout=True,
    gridspec_kw={"height_ratios": [3, 2]},
)

if outage:
    ax_t.axvspan(min(outage) - 0.5, max(outage) + 0.5, color="#888",
                 alpha=0.15, lw=0, label="uplink outage")
for state, colour in (("delivered", DELIVERED), ("buffered", BUFFERED)):
    m = np.array([s == state for s in states])
    if m.any():
        ax_t.scatter(idx[m], np.zeros(m.sum()), c=colour, s=60, zorder=3,
                     edgecolor="#222", linewidth=0.5, label=state)
ax_t.plot(idx, pending, "-", color=FLUSHED, lw=1.8, label="buffered backlog")
ax_t.set(xlabel="station index along the line", ylabel="packets buffered",
         title="Store-and-forward telemetry through an uplink outage")
ax_t.legend(frameon=False, ncol=2, loc="upper left")
style_axis(ax_t)

checks = [
    ("signature\n(correct key)", verify_manifest(signed, key)),
    ("wrong key\nrejected", not verify_manifest(signed, "guessed")),
    ("tamper\ndetected", not verify_manifest(tampered, key)),
    ("QC chain\nintact", verify_hash_chain(chain)),
    ("reorder\ncaught", not verify_hash_chain(reordered)),
    ("no packets\nlost", len(link.delivered) == len(stations)),
]
labels = [c[0] for c in checks]
passed = [1.0 if c[1] else 0.0 for c in checks]
bar_colours = [DELIVERED if p else BUFFERED for p in passed]
ax_i.bar(range(len(checks)), passed, color=bar_colours, alpha=0.9)
ax_i.set_xticks(range(len(checks)))
ax_i.set_xticklabels(labels, fontsize=8)
ax_i.set(ylim=(0, 1.25), ylabel="pass", yticks=[0, 1],
         title="Provenance integrity and delivery guarantees")
for i, p in enumerate(passed):
    ax_i.text(i, p + 0.06, "PASS" if p else "FAIL", ha="center",
              fontsize=8, color="#222", fontweight="bold")
style_axis(ax_i)

plt.show()
