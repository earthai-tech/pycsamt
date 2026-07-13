"""Tests for tamper-evident provenance: HMAC manifest signing and the
QC-decision hash chain in :mod:`pycsamt.iot.provenance`.
"""

from __future__ import annotations

import copy
import json

from pycsamt.iot import (
    ProvenanceRecord,
    build_acquisition_manifest,
    hash_chain,
    sign_mapping,
    verify_hash_chain,
    verify_manifest,
    verify_signature,
)

KEY = "field-secret-2026"


def _manifest():
    m = build_acquisition_manifest(
        "SURV1",
        records=[ProvenanceRecord("S01", lat=6.5, lon=3.4)],
        method="csamt",
        operator="crew",
    )
    m.add_qc_decision(station="S01", decision="accept")
    m.add_qc_decision(station="S01", decision="reject", reasons=["low_snr"])
    return m


# ---------------------------------------------------------------------------
# generic HMAC sign / verify
# ---------------------------------------------------------------------------
def test_sign_and_verify_mapping():
    sig = sign_mapping({"a": 1, "b": 2}, KEY)
    # key ordering does not matter (canonical JSON)
    assert verify_signature({"b": 2, "a": 1}, sig, KEY)


def test_verify_rejects_altered_mapping():
    sig = sign_mapping({"a": 1, "b": 2}, KEY)
    assert not verify_signature({"a": 1, "b": 3}, sig, KEY)


def test_verify_rejects_wrong_key():
    sig = sign_mapping({"a": 1}, KEY)
    assert not verify_signature({"a": 1}, sig, "other-key")


def test_sign_accepts_bytes_key():
    sig = sign_mapping({"a": 1}, b"\x00\x01\x02")
    assert verify_signature({"a": 1}, sig, b"\x00\x01\x02")


# ---------------------------------------------------------------------------
# signed manifest envelope
# ---------------------------------------------------------------------------
def test_manifest_sign_envelope_shape():
    signed = _manifest().sign(KEY)
    assert set(signed) == {"manifest", "signature", "signature_algo"}
    assert signed["signature_algo"] == "hmac-sha256"


def test_verify_manifest_good_and_wrong_key():
    signed = _manifest().sign(KEY)
    assert verify_manifest(signed, KEY)
    assert not verify_manifest(signed, "nope")


def test_verify_manifest_detects_tampering():
    signed = _manifest().sign(KEY)
    tampered = copy.deepcopy(signed)
    tampered["manifest"]["operator"] = "attacker"
    assert not verify_manifest(tampered, KEY)


def test_verify_manifest_malformed_envelope():
    assert not verify_manifest({"signature": "x"}, KEY)  # no manifest
    assert not verify_manifest({"manifest": {}}, KEY)  # no signature


def test_write_signed_roundtrip(tmp_path):
    path = _manifest().write_signed(str(tmp_path / "m.signed.json"), KEY)
    loaded = json.loads(open(path).read())
    assert verify_manifest(loaded, KEY)


# ---------------------------------------------------------------------------
# hash chain
# ---------------------------------------------------------------------------
def test_hash_chain_links_and_verifies():
    entries = [
        {"decision": "accept"},
        {"decision": "reject"},
        {"decision": "accept"},
    ]
    chain = hash_chain(entries)
    assert [c["seq"] for c in chain] == [0, 1, 2]
    assert chain[0]["prev_hash"] == ""
    assert chain[1]["prev_hash"] == chain[0]["entry_hash"]
    assert verify_hash_chain(chain)


def test_hash_chain_detects_entry_tampering():
    chain = hash_chain([{"decision": "accept"}, {"decision": "reject"}])
    bad = copy.deepcopy(chain)
    bad[0]["decision"] = "reject"  # flip a verdict
    assert not verify_hash_chain(bad)


def test_hash_chain_detects_reordering():
    chain = hash_chain([{"n": 1}, {"n": 2}, {"n": 3}])
    assert not verify_hash_chain(list(reversed(chain)))


def test_manifest_chained_qc_decisions():
    chain = _manifest().chained_qc_decisions()
    assert len(chain) == 2
    assert verify_hash_chain(chain)
    # the chain preserves the original decision payload
    assert chain[0]["decision"] == "accept"
    assert chain[1]["decision"] == "reject"


def test_empty_chain_verifies():
    assert verify_hash_chain(hash_chain([]))
