# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for :mod:`pycsamt.format.provenance` — ML-model provenance for a
shared AI/DL inversion result."""

from __future__ import annotations

import hashlib

import pytest

from pycsamt.format.provenance import ModelProvenance, compute_checkpoint_hash


class TestModelProvenance:
    def test_defaults_are_empty_and_free_form(self):
        prov = ModelProvenance()
        d = prov.to_dict()
        assert d["architecture"] == ""
        assert d["hyperparameters"] == {}
        assert d["authors"] == []
        assert d["random_seed"] is None

    def test_to_dict_round_trips_given_fields(self):
        prov = ModelProvenance(
            architecture="GCN",
            framework="pytorch",
            framework_version="2.3.0",
            checkpoint="gcn_v1.pt",
            checkpoint_sha256="deadbeef",
            training_data="synthetic_survey_v2",
            hyperparameters={"lr": 1e-3, "layers": 4},
            random_seed=42,
            git_commit="abc123",
            authors=["A. Researcher"],
            contact="a@example.org",
            notes="trained on synthetic MT data",
        )
        d = prov.to_dict()
        assert d["architecture"] == "GCN"
        assert d["hyperparameters"] == {"lr": 1e-3, "layers": 4}
        assert d["random_seed"] == 42
        assert d["authors"] == ["A. Researcher"]

    def test_random_seed_must_be_int_or_none(self):
        prov = ModelProvenance(random_seed="not-an-int")
        with pytest.raises(TypeError):
            prov.validate()

    def test_mutable_defaults_are_independent_across_instances(self):
        a = ModelProvenance()
        b = ModelProvenance()
        a.hyperparameters["lr"] = 1e-3
        assert b.hyperparameters == {}


class TestComputeCheckpointHash:
    def test_matches_hashlib_reference(self, tmp_path):
        path = tmp_path / "weights.bin"
        payload = b"pretend-model-weights" * 1000
        path.write_bytes(payload)

        expected = hashlib.sha256(payload).hexdigest()
        assert compute_checkpoint_hash(path) == expected

    def test_streams_in_chunks_without_loading_whole_file(self, tmp_path):
        path = tmp_path / "weights.bin"
        payload = b"\x00\x01\x02\x03" * 2048
        path.write_bytes(payload)

        expected = hashlib.sha256(payload).hexdigest()
        assert compute_checkpoint_hash(path, chunk_size=16) == expected
