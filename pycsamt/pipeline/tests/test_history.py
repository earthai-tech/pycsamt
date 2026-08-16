"""Tests for pycsamt.pipeline._history — append_run, load_history.

Public entry point: :mod:`pycsamt.pipeline` (``Pipeline.run(history=...)``,
``load_history``)
Implementation:     :mod:`pycsamt.pipeline._history`
"""

from __future__ import annotations

import json

import pytest

from pycsamt.api.pipe import reset_pipe
from pycsamt.pipeline import Pipeline, Step, load_history, lookup_step
from pycsamt.pipeline._history import append_run, default_history_path


class _CountingSites:
    def __init__(self, n: int = 3, count: int = 0):
        self._n = n
        self.count = count

    def __len__(self) -> int:
        return self._n


class _IdentityStep(Step):
    def transform(self, sites):
        return sites

    def generate_qc_plots(self, sites):
        return []


def _identity(code: str) -> _IdentityStep:
    step = _IdentityStep.__new__(_IdentityStep)
    step.spec = lookup_step(code)
    step.params = dict(step.spec.defaults)
    return step


@pytest.fixture(autouse=True)
def _reset_pipe_cfg():
    reset_pipe()
    yield
    reset_pipe()


@pytest.fixture()
def sites() -> _CountingSites:
    return _CountingSites(n=4, count=0)


@pytest.fixture()
def simple_pipe() -> Pipeline:
    return Pipeline([("notch", _identity("NR001")), ("band", _identity("FREQ001"))])


class TestAppendAndLoadRoundTrip:
    def test_history_false_writes_nothing(self, simple_pipe, sites, tmp_path):
        path = tmp_path / "hist.jsonl"
        simple_pipe.run(sites, outdir=None, save_report=False, save_plots=False)
        assert not path.exists()

    def test_history_true_writes_to_default_location(
        self, simple_pipe, sites, tmp_path, monkeypatch
    ):
        default_path = tmp_path / "default_history.jsonl"
        monkeypatch.setattr(
            "pycsamt.pipeline._history.default_history_path", lambda: default_path
        )
        simple_pipe.run(
            sites, outdir=None, save_report=False, save_plots=False, history=True
        )
        assert default_path.exists()
        assert len(load_history(default_path)) == 1

    def test_history_path_writes_to_explicit_file(self, simple_pipe, sites, tmp_path):
        path = tmp_path / "custom.jsonl"
        simple_pipe.run(
            sites, outdir=None, save_report=False, save_plots=False, history=path
        )
        assert path.exists()
        records = load_history(path)
        assert len(records) == 1

    def test_two_runs_append_two_lines(self, simple_pipe, sites, tmp_path):
        path = tmp_path / "hist.jsonl"
        simple_pipe.run(
            sites, outdir=None, save_report=False, save_plots=False, history=path
        )
        simple_pipe.run(
            sites, outdir=None, save_report=False, save_plots=False, history=path
        )
        assert len(load_history(path)) == 2
        # real append, not overwrite: two distinct JSON lines on disk
        lines = [l for l in path.read_text(encoding="utf-8").splitlines() if l.strip()]
        assert len(lines) == 2
        for line in lines:
            json.loads(line)  # each line individually parses

    def test_record_shape(self, simple_pipe, sites, tmp_path):
        path = tmp_path / "hist.jsonl"
        simple_pipe.run(
            sites, outdir=None, save_report=False, save_plots=False, history=path
        )
        record = load_history(path)[0]
        assert record["pipeline_name"] == simple_pipe.name
        assert record["ok"] is True
        assert record["n_errors"] == 0
        assert record["n_sites_in"] == 4
        assert record["n_sites_out"] == 4
        assert record["outdir"] is None
        assert [s["code"] for s in record["steps"]] == ["NR001", "FREQ001"]
        assert all("cached" in s for s in record["steps"])
        assert "timestamp" in record


class TestLoadHistory:
    def test_missing_file_returns_empty_list(self, tmp_path):
        assert load_history(tmp_path / "does_not_exist.jsonl") == []

    def test_last_limits_to_most_recent(self, tmp_path):
        path = tmp_path / "hist.jsonl"
        for i in range(5):
            path.write_text(
                (path.read_text(encoding="utf-8") if path.exists() else "")
                + json.dumps({"pipeline_name": f"run{i}"})
                + "\n",
                encoding="utf-8",
            )
        records = load_history(path, last=2)
        assert [r["pipeline_name"] for r in records] == ["run3", "run4"]

    def test_corrupt_line_is_skipped_not_raised(self, tmp_path):
        path = tmp_path / "hist.jsonl"
        path.write_text(
            json.dumps({"pipeline_name": "good1"})
            + "\n"
            + "{this is not valid json\n"
            + json.dumps({"pipeline_name": "good2"})
            + "\n",
            encoding="utf-8",
        )
        records = load_history(path)
        assert [r["pipeline_name"] for r in records] == ["good1", "good2"]

    def test_default_path_is_under_home_pycsamt(self):
        assert default_history_path().parts[-2:] == (".pycsamt", "pipeline_history.jsonl")


class TestAppendRunFailureIsolation:
    def test_run_still_succeeds_if_history_write_fails(
        self, simple_pipe, sites, tmp_path, monkeypatch
    ):
        import warnings

        def _boom(path, result):
            raise OSError("disk full (simulated)")

        # append_run is imported locally inside run() (from ._history import
        # append_run, ...), executed fresh each call — patching it on the
        # _history module itself is what that local import actually resolves.
        monkeypatch.setattr("pycsamt.pipeline._history.append_run", _boom)

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = simple_pipe.run(
                sites,
                outdir=None,
                save_report=False,
                save_plots=False,
                history=tmp_path / "hist.jsonl",
            )
        assert result.ok  # the pipeline run itself must still succeed
        assert any("Failed to append run history" in str(w.message) for w in caught)
