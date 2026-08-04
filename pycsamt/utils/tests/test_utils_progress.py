"""Tests for :mod:`pycsamt.utils.progress`."""

from __future__ import annotations

import pytest

from pycsamt.utils.progress import ProgressBar, normalize_verbose, progress


class _FakeTTYStream:
    def isatty(self) -> bool:
        return True


class _FakeNonTTYStream:
    def isatty(self) -> bool:
        return False


# ── normalize_verbose ───────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "value,expected",
    [
        (True, 1),
        (False, 0),
        (0, 0),
        (1, 1),
        (2, 2),
        (5, 2),
        (-3, 0),
        ("silent", 0),
        ("off", 0),
        ("bar", 1),
        ("on", 1),
        ("log", 2),
        ("lines", 2),
    ],
)
def test_normalize_verbose_scalars(value, expected):
    assert normalize_verbose(value) == expected


def test_normalize_verbose_none_uses_default():
    assert normalize_verbose(None, default=2) == 2
    assert normalize_verbose(None) == 1


def test_normalize_verbose_auto_resolves_from_stream():
    assert normalize_verbose("auto", stream=_FakeTTYStream()) == 1
    assert normalize_verbose("auto", stream=_FakeNonTTYStream()) == 2


def test_normalize_verbose_rejects_bad_string():
    with pytest.raises(ValueError):
        normalize_verbose("not-a-mode")


# ── ProgressBar ─────────────────────────────────────────────────────────────


def test_progress_bar_silent_emits_nothing(capsys):
    with ProgressBar(total=5, desc="x", verbose=0) as bar:
        for i in range(5):
            bar.update(1)
    captured = capsys.readouterr()
    assert captured.out == ""


def test_progress_bar_log_mode_prints_lines(capsys):
    with ProgressBar(total=10, desc="samples", unit="sample", verbose=2) as bar:
        for _ in range(10):
            bar.update(1)
    captured = capsys.readouterr()
    lines = [ln for ln in captured.out.splitlines() if ln.strip()]
    assert lines, "expected at least one log line"
    assert "samples" in lines[0]
    assert "10 sample" in lines[-1] or "10/10 sample" in lines[-1]


def test_progress_bar_log_mode_respects_log_every(capsys):
    with ProgressBar(total=20, verbose=2, log_every=5) as bar:
        for _ in range(20):
            bar.update(1)
    captured = capsys.readouterr()
    lines = [ln for ln in captured.out.splitlines() if ln.strip()]
    # updates at 5, 10, 15, 20 -> 4 lines
    assert len(lines) == 4


def test_progress_bar_log_mode_with_metrics(capsys):
    with ProgressBar(total=1, verbose=2, log_every=1) as bar:
        bar.update(1, metrics={"loss": 0.12345})
    captured = capsys.readouterr()
    assert "loss: 0.1234" in captured.out or "loss: 0.1235" in captured.out


def test_progress_bar_bar_mode_uses_tqdm_and_does_not_raise():
    pytest.importorskip("tqdm")
    with ProgressBar(total=3, desc="bar-mode", verbose=1, leave=False) as bar:
        for i in range(3):
            bar.update(1, metrics={"loss": float(i)})
    # No assertion on stderr formatting (tqdm owns that) — just must not raise.


def test_progress_bar_close_is_idempotent():
    bar = ProgressBar(total=1, verbose=0)
    bar.close()
    bar.close()  # must not raise


def test_progress_bar_wrap_iterates_all_items():
    items = list(range(7))
    seen = list(ProgressBar.wrap(items, verbose=0))
    assert seen == items


def test_progress_function_is_progressbar_wrap_shorthand():
    items = list(range(4))
    seen = list(progress(items, verbose=0, desc="x"))
    assert seen == items


def test_progress_bar_exception_inside_context_still_closes():
    bar_holder = {}

    class _Recording(ProgressBar):
        def close(self):
            bar_holder["closed"] = True
            super().close()

    with pytest.raises(RuntimeError):
        with _Recording(total=2, verbose=0) as bar:
            bar.update(1)
            raise RuntimeError("boom")
    assert bar_holder.get("closed") is True
