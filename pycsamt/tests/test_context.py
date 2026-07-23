from __future__ import annotations

import sys

from pycsamt.context import nullify_output


def test_nullify_output_suppresses_both_by_default():
    stdout, stderr = sys.stdout, sys.stderr
    with nullify_output():
        assert sys.stdout is not stdout
        assert sys.stderr is not stderr
        print("this should not raise")
        print("neither should this", file=sys.stderr)
    assert sys.stdout is stdout
    assert sys.stderr is stderr


def test_nullify_output_suppress_stdout_only():
    stdout, stderr = sys.stdout, sys.stderr
    with nullify_output(suppress_stdout=True, suppress_stderr=False):
        assert sys.stdout is not stdout
        assert sys.stderr is stderr
    assert sys.stdout is stdout
    assert sys.stderr is stderr


def test_nullify_output_suppress_stderr_only():
    stdout, stderr = sys.stdout, sys.stderr
    with nullify_output(suppress_stdout=False, suppress_stderr=True):
        assert sys.stdout is stdout
        assert sys.stderr is not stderr
    assert sys.stdout is stdout
    assert sys.stderr is stderr


def test_nullify_output_neither_suppressed():
    stdout, stderr = sys.stdout, sys.stderr
    with nullify_output(suppress_stdout=False, suppress_stderr=False):
        assert sys.stdout is stdout
        assert sys.stderr is stderr


def test_nullify_output_restores_on_exception():
    stdout, stderr = sys.stdout, sys.stderr
    try:
        with nullify_output():
            raise ValueError("boom")
    except ValueError:
        pass
    assert sys.stdout is stdout
    assert sys.stderr is stderr
