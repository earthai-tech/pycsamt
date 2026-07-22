"""
Scratch coverage runner -- NOT a permanent test file.

Works around a reproducible Windows access-violation crash in this conda
env: importing torch's C extension (torch/__init__.py -> torch._C) while
coverage.py's tracer is already active corrupts memory (confirmed to
reproduce even on already-passing, pre-existing test files such as
test_web_callbacks_lines.py -- nothing to do with the new pipeline/profile
tests). Pre-importing torch before coverage.start() avoids the crash since
the C extension is then already resident in sys.modules by the time
tracing begins.
"""
import os
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import torch  # noqa: F401  (pre-load before coverage tracing starts)

import coverage

cov = coverage.Coverage(
    source=[
        "pycsamt.app.web.callbacks.pipeline",
        "pycsamt.app.web.callbacks.profile",
    ],
)
cov.start()

import pytest

rc = pytest.main(
    [
        "pycsamt/app/tests/test_web_callbacks_pipeline.py",
        "pycsamt/app/tests/test_web_callbacks_profile.py",
        "-q",
    ]
)

cov.stop()
cov.save()
cov.report(show_missing=True)

sys.exit(rc)
