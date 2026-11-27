# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Output folders of the pipeline tests that write real results.

These tests keep their results under ``<repo>/pipe_results/`` so they
can be inspected after a run, and each clears its folder before
running. Under pytest-xdist the ``xdist_group`` pins were not enough:
two workers ran the same fixture at once (CI, Python 3.10 workflows
shard), and one worker's ``shutil.rmtree`` deleted the other's plots
mid-run, so the files on disk no longer matched the recorded result.
Giving every worker its own folder removes the race whatever the
scheduling.
"""

from __future__ import annotations

import os
from pathlib import Path

__all__ = ["PROJECT_ROOT", "pipe_results_dir"]

PROJECT_ROOT = Path(__file__).resolve().parents[3]


def pipe_results_dir(name: str) -> Path:
    """Folder ``pipe_results/<name>`` for this test process.

    Parameters
    ----------
    name : str
        Folder name of the test, for example ``"integration_test"``.

    Returns
    -------
    Path
        ``<repo>/pipe_results/<name>`` when tests run in one process,
        ``<repo>/pipe_results/<name>_<worker>`` under pytest-xdist
        (``<worker>`` is ``PYTEST_XDIST_WORKER``, e.g. ``gw0``).

    Examples
    --------
    >>> from pycsamt.pipeline.tests._outdirs import pipe_results_dir
    >>> path = pipe_results_dir("demo")
    >>> path.parent.name, path.name.startswith("demo")
    ('pipe_results', True)
    """
    worker = os.environ.get("PYTEST_XDIST_WORKER")
    folder = f"{name}_{worker}" if worker else name
    return PROJECT_ROOT / "pipe_results" / folder
