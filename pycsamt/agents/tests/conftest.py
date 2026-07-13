# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Shared fixtures for pycsamt.agents tests.

Conventions
-----------
* All agents are instantiated without an api_key → no LLM calls, no cost.
* DL-dependent tests are guarded with ``_has_backend()``.
* ``edi_dir`` points to the bundled 3-EDI dataset; tests skip gracefully
  when the directory is absent (e.g. on minimal CI environments).
"""

from __future__ import annotations

from pathlib import Path

import pytest

# ── paths ─────────────────────────────────────────────────────────────────────

_AGENTS_TESTS = Path(__file__).resolve().parent  # agents/tests/
_PROJECT_ROOT = _AGENTS_TESTS.parents[3]  # repo root
_EDI_DIR = _PROJECT_ROOT / "data" / "3edis"


# ── backend helper ────────────────────────────────────────────────────────────


def has_backend() -> bool:
    """Return True when PyTorch or TensorFlow is installed."""
    try:
        from pycsamt.backends import get_backend

        return get_backend() != "none"
    except Exception:
        return False


# ── fixtures ──────────────────────────────────────────────────────────────────


@pytest.fixture(scope="session")
def project_root() -> Path:
    return _PROJECT_ROOT


@pytest.fixture(scope="session")
def edi_dir() -> Path:
    """Path to the bundled 3-EDI test dataset; skip when absent."""
    if not _EDI_DIR.exists():
        pytest.skip(f"3edis dataset not found: {_EDI_DIR}")
    return _EDI_DIR


@pytest.fixture(scope="session")
def edi_files(edi_dir: Path) -> list[Path]:
    """Sorted list of EDI files in the 3edis directory."""
    files = sorted(_EDI_DIR.glob("*.edi"))
    if not files:
        pytest.skip("No EDI files in 3edis directory.")
    return files


@pytest.fixture(scope="session")
def no_llm_kw() -> dict:
    """Keyword args for agents that disable LLM calls."""
    return {"api_key": None, "model": None, "llm_provider": "claude"}


@pytest.fixture(scope="session")
def loaded_sites(edi_dir: Path):
    """Sites object loaded from the 3edis directory (session-scoped)."""
    try:
        from pycsamt.agents import MTLoaderAgent

        ag = MTLoaderAgent()
        r = ag.execute({"path": str(edi_dir)})
        if r.status == "failed":
            pytest.skip(f"MTLoaderAgent failed: {r.error}")
        return r["sites"]
    except Exception as exc:
        pytest.skip(f"Could not load sites: {exc}")


@pytest.fixture()
def tmp_output(tmp_path: Path) -> Path:
    """A fresh temporary output directory."""
    d = tmp_path / "agent_output"
    d.mkdir()
    return d
