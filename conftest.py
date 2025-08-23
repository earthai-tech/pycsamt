# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
Configuration for the pytest suite.

This file defines shared fixtures for locating and accessing test
data files, ensuring that tests can be run consistently from any
directory.
"""

from __future__ import annotations

from pathlib import Path
import pytest


def get_project_root() -> Path:
    """Find the repo root by locating the 'pycsamt' dir."""
    cur = Path(__file__).resolve()
    for parent in cur.parents:
        if (parent / "pycsamt").is_dir():
            return parent
    raise FileNotFoundError("Project root not found.")


@pytest.fixture(scope="session")
def project_root() -> Path:
    """Project root directory."""
    return get_project_root()


@pytest.fixture(scope="session")
def data_path(project_root: Path) -> Path:
    """Base path to bundled test data."""
    return project_root / "data" / "avg"


@pytest.fixture(scope="session")
def legacy_data_file(data_path: Path) -> Path:
    """
    Path to legacy K1.AVG file; skip if missing.
    """
    p = data_path / "K1.AVG"
    if not p.exists():
        pytest.skip(f"Missing legacy data: {p}")
    return p


@pytest.fixture(scope="session")
def modern_data_file(data_path: Path) -> Path:
    """
    Path to modern K2.AVG file; skip if missing.
    """
    p = data_path / "K2.AVG"
    if not p.exists():
        pytest.skip(f"Missing modern data: {p}")
    return p


@pytest.fixture(scope="session")
def stn_file_k1(data_path: Path) -> Path:
    """
    Path to K1.stn topography file; skip if missing.
    """
    p = data_path / "K1.stn"
    if not p.exists():
        pytest.skip(f"Missing STN file: {p}")
    return p

@pytest.fixture(scope="session")
def stn_file_k2(data_path: Path) -> Path:
    """
    Path to K2.stn topography file; skip if missing.
    """
    p = data_path / "K2.stn"
    if not p.exists():
        pytest.skip(f"Missing STN file: {p}")
    return p


