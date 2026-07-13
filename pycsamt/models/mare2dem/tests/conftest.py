# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Shared pytest fixtures for MARE2DEM tests.

All tests in this package skip gracefully when the
``data/mare2dem/`` example directory is not present (e.g. in CI
environments where large binary data is not checked in).
"""

from __future__ import annotations

from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Root of the example data directory
# ---------------------------------------------------------------------------

_REPO_ROOT = Path(__file__).parents[4]  # pycsamt/
_DATA_DIR = _REPO_ROOT / "data" / "mare2dem"

# Specific example sub-directories
HILL_DIR = _DATA_DIR / "hill"
CSEM_DIR = _DATA_DIR / "demo_csem"
INVERSION_DIR = _DATA_DIR / "demo_mt_inversion"
CSEQ_MT_DIR = _DATA_DIR / "demo_csem_mt"


def _skip_if_missing(path: Path) -> pytest.MarkDecorator:
    return pytest.mark.skipif(
        not path.exists(),
        reason=f"Example data not found: {path}. Copy from mare2dem_examples/.",
    )


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def hill_dir() -> Path:
    if not HILL_DIR.exists():
        pytest.skip(f"Hill example data not found: {HILL_DIR}")
    return HILL_DIR


@pytest.fixture(scope="session")
def csem_dir() -> Path:
    if not CSEM_DIR.exists():
        pytest.skip(f"CSEM example data not found: {CSEM_DIR}")
    return CSEM_DIR


@pytest.fixture(scope="session")
def inversion_dir() -> Path:
    if not INVERSION_DIR.exists():
        pytest.skip(f"MT inversion example not found: {INVERSION_DIR}")
    return INVERSION_DIR


@pytest.fixture(scope="session")
def csem_mt_dir() -> Path:
    if not CSEQ_MT_DIR.exists():
        pytest.skip(f"CSEM+MT example not found: {CSEQ_MT_DIR}")
    return CSEQ_MT_DIR
