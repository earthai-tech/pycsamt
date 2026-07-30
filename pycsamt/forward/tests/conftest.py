# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Test-only fixtures for :mod:`pycsamt.forward.tests`."""

from __future__ import annotations

import os

# The repo-root conftest.py unconditionally pre-imports torch (for an
# unrelated coverage-tracing crash, see its own comment). On this
# Windows/conda environment, torch's bundled Intel OpenMP runtime then
# conflicts with the one scipy.sparse.linalg.spsolve pulls in for a
# large-enough sparse factorization -- e.g. mt3d.py's 3-D curl-curl
# solves -- aborting the whole interpreter ("OMP: Error #15") rather
# than raising a catchable Python exception. Same root cause already
# documented and fixed the same way in pycsamt/app/tests/conftest.py;
# must be set before scipy's OpenMP runtime actually initializes.
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
