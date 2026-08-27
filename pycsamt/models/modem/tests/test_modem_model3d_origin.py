# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for ModEmModel3D.origin/.rotation — the WS-format trailer that
real ModEM output carries (grid centre + rotation) after the
resistivity volume, previously parsed and then silently discarded.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.models.modem.model3d import ModEmModel3D

_DATA_DIR = (
    Path(__file__).parents[4] / "data" / "modem" / "willy_27freq_watex_line02_sample"
)
_REAL_RHO = _DATA_DIR / "27-freq-run-watex01.rho"
_SKIP = pytest.mark.skipif(
    not _REAL_RHO.exists(), reason=f"bundled ModEM data not found: {_REAL_RHO}"
)


@_SKIP
def test_real_file_origin_and_rotation_are_parsed():
    model = ModEmModel3D.read(_REAL_RHO)
    np.testing.assert_allclose(model.origin, [-4509.828, -6752.725, 0.0])
    assert model.rotation == pytest.approx(0.0)


def test_default_origin_and_rotation_are_zeros_not_none():
    # Matches read_mackie3d's own default exactly, so every existing
    # consumer of these attributes (iotools.export.write_meshtools3d,
    # iotools.interpolate) keeps working via their
    # ``getattr(model, "origin", [0.0, 0.0, 0.0])`` fallback.
    model = ModEmModel3D()
    np.testing.assert_array_equal(model.origin, [0.0, 0.0, 0.0])
    assert model.rotation == 0.0


def test_write_meshtools3d_still_works_with_a_real_origin(tmp_path):
    # Regression guard for the actual bug this fix introduced and then
    # fixed: export.py's write_meshtools3d does
    # ``np.asarray(getattr(model, "origin", [0, 0, 0])) * 1000.0`` and
    # unpacks the result as (ox, oy, oz) -- that unpack raised
    # TypeError when origin defaulted to None instead of zeros(3).
    from pycsamt.models.modem.iotools.export import write_meshtools3d

    model = ModEmModel3D()
    model.x_widths = np.array([100.0, 100.0])
    model.y_widths = np.array([100.0, 100.0])
    model.z_widths = np.array([50.0, 50.0])
    model.n_air = 0
    model.rho_loge = np.zeros((2, 2, 2))
    model.origin = np.array([1.5, 2.5, 0.0])

    msh_path, con_path = write_meshtools3d(model, tmp_path / "out")
    assert msh_path.exists() and con_path.exists()
