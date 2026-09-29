"""full_processing on real data must keep the data (2026-09-25).

Found from a desktop run on the Gabbs Valley EMTF-XML files (USGS,
https://doi.org/10.5066/P9GZ9Z56): the preset "finished OK" but wrote
all-NaN impedances.  Causes, each tested below:

* FREQ004 align_grid interpolated on descending frequencies (EDI/XML
  order) -> every target frequency took one of the two end values;
* SK001/SK002 used a Bahr-skew threshold (0.3) on the phase-tensor beta in
  degrees -> 96-100 % of rows masked; unknown skew masked whole stations;
* the skew lookup searched descending periods;
* TZ001 rotated by a NaN strike -> NaN tensor and errors;
* rotated errors used an element-wise (wrong) propagation.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[3]
GV_XML = ROOT / "data" / "gv_data" / "xml"


def test_interp_complex_handles_descending_frequencies():
    from pycsamt.emtools.frequency import _interp_complex

    x = np.log10(np.array([100.0, 10.0, 1.0, 0.1]))  # EDI order
    y = np.array([1, 2, 3, 4], dtype=complex)
    new = np.log10(np.array([0.1, 1.0, 10.0, 100.0]))
    assert np.allclose(_interp_complex(x, y, new, method="nearest"),
                       [4, 3, 2, 1])
    mid = np.log10(np.array([np.sqrt(10.0)]))
    assert np.allclose(_interp_complex(x, y, mid, method="linear"), 2.5)


def test_skew_defaults_are_phase_tensor_degrees():
    from pycsamt.pipeline._registry import lookup_step

    assert lookup_step("SK001").defaults["thresh"] == 3.0
    assert lookup_step("SK002").defaults["thresh"] == 3.0


def test_rotate_nan_angle_is_noop_and_errors_propagate():
    from types import SimpleNamespace

    from pycsamt.site.edit import rotate

    z = np.array([[[1 + 1j, 2 + 0j], [3 + 0j, 4 + 1j]]])
    err = np.array([[[0.1, 0.2], [0.3, 0.4]]])
    ed = SimpleNamespace(Z=SimpleNamespace(z=z.copy(), z_err=err.copy()))
    rotate(ed, float("nan"), inplace=True)
    assert np.array_equal(ed.Z.z, z) and np.array_equal(ed.Z.z_err, err)
    rotate(ed, 90.0, inplace=True)
    # 90 deg swaps x and y: sigma_xx' = sigma_yy etc.
    assert np.allclose(ed.Z.z_err[0], [[0.4, 0.3], [0.2, 0.1]], atol=1e-12)


@pytest.mark.skipif(not GV_XML.is_dir(), reason="gv_data not present")
def test_full_processing_keeps_real_xml_data(tmp_path):
    import matplotlib

    matplotlib.use("Agg")
    from pycsamt.emtools import ensure_sites
    from pycsamt.pipeline import Pipeline

    p = Pipeline.from_preset("full_processing")
    S = ensure_sites(str(GV_XML))
    for item in p._steps:
        st = next(x for x in item if hasattr(x, "transform"))
        S = st.transform(S)
        s = next(iter(S))
        z = np.asarray(s.z)
        assert np.isfinite(z).any(), f"{st.spec.code} wiped the impedance"
        if st.spec.code == "FREQ004":
            assert np.unique(np.round(np.abs(z[:, 0, 1]), 9)).size > 10
        assert np.isfinite(np.asarray(s.z_err)).any(), st.spec.code
