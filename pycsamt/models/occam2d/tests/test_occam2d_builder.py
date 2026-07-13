"""Tests for InputBuilder."""


def test_builder_imports():
    from pycsamt.models.occam2d.builder import InputBuilder

    assert InputBuilder is not None


def test_builder_not_ready_before_build(tmp_path):
    from pycsamt.models.occam2d.builder import InputBuilder

    builder = InputBuilder(source=None, workdir=tmp_path)
    assert not builder.is_ready


def test_builder_summary_not_built(tmp_path):
    from pycsamt.models.occam2d.builder import InputBuilder

    builder = InputBuilder(source=None, workdir=tmp_path)
    assert "not yet built" in builder.summary()


def test_builder_accepts_edi_style_z_container(tmp_path):
    import numpy as np

    from pycsamt.models.occam2d import InputBuilder, OccamConfig

    class ZContainer:
        freq = np.array([10.0, 1.0])
        resistivity = np.ones((2, 2, 2), dtype=float) * 100.0
        resistivity_err = np.ones((2, 2, 2), dtype=float) * 5.0
        phase = np.ones((2, 2, 2), dtype=float) * 45.0
        phase_err = np.ones((2, 2, 2), dtype=float)

    class EDIStyleSite:
        def __init__(self, name, lon):
            self.station = name
            self.coords = (0.0, lon, 0.0)
            self.Z = ZContainer()

    sites = [EDIStyleSite("S001", 0.0), EDIStyleSite("S002", 0.01)]
    cfg = OccamConfig(modes=["TE", "TM"], freq_min=1.0, freq_max=10.0)

    builder = InputBuilder(sites, workdir=tmp_path, config=cfg).build()

    assert builder.is_ready
    assert builder.data.n_sites == 2
    assert builder.data.n_frequencies == 2
    assert builder.data.n_data == 16
