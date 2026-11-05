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


def _phase_sites(xy_phase, yx_phase):
    import numpy as np

    class ZContainer:
        freq = np.array([10.0, 1.0])
        resistivity = np.ones((2, 2, 2), dtype=float) * 100.0
        resistivity_err = np.ones((2, 2, 2), dtype=float) * 5.0
        phase = np.ones((2, 2, 2), dtype=float) * 45.0
        phase_err = np.ones((2, 2, 2), dtype=float)

    z = ZContainer()
    z.phase = z.phase.copy()
    z.phase[:, 0, 1] = xy_phase
    z.phase[:, 1, 0] = yx_phase

    class EDIStyleSite:
        def __init__(self, name, lon):
            self.station = name
            self.coords = (0.0, lon, 0.0)
            self.Z = z

    return [EDIStyleSite("S001", 0.0), EDIStyleSite("S002", 0.01)]


def test_tm_phase_is_wrapped_not_shifted_past_360(tmp_path):
    # A yx phase of +178 deg (wrapped acquisition convention) used to be
    # written as 358 deg; it must be -2 deg (out of quadrant, honestly).
    import numpy as np

    from pycsamt.models.occam2d import InputBuilder, OccamConfig

    cfg = OccamConfig(modes=["TE", "TM"], freq_min=1.0, freq_max=10.0)
    data = InputBuilder(
        _phase_sites(40.0, 178.0), workdir=tmp_path, config=cfg
    ).build().data
    tm_phase = data.data_blocks[data.data_blocks[:, 2] == 6, 3]
    assert np.allclose(tm_phase, -2.0)
    assert data.masked == []


def test_phase_quadrant_drops_rows_and_their_rho_partner(tmp_path):
    import numpy as np

    from pycsamt.models.occam2d import InputBuilder, OccamConfig

    cfg = OccamConfig(
        modes=["TE", "TM"], freq_min=1.0, freq_max=10.0,
        phase_quadrant=(0.0, 90.0),
    )
    # TE (xy) phase -12 deg is out of quadrant; TM (-135+180=45) is kept.
    data = InputBuilder(
        _phase_sites(-12.0, -135.0), workdir=tmp_path, config=cfg
    ).build().data
    types = data.data_blocks[:, 2]
    assert not np.any(np.isin(types, [1, 2]))  # TE rho+phase both dropped
    assert np.count_nonzero(types == 5) == np.count_nonzero(types == 6) == 4
    assert len(data.masked) == 4
    assert {m["mode"] for m in data.masked} == {"TE"}
