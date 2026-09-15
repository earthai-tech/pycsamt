"""Focused branch coverage for the lightweight MARE2DEM I/O helpers."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.models.mare2dem.iotools.emdata import (
    CSEMConfig,
    DCConfig,
    EMDataFile,
    MTConfig,
    UTMOrigin,
    _EMDataReader,
    _parse_line,
    _strtok as em_strtok,
    read_emdata,
    write_emdata,
)
from pycsamt.models.mare2dem.iotools.codes import (
    code_component,
    code_label,
    code_representation,
    is_csem_code,
    is_mt_code,
)
from pycsamt.models.mare2dem.iotools.data_group import (
    DataGroupFile,
    read_data_group,
    write_data_group,
)
from pycsamt.models.mare2dem.iotools.group_rms import read_group_rms_log
from pycsamt.models.mare2dem.iotools.poly import (
    PolyFile,
    _next_data_line,
    _read_float,
    _read_int,
    _triangle_neighbors,
    read_poly,
    read_triangulation,
    write_poly,
    write_triangulation,
)
from pycsamt.models.mare2dem.iotools.resistivity import (
    ResistivityFile,
    _nrho,
    _strtok as rho_strtok,
    read_resistivity,
    write_resistivity,
)


def test_data_code_helpers_cover_known_and_unknown_codes():
    assert code_component(104) == "Zxy (TE)"
    assert code_representation(104) == "Phase"
    assert "Zxy (TE)" in code_label(104)
    assert code_label(999) == "Unknown (code=999)"
    assert code_component(999) == ""
    assert code_representation(999) == ""
    assert is_mt_code(104) and not is_mt_code(1)
    assert is_csem_code(1) and is_csem_code(151) and not is_csem_code(104)


def test_data_group_edge_cases(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_data_group(tmp_path / "missing")

    bad = tmp_path / "bad.group"
    bad.write_text("Format: EMDataGroup_9.9\n")
    with pytest.raises(ValueError, match="Unsupported"):
        read_data_group(bad)

    source = tmp_path / "comments.group"
    source.write_text(
        "% ignored\n! hello\nnot a field\nVersion: EMDataGroup_1.0 % ok\n"
        "# groups: 2 ! inline\n% skip\nA\n! skip\nB\n#data: 3\n1 2 1 2\n"
    )
    dg = read_data_group(source)
    assert dg.comment.strip() == "hello"
    assert dg.n_groups == 2 and dg.n_data == 3
    assert repr(dg) == "DataGroupFile(n_groups=2, n_data=3)"

    with pytest.raises(ValueError, match="must not be empty"):
        write_data_group(DataGroupFile(), tmp_path / "empty")
    with pytest.raises(ValueError, match=">= 1"):
        write_data_group(
            DataGroupFile(group_names=["A"], group_indices=np.array([0])),
            tmp_path / "zero",
        )
    written = write_data_group(
        DataGroupFile(comment=" note", group_names=["A"]), tmp_path / "nested/x"
    )
    assert read_data_group(written).n_data == 0


def test_group_rms_missing_and_empty_inputs(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_group_rms_log(tmp_path / "missing.log")
    empty = tmp_path / "empty.log"
    empty.write_text("")
    log = read_group_rms_log(empty)
    assert log.n_iterations == 0
    assert "n_iterations=0" in repr(log)


def test_poly_rich_roundtrip_and_helpers(tmp_path):
    assert _read_int([], 0, 7) == 7
    assert _read_int(["x"], 0, 8) == 8
    assert _read_float([], 0, 1.5) == 1.5
    assert _read_float(["x"], 0, 2.5) == 2.5
    assert _next_data_line(["", "# comment"], 0) == ([], 2)

    pf = PolyFile()
    pf.nodes = np.array([[0, 0], [1, 0], [0, 1]], dtype=float)
    pf.node_boundary_markers = np.array([1, 2, 3])
    pf.segments = np.array([[1, 2], [2, 3], [3, 1]])
    pf.segment_markers = np.array([4, 5, 6])
    pf.holes = np.array([[0.2, 0.2]])
    pf.regions = np.array([[0.3, 0.3, 7, 0.01]])
    assert (pf.n_nodes, pf.n_segments, pf.n_holes, pf.n_regions) == (3, 3, 1, 1)
    assert "n_regions=1" in repr(pf)
    out = write_poly(pf, tmp_path / "nested/rich.poly")
    got = read_poly(out)
    np.testing.assert_array_equal(got.node_boundary_markers, [1, 2, 3])
    np.testing.assert_array_equal(got.segment_markers, [4, 5, 6])
    np.testing.assert_allclose(got.holes, pf.holes)
    np.testing.assert_allclose(got.regions, pf.regions)
    with pytest.raises(FileNotFoundError):
        read_poly(tmp_path / "absent.poly")


def test_triangulation_roundtrip_defaults_and_validation(tmp_path):
    nodes = np.array([[0, 0], [1, 0], [0, 1], [1, 1]], dtype=float)
    triangles = np.array([[0, 1, 2], [1, 3, 2]])
    path = write_triangulation(nodes, triangles, None, tmp_path / "mesh.node")
    got_nodes, got_triangles, attrs = read_triangulation(path)
    np.testing.assert_allclose(got_nodes, nodes)
    np.testing.assert_array_equal(got_triangles, triangles)
    np.testing.assert_array_equal(attrs, [1, 1])
    neighbors = _triangle_neighbors(triangles)
    assert 1 in neighbors[0] and 0 in neighbors[1]

    with pytest.raises(FileNotFoundError, match="node file"):
        read_triangulation(tmp_path / "none.node")
    lonely = tmp_path / "lonely.node"
    lonely.write_text("0 2 0 0\n")
    with pytest.raises(FileNotFoundError, match="companion"):
        read_triangulation(lonely)
    non_manifold = np.array([[0, 1, 2], [0, 1, 3], [0, 1, 4]])
    with pytest.raises(ValueError, match="non-manifold"):
        _triangle_neighbors(non_manifold)


def test_read_poly_companion_node_and_ele(tmp_path):
    (tmp_path / "done.poly").write_text("0 2 0 0\n0 0\n0\n0\n")
    (tmp_path / "done.node").write_text(
        "3 2 1 1\n1 0 0 10 1\n2 1 0 20 2\n3 0 1 30 3\n"
    )
    (tmp_path / "done.ele").write_text("1 3 1\n1 1 2 3 9\n")
    pf = read_poly(tmp_path / "done.poly")
    assert pf.nodes.shape == (3, 2)
    assert pf.node_attributes.shape == (3, 1)
    np.testing.assert_array_equal(pf.node_boundary_markers, [1, 2, 3])
    np.testing.assert_array_equal(pf.segments, [[1, 2, 3]])


def test_full_emdata_writer_and_reader(tmp_path):
    em = EMDataFile(
        comment=" comprehensive",
        utm=UTMOrigin(grid=30, hemi="N", north0=1, east0=2, theta=3),
        mt=MTConfig(
            frequencies=np.array([1.0]),
            receivers=np.array([[1, 2, 3, 4, 5, 6, 7, 1]], dtype=float),
            receiver_name=["MT1"],
        ),
        csem=CSEMConfig(
            phase_convention="lead",
            reciprocity_used="yes",
            frequencies=np.array([2.0]),
            time_offsets=np.array([0.1, 0.2]),
            transmitters=np.array([[1, 2, 3, 4, 5, 6]], dtype=float),
            transmitter_type=["edipole"],
            transmitter_name=["TX1"],
            receivers=np.array([[1, 2, 3, 4, 5, 6]], dtype=float),
            receiver_name=["RX1"],
        ),
        dc=DCConfig(
            tx_electrodes=np.array([[0, 0, 0], [1, 0, 0]], dtype=float),
            rx_electrodes=np.array([[0, 1, 0], [1, 1, 0]], dtype=float),
            transmitters=np.array([[1, 2]]),
            receivers=np.array([[1, 2]]),
            # The legacy reader consumes only the two integer columns.
            transmitter_name=[],
            receiver_name=[],
        ),
        data=np.array([[1, 1, 1, 1, 2.5, 0.1]], dtype=float),
    )
    path = write_emdata(em, tmp_path / "all/full.emdata")
    got = read_emdata(path)
    assert got.n_data == 1
    assert got.n_mt_frequencies == got.n_mt_receivers == 1
    assert got.n_csem_transmitters == got.n_csem_receivers == 1
    assert got.dc is not None and len(got.dc.transmitters) == 1
    assert "n_data=1" in repr(got)

    em.is_response = True
    em.data = np.array([[1, 1, 1, 1, 2.5, 0.1, 2.4, 1.0]])
    response = read_emdata(write_emdata(em, tmp_path / "full.resp"))
    assert response.is_response and response.data.shape == (1, 8)


def test_emdata_helpers_and_errors(tmp_path):
    assert _parse_line("1D+02 % tail") == ("1e+02", False)
    assert _parse_line("! comment")[1]
    assert em_strtok("plain") == ("plain", "")
    reader = _EMDataReader(["! skip"])
    assert reader._next_data_line() is None
    assert reader._read_float_block(1).size == 0
    with pytest.raises(ValueError, match="reshape"):
        reader._read_int_block(1, 2)
    with pytest.raises(FileNotFoundError):
        read_emdata(tmp_path / "missing")


@pytest.mark.parametrize(
    ("version", "tx", "csem_rx", "mt_rx"),
    [
        ("2.0", "1 2 3 4 5 edipole TX", "1 2 3 4 5 6", "1 2 3 4 5 6"),
        ("2.1", "1 2 3 4 5 edipole TX", "1 2 3 4 5 6 RX", "1 2 3 4 5 6 7 MT"),
        ("2.2", "1 2 3 4 5 6 edipole TX", "1 2 3 4 5 6 7 RX", "1 2 3 4 5 6 7 1 MT"),
    ],
)
def test_emdata_legacy_layouts(tmp_path, version, tx, csem_rx, mt_rx):
    path = tmp_path / f"legacy-{version}.emdata"
    path.write_text(
        f"Format: EMData_{version}\n"
        "# Transmitters: 1\n" + tx + "\n"
        "# CSEM Receivers: 1\n" + csem_rx + "\n"
        "# MT Receivers: 1\n" + mt_rx + "\n"
        "# Data: 1\n! header\n1 1 1 1 2 0.1\n"
    )
    got = read_emdata(path)
    assert got.n_csem_transmitters == got.n_csem_receivers == 1
    assert got.n_mt_receivers == 1


def test_emdata_slow_data_path_and_tdem_writer(tmp_path):
    slow = tmp_path / "slow.emdata"
    slow.write_text(
        "Format: EMData_2.3\n# Data: 2\n"
        "1 1 1 1 2 0.1\n! interspersed\n1 1 1 1 3 0.2\n"
    )
    assert read_emdata(slow).data.shape == (2, 6)

    em = EMDataFile(
        csem=CSEMConfig(
            transmitters=np.array([[0, 0, 0, 0, 0, 1, 0]], dtype=float),
            tdem_waveform=np.array([[0.0, 1.0]]),
        )
    )
    text = write_emdata(em, tmp_path / "tdem.emdata").read_text()
    assert "TDEM waveform" in text


def test_resistivity_full_anisotropic_roundtrip(tmp_path):
    assert _nrho("unknown") == 1
    assert rho_strtok("Key") == ("key", "")
    rf = ResistivityFile(
        resistivity_file=str(tmp_path / "fallback"),
        poly_file="mesh.poly",
        data_file="survey.emdata",
        settings_file="run.settings",
        anisotropy="triaxial",
        data_group_file="Groups.File",
        joint_inv_weight_type="equal",
        penalty_file="penalty.dat",
        roughness_with_prejudice=True,
        anisotropy_penalty_weight=2.5,
        anisotropy_ratio_roughness_weight=3.5,
        fixed_mu_cut=0.2,
        roughness=12.0,
        misfit=1.1,
        date_time="today",
        resistivity=np.array([[10, 20, 30], [40, 50, 60]], dtype=float),
        free_parameter=np.array([[1, 0, 1], [1, 1, 0]], dtype=float),
        bounds=np.tile([1, 100], (2, 3)),
        prejudice=np.tile([10, 0.5], (2, 3)),
    )
    path = write_resistivity(rf, tmp_path / "nested/model")
    assert path.suffix == ".resistivity"
    got = read_resistivity(path)
    assert got.num_regions == 2
    assert got.data_group_file == "groups.file"
    assert got.joint_inv_weight_type == "equal"
    assert got.anisotropy_penalty_weight == 2.5
    assert got.anisotropy_ratio_roughness_weight == 3.5
    assert got.roughness_with_prejudice
    assert "num_regions=2" in repr(got)
    np.testing.assert_allclose(got.resistivity, rf.resistivity)

    rf.free_parameter = np.empty((0, 3))
    rf.bounds = None
    rf.prejudice = None
    assert write_resistivity(rf).exists()
    with pytest.raises(FileNotFoundError):
        read_resistivity(tmp_path / "absent")


def test_resistivity_invalid_optional_values_are_ignored(tmp_path):
    fields = [
        "Maximum Iterations", "Penalty Cut Weight", "Min. Gradient Support Weight",
        "Aniso. Penalty Weight", "Aniso. Ratio Roughness Weight", "Debug Level",
        "Target Misfit", "Iteration", "Lagrange Value", "Fixed Mu Cut",
        "Misfit Decrease Threshold",
    ]
    text = "\n% comment\n" + "\n".join(f"{field}: nope" for field in fields)
    path = tmp_path / "invalid.resistivity"
    path.write_text(text)
    got = read_resistivity(path)
    assert got.max_iterations == 100 and got.target_misfit == 1.0
