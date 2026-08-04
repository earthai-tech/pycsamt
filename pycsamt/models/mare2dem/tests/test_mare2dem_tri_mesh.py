# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for MARE2DEM triangular mesh generation.

Covers :mod:`pycsamt.models.mare2dem.tri_mesh`,
:mod:`pycsamt.models.mare2dem.triangle_exec`,
:meth:`~pycsamt.models.mare2dem.source.SourceManager.resolve_triangle_binary`,
and :func:`~pycsamt.models.mare2dem.iotools.poly.read_triangulation`.

No real Triangle executable is required: the end-to-end path is tested
against a scripted stand-in (mirrors
``pycsamt/forward/tests/test_maxwell_modem3d.py``'s technique for a
real solver binary that isn't installed in this environment either).
"""

from __future__ import annotations

import os
import sys
import textwrap
from pathlib import Path

import numpy as np
import pytest

from pycsamt.models.mare2dem import (
    Mare2DEMConfig,
    MTSurveyConfig,
    build_pslg,
    build_survey_mesh,
    make_data_file,
    read_triangulation,
    survey_points,
    tri_mesh_from_poly,
)
from pycsamt.models.mare2dem.source import SourceManager
from pycsamt.models.mare2dem.triangle_exec import TriangleRunError, run_triangle

# ---------------------------------------------------------------------------
# build_pslg
# ---------------------------------------------------------------------------


def test_build_pslg_closes_the_boundary_polygon():
    box = [[0, 0], [1000, 0], [1000, 500], [0, 500]]
    pf = build_pslg(box)
    assert pf.n_nodes == 4
    assert pf.n_segments == 4
    # closing segment goes from the last node back to the first
    assert pf.segments[-1].tolist() == [4, 1]


def test_build_pslg_with_interior_points_and_region_seeds():
    box = [[0, 0], [1000, 0], [1000, 500], [0, 500]]
    interior = [[200, 0], [800, 0]]
    seeds = [[500, 250]]
    pf = build_pslg(box, interior_points=interior, region_seeds=seeds)
    assert pf.n_nodes == 6
    assert pf.n_segments == 4  # interior points add no segments
    assert pf.n_regions == 1
    np.testing.assert_allclose(pf.regions[0, :2], [500, 250])
    assert pf.regions[0, 2] == 1  # region attribute


def test_build_pslg_rejects_degenerate_boundary():
    with pytest.raises(ValueError, match="shape"):
        build_pslg([[0, 0], [1, 0]])


# ---------------------------------------------------------------------------
# SourceManager.resolve_triangle_binary
# ---------------------------------------------------------------------------


def _make_fake_binary(directory: Path, name: str) -> Path:
    exe_name = f"{name}.exe" if os.name == "nt" else name
    exe = directory / exe_name
    exe.write_text("not a real binary\n")
    if os.name != "nt":
        exe.chmod(0o755)
    return exe


def test_resolve_triangle_binary_finds_source_dir_build_product(tmp_path):
    src = tmp_path / "src"
    src.mkdir()
    exe = _make_fake_binary(src, "triangle")

    mgr = SourceManager(config=Mare2DEMConfig(source_dir=str(src)))
    found = mgr.resolve_triangle_binary()

    assert found is not None
    assert found.resolve() == exe.resolve()


def test_resolve_triangle_binary_returns_none_when_missing(tmp_path):
    mgr = SourceManager(
        config=Mare2DEMConfig(source_dir=str(tmp_path / "empty"))
    )
    assert mgr.resolve_triangle_binary() is None


# ---------------------------------------------------------------------------
# run_triangle: scripted stand-in (no real Triangle binary needed)
# ---------------------------------------------------------------------------

_FAKE_TRIANGLE_SCRIPT = textwrap.dedent(
    """
    import sys
    from pathlib import Path

    import numpy as np
    from scipy.spatial import Delaunay

    # The real subprocess argv is [flags, *extra_args, poly_name] --
    # poly_path is always last regardless of how many extra flags precede
    # it, so index from the end rather than assuming a fixed position.
    poly_path = Path(sys.argv[-1])
    argv_flags = sys.argv[1:-1]
    stem = poly_path.stem

    lines = [
        ln for ln in poly_path.read_text().splitlines() if ln.strip()
    ]
    n_nodes = int(lines[0].split()[0])
    coords = []
    for line in lines[1 : 1 + n_nodes]:
        parts = line.split()
        coords.append((float(parts[1]), float(parts[2])))

    if "-fail" in argv_flags:
        sys.stderr.write("simulated Triangle failure\\n")
        sys.exit(1)

    node_path = poly_path.with_name(f"{stem}.1.node")
    ele_path = poly_path.with_name(f"{stem}.1.ele")
    poly_out = poly_path.with_name(f"{stem}.1.poly")

    # A real constrained-Delaunay quality mesher; not real Triangle, but a
    # genuine (non-fabricated) triangulation of the input points -- robust
    # to the collinear boundary points a flat-topography survey produces,
    # unlike a naive single-hub fan triangulation.
    pts = np.array(coords)
    tris = Delaunay(pts).simplices

    with node_path.open("w") as fh:
        fh.write(f"{len(coords)} 2 0 0\\n")
        for i, (x, z) in enumerate(coords):
            fh.write(f"{i} {x} {z}\\n")

    with ele_path.open("w") as fh:
        fh.write(f"{len(tris)} 3 1\\n")
        for i, tri in enumerate(tris):
            fh.write(f"{i} {tri[0]} {tri[1]} {tri[2]} 1\\n")

    if "-noout" not in argv_flags:
        with poly_out.open("w") as fh:
            fh.write("0 2 0 0\\n0 0\\n0\\n0\\n")
    """
)


def _write_fake_triangle(tmp_path: Path) -> Path:
    script = tmp_path / "fake_triangle.py"
    script.write_text(_FAKE_TRIANGLE_SCRIPT)
    return script


def test_run_triangle_end_to_end_against_scripted_stand_in(tmp_path):
    script = _write_fake_triangle(tmp_path)
    box = [[0, 0], [1000, 0], [1000, 500], [0, 500]]
    from pycsamt.models.mare2dem.iotools.poly import write_poly

    poly_path = write_poly(build_pslg(box), tmp_path / "mesh.poly")

    refined = run_triangle(
        poly_path, executable=[sys.executable, str(script)]
    )

    assert refined.name == "mesh.1.poly"
    nodes, triangles, region_attrs = read_triangulation(
        refined.with_suffix(".node")
    )
    assert nodes.shape == (4, 2)
    assert triangles.shape == (2, 3)
    assert region_attrs.tolist() == [1, 1]
    assert triangles.min() == 0  # normalized 0-based


def test_run_triangle_raises_on_nonzero_exit(tmp_path):
    script = _write_fake_triangle(tmp_path)
    from pycsamt.models.mare2dem.iotools.poly import write_poly

    poly_path = write_poly(
        build_pslg([[0, 0], [1, 0], [1, 1], [0, 1]]), tmp_path / "mesh.poly"
    )

    with pytest.raises(TriangleRunError, match="exited"):
        run_triangle(
            poly_path,
            executable=[sys.executable, str(script)],
            extra_args=("-fail",),
        )


def test_run_triangle_raises_when_output_missing(tmp_path):
    script = _write_fake_triangle(tmp_path)
    from pycsamt.models.mare2dem.iotools.poly import write_poly

    poly_path = write_poly(
        build_pslg([[0, 0], [1, 0], [1, 1], [0, 1]]), tmp_path / "mesh.poly"
    )

    with pytest.raises(TriangleRunError, match="was not produced"):
        run_triangle(
            poly_path,
            executable=[sys.executable, str(script)],
            extra_args=("-noout",),
        )


def test_run_triangle_raises_when_no_executable_resolvable(tmp_path):
    from pycsamt.models.mare2dem.iotools.poly import write_poly

    poly_path = write_poly(
        build_pslg([[0, 0], [1, 0], [1, 1], [0, 1]]), tmp_path / "mesh.poly"
    )
    empty_source = tmp_path / "no_triangle_here"

    with pytest.raises(TriangleRunError, match="no Triangle executable"):
        run_triangle(
            poly_path,
            config=Mare2DEMConfig(source_dir=str(empty_source)),
        )


def test_run_triangle_raises_on_missing_input(tmp_path):
    with pytest.raises(FileNotFoundError):
        run_triangle(tmp_path / "missing.poly")


# ---------------------------------------------------------------------------
# tri_mesh_from_poly / build_survey_mesh
# ---------------------------------------------------------------------------


def test_tri_mesh_from_poly_builds_valid_trimesh(tmp_path):
    script = _write_fake_triangle(tmp_path)
    from pycsamt.models.mare2dem.iotools.poly import write_poly

    box = [[0, 0], [1000, 0], [1000, 500], [0, 500]]
    poly_path = write_poly(build_pslg(box), tmp_path / "mesh.poly")
    refined = run_triangle(
        poly_path, executable=[sys.executable, str(script)]
    )

    mesh = tri_mesh_from_poly(refined)

    assert mesh.n_nodes == 4
    assert mesh.n_triangles == 2
    assert mesh.dimension == 2
    assert mesh.region_ids.tolist() == [1, 1]


def _synthetic_mt_survey():
    mt = MTSurveyConfig(
        frequencies=np.logspace(-2, 2, 5),
        rx_y=np.linspace(-1000, 1000, 5),
        rx_type="land",
        lTE=True,
        lTM=True,
    )
    return mt


def test_survey_points_matches_area_of_interest_geometry(tmp_path):
    mt = _synthetic_mt_survey()
    em = make_data_file(tmp_path / "survey.emdata", topo=0.0, mt=mt)

    points = survey_points(em)

    assert points.shape[1] == 2
    np.testing.assert_allclose(np.sort(points[:, 0]), np.sort(mt.rx_y))


def test_build_survey_mesh_end_to_end_against_scripted_stand_in(tmp_path):
    script = _write_fake_triangle(tmp_path)
    mt = _synthetic_mt_survey()
    em = make_data_file(tmp_path / "survey.emdata", topo=0.0, mt=mt)

    mesh = build_survey_mesh(
        em,
        topo=0.0,
        simplify_tolerance_m=1.0,
        workdir=tmp_path / "run",
        executable=[sys.executable, str(script)],
    )

    assert mesh.n_nodes >= 4  # at least the boundary box
    assert mesh.n_triangles > 0
    # every receiver y-position must appear as a mesh node (required PSLG
    # interior point), since the fake script echoes input nodes back
    receiver_y = np.sort(mt.rx_y)
    mesh_y = np.sort(mesh.nodes_m[:, 0])
    for y in receiver_y:
        assert np.any(np.isclose(mesh_y, y))


def test_build_survey_mesh_rejects_survey_without_geometry(tmp_path):
    from pycsamt.models.mare2dem.iotools.emdata import EMDataFile

    empty_em = EMDataFile()
    with pytest.raises(ValueError, match="geometry"):
        build_survey_mesh(empty_em, workdir=tmp_path)
