"""Tests for pycsamt.models.mare2dem._gnu_compat (GNU source fixes).

Each snippet is the exact upstream construct gfortran/GNU cpp rejected (or
mis-executed) in the 2026-09-25 verification build.
"""

from __future__ import annotations

from pycsamt.models.mare2dem import _gnu_compat as g


def test_cpp_directives_normalised():
    src = ("#IF DEFINED(_WIN32) .OR. DEFINED(_WIN64)\n"
           "#DEFINE a b\n#ENDIF\n"
           "#if .not. (defined (_WIN32) || defined (WIN32))\n"
           "x = .not. y  ! Fortran code untouched\n")
    out = g.fix_cpp_directives(src)
    assert "#if defined(_WIN32) || defined(_WIN64)" in out
    assert "#define a b" in out and "#endif" in out
    assert "#if ! (defined (_WIN32) || defined (WIN32))" in out
    assert "x = .not. y" in out


def test_rank2_initialiser_reshaped():
    src = "    real(8), dimension(3,6) :: V =[1,0,0, 0,1,0, 0,0,1, 1,0,0, 0,1,0, 0,0,1]"
    out = g.fix_rank2_init(src)
    assert "reshape([1,0,0, 0,1,0, 0,0,1, 1,0,0, 0,1,0, 0,0,1], [3,6])" in out


def test_range_constructor_rewritten_with_declared_loop_variable():
    src = ("subroutine s(n)\n    implicit none\n    integer :: n, a(n)\n"
           "    a = [1:n] + 1\n    a(1:2) = [0:1]\nend subroutine\n")
    out = g.fix_range_constructors(src)
    assert "[(j_pyc, j_pyc=1,n)]" in out and "[(j_pyc, j_pyc=0,1)]" in out
    assert out.count("integer :: j_pyc") == 1  # declared once per scope
    assert "a(1:2)" in out  # array sections (parentheses) are untouched
    assert g.fix_range_constructors(out) == out  # idempotent


def test_mt1d_non_short_circuit_allocation_fixed():
    src = ("    if ( (.not.allocated(this%a) ) .or. ( size(this%a) /= "
           "this%nlayer ) )  allocate( this%a(this%nlayer), "
           "this%b(this%nlayer))\n")
    out = g.fix_mt1d_alloc(src)
    assert "size(this%a)" in out and ".or." not in out
    assert "deallocate(this%a, this%b)" in out
    assert "if (.not. allocated(this%a)) allocate(" in out


def test_unallocated_estimate_arrays_get_size_zero():
    src = "    write(*,*) ' '\n\nend subroutine readData\n"
    out = g.fix_unallocated_estimates(src)
    for arr in ("iEstimateTxCorrection", "iEstimateRxCorrection",
                "iEstimateMTStatic"):
        assert f"allocate({arr}(0))" in out
    assert out.index("allocate(") < out.index("end subroutine readData")
    assert g.fix_unallocated_estimates(out) == out


def test_patch_tree_intel_gets_only_memory_fixes(tmp_path):
    (tmp_path / "mt1d.f90").write_text(
        "    if ( (.not.allocated(this%a) ) .or. ( size(this%a) /= "
        "this%nlayer ) )  allocate( this%a(this%nlayer), "
        "this%b(this%nlayer))\n")
    (tmp_path / "call_triangle.f90").write_text("#IF DEFINED(X)\n#ENDIF\n")
    changed = g.patch_source_tree(tmp_path, gnu=False)
    assert changed == ["mt1d.f90"]
    assert "#IF" in (tmp_path / "call_triangle.f90").read_text()
    changed = g.patch_source_tree(tmp_path, gnu=True)
    assert changed == ["call_triangle.f90"]
    assert g.patch_source_tree(tmp_path, gnu=True) == []  # idempotent


def test_generate_inc_gnu_flags_and_versioned_mkl(tmp_path):
    from pycsamt.models.mare2dem import source as s

    mkl = tmp_path / "mkl"
    (mkl / "include").mkdir(parents=True)
    (mkl / "lib").mkdir()
    (mkl / "include" / "mkl_dss.f90").write_text("")
    for lib in ("mkl_core", "mkl_sequential", "mkl_gf_lp64"):
        (mkl / "lib" / f"lib{lib}.so.2").write_bytes(b"")  # pip wheel style
    inc = s._generate_inc("mpifort", "mpicc", str(mkl), tmp_path / "b")
    text = inc.read_text()
    assert "-fallow-argument-mismatch" in text
    assert "-fdec-format-defaults" in text
    assert "TRICOPTS = -O2 -fPIC -std=gnu89" in text
    assert "-l:libmkl_gf_lp64.so.2" in text  # gfortran interface, versioned
    assert "-lmkl_intel_lp64" not in text


def test_runner_launches_wsl_binary_through_wsl(tmp_path):
    """A ``wsl:``-registered MARE2DEM (built in WSL2 by the Solver Builder)
    runs via ``wsl -e bash -lc`` with the work dir mapped to /mnt."""
    from pycsamt.models.mare2dem.config import Mare2DEMConfig
    from pycsamt.models.mare2dem.runner import _wsl_command

    cfg = Mare2DEMConfig(binary="wsl:/home/u/.local/share/pycsamt/mare2dem/"
                                "build/MARE2DEM", n_procs=4)
    wd = tmp_path / "run"
    wd.mkdir()
    cmd = _wsl_command(cfg, wd, "mare2dem.resistivity", True, 4, None)
    assert cmd[:4] == ["wsl", "-e", "bash", "-lc"]
    inner = cmd[4]
    assert "toolchain/mare2dem" in inner  # managed mpirun on PATH
    assert "-np 4" in inner and inner.rstrip().endswith("MARE2DEM mare2dem")
    assert "cd " in inner and "/run" in inner
