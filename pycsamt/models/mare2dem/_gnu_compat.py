# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Source fixes that let MARE2DEM build and run with the GNU toolchain.

Upstream MARE2DEM targets Intel compilers only.  Building it with
gfortran + OpenMPI + the free oneMKL was verified end to end on
2026-09-25 (gfortran 15.2, OpenMPI 5.0, MKL; half-space benchmark in
``pycsamt/forward/tests/test_maxwell_mare2dem.py`` passes) after two
kinds of change, both applied here to the *downloaded* source tree
(MARE2DEM is not vendored), idempotently:

**Memory-safety bugs** (applied for every compiler; the replacements are
standard-conforming, so Intel builds are unaffected):

* ``mt1d.f90`` evaluated ``size(this%a)`` of a possibly unallocated
  array — Fortran's ``.or.`` does not short-circuit — and re-allocated
  an allocated array of the wrong size without deallocating it.
* ``readData`` leaves the CSEM correction arrays unallocated for MT-only
  data (and the MT static-shift array for CSEM-only data), yet
  ``em2d.f90`` calls ``any()`` on them unguarded.  They are now
  allocated with size 0 when absent.

Intel builds happened to survive both; gfortran binaries segfaulted.

**Intel-only syntax** (applied for GNU compilers):

* uppercase preprocessor directives and Fortran logical operators inside
  ``#if`` (Intel ``-fpp`` accepts them, GNU cpp does not);
* rank-1 constructors initialising rank-2 arrays -> ``reshape``;
* the ``[a:b]`` array-constructor range shorthand -> an implied-do
  (with a declared loop variable: gfortran has no typed implied-do).

Compiler flags for the remaining differences (``-fallow-argument-mismatch``,
``-fdec-format-defaults``, ``-std=gnu89``) are set by the include-file
generator in :mod:`pycsamt.models.mare2dem.source`.
"""

from __future__ import annotations

import re
from pathlib import Path

LOOPVAR = "j_pyc"
_MARKER = "pycsamt:"  # all inserted code carries this, for traceability

# ── Intel-only syntax ─────────────────────────────────────────────────────


def fix_cpp_directives(text: str) -> str:
    out = []
    for line in text.splitlines(keepends=True):
        if line.lstrip().startswith("#"):
            line = re.sub(r"^(\s*#\s*)([A-Za-z]+)",
                          lambda m: m.group(1) + m.group(2).lower(), line)
            line = re.sub(r"\bDEFINED\b", "defined", line)
            line = re.sub(r"\.not\.", "!", line, flags=re.I)
            line = re.sub(r"\.and\.", "&&", line, flags=re.I)
            line = re.sub(r"\.or\.", "||", line, flags=re.I)
        out.append(line)
    return "".join(out)


def fix_rank2_init(text: str) -> str:
    return re.sub(
        r"(dimension\((\d+)\s*,\s*(\d+)\)\s*::\s*\w+\s*=\s*)\[([^\]\n]*)\]",
        r"\1reshape([\4], [\2,\3])",
        text,
    )


_RANGE = re.compile(r"\[\s*([-+]?\w+)\s*:\s*([^\]\n,]+?)\s*\]")


def fix_range_constructors(text: str) -> str:
    lines = text.splitlines(keepends=True)
    declare_after: set[int] = set()
    for i, line in enumerate(lines):
        code = line.split("!", 1)[0]
        if not _RANGE.search(code):
            continue
        lines[i] = _RANGE.sub(
            lambda m: f"[({LOOPVAR}, {LOOPVAR}={m.group(1)},{m.group(2)})]",
            code,
        ) + line[len(code):]
        for j in range(i - 1, -1, -1):
            if re.match(r"\s*implicit\s+none\b", lines[j], re.I):
                declare_after.add(j)
                break
    for j in sorted(declare_after, reverse=True):
        nxt = lines[j + 1] if j + 1 < len(lines) else ""
        if LOOPVAR not in nxt:
            indent = re.match(r"\s*", lines[j]).group(0)
            lines.insert(
                j + 1, f"{indent}integer :: {LOOPVAR}  ! {_MARKER} gfortran\n"
            )
    return "".join(lines)


# ── Memory-safety bugs ────────────────────────────────────────────────────

_MT1D_ALLOC = re.compile(
    r"^(\s*)if\s*\(\s*\(\s*\.not\.\s*allocated\(this%a\)\s*\)\s*\.or\.\s*"
    r"\(\s*size\(this%a\)\s*/=\s*this%nlayer\s*\)\s*\)\s*"
    r"allocate\(\s*this%a\(this%nlayer\)\s*,\s*this%b\(this%nlayer\)\s*\)"
    r"[^\n]*$",
    re.I | re.M,
)


def fix_mt1d_alloc(text: str) -> str:
    def repl(m):
        i = m.group(1)
        return (
            f"{i}if (allocated(this%a)) then  ! {_MARKER} conforming form\n"
            f"{i}    if (size(this%a) /= this%nlayer) "
            "deallocate(this%a, this%b)\n"
            f"{i}endif\n"
            f"{i}if (.not. allocated(this%a)) "
            "allocate(this%a(this%nlayer), this%b(this%nlayer))"
        )

    return _MT1D_ALLOC.sub(repl, text)


_READDATA_END = re.compile(r"^(end\s+subroutine\s+readData)\s*$",
                           re.I | re.M)
_ZERO_ALLOC = (
    f"    ! {_MARKER} allocate absent correction arrays with size 0 so the\n"
    "    ! unguarded any(...) tests in em2d.f90 are defined.\n"
    "    if (.not. allocated(iEstimateTxCorrection)) "
    "allocate(iEstimateTxCorrection(0))\n"
    "    if (.not. allocated(iEstimateRxCorrection)) "
    "allocate(iEstimateRxCorrection(0))\n"
    "    if (.not. allocated(iEstimateMTStatic)) "
    "allocate(iEstimateMTStatic(0))\n\n"
)


def fix_unallocated_estimates(text: str) -> str:
    if f"{_MARKER} allocate absent correction arrays" in text:
        return text
    return _READDATA_END.sub(lambda m: _ZERO_ALLOC + m.group(1), text,
                             count=1)


# ── Driver ────────────────────────────────────────────────────────────────


def patch_source_tree(src: str | Path, *, gnu: bool) -> list[str]:
    """Patch the MARE2DEM Fortran sources in *src*; return changed files.

    Memory-safety fixes are always applied; Intel-only syntax is rewritten
    only when ``gnu`` is true.  Idempotent.
    """
    changed = []
    for f in sorted(Path(src).glob("*.f90")):
        old = f.read_text(errors="replace")
        new = fix_mt1d_alloc(old)
        if f.name == "mare2dem_io.f90":
            new = fix_unallocated_estimates(new)
        if gnu:
            new = fix_range_constructors(fix_rank2_init(
                fix_cpp_directives(new)))
        if new != old:
            f.write_text(new)
            changed.append(f.name)
    return changed


__all__ = ["patch_source_tree"]


if __name__ == "__main__":  # used by generated build scripts (stdlib only)
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("source_dir")
    ap.add_argument("--intel", action="store_true",
                    help="Intel compilers: apply only the memory-safety fixes")
    ns = ap.parse_args()
    changed = patch_source_tree(ns.source_dir, gnu=not ns.intel)
    print("pycsamt source fixes applied to:",
          ", ".join(changed) if changed else "(none needed)")
