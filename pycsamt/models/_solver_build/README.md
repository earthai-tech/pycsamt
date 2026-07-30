# pycsamt solver build scripts

Bash scripts that compile pycsamt's external forward/inversion
solvers from source, so using them does not require hand-writing
`make` invocations or knowing each solver's own build quirks.

| Script | Solver | Source | Toolchain |
|---|---|---|---|
| `modem2d.sh` | ModEM `Mod2DMT` (2-D) | vendored, committed | gfortran + make (+ LAPACK/BLAS) |
| `modem3d.sh` | ModEM `Mod3DMT` (3-D) | vendored, committed | gfortran + make (+ LAPACK/BLAS) |
| `occam2d.sh` | Occam2DMT | vendored, committed | gfortran + make |
| `mare2dem.sh` | MARE2DEM | downloaded on demand (not vendored — see below) | Intel `mpiifort`/`mpiicc` + MKL, Linux/macOS/WSL only |

A single dispatcher is also provided:

```bash
pycsamt/models/_solver_build/build.sh modem3d --auto-install -y
pycsamt/models/_solver_build/build.sh occam2d
pycsamt/models/_solver_build/build.sh mare2dem
```

which is equivalent to calling `modem3d.sh`/`occam2d.sh`/`mare2dem.sh`
directly — use whichever is more convenient.

## Preferred entry points

Once pycsamt is installed, you normally don't need this path at all —
the packaged CLI and the repo-root `Makefile` both forward to these
same scripts with the same arguments:

```bash
pycsamt build modem3d --auto-install -y   # works from any directory
make modem3d ARGS="--auto-install -y"     # repo checkout convenience
```

This directory (and the paths above) stay the single source of truth
for what actually runs; see
`docs/source/user_guide/models/compilation.rst` for the full guide.

## Quick start

```bash
# ModEM 3-D, installing a Fortran toolchain automatically if needed
pycsamt/models/_solver_build/modem3d.sh --auto-install

# ModEM 2-D, assuming the toolchain is already installed
pycsamt/models/_solver_build/modem2d.sh

# Occam2D
pycsamt/models/_solver_build/occam2d.sh --auto-install -y

# MARE2DEM: status only (no download/build without --auto-install)
pycsamt/models/_solver_build/mare2dem.sh
```

On success, each script prints the resulting binary's path and how to
point the matching pycsamt config at it (`ModEmConfig(binary_2d=...)`
/ `binary_3d=...`, `Mare2DEMConfig(binary=...)`). If you do nothing,
pycsamt's own adapters already search each solver's `_source/`
directory as a fallback, so a binary built in place is often found
automatically.

## What `--auto-install` actually does, per platform

Never runs silently: every install command is printed before it runs,
and requires a `y` confirmation unless you also pass `-y`/`--yes`.

- **Linux**: `apt-get`/`dnf`/`yum install gfortran` + LAPACK/BLAS dev
  packages (whichever package manager is found first).
- **macOS**: `brew install gcc openblas` (Homebrew's `gcc` formula
  includes `gfortran`).
- **Windows**: gfortran and `make` are not reliably available outside
  a dedicated MinGW-w64 toolchain, so these scripts create a small,
  **isolated** conda environment (`pycsamt-fortran` by default,
  override with `PYCSAMT_FORTRAN_ENV`) with `m2w64-gcc-fortran`,
  `m2w64-openblas`, and `make` from conda-forge — it does not touch
  your existing pycsamt environment. Resolved directly by absolute
  path from the env's own `Library/mingw-w64/bin`/`Library/bin`, so
  you never need to `conda activate` it yourself.

MARE2DEM is the one exception: it needs Intel's commercial compilers
and MKL, which nothing here can install for you (see `mare2dem.sh`'s
own header for details) — `--auto-install` there only runs the
download-then-build step assuming you already have those.

## Windows-specific notes

- A MinGW-w64-built `.exe` needs its own runtime DLLs either on
  `PATH` or sitting next to it (Windows checks the executable's own
  directory first). These scripts copy the handful of DLLs the
  `pycsamt-fortran` conda env provides (`libgcc_s_seh-1.dll`,
  `libgfortran-3.dll`, `libopenblas.dll`, `libquadmath-0.dll`,
  `libwinpthread-1.dll`) alongside the binary automatically, so it
  runs standalone without activating anything.
- MARE2DEM cannot be built natively on Windows at all (its own
  `SourceManager.build()` raises immediately on
  `sys.platform == "win32"`); use WSL2.

## Why these scripts retry `make` with `-k`

The vendored ModEM Makefiles list their object files in an order that
is not strictly a topological sort of Fortran module dependencies (a
pre-existing property of the adapted-from-upstream Makefiles, not
something introduced here). A single-pass `make` can therefore fail
with "module file not found" purely because of build order, not a
real compile error. `lib/common.sh`'s `run_make_with_retries` re-runs
`make -k` (keep-going past errors) for up to a handful of passes,
since each pass can produce `.mod` files a later file depended on —
confirmed empirically while first compiling this source (the 3-D
build needed 4 passes to converge). `-j` (parallel `make`) is
deliberately never used here for the same reason: it would make the
already-imperfect ordering non-deterministic instead of merely
sub-optimal.

## Known model constraint (ModEM 3-D only, not a build issue)

ModEM's own Fortran (not this project's code) hardcodes 10 air layers
for the 3-D solver and its default "mirror" air-layer sizing method
reads that many earth-layer widths from the model with no bounds
check against the actual earth cell count. A model with fewer than 10
earth z-cells will build and launch fine, then crash mid-solve with
`b in QMR contains NaNs`. `pycsamt.forward.maxwell.modem3d.ModEm3DAdapter`
already guards against this (`_MIN_EARTH_Z_CELLS`), but it is worth
knowing if you ever drive `Mod3DMT` directly.

## Bugs these scripts' Makefile fixes correct

Building ModEM for the first time in this project surfaced (and these
now-committed Makefiles fix) two real, platform-independent gaps that
would reproduce identically on Linux/macOS, not just Windows:

1. A missing build rule for `FIELDS/FiniteDiff3D/sg_spherical.f90`
   (3-D only — `use`d directly by `GridCalc.f90` but never in the
   original object list at all).
2. `MPI/Declaration_MPI.f90`, `{2D,3D}_MT/Sub_MPI.f90`, and
   `MPI/Main_MPI.f90` were never compiled, despite being
   unconditionally `use`d by several files regardless of whether MPI
   is enabled (their entire module bodies are wrapped in
   `#ifdef MPI` / `#endif`) — the serial build also never passed
   `-cpp` to strip those guards.

See `pycsamt/forward/maxwell/modem3d.py`'s module docstring and the
`AI-INVERSION-M6-3D-ADR.md`'s 2026-07-29 update for the full account,
including further bugs found in pycsamt's own adapter code (not the
Makefiles) while validating a real compiled binary end-to-end.
