#!/usr/bin/env bash
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
#
# Download and build MARE2DEM (Key, 2016) via pycsamt's own
# pycsamt.models.mare2dem.SourceManager.
#
# Unlike the ModEM/Occam2D scripts in this directory, MARE2DEM's
# source is NOT vendored in this repository (it is ~49 MB and under
# its own license; see pycsamt/models/mare2dem/_source/.gitkeep) and
# its build genuinely needs Intel compilers (mpiifx/mpiicx on current
# oneAPI releases, or the classic mpiifort/mpiicc on older ones) and the
# Intel MKL for ScaLAPACK/BLACS -- there is no generic gfortran/conda
# path for it the way there is for ModEM and Occam2D. pycsamt already
# has a complete Python-side manager for exactly this
# (pycsamt.models.mare2dem.SourceManager: download via git/archive,
# generate a Makefile include, build, and locate the resulting
# binary); this script is a thin, honest wrapper around it, not a
# reimplementation.
#
# Usage:
#   pycsamt/models/_solver_build/mare2dem.sh [options]
#
# Options:
#   --auto-install     Actually run download()+build() (see below for
#                       why this is opt-in). Without it, this script
#                       only reports status and what it would do.
#   -y, --yes          No extra effect beyond --auto-install here
#                       (kept for flag consistency with the other
#                       scripts in this directory); SourceManager
#                       itself does not prompt interactively.
#   --clean             Pass clean_first=True to SourceManager.build().
#   --source-dir DIR    Use an existing source tree instead of
#                       downloading (passed to SourceManager(source_dir=...)).
#   -h, --help          Show this help and exit.
#
# Why --auto-install is required even to check things: downloading
# ~49 MB from Bitbucket and compiling a large Fortran/MPI codebase is
# a heavier, longer-running action than the other scripts in this
# directory perform by default; this script always shows
# SourceManager's own status report first and only proceeds to
# download/build when explicitly asked.
#
# Windows: MARE2DEM cannot be built natively -- SourceManager.build()
# raises immediately on sys.platform == "win32". Run this script from
# inside WSL2 instead; it detects a native Windows/Git-Bash/MSYS
# environment and stops with that same guidance before even trying.
#
# Exit status: 0 on success (or when only reporting status), 1 on any
# failure, matching SourceManager's own exceptions.

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "$SCRIPT_DIR/lib/common.sh"

AUTO_INSTALL=0
ASSUME_YES=0
DO_CLEAN=0
SOURCE_DIR_OVERRIDE=""

print_help() {
    awk '/^set -u/{exit} NR>1 && /^#/{sub(/^# ?/,""); print}' "${BASH_SOURCE[0]}"
}

while [ $# -gt 0 ]; do
    case "$1" in
        --auto-install) AUTO_INSTALL=1 ;;
        -y|--yes) ASSUME_YES=1 ;;
        --clean) DO_CLEAN=1 ;;
        --source-dir) shift; SOURCE_DIR_OVERRIDE="${1:-}" ;;
        -h|--help) print_help; exit 0 ;;
        *) die "Unknown option: $1 (see --help)" ;;
    esac
    shift
done
# ASSUME_YES is accepted for CLI consistency across these scripts but
# has no separate effect here; reference it so shellcheck/set -u
# don't flag it as unused.
: "$ASSUME_YES"

if [ "$(pycsamt_build_os)" = "windows" ]; then
    die "MARE2DEM cannot be built natively on Windows (Intel MKL/ScaLAPACK toolchain is Linux/macOS-oriented). Run this script from inside WSL2 instead: wsl bash pycsamt/models/_solver_build/mare2dem.sh $*"
fi

have_cmd python3 && PY=python3 || PY=python
have_cmd "$PY" || die "No Python interpreter found on PATH."

log_step "MARE2DEM status (via pycsamt.models.mare2dem.SourceManager)"
"$PY" - "$SOURCE_DIR_OVERRIDE" <<'PYEOF'
import sys
from pycsamt.models.mare2dem import SourceManager

source_dir = sys.argv[1] or None
sm = SourceManager(source_dir=source_dir, verbose=1)
sm.print_status()
PYEOF

if [ "$AUTO_INSTALL" != "1" ]; then
    log_info "Status only (pass --auto-install to actually download + build)."
    exit 0
fi

log_step "Downloading (if needed) and building MARE2DEM -- this can take several minutes."
CLEAN_FLAG="False"
[ "$DO_CLEAN" = "1" ] && CLEAN_FLAG="True"

"$PY" - "$SOURCE_DIR_OVERRIDE" "$CLEAN_FLAG" <<'PYEOF'
import sys
from pycsamt.models.mare2dem import SourceManager

source_dir = sys.argv[1] or None
clean_first = sys.argv[2] == "True"

sm = SourceManager(source_dir=source_dir, verbose=1)
sm.download()
binary = sm.build(clean_first=clean_first)
print(f"\nBuilt: {binary}")
print('Point pycsamt at it via Mare2DEMConfig(binary=str(binary)), or leave')
print('the default "MARE2DEM" and rely on SourceManager.resolve_binary().')
PYEOF
