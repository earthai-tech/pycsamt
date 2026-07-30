#!/usr/bin/env bash
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
#
# Build ModEM's serial 2-D forward/inversion solver (Mod2DMT) from the
# vendored source at pycsamt/models/modem/_source/2D/.
#
# Usage:
#   pycsamt/models/_solver_build/modem2d.sh [options]
#
# Options:
#   --auto-install     Install a missing Fortran toolchain for this
#                       platform (see lib/common.sh for exactly what
#                       that means per OS). Off by default.
#   -y, --yes          Skip the confirmation prompt --auto-install
#                       would otherwise show.
#   --clean            Run "make clean" first.
#   --mpi              Build the MPI-parallel variant instead of the
#                       serial default (needs mpif90 on PATH; not
#                       covered by --auto-install on any platform).
#   --intel            Build with the Intel Fortran compiler (ifort)
#                       instead of gfortran. You are responsible for
#                       ifort already being on PATH; --auto-install
#                       never installs a commercial compiler.
#   --prefix DIR        After a successful build, also copy the
#                       binary (and, on Windows, its runtime DLLs)
#                       into DIR.
#   -h, --help          Show this help and exit.
#
# What this actually does (so nothing here is a black box):
#   1. Detects (or, with --auto-install, installs) a Fortran toolchain.
#   2. Runs `make clean` if requested.
#   3. Runs `make` (or `make mpi` / `make intel`), retrying with -k
#      (keep-going) across a few passes -- the vendored Makefile's
#      object list is not in strict Fortran-module dependency order,
#      a pre-existing property fixed just enough (missing files/rules)
#      to make it *convergent* under repeated -k passes, not
#      reordered outright. See lib/common.sh's module docstring.
#   4. On success, prints the resulting binary's path and (Windows
#      only) copies its MinGW-w64 runtime DLLs alongside it so it can
#      run without activating any environment.
#
# Exit status: 0 on a successful build, 1 otherwise (with a specific
# message -- missing toolchain, install declined, or a build failure
# with the tail of the last `make` pass's output).

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "$SCRIPT_DIR/lib/common.sh"

SOURCE_DIR="$(cd "$SCRIPT_DIR/../modem/_source/2D" && pwd)"
OUT_NAME="Mod2DMT"

AUTO_INSTALL=0
ASSUME_YES=0
DO_CLEAN=0
BUILD_TARGET=""       # "" = serial default, "mpi", "intel"
PREFIX_DIR=""

print_help() {
    awk '/^set -u/{exit} NR>1 && /^#/{sub(/^# ?/,""); print}' "${BASH_SOURCE[0]}"
}

while [ $# -gt 0 ]; do
    case "$1" in
        --auto-install) AUTO_INSTALL=1 ;;
        -y|--yes) ASSUME_YES=1 ;;
        --clean) DO_CLEAN=1 ;;
        --mpi) BUILD_TARGET="mpi" ;;
        --intel) BUILD_TARGET="intel" ;;
        --prefix) shift; PREFIX_DIR="${1:-}" ;;
        -h|--help) print_help; exit 0 ;;
        *) die "Unknown option: $1 (see --help)" ;;
    esac
    shift
done

log_step "Building ModEM Mod2DMT from $SOURCE_DIR"

detect_fortran_toolchain
if [ "$TOOLCHAIN_OK" != "1" ]; then
    if [ "$AUTO_INSTALL" = "1" ]; then
        log_step "No usable Fortran toolchain found; installing (--auto-install)."
        install_fortran_toolchain "$ASSUME_YES"
        detect_fortran_toolchain
    fi
fi
if [ "$TOOLCHAIN_OK" != "1" ]; then
    die "No usable Fortran toolchain (gfortran + make) found. Re-run with --auto-install, or install one manually -- see pycsamt/models/modem/_source/README."
fi
log_info "Using F90=$TOOLCHAIN_F90  make=$TOOLCHAIN_MAKE"

if [ -n "$TOOLCHAIN_BIN_DIR" ]; then
    export PATH="$TOOLCHAIN_BIN_DIR:$PATH"
fi

cd "$SOURCE_DIR" || die "Cannot cd into $SOURCE_DIR"

if [ "$DO_CLEAN" = "1" ]; then
    log_step "make clean"
    "$TOOLCHAIN_MAKE" clean
fi

MAKE_ARGS=("$TOOLCHAIN_MAKE" "F90=$TOOLCHAIN_F90")
[ -n "$TOOLCHAIN_LIBS" ] && MAKE_ARGS+=("LIBS=$TOOLCHAIN_LIBS")
if [ "$BUILD_TARGET" = "mpi" ]; then
    have_cmd mpif90 || die "--mpi requested but mpif90 was not found on PATH."
    MAKE_ARGS+=("mpi")
elif [ "$BUILD_TARGET" = "intel" ]; then
    have_cmd ifort || die "--intel requested but ifort was not found on PATH."
    MAKE_ARGS+=("intel")
fi

log_step "Compiling (this can take a minute; a few make -k passes is normal)"
if ! run_make_with_retries 6 "${MAKE_ARGS[@]}"; then
    log_error "Build did not converge after 6 make -k passes. Last pass's output:"
    printf '%s\n' "$MAKE_LAST_OUTPUT" | tail -n 40
    exit 1
fi

BIN_PATH="$SOURCE_DIR/$OUT_NAME"
[ -f "$BIN_PATH.exe" ] && BIN_PATH="$BIN_PATH.exe"
[ -f "$BIN_PATH" ] || die "make reported success but $OUT_NAME was not found in $SOURCE_DIR."

if [ "$(pycsamt_build_os)" = "windows" ]; then
    copy_runtime_dlls "$TOOLCHAIN_BIN_DIR" "$SOURCE_DIR"
fi

if [ -n "$PREFIX_DIR" ]; then
    mkdir -p "$PREFIX_DIR"
    cp -f "$BIN_PATH" "$PREFIX_DIR/"
    [ "$(pycsamt_build_os)" = "windows" ] && cp -f "$SOURCE_DIR"/*.dll "$PREFIX_DIR/" 2>/dev/null
    log_info "Also copied to $PREFIX_DIR"
fi

log_step "Built: $BIN_PATH"
log_info "Point pycsamt at it via ModEmConfig(binary_2d=\"$BIN_PATH\"), or leave the default"
log_info "\"Mod2DMT\" and rely on this directory being pycsamt's own search-path fallback."
