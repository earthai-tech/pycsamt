#!/usr/bin/env bash
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
#
# Build ModEM's serial 3-D forward/inversion solver (Mod3DMT) from the
# vendored source at pycsamt/models/modem/_source/3D/.
#
# Usage:
#   pycsamt/models/_solver_build/modem3d.sh [options]
#
# Options: identical to modem2d.sh -- see that script's header, or
# run with --help.
#
# Known ModEM constraint this script cannot build around (it is in
# the vendored Fortran itself, not in this project's own adapter --
# see pycsamt/forward/maxwell/modem3d.py's module docstring): ModEM's
# WS-format model reader hardcodes 10 air layers and its default
# "mirror" air-layer sizing reads that many earth-layer widths with no
# bounds check against the actual earth cell count. Any model with
# fewer than 10 earth z-cells will compile and run, then crash mid-
# solve with "b in QMR contains NaNs" -- this is a property of the
# *model* you later run, not of the build this script performs, but
# it is worth knowing before your first real solve.
#
# Exit status: 0 on a successful build, 1 otherwise.

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "$SCRIPT_DIR/lib/common.sh"

SOURCE_DIR="$(cd "$SCRIPT_DIR/../modem/_source/3D" && pwd)"
OUT_NAME="Mod3DMT"

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

log_step "Building ModEM Mod3DMT from $SOURCE_DIR"

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
log_info "Point pycsamt at it via ModEmConfig(binary_3d=\"$BIN_PATH\"), or leave the default"
log_info "\"Mod3DMT\" and rely on this directory being pycsamt's own search-path fallback."
log_warn "Remember: models need >= 10 earth z-cells (ModEM's own hardcoded"
log_warn "air-layer count) or the solver will crash -- see this script's header."
