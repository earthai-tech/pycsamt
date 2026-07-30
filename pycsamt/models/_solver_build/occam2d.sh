#!/usr/bin/env bash
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
#
# Build Occam2DMT (deGroot-Hedlin & Constable) from the vendored
# source at pycsamt/models/occam2d/_source/.
#
# Usage:
#   pycsamt/models/_solver_build/occam2d.sh [options]
#
# Options:
#   --auto-install     Install a missing Fortran toolchain (see
#                       lib/common.sh for exactly what that means per
#                       OS). Off by default. Occam2D itself does not
#                       need LAPACK/BLAS (its own linear algebra is
#                       self-contained; a commented-out fast-LAPACK
#                       path exists in Occam.f90 for anyone who wants
#                       to enable it by hand), but the shared
#                       auto-installer provisions the same toolchain
#                       as the ModEM scripts for consistency.
#   -y, --yes          Skip the confirmation prompt --auto-install
#                       would otherwise show.
#   --clean            Remove previous object/module files and the
#                       binary before building.
#   --prefix DIR        After a successful build, also copy the
#                       binary into DIR.
#   -h, --help          Show this help and exit.
#
# The vendored Makefile's own default compiler ("f90") is a legacy
# alias unlikely to exist on a modern system; this script always
# overrides it to the detected gfortran.
#
# Exit status: 0 on a successful build, 1 otherwise.

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "$SCRIPT_DIR/lib/common.sh"

SOURCE_DIR="$(cd "$SCRIPT_DIR/../occam2d/_source" && pwd)"
OUT_NAME="Occam2D"

AUTO_INSTALL=0
ASSUME_YES=0
DO_CLEAN=0
PREFIX_DIR=""

print_help() {
    awk '/^set -u/{exit} NR>1 && /^#/{sub(/^# ?/,""); print}' "${BASH_SOURCE[0]}"
}

while [ $# -gt 0 ]; do
    case "$1" in
        --auto-install) AUTO_INSTALL=1 ;;
        -y|--yes) ASSUME_YES=1 ;;
        --clean) DO_CLEAN=1 ;;
        --prefix) shift; PREFIX_DIR="${1:-}" ;;
        -h|--help) print_help; exit 0 ;;
        *) die "Unknown option: $1 (see --help)" ;;
    esac
    shift
done

log_step "Building Occam2D from $SOURCE_DIR"

detect_fortran_toolchain
if [ "$TOOLCHAIN_OK" != "1" ]; then
    if [ "$AUTO_INSTALL" = "1" ]; then
        log_step "No usable Fortran toolchain found; installing (--auto-install)."
        install_fortran_toolchain "$ASSUME_YES"
        detect_fortran_toolchain
    fi
fi
if [ "$TOOLCHAIN_OK" != "1" ]; then
    die "No usable Fortran toolchain (gfortran + make) found. Re-run with --auto-install, or install one manually."
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

log_step "Compiling"
if ! run_make_with_retries 3 "$TOOLCHAIN_MAKE" \
        "FC90=$TOOLCHAIN_F90" "FCFLAGS=-O2 -ffree-line-length-none"; then
    log_error "Build failed. Last pass's output:"
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
log_info "Point pycsamt.models.occam2d.runner at this binary's path."
