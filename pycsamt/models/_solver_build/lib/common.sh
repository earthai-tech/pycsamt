#!/usr/bin/env bash
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
#
# Shared helpers for pycsamt/models/_solver_build/*.sh (ModEM, Occam2D,
# MARE2DEM compile scripts). Not meant to be run directly.
#
# Sourced with:
#   COMMON_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
#   source "$COMMON_DIR/lib/common.sh"
#
# What this provides
# -------------------
# * coloured logging (info/warn/error/step), auto-disabled when stdout
#   is not a terminal or $NO_COLOR is set;
# * OS detection (linux / macos / windows-via-git-bash-or-msys / wsl);
# * dependency probing for a Fortran toolchain (gfortran, make, and a
#   linkable LAPACK/BLAS) without assuming any one platform's package
#   manager;
# * an *opt-in*, explicit installer for that toolchain per platform
#   (apt/dnf/yum on Linux, Homebrew on macOS, a dedicated conda
#   environment on Windows -- gfortran/make are not otherwise
#   reliably available there without MSYS2/WSL);
# * a resilient `make` runner that retries with `-k` (keep-going)
#   across a handful of passes -- the vendored ModEM Makefiles' object
#   lists are not in strict Fortran-module dependency order (a
#   pre-existing property of the upstream/adapted Makefiles, not
#   something these scripts can safely reorder), so a single-pass
#   `make` can fail with "module file not found" purely because of
#   build order, not a real error. Confirmed empirically while first
#   compiling this source in 2026-07-29: the 3-D Makefile needed 4
#   `make -k` passes to converge even after every missing file/rule
#   was fixed. Never pass `-j` (parallel) for this reason.

set -u

# --------------------------------------------------------------------------
# Logging
# --------------------------------------------------------------------------

if [ -t 1 ] && [ -z "${NO_COLOR:-}" ]; then
    _C_RESET=$'\033[0m'
    _C_INFO=$'\033[36m'
    _C_WARN=$'\033[33m'
    _C_ERR=$'\033[31m'
    _C_STEP=$'\033[1;32m'
else
    _C_RESET=""; _C_INFO=""; _C_WARN=""; _C_ERR=""; _C_STEP=""
fi

log_info()  { printf '%s[info]%s  %s\n'  "$_C_INFO" "$_C_RESET" "$*"; }
log_warn()  { printf '%s[warn]%s  %s\n'  "$_C_WARN" "$_C_RESET" "$*" >&2; }
log_error() { printf '%s[error]%s %s\n'  "$_C_ERR"  "$_C_RESET" "$*" >&2; }
log_step()  { printf '%s==>%s %s\n'      "$_C_STEP" "$_C_RESET" "$*"; }

die() { log_error "$*"; exit 1; }

# --------------------------------------------------------------------------
# OS / environment detection
# --------------------------------------------------------------------------

# Echoes one of: linux, macos, windows, unknown
pycsamt_build_os() {
    case "$(uname -s 2>/dev/null || echo unknown)" in
        Linux*)   echo "linux" ;;
        Darwin*)  echo "macos" ;;
        MINGW*|MSYS*|CYGWIN*) echo "windows" ;;
        *)        echo "unknown" ;;
    esac
}

have_cmd() { command -v "$1" >/dev/null 2>&1; }

# Best-effort path to an activated (or discoverable) conda installation's
# `conda` executable; empty string if none is found.
_conda_exe() {
    if have_cmd conda; then
        command -v conda
        return 0
    fi
    for candidate in "$HOME/anaconda3/condabin/conda" \
                      "$HOME/miniconda3/condabin/conda" \
                      "/c/Users/$USER/anaconda3/condabin/conda"; do
        if [ -x "$candidate" ]; then
            echo "$candidate"
            return 0
        fi
    done
    return 1
}

# Convert a Windows-style path (C:\Users\...) to the POSIX form Git
# Bash/MSYS expects (/c/Users/...) -- passes non-Windows paths through
# unchanged. This matters because `conda env list` always prints
# Windows-style paths on Windows regardless of the shell, and splicing
# "C:\foo\bin" directly into a colon-separated PATH string is invalid:
# the colon after the drive letter is read as a *second* PATH entry
# separator, silently corrupting PATH rather than raising an error.
# That corruption was the actual root cause of an intermittent-looking
# "internal compiler error: Segmentation fault" from gfortran while
# first writing this script: gfortran internally re-execs its own
# cc1/collect2/as/ld helpers by searching PATH, and a corrupted PATH
# let it pick up mismatched helpers from elsewhere on the system.
_to_posix_path() {
    local path="$1"
    # Deliberately not using `cygpath`: it is not guaranteed to be the
    # MSYS variant (Git Bash's own), and a Cygwin-flavoured `cygpath`
    # emits "/cygdrive/c/..." rather than the "/c/..." this shell
    # actually needs -- confirmed to silently produce an unusable path
    # while first writing this script. Plain parameter expansion is
    # unambiguous and has no such variant to get wrong.
    case "$path" in
        [A-Za-z]:\\*|[A-Za-z]:/*)
            local drive="${path:0:1}"
            local rest="${path:2}"
            rest="${rest//\\//}"
            printf '/%s%s\n' "$(printf '%s' "$drive" | tr 'A-Z' 'a-z')" "$rest"
            ;;
        *)
            printf '%s\n' "$path"
            ;;
    esac
}

# Prefix directory of a named conda environment, in POSIX form (see
# _to_posix_path), or empty if it does not exist. Requires `conda` to
# be resolvable via _conda_exe.
conda_env_prefix() {
    local env_name="$1" conda_bin raw
    conda_bin="$(_conda_exe)" || return 1
    raw="$("$conda_bin" env list 2>/dev/null | awk -v name="$env_name" \
        '$1 == name { print $NF; found=1 } END { exit !found }')" || return 1
    [ -n "$raw" ] && _to_posix_path "$raw"
}

# --------------------------------------------------------------------------
# Fortran toolchain detection
# --------------------------------------------------------------------------

# Default name for the dedicated Windows conda environment these
# scripts create/use for gfortran+make+OpenBLAS. Override with
# PYCSAMT_FORTRAN_ENV.
PYCSAMT_FORTRAN_ENV="${PYCSAMT_FORTRAN_ENV:-pycsamt-fortran}"

# Populates the following globals for the caller (not returned, since
# bash has no clean multi-value return):
#   TOOLCHAIN_OK        1 if a usable toolchain was found, else 0
#   TOOLCHAIN_F90        path/name to invoke as F90 (e.g. "gfortran")
#   TOOLCHAIN_MAKE       path/name to invoke as make
#   TOOLCHAIN_LIBS       LIBS= override to pass to `make` (may be empty
#                        to keep the Makefile's own default)
#   TOOLCHAIN_BIN_DIR    directory to prepend to PATH before running,
#                        or empty if nothing extra is needed
detect_fortran_toolchain() {
    TOOLCHAIN_OK=0
    TOOLCHAIN_F90=""
    TOOLCHAIN_MAKE=""
    TOOLCHAIN_LIBS=""
    TOOLCHAIN_BIN_DIR=""

    local os
    os="$(pycsamt_build_os)"

    if [ "$os" = "windows" ]; then
        # gfortran/make are not reliably available on Windows outside
        # a dedicated MinGW-w64 toolchain; look for our conda env
        # first, then fall back to whatever is already on PATH (e.g.
        # a user-installed MSYS2).
        local prefix
        if prefix="$(conda_env_prefix "$PYCSAMT_FORTRAN_ENV")" && [ -n "$prefix" ]; then
            local gfortran_bin="$prefix/Library/mingw-w64/bin/gfortran.exe"
            local make_bin="$prefix/Library/bin/make.exe"
            local openblas_lib="$prefix/Library/mingw-w64/lib/libopenblas.a"
            if [ -x "$gfortran_bin" ] && [ -x "$make_bin" ] && [ -f "$openblas_lib" ]; then
                TOOLCHAIN_OK=1
                TOOLCHAIN_F90="$gfortran_bin"
                TOOLCHAIN_MAKE="$make_bin"
                TOOLCHAIN_LIBS="-L$prefix/Library/mingw-w64/lib -lopenblas"
                TOOLCHAIN_BIN_DIR="$prefix/Library/mingw-w64/bin"
                return 0
            fi
        fi
        if have_cmd gfortran && have_cmd make; then
            # A non-conda MinGW-w64/MSYS2 toolchain already on PATH.
            # We cannot know its LAPACK/BLAS story, so leave LIBS at
            # the Makefile's own default and let a link failure speak
            # for itself with a clear message.
            TOOLCHAIN_OK=1
            TOOLCHAIN_F90="gfortran"
            TOOLCHAIN_MAKE="make"
            return 0
        fi
        return 0
    fi

    # Linux / macOS: trust PATH. LAPACK/BLAS are checked at link time
    # (there is no single reliable cross-distro way to probe for them
    # up front); the Makefile's own `-llapack -lblas` default is left
    # untouched.
    have_cmd make || return 0
    TOOLCHAIN_MAKE="make"

    if have_cmd gfortran; then
        TOOLCHAIN_OK=1
        TOOLCHAIN_F90="gfortran"
        return 0
    fi
    # Homebrew's gcc formula on macOS does not always leave an
    # unversioned "gfortran" symlink on PATH (in particular when more
    # than one GCC version is installed side by side) -- only
    # "gfortran-14", "gfortran-13", etc. Search recent versions before
    # giving up, newest first.
    local ver candidate
    for ver in 15 14 13 12 11 10 9; do
        candidate="gfortran-$ver"
        if have_cmd "$candidate"; then
            TOOLCHAIN_OK=1
            TOOLCHAIN_F90="$candidate"
            return 0
        fi
    done
}

# --------------------------------------------------------------------------
# Opt-in toolchain installation
# --------------------------------------------------------------------------

# Prefix a package-manager command with "sudo" unless already root or
# sudo does not exist (e.g. minimal containers that run as root only).
_sudo_prefix() {
    if [ "$(id -u 2>/dev/null || echo 1000)" = "0" ] || ! have_cmd sudo; then
        echo ""
    else
        echo "sudo "
    fi
}

# install_fortran_toolchain <assume_yes:0|1>
# Prints what it is about to do and, unless assume_yes=1, asks for
# confirmation before running anything. Never called unless the
# caller script was given --auto-install.
install_fortran_toolchain() {
    local assume_yes="$1"
    local os sudo_
    os="$(pycsamt_build_os)"
    sudo_="$(_sudo_prefix)"

    case "$os" in
        linux)
            if have_cmd apt-get; then
                _confirm_and_run "$assume_yes" \
                    "${sudo_}apt-get update && ${sudo_}apt-get install -y gfortran liblapack-dev libblas-dev"
            elif have_cmd dnf; then
                _confirm_and_run "$assume_yes" \
                    "${sudo_}dnf install -y gcc-gfortran lapack-devel blas-devel"
            elif have_cmd yum; then
                _confirm_and_run "$assume_yes" \
                    "${sudo_}yum install -y gcc-gfortran lapack-devel blas-devel"
            else
                die "No supported package manager found (looked for apt-get, dnf, yum). Install gfortran + LAPACK + BLAS manually, then re-run."
            fi
            ;;
        macos)
            have_cmd brew || die "Homebrew not found (https://brew.sh). Install gfortran + OpenBLAS manually, then re-run."
            _confirm_and_run "$assume_yes" "brew install gcc openblas"
            ;;
        windows)
            local conda_bin
            conda_bin="$(_conda_exe)" || die "conda was not found on PATH. Install Miniconda/Anaconda, or gfortran+make via MSYS2, then re-run."
            _confirm_and_run "$assume_yes" \
                "\"$conda_bin\" create -y -n $PYCSAMT_FORTRAN_ENV -c conda-forge m2w64-gcc-fortran m2w64-openblas make"
            ;;
        *)
            die "Unrecognized platform; install gfortran + make + LAPACK/BLAS manually, then re-run."
            ;;
    esac
}

_confirm_and_run() {
    local assume_yes="$1" cmd="$2"
    log_step "About to run:"
    printf '    %s\n' "$cmd"
    if [ "$assume_yes" != "1" ]; then
        printf 'Proceed? [y/N] '
        read -r reply
        case "$reply" in
            y|Y|yes|YES) ;;
            *) die "Aborted by user." ;;
        esac
    fi
    eval "$cmd" || die "Installation command failed: $cmd"
}

# --------------------------------------------------------------------------
# Resilient make runner
# --------------------------------------------------------------------------

# run_make_with_retries <max_passes> <make_args...>
# Runs `make` (never with -j: see the module-header note on why) up
# to <max_passes> times with -k (keep-going past errors), stopping as
# soon as a pass has zero errors. Each pass can fix a subset of a
# previous pass's "module file not found" failures once the .mod files
# that pass produced exist. Returns 0 on success, 1 if it never
# converges (last pass's output is left in $MAKE_LAST_OUTPUT for the
# caller to show).
run_make_with_retries() {
    local max_passes="$1"
    shift
    local pass
    for ((pass = 1; pass <= max_passes; pass++)); do
        MAKE_LAST_OUTPUT="$("$@" -k 2>&1)"
        local errors
        errors=$(printf '%s\n' "$MAKE_LAST_OUTPUT" | grep -c "Fatal Error\|Error [0-9]" || true)
        log_info "make pass $pass/$max_passes: $errors error line(s) remaining"
        if [ "$errors" -eq 0 ]; then
            return 0
        fi
    done
    return 1
}

# --------------------------------------------------------------------------
# Self-contained Windows deployment
# --------------------------------------------------------------------------

# copy_runtime_dlls <bin_dir> <dest_dir>
# On Windows, a MinGW-w64-built executable needs its runtime DLLs
# either on PATH or next to itself (Windows checks the executable's
# own directory first). Copies the ones this toolchain actually needs
# if found in <bin_dir>; a no-op (with a warning) if <bin_dir> is empty
# or the DLLs are not there.
copy_runtime_dlls() {
    local bin_dir="$1" dest_dir="$2"
    if [ -z "$bin_dir" ] || [ ! -d "$bin_dir" ]; then
        log_warn "No known MinGW-w64 bin directory to copy runtime DLLs from; the built binary may fail to launch unless its DLLs are already on PATH."
        return 0
    fi
    local dll copied=0
    for dll in libgcc_s_seh-1.dll libgfortran-3.dll libopenblas.dll \
               libquadmath-0.dll libwinpthread-1.dll; do
        if [ -f "$bin_dir/$dll" ]; then
            cp -f "$bin_dir/$dll" "$dest_dir/"
            copied=$((copied + 1))
        fi
    done
    log_info "Copied $copied runtime DLL(s) into $dest_dir for a self-contained build."
}

# --------------------------------------------------------------------------
# Argument helpers shared across the per-solver scripts
# --------------------------------------------------------------------------

# print_common_help <solver_name>
print_common_help() {
    cat <<EOF
Common options (all $1 build scripts):
  --auto-install     Attempt to install a missing Fortran toolchain
                      (gfortran, make, LAPACK/BLAS) for this platform.
                      Always prints the exact command first and asks
                      for confirmation unless -y/--yes is also given.
  -y, --yes           Assume "yes" to any confirmation prompt (only
                      relevant together with --auto-install).
  --clean             Run "make clean" before building.
  -h, --help          Show this help and exit.
EOF
}
