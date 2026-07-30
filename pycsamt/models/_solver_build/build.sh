#!/usr/bin/env bash
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
#
# Single entry point for pycsamt/models/_solver_build/*.sh.
#
# Usage:
#   pycsamt/models/_solver_build/build.sh <solver> [options...]
#
# <solver> is one of: modem2d, modem3d, occam2d, mare2dem
# [options...] are forwarded verbatim to that solver's own script --
# see modem2d.sh/modem3d.sh/occam2d.sh/mare2dem.sh (or run
# `build.sh <solver> --help`) for what each accepts.
#
# This script does nothing besides dispatch; each solver script is
# fully usable on its own and this wrapper exists only so a new user
# has one obvious place to start.

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat <<EOF
Usage: $(basename "${BASH_SOURCE[0]}") <solver> [options...]

  solver    modem2d | modem3d | occam2d | mare2dem

Run '$(basename "${BASH_SOURCE[0]}") <solver> --help' for that
solver's own options.

Examples:
  $(basename "${BASH_SOURCE[0]}") modem3d --auto-install -y
  $(basename "${BASH_SOURCE[0]}") occam2d --clean
  $(basename "${BASH_SOURCE[0]}") mare2dem            # status only
EOF
}

if [ $# -lt 1 ]; then
    usage
    exit 1
fi

solver="$1"
shift

case "$solver" in
    modem2d|modem3d|occam2d|mare2dem)
        exec "$SCRIPT_DIR/$solver.sh" "$@"
        ;;
    -h|--help)
        usage
        exit 0
        ;;
    *)
        echo "Unknown solver: $solver" >&2
        usage
        exit 1
        ;;
esac
