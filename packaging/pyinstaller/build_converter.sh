#!/usr/bin/env bash
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
#
# Build the pyCSAMT Format Studio standalone binary (Linux/macOS).
#
# Usage (from anywhere):
#   bash packaging/pyinstaller/build_converter.sh
#
# Requires: the environment that has pycsamt + PySide6 installed also
# has `pyinstaller` (pip install pyinstaller). See
# packaging/pyinstaller/README.md for details and troubleshooting.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
SPEC_FILE="$SCRIPT_DIR/pycsamt_converter.spec"

cd "$REPO_ROOT"

if ! python -c "import PyInstaller" >/dev/null 2>&1; then
    echo "PyInstaller is not installed in this Python environment. Run: pip install pyinstaller" >&2
    exit 1
fi

echo "Cleaning previous build artifacts..."
rm -rf "$REPO_ROOT/build/pycsamt-converter" "$REPO_ROOT/dist/pycsamt-converter"

echo "Running PyInstaller..."
pyinstaller --noconfirm --clean "$SPEC_FILE"

EXE_PATH="$REPO_ROOT/dist/pycsamt-converter/pycsamt-converter"
if [ -e "$EXE_PATH" ]; then
    echo
    echo "Build succeeded: $EXE_PATH"
else
    echo "Build finished but $EXE_PATH was not found." >&2
    exit 1
fi
