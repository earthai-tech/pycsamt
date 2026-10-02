#!/usr/bin/env bash
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
#
# Build the pycsamt-desktop standalone binary (Linux/macOS).
#
# Usage (from anywhere):
#   bash packaging/pyinstaller/build_desktop.sh
#
# Requires: the environment that has pycsamt installed together with its
# `desktop` AND `agents` extras (PySide6, pyqtgraph, contextily, torch,
# scikit-learn, the LLM provider clients) -- this build cannot exclude
# torch/tensorflow the way the converter build does, see
# pycsamt_desktop.spec's own docstring -- plus `pyinstaller` itself
# (pip install pyinstaller). See packaging/pyinstaller/README.md for
# details and troubleshooting.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
SPEC_FILE="$SCRIPT_DIR/pycsamt_desktop.spec"

cd "$REPO_ROOT"

if ! python -c "import PyInstaller" >/dev/null 2>&1; then
    echo "PyInstaller is not installed in this Python environment. Run: pip install pyinstaller" >&2
    exit 1
fi

echo "Cleaning previous build artifacts..."
rm -rf "$REPO_ROOT/build/pycsamt-desktop" "$REPO_ROOT/dist/pycsamt-desktop"

echo "Running PyInstaller (this build is substantially larger than the converter's -- expect several minutes)..."
pyinstaller --noconfirm --clean "$SPEC_FILE"

EXE_PATH="$REPO_ROOT/dist/pycsamt-desktop/pycsamt-desktop"
if [ -e "$EXE_PATH" ]; then
    echo
    echo "Build succeeded: $EXE_PATH"
    echo "Next: wrap it with packaging/linux/build_installer.sh for a real self-extracting installer."
else
    echo "Build finished but $EXE_PATH was not found." >&2
    exit 1
fi
