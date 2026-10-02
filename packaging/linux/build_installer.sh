#!/usr/bin/env bash
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
#
# Build the pycsamt-desktop self-extracting Linux installer from the
# PyInstaller onedir build (see installer_template.sh for what the
# finished installer actually does).
#
# Usage (from anywhere):
#   bash packaging/linux/build_installer.sh
#
# Requires: dist/pycsamt-desktop/ already built (this script builds it
# first via ../pyinstaller/build_desktop.sh if it's missing -- but does
# NOT rebuild it if it's already there, since that build takes several
# minutes; delete dist/pycsamt-desktop/ yourself first to force a
# rebuild).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
DIST_DIR="$REPO_ROOT/dist/pycsamt-desktop"
TEMPLATE="$SCRIPT_DIR/installer_template.sh"
OUT_DIR="$REPO_ROOT/dist/installer"

cd "$REPO_ROOT"

if [ ! -e "$DIST_DIR/pycsamt-desktop" ]; then
    echo "dist/pycsamt-desktop/pycsamt-desktop not found -- building it first..."
    bash "$REPO_ROOT/packaging/pyinstaller/build_desktop.sh"
fi

VERSION="$(python -c "import pycsamt; print(pycsamt.__version__)")"
if [ -z "$VERSION" ]; then
    echo "Could not read pycsamt.__version__ from the active Python environment." >&2
    exit 1
fi
echo "Building installer for version $VERSION..."

mkdir -p "$OUT_DIR"
OUT_FILE="$OUT_DIR/pycsamt-desktop-$VERSION-linux-x86_64.sh"

cat "$TEMPLATE" > "$OUT_FILE"
# --strip-components=1 (in installer_template.sh) expects a single
# top-level directory inside the archive -- -C dist pycsamt-desktop
# gives it exactly that, matching dist/pycsamt-desktop/'s own onedir
# layout (the pycsamt-desktop executable + _internal/).
tar czf - -C "$REPO_ROOT/dist" pycsamt-desktop >> "$OUT_FILE"
chmod +x "$OUT_FILE"

echo ""
echo "Installer built: $OUT_FILE"
echo "Install with:   bash $OUT_FILE"
echo "Uninstall with: bash $OUT_FILE --uninstall"
