#!/usr/bin/env bash
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
#
# Self-extracting pycsamt-desktop installer (Linux). This file is a
# template: build_installer.sh appends a tar.gz of dist/pycsamt-desktop/
# after the __PYCSAMT_ARCHIVE_BELOW__ marker line and writes the result to
# dist/installer/ -- do not run this template file directly, it has no
# payload attached.
#
# What running the finished installer does:
#   - Extracts the onedir PyInstaller build to
#     ${PYCSAMT_PREFIX:-$HOME/.local/share/pycsamt-desktop}
#   - Installs a launcher script to ~/.local/bin/pycsamt-desktop
#   - Installs a .desktop entry + icon so the app shows up in a normal
#     desktop-environment application menu, not just a terminal command
#   - Supports --uninstall to remove exactly what it installed
#
# Trial/license state (~/.pycsamt, QSettings' Linux equivalent under
# ~/.config/earthai-tech/) is deliberately left alone by --uninstall --
# same reasoning as the Windows Inno Setup script's [UninstallDelete]
# comment: removing it would silently reset a user's trial or drop their
# license key the moment they reinstall.

set -euo pipefail

APP_NAME="pycsamt-desktop"
PREFIX="${PYCSAMT_PREFIX:-$HOME/.local/share/$APP_NAME}"
BIN_DIR="$HOME/.local/bin"
DESKTOP_FILE_DIR="$HOME/.local/share/applications"
ICON_DIR="$HOME/.local/share/icons/hicolor/256x256/apps"

_usage() {
    echo "Usage: $0 [--prefix DIR] [--uninstall]"
    echo "  --prefix DIR   Install location (default: \$HOME/.local/share/$APP_NAME,"
    echo "                 or \$PYCSAMT_PREFIX if set)"
    echo "  --uninstall    Remove a previous install (launcher, .desktop entry,"
    echo "                 icon, and the extracted app directory). Trial/license"
    echo "                 state under ~/.pycsamt and ~/.config/earthai-tech is"
    echo "                 left untouched."
}

_uninstall() {
    echo "Removing $APP_NAME..."
    rm -rf "$PREFIX"
    rm -f "$BIN_DIR/$APP_NAME"
    rm -f "$DESKTOP_FILE_DIR/$APP_NAME.desktop"
    rm -f "$ICON_DIR/$APP_NAME.png"
    echo "Done. Trial/license state under ~/.pycsamt and ~/.config/earthai-tech was left in place."
    exit 0
}

DO_UNINSTALL=0
while [ $# -gt 0 ]; do
    case "$1" in
        --prefix)
            PREFIX="$2"
            shift 2
            ;;
        --uninstall)
            DO_UNINSTALL=1
            shift
            ;;
        -h|--help)
            _usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            _usage
            exit 1
            ;;
    esac
done

if [ "$DO_UNINSTALL" -eq 1 ]; then
    _uninstall
fi

ARCHIVE_MARKER_LINE=$(awk '/^__PYCSAMT_ARCHIVE_BELOW__$/{print NR + 1; exit}' "$0")
if [ -z "$ARCHIVE_MARKER_LINE" ]; then
    echo "This installer has no archive payload attached -- it looks like" >&2
    echo "the unmodified template, not a build produced by build_installer.sh." >&2
    exit 1
fi

echo "Installing $APP_NAME to $PREFIX ..."
mkdir -p "$PREFIX" "$BIN_DIR" "$DESKTOP_FILE_DIR" "$ICON_DIR"
tail -n +"$ARCHIVE_MARKER_LINE" "$0" | tar xzf - -C "$PREFIX" --strip-components=1

cat > "$BIN_DIR/$APP_NAME" <<EOF
#!/usr/bin/env bash
exec "$PREFIX/$APP_NAME" "\$@"
EOF
chmod +x "$BIN_DIR/$APP_NAME"

# PyInstaller 6.x onedir layout nests bundled `datas` under _internal/,
# preserving the destination path given in pycsamt_desktop.spec's `datas`
# tuple -- not at the extracted tree's root.
_BUNDLED_ICON="$PREFIX/_internal/pycsamt/app/desktop/resources/icons/pycsamt_256.png"
if [ -f "$_BUNDLED_ICON" ]; then
    cp "$_BUNDLED_ICON" "$ICON_DIR/$APP_NAME.png"
fi

cat > "$DESKTOP_FILE_DIR/$APP_NAME.desktop" <<EOF
[Desktop Entry]
Type=Application
Name=pyCSAMT
Comment=Geophysical processing suite for MT/AMT/CSAMT/CSEM
Exec=$BIN_DIR/$APP_NAME
Icon=$APP_NAME
Terminal=false
Categories=Science;Education;
EOF

echo ""
echo "Installed. Launch with '$APP_NAME' (if ~/.local/bin is on your PATH) or from your"
echo "application menu. Uninstall with: $0 --uninstall"

exit 0
__PYCSAMT_ARCHIVE_BELOW__
