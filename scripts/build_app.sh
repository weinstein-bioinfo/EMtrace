#!/usr/bin/env bash
# Build the EMtrace standalone macOS .app bundle using PyInstaller.
#
# Usage:
#   bash scripts/build_app.sh            # build for the current arch (arm64)
#   bash scripts/build_app.sh x86_64     # build x86_64 via Rosetta 2
#
# Output:
#   dist/EMtrace-{arch}.app
#   dist/EMtrace-macos-{arch}-{version}.dmg
#
# Prerequisites:
#   - Python 3.10-3.12 (3.12 recommended: widest wheel coverage for
#     PyQt6 / pysam / kaleido)
#   - For x86_64: Rosetta 2 and an x86_64 Homebrew Python, in .venv-x86
#
# The venv is created and populated automatically if missing.

set -euo pipefail
cd "$(dirname "$0")/.."

TARGET="${1:-}"   # optional: "x86_64" or "arm64"
PYTHON_BIN="${PYTHON_BIN:-python3.12}"

echo "=== EMtrace macOS App Bundle Builder ==="
echo ""

# ------------------------------------------------------------------
# Select venv and arch prefix
# ------------------------------------------------------------------
if [ "$TARGET" = "x86_64" ]; then
    VENV=".venv-x86"
    ARCH_CMD="arch -x86_64"
    ARCH="x86_64"
elif [ "$TARGET" = "arm64" ] || [ -z "$TARGET" ]; then
    VENV=".venv"
    ARCH_CMD=""
    ARCH="arm64"
else
    echo "ERROR: unknown target '$TARGET'. Use arm64 or x86_64."
    exit 1
fi

PYTHON="$VENV/bin/python3"
PIP="$VENV/bin/pip"
PYINSTALLER="$VENV/bin/pyinstaller"

# ------------------------------------------------------------------
# Create the venv on first run
# ------------------------------------------------------------------
if [ ! -f "$PYTHON" ]; then
    if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
        echo "ERROR: $PYTHON_BIN not found. Install it, or set"
        echo "       PYTHON_BIN=python3.11 bash scripts/build_app.sh"
        exit 1
    fi
    echo "Creating virtual environment in $VENV ($($PYTHON_BIN --version)) ..."
    $ARCH_CMD "$PYTHON_BIN" -m venv "$VENV"
fi

echo "Target arch: $ARCH"
echo "Venv: $VENV ($("$PYTHON" --version))"
echo ""

# ------------------------------------------------------------------
# Ensure dependencies
# ------------------------------------------------------------------
echo "Installing / updating dependencies ..."
$ARCH_CMD "$PIP" install --upgrade pip --quiet
$ARCH_CMD "$PIP" install -e ".[gui,plots,excel]" --quiet
$ARCH_CMD "$PIP" install pyinstaller pillow --quiet
echo "Dependencies ready."
echo ""

# ------------------------------------------------------------------
# Ensure the app icon exists
# ------------------------------------------------------------------
ICON="src/emtrace/gui/resources/icons/app.icns"
if [ ! -f "$ICON" ]; then
    echo "App icon not found — generating ..."
    $ARCH_CMD "$PYTHON" scripts/generate_icon.py
    echo ""
fi

# ------------------------------------------------------------------
# Build
# ------------------------------------------------------------------
echo "Building EMtrace-${ARCH}.app ..."
echo ""
TARGET_ARCH="$ARCH" $ARCH_CMD "$PYINSTALLER" EMtrace.spec --clean --noconfirm

echo ""
echo "=== Build complete ==="

APP="dist/EMtrace-${ARCH}.app"
if [ ! -d "$APP" ]; then
    echo "ERROR: $APP was not created."
    exit 1
fi

# Ad-hoc code signing for local use (no Developer ID / notarisation)
echo "Signing $APP (ad-hoc) ..."
xattr -cr "$APP"
codesign --sign - --force --deep --timestamp=none "$APP"

SIZE=$(du -sh "$APP" | cut -f1)
echo ""
echo "  APP: $APP ($SIZE)"

VERSION=$(grep '^version' pyproject.toml | head -1 | sed 's/.*"\(.*\)"/\1/')
DMG="dist/EMtrace-macos-${ARCH}-${VERSION}.dmg"
if [ -f "$DMG" ]; then
    echo "  DMG: $DMG ($(du -sh "$DMG" | cut -f1))"
fi

echo ""
echo "To launch:"
echo "  First launch: unsigned build, see the User Guide, 'Opening the app'"
echo "  (macOS 15+: System Settings → Privacy & Security → Open Anyway)"
echo "  open $APP            (subsequent launches)"
echo ""
echo "CLI inside the bundle:"
echo "  $APP/Contents/MacOS/EMtrace --help"
