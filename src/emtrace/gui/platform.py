"""Cross-platform helpers for the EMtrace GUI.

Handles per-OS tweaks (PATH augmentation for Finder-launched apps on
macOS, taskbar grouping on Windows), a forced light Fusion theme so
the app looks identical under system dark-mode, and monospace font
selection for the log viewer.
"""

from __future__ import annotations

import json
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import List

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor, QFont, QFontDatabase, QPalette
from PyQt6.QtWidgets import QApplication, QStyleFactory


# ---------------------------------------------------------------------------
# Theme
# ---------------------------------------------------------------------------

def _light_palette() -> QPalette:
    """Light-theme QPalette matching the app's QSS colors."""
    p = QPalette()
    white = QColor("#ffffff")
    bg = QColor("#f5f5f5")
    mid = QColor("#e0e0e0")
    dark = QColor("#bbbbbb")
    text = QColor("#222222")
    distext = QColor("#aaaaaa")
    hilite = QColor("#cce5ff")
    hiltxt = QColor("#000000")

    for role, colour in (
        (QPalette.ColorRole.Window, bg),
        (QPalette.ColorRole.WindowText, text),
        (QPalette.ColorRole.Base, white),
        (QPalette.ColorRole.AlternateBase, mid),
        (QPalette.ColorRole.Text, text),
        (QPalette.ColorRole.BrightText, text),
        (QPalette.ColorRole.Button, mid),
        (QPalette.ColorRole.ButtonText, text),
        (QPalette.ColorRole.Highlight, hilite),
        (QPalette.ColorRole.HighlightedText, hiltxt),
        (QPalette.ColorRole.ToolTipBase, white),
        (QPalette.ColorRole.ToolTipText, text),
        (QPalette.ColorRole.PlaceholderText, distext),
        (QPalette.ColorRole.Mid, dark),
        (QPalette.ColorRole.Dark, dark),
        (QPalette.ColorRole.Shadow, QColor("#888888")),
        (QPalette.ColorRole.Link, QColor("#0057ae")),
        (QPalette.ColorRole.LinkVisited, QColor("#7f007f")),
    ):
        p.setColor(QPalette.ColorGroup.All, role, colour)

    # Disabled group — dim text
    p.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.WindowText, distext)
    p.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.Text, distext)
    p.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.ButtonText, distext)
    return p


# ---------------------------------------------------------------------------
# Bundle / resource paths
# ---------------------------------------------------------------------------

def bundle_root() -> Path:
    """Return the project root data directory.

    In a PyInstaller frozen bundle this is ``sys._MEIPASS``.
    In a development checkout it is the repository root (three levels
    above this file: ``gui/platform.py`` → ``emtrace`` → ``src`` →
    repo root).
    """
    if getattr(sys, "frozen", False):
        return Path(sys._MEIPASS)  # type: ignore[attr-defined]
    return Path(__file__).resolve().parent.parent.parent.parent


# ---------------------------------------------------------------------------
# PATH augmentation (macOS Finder-launched apps lose Homebrew/conda)
# ---------------------------------------------------------------------------

def _augment_path_macos() -> None:
    """Prepend common tool directories to PATH on macOS.

    When a bundled app launches from Finder the process inherits a
    minimal PATH (``/usr/bin:/bin:/usr/sbin:/sbin``) that omits
    Homebrew, conda, and MacPorts.  This ensures ``shutil.which`` can
    locate minimap2, samtools, BLAST+, diamond, etc.
    """
    extra = [
        "/opt/homebrew/bin",           # Apple Silicon Homebrew
        "/usr/local/bin",              # Intel Homebrew / manual installs
        "/opt/local/bin",              # MacPorts
        "/usr/local/ncbi/blast/bin",   # NCBI standalone installer
        str(Path.home() / "bin"),
        str(Path.home() / ".local" / "bin"),
        str(Path.home() / "miniconda3" / "bin"),
        str(Path.home() / "mambaforge" / "bin"),
    ]
    current = os.environ.get("PATH", "")
    existing = set(current.split(":"))
    additions = [p for p in extra if p not in existing and Path(p).is_dir()]
    if additions:
        os.environ["PATH"] = ":".join(additions) + (":" + current if current else "")


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def configure_platform(app: QApplication) -> None:
    """Apply per-OS tweaks after QApplication creation.

    - Forces the Fusion style so Qt draws every widget itself, overriding
      the OS dark-mode theme (keeps the app's look consistent on
      macOS/Windows/Linux).
    - Applies the light palette.
    - macOS: augments PATH so Finder-launched bundles can find CLI tools.
    - Windows: sets an explicit AppUserModelID for correct taskbar grouping.
    """
    if "Fusion" in QStyleFactory.keys():
        app.setStyle("Fusion")
    app.setPalette(_light_palette())

    if sys.platform == "darwin":
        _augment_path_macos()
        app.setAttribute(
            Qt.ApplicationAttribute.AA_DontShowIconsInMenus, False,
        )
    elif sys.platform == "win32":
        try:
            import ctypes

            ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(
                "ngs.emtrace.gui",
            )
        except Exception:
            # Not fatal — taskbar grouping just falls back to default.
            pass


# ---------------------------------------------------------------------------
# Fonts
# ---------------------------------------------------------------------------

def monospace_font(size: int = 10) -> QFont:
    """Return the system fixed-pitch font at *size* pt.

    Uses ``QFontDatabase.systemFont(FixedFont)`` so no font-family string
    lookup is needed — avoids the slow alias-population warning on macOS.
    """
    font = QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont)
    font.setPointSize(size)
    return font


# ---------------------------------------------------------------------------
# Per-user data directory (recent samples, session state)
# ---------------------------------------------------------------------------

def _data_dir() -> Path:
    """Return the platform-appropriate per-user application data directory."""
    if sys.platform == "win32":
        base = Path(os.environ.get("APPDATA", str(Path.home())))
    elif sys.platform == "darwin":
        base = Path.home() / "Library" / "Application Support"
    else:
        base = Path(
            os.environ.get("XDG_DATA_HOME", str(Path.home() / ".local" / "share"))
        )
    d = base / "EMtrace"
    d.mkdir(parents=True, exist_ok=True)
    return d


def recent_file() -> Path:
    """Path to the recent-samples JSON file."""
    return _data_dir() / "recent.json"


def load_recent() -> List[dict]:
    """Load recent-sample entries from disk (most recent first)."""
    rf = recent_file()
    if rf.exists():
        try:
            return json.loads(rf.read_text())
        except (json.JSONDecodeError, OSError):
            return []
    return []


def add_recent(
    path: str,
    name: str,
    input_reads: str = "",
    max_entries: int = 10,
) -> None:
    """Add a sample to the recent list (most recent first, capped)."""
    entries = load_recent()
    entries = [e for e in entries if e.get("path") != path]
    entries.insert(0, {
        "path": path,
        "name": name,
        "input_reads": input_reads,
        "last_opened": datetime.now().isoformat(timespec="seconds"),
    })
    entries = entries[:max_entries]
    recent_file().write_text(
        json.dumps(entries, indent=2), encoding="utf-8", newline="\n",
    )
