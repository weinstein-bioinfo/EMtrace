#!/usr/bin/env python3
"""Minimal entry point for PyInstaller builds of EMtrace.

Dispatches to CLI or GUI based on sys.argv:
  - No arguments, or launched by macOS as a .app (Finder passes no
    subcommand) -> launch the GUI
  - Any CLI subcommand present (run, databases, watch, batch-summary,
    batch-taxa, batch-abundance, gui, --version, --help) -> run the CLI
"""
import multiprocessing
import os
import sys

# Required for PyInstaller --windowed on Windows: prevents child processes from
# re-executing the entry point and re-launching the GUI (window flash/reopen).
multiprocessing.freeze_support()

# Point matplotlib to the pre-built font cache bundled inside _MEIPASS, so the
# first launch does not pause to scan every system font.
#
# This must be a plain assignment, not setdefault(): PyInstaller's own
# pyi_rth_mplconfig runtime hook runs before this script and unconditionally
# points MPLCONFIGDIR at a fresh temporary directory, which would leave
# setdefault() a no-op and rebuild the cache on every single launch.
if hasattr(sys, '_MEIPASS'):
    _bundled_cache = os.path.join(sys._MEIPASS, 'matplotlib-fontcache')
    if os.path.isdir(_bundled_cache):
        os.environ['MPLCONFIGDIR'] = _bundled_cache


def _augment_path() -> None:
    """Prepend common tool directories to PATH.

    When launched from Finder / a desktop shortcut the process inherits a
    minimal PATH that omits Homebrew, conda, and MacPorts.  EMtrace shells
    out to minimap2, samtools, diamond, porechop, chopper, medaka, etc., so
    shutil.which() must be able to find them.

    Mirrors emtrace.gui.platform._augment_path_macos(), which only runs for
    the GUI; this covers CLI invocations of the frozen binary too.
    """
    from pathlib import Path

    if sys.platform == 'darwin':
        extra = [
            '/opt/homebrew/bin',                  # Apple Silicon Homebrew
            '/usr/local/bin',                     # Intel Homebrew / manual
            '/opt/local/bin',                     # MacPorts
            '/usr/local/ncbi/blast/bin',          # NCBI standalone installer
            str(Path.home() / 'bin'),
            str(Path.home() / '.local' / 'bin'),
            str(Path.home() / 'miniconda3' / 'bin'),
            str(Path.home() / 'mambaforge' / 'bin'),
            str(Path.home() / 'anaconda3' / 'bin'),
        ]
        sep = ':'
    elif sys.platform == 'win32':
        import glob as _glob
        # NCBI BLAST+ installs into a versioned folder: blast-2.16.0+\bin
        _blast = sorted(
            _glob.glob(r'C:\Program Files\NCBI\blast-*\bin')
            + _glob.glob(r'C:\Program Files (x86)\NCBI\blast-*\bin')
        )
        extra = _blast + [
            str(Path.home() / 'miniconda3' / 'Scripts'),
            str(Path.home() / 'mambaforge' / 'Scripts'),
            str(Path.home() / 'anaconda3' / 'Scripts'),
        ]
        sep = ';'
    else:
        return

    current = os.environ.get('PATH', '')
    existing = set(current.split(sep))
    additions = [p for p in extra if p not in existing and os.path.isdir(p)]
    if additions:
        os.environ['PATH'] = sep.join(additions) + (sep + current if current else '')


_augment_path()


def _ensure_home() -> None:
    """Set HOME on Windows if missing.

    Several Unix-originated libraries read os.environ['HOME'] directly, which
    raises KeyError on Windows.  Map it to USERPROFILE.
    """
    if sys.platform == 'win32' and 'HOME' not in os.environ:
        os.environ['HOME'] = os.environ.get('USERPROFILE') or os.path.expanduser('~')


_ensure_home()


def _wants_cli() -> bool:
    """Return True if the invocation looks like a CLI call."""
    args = sys.argv[1:]
    if not args:
        return False
    cli_commands = {
        'run', 'databases', 'watch', 'gui',
        'batch-summary', 'batch-taxa', 'batch-abundance',
        '--version', '--help', '-h',
    }
    return args[0] in cli_commands


if _wants_cli():
    from emtrace.__main__ import main
    main()
else:
    from emtrace.gui.main_window import launch_gui
    sys.exit(launch_gui())
