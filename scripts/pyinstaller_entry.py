#!/usr/bin/env python3
"""Minimal entry point for PyInstaller builds of EMtrace.

Dispatches to CLI or GUI based on sys.argv:
  - choreographer's Chrome wrapper script as the first argument -> run that
    wrapper (see below), never the GUI
  - No arguments, or launched by macOS as a .app (Finder passes no
    subcommand, or a legacy -psn_* argument) -> launch the GUI
  - Anything else -> run the CLI (which reports unknown commands), so a
    helper call of the frozen binary can never open a second GUI window
"""
import multiprocessing
import os
import sys

# kaleido (Plotly PNG/PDF export) starts Chrome through choreographer, which on
# macOS/Linux runs `sys.executable _unix_pipe_chromium_wrapper.py <chrome> ...`.
# In a frozen bundle sys.executable is this binary, so run the wrapper module
# (bundled as a hidden import in EMtrace.spec) instead of starting EMtrace.
# Must stay first: the wrapper dup2()s fds 3/4 before anything else opens files.
_CHROME_WRAPPER = '_unix_pipe_chromium_wrapper.py'
if len(sys.argv) > 1 and os.path.basename(sys.argv[1]) == _CHROME_WRAPPER:
    import runpy
    sys.argv = sys.argv[1:]
    runpy.run_module(
        'choreographer.browsers._unix_pipe_chromium_wrapper',
        run_name='__main__', alter_sys=True,
    )
    sys.exit(0)

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
    """Return True unless the invocation is a plain app launch.

    Only no arguments (or macOS's legacy -psn_* process serial number) opens
    the GUI; any other argument goes to the CLI, which handles subcommands
    and rejects unknown ones instead of opening a second GUI window.
    """
    args = sys.argv[1:]
    return bool(args) and not args[0].startswith('-psn_')


if _wants_cli():
    from emtrace.__main__ import main
    main()
else:
    from emtrace.gui.main_window import launch_gui
    sys.exit(launch_gui())
