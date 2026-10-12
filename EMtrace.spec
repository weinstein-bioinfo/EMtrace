# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller spec file for EMtrace.

Build with:
    pyinstaller EMtrace.spec --clean --noconfirm

macOS:   produces  dist/EMtrace-{arch}.app  +  dist/EMtrace-macos-{arch}-{ver}.dmg
Windows: produces  dist/EMtrace-x86_64/EMtrace.exe  +  dist/EMtrace-win-x86_64-{ver}.zip

Prefer the wrapper scripts, which set the venv up first:
    bash scripts/build_app.sh            (macOS, current arch)
    .\scripts\build_app.ps1              (Windows)
"""

import os
import platform
import sys
from pathlib import Path

block_cipher = None

# ------------------------------------------------------------------
# Paths  (EMtrace uses a src/ layout: src/emtrace/...)
# ------------------------------------------------------------------
PROJECT_ROOT = Path(SPECPATH)  # noqa: F821  (SPECPATH injected by PyInstaller)
SRC_ROOT = PROJECT_ROOT / 'src'
PKG_ROOT = SRC_ROOT / 'emtrace'
ENTRY_POINT = str(PROJECT_ROOT / 'scripts' / 'pyinstaller_entry.py')
ICON_ICNS = str(PKG_ROOT / 'gui' / 'resources' / 'icons' / 'app.icns')
ICON_ICO = str(PKG_ROOT / 'gui' / 'resources' / 'icons' / 'app.ico')

# ------------------------------------------------------------------
# Architecture detection (macOS: arm64 vs x86_64; Windows: x86_64)
# ------------------------------------------------------------------
IS_DARWIN = sys.platform == 'darwin'
IS_WIN = sys.platform == 'win32'

_machine = platform.machine()          # 'arm64', 'x86_64', or 'AMD64'
_override = os.environ.get('TARGET_ARCH')
BUILD_ARCH = _override if _override in ('arm64', 'x86_64') else _machine

# Normalise the Windows arch tag
if BUILD_ARCH.upper() == 'AMD64':
    BUILD_ARCH = 'x86_64'

# Resolve the venv from the active Python executable
_venv = Path(sys.executable).parent.parent
if IS_WIN:
    VENV = _venv
elif '.venv-x86' in str(_venv):
    BUILD_ARCH = 'x86_64'
    VENV = _venv
elif '.venv' in str(_venv):
    BUILD_ARCH = 'arm64'
    VENV = _venv
else:
    VENV = PROJECT_ROOT / ('.venv-x86' if BUILD_ARCH == 'x86_64' else '.venv')

if IS_WIN:
    SITE_PACKAGES = VENV / 'Lib' / 'site-packages'
else:
    SITE_PACKAGES = (
        VENV / 'lib' / f'python3.{sys.version_info.minor}' / 'site-packages'
    )

# ------------------------------------------------------------------
# Version (read from pyproject.toml)
# ------------------------------------------------------------------
import tomllib

with open(PROJECT_ROOT / 'pyproject.toml', 'rb') as _f:
    APP_VERSION = tomllib.load(_f)['project']['version']

print(f"[spec] BUILD_ARCH={BUILD_ARCH}  VENV={VENV}  SITE_PACKAGES={SITE_PACKAGES}")

# ------------------------------------------------------------------
# Data files to bundle
# ------------------------------------------------------------------
datas = [
    # Packaged default configuration — EMtraceConfig.create_default() reads
    # these through importlib.resources, so they must keep their package path.
    (str(PKG_ROOT / 'config' / 'default.yaml'),
     os.path.join('emtrace', 'config')),
    (str(PKG_ROOT / 'config' / 'default.hostdbs.yaml'),
     os.path.join('emtrace', 'config')),
]

# GUI resources (icons); created by scripts/generate_icon.py
_gui_resources = PKG_ROOT / 'gui' / 'resources'
if _gui_resources.is_dir():
    datas.append((str(_gui_resources), os.path.join('emtrace', 'gui', 'resources')))

# Optional package data: report.py looks up emtrace/data/ICTV.bphage.filter.tsv
_pkg_data = PKG_ROOT / 'data'
if _pkg_data.is_dir():
    datas.append((str(_pkg_data), os.path.join('emtrace', 'data')))

# User manual — shipped so the GUI's help can open it offline
for _doc in ('UserGuide.v1.0.md', 'UserGuide.v1.0.pdf'):
    _p = PROJECT_ROOT / 'docs' / _doc
    if _p.is_file():
        datas.append((str(_p), 'docs'))

# ------------------------------------------------------------------
# Pre-build the matplotlib font cache and ship it as bundled data, so the
# first launch does not pause to scan system fonts.
#
# This must happen BEFORE Analysis: BUNDLE copies the files listed in
# COLLECT's table of contents, so a directory dropped into the collect tree
# afterwards never reaches the .app.  Routing it through `datas` is what makes
# it land in Contents/Frameworks (= sys._MEIPASS) as well as in the Windows
# collect directory.
#
# The cache is portable: matplotlib's own faces are stored as paths relative
# to mpl-data, and the absolute entries point at /System/Library/Fonts, which
# exists on every macOS install.
# ------------------------------------------------------------------
import subprocess as _subprocess

_fontcache_dir = PROJECT_ROOT / 'build' / 'matplotlib-fontcache'
_fontcache_dir.mkdir(parents=True, exist_ok=True)
_fc = _subprocess.run(
    [sys.executable, '-c',
     'import matplotlib.font_manager; '
     'matplotlib.font_manager._load_fontmanager(try_read_cache=False)'],
    env={**os.environ, 'MPLCONFIGDIR': str(_fontcache_dir)},
    check=False, capture_output=True,
)
if any(_fontcache_dir.glob('fontlist-*.json')):
    datas.append((str(_fontcache_dir), 'matplotlib-fontcache'))
    print(f"[spec] matplotlib font cache pre-built in {_fontcache_dir}")
else:
    print(f"[spec] WARNING: font cache not built; first launch will be slow. "
          f"{_fc.stderr.decode(errors='replace')[-400:]}")

# ------------------------------------------------------------------
# Collect packages that use dynamic imports / ship data PyInstaller misses
# ------------------------------------------------------------------
from PyInstaller.utils.hooks import collect_all, collect_submodules

binaries = []
hidden_imports = []

# kaleido is deliberately absent: pyinstaller-hooks-contrib ships a working
# hook-kaleido.py, and collect_all() cannot import kaleido.mocker (its _args
# module calls parser.parse_args() at import time and exits).
for _pkg in ('plotly', 'PyComplexHeatmap', 'pysam'):
    try:
        _d, _b, _h = collect_all(_pkg)
    except Exception as _exc:                       # package not installed
        print(f"[spec] skipping optional package {_pkg!r}: {_exc}")
        continue
    datas += _d
    binaries += _b
    hidden_imports += _h

# ------------------------------------------------------------------
# Hidden imports — modules PyInstaller cannot detect statically
# ------------------------------------------------------------------
# Every emtrace submodule: steps register themselves via the
# @register_step decorator at import time, and the GUI widgets/workers are
# resolved lazily, so a static scan misses most of them.
hidden_imports += collect_submodules('emtrace')

hidden_imports += [
    # ---- third-party libraries often missed by PyInstaller ----
    'pydantic',
    'pydantic.deprecated',
    'pydantic.deprecated.decorator',
    'pydantic_core',
    'pkg_resources',
    'Bio',
    'Bio.SeqIO',
    'Bio.Seq',
    'Bio.Align',
    'matplotlib',
    'matplotlib.backends.backend_agg',
    'matplotlib.backends.backend_pdf',
    'matplotlib.backends.backend_svg',
    'matplotlib.backends.backend_qtagg',
    'seaborn',
    'openpyxl',
    'watchdog',
    'watchdog.observers',
    'watchdog.observers.fsevents',
    'watchdog.observers.polling',
    'watchdog.events',
    'tqdm',
    'rich',
    'yaml',
    'click',
    'colorlog',
    'numpy',
    'pandas',
]

# ------------------------------------------------------------------
# Analysis
# ------------------------------------------------------------------
a = Analysis(
    [ENTRY_POINT],
    pathex=[str(SRC_ROOT), str(PROJECT_ROOT)],
    binaries=binaries,
    datas=datas,
    hiddenimports=hidden_imports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        'tkinter',
        # ---- Qt modules EMtrace never instantiates ----
        # QtWebEngine is ~200 MB; figure_viewer.py already falls back to the
        # system browser when QWebEngineView is unavailable.  Drop the exclude
        # (and pip install PyQt6-WebEngine) if you want in-app HTML plots.
        'PyQt6.QtWebEngineWidgets',
        'PyQt6.QtWebEngineCore',
        'PyQt6.QtWebEngineQuick',
        'PyQt6.QtQuick',
        'PyQt6.QtQml',
        'PyQt6.Qt3DCore',
        'PyQt6.QtBluetooth',
        'PyQt6.QtNfc',
        'PyQt6.QtMultimedia',
        'PyQt6.QtDesigner',
        'PyQt6.QtTest',
        # ---- Development tools (not runtime deps) ----
        'pytest',
        'IPython',
    ],
    noarchive=False,
    optimize=0,
    cipher=block_cipher,
)

# ------------------------------------------------------------------
# PYZ (Python bytecode archive)
# ------------------------------------------------------------------
pyz = PYZ(a.pure, cipher=block_cipher)

# ------------------------------------------------------------------
# EXE
# ------------------------------------------------------------------
exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='EMtrace',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,           # --windowed (no terminal on macOS/Windows)
    disable_windowed_traceback=False,
    argv_emulation=IS_DARWIN,          # macOS only: accept file drops via argv
    target_arch=BUILD_ARCH if IS_DARWIN else None,
    codesign_identity=None,
    entitlements_file=None,
    icon=ICON_ICNS if IS_DARWIN else ICON_ICO,
)

# ------------------------------------------------------------------
# COLLECT (directory-based bundle)
# ------------------------------------------------------------------
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name=f'EMtrace-{BUILD_ARCH}',
)

# ------------------------------------------------------------------
# Post-COLLECT cleanup — trim unused data from the collected directory
# ------------------------------------------------------------------
import glob as _glob
import shutil

_collect_dir = os.path.join(DISTPATH, f'EMtrace-{BUILD_ARCH}')  # noqa: F821

# Remove matplotlib sample_data — never accessed at runtime
_sample_data = os.path.join(_collect_dir, 'matplotlib', 'mpl-data', 'sample_data')
if os.path.isdir(_sample_data):
    shutil.rmtree(_sample_data)

# Remove dist-info directories (~1 MB) — pip metadata, not needed in the bundle
for _d in _glob.glob(os.path.join(_collect_dir, '*.dist-info')):
    shutil.rmtree(_d)

# ------------------------------------------------------------------
# macOS .app BUNDLE
# ------------------------------------------------------------------
if IS_DARWIN:
    app = BUNDLE(
        coll,
        name=f'EMtrace-{BUILD_ARCH}.app',
        icon=ICON_ICNS,
        bundle_identifier='fi.helsinki.emtrace.app',
        info_plist={
            'CFBundleName': 'EMtrace',
            'CFBundleDisplayName': 'EMtrace',
            'CFBundleShortVersionString': APP_VERSION,
            'CFBundleVersion': APP_VERSION,
            'NSHighResolutionCapable': True,
            'NSPrincipalClass': 'NSApplication',
            'NSAppleScriptEnabled': False,
            'LSMinimumSystemVersion': '11.0',
        },
    )

    # --------------------------------------------------------------
    # Post-BUNDLE: create the versioned .dmg
    # e.g. dist/EMtrace-macos-arm64-1.0.0.dmg
    # --------------------------------------------------------------
    _app_path = os.path.join(DISTPATH, f'EMtrace-{BUILD_ARCH}.app')  # noqa: F821
    _dmg_path = os.path.join(                                            # noqa: F821
        DISTPATH, f'EMtrace-macos-{BUILD_ARCH}-{APP_VERSION}.dmg',
    )
    _subprocess.run(
        ['hdiutil', 'create', '-volname', 'EMtrace',
         '-srcfolder', _app_path, '-ov', '-format', 'UDZO', _dmg_path],
        check=True,
    )
    print(f"DMG created: {_dmg_path}")

# ------------------------------------------------------------------
# Windows: create the versioned .zip
# e.g. dist/EMtrace-win-x86_64-1.0.0.zip
# ------------------------------------------------------------------
if IS_WIN:
    import zipfile as _zipfile

    _zip_path = os.path.join(                                            # noqa: F821
        DISTPATH, f'EMtrace-win-{BUILD_ARCH}-{APP_VERSION}.zip',
    )
    _src_dir = os.path.join(DISTPATH, f'EMtrace-{BUILD_ARCH}')        # noqa: F821
    with _zipfile.ZipFile(_zip_path, 'w', _zipfile.ZIP_DEFLATED) as _zf:
        for _root, _dirs, _files in os.walk(_src_dir):
            for _file in _files:
                _abs = os.path.join(_root, _file)
                _arc = os.path.join(
                    f'EMtrace-{BUILD_ARCH}', os.path.relpath(_abs, _src_dir),
                )
                _zf.write(_abs, _arc)
    print(f"ZIP created: {_zip_path} "
          f"({os.path.getsize(_zip_path) / (1024 * 1024):.1f} MB)")
