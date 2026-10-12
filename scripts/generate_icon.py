#!/usr/bin/env python3
"""Generate the EMtrace app icon (.icns for macOS, .ico for Windows).

Uses matplotlib to render a simple branded icon, then converts it with
macOS ``iconutil`` (.icns) or Pillow (.ico).

Usage:
    python scripts/generate_icon.py
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.patches as mpatches       # noqa: E402
import matplotlib.pyplot as plt             # noqa: E402
import numpy as np                          # noqa: E402


ICON_DIR = (
    Path(__file__).resolve().parent.parent
    / 'src' / 'emtrace' / 'gui' / 'resources' / 'icons'
)

# Required sizes for a macOS .iconset (filename -> pixels)
ICONSET_SIZES = {
    'icon_16x16.png': 16,
    'icon_16x16@2x.png': 32,
    'icon_32x32.png': 32,
    'icon_32x32@2x.png': 64,
    'icon_128x128.png': 128,
    'icon_128x128@2x.png': 256,
    'icon_256x256.png': 256,
    'icon_256x256@2x.png': 512,
    'icon_512x512.png': 512,
    'icon_512x512@2x.png': 1024,
}


def render_icon(size: int, path: Path) -> None:
    """Render the EMtrace icon at the given pixel size."""
    dpi = 100
    fig, ax = plt.subplots(1, 1, figsize=(size / dpi, size / dpi), dpi=dpi)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    ax.axis('off')

    # Background: rounded square in the EMtrace blue
    bg = mpatches.FancyBboxPatch(
        (0.02, 0.02), 0.96, 0.96,
        boxstyle='round,pad=0.08',
        facecolor='#1f4e79',
        edgecolor='#154360',
        linewidth=max(1, size / 128),
    )
    ax.add_patch(bg)

    # A nanopore read squiggle — the raw-signal motif, drawn as a trace
    rng = np.random.default_rng(7)
    t = np.linspace(0.12, 0.88, 240)
    squiggle = (
        0.5
        + 0.16 * np.sin(2 * np.pi * 2.1 * t)
        + 0.06 * np.sin(2 * np.pi * 6.7 * t + 1.1)
        + 0.02 * rng.standard_normal(t.size)
    )
    ax.plot(
        t, squiggle,
        color='#5dade2', linewidth=max(1, size / 40),
        alpha=0.85, solid_capstyle='round',
    )

    # Text: "EMt" centred over the squiggle
    ax.text(
        0.50, 0.50, 'EMt',
        fontsize=size * 0.23, fontweight='bold',
        color='white', ha='center', va='center',
        fontfamily='sans-serif',
    )

    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    fig.savefig(path, dpi=dpi, transparent=True, pad_inches=0)
    plt.close(fig)


def _create_ico(png_path: Path, ico_path: Path) -> None:
    """Create a Windows .ico file from a PNG using Pillow."""
    from PIL import Image

    img = Image.open(png_path)
    img.save(ico_path, format='ICO',
             sizes=[(s, s) for s in (16, 32, 48, 64, 128, 256)])
    print(f'Created: {ico_path}')


def main() -> None:
    ICON_DIR.mkdir(parents=True, exist_ok=True)

    # Always save a 512px PNG for other uses and as the ICO source
    png_path = ICON_DIR / 'app.png'
    render_icon(512, png_path)
    print(f'Created: {png_path}')

    if sys.platform == 'win32':
        _create_ico(png_path, ICON_DIR / 'app.ico')
    else:
        # macOS: build an .iconset directory and convert it to .icns
        with tempfile.TemporaryDirectory() as tmp:
            iconset = Path(tmp) / 'EMtrace.iconset'
            iconset.mkdir()
            for filename, size in ICONSET_SIZES.items():
                render_icon(size, iconset / filename)
                print(f'  {filename} ({size}x{size})')

            icns_path = ICON_DIR / 'app.icns'
            subprocess.run(
                ['iconutil', '-c', 'icns', str(iconset), '-o', str(icns_path)],
                check=True,
            )
            print(f'\nCreated: {icns_path}')

        # Also emit an .ico so a Windows build can reuse it if Pillow is
        # unavailable there; harmless on macOS.
        try:
            _create_ico(png_path, ICON_DIR / 'app.ico')
        except ImportError:
            print('Pillow not installed - skipping app.ico')


if __name__ == '__main__':
    main()
