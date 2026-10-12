"""Path resolution utilities with environment variable expansion."""

import os
from pathlib import Path
from typing import Optional


def expand_path(path: str | Path) -> Path:
    """Expand environment variables ($VAR) and user home (~) in a path."""
    return Path(os.path.expandvars(os.path.expanduser(str(path))))


def resolve_path(path: str | Path, base: Optional[Path] = None) -> Path:
    """Expand and resolve a path, optionally relative to *base*."""
    expanded = expand_path(path)
    if not expanded.is_absolute() and base is not None:
        expanded = base / expanded
    return expanded.resolve()


def ensure_dir(path: Path) -> Path:
    """Create directory (and parents) if it doesn't exist; return the path."""
    path.mkdir(parents=True, exist_ok=True)
    return path
