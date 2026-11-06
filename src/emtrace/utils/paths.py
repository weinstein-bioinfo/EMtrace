"""Path resolution utilities with environment variable expansion."""

import os
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, Optional


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


@contextmanager
def working_dir(path: Path) -> Iterator[Path]:
    """Temporarily change the process working directory to *path*.

    Used by the GUI workers: an app launched from Finder/Spotlight starts in
    the read-only ``/``, and external tools that write temp files to the
    working directory (porechop_abi's ``TEMP_<pid>.fastq``) fail there.
    The working directory is process-wide, so only one worker should use
    this at a time.
    """
    path = ensure_dir(Path(path).resolve())
    previous = os.getcwd()
    os.chdir(path)
    try:
        yield path
    finally:
        os.chdir(previous)
