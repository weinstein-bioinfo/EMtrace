"""Input validation utilities."""

import shutil
from pathlib import Path

from emtrace.core.exceptions import DependencyError, SampleFileNotFoundError


def check_tool(name: str, install_hint: str = "") -> str:
    """Verify an external tool is on PATH; return its full path."""
    path = shutil.which(name)
    if path is None:
        raise DependencyError(name, install_hint)
    return path


def check_file(path: Path | str, step: str = "") -> Path:
    """Verify a file exists; return its Path."""
    p = Path(path)
    if not p.exists():
        raise SampleFileNotFoundError(str(p), step)
    return p
