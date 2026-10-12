"""Per-user persistent config for EMtrace.

Stores GUI/session preferences (e.g. db_dir) in a human-readable YAML file
under the platform-appropriate per-user data directory:

  macOS:   ~/Library/Application Support/EMtrace/user_config.yaml
  Windows: %APPDATA%/EMtrace/user_config.yaml
  Linux:   $XDG_DATA_HOME/EMtrace/user_config.yaml  (default ~/.local/share)

This module has no Qt dependency and can be imported from core/, utils/, or gui/.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path


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


def user_config_file() -> Path:
    """Path to the per-user config YAML file."""
    return _data_dir() / "user_config.yaml"


def load_user_config() -> dict:
    """Load persisted user config.  Returns ``{}`` on missing or corrupt file."""
    cf = user_config_file()
    if cf.exists():
        try:
            import yaml  # PyYAML — always available (it's a core dependency)

            data = yaml.safe_load(cf.read_text(encoding="utf-8")) or {}
            return data if isinstance(data, dict) else {}
        except Exception:  # noqa: BLE001
            return {}
    return {}


def save_user_config(config: dict) -> None:
    """Persist *config* dict to disk as YAML."""
    import yaml

    user_config_file().write_text(
        yaml.dump(config, default_flow_style=False, allow_unicode=True),
        encoding="utf-8",
    )
