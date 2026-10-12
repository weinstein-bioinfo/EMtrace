"""Tests for the best-effort advisory file lock."""

from __future__ import annotations

from emtrace.utils import filelock
from emtrace.utils.filelock import file_lock


def test_file_lock_acquires(tmp_path):
    lock = tmp_path / "x.lock"
    with file_lock(lock) as held:
        assert held is True
    # Re-acquirable after release.
    with file_lock(lock) as held:
        assert held is True


def test_file_lock_falls_back_without_fcntl(tmp_path, monkeypatch):
    """When fcntl is unavailable the lock yields False (caller stays correct)."""
    monkeypatch.setattr(filelock, "_HAVE_FCNTL", False)
    with file_lock(tmp_path / "y.lock") as held:
        assert held is False
