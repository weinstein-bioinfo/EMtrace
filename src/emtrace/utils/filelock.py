"""Best-effort advisory file locking for cross-process coordination.

Used to serialize work that several processes might attempt at once against a
shared file — most importantly building a host-genome ``.mmi`` index in a shared
database directory when many HPC array tasks start together.

The lock is *advisory* and *best-effort*: it relies on ``fcntl.flock`` and
degrades gracefully on platforms / filesystems where that is unavailable
(``fcntl`` missing on Windows, or a network filesystem mounted without flock
support). In that case the context manager yields ``False`` so callers can fall
back to a still-correct (if less efficient) path rather than failing.
"""

from __future__ import annotations

import logging
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, Union

logger = logging.getLogger("emtrace.utils.filelock")

try:
    import fcntl  # POSIX only
    _HAVE_FCNTL = True
except ImportError:  # pragma: no cover - platform dependent
    fcntl = None  # type: ignore[assignment]
    _HAVE_FCNTL = False


@contextmanager
def file_lock(lock_path: Union[str, Path]) -> Iterator[bool]:
    """Acquire an exclusive advisory lock on *lock_path* for the block's duration.

    Yields ``True`` when the lock is held and ``False`` when locking is
    unsupported (so the caller can proceed without serialization). The lock is
    blocking: concurrent holders wait until the current holder exits the block.
    The lock file is created if missing and is left on disk (its presence is
    harmless; the lock state lives in the kernel, not the file contents).
    """
    lock_path = Path(lock_path)

    if not _HAVE_FCNTL:
        logger.debug("fcntl unavailable — proceeding without file lock (%s)", lock_path)
        yield False
        return

    fd = None
    try:
        lock_path.parent.mkdir(parents=True, exist_ok=True)
        fd = open(lock_path, "w")
        try:
            fcntl.flock(fd.fileno(), fcntl.LOCK_EX)
        except OSError as exc:
            # Filesystem doesn't support flock (e.g. some network mounts).
            logger.debug("flock unsupported on %s (%s) — proceeding unlocked", lock_path, exc)
            fd.close()
            fd = None
            yield False
            return
        yield True
    finally:
        if fd is not None:
            try:
                fcntl.flock(fd.fileno(), fcntl.LOCK_UN)
            except OSError:
                pass
            fd.close()
