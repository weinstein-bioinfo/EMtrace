"""Directory watcher for MinKNOW output.

Monitors a directory for new ``.fastq`` / ``.fastq.gz`` files using
the ``watchdog`` library.  Debounces file detection (MinKNOW writes
incrementally) and pushes stable files into a thread-safe queue for
the batch processor.

Files count whether they are written, copied or moved/renamed into the
directory (Finder drag, ``mv``, rsync's temp-then-rename). A periodic
rescan catches files whose events were missed (network shares, coalesced
FSEvents). A file already processed is recognised after a rename by its
identity (device, inode, size, mtime) and not processed twice.

Usage::

    from emtrace.realtime.watcher import FastqWatcher

    watcher = FastqWatcher("/path/to/minknow/output")
    watcher.start()
    while True:
        path = watcher.get()  # blocks until a file is ready
        process(path)
    watcher.stop()
"""

from __future__ import annotations

import logging
import threading
import time
from typing import cast
from pathlib import Path
from queue import Empty, Queue
from typing import Optional, Set

logger = logging.getLogger("emtrace.realtime.watcher")

# Glob patterns for FASTQ files
_FASTQ_SUFFIXES = {".fastq", ".fq"}
_GZIP_EXT = ".gz"


def _is_fastq(path: Path) -> bool:
    """Check if a path looks like a FASTQ file."""
    name = path.name
    if name.endswith(_GZIP_EXT):
        stem = Path(name[: -len(_GZIP_EXT)])
        return stem.suffix in _FASTQ_SUFFIXES
    return path.suffix in _FASTQ_SUFFIXES


class FastqWatcher:
    """Watch a directory for new FASTQ files.

    Parameters
    ----------
    watch_dir : Path or str
        Directory to monitor.
    debounce_seconds : float
        Wait this long after a file stops changing before
        enqueuing it (MinKNOW writes incrementally).
    recursive : bool
        Watch subdirectories as well.
    poll_interval : float
        Fallback polling interval in seconds (used when
        watchdog is not available or for the debounce loop).
    scan_existing : bool
        If True, enqueue FASTQ files already present in the
        directory at startup. If False they are ignored, also by
        the rescan and the polling fallback.
    rescan_interval : float
        With watchdog, also rescan the directory every this many
        seconds for files whose events were missed; 0 = off.
    """

    def __init__(
        self,
        watch_dir: str | Path,
        *,
        debounce_seconds: float = 2.0,
        recursive: bool = True,
        poll_interval: float = 1.0,
        scan_existing: bool = False,
        rescan_interval: float = 30.0,
    ) -> None:
        self.watch_dir = Path(watch_dir)
        self.debounce_seconds = debounce_seconds
        self.recursive = recursive
        self.poll_interval = poll_interval
        self.scan_existing = scan_existing
        self.rescan_interval = rescan_interval

        self._queue: Queue[Path] = Queue()
        self._seen: Set[Path] = set()
        # Identities of processed files: a renamed file keeps its identity
        self._seen_ids: Set[tuple] = set()
        # Files present at start with scan_existing=False: skipped by scans
        self._preexisting: Set[Path] = set()
        self._watchdog_active = False
        # Size at the last scan, so scans restart the debounce only on growth
        self._scanned_size: dict[Path, int] = {}
        self._stop_event = threading.Event()
        self._observer = None
        self._debounce_thread: Optional[threading.Thread] = None

        # Pending files: path → last-modified-time
        self._pending: dict[Path, float] = {}
        self._pending_lock = threading.Lock()

    # ----------------------------------------------------------
    # Public API
    # ----------------------------------------------------------

    def start(self) -> None:
        """Begin watching. Non-blocking."""
        if not self.watch_dir.is_dir():
            raise FileNotFoundError(
                f"Watch directory does not exist: "
                f"{self.watch_dir}"
            )

        logger.info(
            "Watching %s for FASTQ files "
            "(debounce=%.1fs, recursive=%s)",
            self.watch_dir,
            self.debounce_seconds,
            self.recursive,
        )

        # Scan existing files if requested; otherwise keep scans off them
        if self.scan_existing:
            self._scan_existing()
        else:
            self._preexisting = set(self._list_fastq())

        # Start watchdog observer
        self._start_observer()

        # Start debounce thread
        self._stop_event.clear()
        self._debounce_thread = threading.Thread(
            target=self._debounce_loop,
            name="fastq-debounce",
            daemon=True,
        )
        self._debounce_thread.start()

    def stop(self) -> None:
        """Stop watching and clean up."""
        self._stop_event.set()
        if self._observer is not None:
            self._observer.stop()
            self._observer.join(timeout=5)
            self._observer = None
        if self._debounce_thread is not None:
            self._debounce_thread.join(timeout=5)
            self._debounce_thread = None
        logger.info("Watcher stopped")

    def get(
        self, timeout: Optional[float] = None,
    ) -> Optional[Path]:
        """Get the next stable FASTQ file path.

        Parameters
        ----------
        timeout : float, optional
            Seconds to wait. None = block forever.

        Returns
        -------
        Path or None
            Path to a FASTQ file, or None on timeout.
        """
        try:
            return self._queue.get(timeout=timeout)
        except Empty:
            return None

    def get_nowait(self) -> Optional[Path]:
        """Non-blocking get."""
        try:
            return self._queue.get_nowait()
        except Empty:
            return None

    @property
    def pending_count(self) -> int:
        """Number of files waiting for debounce."""
        with self._pending_lock:
            return len(self._pending)

    @property
    def queued_count(self) -> int:
        """Number of files ready in the queue."""
        return self._queue.qsize()

    @property
    def seen_count(self) -> int:
        """Total files seen (queued + processed)."""
        return len(self._seen)

    @property
    def is_running(self) -> bool:
        """Whether the watcher is active."""
        return not self._stop_event.is_set()

    # ----------------------------------------------------------
    # Internal: watchdog observer
    # ----------------------------------------------------------

    def _start_observer(self) -> None:
        """Start the watchdog filesystem observer."""
        try:
            from watchdog.events import (
                FileSystemEventHandler,
            )
            from watchdog.observers import Observer
        except ImportError:
            logger.warning(
                "watchdog not installed — falling back "
                "to polling mode "
                "(pip install watchdog)"
            )
            self._start_poller()
            return

        watcher = self

        class _Handler(FileSystemEventHandler):
            def on_created(self, event):
                if not event.is_directory:
                    watcher._on_file_event(
                        Path(cast(str,event.src_path)),
                    )

            def on_modified(self, event):
                if not event.is_directory:
                    watcher._on_file_event(
                        Path(cast(str,event.src_path)),
                    )

            def on_moved(self, event):
                if not event.is_directory:
                    watcher._on_file_moved(
                        Path(cast(str,event.src_path)),
                        Path(cast(str,event.dest_path)),
                    )

        self._observer = Observer()
        self._observer.schedule(
            _Handler(),
            str(self.watch_dir),
            recursive=self.recursive,
        )
        self._observer.daemon = True
        self._observer.start()
        self._watchdog_active = True

    def _start_poller(self) -> None:
        """Fallback: poll directory for new files."""
        t = threading.Thread(
            target=self._poll_loop,
            name="fastq-poller",
            daemon=True,
        )
        t.start()

    def _poll_loop(self) -> None:
        """Poll for new FASTQ files (fallback mode)."""
        while not self._stop_event.is_set():
            self._scan_dir()
            self._stop_event.wait(self.poll_interval)

    def _list_fastq(self) -> list[Path]:
        """FASTQ files now in the watched directory (resolved paths)."""
        pattern = "**/*" if self.recursive else "*"
        return [
            p.resolve() for p in self.watch_dir.glob(pattern)
            if p.is_file() and _is_fastq(p)
        ]

    def _scan_dir(self) -> None:
        """Scan directory for new FASTQ files (poller, rescan)."""
        for p in self._list_fastq():
            if p in self._preexisting or p in self._seen:
                continue
            try:
                size = p.stat().st_size
            except OSError:
                continue
            # Restart the debounce only when the file grew since the last scan
            if self._scanned_size.get(p) != size:
                self._scanned_size[p] = size
                self._on_file_event(p)

    def _scan_existing(self) -> None:
        """Enqueue FASTQ files already present."""
        pattern = "**/*" if self.recursive else "*"
        count = 0
        for p in sorted(self.watch_dir.glob(pattern)):
            if p.is_file() and _is_fastq(p):
                self._on_file_event(p)
                count += 1
        if count:
            logger.info(
                "Found %d existing FASTQ files", count,
            )

    # ----------------------------------------------------------
    # Internal: debouncing
    # ----------------------------------------------------------

    def _on_file_event(self, path: Path) -> None:
        """Called when a FASTQ file is created/modified."""
        path = path.resolve()
        if path in self._seen:
            return
        if not _is_fastq(path):
            return
        with self._pending_lock:
            self._pending[path] = time.monotonic()

    def _on_file_moved(self, src: Path, dest: Path) -> None:
        """A file was moved or renamed: track it under its new path."""
        with self._pending_lock:
            self._pending.pop(src.resolve(), None)
        if self._in_scope(dest):
            self._on_file_event(dest)

    def _in_scope(self, path: Path) -> bool:
        """True if *path* lies in the watched directory (or below, if recursive)."""
        path, root = path.resolve(), self.watch_dir.resolve()
        return path.is_relative_to(root) if self.recursive else path.parent == root

    @staticmethod
    def _file_id(path: Path) -> Optional[tuple]:
        """Identity that survives a rename: (device, inode, size, mtime)."""
        try:
            st = path.stat()
        except OSError:
            return None
        return (st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns)

    def _debounce_loop(self) -> None:
        """Promote pending files to queue once stable."""
        last_rescan = time.monotonic()
        while not self._stop_event.is_set():
            now = time.monotonic()
            ready = []

            # Safety net for missed events (watchdog mode; the poller scans anyway)
            if (self._watchdog_active and self.rescan_interval > 0
                    and now - last_rescan >= self.rescan_interval):
                last_rescan = now
                self._scan_dir()

            with self._pending_lock:
                for path, last_mod in list(
                    self._pending.items()
                ):
                    age = now - last_mod
                    if age >= self.debounce_seconds:
                        # Verify file still exists and is
                        # not actively being written
                        if self._file_is_stable(path):
                            ready.append(path)
                            del self._pending[path]

            for path in ready:
                if path in self._seen:
                    continue
                self._seen.add(path)
                fid = self._file_id(path)
                if fid is not None and fid in self._seen_ids:
                    logger.info(
                        "Skipping %s: renamed copy of a file "
                        "already processed", path.name,
                    )
                    continue
                if fid is not None:
                    self._seen_ids.add(fid)
                self._queue.put(path)
                logger.info(
                    "New FASTQ: %s", path.name,
                )

            self._stop_event.wait(
                min(self.poll_interval, 0.5),
            )

    @staticmethod
    def _file_is_stable(path: Path) -> bool:
        """Check that a file exists and is not zero-length."""
        try:
            st = path.stat()
            return st.st_size > 0
        except OSError:
            return False
