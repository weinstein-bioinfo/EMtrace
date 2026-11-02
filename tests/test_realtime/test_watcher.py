"""Tests for realtime/watcher.py — directory watcher."""

import time
from pathlib import Path

import pytest

from emtrace.realtime.watcher import (
    FastqWatcher,
    _is_fastq,
)


class TestIsFastq:
    """Test FASTQ file detection."""

    def test_plain_fastq(self):
        assert _is_fastq(Path("reads.fastq")) is True
        assert _is_fastq(Path("reads.fq")) is True

    def test_gzipped_fastq(self):
        assert _is_fastq(Path("reads.fastq.gz")) is True
        assert _is_fastq(Path("reads.fq.gz")) is True

    def test_non_fastq(self):
        assert _is_fastq(Path("reads.fa")) is False
        assert _is_fastq(Path("reads.bam")) is False
        assert _is_fastq(Path("reads.txt")) is False


class TestFastqWatcher:
    """Test FastqWatcher lifecycle."""

    def test_missing_dir_raises(self, tmp_dir: Path):
        w = FastqWatcher(tmp_dir / "nonexistent")
        with pytest.raises(FileNotFoundError):
            w.start()

    def test_start_stop(self, tmp_dir: Path):
        w = FastqWatcher(
            tmp_dir, debounce_seconds=0.1,
        )
        w.start()
        assert w.is_running
        w.stop()
        assert not w.is_running

    def test_scan_existing(self, tmp_dir: Path):
        # Create a FASTQ file before starting
        fq = tmp_dir / "existing.fq"
        fq.write_text(
            "@r1\nATGC\n+\nIIII\n",
        )

        w = FastqWatcher(
            tmp_dir,
            debounce_seconds=0.1,
            poll_interval=0.1,
            scan_existing=True,
        )
        w.start()
        # Use blocking get with timeout
        path = w.get(timeout=3.0)
        w.stop()

        assert path is not None
        assert path.name == "existing.fq"

    def test_detects_new_file(self, tmp_dir: Path):
        w = FastqWatcher(
            tmp_dir,
            debounce_seconds=0.2,
            poll_interval=0.1,
        )
        w.start()

        # Create a new file
        fq = tmp_dir / "new.fastq"
        fq.write_text("@r1\nATGC\n+\nIIII\n")

        # Wait for debounce
        path = w.get(timeout=3.0)
        w.stop()

        assert path is not None
        assert path.name == "new.fastq"

    def test_ignores_non_fastq(self, tmp_dir: Path):
        w = FastqWatcher(
            tmp_dir,
            debounce_seconds=0.1,
            poll_interval=0.1,
        )
        w.start()

        # Create a non-FASTQ file
        (tmp_dir / "data.txt").write_text("hello")
        time.sleep(0.5)

        path = w.get_nowait()
        w.stop()
        assert path is None

    def test_get_timeout(self, tmp_dir: Path):
        w = FastqWatcher(tmp_dir, debounce_seconds=0.1)
        w.start()
        path = w.get(timeout=0.2)
        w.stop()
        assert path is None

    def test_seen_count(self, tmp_dir: Path):
        fq = tmp_dir / "test.fq"
        fq.write_text("@r1\nATGC\n+\nIIII\n")

        w = FastqWatcher(
            tmp_dir,
            debounce_seconds=0.1,
            poll_interval=0.1,
            scan_existing=True,
        )
        w.start()
        # Use blocking get to ensure file is processed
        w.get(timeout=3.0)
        w.stop()
        assert w.seen_count >= 1


class TestMovedAndMissedFiles:
    """Files moved/renamed in, missed events, and renames of processed files."""

    FQ = "@r1\nATGC\n+\nIIII\n"

    def _watcher(self, d: Path, **kw) -> FastqWatcher:
        return FastqWatcher(d, debounce_seconds=0.2, poll_interval=0.1, **kw)

    def test_detects_file_moved_in(self, tmp_dir: Path):
        watched, outside = tmp_dir / "watch", tmp_dir / "outside"
        watched.mkdir(); outside.mkdir()
        src = outside / "moved.fastq"
        src.write_text(self.FQ)
        w = self._watcher(watched)
        w.start()
        src.rename(watched / "moved.fastq")              # Finder drag / mv
        path = w.get(timeout=5.0)
        w.stop()
        assert path is not None and path.name == "moved.fastq"

    def test_detects_temp_then_rename(self, tmp_dir: Path):
        w = self._watcher(tmp_dir)
        w.start()
        tmp = tmp_dir / ".reads.fastq.gz.XYZ123"           # rsync-style temp name
        tmp.write_text(self.FQ)
        tmp.rename(tmp_dir / "reads.fastq")
        path = w.get(timeout=5.0)
        w.stop()
        assert path is not None and path.name == "reads.fastq"

    def test_renamed_processed_file_is_not_processed_again(self, tmp_dir: Path):
        w = self._watcher(tmp_dir)
        w.start()
        fq = tmp_dir / "a.fastq"
        fq.write_text(self.FQ)
        assert w.get(timeout=5.0) is not None
        fq.rename(tmp_dir / "a_old.fastq")
        time.sleep(1.0)
        assert w.get_nowait() is None
        w.stop()

    def test_on_moved_out_of_scope_is_ignored(self, tmp_dir: Path):
        w = self._watcher(tmp_dir / "watch")
        (tmp_dir / "watch").mkdir()
        w._on_file_moved(tmp_dir / "watch" / "x.fastq", tmp_dir / "elsewhere" / "x.fastq")
        assert w._pending == {}

    def test_rescan_catches_missed_events(self, tmp_dir: Path):
        w = self._watcher(tmp_dir, rescan_interval=0.3)
        w.start()
        w._watchdog_active = True                          # as if watchdog ran
        if w._observer is not None:                        # ...but drops every event
            w._observer.stop(); w._observer.join(timeout=5); w._observer = None
        (tmp_dir / "missed.fastq").write_text(self.FQ)
        path = w.get(timeout=5.0)
        w.stop()
        assert path is not None and path.name == "missed.fastq"

    def test_preexisting_files_stay_ignored_by_scans(self, tmp_dir: Path):
        (tmp_dir / "old.fastq").write_text(self.FQ)
        w = self._watcher(tmp_dir, rescan_interval=0.2)
        w._start_observer = w._start_poller                # polling scans every 0.1 s
        w.start()
        time.sleep(0.8)
        assert w.get_nowait() is None
        w.stop()

    def test_poller_promotes_a_file_despite_frequent_scans(self, tmp_dir: Path):
        # Scans restart the debounce only when the file grew, so poll < debounce works
        w = FastqWatcher(tmp_dir, debounce_seconds=0.5, poll_interval=0.1)
        w._start_observer = w._start_poller
        w.start()
        (tmp_dir / "p.fastq").write_text(self.FQ)
        path = w.get(timeout=5.0)
        w.stop()
        assert path is not None and path.name == "p.fastq"
