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
