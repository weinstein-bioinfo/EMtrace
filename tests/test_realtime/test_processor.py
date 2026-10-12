"""Tests for realtime/processor.py — batch processor."""

from pathlib import Path

import pytest

from emtrace.core.config import EMtraceConfig
from emtrace.realtime.processor import (
    BatchProcessor,
    RealtimeState,
)


class TestRealtimeState:
    """Test thread-safe state accumulation."""

    def test_initial_state(self):
        state = RealtimeState()
        assert state.batch_count == 0
        assert state.total_reads == 0
        assert state.total_files == 0
        assert state.F == {}

    def test_snapshot(self):
        state = RealtimeState()
        snap = state.get_snapshot()
        assert snap["batch_count"] == 0
        assert snap["total_reads"] == 0
        assert snap["n_taxa"] == 0
        assert snap["top_taxa"] == []

    def test_update(self):
        import pandas as pd
        state = RealtimeState()
        df = pd.DataFrame({
            "qseqid": ["r1"],
            "staxid": ["A"],
            "bitscore": [100.0],
        })
        state.update(df, None, 10, 1)
        assert state.batch_count == 1
        assert state.total_reads == 10
        assert state.total_files == 1
        assert len(state.annotations) == 1


class TestBatchProcessorInit:
    """Test BatchProcessor construction."""

    def test_creates_directories(self, tmp_dir: Path):
        cfg = EMtraceConfig()
        out = tmp_dir / "rt_output"
        proc = BatchProcessor(
            config=cfg,
            output_dir=out,
            job_name="test",
        )
        assert (out / "batches").is_dir()
        assert (out / "results").is_dir()

    def test_concat_fastq(self, tmp_dir: Path):
        """Test FASTQ concatenation helper."""
        # Create two small FASTQs
        fq1 = tmp_dir / "a.fq"
        fq2 = tmp_dir / "b.fq"
        fq1.write_text("@r1\nATGC\n+\nIIII\n")
        fq2.write_text("@r2\nGCTA\n+\nIIII\n")

        out = tmp_dir / "merged.fq.gz"
        BatchProcessor._concat_fastq([fq1, fq2], out)
        assert out.exists()

        # Read back
        import gzip
        with gzip.open(out, "rt") as fh:
            content = fh.read()
        assert "@r1" in content
        assert "@r2" in content

    def test_count_reads(self, sample_fastq: Path):
        n = BatchProcessor._count_reads(sample_fastq)
        assert n == 4

    def test_write_em_tsv(self, tmp_dir: Path):
        import pandas as pd
        p = tmp_dir / "em.tsv"
        BatchProcessor._write_em_tsv(
            {"A": 0.6, "B": 0.4}, p,
        )
        assert p.exists()
        df = pd.read_csv(p, sep="\t")
        assert len(df) == 2
        assert df.iloc[0]["taxid"] == "A"
        assert df.iloc[0]["prob"] == 0.6
