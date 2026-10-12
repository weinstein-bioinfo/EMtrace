"""Tests for steps/stats.py — statistics step helpers."""

from pathlib import Path

import pandas as pd

from emtrace.steps.stats import (
    _division_counts,
    _fastq_stats,
    _read_lengths,
)


class TestFastqStats:
    """Test FASTQ statistics computation."""

    def test_basic_stats(self, sample_fastq: Path):
        stats = _fastq_stats(sample_fastq)
        assert stats["read_count"] == 4
        assert stats["min_len"] > 0
        assert stats["max_len"] >= stats["min_len"]
        assert stats["total_bp"] > 0
        assert stats["n50"] > 0

    def test_specific_values(self, sample_fastq: Path):
        stats = _fastq_stats(sample_fastq)
        # Reads: 20, 30, 15, 23 (from conftest)
        assert stats["read_count"] == 4
        assert stats["min_len"] == 15
        assert stats["max_len"] == 30

    def test_empty_fastq(self, tmp_dir: Path):
        empty = tmp_dir / "empty.fq"
        empty.write_text("")
        stats = _fastq_stats(empty)
        assert stats["read_count"] == 0
        assert stats["total_bp"] == 0
        assert stats["n50"] == 0

    def test_gzipped(self, sample_fastq_gz: Path):
        stats = _fastq_stats(sample_fastq_gz)
        assert stats["read_count"] == 4


class TestReadLengths:
    """Test read length extraction."""

    def test_lengths(self, sample_fastq: Path):
        lengths = _read_lengths(sample_fastq)
        assert len(lengths) == 4
        assert 15 in lengths
        assert 20 in lengths


class TestDivisionCounts:
    """Test division read counts from abundance table."""

    def test_counts(self, tmp_dir: Path):
        abund = tmp_dir / "abund_table.tsv"
        df = pd.DataFrame({
            "division": [
                "Viruses", "Viruses", "Bacteria",
                "Viruses",
            ],
            "readn": ["10", "20", "5", "15"],
        })
        df.to_csv(abund, sep="\t", index=False)

        counts = _division_counts(abund)
        assert len(counts) == 2
        # Viruses: 10+20+15=45
        virus_entry = next(
            c for c in counts
            if c["division"] == "Viruses"
        )
        assert virus_entry["read_count"] == 45

    def test_missing_file(self, tmp_dir: Path):
        counts = _division_counts(
            tmp_dir / "nonexistent.tsv",
        )
        assert counts == []
