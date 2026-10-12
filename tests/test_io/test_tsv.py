"""Tests for io/tsv.py — TSV I/O and filter_tophits."""

from pathlib import Path

import pandas as pd

from emtrace.io.tsv import (
    filter_tophits,
    read_tsv,
    write_tsv,
)


class TestReadWriteTsv:
    """Test TSV read/write roundtrip."""

    def test_roundtrip(self, tmp_dir: Path):
        df = pd.DataFrame({
            "a": ["x", "y", "z"],
            "b": [1, 2, 3],
        })
        p = tmp_dir / "test.tsv"
        write_tsv(df, p)
        df2 = read_tsv(p)
        assert len(df2) == 3
        assert list(df2.columns) == ["a", "b"]


class TestFilterTophits:
    """Test top-hit filtering."""

    def test_single_top_hit(self):
        df = pd.DataFrame({
            "qseqid": ["r1", "r1", "r2", "r2"],
            "bitscore": [100, 200, 300, 150],
        })
        result = filter_tophits(df)
        assert len(result) == 2
        # r1's top hit = 200
        r1 = result[result["qseqid"] == "r1"]
        assert r1.iloc[0]["bitscore"] == 200
        # r2's top hit = 300
        r2 = result[result["qseqid"] == "r2"]
        assert r2.iloc[0]["bitscore"] == 300

    def test_retain_ties(self):
        df = pd.DataFrame({
            "qseqid": ["r1", "r1", "r1"],
            "bitscore": [100, 200, 200],
        })
        result = filter_tophits(df, retain_ties=True)
        assert len(result) == 2  # both 200s kept

    def test_unsorted_input(self):
        """filter_tophits should handle unsorted input."""
        df = pd.DataFrame({
            "qseqid": ["r2", "r1", "r2", "r1"],
            "bitscore": [50, 300, 400, 100],
        })
        result = filter_tophits(df)
        assert len(result) == 2
        r1 = result[result["qseqid"] == "r1"]
        assert r1.iloc[0]["bitscore"] == 300

    def test_string_bitscores(self):
        """Bitscores as strings should be coerced."""
        df = pd.DataFrame({
            "qseqid": ["r1", "r1"],
            "bitscore": ["100", "200"],
        })
        result = filter_tophits(df)
        assert len(result) == 1
        assert result.iloc[0]["bitscore"] == 200

    def test_empty_dataframe(self):
        df = pd.DataFrame(
            columns=["qseqid", "bitscore"],
        )
        result = filter_tophits(df)
        assert len(result) == 0
