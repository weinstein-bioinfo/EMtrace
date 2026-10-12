"""Tests for models/abundance.py — Abundance table schema and I/O."""

from pathlib import Path

import pandas as pd

from emtrace.models.abundance import (
    ABUNDANCE_FIELDS,
    AbundanceRow,
    empty_abundance_df,
    read_abundance_csv,
    read_abundance_tsv,
    write_abundance_csv,
    write_abundance_tsv,
)


class TestAbundanceRow:
    """Test Pydantic AbundanceRow model."""

    def test_valid_row(self):
        row = AbundanceRow(taxid="12345", readn=100)
        assert row.taxid == "12345"
        assert row.readn == 100
        assert row.prob == 0.0


class TestEmptyAbundanceDf:
    """Test empty abundance DataFrame."""

    def test_has_all_columns(self):
        df = empty_abundance_df()
        for col in ABUNDANCE_FIELDS:
            assert col in df.columns
        assert len(df) == 0


def _make_abundance_df():
    return pd.DataFrame({
        "taxid": ["10001", "10002"],
        "readn": [50, 30],
        "bitscore": [500.0, 300.0],
        "prob": [0.6, 0.4],
        "species": ["Virus A", "Virus, B"],
        "genus": ["GenusA", "GenusB"],
        "family": ["FamA", "FamB"],
        "order": ["OrdA", "OrdB"],
        "class": ["ClassA", "ClassB"],
        "division": ["Viruses", "Viruses"],
    })


class TestAbundanceTsvIO:
    """Test TSV roundtrip."""

    def test_roundtrip(self, tmp_dir: Path):
        df = _make_abundance_df()
        p = tmp_dir / "abund.tsv"
        write_abundance_tsv(df, p)
        df2 = read_abundance_tsv(p)
        assert len(df2) == 2
        assert "taxid" in df2.columns


class TestAbundanceCsvIO:
    """Test CSV roundtrip with quoting."""

    def test_roundtrip(self, tmp_dir: Path):
        df = _make_abundance_df()
        p = tmp_dir / "abund.csv"
        write_abundance_csv(df, p)
        df2 = read_abundance_csv(p)
        assert len(df2) == 2
        # Comma in species name should survive
        assert df2.iloc[1]["species"] == "Virus, B"
