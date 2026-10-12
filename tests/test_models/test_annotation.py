"""Tests for models/annotation.py — Annotation table schema and I/O."""

from pathlib import Path

import pandas as pd

from emtrace.models.annotation import (
    ANNOT_FIELDS,
    AnnotationRow,
    empty_annotation_df,
    read_annotation_csv,
    read_annotation_tsv,
    validate_annotation_df,
    write_annotation_csv,
    write_annotation_tsv,
)


class TestAnnotationRow:
    """Test Pydantic AnnotationRow model."""

    def test_valid_row(self):
        row = AnnotationRow(
            search="minimap2",
            db="nt",
            qseqid="read1",
            sseqid="ref1",
            bitscore=500.0,
        )
        assert row.search == "minimap2"
        assert row.orf == "NA"
        assert row.division == ""

    def test_evalue_validation(self):
        import pytest
        with pytest.raises(Exception):
            AnnotationRow(
                search="x", db="x", qseqid="x",
                sseqid="x", bitscore=1.0,
                evalue=-1.0,
            )


class TestEmptyAnnotationDf:
    """Test empty DataFrame creation."""

    def test_has_all_columns(self):
        df = empty_annotation_df()
        for col in ANNOT_FIELDS:
            assert col in df.columns
        assert len(df) == 0


class TestAnnotationTsvIO:
    """Test TSV read/write roundtrip."""

    def test_roundtrip(
        self, tmp_dir: Path, sample_annotation_df,
    ):
        p = tmp_dir / "annot.tsv"
        write_annotation_tsv(sample_annotation_df, p)
        df2 = read_annotation_tsv(p)
        assert len(df2) == len(sample_annotation_df)
        assert "qseqid" in df2.columns
        assert df2["bitscore"].dtype == float


class TestAnnotationCsvIO:
    """Test CSV read/write with quoting."""

    def test_roundtrip(
        self, tmp_dir: Path, sample_annotation_df,
    ):
        p = tmp_dir / "annot.csv"
        write_annotation_csv(sample_annotation_df, p)
        df2 = read_annotation_csv(p)
        assert len(df2) == len(sample_annotation_df)

    def test_special_chars_preserved(self, tmp_dir: Path):
        """Commas in sname should survive CSV roundtrip."""
        df = empty_annotation_df()
        row = {
            "search": "minimap2", "db": "nt",
            "qseqid": "r1", "sseqid": "s1",
            "bitscore": 100.0,
            "sname": "Virus, subtype A",
            "division": "Viruses",
        }
        df = pd.concat(
            [df, pd.DataFrame([row])],
            ignore_index=True,
        )
        p = tmp_dir / "annot.csv"
        write_annotation_csv(df, p)
        df2 = read_annotation_csv(p)
        assert df2.iloc[0]["sname"] == "Virus, subtype A"


class TestValidateAnnotationDf:
    """Test annotation DataFrame validation."""

    def test_valid(self, sample_annotation_df):
        result = validate_annotation_df(sample_annotation_df)
        assert len(result) == 5

    def test_missing_column_raises(self):
        df = pd.DataFrame({"qseqid": ["r1"]})
        import pytest
        with pytest.raises(ValueError, match="missing columns"):
            validate_annotation_df(df)
