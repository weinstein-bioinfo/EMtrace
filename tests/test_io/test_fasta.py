"""Tests for io/fasta.py — FASTA I/O utilities."""

from pathlib import Path

from emtrace.io.fasta import (
    exclude_fasta_by_ids,
    extract_fasta_by_ids,
    fasta_id_from_header,
    fasta_lengths,
    fastq_to_fasta,
    iter_fasta,
    reads_info_tsv,
)


class TestIterFasta:
    """Test FASTA iteration."""

    def test_read_fasta(self, sample_fasta: Path):
        records = list(iter_fasta(sample_fasta))
        assert len(records) == 3
        assert records[0][0] == ">seq1 some description"
        assert records[0][1] == "ATGCATGCATGCATGCATGC"


class TestFastaIdFromHeader:
    """Test FASTA header ID extraction."""

    def test_simple(self):
        assert fasta_id_from_header(">seq1") == "seq1"

    def test_with_description(self):
        assert fasta_id_from_header(
            ">seq1 some description",
        ) == "seq1"


class TestFastqToFasta:
    """Test FASTQ to FASTA conversion."""

    def test_conversion(
        self, sample_fastq: Path, tmp_dir: Path,
    ):
        fa = fastq_to_fasta(
            sample_fastq, tmp_dir / "out.fa",
        )
        records = list(iter_fasta(fa))
        assert len(records) == 4
        # Headers should start with >
        assert records[0][0].startswith(">")
        assert "read1" in records[0][0]


class TestReadsInfoTsv:
    """Test reads info TSV generation."""

    def test_info_tsv(
        self, sample_fastq: Path, tmp_dir: Path,
    ):
        import pandas as pd
        out = reads_info_tsv(
            sample_fastq, tmp_dir / "info.tsv",
        )
        df = pd.read_csv(out, sep="\t")
        assert len(df) == 4
        assert "seqid" in df.columns
        assert "length" in df.columns
        assert df.loc[0, "seqid"] == "read1"
        assert df.loc[0, "length"] == 20


class TestExtractFasta:
    """Test FASTA extraction by ID."""

    def test_extract(
        self, sample_fasta: Path, tmp_dir: Path,
    ):
        out = tmp_dir / "extracted.fa"
        n = extract_fasta_by_ids(
            sample_fasta, out, keep_ids={"seq1", "seq3"},
        )
        assert n == 2
        records = list(iter_fasta(out))
        ids = {
            fasta_id_from_header(h) for h, _ in records
        }
        assert ids == {"seq1", "seq3"}

    def test_exclude(
        self, sample_fasta: Path, tmp_dir: Path,
    ):
        out = tmp_dir / "excluded.fa"
        n = exclude_fasta_by_ids(
            sample_fasta, out,
            exclude_ids={"seq2"},
        )
        assert n == 2


class TestFastaLengths:
    """Test FASTA length extraction."""

    def test_lengths(self, sample_fasta: Path):
        lengths = fasta_lengths(sample_fasta)
        assert lengths["seq1"] == 20
        assert lengths["seq2"] == 20
        assert lengths["seq3"] == 15
