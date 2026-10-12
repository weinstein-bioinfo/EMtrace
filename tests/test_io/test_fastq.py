"""Tests for io/fastq.py — FASTQ reading and partitioning."""

from pathlib import Path

from emtrace.io.fastq import (
    count_reads,
    extract_fastq,
    iter_fastq,
    partition_fastq,
    read_id_from_header,
)


class TestIterFastq:
    """Test FASTQ iteration."""

    def test_plain_fastq(self, sample_fastq: Path):
        records = list(iter_fastq(sample_fastq))
        assert len(records) == 4
        assert records[0][0] == "@read1 length=20"
        assert len(records[0][1]) == 20

    def test_gzipped_fastq(self, sample_fastq_gz: Path):
        records = list(iter_fastq(sample_fastq_gz))
        assert len(records) == 4

    def test_empty_file(self, tmp_dir: Path):
        empty = tmp_dir / "empty.fq"
        empty.write_text("")
        assert list(iter_fastq(empty)) == []


class TestReadIdFromHeader:
    """Test read ID extraction."""

    def test_simple(self):
        assert read_id_from_header("@read1") == "read1"

    def test_with_description(self):
        assert read_id_from_header(
            "@read1 length=100",
        ) == "read1"

    def test_with_at_sign(self):
        # lstrip("@") strips all leading @s
        assert read_id_from_header(
            "@@complex_id desc",
        ) == "complex_id"


class TestCountReads:
    """Test read counting."""

    def test_plain(self, sample_fastq: Path):
        assert count_reads(sample_fastq) == 4

    def test_gzipped(self, sample_fastq_gz: Path):
        assert count_reads(sample_fastq_gz) == 4


class TestPartitionFastq:
    """Test FASTQ partitioning by ID set."""

    def test_partition(
        self, sample_fastq: Path, tmp_dir: Path,
    ):
        pass_fq = tmp_dir / "pass.fq.gz"
        flt_fq = tmp_dir / "flt.fq.gz"
        n_pass, n_flt = partition_fastq(
            sample_fastq, pass_fq, flt_fq,
            filter_ids={"read1", "read3"},
        )
        assert n_flt == 2
        assert n_pass == 2
        assert list(iter_fastq(flt_fq))
        assert list(iter_fastq(pass_fq))

    def test_all_pass(
        self, sample_fastq: Path, tmp_dir: Path,
    ):
        pass_fq = tmp_dir / "pass.fq.gz"
        flt_fq = tmp_dir / "flt.fq.gz"
        n_pass, n_flt = partition_fastq(
            sample_fastq, pass_fq, flt_fq,
            filter_ids=set(),
        )
        assert n_pass == 4
        assert n_flt == 0


class TestExtractFastq:
    """Test FASTQ extraction by ID."""

    def test_extract(
        self, sample_fastq: Path, tmp_dir: Path,
    ):
        out = tmp_dir / "extracted.fq.gz"
        n = extract_fastq(
            sample_fastq, out,
            keep_ids={"read2", "read4"},
        )
        assert n == 2
        records = list(iter_fastq(out))
        ids = {
            read_id_from_header(r[0]) for r in records
        }
        assert ids == {"read2", "read4"}
