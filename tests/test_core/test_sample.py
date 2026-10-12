"""Tests for core/sample.py — Sample state object."""

from pathlib import Path

from emtrace.core.pipelinejob import AnnotationRound, PipelineJob as Sample


class TestSample:
    """Test Sample dataclass."""

    def test_default_construction(self):
        s = Sample(name="s1")
        assert s.name == "s1"
        assert s.input_reads is None
        assert s.num_threads == 16

    def test_directory_properties(self, tmp_dir: Path):
        s = Sample(
            name="s1",
            output_dir=tmp_dir / "out",
        )
        assert s.reads_dir == tmp_dir / "out" / "reads"
        assert s.annotation_dir == tmp_dir / "out" / "annotation"
        assert s.assembly_dir == tmp_dir / "out" / "assembly"
        assert s.reports_dir == tmp_dir / "out" / "results"

    def test_init_dirs_creates_tree(self, tmp_dir: Path):
        s = Sample(
            name="s1",
            output_dir=tmp_dir / "out",
        )
        s.init_dirs()
        assert s.reads_dir.is_dir()
        assert s.annotation_dir.is_dir()
        assert s.assembly_dir.is_dir()
        assert s.reports_dir.is_dir()
        assert s.log_path.is_dir()

    def test_file_registry(self, tmp_dir: Path):
        s = Sample(name="s1")
        p = tmp_dir / "test.txt"
        p.write_text("hello")

        s.set_file("test", p)
        assert s.get_file("test") == p
        assert s.has_file("test") is True
        assert s.get_file("missing") is None
        assert s.has_file("missing") is False

    def test_current_reads_fallback(self, sample_fastq: Path):
        s = Sample(
            name="s1",
            input_reads=sample_fastq,
        )
        # Falls back to input_reads
        assert s.current_reads == sample_fastq

        # After setting current_reads
        new = Path("/tmp/new.fq")
        s.set_file("current_reads", new)
        assert s.current_reads == new


class TestAnnotationRound:
    """Test AnnotationRound."""

    def test_empty(self):
        r = AnnotationRound()
        assert r.databases == []

    def test_with_databases(self):
        r = AnnotationRound(
            databases=[
                {"db": "nt", "search": "minimap"},
            ],
        )
        assert len(r.databases) == 1
