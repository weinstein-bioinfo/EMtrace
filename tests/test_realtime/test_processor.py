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


class TestTopTaxaNames:
    """Names for the live top-taxa list (GUI Realtime panel, `watch` console)."""

    @pytest.fixture
    def proc(self, tmp_dir: Path) -> BatchProcessor:
        return BatchProcessor(config=EMtraceConfig(), output_dir=tmp_dir / "rt", job_name="t")

    def test_snapshot_names_top_taxa(self):
        state = RealtimeState()
        state.F = {"11084": 0.9, "999": 0.1}
        state.taxa_names = {"11084": "Orthoflavivirus encephalitidis"}
        snap = state.get_snapshot()
        assert snap["top_taxa"] == [("11084", 0.9), ("999", 0.1)]
        assert snap["taxa_names"] == {"11084": "Orthoflavivirus encephalitidis", "999": ""}

    def test_names_come_from_the_abundance_table(self, proc):
        import pandas as pd
        abund = pd.DataFrame({
            "taxid": [11084, 1280, 42],
            "species": ["Orthoflavivirus encephalitidis", "NA", None],
            "organism": ["Tick-borne encephalitis virus", "Staphylococcus aureus", ""],
            "genus": ["Orthoflavivirus", "Staphylococcus", "Bacillus"],
        })
        proc._remember_names(abund)
        assert proc.state.taxa_names == {
            "11084": "Orthoflavivirus encephalitidis",   # species first
            "1280": "Staphylococcus aureus",             # organism when species is NA
            "42": "Bacillus",                            # then genus
        }

    def test_missing_names_are_looked_up_once(self, proc, monkeypatch):
        calls = []

        class FakeNCBI:
            def __init__(self, db_path):
                pass

            def get_organism(self, taxid):
                calls.append(taxid)
                return {11084: "Tick-borne encephalitis virus"}.get(taxid, "")

        monkeypatch.setattr("emtrace.models.taxonomy.NCBITaxonomy", FakeNCBI)
        db = proc.output_dir / "taxonomy.db"
        db.touch()
        proc.config.taxonomy.NCBI.db = str(db)
        proc.state.F = {"11084": 0.9, "7": 0.1}
        proc._name_top_taxa()
        assert proc.state.taxa_names == {"11084": "Tick-borne encephalitis virus", "7": ""}
        proc._name_top_taxa()
        assert sorted(calls) == [7, 11084]                # not again on the next batch

    def test_no_lookup_without_a_taxonomy_db(self, proc, tmp_dir):
        proc.config.taxonomy.NCBI.db = str(tmp_dir / "absent.db")
        proc.state.F = {"11084": 1.0}
        proc._name_top_taxa()
        assert proc.state.taxa_names == {}               # retried once a DB exists
