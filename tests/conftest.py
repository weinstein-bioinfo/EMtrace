"""Pytest shared fixtures for emtrace tests."""

from __future__ import annotations

import gzip
import logging
import textwrap
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Generator

import pandas as pd
import pytest

from emtrace.core.config import EMtraceConfig
from emtrace.core.pipelinejob import PipelineJob as Sample


# ------------------------------------------------------------------
# Logging isolation
# ------------------------------------------------------------------

@pytest.fixture(autouse=True)
def _restore_logging():
    """Undo global logging changes a test may make.

    Invoking the CLI (`run`, `watch`) calls `setup_logging`, which sets
    ``logging.getLogger("emtrace").propagate = False``. That leaks into
    every later test in the session: `caplog` attaches to the root logger, so
    it silently captures nothing, and the failure surfaces in an unrelated
    file, ordered by collection. Snapshot and restore instead, so the suite is
    order-independent.

    Loggers that do not exist yet are the subtle part. Before the first
    `setup_logging`, "emtrace" is only a `PlaceHolder` in `loggerDict` --
    implicitly created as the parent of "emtrace.steps.report" -- with no
    state to save. `setup_logging` then promotes it to a real Logger. Anything
    promoted or created during a test is therefore reset to a pristine,
    propagating state rather than restored.
    """
    root = logging.getLogger()
    manager = root.manager
    saved_root = (root.level, list(root.handlers))
    saved_disable = manager.disable
    saved = {
        name: (lg.level, lg.propagate, lg.disabled, list(lg.handlers))
        for name, lg in manager.loggerDict.items()
        if isinstance(lg, logging.Logger)
    }
    try:
        yield
    finally:
        root.setLevel(saved_root[0])
        root.handlers[:] = saved_root[1]
        manager.disable = saved_disable
        for name, lg in list(manager.loggerDict.items()):
            if not isinstance(lg, logging.Logger):
                continue
            if name in saved:
                level, propagate, disabled, handlers = saved[name]
                lg.setLevel(level)
                lg.propagate = propagate
                lg.disabled = disabled
                lg.handlers[:] = handlers
            else:
                lg.setLevel(logging.NOTSET)
                lg.propagate = True
                lg.disabled = False
                lg.handlers.clear()


# ------------------------------------------------------------------
# Directories
# ------------------------------------------------------------------

@pytest.fixture
def tmp_dir() -> Generator[Path, None, None]:
    """Create a temporary directory for test files."""
    with TemporaryDirectory() as d:
        yield Path(d)


# ------------------------------------------------------------------
# FASTQ fixtures
# ------------------------------------------------------------------

_FASTQ_RECORDS = [
    ("@read1 length=20", "ATGCATGCATGCATGCATGC", "+", "IIIIIIIIIIIIIIIIIIII"),
    ("@read2 length=30", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGC", "+", "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII"),
    ("@read3 length=15", "AAACCCGGGTTTAAA", "+", "IIIIIIIIIIIIIII"),
    ("@read4 length=25", "TTTTCCCCCAAAAAGGGGGTTTTC", "+", "IIIIIIIIIIIIIIIIIIIIIIIII"),
]


@pytest.fixture
def sample_fastq(tmp_dir: Path) -> Path:
    """Create a small plain FASTQ file with 4 reads."""
    fq = tmp_dir / "reads.fq"
    lines = []
    for h, s, p, q in _FASTQ_RECORDS:
        lines.extend([h, s, p, q])
    fq.write_text("\n".join(lines) + "\n")
    return fq


@pytest.fixture
def sample_fastq_gz(tmp_dir: Path) -> Path:
    """Create a small gzipped FASTQ file with 4 reads."""
    fq = tmp_dir / "reads.fq.gz"
    lines = []
    for h, s, p, q in _FASTQ_RECORDS:
        lines.extend([h, s, p, q])
    with gzip.open(fq, "wt") as fh:
        fh.write("\n".join(lines) + "\n")
    return fq


# ------------------------------------------------------------------
# FASTA fixtures
# ------------------------------------------------------------------

@pytest.fixture
def sample_fasta(tmp_dir: Path) -> Path:
    """Create a small FASTA file."""
    fa = tmp_dir / "contigs.fa"
    fa.write_text(
        ">seq1 some description\n"
        "ATGCATGCATGCATGCATGC\n"
        ">seq2\n"
        "GCTAGCTAGCTAGCTAGCTA\n"
        ">seq3\n"
        "AAACCCGGGTTTAAA\n"
    )
    return fa


# ------------------------------------------------------------------
# Annotation table fixtures
# ------------------------------------------------------------------

@pytest.fixture
def sample_annotation_df() -> pd.DataFrame:
    """Create a small annotation DataFrame."""
    return pd.DataFrame({
        "search": ["minimap2"] * 5,
        "db": ["nt.abv"] * 5,
        "qseqid": ["read1", "read1", "read2", "read3", "read3"],
        "orf": ["NA"] * 5,
        "qseqlen": [200, 200, 300, 150, 150],
        "sseqid": ["ref1", "ref2", "ref3", "ref4", "ref5"],
        "bitscore": [500.0, 300.0, 450.0, 200.0, 180.0],
        "evalue": [0.0] * 5,
        "alen": [180, 150, 280, 140, 130],
        "pident": [95.0, 80.0, 90.0, 85.0, 82.0],
        "qlen": [190, 180, 290, 145, 140],
        "qcov": [0.95, 0.90, 0.97, 0.97, 0.93],
        "slen": [1000, 2000, 1500, 800, 900],
        "scov": [0.19, 0.09, 0.19, 0.18, 0.16],
        "staxid": ["10001", "10002", "10003", "10001", "10004"],
        "sname": ["Virus A", "Virus B", "Bact C", "Virus A", "Virus D"],
        "division": ["Viruses", "Viruses", "Bacteria", "Viruses", "Viruses"],
    })


@pytest.fixture
def sample_annotation_tsv(
    tmp_dir: Path, sample_annotation_df: pd.DataFrame,
) -> Path:
    """Write annotation DataFrame to a TSV file."""
    p = tmp_dir / "annot1.tsv"
    sample_annotation_df.to_csv(p, sep="\t", index=False)
    return p


# ------------------------------------------------------------------
# Config fixtures
# ------------------------------------------------------------------

@pytest.fixture
def sample_config(tmp_dir: Path) -> EMtraceConfig:
    """Create a minimal config for testing."""
    return EMtraceConfig.create_default()


@pytest.fixture
def sample_config_yaml(tmp_dir: Path) -> Path:
    """Write a minimal config YAML."""
    p = tmp_dir / "config.yaml"
    p.write_text(textwrap.dedent("""\
        runtime:
          num_threads: 4
        pipe: main
        pre: chopper+porechop
        hostgen: "0"
        ann1: minimap.nt.abv
        ann2: ""
        anns: "0"
        ass: refbased
        gen: mga
        annotation:
          blastn:
            bitscore: 200
          minimap:
            score: 200
          diamond_blastp:
            bitscore: 200
          orfipy:
            min_len: 300
        general.parameters:
          min_qcov: 0.40
        em:
          temp: 1.0
          max_iter: 500
          log_L_threshold: 0.01
        taxonomy:
          selected: NCBI
          NCBI:
            db: /tmp/taxonomy
          SILVA:
            db: /tmp/SILVA
        ann.databases: {}
        host.databases: {}
        ann.strategies: {}
    """))
    return p


# ------------------------------------------------------------------
# Sample fixtures
# ------------------------------------------------------------------

@pytest.fixture
def sample_obj(tmp_dir: Path, sample_fastq: Path) -> Sample:
    """Create a Sample with initialized directories."""
    s = Sample(
        name="test_sample",
        input_reads=sample_fastq,
        output_dir=tmp_dir / "output",
        num_threads=2,
    )
    s.init_dirs()
    return s
