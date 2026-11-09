"""`watch` applies `report.tail` the same way `run` does.

`_write_abundance_tsv` had no test at all before this; it is the half of the
tail trim a live run exercises. The trim rule itself belongs to ReportStep and
is covered in tests/test_steps/test_report_tail.py — what is pinned here is
that `watch` reaches it, passes `self.config.report`, and writes what it gets
back.

These tests deliberately point the taxonomy DB at a path that does not exist.
That is not an evasion: `watch` links lineage best-effort inside a try/except
and logs "Taxonomy linkage skipped" when the DB is unavailable, so a table
with no `division` column is a real and common shape on this path — and the
one where the trim must refuse rather than drop taxa blind. Using the real DB
instead would tie the test to $EMTRACE_DB_DIR and cost ~40 s.
"""

from __future__ import annotations

import logging

import pandas as pd
import pytest

from emtrace.core.config import EMtraceConfig
from emtrace.realtime.processor import BatchProcessor


@pytest.fixture
def cfg(tmp_path):
    c = EMtraceConfig.create_default()
    c.taxonomy.NCBI.db = str(tmp_path / "no-such-taxonomy.db")
    return c


def processor(cfg, job_dir):
    return BatchProcessor(config=cfg, job_dir=job_dir, job_name="rt")


# Six reads over three taxa: taxid 1 dominates, 2 and 3 are single-read tails.
# Columns follow the annotation schema _build_abundance_table consumes.
ANNOT = pd.DataFrame({
    "search": ["Minimap2"] * 6,
    "db": ["db"] * 6,
    "qseqid": ["r1", "r2", "r3", "r4", "r5", "r6"],
    "orf": ["NA"] * 6,
    "staxid": ["1", "1", "1", "1", "2", "3"],
    "bitscore": [100.0] * 6,
    "qcov": [1.0] * 6,
    "species": ["Big bug"] * 4 + ["Rare virus", "Rare bug"],
})
F = {"1": 0.98, "2": 0.01, "3": 0.01}


def _write(proc):
    """Drive _write_abundance_tsv and return the table it wrote."""
    proc.state.annotations = ANNOT.copy()
    proc.state.F = dict(F)
    proc._write_em_tsv(proc.state.F, proc._em_prob_path, n_queries=len(ANNOT))
    proc._write_abundance_tsv()
    out = proc.output_dir / "results" / "abund_table.tsv"
    return pd.read_csv(out, sep="\t", dtype={"taxid": str}) if out.exists() else None


def test_watch_writes_an_abundance_table(cfg, tmp_path):
    got = _write(processor(cfg, tmp_path / "w"))
    assert got is not None
    assert set(got["taxid"]) == {"1", "2", "3"}


def test_watch_does_not_trim_by_default(cfg, tmp_path):
    assert cfg.report.tail == 0.0
    assert len(_write(processor(cfg, tmp_path / "w"))) == 3


def test_watch_refuses_to_trim_without_a_division_column(cfg, tmp_path, caplog):
    # The realistic `watch` shape when taxonomy linkage was skipped. Reaching
    # this warning proves _write_abundance_tsv calls ReportStep._apply_tail_trim
    # with the live config: the message is that helper's.
    cfg.report.tail = 0.01
    with caplog.at_level(logging.WARNING):
        got = _write(processor(cfg, tmp_path / "w"))
    assert len(got) == 3, "refusal must leave every taxon in place"
    assert "NOT trimming" in caplog.text
    assert "report.tail=0.01" in caplog.text


def test_watch_reads_tail_from_the_live_config(cfg, tmp_path, caplog):
    # The config must be the one the processor was built with, not a default.
    cfg.report.tail = 0.05
    with caplog.at_level(logging.WARNING):
        _write(processor(cfg, tmp_path / "w"))
    assert "report.tail=0.05" in caplog.text
