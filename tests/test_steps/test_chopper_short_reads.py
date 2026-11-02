"""Chopper step: warning when many reads are shorter than chopper.minlength."""

import logging
from pathlib import Path

import pytest

from emtrace.core.config import EMtraceConfig
from emtrace.steps.chopper import ChopperStep


def _fastq(path: Path, lengths: list[int]) -> Path:
    with open(path, "w") as fh:
        for i, n in enumerate(lengths):
            fh.write(f"@r{i}\n{'A' * n}\n+\n{'I' * n}\n")
    return path


@pytest.fixture
def step():
    return ChopperStep(config=EMtraceConfig())


def test_warns_when_quarter_of_reads_are_short(step, tmp_dir, caplog):
    # 4 of 10 reads below 400 bp (40 %)
    fq = _fastq(tmp_dir / "r.fq", [200] * 4 + [800] * 6)
    with caplog.at_level(logging.WARNING):
        step._warn_short_reads(fq, 400)
    assert "40% of reads < minlength 400 bp (median 800 bp)" in caplog.text


def test_silent_below_threshold(step, tmp_dir, caplog):
    # 2 of 10 reads below 400 bp (20 %)
    fq = _fastq(tmp_dir / "r.fq", [200] * 2 + [800] * 8)
    with caplog.at_level(logging.WARNING):
        step._warn_short_reads(fq, 400)
    assert "minlength" not in caplog.text


def test_silent_without_minlength_or_reads(step, tmp_dir, caplog):
    with caplog.at_level(logging.WARNING):
        step._warn_short_reads(_fastq(tmp_dir / "a.fq", [50] * 10), 0)
        step._warn_short_reads(_fastq(tmp_dir / "b.fq", []), 400)
    assert caplog.text == ""
