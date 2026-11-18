"""The clean step keeps the final assembly BAMs and removes intermediates."""

from __future__ import annotations

import pytest

from emtrace.core.config import EMtraceConfig
from emtrace.core.pipelinejob import PipelineJob
from emtrace.steps.assembly import FINAL_BAMS
from emtrace.steps.pack import CleanStep

KEPT = [
    *FINAL_BAMS,
    *(f"{b}.bai" for b in FINAL_BAMS),
    "consensus.fa",
    "vcf.gz",
    "reads.refgen.lowcovmask.bed",
]
REMOVED = [
    "reads.refgen.sam",
    "reads.refgen.bam",
    "reads.draft.md.bam",        # left behind by a failed medaka polish
    "reads.draft.md.bam.bai",
    "vcf.raw.gz",
]


@pytest.fixture
def cfg(tmp_path) -> EMtraceConfig:
    c = EMtraceConfig.create_default()
    c.taxonomy.NCBI.db = str(tmp_path / "no-such-taxonomy.db")
    return c


@pytest.fixture
def job(tmp_path) -> PipelineJob:
    j = PipelineJob(name="j", output_dir=tmp_path / "j")
    j.init_dirs()
    return j


def test_keeps_final_assembly_bams(cfg, job):
    asm = job.assembly_dir
    asm.mkdir(parents=True, exist_ok=True)
    for name in KEPT + REMOVED:
        (asm / name).write_text("x")

    CleanStep(cfg).execute(job)

    assert sorted(p.name for p in asm.iterdir()) == sorted(KEPT)


def test_no_assembly_dir(cfg, job):
    if job.assembly_dir.exists():
        job.assembly_dir.rmdir()
    CleanStep(cfg).execute(job)
    assert not job.assembly_dir.exists()
