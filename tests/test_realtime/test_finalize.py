"""Watch results reach the job's standard layout, and a watch can resume.

The per-batch pipeline needs external tools and databases, so these tests stub
it with fake batch outputs. Everything after a batch runs for real: the merge,
EM, the report, resume and clean.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import List, Optional

import pandas as pd
import pytest

from emtrace.analysis.ref_select import candidates_are_current
from emtrace.core.config import EMtraceConfig
from emtrace.core.pipelinejob import PipelineJob
from emtrace.realtime.finalize import (
    BATCH_INPUT,
    MERGED_INPUT,
    completed_batches,
    dirs_overlap,
    has_unmerged_batches,
    is_watch_job,
    merge_batches_into_job,
    merged_input,
    next_batch_number,
    watch_conflict,
)
from emtrace.realtime.processor import BatchProcessor, run_realtime
from emtrace.realtime.watcher import FastqWatcher
from emtrace.steps import StepName
from emtrace.steps.pack import CleanStep
from emtrace.steps.report import ReportStep


def fastq(*ids: str) -> str:
    return "".join(f"@{i}\nACGTACGT\n+\nIIIIIIII\n" for i in ids)


def annot_table(*rows) -> pd.DataFrame:
    """Annotation rows of (read id, taxid)."""
    return pd.DataFrame({
        "search": ["minimap2"] * len(rows),
        "db": ["db"] * len(rows),
        "qseqid": [r[0] for r in rows],
        "sseqid": [f"acc{r[1]}" for r in rows],
        "orf": ["NA"] * len(rows),
        "staxid": [str(r[1]) for r in rows],
        "bitscore": [100.0] * len(rows),
        "qcov": [1.0] * len(rows),
        "species": [f"Species {r[1]}" for r in rows],
    })


def fake_batch_outputs(job: PipelineJob, reads: List[str], taxid: int = 1,
                       annot2: bool = False, candidates: bool = False) -> None:
    """What PORE, CHOP and ANN leave in a batch directory."""
    for fname in ("reads.porechop.fq.gz", "reads.chopped.fq.gz"):
        with gzip.open(job.reads_dir / fname, "wt") as fh:
            fh.write(fastq(*reads))
    (job.reads_dir / "reads.ann1.vi.fq").write_text(fastq(*reads))
    (job.reads_dir / "reads.adaptertrim.tmp.fq.gz").write_text("temporary")
    table = annot_table(*[(r, taxid) for r in reads])
    table.to_csv(job.annotation_dir / "annot1.tsv", sep="\t", index=False)
    if annot2:
        table.to_csv(job.annotation_dir / "annot2.tsv", sep="\t", index=False)
    else:
        annot_table().to_csv(job.annotation_dir / "annot2.tsv", sep="\t", index=False)
    if candidates:
        table.to_csv(job.annotation_dir / "assembly_candidates.tsv", sep="\t", index=False)


@pytest.fixture
def cfg(tmp_path) -> EMtraceConfig:
    c = EMtraceConfig.create_default()
    c.taxonomy.NCBI.db = str(tmp_path / "no-such-taxonomy.db")
    return c


def make_processor(cfg, job_dir: Path, monkeypatch, *, fail: Optional[set] = None,
                   annot2: bool = False, candidates: bool = False,
                   no_hits: bool = False) -> BatchProcessor:
    """A processor whose batch steps are stubbed; EM and REP stay real."""
    proc = BatchProcessor(config=cfg, job_dir=job_dir, job_name="job")
    monkeypatch.setattr(proc, "_batch_steps", lambda: ["stub"])
    real_run = proc.pipeline.run

    def run(job, steps, **kwargs):
        if steps != ["stub"]:
            return real_run(job, steps, **kwargs)
        if job.name in (fail or set()):
            raise RuntimeError(f"{job.name} failed")
        with gzip.open(job.input_reads, "rt") as fh:
            reads = [line[1:].strip() for line in fh.readlines()[::4]]
        fake_batch_outputs(job, reads, annot2=annot2, candidates=candidates)
        if no_hits:
            for name in ("annot1.tsv", "annot2.tsv"):
                annot_table().to_csv(job.annotation_dir / name, sep="\t", index=False)
        return job

    monkeypatch.setattr(proc.pipeline, "run", run)
    return proc


def source_fastq(directory: Path, name: str, *reads: str) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / name
    path.write_text(fastq(*reads))
    return path


def read_ids(path: Path) -> List[str]:
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt") as fh:
        return [line[1:].strip() for line in fh.readlines()[::4]]


@pytest.fixture
def stopped_job(cfg, tmp_path, monkeypatch):
    """A job with two batches, not merged yet."""
    proc = make_processor(cfg, tmp_path / "job", monkeypatch, candidates=True)
    src = tmp_path / "minknow"
    proc.process_files([source_fastq(src, "a.fq", "r1", "r2")])
    proc.process_files([source_fastq(src, "b.fq", "r3")])
    return proc


class TestMerge:
    def test_reads_are_concatenated_and_batch_copies_removed(self, stopped_job):
        proc = stopped_job
        merge_batches_into_job(proc.job, proc.output_dir)
        reads = proc.job.reads_dir

        assert read_ids(reads / MERGED_INPUT) == ["r1", "r2", "r3"]
        assert read_ids(reads / "reads.chopped.fq.gz") == ["r1", "r2", "r3"]
        assert read_ids(reads / "reads.ann1.vi.fq") == ["r1", "r2", "r3"]
        for batch in completed_batches(proc.output_dir):
            assert not (batch / BATCH_INPUT).exists()
            assert not list((batch / "reads").glob("reads.chopped*"))

    def test_temporary_files_are_not_merged(self, stopped_job):
        merge_batches_into_job(stopped_job.job, stopped_job.output_dir)
        assert not list(stopped_job.job.reads_dir.glob("*tmp*"))

    def test_input_reads_recorded_for_open_job(self, stopped_job):
        job = stopped_job.job
        merge_batches_into_job(job, stopped_job.output_dir)
        assert job.input_reads == job.reads_dir / MERGED_INPUT
        assert PipelineJob.read_job_manifest(job.output_dir)["input_reads"] == str(job.input_reads)
        assert merged_input(job.output_dir) == job.input_reads
        assert job.current_reads == job.reads_dir / "reads.chopped.fq.gz"

    def test_annotation_rounds_stay_separate(self, cfg, tmp_path, monkeypatch):
        proc = make_processor(cfg, tmp_path / "job", monkeypatch, annot2=True)
        proc.process_files([source_fastq(tmp_path / "src", "a.fq", "r1", "r2")])
        merge_batches_into_job(proc.job, proc.output_dir)
        for name in ("annot1.tsv", "annot2.tsv"):
            table = pd.read_csv(proc.job.annotation_dir / name, sep="\t")
            assert list(table["qseqid"]) == ["r1", "r2"], name

    def test_empty_round_leaves_no_table(self, stopped_job):
        # A header-only annot2.tsv would make report and EM ignore annot1.tsv.
        merge_batches_into_job(stopped_job.job, stopped_job.output_dir)
        annotation = stopped_job.job.annotation_dir
        assert (annotation / "annot1.tsv").exists()
        assert not (annotation / "annot2.tsv").exists()
        assert ReportStep(stopped_job.config)._resolve_annotation_tsv(stopped_job.job).name == "annot1.tsv"

    def test_failed_batch_is_left_out(self, cfg, tmp_path, monkeypatch):
        proc = make_processor(cfg, tmp_path / "job", monkeypatch, fail={"job_b0002"})
        src = tmp_path / "src"
        proc.process_files([source_fastq(src, "a.fq", "r1")])
        with pytest.raises(RuntimeError):
            proc.process_files([source_fastq(src, "b.fq", "bad")])
        proc.process_files([source_fastq(src, "c.fq", "r3")])

        merge_batches_into_job(proc.job, proc.output_dir)
        assert read_ids(proc.job.reads_dir / MERGED_INPUT) == ["r1", "r3"]
        assert proc.state.total_files == 2

    def test_repeating_a_merge_does_not_duplicate(self, stopped_job):
        job = stopped_job.job
        merge_batches_into_job(job, stopped_job.output_dir)
        merge_batches_into_job(job, stopped_job.output_dir)
        assert read_ids(job.reads_dir / "reads.chopped.fq.gz") == ["r1", "r2", "r3"]
        assert list(pd.read_csv(job.annotation_dir / "annot1.tsv", sep="\t")["qseqid"]) == ["r1", "r2", "r3"]

    def test_candidates_stay_current_for_assembly(self, stopped_job):
        job = stopped_job.job
        merge_batches_into_job(job, stopped_job.output_dir)
        assert candidates_are_current(
            job.annotation_dir / "assembly_candidates.tsv",
            job.reads_dir / "reads.chopped.fq.gz",
        )


class TestFinalize:
    def test_report_has_tables_and_fastq_stats(self, stopped_job):
        proc = stopped_job
        assert proc.finalize() is True
        reports = proc.job.reports_dir / "reports"

        abund = pd.read_csv(reports / "abund_table.tsv", sep="\t")
        assert abund["readn"].sum() == 3
        assert (reports / "annot_table.tsv").exists()
        assert (proc.job.reports_dir / "EM.taxid_prob.tsv").exists()

        stages = pd.read_csv(reports / "stats_fastq.csv")
        assert {"input", "chopped"} <= set(stages["stage"])
        assert "chopped" in (reports / "stats_report.html").read_text()

    def test_assembly_finds_the_merged_reads(self, stopped_job):
        from emtrace.steps.assembly import AssemblyStep

        proc = stopped_job
        proc.finalize()
        step = AssemblyStep(proc.config)
        assert step._find_assembly_reads(proc.job) == proc.job.reads_dir / "reads.chopped.fq.gz"
        proc.config.assembly.input_reads = "viral"
        assert step._find_assembly_reads(proc.job) == proc.job.reads_dir / "reads.ann1.vi.fq"

    def test_no_hits_still_gives_a_qc_report(self, cfg, tmp_path, monkeypatch):
        proc = make_processor(cfg, tmp_path / "job", monkeypatch, no_hits=True)
        proc.process_files([source_fastq(tmp_path / "src", "a.fq", "r1")])
        assert proc.finalize() is True
        reports = proc.job.reports_dir / "reports"
        assert (reports / "stats_report.html").exists()
        assert (reports / "stats_fastq.csv").exists()
        assert not (reports / "abund_table.tsv").exists()

    def test_no_completed_batch_means_no_report(self, cfg, tmp_path, monkeypatch):
        proc = make_processor(cfg, tmp_path / "job", monkeypatch)
        assert proc.finalize() is False
        assert not (proc.job.reports_dir / "reports" / "stats_report.html").exists()

    def test_em_failure_still_reports(self, stopped_job, monkeypatch):
        proc = stopped_job
        stale = proc.job.reports_dir / "EM.taxid_prob.tsv"
        stale.write_text("taxid\tprob\treadn_em\n999\t1.0\t1\n")
        stubbed_run = proc.pipeline.run

        def run(job, steps, **kwargs):
            if steps == [StepName.EM]:
                raise RuntimeError("EM failed")
            return stubbed_run(job, steps, **kwargs)

        monkeypatch.setattr(proc.pipeline, "run", run)
        assert proc.finalize() is True
        assert not stale.exists()
        assert (proc.job.reports_dir / "reports" / "stats_report.html").exists()


class TestResume:
    def test_state_is_restored(self, cfg, stopped_job, monkeypatch):
        job_dir = stopped_job.job.output_dir
        stopped_job.finalize()

        proc = make_processor(cfg, job_dir, monkeypatch)
        processed = proc.resume()

        assert proc.state.batch_count == 2
        assert proc.state.total_reads == 3
        assert proc.state.total_files == 2
        assert len(proc.state.annotations) == 3
        assert proc.state.F  # warm start from the live estimate
        assert sorted(p.name for p, _ in processed) == ["a.fq", "b.fq"]
        assert proc.job.input_reads == job_dir / "reads" / MERGED_INPUT

    def test_a_new_session_adds_to_the_results(self, cfg, tmp_path, stopped_job, monkeypatch):
        job_dir = stopped_job.job.output_dir
        stopped_job.finalize()

        proc = make_processor(cfg, job_dir, monkeypatch)
        proc.resume()
        assert proc.next_batch == 3
        proc.process_files([source_fastq(tmp_path / "minknow", "c.fq", "r4")])
        assert (proc.output_dir / "batches" / "batch_0001" / "batch.json").exists()
        assert proc.finalize() is True

        assert read_ids(job_dir / "reads" / MERGED_INPUT) == ["r1", "r2", "r3", "r4"]
        assert list(pd.read_csv(job_dir / "annotation" / "annot1.tsv", sep="\t")["qseqid"]) == ["r1", "r2", "r3", "r4"]
        abund = pd.read_csv(job_dir / "results" / "reports" / "abund_table.tsv", sep="\t")
        assert abund["readn"].sum() == 4

    def test_nothing_to_resume(self, cfg, tmp_path, monkeypatch):
        proc = make_processor(cfg, tmp_path / "job", monkeypatch)
        assert proc.resume() == []
        assert proc.next_batch == 1

    def test_failed_batch_keeps_its_directory_and_is_redone(self, cfg, tmp_path, monkeypatch):
        proc = make_processor(cfg, tmp_path / "job", monkeypatch, fail={"job_b0001"})
        bad = source_fastq(tmp_path / "src", "bad.fq", "r1")
        with pytest.raises(RuntimeError):
            proc.process_files([bad])
        assert proc.next_batch == 2

        again = make_processor(cfg, tmp_path / "job", monkeypatch)
        assert again.resume() == []                  # not counted as processed
        assert next_batch_number(again.output_dir) == 2   # directory not reused

    def test_processed_files_are_skipped_even_after_a_rename(self, cfg, tmp_path, stopped_job, monkeypatch):
        job_dir = stopped_job.job.output_dir
        src = tmp_path / "minknow"
        stopped_job.finalize()
        (src / "b.fq").rename(src / "renamed.fq")
        source_fastq(src, "new.fq", "r9")

        processed = make_processor(cfg, job_dir, monkeypatch).resume()
        watcher = FastqWatcher(
            src, debounce_seconds=0.1, poll_interval=0.05,
            scan_existing=True, already_processed=processed,
        )
        watcher.start()
        try:
            first = watcher.get(timeout=5)
            second = watcher.get(timeout=1)
        finally:
            watcher.stop()
        assert first is not None and first.name == "new.fq"
        assert second is None

    def test_run_realtime_resumes_and_finalizes(self, cfg, tmp_path, stopped_job, monkeypatch):
        import threading

        job_dir = stopped_job.job.output_dir
        stopped_job.finalize()
        src = tmp_path / "minknow"
        source_fastq(src, "c.fq", "r4")
        stop = threading.Event()

        def batch_steps(self):
            return ["stub"]

        real_run = stopped_job.pipeline.__class__.run

        def run(self, job, steps, **kwargs):
            if steps != ["stub"]:
                return real_run(self, job, steps, **kwargs)
            fake_batch_outputs(job, read_ids(job.input_reads))
            stop.set()
            return job

        monkeypatch.setattr(BatchProcessor, "_batch_steps", batch_steps)
        monkeypatch.setattr(stopped_job.pipeline.__class__, "run", run)

        state = run_realtime(src, job_dir, cfg, job_name="job", debounce=0.1, stop_event=stop)

        assert state.total_files == 3                 # a.fq and b.fq not redone
        assert read_ids(job_dir / "reads" / MERGED_INPUT) == ["r1", "r2", "r3", "r4"]
        assert (job_dir / "results" / "reports" / "stats_report.html").exists()


class TestWatchJobRules:
    def test_a_new_directory_can_host_a_watch(self, tmp_path):
        assert watch_conflict(tmp_path / "new") is None

    def test_job_with_its_own_reads_is_refused(self, tmp_path):
        PipelineJob(
            name="j", input_reads=tmp_path / "reads.fq", output_dir=tmp_path / "j",
        ).init_dirs()
        assert "input reads" in watch_conflict(tmp_path / "j")
        assert not is_watch_job(tmp_path / "j")

    def test_job_with_results_is_refused(self, tmp_path):
        (tmp_path / "j" / "annotation").mkdir(parents=True)
        (tmp_path / "j" / "annotation" / "annot1.tsv").write_text("qseqid\n")
        assert "annotation" in watch_conflict(tmp_path / "j")

    def test_watched_job_can_continue(self, stopped_job):
        job_dir = stopped_job.job.output_dir
        stopped_job.finalize()
        assert watch_conflict(job_dir) is None
        assert is_watch_job(job_dir)

    def test_job_created_for_watching_counts_before_its_first_start(self, tmp_path):
        PipelineJob(name="j", output_dir=tmp_path / "j").init_dirs()
        assert is_watch_job(tmp_path / "j")
        assert watch_conflict(tmp_path / "j") is None

    def test_overlap_in_either_direction(self, tmp_path):
        assert dirs_overlap(tmp_path / "a", tmp_path / "a" / "job")
        assert dirs_overlap(tmp_path / "a" / "minknow", tmp_path / "a")
        assert not dirs_overlap(tmp_path / "a", tmp_path / "b")

    def test_unmerged_batches_are_detected(self, stopped_job):
        assert has_unmerged_batches(stopped_job.output_dir)
        merge_batches_into_job(stopped_job.job, stopped_job.output_dir)
        assert not has_unmerged_batches(stopped_job.output_dir)


class TestClean:
    def test_keeps_what_resume_needs(self, cfg, stopped_job, monkeypatch):
        job_dir = stopped_job.job.output_dir
        stopped_job.finalize()
        batch = completed_batches(stopped_job.output_dir)[0]
        (batch / "annotation" / "ann1.db.paf").write_text("raw search output")
        (batch / "logs" / "chopper.log").write_text("log")

        CleanStep(cfg).execute(stopped_job.job)

        assert (batch / "batch.json").exists()
        assert (batch / "annotation" / "annot1.tsv").exists()
        assert (batch / "annotation" / "assembly_candidates.tsv").exists()
        assert not (batch / "annotation" / "ann1.db.paf").exists()
        assert not (batch / "logs").exists()
        assert not (batch / "job.yaml").exists()

        proc = make_processor(cfg, job_dir, monkeypatch)
        proc.resume()
        assert proc.state.total_reads == 3
        assert len(proc.state.annotations) == 3

    def test_unmerged_batches_are_not_touched(self, cfg, stopped_job):
        batch = completed_batches(stopped_job.output_dir)[0]
        CleanStep(cfg).execute(stopped_job.job)
        assert (batch / BATCH_INPUT).exists()
        assert read_ids(batch / "reads" / "reads.chopped.fq.gz") == ["r1", "r2"]

    def test_failed_batches_are_removed(self, cfg, tmp_path, monkeypatch):
        proc = make_processor(cfg, tmp_path / "job", monkeypatch, fail={"job_b0001"})
        with pytest.raises(RuntimeError):
            proc.process_files([source_fastq(tmp_path / "src", "a.fq", "r1")])
        CleanStep(cfg).execute(proc.job)
        assert not (proc.output_dir / "batches" / "batch_0001").exists()

    def test_other_jobs_are_unaffected(self, cfg, tmp_path):
        job = PipelineJob(name="j", output_dir=tmp_path / "j")
        job.init_dirs()
        CleanStep(cfg).execute(job)   # no watch directory: nothing to do


class TestCli:
    """`emtrace watch` works on the job directory RES/SAMPLE."""

    @staticmethod
    def invoke(tmp_path, *args):
        from click.testing import CliRunner

        from emtrace.__main__ import cli

        return CliRunner().invoke(cli, [
            "watch", "--dir", str(tmp_path / "minknow"),
            "-r", str(tmp_path / "res"), "-s", "run01", *args,
        ])

    def test_passes_the_job_directory(self, tmp_path, monkeypatch):
        from emtrace.realtime.processor import RealtimeState

        (tmp_path / "minknow").mkdir()
        seen = {}

        def fake_run(**kwargs):
            seen.update(kwargs)
            return RealtimeState()

        monkeypatch.setattr("emtrace.realtime.processor.run_realtime", fake_run)
        result = self.invoke(tmp_path)
        assert result.exit_code == 0, result.output
        assert Path(seen["job_dir"]) == tmp_path / "res" / "run01"

    def test_refuses_a_job_with_its_own_reads(self, tmp_path):
        (tmp_path / "minknow").mkdir()
        PipelineJob(
            name="run01", input_reads=tmp_path / "reads.fq",
            output_dir=tmp_path / "res" / "run01",
        ).init_dirs()
        assert self.invoke(tmp_path).exit_code == 1

    def test_refuses_a_job_inside_the_watch_directory(self, tmp_path):
        from click.testing import CliRunner

        from emtrace.__main__ import cli

        (tmp_path / "minknow").mkdir()
        result = CliRunner().invoke(cli, [
            "watch", "--dir", str(tmp_path / "minknow"),
            "-r", str(tmp_path / "minknow"), "-s", "run01",
        ])
        assert result.exit_code == 1


def test_running_the_report_twice_does_not_duplicate_sorted_reads(cfg, tmp_path):
    reads = tmp_path / "reads.fq"
    reads.write_text(fastq("r1"))
    table = annot_table(("r1", 1)).assign(division="Viruses", family="F", genus="G")
    out = tmp_path / "sorted"
    step = ReportStep(cfg)
    step._sort_reads_by_taxonomy(reads, table, out)
    step._sort_reads_by_taxonomy(reads, table, out)
    assert read_ids(next(out.rglob("*.fq"))) == ["r1"]
