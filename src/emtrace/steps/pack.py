"""Pack and Clean steps.

PackStep: creates a ``{sample_name}.tar.gz`` archive of the shareable
output files (reports, tables, figures, assemblies, viral reads) in the
parent directory of ``sample.reports_dir``.

CleanStep: removes intermediate BAM, SAM and temp files that are no
longer needed after packing, and the intermediates of watch-mode batches.

Configuration (from ``config.output``):
    pack_reads    — include all reads/ directory (default 0 / False)
    pack_reads_vi — include viral FASTQ reads (default 1 / True)
"""

from __future__ import annotations

import shutil
import tarfile
from pathlib import Path
from typing import TYPE_CHECKING, List

from emtrace.realtime.finalize import (
    BATCH_INPUT,
    BATCH_MANIFEST,
    batch_dirs,
    completed_batches,
    resume_tables,
)
from emtrace.steps import StepName, register_step
from emtrace.steps.base import PipelineStep

if TYPE_CHECKING:
    from emtrace.core.pipelinejob import PipelineJob


# ------------------------------------------------------------------
# Pack step
# ------------------------------------------------------------------

@register_step(StepName.PACK)
class PackStep(PipelineStep):
    """Create a shareable tar.gz archive of pipeline outputs."""

    step_name = "pack"
    step_description = "Pack output files for sharing"

    def validate_inputs(self, job: PipelineJob, **kwargs) -> bool:
        if not job.reports_dir.exists():
            self.logger.warning(
                "reports_dir not found — nothing to pack"
            )
            return False
        return True

    def run(self, job: PipelineJob, **kwargs) -> PipelineJob:
        res = job.reports_dir             # output_dir/results
        parent = res.parent               # output_dir
        basename = res.name               # "results"
        archive = parent / f"{job.name or basename}.tar.gz"

        out_cfg = self.config.output
        pack_reads        = bool(getattr(out_cfg, "pack_reads",        False))
        pack_reads_vi     = bool(getattr(out_cfg, "pack_reads_vi",     True))
        pack_assembly_bam = bool(getattr(out_cfg, "pack_assembly_bam", False))

        reports_subdir = res / "reports"
        plots_subdir   = res / "plots"
        asm_dir        = job.assembly_dir

        # Collect candidate paths
        candidates: List[Path] = [
            # Result tables (written by ReportStep into reports/)
            reports_subdir / "abund_table.tsv",
            reports_subdir / "annot_table.tsv",
            reports_subdir / "abund_table.xlsx",
            reports_subdir / "annot_table.xlsx",
            # Statistics (ReportStep → reports/)
            reports_subdir / "stats_fastq.csv",
            reports_subdir / "stats_div.csv",
            reports_subdir / "stats_assembly.csv",
            reports_subdir / "stats_report.html",
            # EM probability table (ReportStep → results/ root)
            res / "EM.taxid_prob.tsv",
            # Plots directory (sankey, sunburst, heatmaps — PlotStep)
            plots_subdir,
            # Assembly: consensus, VCF, reference, low-coverage mask
            asm_dir / "consensus.fa",
            asm_dir / "vcf.gz",
            asm_dir / "combined_ref.fa",
            asm_dir / "reads.refgen.lowcovmask.bed",
            asm_dir / "ref_selection.tsv",           # ref_auto choice (when used)
            asm_dir / "consensus.bcftools.fa",       # unpolished (medaka runs only)
            asm_dir / "reads.polished.lowcovmask.bed",
        ]

        # Assembly BAM (optional — large file)
        if pack_assembly_bam:
            candidates.append(asm_dir / "reads.refgen.md.bam")
            candidates.append(asm_dir / "reads.refgen.md.bam.bai")

        # Reads FASTQ files
        reads_dir = job.reads_dir
        if pack_reads:
            candidates.append(reads_dir)
        elif pack_reads_vi:
            # Prefer ann2 viral reads; fall back to ann1 if ann2 was not run
            vi_ann2 = reads_dir / "reads.ann2.vi.fq"
            vi_ann2gz = reads_dir / "reads.ann2.vi.fq.gz"
            vi_ann1 = reads_dir / "reads.ann1.vi.fq"
            vi_ann1gz = reads_dir / "reads.ann1.vi.fq.gz"
            if vi_ann2.exists() or vi_ann2gz.exists():
                candidates.extend([vi_ann2, vi_ann2gz])
            else:
                candidates.extend([vi_ann1, vi_ann1gz])

        existing = [p for p in candidates if p.exists()]
        if not existing:
            self.logger.warning("No output files found to pack")
            return job

        self.logger.info( "Packing %d paths -> %s", len(existing), archive )

        with tarfile.open(archive, "w:gz") as tf:
            for path in existing:
                # Archive relative to output_dir so all paths (results/, assembly/)
                # are preserved with their natural layout inside the tar.
                arcname = Path(job.name or basename) / path.relative_to(parent)
                tf.add(path, arcname=str(arcname), recursive=True)

        self.logger.info("Archive created: %s", archive)
        job.set_file("pack_archive", archive)
        return job


# ------------------------------------------------------------------
# Clean step
# ------------------------------------------------------------------

@register_step(StepName.CLEAN)
class CleanStep(PipelineStep):
    """Remove intermediate files no longer needed after packing."""

    step_name = "clean"
    step_description = "Remove intermediate pipeline files"

    def validate_inputs(self, job: PipelineJob, **kwargs) -> bool:
        return True

    def run(self, job: PipelineJob, **kwargs) -> PipelineJob:
        removed = 0

        def _try_unlink(p: Path) -> None:
            nonlocal removed
            try:
                p.unlink()
                removed += 1
                self.logger.debug("Removed: %s", p)
            except OSError as exc:
                self.logger.warning("Could not remove %s: %s", p, exc)

        def _try_rmtree(d: Path) -> None:
            nonlocal removed
            if d.exists():
                try:
                    shutil.rmtree(d)
                    removed += 1
                    self.logger.debug("Removed dir: %s", d)
                except OSError as exc:
                    self.logger.warning("Could not remove %s: %s", d, exc)

        # Intermediate alignment files in assembly_dir
        asm_dir = job.assembly_dir
        if asm_dir.exists():
            for pattern in ("*.bam", "*.bai", "*.sam", "*.raw.gz"):
                for p in asm_dir.glob(pattern):
                    _try_unlink(p)

        # Temporary reads files
        reads_dir = job.reads_dir
        if reads_dir.exists():
            for pattern in ("*tmp*",):
                for p in reads_dir.glob(pattern):
                    _try_unlink(p)

        # Stale assembler working directories (legacy)
        for dname in ("Roadmaps", "assembler_out"):
            _try_rmtree(job.output_dir / dname)

        # Watch batches: only what a later watch of the job reloads is kept.
        # Batches not merged yet still hold the job's only copy of their reads;
        # batches that failed are redone by a later watch.
        watch_dir = job.watch_dir
        completed = set(completed_batches(watch_dir))
        for batch in batch_dirs(watch_dir):
            if batch not in completed:
                _try_rmtree(batch)
            elif not (batch / BATCH_INPUT).exists():
                keep = {batch / BATCH_MANIFEST, *resume_tables(batch)}
                for p in batch.rglob("*"):
                    if p.is_file() and p not in keep:
                        _try_unlink(p)
                # Deepest first, so a parent is empty once its children are gone
                for d in sorted(
                    (p for p in batch.rglob("*") if p.is_dir()),
                    key=lambda p: len(p.parts), reverse=True,
                ):
                    if not any(d.iterdir()):
                        d.rmdir()

        self.logger.info("Clean: removed %d items", removed)
        return job
