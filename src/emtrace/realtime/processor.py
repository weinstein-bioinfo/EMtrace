"""Batch processor for real-time Nanopore analysis.

Processes incoming FASTQ files through the pipeline in batches:
  1. Preprocess (chopper/porechop)
  2. Host filter
  3. Annotate (rounds 1 & 2)
  4. EM abundance (warm-start from previous F(t) estimates)
  5. Report (optional, on demand or periodic)

Thread-safe: accumulates annotation tables across batches and
updates EM abundance estimates incrementally.

Option B strategy (from ROADMAP):
  Warm-start EM from previous F(t) estimates rather than
  re-running from scratch each time.

Batches and live tables live in ``<job>/watch``. :meth:`BatchProcessor.finalize`
merges them into the job's standard layout (see ``realtime/finalize.py``), so
the other pipeline steps can use the results, and a later watch of the same job
continues from them (:meth:`BatchProcessor.resume`).
"""

from __future__ import annotations

import logging
import os
import shutil
import threading
from dataclasses import dataclass, field
from pathlib import Path
from typing import (Callable, Dict, List, Optional, Sequence, Tuple)

import pandas as pd

from emtrace.analysis.em import EMResult, em_loop, log_hit_structure
from emtrace.core.config import EMtraceConfig
from emtrace.core.pipeline import Pipeline
from emtrace.core.pipelinejob import PipelineJob
from emtrace.models.annotation import (read_annotation_tsv, write_annotation_tsv)
from emtrace.realtime.finalize import (
    MERGED_INPUT,
    FileId,
    completed_batches,
    merge_batches_into_job,
    next_batch_number,
    publish_annotations,
    read_batch_manifest,
    write_batch_manifest,
)
from emtrace.realtime.watcher import FastqWatcher
from emtrace.steps import StepName
from emtrace.steps.report import ReportStep
from emtrace.utils.paths import ensure_dir

logger = logging.getLogger("emtrace.realtime.processor")


# ------------------------------------------------------------------
# Accumulated state
# ------------------------------------------------------------------

# Taxa listed in a snapshot (and named, see BatchProcessor._name_top_taxa)
TOP_TAXA_N = 20


def _taxon_name(values: Sequence[object]) -> str:
    """First real name among lineage values (species, organism, genus); "" if none."""
    for v in values:
        if isinstance(v, str) and v.strip() and v.strip() not in ("NA", "nan", "None"):
            return v.strip()
    return ""


@dataclass
class RealtimeState:
    """Thread-safe accumulated state across batches."""

    # Annotation table (all batches merged)
    annotations: pd.DataFrame = field( default_factory=lambda: pd.DataFrame() )
    # Current EM abundance estimates F(t)
    F: Dict[str, float] = field(default_factory=dict)
    # Latest EM result
    em_result: Optional[EMResult] = None
    # Batch counter
    batch_count: int = 0
    # Total reads processed
    total_reads: int = 0
    # Total files processed
    total_files: int = 0
    # taxid → display name ("" = looked up, none found), for the top-taxa list
    taxa_names: Dict[str, str] = field(default_factory=dict)
    # Lock for thread-safe updates
    _lock: threading.Lock = field( default_factory=threading.Lock )

    def update(
        self,
        new_annot: pd.DataFrame,
        em_result: Optional[EMResult],
        n_reads: int,
        n_files: int,
    ) -> None:
        """Merge new annotations and EM result."""
        with self._lock:
            self.batch_count += 1
            self.total_reads += n_reads
            self.total_files += n_files

            if not new_annot.empty:
                if self.annotations.empty:
                    self.annotations = new_annot.copy()
                else:
                    self.annotations = pd.concat(
                        [self.annotations, new_annot], ignore_index=True)

            if em_result is not None:
                self.em_result = em_result
                self.F = dict(em_result.F)

    def get_snapshot(self) -> dict:
        """Return a summary dict (thread-safe)."""
        with self._lock:
            top_taxa = sorted( self.F.items(), key=lambda x: x[1], reverse=True)[:TOP_TAXA_N]
            return {
                "batch_count": self.batch_count,
                "total_reads": self.total_reads,
                "total_files": self.total_files,
                "n_annotations": len(self.annotations),
                "n_taxa": len(self.F),
                "log_L": (self.em_result.log_L if self.em_result else None),
                "top_taxa": top_taxa,
                # taxid → name for top_taxa; "" when unknown (show the taxid)
                "taxa_names": {t: self.taxa_names.get(t, "") for t, _ in top_taxa},
            }


# ------------------------------------------------------------------
# Callbacks
# ------------------------------------------------------------------

# Type alias for callbacks
BatchCallback = Callable[ [int, RealtimeState], None]

# ------------------------------------------------------------------
# Batch processor
# ------------------------------------------------------------------

class BatchProcessor:
    """Process FASTQ batches through the pipeline.

    Accumulates annotations across batches and runs warm-start
    EM after each batch.

    Parameters
    ----------
    config : EMtraceConfig
        Pipeline configuration.
    job_dir : Path
        The job's output directory. Batches and live tables go to
        ``job_dir/watch`` (``self.output_dir``); finalize() merges them
        into ``job_dir``.
    job_name : str
        PipelineJob name / label.
    ann1 : str
        Annotation database string for round 1.
    ann2 : str
        Annotation database string for round 2.
    anns : str
        Named annotation strategy.
    hostgen : str
        Host genome key(s).
    batch_size : int
        Number of files per batch.
    on_batch_complete : callback, optional
        Called after each batch with (batch_num, state).
    run_report : bool
        Run report step after each batch.
    """

    def __init__(
        self,
        config: EMtraceConfig,
        job_dir: Path,
        job_name: str = "realtime",
        *,
        ann1: str = "",
        ann2: str = "",
        anns: str = "",
        hostgen: str = "",
        batch_size: int = 1,
        on_batch_complete: Optional[BatchCallback] = None,
        run_report: bool = False,
    ) -> None:
        self.config = config
        job_dir = Path(job_dir)
        merged_input = job_dir / "reads" / MERGED_INPUT
        self.job = PipelineJob(
            name=job_name,
            input_reads=merged_input if merged_input.exists() else None,
            output_dir=job_dir,
            results_dir="results",
            num_threads=config.runtime.num_threads,
            profile=config.active_profile,
        )
        self.job.init_dirs()
        self.output_dir = self.job.watch_dir
        self.job_name = job_name
        self.ann1 = ann1
        self.ann2 = ann2
        self.anns = anns
        self.hostgen = hostgen
        self.batch_size = batch_size
        self.on_batch_complete = on_batch_complete
        self.run_report = run_report
        self.state = RealtimeState()
        self.pipeline = Pipeline(config)

        # Accumulated annotation table path
        self._accum_annot_path = self.output_dir / "realtime_annot.tsv"

        # EM output path
        self._em_prob_path = self.output_dir / "results" / "EM.taxid_prob.tsv"
        # Hit-structure warning is logged once per run, not once per batch.
        self._warned_hit_structure = False

        # Directories
        ensure_dir(self.output_dir)
        ensure_dir(self.output_dir / "batches")
        ensure_dir(self.output_dir / "results")

        # Failed batches keep their directory, so numbers are never reused
        self.next_batch = next_batch_number(self.output_dir)

    # ----------------------------------------------------------
    # Public API
    # ----------------------------------------------------------
    def resume(self) -> List[Tuple[Path, Optional[FileId]]]:
        """Reload state from earlier watch sessions of this job.

        Returns the source files of the completed batches, as
        ``(path, identity)`` for ``FastqWatcher(already_processed=...)``.
        Empty when there is nothing to resume.
        """
        batches = completed_batches(self.output_dir)
        if not batches:
            return []

        frames: List[pd.DataFrame] = []
        processed: List[Tuple[Path, Optional[FileId]]] = []
        n_reads = 0
        for batch in batches:
            meta = read_batch_manifest(batch)
            n_reads += int(meta.get("n_reads", 0))
            for entry in meta.get("files", []):
                fid = entry.get("id")
                processed.append((Path(entry["path"]), tuple(fid) if fid else None))
            annot = self._read_batch_annotations(PipelineJob(output_dir=batch))
            if not annot.empty:
                frames.append(annot)

        em_table = ReportStep._read_em_table(self._em_prob_path)
        with self.state._lock:
            self.state.batch_count = len(batches)
            self.state.total_reads = n_reads
            self.state.total_files = len(processed)
            if frames:
                self.state.annotations = pd.concat(frames, ignore_index=True)
            if em_table is not None:
                self.state.F = dict(zip(em_table["taxid"], em_table["prob"]))
        self._name_top_taxa()

        logger.info(
            "Resuming watch: %d earlier batch(es), %d file(s), %d read(s)",
            len(batches), len(processed), n_reads,
        )
        return processed

    def process_files( self, fastq_files: Sequence[Path]) -> RealtimeState:
        """Process one or more FASTQ files as a batch.

        Steps executed per batch:
          1. Concatenate input files
          2. Preprocess
          3. Host filter (if configured)
          4. Annotate
          5. Merge annotations with accumulated table
          6. Warm-start EM on full accumulated annotations
          7. Write updated outputs
          8. Report (if enabled)

        Returns the updated state. Raises if the batch fails: it then has no
        ``batch.json``, so it is left out of the results and its files are
        processed again by a later watch of this job.
        """
        batch_num = self.next_batch
        self.next_batch += 1
        batch_dir = self.output_dir / "batches" / f"batch_{batch_num:04d}"
        ensure_dir(batch_dir)

        logger.info("=== BATCH %d: %d file(s) ===", batch_num, len(fastq_files))

        # Identity before processing, as recorded by the watcher
        sources = [
            (Path(fq).resolve(), FastqWatcher._file_id(Path(fq)))
            for fq in fastq_files
        ]

        # --- Concatenate input files ---
        batch_fq = batch_dir / "batch_input.fq.gz"
        self._concat_fastq(fastq_files, batch_fq)

        # --- Build per-batch PipelineJob ---
        job = self._make_batch_job( batch_fq, batch_dir, batch_num )

        # --- Run pre + flt + ann steps ---
        self.pipeline.run(job, self._batch_steps())

        # --- Count reads processed ---
        n_reads = self._count_reads(batch_fq)

        # --- The batch is complete from here on ---
        write_batch_manifest(batch_dir, batch_num, sources, n_reads)

        # --- Merge batch annotations ---
        batch_annot = self._read_batch_annotations(job)

        # --- EM with warm-start ---
        em_result = self._run_em(batch_annot)

        # --- Update state ---
        self.state.update( batch_annot, em_result, n_reads, len(fastq_files))

        # --- Write accumulated outputs ---
        self._write_outputs()

        # --- Report (optional) ---
        if self.run_report:
            self._run_report_step()

        # --- Callback ---
        if self.on_batch_complete:
            try:
                self.on_batch_complete( batch_num, self.state)
            except Exception as exc:
                logger.warning( "Batch callback error: %s", exc )

        logger.info(
            "Batch %d done: %d reads, %d total, %d taxa, logL=%.2f",
            batch_num,
            n_reads,
            self.state.total_reads,
            len(self.state.F),
            (em_result.log_L if em_result else 0.0),
        )

        return self.state

    def process_single(self, fastq: Path) -> RealtimeState:
        """Process a single FASTQ file."""
        return self.process_files([fastq])

    def finalize(self) -> bool:
        """Merge the batches into the job, then run EM and the report on it.

        EM runs from scratch on the merged tables: the live estimate is
        warm-started, skips ``prune_min_reads`` and mixes annot2/annot1 per
        batch, whereas this one matches ``run -p em``. The final numbers can
        therefore differ slightly from the live ones. The report runs even if
        EM fails, and is QC-only when there are no hits.

        Returns False if no batch completed (nothing to report). Raises if the
        merge or the report fails.
        """
        if not completed_batches(self.output_dir):
            logger.warning("No completed batches — no report generated")
            return False

        merge_batches_into_job(self.job, self.output_dir)
        try:
            self.job = self.pipeline.run(self.job, [StepName.EM])
        except Exception as exc:  # noqa: BLE001
            logger.error("Final EM failed, reporting without it: %s", exc)
            # An EM table from an earlier session would no longer match the reads
            for key, name in (
                ("em_taxid_prob", "EM.taxid_prob.tsv"),
                ("em_taxid_prob_byiter", "EM.taxid_prob.byiter.tsv"),
            ):
                (self.job.reports_dir / name).unlink(missing_ok=True)
                self.job.files.pop(key, None)
        self.job = self.pipeline.run(self.job, [StepName.REP])
        logger.info(
            "Final report -> %s",
            self.job.reports_dir / "reports" / "stats_report.html",
        )
        return True

    # ----------------------------------------------------------
    # Internal: PipelineJob construction
    # ----------------------------------------------------------

    def _make_batch_job(self, batch_fq: Path, batch_dir: Path, batch_num: int
                        ) -> PipelineJob:
        """Build a PipelineJob for one batch."""
        job: PipelineJob = PipelineJob(
            name=f"{self.job_name}_b{batch_num:04d}",
            input_reads=batch_fq,
            output_dir=batch_dir,
            results_dir="results",
            num_threads=self.config.runtime.num_threads,
        )
        job.init_dirs()
        self.pipeline.setup_job(
            job,
            ann1=self.ann1,
            ann2=self.ann2,
            anns=self.anns,
            bgfilter=self.hostgen,
        )
        return job

    # ----------------------------------------------------------
    # Internal: pipeline steps
    # ----------------------------------------------------------

    def _batch_steps(self) -> List[StepName]:
        """Steps to run per batch."""
        steps = [StepName.PORE, StepName.CHOP, StepName.ANN]
        # Insert host filtering if configured
        if self.hostgen and self.hostgen != "0":
            steps.insert(2, StepName.FLT)
        return steps

    def _run_report_step(self) -> None:
        """Run the report on the job's data so far (``--report`` after a batch).

        Reads are merged only at finalize(), so this report's tables are
        current but it has no FASTQ statistics yet.
        """
        try:
            publish_annotations(self.job, self.output_dir)
            if self._em_prob_path.exists():
                shutil.copyfile(
                    self._em_prob_path, self.job.reports_dir / "EM.taxid_prob.tsv",
                )
            self.job = self.pipeline.run(self.job, [StepName.REP])
        except Exception as exc:
            logger.warning(
                "Report step failed: %s", exc,
            )

    # ----------------------------------------------------------
    # Internal: annotation merging
    # ----------------------------------------------------------
    def _read_batch_annotations( self, sample: PipelineJob ) -> pd.DataFrame:
        """Table of a batch for the live estimate: annot2, else annot1."""
        for name in ("annot2.tsv", "annot1.tsv"):
            p = sample.annotation_dir / name
            if p.exists() and p.stat().st_size > 0:
                try:
                    df = read_annotation_tsv(p)
                except Exception as exc:  # noqa: BLE001
                    logger.warning("Could not read %s: %s", p, exc)
                    continue
                if not df.empty:
                    return df
        return pd.DataFrame()

    # ----------------------------------------------------------
    # Internal: EM with warm-start
    # ----------------------------------------------------------
    def _run_em( self, new_annot: pd.DataFrame ) -> Optional[EMResult]:
        """Run EM on accumulated annotations with warm-start.

        Uses the full accumulated annotation table (not just
        the new batch) and warm-starts from previous F(t).
        """
        # Merge new with accumulated
        with self.state._lock:
            if not self.state.annotations.empty:
                full_annot = pd.concat(
                    [self.state.annotations, new_annot],
                    ignore_index=True,
                )
            else:
                full_annot = new_annot

        if full_annot.empty:
            return None

        # Check required columns
        required = {"qseqid", "staxid", "bitscore"}
        if not required.issubset(full_annot.columns):
            logger.warning("Annotation missing columns for EM: %s",
                           required - set(full_annot.columns))
            return None

        # EM parameters from config
        em_cfg = self.config.em
        F_init = self.state.F or None

        try:
            result = em_loop(
                full_annot,
                temp=em_cfg.temp,
                max_iter=em_cfg.max_iter,
                log_L_threshold=em_cfg.log_L_threshold,
                exclude_no_species=(
                    em_cfg.exclude_no_species
                ),
                F_init=F_init,
            )
            if not self._warned_hit_structure:
                self._warned_hit_structure = log_hit_structure(
                    result, logger, self.config.annotation.filter_tophits,
                )
            return result
        except Exception as exc:
            logger.error("EM failed: %s", exc)
            return None

    # ----------------------------------------------------------
    # Internal: I/O helpers
    # ----------------------------------------------------------

    @staticmethod
    def _concat_fastq(
        files: Sequence[Path],
        output: Path,
    ) -> None:
        """Concatenate FASTQ files (gz-aware)."""
        import gzip

        output.parent.mkdir(parents=True, exist_ok=True)

        with gzip.open(output, "wb") as out_fh:
            for fq in files:
                opener = (
                    gzip.open if fq.name.endswith(".gz")
                    else open
                )
                with opener(fq, "rb") as in_fh:
                    shutil.copyfileobj(in_fh, out_fh)

    @staticmethod
    def _count_reads(fastq_path: Path) -> int:
        """Count reads in a FASTQ file (4 lines each)."""
        import gzip

        count = 0
        opener = (
            gzip.open if fastq_path.name.endswith(".gz")
            else open
        )
        try:
            with opener(fastq_path, "rt") as fh:
                for _ in fh:
                    count += 1
        except Exception:
            return 0
        return count // 4

    def _em_n_queries(self) -> int:
        """Reads in the latest EM input (0 if EM has not run)."""
        result = self.state.em_result
        return result.n_queries if result is not None else 0

    @staticmethod
    def _write_em_tsv(
        F: Dict[str, float], path: Path, n_queries: int = 0,
    ) -> None:
        """Write taxid, prob, readn_em TSV (readn_em "NA" if reads unknown)."""
        path.parent.mkdir(parents=True, exist_ok=True)
        sorted_taxa = sorted(
            F.items(),
            key=lambda x: x[1],
            reverse=True,
        )
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("taxid\tprob\treadn_em\n")
            for taxid, prob in sorted_taxa:
                readn_em = prob * n_queries if n_queries > 0 else "NA"
                fh.write(f"{taxid}\t{prob}\t{readn_em}\n")

    def _write_outputs(self) -> None:
        """Write accumulated annotations and EM probs."""
        with self.state._lock:
            if not self.state.annotations.empty:
                write_annotation_tsv(
                    self.state.annotations,
                    self._accum_annot_path,
                )
            if self.state.F:
                self._write_em_tsv(
                    self.state.F, self._em_prob_path, self._em_n_queries(),
                )
        if self.config.output.realtime_abund_table:
            self._write_abundance_tsv()
        self._name_top_taxa()

    def _remember_names(self, abund_df: pd.DataFrame) -> None:
        """Cache taxid → name from a lineage-linked abundance table.

        Taken before the tail trim and other filters, so taxa that drop out
        of abund_table.tsv keep their names in the live top-taxa list.
        """
        cols = [c for c in ("species", "organism", "genus") if c in abund_df.columns]
        if "taxid" not in abund_df.columns or not cols:
            return
        names = {
            str(row[0]): name
            for row in abund_df[["taxid", *cols]].itertuples(index=False)
            if (name := _taxon_name(row[1:]))
        }
        with self.state._lock:
            self.state.taxa_names.update(names)

    def _name_top_taxa(self) -> None:
        """Look up names of top taxa the abundance table did not provide.

        Covers ``output.realtime_abund_table: false`` and failed linkage. NCBI
        taxonomy only (SILVA taxids are named by the table's SILVA linkage);
        best-effort, and taxids without a name are cached as "" so each is
        looked up once.
        """
        with self.state._lock:
            top = sorted(self.state.F.items(), key=lambda x: x[1], reverse=True)[:TOP_TAXA_N]
            missing = [t for t, _ in top if t not in self.state.taxa_names]
            annot = self.state.annotations
        if not missing:
            return
        try:
            from emtrace.analysis.taxonomy_utils import is_silva_annotation
            from emtrace.core.config import TaxonomyType
            from emtrace.models.taxonomy import NCBITaxonomy

            if is_silva_annotation(annot, self.config):
                return
            db = Path(os.path.expandvars(self.config.taxonomy.source(TaxonomyType.NCBI).db))
            if not db.exists():
                return                      # no taxonomy installed: keep taxids, retry later
            tax = NCBITaxonomy(db)
            names = {t: (tax.get_organism(int(t)) if t.isdigit() else "") for t in missing}
        except Exception as exc:  # noqa: BLE001
            logger.debug("Top-taxa name lookup skipped: %s", exc)
            return
        with self.state._lock:
            for t, name in names.items():
                self.state.taxa_names.setdefault(t, name or "")

    def _write_abundance_tsv(self) -> None:
        """Build and write abund_table.tsv from accumulated annotations + EM probs.

        Reuses ReportStep._build_abundance_table() to stay in sync with the
        batch pipeline. Skips gracefully if state is empty or build fails.
        """
        with self.state._lock:
            annot = self.state.annotations.copy()
            has_f = bool(self.state.F)
        if annot.empty or not has_f:
            return

        from emtrace.models.abundance import write_abundance_tsv

        report_step: ReportStep = ReportStep(config=self.config)
        try:
            abund_df = report_step._build_abundance_table(  # pyright: ignore [reportAttributeAccessIssue]
                annot=annot,
                min_qcov=0,
                min_bits=0,
                em_prob_path=self._em_prob_path,
            )
        except Exception as exc:
            logger.warning("Could not build abundance table: %s", exc)
            return

        if abund_df.empty:
            return

        # Link taxonomy lineage (NCBI or SILVA).
        #
        # The taxonomy follows the annotation database (the 'db' column),
        # exactly as ReportStep does: SILVA accessions only resolve in SILVA
        # taxonomy, everything else in NCBI taxonomy.
        try:
            from emtrace.analysis.taxonomy_utils import is_silva_annotation
            from emtrace.core.config import TaxonomyType

            tax_type = (
                TaxonomyType.SILVA if is_silva_annotation(annot, self.config)
                else TaxonomyType.NCBI
            )
            tax_cfg = self.config.taxonomy.source(tax_type)
            tax_db = Path(tax_cfg.db)

            if tax_type == TaxonomyType.SILVA:
                report_step._link_silva_lineage_abund(abund_df, tax_cfg) # pyright: ignore [reportAttributeAccessIssue]
            else:
                report_step._link_ncbi_lineage_abund(abund_df, tax_db) # pyright: ignore [reportAttributeAccessIssue]
                from emtrace.analysis.taxonomy_utils import add_division_field
                add_division_field(abund_df, tax_db, taxid_col="taxid")
        except Exception as exc:
            logger.warning("Taxonomy linkage skipped: %s", exc)
        self._remember_names(abund_df)

        # Same opt-in filter the batch report step applies, so `watch`
        # output matches `run` output.
        if getattr(self.config.report, "exclude_no_species", False):
            from emtrace.analysis.taxonomy_utils import filter_no_species
            abund_df = filter_no_species(abund_df, context="realtime abund_table")

        # Same cumulative-share tail trim `run` applies (report.tail), so the
        # live table means the same thing as the batch one. ReportStep owns
        # the rule, including its refusal to trim without a `division` column
        # -- which is likelier here, since taxonomy linkage above is
        # best-effort and may have been skipped.
        abund_df = report_step._apply_tail_trim(abund_df, self.config.report)  # pyright: ignore [reportAttributeAccessIssue]

        abund_path = self.output_dir / "results" / "abund_table.tsv"
        write_abundance_tsv(abund_df, abund_path)
        logger.info("Abundance table -> %s (%d taxa)", abund_path, len(abund_df))


# ------------------------------------------------------------------
# Convenience: run watcher + processor together
# ------------------------------------------------------------------

def run_realtime(
    watch_dir: str | Path,
    job_dir: str | Path,
    config: EMtraceConfig,
    job_name: str = "realtime",
    *,
    ann1: str = "",
    ann2: str = "",
    anns: str = "",
    hostgen: str = "",
    batch_size: int = 1,
    debounce: float = 2.0,
    scan_existing: bool = False,
    on_batch_complete: Optional[BatchCallback] = None,
    stop_event: Optional[threading.Event] = None,
    run_report: bool = False,
) -> RealtimeState:
    """Run the full real-time pipeline: watch + process.

    Blocks until *stop_event* is set (e.g. via Ctrl+C handler).

    Parameters
    ----------
    watch_dir : Path
        Directory to monitor for FASTQ files.
    job_dir : Path
        The job's output directory; see :class:`BatchProcessor`.
    config : EMtraceConfig
        Pipeline configuration.
    job_name : str
        PipelineJob label.
    ann1, ann2, anns, hostgen
        Pipeline configuration (same as CLI options).
    batch_size : int
        Files per batch (1 = process each file immediately).
    debounce : float
        Seconds to wait after file stops changing.
    scan_existing : bool
        Process FASTQ files already in the directory. Always on when
        resuming an earlier watch of this job, to catch up on files that
        arrived while stopped.
    on_batch_complete : callback
        Called after each batch.
    stop_event : threading.Event
        Set this to stop the loop.
    run_report : bool
        Run report step after each batch.

    Returns
    -------
    RealtimeState
        Final accumulated state.
    """
    if stop_event is None:
        stop_event = threading.Event()

    processor = BatchProcessor(
        config=config,
        job_dir=Path(job_dir),
        job_name=job_name,
        ann1=ann1,
        ann2=ann2,
        anns=anns,
        hostgen=hostgen,
        batch_size=batch_size,
        on_batch_complete=on_batch_complete,
        run_report=run_report,
    )
    processed = processor.resume()

    watcher = FastqWatcher(
        watch_dir,
        debounce_seconds=debounce,
        scan_existing=scan_existing or bool(processed),
        already_processed=processed,
    )

    def _process(files: List[Path]) -> None:
        """Process a batch; a failed batch does not end the watch."""
        try:
            processor.process_files(files)
        except Exception as exc:  # noqa: BLE001
            logger.error("Batch failed: %s", exc)

    watcher.start()
    batch_buffer: List[Path] = []

    try:
        while not stop_event.is_set():
            path = watcher.get(timeout=1.0)
            if path is None:
                continue

            batch_buffer.append(path)

            if len(batch_buffer) >= batch_size:
                _process(batch_buffer)
                batch_buffer.clear()
    except KeyboardInterrupt:
        logger.info("Interrupted — shutting down")
    finally:
        watcher.stop()

        # Process remaining files
        if batch_buffer:
            logger.info(  "Processing %d remaining files", len(batch_buffer) )
            _process(batch_buffer)
            batch_buffer.clear()
        if watcher.queued_count:
            logger.info(
                "%d queued file(s) not processed; a later watch of this job "
                "picks them up", watcher.queued_count,
            )

        try:
            processor.finalize()
        except Exception as exc:  # noqa: BLE001
            logger.error("Final report failed: %s", exc)

    return processor.state
