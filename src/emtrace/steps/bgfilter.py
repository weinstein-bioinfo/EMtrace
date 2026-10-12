"""Background (host) filtering step — port of ``filter_host_reads()`` from SeqAn.pm.

Depletes background/host reads from a FASTQ file by aligning against one or
more host genome databases.  Two aligners are supported:

``bwa`` (default)
    Uses ``bwa mem`` + pure-Python SAM parsing (replaces Perl's sambamba +
    seqkit).  BWA index is built automatically when missing.
    Host DB lookup: any file matching ``<db>.*`` (index files).

``minimap2``
    Uses ``minimap2 -x map-ont`` + pure-Python PAF parsing.
    Host DB lookup: the exact ``.fna.gz`` path stored in
    ``host.databases.{key}.db``; minimap2 indexes on the fly.

Supports chained filtering over multiple host genomes.

Output per run
--------------
``reads/reads.trim.hflt.fq.gz``
    Reads that passed all host filters.
``reads/hflt.{host_label}.reads.fq.gz``
    Reads filtered by each host genome.

Sample file keys set
--------------------
``"host_filtered_reads"`` / ``"current_reads"``
    Final passing reads.
``"hflt.{label}"``
    Per-host filtered reads.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Dict, List, Set

from emtrace.core.exceptions import StepValidationError
from emtrace.io.fastq import (
    copy_fastq_gz,
    mapped_ids_from_paf,
    partition_fastq,
)
from emtrace.steps import StepName, register_step
from emtrace.steps.base import PipelineStep, resolve_preprocess_input
from emtrace.utils.paths import ensure_dir

if TYPE_CHECKING:
    from emtrace.core.pipelinejob import PipelineJob


@register_step(StepName.FLT)
class BgFilterStep(PipelineStep):
    """Deplete background/host reads using BWA (default) or minimap2.

    Mirrors ``filter_host_reads()`` in ``Lazypipe::SeqAn``.

    For each host genome in ``sample.bgfilter_databases`` the step:

    1. Aligns current reads against the host database.
    2. Collects mapped read IDs (from SAM or PAF, depending on aligner).
    3. Partitions the FASTQ: host-matched reads → per-host file;
       remaining reads → input for the next host filter or final output.

    If no host databases are configured the step is a no-op.
    """

    step_name = "bgfilter"
    step_description = "Background/host read depletion (BWA or minimap2)"

    def validate_inputs(self, sample: PipelineJob, **kwargs) -> bool:
        if sample.input_reads is None or not Path(sample.input_reads).exists():
            raise StepValidationError(
                self.step_name,
                "No input reads available",
            )
        return True

    def run(self, sample: PipelineJob, **kwargs) -> "PipelineJob":
        host_dbs: List[Dict] = sample.bgfilter_databases

        ensure_dir(sample.reads_dir)
        ensure_dir(sample.log_path)
        reads_in = resolve_preprocess_input(
            sample, self.config, self.logger,
            predecessors=("reads.chopped.fq.gz", "reads.porechop.fq.gz"),
        )

        if not host_dbs:
            self.logger.info(
                "No host databases configured — skipping background filtering"
            )
            sample.set_file("host_filtered_reads", str(reads_in))
            return sample

        log_path = sample.log_path / f"{sample.name}.bgfilter.log"
        log_path.write_text("")
        final_pass = sample.reads_dir / "reads.trim.hflt.fq.gz"
        reads_start = reads_in  # remember starting point before loop updates reads_in
        bgf_cfg = self.config.bgfilter
        default_aligner = bgf_cfg.aligner.lower()
        save_mmi = bgf_cfg.minimap.save_mmi_index
        mm2_c: bool = bgf_cfg.minimap.c
        mm2_m: int = bgf_cfg.minimap.m
        mm2_s: int = bgf_cfg.minimap.s
        mm2_k: int = bgf_cfg.minimap.k
        mm2_w: int = bgf_cfg.minimap.w
        min_T: int = bgf_cfg.bwa_mem.T
        min_score: int = bgf_cfg.bwa_mem.min_score
        min_mapq: int = bgf_cfg.bwa_mem.min_mapq

        for host_entry in host_dbs:
            host_db = host_entry.get("db", "")
            if not host_db:
                self.logger.warning(
                    "Host entry missing 'db' field, skipping: %s", host_entry
                )
                continue

            # Per-entry aligner overrides the global default when "search" is set
            aligner = host_entry.get("search", default_aligner).lower()

            label = self._make_label(host_entry)
            self.logger.info(
                "Filtering against host '%s' using %s (%s)",
                label, aligner, host_db,
            )

            reads_flt = sample.reads_dir / f"hflt.{label}.reads.fq.gz"
            reads_pass_tmp = sample.reads_dir / f"reads.hflt.{label}.tmp.fq.gz"

            if aligner == "minimap2":
                mapped_ids = self._align_minimap2(
                    host_db=host_db,
                    reads_in=reads_in,
                    mm2_c=mm2_c,
                    mm2_m=mm2_m,
                    mm2_s=mm2_s,
                    mm2_k=mm2_k,
                    mm2_w=mm2_w,
                    sample=sample,
                    save_mmi=save_mmi,
                )
            else:
                mapped_ids = self._align_bwa(
                    host_db=host_db,
                    reads_in=reads_in,
                    min_T=min_T,
                    min_score=min_score,
                    min_mapq=min_mapq,
                    sample=sample,
                    log_file=log_path,
                )

            self.logger.info("Host '%s': %d reads mapped", label, len(mapped_ids))

            if mapped_ids:
                n_pass, n_flt = partition_fastq(
                    input_fq=reads_in,
                    pass_fq=reads_pass_tmp,
                    flt_fq=reads_flt,
                    filter_ids=mapped_ids,
                    gzip_output=True,
                )
                self.logger.info(
                    "Host '%s': %d passed, %d filtered", label, n_pass, n_flt
                )
                sample.set_file(f"hflt.{label}", reads_flt)
                reads_in = reads_pass_tmp
            else:
                self.logger.info(
                    "Host '%s': no reads mapped, forwarding unchanged", label
                )
                reads_flt.write_bytes(b"")
                sample.set_file(f"hflt.{label}", reads_flt)

        # Rename the last pass-through to the canonical final output
        if reads_in is not reads_start:
            # reads_in was updated by the loop — it points to a temp file
            if reads_in.exists():
                reads_in.rename(final_pass)
            else:
                copy_fastq_gz(reads_start, final_pass)
        else:
            copy_fastq_gz(reads_in, final_pass)

        self.logger.info("Background-filtered reads -> %s", final_pass)
        sample.set_file("host_filtered_reads", final_pass)
        sample.set_file("current_reads", final_pass)
        return sample

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _make_label(self, host_entry: Dict) -> str:
        """Derive a filesystem-safe label from a host database entry."""
        label = (
            host_entry.get("latinName")
            or host_entry.get("commonName")
            or host_entry.get("name")
            or host_entry.get("db", "host")
        )
        label = str(label).split(".")[0].replace(" ", "_")
        return label or "host"

    def _align_bwa(
        self,
        host_db: str,
        reads_in: Path,
        min_T: int,
        min_score: int,
        min_mapq: int,
        sample: "PipelineJob",
        log_file: Path,
    ) -> Set[str]:
        """Run BWA mem and return a set of mapped read IDs.

        The BWA index is built automatically when ``<db>.amb`` / ``<db>.sa``
        are missing.  SAM is parsed in pure Python; no sambamba required.
        """
        from emtrace.wrappers.bwa import BWA

        bwa = BWA(num_threads=sample.num_threads)
        sam_tmp = sample.reads_dir / "hostgen_bwa.sam"
        return bwa.host_filter(
            host_db=host_db,
            query=reads_in,
            score_threshold=min_T,
            min_score=min_score,
            min_mapq=min_mapq,
            sam_tmp=sam_tmp,
            log_file=log_file,
        )

    def _align_minimap2(
        self,
        host_db: str,
        reads_in: Path,
        mm2_c: bool,
        mm2_m: int,
        mm2_s: int,
        mm2_k: int,
        mm2_w: int,
        sample: "PipelineJob",
        save_mmi: bool = False,
    ) -> Set[str]:
        """Run minimap2 (map-ont) and return a set of mapped read IDs.

        All filtering is performed by minimap2 itself via ``-m`` / ``-s``.
        Any read that appears in the PAF output is treated as host-mapped.
        When *save_mmi* is True, a prebuilt ``.mmi`` index is created next to
        the FASTA on first use and reused on subsequent runs.

        Non-default ``-k`` / ``-w`` get their own cached index
        (``<db>.k{k}w{w}.mmi``) so the default and k/w-tuned indexes never
        clobber or get cross-used with one another.
        """
        from emtrace.wrappers.minimap2 import Minimap2

        paf_tmp = sample.reads_dir / "hostgen_mm2.paf"
        # k/w-aware index path: keep tuned indexes separate from the default one.
        if mm2_k > 0 or mm2_w > 0:
            mmi_suffix = f".k{mm2_k}w{mm2_w}.mmi"
        else:
            mmi_suffix = ".mmi"
        mmi_path = Path(str(host_db) + mmi_suffix)  # built only when save_mmi=True
        try:
            mm2 = Minimap2(num_threads=sample.num_threads)
            mm2.host_filter(
                host_db=host_db,
                query=reads_in,
                output=paf_tmp,
                c=mm2_c,
                m=mm2_m,
                s=mm2_s,
                k=mm2_k,
                w=mm2_w,
                mmi_path=mmi_path,
                save_mmi=save_mmi,
            )
            # No post-filter: minimap2 already applied -m/-s thresholds
            return mapped_ids_from_paf(paf_tmp)
        finally:
            if paf_tmp.exists():
                paf_tmp.unlink()
