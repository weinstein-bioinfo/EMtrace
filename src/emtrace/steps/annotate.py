"""Annotation step — port of annotation rounds 1 & 2 from Lazypipe

Runs configurable search tools (minimap2, blastn, diamond blastp) against
annotation databases, converts hits to the 18-field annotation table,
assigns NCBI division/SILVA domain and bacteriophage flags, and partitions reads by
taxonomy (vi/ph/ab/un).

Mirrors the Lazypipe logic:
  - Round 1: iterate ``sample.ann1.databases``, first search targets all
    reads, subsequent searches target unannotated reads (target='un').
  - After all round-1 searches: assign division + bphage, partition reads.
  - Round 2: iterate ``sample.ann2.databases``, each targeting a specific
    read partition (vi/ph/ab) from round 1.
  - After all round-2 searches: assign division + bphage, partition reads.

Output files:
    ``annotation/annot1.tsv``        — round 1 annotation table (merged)
    ``annotation/annot2.tsv``        — round 2 annotation table (merged)
    ``annotation/annot1.{db}.tsv``   — per-database round 1 table, written
                                       immediately after each successful
                                       subsearch so a later subsearch crash in a
                                       chain cannot discard earlier results
    ``annotation/annot2.{db}.tsv``   — per-database round 2 table (same)
    ``annotation/ann1.*.paf``        — raw PAF hits (when save_hits=true)
    ``annotation/ann1.*.tsv``        — raw BLAST/Diamond hits (when save_hits=true)
    ``reads/reads.ann1.vi.fq``  — viral reads (round 1)
    ``reads/reads.ann1.ph.fq``  — phage reads (round 1)
    ``reads/reads.ann1.ab.fq``  — bacterial reads (round 1)
    ``reads/reads.ann1.un.fq``  — unannotated reads (round 1)
    ``reads/reads.ann2.*``      — same for round 2
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, cast, Dict, List, Optional, Set

import pandas as pd

from emtrace.analysis.ref_select import search_reference_db
from emtrace.analysis.taxonomy_utils import (
    add_bphage_field,
    add_division_field,
    add_domain_field,
    add_species_field,
)
from emtrace.core.exceptions import StepValidationError
from emtrace.core.config import MinimapConfig, BlastnConfig, DiamondBlastpConfig, OrfipyConfig, TaxonomyType
from emtrace.io.fasta import (
    fastq_to_fasta,
    reads_info_tsv,
)
from emtrace.io.fastq import (
    iter_fastq,
    read_id_from_header,
)
from emtrace.io.tsv import filter_tophits as filter_tophits_fn, top_hits_per_search
from emtrace.models.annotation import (
    blast_to_annotation,
    empty_annotation_df,
    paf_to_annotation,
    write_annotation_tsv,
)
from emtrace.steps import StepName, register_step
from emtrace.steps.base import PipelineStep, apply_preprocess_subsample
from emtrace.utils.paths import ensure_dir
from emtrace.wrappers.blast import Blastn
from emtrace.wrappers.minimap2 import Minimap2, parse_paf

if TYPE_CHECKING:
    from emtrace.core.pipelinejob import PipelineJob

logger = logging.getLogger(__name__)


# -------------------------------------------------------------------
# Read partitioning helpers (filter_seqs / filter_unseqs ports)
# -------------------------------------------------------------------

def _annotated_ids(annot: pd.DataFrame, col: str = "qseqid") -> Set[str]:
    """Return set of query IDs present in an annotation table."""
    if annot.empty or col not in annot.columns:
        return set()
    return set(annot[col].dropna().astype(str).unique())


def filter_seqs(
    input_fq: Path,
    output_fq: Path,
    annot: pd.DataFrame,
    condition_fn,
    seqid_col: str = "qseqid",
) -> int:
    """Extract reads matching a taxonomy condition from a FASTQ file.

    Filter the annotation DataFrame with
    *condition_fn* and write matching reads from *input_fq*.

    Parameters
    ----------
    input_fq : Path
        Input FASTQ (plain or gzipped).
    output_fq : Path
        Output FASTQ (plain, uncompressed).
    annot : pd.DataFrame
        Annotation table with division/bphage/qseqid columns.
    condition_fn : callable
        A function that takes a DataFrame and returns a boolean mask.
    seqid_col : str
        Column in *annot* containing sequence IDs.

    Returns
    -------
    int
        Number of reads written.
    """
    if annot.empty:
        output_fq.write_text("")
        return 0

    mask = condition_fn(annot)
    series_results = cast(pd.Series, annot.loc[mask, seqid_col])
    keep_ids: Set[str] = set(series_results.dropna().astype(str).unique())

    n = 0
    with open(output_fq, "w") as out:
        for header, seq, plus, qual in iter_fastq(input_fq):
            rid = read_id_from_header(header)
            if rid in keep_ids:
                out.write(f"{header}\n{seq}\n{plus}\n{qual}\n")
                n += 1
    return n


def filter_unseqs(
    input_fq: Path,
    output_fq: Path,
    annot: pd.DataFrame,
    seqid_col: str = "qseqid",
) -> int:
    """Extract reads NOT present in the annotation table.

    Port of ``filter_unseqs()`` from SeqAn.pm.

    Returns
    -------
    int
        Number of unannotated reads written.
    """
    annotated = _annotated_ids(annot, seqid_col)

    if not annotated:
        # No annotations — every read is "unannotated". Write them out as plain
        # (decompressed) FASTQ via iter_fastq. A raw file copy would preserve the
        # input's gzip compression under a ``.fq`` name, which the next
        # subsearch's FASTQ→FASTA step would then misread as text (failing on the
        # gzip magic bytes).
        n = 0
        with open(output_fq, "w") as out:
            for header, seq, plus, qual in iter_fastq(input_fq):
                out.write(f"{header}\n{seq}\n{plus}\n{qual}\n")
                n += 1
        return n

    n = 0
    with open(output_fq, "w") as out:
        for header, seq, plus, qual in iter_fastq(input_fq):
            rid = read_id_from_header(header)
            if rid not in annotated:
                out.write(f"{header}\n{seq}\n{plus}\n{qual}\n")
                n += 1
    return n


# -------------------------------------------------------------------
# Taxonomy condition functions (mirrors Perl $cond_vi, $cond_ph, etc.)
# -------------------------------------------------------------------

# -- NCBI conditions (use 'division' field) --

def _cond_vi(df: pd.DataFrame) -> pd.Series:
    """NCBI — viral reads: division=='Viruses' AND bphage=='no'."""
    div_ok = df["division"].str.lower() == "viruses"
    bph_ok = df["bphage"].str.lower() == "no"
    return div_ok & bph_ok


def _cond_ph(df: pd.DataFrame) -> pd.Series:
    """NCBI — phage reads: division=='Phages' OR bphage=='yes'."""
    div_ph = df["division"].str.lower() == "phages"
    bph_yes = df["bphage"].str.lower() == "yes"
    return div_ph | bph_yes


def _cond_ab(df: pd.DataFrame) -> pd.Series:
    """NCBI — bacterial reads: division=='Bacteria'."""
    return df["division"].str.lower() == "bacteria"


# -- SILVA conditions (use 'domain' field; no viruses, no phages) --

def _cond_silva_ba(df: pd.DataFrame) -> pd.Series:
    """SILVA — bacterial/archaeal reads: domain in {'Archaea', 'Bacteria'}."""
    dom = df["domain"].str.lower()
    return dom.isin({"archaea", "bacteria"})


def _cond_silva_euk(df: pd.DataFrame) -> pd.Series:
    """SILVA — eukaryotic reads: domain=='Eukaryota'."""
    return df["domain"].str.lower() == "eukaryota"


# -------------------------------------------------------------------
# Individual search runners
# -------------------------------------------------------------------

def _run_minimap_search(
    query: Path,
    db_entry: Dict[str, str],
    config: MinimapConfig,
    sample,
    acc2taxid_path: Optional[Path] = None,
    output: Optional[Path] = None,
    filter_tophits: bool = False,
) -> pd.DataFrame:
    """Run minimap2 and return annotation DataFrame.

    Port of ``annotate_minimap()`` from Lazypipe
    """
    db_path = db_entry.get("db", '')
    db_name = db_entry.get("name", 'NA')
    minimap2 = Minimap2(num_threads=sample.num_threads)

    # Resolve the database path to pass to minimap2.
    # Priority:
    #   1. db_path ends with .mmi and exists → use as-is
    #   2. db_path + ".mmi" exists → use pre-built index
    #   3. Otherwise → pass db_path directly (FASTA; minimap2 indexes on the fly)
    db_p = Path(db_path)
    mmi_candidate = Path(db_path + ".mmi")
    if str(db_path).endswith(".mmi") and db_p.exists():
        effective_db = db_p
    elif mmi_candidate.exists():
        effective_db = mmi_candidate
    else:
        effective_db = db_p   # FASTA — minimap2 handles on-the-fly indexing

    # Run alignment.
    # Thresholding is done inside minimap2 via -m (chain score) and -s (DP score, requires -c).
    # No post-hoc score filtering: all PAF output already passed minimap2's own filters.
    # c=True: full DP alignment; AS:i: score; accurate qcov/scov; 2-5x slower.
    # c=False (default): chain-only; s1:i: score; qcov/scov set to NA; fast.
    c = config.c
    log_file = sample.log_path / f"minimap2.{db_name}.log"
    sample.log_path.mkdir(parents=True, exist_ok=True)
    hits = minimap2.align(
        db=effective_db,
        query=query,
        preset=config.preset,
        s=config.s,
        m=config.m,
        c=c,
        f=config.f,
        secondary=config.secondary,
        n=config.max_target_seqs,
        k=config.k,
        w=config.w,
        K=config.K,
        output=output,
        log_file=log_file,
    )

    _stderr = minimap2.last_stderr or ""
    logger.info(
        "[minimap2/%s] stderr captured: %d bytes, log_file=%s",
        db_name, len(_stderr), log_file,
    )
    if _stderr:
        for _line in _stderr.splitlines():
            _line = _line.strip()
            if _line and ("[WARNING]" in _line or "warning" in _line.lower()):
                logger.warning("[minimap2/%s] %s", db_name, _line)

    if isinstance(hits, Path):
        hits = parse_paf(hits, c=c)

    if hits.empty:
        return empty_annotation_df()

    # Keep only top hits when annotation.filter_tophits is on; otherwise every
    # hit reaches EM (downstream top-hit consumers select their own).
    if filter_tophits:
        hits = filter_tophits_fn(hits, query_col="qseqid", score_col="bitscore", retain_ties=True)

    # Convert to annotation table format; pass c so coverage fields are set correctly
    annot = paf_to_annotation(hits, search="Minimap2", db=db_name, c=c)

    # Resolve staxid from acc2taxid map if available.
    # Check alongside the effective db path (strip .mmi if present).
    if acc2taxid_path is None:
        base = str(effective_db)
        if base.endswith(".mmi"):
            base = base[:-4]
        candidate = Path(base + ".acc2taxid")
        if candidate.exists():
            acc2taxid_path = candidate

    if acc2taxid_path and acc2taxid_path.exists():
        _apply_acc2taxid(annot, acc2taxid_path)

    return annot


def _run_blastn_search(
    query_fa: Path,
    db_entry: Dict[str, str],
    config: BlastnConfig,
    sample,
    output: Optional[Path] = None,
    filter_tophits: bool = False,
) -> pd.DataFrame:
    """Run blastn and return annotation DataFrame.

    Port of ``annotate_blastn()`` from SeqAn.pm.
    """
    db_path = db_entry.get("db", '')
    db_name = db_entry.get("name", 'NA')

    blastn: Blastn = Blastn(num_threads=sample.num_threads)
    hits = blastn.search(
        query=query_fa,
        db=db_path,
        evalue=config.evalue,
        max_target_seqs=config.max_target_seqs,
        output=output,
    )

    if isinstance(hits, Path):
        from emtrace.wrappers.blast import parse_blast_file
        hits = parse_blast_file(hits)

    if hits.empty:
        return empty_annotation_df()

    # Filter by min bitscore
    hits["bitscore"] = pd.to_numeric(
        hits["bitscore"], errors="coerce"
    )
    hits = hits[hits["bitscore"] >= config.bitscore].copy()

    if hits.empty:
        return empty_annotation_df()

    if filter_tophits:
        hits = filter_tophits_fn(hits, query_col="qseqid", score_col="bitscore", retain_ties=True)

    return blast_to_annotation(hits, search="Blastn", db=db_name, dbtype="nucl")


def _run_blastp_search(
    query_fa: Path,
    db_entry: Dict[str, str],
    config: BlastnConfig,
    sample,
    output: Optional[Path] = None,
    filter_tophits: bool = False,
) -> pd.DataFrame:
    """Run blastp (protein DB) and return annotation DataFrame."""
    from emtrace.wrappers.blast import Blastp

    db_path = db_entry.get("db", "")
    db_name = db_entry.get("name", "NA")

    blastp: Blastp = Blastp(num_threads=sample.num_threads)
    hits = blastp.search(
        query=query_fa,
        db=db_path,
        evalue=config.evalue,
        max_target_seqs=config.max_target_seqs,
        output=output,
    )

    if isinstance(hits, Path):
        from emtrace.wrappers.blast import parse_blast_file
        hits = parse_blast_file(hits)

    if hits.empty:
        return empty_annotation_df()

    hits["bitscore"] = pd.to_numeric(hits["bitscore"], errors="coerce")
    hits = hits[hits["bitscore"] >= config.bitscore].copy()

    if hits.empty:
        return empty_annotation_df()

    if filter_tophits:
        hits = filter_tophits_fn(hits, query_col="qseqid", score_col="bitscore", retain_ties=True)
    return blast_to_annotation(hits, search="Blastp", db=db_name, dbtype="prot")


def _run_diamondx_search(
    query_fa: Path,
    db_entry: Dict[str, str],
    config: DiamondBlastpConfig,
    sample,
    output: Optional[Path] = None,
    filter_tophits: bool = False,
) -> pd.DataFrame:
    """Run Diamond blastx and return annotation DataFrame.

    Translates nucleotide reads on-the-fly against a protein DB.
    No separate ORF detection step — Diamond handles translation internally.
    """
    from emtrace.wrappers.diamond import DiamondBlastx, parse_diamond_file

    db_path = db_entry.get("db", "")
    db_name = db_entry.get("name", "NA")

    diamond = DiamondBlastx(num_threads=sample.num_threads)
    hits = diamond.search(
        query=query_fa,
        db=db_path,
        evalue=config.evalue,
        max_target_seqs=config.max_target_seqs,
        min_score=config.bitscore,
        sensitivity=config.mode,
        output=output,
    )

    if isinstance(hits, Path):
        hits = parse_diamond_file(hits)

    if hits.empty:
        return empty_annotation_df()

    hits["bitscore"] = pd.to_numeric(hits["bitscore"], errors="coerce")
    hits = hits[hits["bitscore"] >= config.bitscore].copy()
    if hits.empty:
        return empty_annotation_df()

    if filter_tophits:
        hits = filter_tophits_fn(
            hits, query_col="qseqid", score_col="bitscore", retain_ties=True,
        )

    return blast_to_annotation(hits, search="Diamond blastx", db=db_name, dbtype="prot")


def _run_diamondp_search(
    query_fa: Path,
    db_entry: Dict[str, str],
    config: DiamondBlastpConfig,
    orfipy_config: OrfipyConfig,
    sample,
    output: Optional[Path] = None,
    filter_tophits: bool = False,
) -> pd.DataFrame:
    """Detect ORFs with orfipy then search with Diamond blastp.

    Port of ``annotate_diamondp()`` from Lazypipe:
      1. Run orfipy on the nucleotide query FASTA to produce a
         peptide FASTA and a BED-derived orf_map.
      2. Run Diamond blastp (protein vs protein DB) on the peptides.
      3. Map hits back to original read IDs using the orf_map and
         store ORF coordinates in the ``orf`` annotation column.
    """
    from emtrace.wrappers.diamond import DiamondBlastp, parse_diamond_file
    from emtrace.wrappers.orfipy import Orfipy

    db_path = db_entry.get("db", "")
    db_name = db_entry.get("name", "NA")

    # ------------------------------------------------------------------
    # 1. ORF detection
    # ------------------------------------------------------------------
    orf_dir = query_fa.parent / f"{query_fa.stem}.orfs"
    orfipy = Orfipy(config=orfipy_config)
    pep_fa, orf_map = orfipy.detect_orfs(query_fa, orf_dir)

    if pep_fa is None or orf_map.empty:
        return empty_annotation_df()

    # ------------------------------------------------------------------
    # 2. Diamond blastp on peptide sequences
    # ------------------------------------------------------------------
    diamond = DiamondBlastp(num_threads=sample.num_threads)
    hits = diamond.search(
        query=pep_fa,
        db=db_path,
        evalue=config.evalue,
        max_target_seqs=config.max_target_seqs,
        min_score=config.bitscore,
        sensitivity=config.mode,
        output=output,
    )

    if isinstance(hits, Path):
        hits = parse_diamond_file(hits)

    if hits.empty:
        return empty_annotation_df()

    # ------------------------------------------------------------------
    # 3. Map ORF IDs back to original read IDs
    #
    # orf_map columns: orf_id, read_id, orf_start, orf_stop, strand
    # hits["qseqid"] == orf_id produced by orfipy
    # ------------------------------------------------------------------
    hits = hits.merge(
        orf_map[["orf_id", "read_id", "orf_start", "orf_stop", "strand"]],
        left_on="qseqid",
        right_on="orf_id",
        how="left",
    )

    # Restore original read ID in qseqid; build orf coordinate string
    hits["qseqid"] = hits["read_id"].fillna(hits["qseqid"])
    hits["orf_tag"] = (
        hits["orf_start"].astype(str)
        + "-"
        + hits["orf_stop"].astype(str)
        + ":"
        + hits["strand"].fillna("+")
    )

    # ------------------------------------------------------------------
    # 4. Filter by bitscore and keep top hits per read
    # ------------------------------------------------------------------
    hits["bitscore"] = pd.to_numeric(hits["bitscore"], errors="coerce")
    hits = hits[hits["bitscore"] >= config.bitscore].copy()
    if hits.empty:
        return empty_annotation_df()

    if filter_tophits:
        hits = filter_tophits_fn(
            hits, query_col="qseqid", score_col="bitscore", retain_ties=True,
        )
    hits = hits.reset_index(drop=True)

    # ------------------------------------------------------------------
    # 5. Convert to annotation table format
    # ------------------------------------------------------------------
    annot = blast_to_annotation(hits, search="Diamond blastp", db=db_name, dbtype="prot")

    # blast_to_annotation looks for "ORF=" in qseqid (won't find it since
    # we restored the original read ID).  Overwrite with the actual orf coords.
    annot["orf"] = hits["orf_tag"].values

    return annot


# -------------------------------------------------------------------
# Acc2taxid resolution
# -------------------------------------------------------------------

def _apply_acc2taxid(
    annot: pd.DataFrame,
    acc2taxid_path: Path,
) -> None:
    """Map sseqid to staxid using a TSV acc2taxid file (in-place).

    The acc2taxid file is a two-column headerless TSV: accession, taxid.
    """
    acc_map: Dict[str, str] = {}
    with open(acc2taxid_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                acc_map[parts[0]] = parts[1]

    if acc_map and "sseqid" in annot.columns:
        annot["staxid"] = annot["sseqid"].map(
            lambda s: acc_map.get(str(s), "NA")
        )


# -------------------------------------------------------------------
# Main annotation step
# -------------------------------------------------------------------

@register_step(StepName.ANN)
class AnnotateStep(PipelineStep):
    """Annotation rounds 1 & 2.

    Port of the annotation logic from emtrace.pl lines 192-384
    and SeqAn.pm annotate_* functions.

    Configuration is driven by ``sample.ann1`` and ``sample.ann2``
    (``AnnotationRound`` dataclasses), each containing a list of
    database entries with keys: db, name, search, target.
    """

    step_name = "annotate"
    step_description = "Taxonomic annotation (rounds 1 & 2)"

    def validate_inputs(self, sample: PipelineJob, **kwargs) -> bool:
        current = sample.current_reads
        if current is None:
            raise StepValidationError(
                self.step_name,
                "No current reads; run preprocess/bgfilter first",
            )
        if not Path(current).exists():
            raise StepValidationError(
                self.step_name,
                f"Input reads not found: {current}",
            )
        if (
            not sample.ann1.databases
            and not sample.ann2.databases
        ):
            self.logger.warning("No annotation databases configured")
        return True

    def run(self, sample: PipelineJob, **kwargs) -> PipelineJob:
        ensure_dir(sample.reads_dir)
        ensure_dir(sample.annotation_dir)
        ensure_dir(sample.output_dir)
        ensure_dir(sample.log_path)

        # Apply subsampling if not already done in preprocessing steps.
        # Idempotent: reuses input.subsample.fq when preprocessing already ran.
        apply_preprocess_subsample(sample, self.config, self.logger)

        reads: Path = self._select_query_reads(sample)  # type: ignore[assignment]

        # Create reads info TSV (seqid + length)
        info_path = Path(str(reads) + ".info.tsv")
        if not info_path.exists():
            reads_info_tsv(reads, info_path)
            sample.set_file("reads_info", info_path)

        # -- Round 1 --
        annot1_path = sample.annotation_dir / "annot1.tsv"
        if sample.ann1.databases:
            self.logger.info("=== ANNOTATION ROUND 1 ===")
            # _run_round assigns taxonomy and persists annot1.tsv incrementally
            # (complete with division/bphage), so a chain that fails midway still
            # leaves a usable annot1.tsv on disk.
            annot1 = self._run_round(
                databases=sample.ann1.databases,
                reads=reads,
                annot_path=annot1_path,
                sample=sample,
                round_num=1,
            )

            tax_type1 = self._detect_tax_type(sample.ann1.databases)
            write_annotation_tsv(annot1, annot1_path)
            sample.set_file("annot1", annot1_path)

            # Partition reads by taxonomy
            self._partition_reads(
                reads=reads,
                annot=annot1,
                prefix="ann1",
                sample=sample,
                tax_type=tax_type1,
            )
        else:
            annot1 = empty_annotation_df()

        # -- Round 2 --
        annot2_path = sample.annotation_dir / "annot2.tsv"
        if sample.ann2.databases:
            self.logger.info("=== ANNOTATION ROUND 2 ===")
            annot2 = self._run_round(
                databases=sample.ann2.databases,
                reads=reads,
                annot_path=annot2_path,
                sample=sample,
                round_num=2,
            )

            tax_type2 = self._detect_tax_type(sample.ann2.databases)
            write_annotation_tsv(annot2, annot2_path)
            sample.set_file("annot2", annot2_path)

            # Partition reads (round 2)
            self._partition_reads(
                reads=reads,
                annot=annot2,
                prefix="ann2",
                sample=sample,
                tax_type=tax_type2,
            )

        # Offline mode: re-search viral reads against the assembly DB so the
        # GUI's assembly-candidate picker has a curated viral-only hit list.
        self._run_assembly_candidates(sample)

        return sample

    # ---------------------------------------------------------------
    # Assembly candidate search (offline mode)
    # ---------------------------------------------------------------

    def _run_assembly_candidates(self, sample: PipelineJob) -> None:
        """Re-search viral reads against ``assembly.database``.

        Writes ``annotation/assembly_candidates.tsv`` so the GUI's
        "Pick from viral hits" combo can offer curated viral accessions; the
        assembly step's ref_auto reuses it as its primary hits.
        No-ops when ``assembly.database`` is not configured.
        """
        db_name = self.config.assembly.database
        if not db_name:
            return

        db_entry = self.config.databases.get(db_name)
        if db_entry is None:
            self.logger.warning(
                "assembly.database '%s' not in config.databases; "
                "skipping assembly_candidates search",
                db_name,
            )
            return

        # Pick the most-recent viral reads (round 2 preferred, fall back to round 1)
        viral_reads = sample.get_file("reads_ann2_vi") or sample.get_file("reads_ann1_vi")
        if not viral_reads:
            self.logger.info(
                "Offline mode: no viral reads partition registered; skipping assembly_candidates"
            )
            return
        viral_path = Path(viral_reads)
        if not viral_path.exists() or viral_path.stat().st_size == 0:
            self.logger.info(
                "Offline mode: viral reads file empty (%s); skipping assembly_candidates",
                viral_path,
            )
            return

        self.logger.info(
            "Searching viral reads against %s for assembly candidates", db_name,
        )
        # Same search as assembly's ref_auto (one path for GUI and CLI):
        # base-level alignment, secondary hits kept, staxid from acc2taxid.
        ra_mm = self.config.assembly.ref_auto.minimap
        sample.log_path.mkdir(parents=True, exist_ok=True)
        try:
            annot = search_reference_db(
                viral_path, db_entry, sample.num_threads, k=ra_mm.k, w=ra_mm.w,
                log_file=sample.log_path / f"minimap2.assembly_candidates.{db_name}.log",
            )
        except Exception as exc:  # noqa: BLE001
            self.logger.warning("Assembly candidate search failed: %s", exc)
            return

        out_path = sample.annotation_dir / "assembly_candidates.tsv"
        write_annotation_tsv(annot, out_path)
        sample.set_file("assembly_candidates", out_path)
        self.logger.info(
            "Wrote %d assembly candidate hit(s) to %s", len(annot), out_path,
        )

    # ---------------------------------------------------------------
    # Round execution
    # ---------------------------------------------------------------

    def _run_round(
        self,
        databases: List[Dict[str, str]],
        reads: Path,
        annot_path: Path,
        sample: PipelineJob,
        round_num: int,
    ) -> pd.DataFrame:
        """Execute one annotation round over a list of databases.

        For round 1: first search targets all reads; subsequent
        searches target unannotated reads (target='un').

        For round 2: each search targets a specific partition
        from round 1 (vi, ph, ab) as specified by target key.

        Returns the merged annotation DataFrame, with taxonomy fields
        (division/bphage for NCBI, domain for SILVA) already assigned.
        """
        accumulated = empty_annotation_df()
        self.logger.info(
            "Round %d hits: %s", round_num,
            "top hits only (annotation.filter_tophits)"
            if self.config.annotation.filter_tophits
            else "all hits (EM weighs alternatives; tables use top hits)",
        )

        # Resolve taxonomy source once for the whole round; used to classify
        # each subsearch's hits as they complete so the incrementally-persisted
        # tables are complete (see below).
        tax_type = self._detect_tax_type(databases)
        tax_db = self._taxonomy_db_path(tax_type)
        warned_no_taxdb = False

        for i, db_entry in enumerate(databases):
            search = db_entry.get("search", "")
            target = db_entry.get("target", "all")
            db_name = db_entry.get("name", db_entry.get("db", ""))

            self.logger.info(
                "Round %d, search %d: %s on %s (target=%s)",
                round_num, i + 1, search, db_name, target,
            )

            # Determine input reads for this search
            query_fq = self._resolve_target_reads(
                target=target,
                reads=reads,
                sample=sample,
                round_num=round_num,
            )

            if query_fq is None or not query_fq.exists():
                self.logger.warning(
                    "Target reads not found for target='%s', skipping",
                    target,
                )
                continue

            # Run the search
            annot_df = self._dispatch_search(
                search=search,
                query_fq=query_fq,
                db_entry=db_entry,
                sample=sample,
                round_num=round_num,
            )

            if annot_df is not None and not annot_df.empty:
                # Assign taxonomy now (per-row by staxid, so per-subsearch is
                # equivalent to classifying the merged table) so the tables we
                # persist below are complete — division/bphage (NCBI) or domain
                # (SILVA).  This lets a partial round left behind by a later
                # subsearch failure feed em/report directly (e.g. via
                # `--pipe em,rep,pack`) without rerunning annotation.
                if tax_db:
                    self._assign_taxonomy_fields(annot_df, tax_type, tax_db)
                elif not warned_no_taxdb:
                    self.logger.warning(
                        "Taxonomy DB not found for %s — division field will be "
                        "empty. Check config.taxonomy.%s.db",
                        tax_type.value, tax_type.value,
                    )
                    warned_no_taxdb = True

                # Persist this subsearch's table immediately (persist-then-fail):
                # a later subsearch crash in the chain must not discard results
                # already obtained from an earlier successful subsearch.
                per_db_path = (
                    sample.annotation_dir
                    / f"annot{round_num}.{db_name}.tsv"
                )
                write_annotation_tsv(annot_df, per_db_path)

                accumulated = pd.concat(
                    [accumulated, annot_df],
                    ignore_index=True,
                )
                # Incrementally refresh the merged round table so it always
                # reflects completed (and fully classified) subsearches, even if
                # a later one fails.
                write_annotation_tsv(accumulated, annot_path)

                self.logger.info(
                    "  -> %d hits (total: %d); persisted %s",
                    len(annot_df),
                    len(accumulated),
                    per_db_path.name,
                )
            else:
                self.logger.info("  -> 0 hits")

            # In round 1, after each search, extract unannotated
            # reads for the next search (target='un')
            if round_num == 1:
                un_path = sample.reads_dir / "reads.ann1.un.fq"
                filter_unseqs(reads, un_path, accumulated)
                sample.set_file("reads_ann1_un", un_path)

        return accumulated

    def _assign_taxonomy_fields(
        self,
        annot: pd.DataFrame,
        tax_type: TaxonomyType,
        tax_db: Path,
    ) -> None:
        """Assign taxonomy fields to *annot* in place.

        NCBI: ``division`` + ``bphage`` + ``species``.  SILVA: ``domain``
        (plus ``species`` and the full lineage, added by ``add_domain_field``).
        No-op for an empty table or a missing taxonomy DB.

        ``species`` is assigned here, rather than only in the reporting step,
        because the ``em`` step runs before ``rep`` and needs it to honour
        ``em.exclude_no_species``.
        """
        if annot.empty or not tax_db:
            return
        if tax_type == TaxonomyType.SILVA:
            add_domain_field(annot, tax_db, taxid_col="staxid")
        else:
            add_division_field(annot, tax_db, taxid_col="staxid")
            add_species_field(annot, tax_db, taxid_col="staxid")
            add_bphage_field(
                annot, tax_db,
                phfilter_path=self._phage_filter_path(),
                taxid_col="staxid",
            )

    def _resolve_target_reads(
        self,
        target: str,
        reads: Path,
        sample: PipelineJob,
        round_num: int,
    ) -> Optional[Path]:
        """Resolve which FASTQ to use as query for a given target.

        Round 1:
          - 'all': use all reads
          - 'un': use unannotated reads from previous search

        Round 2:
          - 'vi'/'ph'/'ab': use partitioned reads from round 1
        """
        if target == "all":
            return reads

        if round_num == 1 and target == "un":
            un = sample.get_file("reads_ann1_un")
            return un if un and un.exists() else reads

        if round_num == 2:
            key = f"reads_ann1_{target}"
            f = sample.get_file(key)
            if f and f.exists():
                return f
            self.logger.warning(
                "Round 2 target '%s' file not found (key=%s)",
                target, key,
            )
            return None

        return reads

    def _dispatch_search(
        self,
        search: str,
        query_fq: Path,
        db_entry: Dict[str, str],
        sample: PipelineJob,
        round_num: int = 1,
    ) -> Optional[pd.DataFrame]:
        """Dispatch to the appropriate search function."""
        search_lower = search.lower()
        db_name = db_entry.get("name", db_entry.get("db", "search"))
        ann = self.config.annotation

        def _out(ext: str, cfg_save: bool) -> Optional[Path]:
            if not cfg_save:
                return None
            return sample.annotation_dir / f"ann{round_num}.{db_name}.{ext}"

        if search_lower in ("minimap", "minimap2"):
            return _run_minimap_search(
                query=query_fq,
                db_entry=db_entry,
                config=ann.minimap,
                sample=sample,
                output=_out("paf", ann.minimap.save_hits),
                filter_tophits=ann.filter_tophits,
            )

        # For blastn/blastp/diamond, convert FASTQ -> FASTA first
        query_fa = Path(str(query_fq).removesuffix(".gz") + ".fa")
        fastq_to_fasta(query_fq, query_fa)
        self.logger.debug("converted fastq to fasta: %s", query_fa)
        try:
            if search_lower in ("blastn",):
                return _run_blastn_search(
                    query_fa=query_fa,
                    db_entry=db_entry,
                    config=ann.blastn,
                    sample=sample,
                    output=_out("blast.tsv", ann.blastn.save_hits),
                    filter_tophits=ann.filter_tophits,
                )
            elif search_lower in ("diamond", "diamondx", "diamond.blastx"):
                return _run_diamondx_search(
                    query_fa=query_fa,
                    db_entry=db_entry,
                    config=ann.diamond_blastp,
                    sample=sample,
                    output=_out("diamond.tsv", ann.diamond_blastp.save_hits),
                    filter_tophits=ann.filter_tophits,
                )
            elif search_lower in ("diamondp", "diamond.blastp"):
                return _run_diamondp_search(
                    query_fa=query_fa,
                    db_entry=db_entry,
                    config=ann.diamond_blastp,
                    orfipy_config=ann.orfipy,
                    sample=sample,
                    output=_out("diamond.tsv", ann.diamond_blastp.save_hits),
                    filter_tophits=ann.filter_tophits,
                )
            else:
                self.logger.warning("Unknown search type: '%s', skipping",search)
                return None
        finally:
            if query_fa.exists():
                query_fa.unlink()

    # ---------------------------------------------------------------
    # Read partitioning
    # ---------------------------------------------------------------

    def _partition_reads(
        self,
        reads: Path,
        annot: pd.DataFrame,
        prefix: str,
        sample: PipelineJob,
        tax_type: TaxonomyType = TaxonomyType.NCBI,
    ) -> None:
        """Partition reads by taxonomy after an annotation round.

        NCBI path (division field): vi / ph / ab / un
        SILVA path (domain field):  ab (archaea+bacteria) / euk / un
        No bphage assignment or reads.ph produced for SILVA.
        """
        un_path = sample.reads_dir / f"reads.{prefix}.un.fq"

        # Partition on each read's top hits only: with all hits kept, a weak
        # sub-best hit (e.g. viral) must not pull a read into another
        # partition. A no-op when annotation.filter_tophits is on. The
        # unannotated set is unaffected (any hit means annotated).
        annot = top_hits_per_search(annot).copy()

        if tax_type == TaxonomyType.SILVA:
            # Ensure domain column exists
            if "domain" not in annot.columns:
                annot["domain"] = ""

            ab_path = sample.reads_dir / f"reads.{prefix}.ab.fq"
            euk_path = sample.reads_dir / f"reads.{prefix}.euk.fq"

            n_ab = filter_seqs(reads, ab_path, annot, _cond_silva_ba)
            n_euk = filter_seqs(reads, euk_path, annot, _cond_silva_euk)
            n_un = filter_unseqs(reads, un_path, annot)

            self.logger.info(
                "Partition '%s' (SILVA): ab=%d euk=%d un=%d",
                prefix, n_ab, n_euk, n_un,
            )

            sample.set_file(f"reads_{prefix}_ab", ab_path)
            sample.set_file(f"reads_{prefix}_euk", euk_path)
            sample.set_file(f"reads_{prefix}_un", un_path)

        else:
            # NCBI: ensure division and bphage columns exist
            for col, default in [("division", ""), ("bphage", "NA")]:
                if col not in annot.columns:
                    annot[col] = default

            vi_path = sample.reads_dir / f"reads.{prefix}.vi.fq"
            ph_path = sample.reads_dir / f"reads.{prefix}.ph.fq"
            ab_path = sample.reads_dir / f"reads.{prefix}.ab.fq"

            n_vi = filter_seqs(reads, vi_path, annot, _cond_vi)
            n_ph = filter_seqs(reads, ph_path, annot, _cond_ph)
            n_ab = filter_seqs(reads, ab_path, annot, _cond_ab)
            n_un = filter_unseqs(reads, un_path, annot)

            self.logger.info(
                "Partition '%s' (NCBI): vi=%d ph=%d ab=%d un=%d",
                prefix, n_vi, n_ph, n_ab, n_un,
            )

            sample.set_file(f"reads_{prefix}_vi", vi_path)
            sample.set_file(f"reads_{prefix}_ph", ph_path)
            sample.set_file(f"reads_{prefix}_ab", ab_path)
            sample.set_file(f"reads_{prefix}_un", un_path)

    # ---------------------------------------------------------------
    # Config helpers
    # ---------------------------------------------------------------

    def _select_query_reads(self, sample: PipelineJob) -> Optional[Path]:
        """Return the best available query FASTQ for annotation.

        Priority:
          1. reads/reads.trim.hflt.fq.gz  — bgfiltered reads
          2. reads/reads.chopped.fq.gz     — chopper output
          3. reads/reads.porechop.fq.gz    — porechop output
          4. reads/input.subsample.fq      — pre-existing subsample
          5. sample.input_reads (subsampled on-the-fly if config active)
        """
        candidates = [
            sample.reads_dir / "reads.trim.hflt.fq.gz",
            sample.reads_dir / "reads.chopped.fq.gz",
            sample.reads_dir / "reads.porechop.fq.gz",
            sample.reads_dir / "input.subsample.fq",
        ]
        for path in candidates:
            if path.exists():
                self.logger.info("Annotation query reads: %s", path)
                return path

        # Fall back to raw input_reads, subsampling if configured
        if sample.input_reads:
            sub = self.config.preprocess.subsample
            if sub.active and sub.size > 0:
                sub_path = sample.reads_dir / "input.subsample.fq"
                from emtrace.io.fastq import subsample_fastq
                from emtrace.utils.paths import ensure_dir
                ensure_dir(sample.reads_dir)
                n = subsample_fastq(Path(sample.input_reads), sub_path, size=sub.size)
                self.logger.info(
                    "Subsampled %d reads (size=%s) for annotation: %s",
                    n, sub.size, sub_path,
                )
                sample.set_file("current_reads", sub_path)
                return sub_path
            self.logger.info("Annotation query reads (raw input): %s", sample.input_reads)
            return Path(sample.input_reads)

        self.logger.warning("No input reads found for annotation")
        return None

    def _detect_tax_type(self, databases: List[Dict[str, str]]) -> TaxonomyType:
        """Infer the taxonomy a round reports through from its databases.

        SILVA if any database resolves to SILVA taxonomy
        (``EMtraceConfig.taxonomy_for_database``), otherwise NCBI.
        """
        for db_entry in databases:
            name = db_entry.get("name", "") or db_entry.get("db", "")
            if self.config.taxonomy_for_database(name) == TaxonomyType.SILVA:
                return TaxonomyType.SILVA
        return TaxonomyType.NCBI

    def _taxonomy_db_path(self, tax_type: TaxonomyType = TaxonomyType.NCBI) -> Optional[Path]:
        """Return the taxonomy SQLite file for *tax_type* (NCBI by default)."""
        tax_cfg = self.config.taxonomy.source(tax_type)
        return Path(tax_cfg.db) if tax_cfg.db else None

    def _phage_filter_path(self) -> Optional[Path]:
        """Locate the ICTV phage filter TSV."""
        # Check common locations
        candidates = [
            Path(__file__).parent.parent
            / "data" / "ICTV.bphage.filter.tsv",
            Path(__file__).parent.parent.parent.parent
            / "emtrace" / "perl" / "ICTV.bphage.filter.tsv",
        ]
        for c in candidates:
            if c.exists():
                return c
        return None
