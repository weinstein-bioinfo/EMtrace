"""Reference-based assembly step — port of ``assemble_refbased()``.

Pipeline:
  1. Resolve reference genomes from config.assembly.refseqs (FASTA paths)
     and config.assembly.refaccs (NCBI accessions fetched via Entrez).
     When both are empty and assembly.ref_auto.enabled, accessions are
     selected from the reads' hits (analysis/ref_select.py) and the choice
     is written to ``ref_selection.tsv``.
     All references are concatenated into a single combined FASTA.
  2. Align reads to combined reference with minimap2. Input: all
     host-filtered reads by default (assembly.input_reads); minimap2 with
     -k 11 -w 5 and minimap2's default -m/-s (assembly.minimap).
  3. samtools sort → fixmate → sort → markdup → index
  4. samtools depth -J → low-coverage mask → maskfasta (deleted positions
     count as covered, so called deletions are applied rather than masked)
  5. bcftools mpileup -X ont-sup -B + call + norm -m -any → VCF → filter (report tier) → tabix
  6. bcftools consensus -i (dominant-allele tier) → final consensus FASTA
  7. Optional (assembly.medaka.enabled, ONT only): medaka polishes a draft
     (VCF applied to the unmasked reference) using the same reads mapped to
     it; the polished sequence is re-masked below min_coverage.

Output (under ``{sample.output_dir}/assembly/``):
    ``reads.refgen.md.bam``          — deduplicated, indexed BAM
    ``consensus.fa``                 — consensus sequence
    ``vcf.gz``                       — variant calls at AF >= report_min_af
    ``reads.refgen.lowcovmask.bed``  — low-coverage regions
    ``ref_selection.tsv``            — automatic reference selection (ref_auto)
  with medaka polishing, ``consensus.fa`` is the polished consensus and also:
    ``consensus.bcftools.fa``        — the unpolished bcftools consensus
    ``consensus.medaka.unmasked.fa`` — polished sequence before masking
    ``reads.polished.md.bam``        — reads mapped to the polished sequence
    ``reads.polished.lowcovmask.bed`` — its low-coverage mask (polished coords)
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, List, Optional

import pandas as pd

from emtrace.analysis.ref_select import (
    candidates_are_current, search_reference_db, select_references,
)
from emtrace.analysis.taxonomy_utils import add_species_field
from emtrace.core.config import AssemblyConfig, DatabaseEntry, RuntimeMode, TaxonomyType
from emtrace.core.exceptions import DependencyError, StepValidationError
from emtrace.io.fasta import extract_fasta_by_ids, fasta_id_from_header, iter_fasta
from emtrace.models.annotation import read_annotation_tsv, write_annotation_tsv
from emtrace.steps import StepName, register_step
from emtrace.steps.base import PipelineStep
from emtrace.utils.entrez import EntrezFetch


# Allele frequency of the (single) ALT of a record, as a bcftools expression.
# The denominator is INFO/DP — the depth the caller actually used — not
# AD[0]+AD[1]: after ``bcftools norm -m -any`` splits a multiallelic site, the
# other ALT's reads are no longer in AD, so AD[0]+AD[1] overstates the frequency
# of every allele at that position.
_AF = "(INFO/AD[1])/(INFO/DP)"


def build_vcf_filters(config: AssemblyConfig) -> tuple[str, str]:
    """Build the (report, consensus) bcftools filter expressions.

    Two different products, so two different thresholds:

    * the **variant table** lists every allele reaching ``report_min_af``,
      which may include sub-dominant alleles (co-infection, intra-host
      variation);
    * the **consensus FASTA** holds one base per position, so only an allele
      strictly exceeding ``min_af`` — the dominant allele — is applied.

    Both require ``min_ac`` supporting reads. ``vcf_flt``, if set, is ANDed into
    both as an escape hatch for advanced use.

    Records are assumed to be biallelic (``bcftools norm -m -any``), so
    ``AD[1]`` is unambiguously *the* ALT of the record.
    """
    base = f"INFO/AD[1]>={config.min_ac}"
    report = f"{base} && {_AF}>={config.report_min_af}"
    consensus = f"{base} && {_AF}>{config.min_af}"
    if config.vcf_flt:
        report = f"({report}) && ({config.vcf_flt})"
        consensus = f"({consensus}) && ({config.vcf_flt})"
    return report, consensus
from emtrace.utils.paths import ensure_dir
from emtrace.utils.shell import run_shell, shell_path
from emtrace.wrappers.bcftools import Bcftools
from emtrace.wrappers.medaka import Medaka
from emtrace.wrappers.minimap2 import Minimap2
from emtrace.wrappers.samtools import Samtools

if TYPE_CHECKING:
    from emtrace.core.pipelinejob import PipelineJob


@register_step(StepName.ASS)
class AssemblyStep(PipelineStep):
    """Reference-based assembly.

    Reference genomes are sourced from ``config.assembly.refseqs`` (CSV of
    FASTA paths) and ``config.assembly.refaccs`` (CSV of NCBI accessions
    fetched via Entrez). All are concatenated into a single combined reference.
    If neither is set, the step is skipped.
    """

    step_name = "assembly"
    step_description = "Reference-based consensus assembly"

    def validate_inputs(self, sample: PipelineJob, **kwargs) -> bool:
        refseqs, refaccs = self._parse_refs()
        if not refseqs and not refaccs and not self.config.assembly.ref_auto.enabled:
            self.logger.warning("No refseqs/refaccs configured — assembly will be skipped")
            return True  # run() handles gracefully

        for path_str in refseqs:
            p = Path(path_str)
            if not p.exists():
                raise StepValidationError(self.step_name, f"Reference FASTA not found: {p}")

        # Polishing is on by default, so a missing medaka must not block assembly:
        # warn here; _polish_medaka then keeps the unpolished consensus.
        if self.config.assembly.medaka.enabled and not sample.is_paired:
            try:
                Medaka(path=self.config.assembly.medaka.path)
            except DependencyError as exc:
                self.logger.warning(
                    "medaka polishing is enabled but medaka is unavailable (%s) — the consensus "
                    "will not be polished. Install medaka, set assembly.medaka.path, or use "
                    "--no-medaka.", exc,
                )

        if sample.is_paired:
            r2 = sample.input_reads_r2
            if r2 is None or not Path(r2).exists():
                raise StepValidationError(
                    self.step_name,
                    f"Paired-end job missing R2 mate: {r2}",
                )
            return True

        reads = self._find_assembly_reads(sample)
        if reads is None:
            raise StepValidationError(
                self.step_name,
                "No reads found: no viral reads from annotation and no input reads available"
            )
        return True

    def run(self, sample: PipelineJob, **kwargs) -> PipelineJob:
        refseqs, refaccs = self._parse_refs()
        auto = self.config.assembly.ref_auto.enabled and not refseqs and not refaccs
        if self.config.assembly.ref_auto.enabled and not auto:
            self.logger.info("ref_auto ignored: explicit refseqs/refaccs are set")
        if not refseqs and not refaccs and not auto:
            self.logger.info("No refseqs/refaccs — skipping assembly")
            return sample

        reads, reads_r2 = self._resolve_assembly_reads(sample)
        if reads is None:
            self.logger.warning("No reads found for assembly — skipping")
            return sample

        asm_dir = sample.assembly_dir
        ensure_dir(asm_dir)
        ensure_dir(sample.log_path)

        self._ref_dbs: List[str] = []
        if auto:
            # candidates come from the viral partition when there is one: it is
            # what the annotation-time search (assembly_candidates.tsv) used
            refaccs = self._select_refs(sample, self._find_viral_reads(sample) or reads, asm_dir)
            if not refaccs:
                self.logger.warning("ref_auto: no reference selected — skipping assembly")
                return sample

        log_path = sample.log_path / f"{sample.name}.assemble_refbased.log"
        log_path.write_text("")

        # Build combined reference FASTA
        refgen = self._build_combined_ref(refseqs, refaccs, asm_dir, sample)

        config: AssemblyConfig = self.config.assembly
        num_threads = sample.num_threads

        # Output paths
        sam_file = asm_dir / "reads.refgen.sam"
        bam_file = asm_dir / "reads.refgen.bam"
        md_bam = asm_dir / "reads.refgen.md.bam"
        lowcov_bed = asm_dir / "reads.refgen.lowcovmask.bed"
        refgen_masked = asm_dir / "refgen.lowcovmasked.fa"
        vcf_file = asm_dir / "vcf.gz"
        consensus = asm_dir / "consensus.fa"

        # medaka outputs too, so a run without polishing leaves no stale ones
        polish_outputs = [asm_dir / n for n in (
            "consensus.bcftools.fa", "reads.polished.md.bam", "reads.polished.md.bam.bai",
            "reads.polished.lowcovmask.bed", "consensus.medaka.unmasked.fa")]
        for f in [sam_file, bam_file, md_bam, lowcov_bed, refgen_masked, vcf_file, consensus,
                  *polish_outputs]:
            if f.exists():
                f.unlink()

        # 1. Align reads to reference
        self.logger.info(
            "Aligning reads to reference (minimap2, %s)...",
            "paired-end" if reads_r2 is not None else "single-end",
        )
        self._align_to_ref(reads, refgen, sam_file, config, num_threads, reads_r2=reads_r2)

        # 2. SAM → sorted BAM → fixmate → markdup → index
        self.logger.info("Processing BAM (sort/fixmate/markdup)...")
        st = Samtools(num_threads=num_threads)
        st.sort_fixmate_markdup(sam_file, bam_file)
        if not md_bam.exists():
            candidate = bam_file.with_suffix(".md.bam")
            if candidate.exists() and candidate != md_bam:
                candidate.rename(md_bam)
                idx = candidate.with_suffix(".md.bam.bai")
                if idx.exists():
                    idx.rename(md_bam.with_suffix(".md.bam.bai"))
        if sam_file.exists():
            sam_file.unlink()

        # 2b. Guard: no reads mapped to the reference -> nothing to reconstruct.
        # Emit an all-N consensus (0% recovery) instead of falling through the
        # pileup, where an empty depth profile would leave the reference unmasked
        # and be mis-scored as 100% coverage. (Belt-and-suspenders with the
        # `samtools depth -aa` fix in depth_low_coverage_mask.)
        st = Samtools(num_threads=num_threads)
        if st.count_mapped(md_bam) == 0:
            self.logger.warning(
                "No reads mapped to reference — emitting all-N consensus (no recovery)"
            )
            self._write_all_n_consensus(refgen, consensus)
            lowcov_bed.write_text("")  # empty mask for provenance
            sample.set_file("assembly_bam", md_bam)
            sample.set_file("assembly_consensus", consensus)
            sample.set_file("assembly_lowcov_bed", lowcov_bed)
            return sample

        # 3. Low-coverage masking
        self.logger.info("Masking low-coverage regions (cov < %d)...", config.min_coverage)
        st = Samtools(num_threads=num_threads)
        st.depth_low_coverage_mask(bam=md_bam, output_bed=lowcov_bed, min_coverage=config.min_coverage)
        Samtools.maskfasta(fasta=refgen, bed=lowcov_bed, output=refgen_masked)

        # 4. Variant calling
        report_expr, consensus_expr = build_vcf_filters(config)
        self.logger.info(
            "Calling variants (bcftools); reporting alleles at AF>=%.2f, "
            "consensus at AF>%.2f, both with >=%d supporting reads",
            config.report_min_af, config.min_af, config.min_ac,
        )
        bc = Bcftools(num_threads=num_threads)
        vcf_raw = asm_dir / "vcf.raw.gz"
        preset = config.mpileup_preset
        if preset == "auto":
            preset = "ont-sup" if sample.platform == "nanopore" and not sample.is_paired else ""
        if preset and not bc.has_mpileup_preset(preset):
            self.logger.warning(
                "bcftools %s has no mpileup profile '%s' (needs bcftools >= 1.20); "
                "calling without it — indels will be largely missed on ONT data",
                bc.version(), preset,
            )
            preset = ""
        bc.mpileup_call(bam=md_bam, ref=refgen_masked, output_vcf=vcf_raw, preset=preset)
        bc.filter(input_vcf=vcf_raw, output_vcf=vcf_file, expression=report_expr, output_type="z")
        if vcf_raw.exists():
            vcf_raw.unlink()
        self._tabix_index(vcf_file)

        # 5. Consensus — only the dominant allele is applied, so the reported
        #    VCF may legitimately contain records absent from the FASTA.
        self.logger.info("Generating consensus sequence...")
        bc.consensus(
            vcf=vcf_file, ref=refgen_masked, output_fa=consensus, include=consensus_expr
        )

        if refgen_masked.exists():
            refgen_masked.unlink()
        fai = Path(str(refgen_masked) + ".fai")
        if fai.exists():
            fai.unlink()

        # 6. Optional medaka polishing (ONT single-end)
        if config.medaka.enabled:
            self._polish_medaka(
                sample, reads, reads_r2, refgen, vcf_file, consensus_expr,
                consensus, config, num_threads,
            )

        self.logger.info("Assembly complete -> %s", consensus)
        sample.set_file("assembly_bam", md_bam)
        sample.set_file("assembly_consensus", consensus)
        sample.set_file("assembly_vcf", vcf_file)
        sample.set_file("assembly_lowcov_bed", lowcov_bed)
        return sample

    def _polish_medaka(
        self,
        sample: PipelineJob,
        reads: Path,
        reads_r2: Optional[Path],
        refgen: Path,
        vcf_file: Path,
        consensus_expr: str,
        consensus: Path,
        config: AssemblyConfig,
        num_threads: int,
    ) -> None:
        """Polish *consensus* in place with medaka; keep the bcftools one aside.

        1. Draft = the dominant-allele VCF applied to the UNMASKED reference.
           The masked consensus is a poor draft: reads cannot align across its
           N runs, so the flanks of every masked region would lose support.
        2. The assembly reads (``assembly.input_reads``, default all host-
           filtered reads) are mapped to the draft with the assembly aligner
           settings (small seeds), so medaka sees the same reads as bcftools.
        3. ``medaka inference`` on the draft regions with depth >=
           ``min_coverage`` only, and ``medaka sequence --fill_char N``. The
           rest of the draft is still reference sequence (no variants were
           called there); polishing it lets thin, poorly aligned reads carry
           it into the consensus, so it is written as N instead.
        4. The reads are mapped to the polished sequence and positions below
           ``min_coverage`` are masked with N — in the polished coordinates,
           since medaka's indels shift them relative to the reference.

        On failure the bcftools consensus is kept and a warning is logged.
        """
        mcfg = config.medaka
        if reads_r2 is not None or sample.platform != "nanopore":
            self.logger.warning("medaka polishing is for ONT reads only — skipped")
            return

        asm_dir = consensus.parent
        bcf_consensus = asm_dir / "consensus.bcftools.fa"
        draft = asm_dir / "consensus.draft.fa"
        hdf = asm_dir / "medaka.hdf"
        polished = asm_dir / "consensus.medaka.unmasked.fa"
        pol_bam = asm_dir / "reads.polished.md.bam"
        pol_bed = asm_dir / "reads.polished.lowcovmask.bed"
        log = sample.log_path / f"{sample.name}.medaka.log"
        for f in (bcf_consensus, draft, hdf, polished, pol_bam, pol_bed):
            f.unlink(missing_ok=True)

        try:
            md = Medaka(num_threads=num_threads, path=mcfg.path)
            model = md.resolve_model(mcfg.model, reads, mcfg.fallback_model)

            # The VCF was called on the masked reference: a record overlapping a
            # masked position has N in REF and cannot be applied to the unmasked
            # one ("fasta sequence does not match the REF allele").
            Bcftools(num_threads=num_threads).consensus(
                vcf=vcf_file, ref=refgen, output_fa=draft,
                include=f'({consensus_expr}) && REF!~"N"',
            )
            self.logger.info("medaka: mapping reads to the draft...")
            self._map_sorted(reads, draft, asm_dir / "reads.draft.md.bam", config, num_threads)

            draft_bed = asm_dir / "reads.draft.lowcovmask.bed"
            Samtools(num_threads=num_threads).depth_low_coverage_mask(
                bam=asm_dir / "reads.draft.md.bam", output_bed=draft_bed,
                min_coverage=config.min_coverage,
            )
            covered = asm_dir / "reads.draft.covered.bed"
            if not self._complement_bed(draft_bed, draft, covered):
                self.logger.warning("medaka: no draft region reaches min_coverage — not polished")
                return

            self.logger.info("medaka: polishing (model %s)...", Path(model).name or "medaka default")
            md.inference(asm_dir / "reads.draft.md.bam", hdf, model, regions=covered,
                         batch_size=mcfg.batch_size, cpu=mcfg.cpu, log_file=log)
            md.sequence(hdf, draft, polished, fill_char="N")

            self.logger.info("medaka: re-masking positions with coverage < %d...", config.min_coverage)
            self._map_sorted(reads, polished, pol_bam, config, num_threads)
            Samtools(num_threads=num_threads).depth_low_coverage_mask(
                bam=pol_bam, output_bed=pol_bed, min_coverage=config.min_coverage
            )
            consensus.rename(bcf_consensus)
            Samtools.maskfasta(fasta=polished, bed=pol_bed, output=consensus)
        except Exception as exc:  # noqa: BLE001
            self.logger.warning("medaka polishing failed — keeping the bcftools consensus: %s", exc)
            if bcf_consensus.exists():
                bcf_consensus.replace(consensus)
            return
        finally:
            for f in (hdf, asm_dir / "reads.draft.md.bam", asm_dir / "reads.draft.md.bam.bai",
                      asm_dir / "reads.draft.lowcovmask.bed", asm_dir / "reads.draft.covered.bed",
                      draft, Path(f"{draft}.fai"), Path(f"{polished}.fai"),
                      Path(f"{polished}.gaps_in_draft_coords.bed")):
                f.unlink(missing_ok=True)

        self.logger.info("medaka: polished consensus -> %s (bcftools: %s)", consensus, bcf_consensus)
        sample.set_file("assembly_consensus_bcftools", bcf_consensus)
        sample.set_file("assembly_polish_bam", pol_bam)
        sample.set_file("assembly_polish_lowcov_bed", pol_bed)

    @staticmethod
    def _complement_bed(bed: Path, fasta: Path, out: Path) -> int:
        """Write the regions of *fasta* NOT in *bed* to *out*; return their count."""
        masked: dict[str, list[tuple[int, int]]] = {}
        for line in bed.read_text().splitlines():
            chrom, start, end = line.split("\t")[:3]
            masked.setdefault(chrom, []).append((int(start), int(end)))
        n = 0
        with open(out, "w", newline="\n") as fh:
            for header, seq in iter_fasta(fasta):
                chrom, pos = fasta_id_from_header(header), 0
                for start, end in sorted(masked.get(chrom, [])) + [(len(seq), len(seq))]:
                    if start > pos:
                        fh.write(f"{chrom}\t{pos}\t{start}\n")
                        n += 1
                    pos = max(pos, end)
        return n

    def _map_sorted(self, reads: Path, ref: Path, out_bam: Path,
                    config: AssemblyConfig, num_threads: int) -> Path:
        """Map *reads* to *ref* → sorted, deduplicated, indexed *out_bam*.

        *out_bam* must end in ``.md.bam``: sort_fixmate_markdup writes the
        sorted ``X.bam`` and derives ``X.md.bam`` from it.
        """
        bam = out_bam.parent / out_bam.name.replace(".md.bam", ".bam")
        sam = bam.with_suffix(".sam")
        self._align_to_ref(reads, ref, sam, config, num_threads)
        Samtools(num_threads=num_threads).sort_fixmate_markdup(sam, bam)
        for f in (sam, bam):
            f.unlink(missing_ok=True)
        return out_bam

    @staticmethod
    def _write_all_n_consensus(refgen: Path, output: Path) -> None:
        """Write an all-N FASTA mirroring the reference records (0% recovery).

        Used when no reads map to the reference, so the "consensus" carries the
        reference headers/lengths but no called bases.
        """
        records: list[tuple[str, int]] = []
        header: str | None = None
        seqlen = 0
        with open(refgen) as fh:
            for line in fh:
                if line.startswith(">"):
                    if header is not None:
                        records.append((header, seqlen))
                    header = line.rstrip("\n")
                    seqlen = 0
                else:
                    seqlen += len(line.strip())
            if header is not None:
                records.append((header, seqlen))
        with open(output, "w", newline="\n") as out:
            for hdr, length in records:
                out.write(hdr + "\n")
                out.write("N" * length + "\n")

    # ---------------------------------------------------------------
    # Reference resolution
    # ---------------------------------------------------------------

    def _select_refs(self, sample: PipelineJob, reads: Path, asm_dir: Path) -> List[str]:
        """Pick reference accessions automatically (``assembly.ref_auto``).

        Primary hits: *reads* searched against ``assembly.database``
        (``annotation/assembly_candidates.tsv``, reused when current).
        Secondary hits: annotation (``annot2.tsv``, else ``annot1.tsv``); also
        the only source when ``assembly.database`` is unset.  Writes
        ``ref_selection.tsv`` and remembers the annotation DBs of secondary
        picks, so their sequences can be fetched from those DBs.
        """
        ra = self.config.assembly.ref_auto
        db_name = self.config.assembly.database
        primary = None
        if ra.use_assembly_db:
            if db_name:
                primary = self._assembly_candidates(sample, reads, db_name)
            else:
                self.logger.warning("ref_auto: assembly.database is not set — using annotation hits only")
        secondary = None
        if ra.use_annotation or primary is None:
            secondary = self._annotation_hits(sample)
        for df in (primary, secondary):
            self._ensure_species(df)

        picks, table = select_references(primary, secondary, ra)
        out = asm_dir / "ref_selection.tsv"
        table.to_csv(out, sep="\t", index=False)
        sample.set_file("ref_selection", out)
        self._ref_dbs = sorted({p.db for p in picks if p.source == "secondary" and p.db})
        self.logger.info("ref_auto: selected %s (details: %s)",
                         ", ".join(p.accession for p in picks) or "none", out)
        return [p.accession for p in picks]

    def _assembly_candidates(self, sample: PipelineJob, reads: Path, db_name: str) -> Optional[pd.DataFrame]:
        """Hits of *reads* against the assembly DB, from cache or a fresh search."""
        path = sample.annotation_dir / "assembly_candidates.tsv"
        if candidates_are_current(path, reads):
            self.logger.info("ref_auto: reusing %s", path)
            return read_annotation_tsv(path)
        db_entry = self.config.databases.get(db_name)
        if db_entry is None:
            self.logger.warning("ref_auto: assembly.database '%s' not in config.databases", db_name)
            return None
        ra = self.config.assembly.ref_auto
        self.logger.info("ref_auto: searching %s against %s", reads.name, db_name)
        try:
            hits = search_reference_db(
                reads, db_entry, sample.num_threads, k=ra.minimap.k, w=ra.minimap.w,
                log_file=sample.log_path / f"minimap2.assembly_candidates.{db_name}.log",
            )
        except Exception as exc:  # noqa: BLE001
            self.logger.warning("ref_auto: assembly DB search failed: %s", exc)
            return None
        ensure_dir(sample.annotation_dir)
        write_annotation_tsv(hits, path)
        sample.set_file("assembly_candidates", path)
        return hits

    def _annotation_hits(self, sample: PipelineJob) -> Optional[pd.DataFrame]:
        """Annotation hits of the latest round (annot2, else annot1), or None."""
        for name in ("annot2.tsv", "annot1.tsv"):
            p = sample.annotation_dir / name
            if p.exists() and p.stat().st_size > 0:
                return read_annotation_tsv(p)
        self.logger.info("ref_auto: no annotation tables in %s", sample.annotation_dir)
        return None

    def _ensure_species(self, df: Optional[pd.DataFrame]) -> None:
        """Fill ``species`` from the NCBI taxonomy when the table lacks it."""
        if df is None or df.empty:
            return
        if "species" in df.columns and (df["species"].astype(str).str.strip() != "").any():
            return
        tax_db = self.config.taxonomy.source(TaxonomyType.NCBI).db
        if tax_db and Path(tax_db).exists():
            add_species_field(df, Path(tax_db), taxid_col="staxid")

    def _parse_refs(self) -> tuple[List[str], List[str]]:
        """Return (refseqs_list, refaccs_list) from config."""
        cfg = self.config.assembly
        refseqs = [s.strip() for s in cfg.refseqs.split(",") if s.strip()]
        refaccs = [a.strip() for a in cfg.refaccs.split(",") if a.strip()]
        return refseqs, refaccs

    def _build_combined_ref(
        self,
        refseqs: List[str],
        refaccs: List[str],
        asm_dir: Path,
        sample: PipelineJob,
    ) -> Path:
        """Concatenate refseqs + fetched refaccs into a single FASTA.

        Blank lines are dropped: NCBI efetch writes one after every record, and
        ``bcftools consensus`` stops at the first one, so the unmasked-reference
        draft of the medaka step kept only the first record of a multi-record
        reference (segmented genomes).
        """
        combined = asm_dir / "combined_ref.fa"

        def _write(text: str) -> None:
            for line in text.splitlines():
                if line.strip():
                    out.write(line.rstrip() + "\n")

        with open(combined, "w", newline="\n") as out:
            for path_str in refseqs:
                p = Path(path_str)
                self.logger.debug("Adding refseq: %s", p)
                _write(p.read_text())

            if refaccs:
                fetched = self._fetch_accessions(refaccs, asm_dir, sample)
                if fetched and fetched.exists():
                    _write(fetched.read_text())

        self.logger.info("Combined reference: %s (%d refseqs, %d refaccs)", combined, len(refseqs), len(refaccs))
        return combined

    def _fetch_accessions(
        self,
        refaccs: List[str],
        asm_dir: Path,
        sample: PipelineJob,
    ) -> Optional[Path]:
        """Resolve accessions to a FASTA file.

        Always tries the local ``assembly.database`` first (when configured).
        Behaviour when the accession is absent from the local DB:
          online  → fall back to NCBI Entrez with a warning.
          offline → log WARNING and abort (no network calls).
        When ``assembly.database`` is not configured, go straight to Entrez.
        """
        db_name = self.config.assembly.database
        offline = self.config.runtime.mode == RuntimeMode.OFFLINE
        extra_dbs = getattr(self, "_ref_dbs", [])

        if not db_name and not extra_dbs and offline:
            self.logger.warning(
                "Offline mode: assembly.database is not configured. "
                "Set it to a Minimap2 database with a FASTA sidecar (e.g. 'minimap.assembly')."
            )
            raise StepValidationError(
                self.step_name,
                "Offline mode: assembly.database is not set.\n"
                "Configure an assembly reference database (e.g. 'minimap.assembly') "
                "in the Assembly section of the Config editor, or switch to online mode.",
            )

        if db_name:
            local_result = self._try_fetch_from_local_db(refaccs, asm_dir, db_name)
            if local_result is not None:
                return local_result

        # ref_auto secondary picks: also try the annotation DBs that produced them
        if extra_dbs:
            got = self._try_fetch_from_dbs(refaccs, asm_dir, ([db_name] if db_name else []) + extra_dbs)
            if got is not None:
                return got

        if db_name or extra_dbs:
            db_label = ", ".join(([db_name] if db_name else []) + extra_dbs)
            # Accession(s) absent from local DB(s)
            missing = self._missing_accessions(refaccs, asm_dir, db_label)
            if offline:
                self.logger.warning(
                    "Offline mode: accession(s) %s not found in assembly DB '%s'. "
                    "Assembly aborted — select accessions listed in 'Pick from viral hits' "
                    "or switch to online mode.",
                    missing, db_label,
                )
                raise StepValidationError(
                    self.step_name,
                    f"Offline mode: {len(missing)} accession(s) not in '{db_label}': "
                    f"{missing}.\nUse 'Pick from viral hits' to select an accession "
                    "that is present in the local assembly database, or switch to online mode.",
                )
            # Online: fall back to Entrez
            self.logger.warning(
                "Accession(s) %s not in local assembly DB '%s'; falling back to Entrez.",
                missing, db_label,
            )

        return self._fetch_via_entrez(refaccs, asm_dir)

    def _try_fetch_from_local_db(
        self,
        refaccs: List[str],
        asm_dir: Path,
        db_name: str,
    ) -> Optional[Path]:
        """Try to extract *refaccs* from the local assembly DB FASTA.

        Returns the output Path if ALL accessions were found, otherwise None
        (so the caller can decide on the fallback strategy).
        """
        db_entry = self.config.databases.get(db_name)
        if db_entry is None:
            self.logger.warning("assembly.database '%s' not in config.databases.", db_name)
            return None

        fasta = self._resolve_db_fasta(db_entry)
        if fasta is None or not fasta.exists():
            self.logger.warning(
                "No FASTA sidecar found for assembly DB '%s' (path=%s). "
                "Ensure the .fa/.fasta file is installed beside the .mmi index.",
                db_name, db_entry.db,
            )
            return None

        out = asm_dir / "refaccs.fetched.fa"
        wanted = set(refaccs)
        n = extract_fasta_by_ids(fasta, out, wanted)
        if n < len(wanted):
            return None  # caller decides what to do with the miss

        self.logger.info(
            "Extracted %d/%d accession(s) from local DB '%s'.", n, len(wanted), db_name,
        )
        return out

    def _try_fetch_from_dbs(
        self,
        refaccs: List[str],
        asm_dir: Path,
        db_names: List[str],
    ) -> Optional[Path]:
        """Extract *refaccs* from the FASTA of several DBs, first DB first.

        Used for ref_auto picks from annotation hits, whose accessions may be
        missing from ``assembly.database``.  Returns the output Path if ALL
        accessions were found, otherwise None (the file keeps what was found).
        """
        out = asm_dir / "refaccs.fetched.fa"
        tmp = asm_dir / "refaccs.fetch.tmp.fa"
        wanted, found = set(refaccs), set()
        with open(out, "w", newline="\n") as fh:
            for name in db_names:
                need = wanted - found
                if not need:
                    break
                entry = self.config.databases.get(name)
                fasta = self._resolve_db_fasta(entry) if entry is not None else None
                if fasta is None:
                    self.logger.info("No FASTA for DB '%s'; cannot extract %s from it", name, sorted(need))
                    continue
                extract_fasta_by_ids(fasta, tmp, need)
                for header, seq in iter_fasta(tmp):
                    acc = fasta_id_from_header(header)
                    if acc in need and acc not in found:
                        fh.write(f"{header.rstrip()}\n{seq}\n")
                        found.add(acc)
        tmp.unlink(missing_ok=True)
        self.logger.info("Extracted %d/%d accession(s) from DB(s) %s.", len(found), len(wanted), db_names)
        return out if found == wanted else None

    def _missing_accessions(
        self,
        refaccs: List[str],
        asm_dir: Path,
        db_name: str,
    ) -> List[str]:
        """Return accessions from *refaccs* that are absent from the local DB."""
        out = asm_dir / "refaccs.fetched.fa"
        if not out.exists():
            return list(refaccs)
        present: set[str] = set()
        for header, _ in iter_fasta(out):
            present.add(fasta_id_from_header(header))
        return sorted(set(refaccs) - present)

    def _fetch_via_entrez(self, refaccs: List[str], asm_dir: Path) -> Optional[Path]:
        """Fetch accessions from NCBI Entrez."""
        rt = self.config.runtime
        if not rt.entrez_email:
            self.logger.error(
                "Cannot fetch accessions via Entrez: runtime.entrez_email is not set. "
                "Add 'entrez_email: your@email.com' under 'runtime:' in your YAML."
            )
            return None
        try:
            fetcher = EntrezFetch(email=rt.entrez_email, api_key=rt.entrez_api_key)
            out = asm_dir / "refaccs.fetched.fa"
            return fetcher.fetch_fasta(refaccs, out)
        except Exception as exc:
            self.logger.error("Entrez fetch failed: %s", exc)
            return None

    @staticmethod
    def _resolve_db_fasta(db_entry: DatabaseEntry) -> Optional[Path]:
        """Return the local FASTA path beside *db_entry*, or None if absent.

        Tries the configured path first; if it ends in ``.mmi``, also tries
        the ``.mmi``-stripped base, ``.fa``, and ``.fasta`` siblings.
        """
        raw = Path(db_entry.db)
        candidates: List[Path] = [raw]
        if raw.suffix == ".mmi":
            base = raw.with_suffix("")
            candidates.extend([base, base.with_suffix(".fa"), base.with_suffix(".fasta")])
        else:
            candidates.extend([raw.with_suffix(".fa"), raw.with_suffix(".fasta")])

        for c in candidates:
            if c.exists() and c.suffix.lower() in (".fa", ".fasta"):
                return c
        return None

    # ---------------------------------------------------------------
    # Internal helpers
    # ---------------------------------------------------------------

    def _resolve_assembly_reads(
        self, sample: PipelineJob,
    ) -> tuple[Optional[Path], Optional[Path]]:
        """Return ``(reads_r1, reads_r2)`` for the assembly step.

        Paired-end Illumina jobs always assemble directly from the raw
        ``input_reads`` / ``input_reads_r2`` pair (annotation steps do not
        run for PE input). Nanopore / single-end jobs use
        :meth:`_find_assembly_reads` (``assembly.input_reads``).
        """
        if sample.is_paired:
            r1 = Path(sample.input_reads) if sample.input_reads else None
            r2 = Path(sample.input_reads_r2) if sample.input_reads_r2 else None
            return r1, r2
        return self._find_assembly_reads(sample), None

    def _find_assembly_reads(self, sample: PipelineJob) -> Optional[Path]:
        """Reads to map, per ``assembly.input_reads``.

        ``bgfiltered`` (default): all host-filtered reads — registry
        ``host_filtered_reads``, else ``reads/reads.trim.hflt.fq.gz`` on disk,
        else the viral-reads chain's last resort (current / input reads).
        Mapping everything to the (small) reference makes the assembly
        independent of annotation sensitivity: reads too divergent for the
        annotation index still map with small seeds.
        ``viral``: the annotation's viral partition (:meth:`_find_viral_reads`).
        """
        if self.config.assembly.input_reads == "viral":
            return self._find_viral_reads(sample)
        candidates = [sample.get_file("host_filtered_reads"),
                      sample.reads_dir / "reads.trim.hflt.fq.gz"]
        for f in candidates:
            if f and Path(f).exists() and Path(f).stat().st_size > 0:
                self.logger.info("Assembling from host-filtered reads: %s", f)
                return Path(f)
        cr = sample.current_reads
        if cr and Path(cr).exists() and Path(cr).stat().st_size > 0:
            self.logger.info("No host-filtered reads — assembling from current reads: %s", cr)
            return Path(cr)
        return None

    def _find_viral_reads(self, sample: PipelineJob) -> Optional[Path]:
        """Locate reads for assembly.

        Priority:
          1. Viral reads from annotation (ann2 preferred, then ann1) via registry
          2. Viral reads from annotation found on disk (registry not populated)
          3. current_reads / input_reads — allows assembly without prior annotation
        """
        # Registry keys set by AnnotateStep at runtime or by _scan_existing_files on open
        registry_keys = ("reads_ann2_vi", "ann2_vi_reads", "reads_ann1_vi", "ann1_vi_reads")
        for key in registry_keys:
            f = sample.get_file(key)
            if f and f.exists() and f.stat().st_size > 0:
                return f
        # Disk fallback — handles cases where the registry was not populated
        rd = sample.reads_dir
        for rnd in ("ann2", "ann1"):
            for suffix in (".fq.gz", ".fq"):
                p = rd / f"reads.{rnd}.vi{suffix}"
                if p.exists() and p.stat().st_size > 0:
                    return p
        # Last resort: use whatever reads the pipeline currently has (or raw input)
        cr = sample.current_reads
        if cr and Path(cr).exists() and Path(cr).stat().st_size > 0:
            self.logger.info(
                "No viral reads found — assembling from current reads: %s", cr
            )
            return Path(cr)
        return None

    def _align_to_ref(self, reads: Path, refgen: Path, sam_out: Path,
                      config: AssemblyConfig, num_threads: int,
                      reads_r2: Optional[Path] = None) -> None:
        self._align_minimap2(reads, refgen, sam_out, config.minimap, num_threads, reads_r2=reads_r2)

    def _align_minimap2(self, reads: Path, refgen: Path, sam_out: Path,
                        mm_cfg, num_threads: int,
                        reads_r2: Optional[Path] = None) -> None:
        mm2 = Minimap2(num_threads=num_threads)
        # Short-read preset for paired-end Illumina; map-ont otherwise.
        preset = "sr" if reads_r2 is not None else "map-ont"
        stdout = mm2.align(
            db=refgen,
            query=reads,
            query2=reads_r2,
            preset=preset,
            m=mm_cfg.m,
            s=mm_cfg.s,
            c=mm_cfg.c,
            k=mm_cfg.k,
            w=mm_cfg.w,
            secondary="yes" if mm_cfg.secondary else "no",
            n=mm_cfg.n,
            sam_output=True,
        )
        sam_out.write_text(stdout, newline="\n")

    def _tabix_index(self, vcf: Path) -> None:
        run_shell(f"tabix -f {shell_path(vcf)}", tool_name="tabix")
