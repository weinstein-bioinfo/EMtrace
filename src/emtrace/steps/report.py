"""Reporting step — port of ``generate_reports()`` + merged statistics.

Pipeline:
  1. Build annotation table: nucl/prot split, dedup, union
  2. Link taxonomy (NCBI or SILVA lineage columns)
  3. Assign division + bphage fields
  4. Export annotation table to TSV + Excel
  5. Sort reads into taxonomy directory hierarchy
  6. Build abundance table (taxid × readn × bitscore)
  7. Merge EM probabilities (if available)
  8. Export abundance table to TSV + Excel
  9. FASTQ stage statistics + read-length histogram
 10. Per-division read counts, assembly BAM stats
 11. HTML stats report

Output (under ``{sample.reports_dir}/``):
    ``reports/annot_table.tsv``       — annotation table with lineage
    ``reports/annot_table.xlsx``      — annotation table Excel
    ``reports/abund_table.tsv``       — abundance table with lineage
    ``reports/abund_table.xlsx``      — abundance table Excel
    ``reports/stats_fastq.csv``       — per-stage FASTQ statistics
    ``reports/stats_div.csv``         — per-division read counts
    ``reports/stats_assembly.csv``    — per-BAM assembly statistics
    ``reports/stats_report.html``     — HTML report with tables + figures
    ``plots/qc/read_lengths.png``     — multi-panel read-length histogram

All abundance visualizations are produced by the separate PlotStep
(``steps/plots.py``).
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import TYPE_CHECKING, Dict, List, Optional, Set

import pandas as pd

from emtrace.analysis.tailtrim import resolve_count_col, tail_trim

from emtrace.analysis.taxonomy_utils import (
    add_bphage_field,
    add_division_field,
    filter_no_species,
    is_silva_annotation,
)
from emtrace.models.annotation import (
    ANNOT_LINEAGE_FIELDS,
    read_annotation_tsv,
    write_annotation_tsv,
)
from emtrace.models.abundance import (
    ABUNDANCE_FIELDS,
    write_abundance_tsv,
)
from emtrace.core.config import ReportConfig, TaxonomyType, TaxonomySourceConfig
from emtrace.steps import StepName, register_step
from emtrace.steps.base import PipelineStep
from emtrace.steps.stats import (
    _FASTQ_HEADERS, _DIV_HEADERS, _ASM_HEADERS,
    _fastq_stats, _read_lengths,
    _bam_stats,
    _division_counts,
    _plot_histograms,
    _generate_html_report,
    _YAML_AVAILABLE,
)
try:
    import yaml as _yaml
except ImportError:
    _yaml = None  # type: ignore[assignment]
from emtrace.io.tsv import top_hits_per_search
from emtrace.utils.paths import ensure_dir

#: NCBI division of viral taxa, as written by ``add_division_field``. SILVA
#: tables carry the domain in the same column and never hold this value, so
#: ``report.tail_viruses: false`` is vacuous under SILVA taxonomy.
VIRUSES_DIVISION = "Viruses"

if TYPE_CHECKING:
    from emtrace.core.pipelinejob import PipelineJob

# Search programs that produce nucleotide-level annotations
_NUCL_SEARCHES = {"minimap2", "blastn"}
# Search programs that produce protein-level annotations
_PROT_SEARCHES = {
    "diamond.blastx", "diamond.blastp", "diamondx", "diamondp", "blastx", "blastp"}


@register_step(StepName.REP)
class ReportStep(PipelineStep):
    """Generate annotation/abundance tables and Krona data.

    Port of ``generate_reports()`` from Lazypipe
    Replaces csvtk, taxonkit, and R scripts with pure Python/Pandas.
    """

    step_name = "report"
    step_description = "Reporting and abundance estimation"

    def _detect_tax_type_from_df(self, annot: pd.DataFrame) -> TaxonomyType:
        """Infer taxonomy type from the 'db' column of a loaded annotation TSV."""
        return (
            TaxonomyType.SILVA if is_silva_annotation(annot, self.config)
            else TaxonomyType.NCBI
        )

    def _resolve_annotation_tsv(self, job: PipelineJob) -> Optional[Path]:
        """Return the annotation TSV to report on, or None if none exists.

        Prefers annot2.tsv (round-2 results) when it is present and non-empty;
        falls back to annot1.tsv; returns None if neither is usable.
        """
        ann_dir = job.annotation_dir
        for name in ("annot2.tsv", "annot1.tsv"):
            p = ann_dir / name
            if p.exists() and p.stat().st_size > 0:
                return p
        return None

    def validate_inputs(
        self, sample: PipelineJob, **kwargs,
    ) -> bool:
        if self._resolve_annotation_tsv(sample) is None:
            self.logger.warning(
                "No annotation data — reporting will produce empty results"
            )
        return True

    def run(self, job: PipelineJob, **kwargs) -> PipelineJob:

        rep_cfg: ReportConfig = self.config.report

        reports_dir = job.reports_dir
        reports_subdir = reports_dir / "reports"
        ensure_dir(reports_dir)
        ensure_dir(reports_subdir)

        # When there are no annotations we still emit FASTQ QC stats
        # (stats_fastq.csv) — read QC is independent of database hits, and a
        # QC-only report keeps the sample visible to the cross-sample summary.
        ann_tsv = self._resolve_annotation_tsv(job)
        if ann_tsv is None:
            self.logger.warning("No annotation data — emitting QC stats only")
            self._run_stats(job, reports_subdir)
            return job

        # Read the raw annotation
        annot = read_annotation_tsv(ann_tsv)
        if annot.empty:
            self.logger.warning("Annotation file is empty — emitting QC stats only")
            self._run_stats(job, reports_subdir)
            return job

        # Detect taxonomy from the annotation file itself, not from global config.
        tax_type = self._detect_tax_type_from_df(annot)
        tax_cfg: TaxonomySourceConfig = (
            self.config.taxonomy.SILVA if tax_type == TaxonomyType.SILVA
            else self.config.taxonomy.NCBI
        )
        tax_db = Path(tax_cfg.db)
        self.logger.info(
            "Detected %s taxonomy from %s (db=%s, toptaxrank=%s)",
            tax_type.value, ann_tsv.name, tax_db, tax_cfg.toptaxrank,
        )

        # ----------------------------------------------------------
        # 1. Build annotation table (nucl/prot dedup + union)
        # ----------------------------------------------------------
        self.logger.info("Building annotation table (nucl/prot split + dedup)...")
        annot_table = self._build_annotation_table(annot, rep_cfg.min_qcov, rep_cfg.min_bitscore)
        if annot_table.empty:
            self.logger.warning("No annotations after filtering — emitting QC stats only")
            self._run_stats(job, reports_subdir)
            return job

        # ----------------------------------------------------------
        # 2. Link taxonomy
        # ----------------------------------------------------------
        if tax_type == TaxonomyType.NCBI:
            self._link_ncbi_taxonomy(annot_table, tax_db)
        elif tax_type == TaxonomyType.SILVA:
            self._link_silva_taxonomy(annot_table, tax_cfg)
            # "species" is the cleaned binomial set by add_lineage_to_df (empty
            # when organism_name names no species); the raw name stays in
            # "organism" so both taxonomies expose a uniform column.
            annot_table["species"] = annot_table.get("species", pd.Series("", index=annot_table.index)).fillna("")
            annot_table["organism"] = annot_table["organism_name"]

        # ----------------------------------------------------------
        # 3. Division + bphage
        # ----------------------------------------------------------
        if tax_type == TaxonomyType.NCBI:
            self.logger.info("Assigning division + bphage...")
            add_division_field(annot_table, tax_db)
            phfilter = self._find_phage_filter()
            add_bphage_field(annot_table, tax_db, phfilter_path=phfilter)

        # ----------------------------------------------------------
        # 3b. Drop taxa with no resolved binomial species (opt-in).
        #     Applied here, after lineage linking, so `species` exists.
        #     Only the reported table is affected — the raw annotation
        #     tables written by the annotate step keep every row.
        # ----------------------------------------------------------
        if rep_cfg.exclude_no_species:
            annot_table = filter_no_species(
                annot_table, context="annot_table",
            )
            if annot_table.empty:
                self.logger.warning(
                    "No annotations left after exclude_no_species — "
                    "emitting QC stats only"
                )
                self._run_stats(job, reports_subdir)
                return job

        # ----------------------------------------------------------
        # 4. Abundance table
        #
        #    Built BEFORE the annotation table is written, because
        #    `report.tail` is computed here and every other report follows
        #    the taxon set it leaves: annot_table, the sorted reads, the
        #    division stats and the plots all agree on what was detected.
        # ----------------------------------------------------------
        self.logger.info("Generating abundance table...")
        abund = self._build_abundance_table(annot, rep_cfg.min_qcov, rep_cfg.min_bitscore, job)

        if not abund.empty:
            # Link lineage
            if tax_type == TaxonomyType.NCBI:
                self._link_ncbi_lineage_abund( abund, tax_db )
                add_division_field( abund, tax_db, taxid_col="taxid" )
            elif tax_type == TaxonomyType.SILVA:
                self._link_silva_lineage_abund( abund, tax_cfg )

            # Independent of the annot_table filter above: the abundance
            # table is built from the raw `annot`, not from `annot_table`.
            if rep_cfg.exclude_no_species:
                abund = filter_no_species(abund, context="abund_table")

            # 4b. Cumulative-share tail trim, and its propagation. The raw
            #     annotation/annot*.tsv tables keep every row either way.
            trimmed = self._apply_tail_trim(abund, rep_cfg)
            if len(trimmed) != len(abund):
                annot_table = self._restrict_annot_to_taxa(annot_table, trimmed)
            abund = trimmed

            abund_out = reports_subdir / "abund_table.tsv"
            write_abundance_tsv(abund, abund_out)
            self.logger.info( "Abundance table: %s (%d taxa)", abund_out, len(abund) )
            job.set_file("abund_table", abund_out)

            if rep_cfg.export_excel:
                abund_xlsx = reports_subdir / "abund_table.xlsx"
                self._to_excel(abund, abund_xlsx)
                job.set_file("abund_table_xlsx", abund_xlsx)

        # ----------------------------------------------------------
        # 5. Write annotation table TSV + Excel
        # ----------------------------------------------------------
        annot_out = reports_subdir / "annot_table.tsv"
        write_annotation_tsv(annot_table, annot_out)
        self.logger.info("Annotation table: %s (%d rows)",  annot_out, len(annot_table))
        job.set_file("annot_table", annot_out)

        if rep_cfg.export_excel:
            annot_xlsx = reports_subdir / "annot_table.xlsx"
            self._to_excel(annot_table, annot_xlsx)
            job.set_file("annot_table_xlsx", annot_xlsx)

        # ----------------------------------------------------------
        # 6. Sort reads by taxonomy
        # ----------------------------------------------------------
        if rep_cfg.sort_reads:
            reads_gz = self._find_reads_for_sorting(job)
            if reads_gz and reads_gz.exists():
                requested = rep_cfg.toptaxrank
                if requested == "auto" or requested not in annot_table.columns:
                    if requested not in ("auto",) and requested not in annot_table.columns:
                        self.logger.warning(
                            "Requested toptaxrank '%s' not in annotation table — "
                            "falling back to auto-detected rank '%s'",
                            requested, tax_cfg.toptaxrank,
                        )
                    toptaxrank = tax_cfg.toptaxrank
                else:
                    toptaxrank = requested
                self.logger.info("Sorting reads by taxonomy (rank: %s)...", toptaxrank)
                self._sort_reads_by_taxonomy(
                    reads_gz, annot_table, job.reads_dir, toptaxrank)

        self.logger.info("Reporting complete")

        # ----------------------------------------------------------
        # 7-9. Statistics (merged from former StatsStep; sub-steps 9a-9d below)
        # ----------------------------------------------------------
        self._run_stats(job, reports_subdir)

        return job

    def _apply_tail_trim(
        self, abund: pd.DataFrame, rep_cfg: ReportConfig,
    ) -> pd.DataFrame:
        """Drop the cumulative tail from the abundance table (``report.tail``).

        Returns `abund` unchanged when the trim is off. Refuses to trim, with
        a warning, when the table carries no usable ``division`` column: a
        blind trim would remove exactly the rare viral taxa that
        ``tail_viruses: false`` exists to protect, and a silent one would be
        worse than none.

        Note the share is computed over every taxon, viral ones included --
        protection changes a group's threshold, not the denominator -- so on a
        virus-dominated library even a small `tail` can remove much of the
        non-viral background.
        """
        tail = float(getattr(rep_cfg, "tail", 0.0) or 0.0)
        if tail <= 0 or abund.empty:
            return abund

        divisions = (
            abund["division"].astype("string").fillna("")
            if "division" in abund.columns
            else pd.Series("", index=abund.index, dtype="string")
        )
        if not divisions.str.len().gt(0).any():
            self.logger.warning(
                "report.tail=%.4g requested, but the abundance table has no "
                "usable 'division' column (taxonomy linkage did not run?) — "
                "NOT trimming, to avoid dropping taxa blind.", tail,
            )
            return abund

        tail_viruses = bool(getattr(rep_cfg, "tail_viruses", False))
        count_col = resolve_count_col(abund)
        kept = tail_trim(
            abund, tail, count_col,
            group_col="division",
            group_tails=None if tail_viruses else {VIRUSES_DIVISION: 0.0},
        )
        dropped = len(abund) - len(kept)
        if dropped:
            counts = pd.to_numeric(abund[count_col], errors="coerce").fillna(0)
            lost = counts.sum() - pd.to_numeric(
                kept[count_col], errors="coerce").fillna(0).sum()
            self.logger.warning(
                "report.tail=%.4g dropped %d of %d taxa (%.0f %s reads, %.2f%% "
                "of the table)%s. The raw annotation/annot*.tsv tables keep "
                "every row.",
                tail, dropped, len(abund), lost, count_col,
                100.0 * lost / counts.sum() if counts.sum() else 0.0,
                "" if tail_viruses else "; viral taxa exempt",
            )
        return kept

    @staticmethod
    def _restrict_annot_to_taxa(
        annot_table: pd.DataFrame, abund: pd.DataFrame,
    ) -> pd.DataFrame:
        """Cut the annotation table down to the taxa the abundance table kept.

        This is what makes the reports agree after a tail trim: one taxon set
        for annot_table, the sorted reads and every plot. Joined on
        ``annot_table.staxid`` -> ``abund.taxid``, compared as strings because
        the two tables do not always agree on dtype.
        """
        if "staxid" not in annot_table.columns or "taxid" not in abund.columns:
            return annot_table
        keep = set(abund["taxid"].astype(str))
        return annot_table[annot_table["staxid"].astype(str).isin(keep)]

    def _run_stats(self, job: "PipelineJob", reports_subdir: Path) -> None:
        """Compute FASTQ/division/assembly stats and write reports."""
        sta_cfg = self.config.stats
        qc_dir = job.reports_dir / "plots" / "qc"
        ensure_dir(reports_subdir)
        ensure_dir(qc_dir)

        fastq_rows: List[Dict[str, object]] = []
        hist_panels: List[Tuple[str, List[int]]] = []

        # 9a. FASTQ stats at each pipeline stage
        for label, fq_path in self._collect_fastq_stages(job):
            self.logger.info("Computing stats: %s (%s)", label, fq_path.name)
            row = _fastq_stats(fq_path)
            row["stage"] = label
            fastq_rows.append(row)
            if sta_cfg.plot_histograms:
                lengths = _read_lengths(fq_path)
                if lengths:
                    hist_panels.append((label, lengths))

        # 9b. Division read counts. Both taxonomies write "division" into the
        # abundance table (SILVA's domain is copied into it), so one reader
        # serves both.
        abund_path = reports_subdir / "abund_table.tsv"
        div_rows = _division_counts(abund_path)
        if div_rows:
            self.logger.info("Division counts: %d divisions", len(div_rows))

        # 9c. Assembly BAM stats
        asm_rows = self._collect_assembly_stats(job)

        # 9d. Write CSV or YAML
        fmt = sta_cfg.output_format.lower()
        if fmt == "yaml" and _YAML_AVAILABLE and _yaml is not None:
            stats_path = job.reports_dir / "stats.yaml"
            data: Dict[str, object] = {}
            if fastq_rows:
                data["fastq_statistics"] = list(fastq_rows)
            if div_rows:
                data["division_read_counts"] = list(div_rows)
            if asm_rows:
                data["assembly_statistics"] = list(asm_rows)
            stats_path.parent.mkdir(parents=True, exist_ok=True)
            with open(stats_path, "w", encoding="utf-8") as fh:
                _yaml.dump(data, fh, default_flow_style=False, allow_unicode=True)
            self.logger.info("Stats YAML -> %s", stats_path)
            job.set_file("stats_yaml", stats_path)
        else:
            self._write_section_csv(fastq_rows, _FASTQ_HEADERS, reports_subdir / "stats_fastq.csv")
            self._write_section_csv(div_rows,   _DIV_HEADERS,   reports_subdir / "stats_div.csv")
            self._write_section_csv(asm_rows,   _ASM_HEADERS,   reports_subdir / "stats_assembly.csv")

        # 9e. Read-length histogram → plots/qc/
        fig_path: Optional[Path] = None
        if sta_cfg.plot_histograms and hist_panels:
            fig_path = qc_dir / "read_lengths.png"
            ok = _plot_histograms(hist_panels, fig_path)
            if ok:
                self.logger.info("Histogram -> %s", fig_path)
                job.set_file("stats_histogram", fig_path)
            else:
                self.logger.warning("matplotlib not available — skipping histogram")
                fig_path = None

        # 9e2. Stacked-bar composition plot → plots/ (canonical location; the
        # PLOT step regenerates the same file later). Embedded in the HTML
        # report so it stays self-contained even when run as `-p rep` alone.
        stacked_bar_path: Optional[Path] = self._stacked_bar_for_report(job, abund_path, qc_dir.parent)

        # 9f. HTML stats report → reports/
        html_path = reports_subdir / "stats_report.html"
        _generate_html_report(
            fastq_rows, div_rows, asm_rows, fig_path, html_path,
            stacked_bar_path=stacked_bar_path,
        )
        job.set_file("stats_html", html_path)
        self.logger.info("Stats report -> %s", html_path)

    def _stacked_bar_for_report(
        self, job: "PipelineJob", abund_path: Path, plots_dir: Path,
    ) -> Optional[Path]:
        """Generate the stacked-bar composition PNG for embedding in the HTML report.

        Returns the PNG path if produced, else None. Gated on
        ``config.report.stacked_bar``; matplotlib absence is handled inside
        ``plot_stacked_bar`` (PNG simply won't exist).
        """
        if not self.config.report.stacked_bar or not abund_path.exists():
            return None
        from emtrace.steps.plots import plot_stacked_bar

        abund = pd.read_csv(abund_path, sep="\t", dtype=str, keep_default_na=False)
        if abund.empty:
            return None
        for col in ("readn", "bitscore", "prob"):
            if col in abund.columns:
                abund[col] = pd.to_numeric(abund[col], errors="coerce").fillna(0)

        # "division" under either taxonomy: SILVA's domain is copied into it,
        # and the abundance table carries no "domain" column.
        plot_stacked_bar(
            abund, plots_dir, toptaxrank="division", label=job.name or "Sample",
        )
        png = plots_dir / "stacked_bar.png"
        return png if png.exists() else None

    # ==============================================================
    # Annotation table construction
    # ==============================================================

    def _infer_dbtype(self, search: str) -> str:
        """Map search program to dbtype (nucl/prot)."""
        s = search.lower().strip()
        if s in _NUCL_SEARCHES:
            return "nucl"
        if s in _PROT_SEARCHES:
            return "prot"
        # Default heuristic
        if "blast" in s and "p" in s:
            return "prot"
        return "nucl"

    def _build_annotation_table(
        self,
        annot: pd.DataFrame,
        min_qcov: float,
        min_bits: float,
    ) -> pd.DataFrame:
        """Nucl/prot split, dedup, union — port of Lazypipe lines 740-773.

        Perl logic:
          nucl: dbtype=="nucl" && qcov>=min_qcov && bitscore>=min_bits
                → dedup by (qseqid, staxid), keep top bitscore
          prot: dbtype=="prot" && qcov>=min_qcov && bitscore>=min_bits
                → dedup by (qseqid, orf, staxid), keep top bitscore
          union: nucl + prot, sorted by qseqid

        Rows are first reduced to each read's top-scoring hits per search
        (``top_hits_per_search``, ties kept), so the table lists per-read top
        classifications whether or not ``annotation.filter_tophits`` dropped
        sub-best hits at annotation. Selection precedes the qcov/bitscore
        filters, as when annotate filtered: a read whose top hit fails them
        drops out rather than being represented by a sub-best hit.
        """
        annot = top_hits_per_search(annot).copy()

        # Infer dbtype from 'search' column
        annot["_dbtype"] = annot["search"].map(self._infer_dbtype)

        # Ensure numeric columns
        # qcov stays NaN where unknown (e.g. minimap2 without -c): see _qcov_ok.
        annot["qcov"] = pd.to_numeric(annot["qcov"], errors="coerce")
        annot["bitscore"] = pd.to_numeric(annot["bitscore"], errors="coerce").fillna(0)
        qcov_ok = self._qcov_ok(annot, min_qcov)

        # --- Nucleotide annotations ---
        mask_n = (
            (annot["_dbtype"] == "nucl")
            & qcov_ok
            & (annot["bitscore"] >= min_bits)
        )
        nucl = annot.loc[mask_n].copy()
        if not nucl.empty:
            nucl = (
                nucl.sort_values(
                    ["qseqid", "staxid", "bitscore"],
                    ascending=[True, True, False],
                )
                .drop_duplicates( subset=["qseqid", "staxid"], keep="first" )
            )

        # --- Protein annotations ---
        mask_p = (
            (annot["_dbtype"] == "prot")
            & qcov_ok
            & (annot["bitscore"] >= min_bits)
        )
        prot = annot.loc[mask_p].copy()
        if not prot.empty:
            prot = (
                prot.sort_values(
                    ["qseqid", "orf", "staxid", "bitscore"],
                    ascending=[True, True, True, False],
                )
                .drop_duplicates( subset=["qseqid", "orf", "staxid"], keep="first" )
            )

        # --- Union ---
        combined = pd.concat( [nucl, prot], ignore_index=True )
        if not combined.empty:
            combined = combined.sort_values(["qseqid", "staxid"]).reset_index(drop=True)

        # Drop internal column
        combined.drop( columns=["_dbtype"], inplace=True, errors="ignore" )

        return combined

    # ==============================================================
    # Taxonomy linking
    # ==============================================================

    def _link_ncbi_taxonomy( self, df: pd.DataFrame, tax_db: Path) -> None:
        """Add NCBI lineage columns via bulk SQLite lookup.

        Replaces taxonkit reformat -f '{s}\t{g}\t{f}\t{o}\t{c}\t{k}'.
        """
        self.logger.info( "Linking NCBI taxonomy lineage...")
        from emtrace.models.taxonomy import NCBITaxonomy

        tax = NCBITaxonomy(tax_db)
        if not tax.db_path.exists():
            self.logger.warning( "NCBI taxonomy DB not found: %s", tax_db )
            for col in ANNOT_LINEAGE_FIELDS:
                df[col] = ""
            return

        # Bulk lookup
        import sqlite3
        with sqlite3.connect(tax.db_path) as conn:
            merged_map = dict(
                conn.execute(
                    "SELECT old_taxid, new_taxid "
                    "FROM merged"
                ).fetchall()
            )
            try:
                tax_rows = conn.execute(
                    "SELECT taxid, organism, species, genus, family, "
                    '"order", class, kingdom '
                    "FROM taxonomy"
                ).fetchall()
                _has_organism = True
            except sqlite3.OperationalError:
                self.logger.warning(
                    "NCBI taxonomy DB lacks 'organism' column — rebuild the DB "
                    "to enable strain/isolate-level annotation"
                )
                tax_rows = conn.execute(
                    "SELECT taxid, species, genus, family, "
                    '"order", class, kingdom '
                    "FROM taxonomy"
                ).fetchall()
                _has_organism = False
        lineage_map: Dict[int, Dict[str, str]] = {}
        for row in tax_rows:
            if _has_organism:
                lineage_map[row[0]] = {
                    "organism": row[1] or "",
                    "species": row[2] or "",
                    "genus": row[3] or "",
                    "family": row[4] or "",
                    "order": row[5] or "",
                    "class": row[6] or "",
                    "kingdom": row[7] or "",
                }
            else:
                lineage_map[row[0]] = {
                    "organism": "",
                    "species": row[1] or "",
                    "genus": row[2] or "",
                    "family": row[3] or "",
                    "order": row[4] or "",
                    "class": row[5] or "",
                    "kingdom": row[6] or "",
                }

        def _lookup_lineage(
            taxid_str: str,
        ) -> Dict[str, str]:
            empty = {r: "" for r in ANNOT_LINEAGE_FIELDS}
            try:
                tid = int(str(taxid_str).strip())
            except (ValueError, TypeError):
                return empty
            resolved = merged_map.get(tid, tid)
            return lineage_map.get(resolved, empty)

        lookups = df["staxid"].map(_lookup_lineage)
        for col in ANNOT_LINEAGE_FIELDS:
            df[col] = lookups.map(
                lambda d, c=col: d.get(c, "") if pd.notna(d) else "",
            )

    def _link_silva_taxonomy( self, df: pd.DataFrame, tax_cfg) -> None:
        """Add SILVA lineage columns via accession lookup."""
        self.logger.info( "Linking SILVA taxonomy lineage..." )
        from emtrace.models.taxonomy import (SILVATaxonomy)

        silva_db = Path(tax_cfg.db)
        silva = SILVATaxonomy(silva_db)
        silva.add_lineage_to_df(df)

    # ==============================================================
    # Abundance table
    # ==============================================================

    def _qcov_ok(self, annot: pd.DataFrame, min_qcov: float) -> pd.Series:
        """Rows passing ``report.min_qcov``; unknown qcov is not filtered.

        qcov is empty/NA when it cannot be computed, e.g. minimap2 without
        ``-c`` reports only approximate chain coordinates. Treating that as 0
        would drop every such hit whenever ``min_qcov > 0``, so those rows
        skip the coverage test instead.
        """
        unknown = annot["qcov"].isna()
        if min_qcov > 0 and unknown.any():
            self.logger.info(
                "min_qcov=%.2f not applied to %d hit(s) with unknown query "
                "coverage (e.g. minimap2 without -c)",
                min_qcov, int(unknown.sum()),
            )
        return unknown | (annot["qcov"] >= min_qcov)

    @staticmethod
    def _read_em_table(em_path: Optional[Path]) -> Optional[pd.DataFrame]:
        """Read ``EM.taxid_prob.tsv`` as ``taxid`` (str), ``prob`` and, when
        present, ``readn_em``; None if the file or its ``prob`` column is missing."""
        if not em_path or not Path(em_path).exists():
            return None
        em_df = pd.read_csv(em_path, sep="\t", dtype=str, keep_default_na=False)
        if "prob" not in em_df.columns:
            return None
        cols = ["taxid", "prob"] + (["readn_em"] if "readn_em" in em_df.columns else [])
        em_df = em_df[cols].copy()
        for col in cols[1:]:
            em_df[col] = pd.to_numeric(em_df[col], errors="coerce")
        return em_df

    def _build_abundance_table(
        self,
        annot: pd.DataFrame,
        min_qcov: float,
        min_bits: float,
        sample: Optional["PipelineJob"] = None,
        em_prob_path: Optional[Path] = None,
    ) -> pd.DataFrame:
        """Build abundance table from raw annotation.

        Port of Lazypipe lines 863-910, one taxon per read:
          1. Reduce to each read's top-scoring hits per search
             (``top_hits_per_search``, ties kept), then apply the
             qcov/bitscore filters, in the same order as the annotation table
          2. nucl: one row per read (``_one_taxon_per_read``)
          3. prot: sum bitscore over ORFs per qseqid+staxid, one taxon per
             read, excluding qseqids already covered by nucl
          4. Union → group by staxid: readn = count(qseqid), sum(bitscore)
          5. Merge EM ``prob`` and ``readn_em`` (expected reads) if available

        A read whose top hits tie across taxa is credited to the tied taxon
        with the highest EM posterior, then the lowest taxid, so ``readn`` is
        reproducible and follows abundance rather than row order.
        """
        em_path: Optional[Path]
        if em_prob_path is not None:
            em_path = em_prob_path
        elif sample is not None:
            em_path = sample.reports_dir / "EM.taxid_prob.tsv"
        else:
            em_path = None
        em_df = self._read_em_table(em_path)
        em_prob: Dict[str, float] = (
            dict(zip(em_df["taxid"], em_df["prob"])) if em_df is not None else {}
        )

        annot = top_hits_per_search(annot).copy()
        annot["_dbtype"] = annot["search"].map(self._infer_dbtype)
        # qcov stays NaN where unknown (e.g. minimap2 without -c): see _qcov_ok.
        annot["qcov"] = pd.to_numeric(annot["qcov"], errors="coerce")
        annot["bitscore"] = pd.to_numeric(annot["bitscore"], errors="coerce").fillna(0)
        qcov_ok = self._qcov_ok(annot, min_qcov)

        def _one_taxon_per_read(df: pd.DataFrame) -> pd.DataFrame:
            if df.empty:
                return df
            df = df.assign(
                _taxkey=df["staxid"].astype(str),
                _em_prob=df["staxid"].astype(str).map(em_prob).fillna(0.0),
            )
            df = df.sort_values(
                ["qseqid", "bitscore", "_em_prob", "_taxkey"],
                ascending=[True, False, False, True],
            )
            return (
                df.drop_duplicates(subset=["qseqid"], keep="first")
                .drop(columns=["_taxkey", "_em_prob"])
            )

        # --- Nucl subset ---
        mask_n = (
            (annot["_dbtype"] == "nucl")
            & qcov_ok
            & (annot["bitscore"] >= min_bits)
        )
        nucl = _one_taxon_per_read(annot.loc[mask_n, ["qseqid", "staxid", "bitscore"]].copy())
        nucl_qseqids: Set[str] = set(nucl["qseqid"]) if not nucl.empty else set()

        # --- Prot subset ---
        mask_p = (
            (annot["_dbtype"] == "prot")
            & qcov_ok
            & (annot["bitscore"] >= min_bits)
        )
        prot = annot.loc[mask_p, ["qseqid", "staxid", "bitscore"]].copy()
        if not prot.empty:
            # Sum bitscore over ORFs per qseqid+staxid, then one taxon per read
            prot = prot.groupby(["qseqid", "staxid"], as_index=False).agg({"bitscore": "sum"})
            prot = _one_taxon_per_read(prot)
            # Exclude qseqids covered by nucl
            prot = prot[~prot["qseqid"].isin(nucl_qseqids)]

        if (not nucl.empty or not prot.empty) and not em_prob:
            self.logger.info(
                "No EM probabilities: reads with tied top hits are credited to "
                "the lowest taxid"
            )

        # --- Union ---
        union = pd.concat(
            [nucl, prot], ignore_index=True,
        )
        if union.empty:
            return pd.DataFrame(columns=ABUNDANCE_FIELDS)

        # Group by taxid: count reads, sum bitscores
        abund = (
            union.groupby("staxid", as_index=False)
            .agg(
                readn=("qseqid", "nunique"),
                bitscore=("bitscore", "sum"),
            )
            .rename(columns={"staxid": "taxid"})
            .sort_values("readn", ascending=False)
            .reset_index(drop=True)
        )

        # --- Merge EM probabilities and expected read counts ---
        if em_df is not None:
            abund["taxid"] = abund["taxid"].astype(str)
            abund = abund.merge(em_df, on="taxid", how="left")
            abund["prob"] = abund["prob"].fillna(0)
            abund["readn_em"] = (
                abund["readn_em"].fillna(0) if "readn_em" in em_df.columns else "NA"
            )
        else:
            abund["prob"] = "NA"
            abund["readn_em"] = "NA"

        # Add placeholder lineage columns
        for col in [
            "species", "genus", "family",
            "order", "class", "division",
        ]:
            if col not in abund.columns:
                abund[col] = ""

        return abund

    def _link_ncbi_lineage_abund(
        self, df: pd.DataFrame, tax_db: Path,
    ) -> None:
        """Add NCBI lineage columns to abundance table."""
        from emtrace.models.taxonomy import NCBITaxonomy
        import sqlite3

        tax = NCBITaxonomy(tax_db)
        if not tax.db_path.exists():
            return

        with sqlite3.connect(tax.db_path) as conn:
            merged_map = dict(
                conn.execute(
                    "SELECT old_taxid, new_taxid "
                    "FROM merged"
                ).fetchall()
            )
            try:
                tax_rows = conn.execute(
                    'SELECT taxid, organism, species, genus, family, "order", class FROM taxonomy'
                ).fetchall()
                _has_organism = True
            except sqlite3.OperationalError:
                self.logger.warning(
                    "NCBI taxonomy DB lacks 'organism' column — rebuild the DB "
                    "to enable strain/isolate-level annotation"
                )
                tax_rows = conn.execute(
                    'SELECT taxid, species, genus, family, "order", class FROM taxonomy'
                ).fetchall()
                _has_organism = False
        lineage_map: Dict[int, Dict[str, str]] = {}
        for row in tax_rows:
            if _has_organism:
                lineage_map[row[0]] = {
                    "organism": row[1] or "",
                    "species": row[2] or "",
                    "genus": row[3] or "",
                    "family": row[4] or "",
                    "order": row[5] or "",
                    "class": row[6] or "",
                }
            else:
                lineage_map[row[0]] = {
                    "organism": "",
                    "species": row[1] or "",
                    "genus": row[2] or "",
                    "family": row[3] or "",
                    "order": row[4] or "",
                    "class": row[5] or "",
                }

        def _lookup(taxid_str: str) -> Dict[str, str]:
            empty = {
                "organism": "", "species": "", "genus": "", "family": "", "order": "", "class": "",
            }
            try:
                tid = int(str(taxid_str).strip())
            except (ValueError, TypeError):
                return empty
            resolved = merged_map.get(tid, tid)
            return lineage_map.get(resolved, empty)

        lookups = df["taxid"].map(_lookup)
        for col in [
            "organism", "species", "genus", "family", "order", "class",
        ]:
            df[col] = lookups.map(
                lambda d, c=col: d.get(c, "") if pd.notna(d) else "",
            )

    def _link_silva_lineage_abund(
        self, df: pd.DataFrame, tax_cfg,
    ) -> None:
        """Add SILVA lineage columns to abundance table."""
        from emtrace.models.taxonomy import (
            SILVATaxonomy,
            SILVA_TAXMAP_RANKS,
        )

        silva = SILVATaxonomy(Path(tax_cfg.db))
        if not silva.db_path.exists():
            return

        # A SILVA lineage taxid carries thousands of organism names, so no single
        # organism_name can label it. Species comes only from the species key the
        # abundance rows are grouped by; lineage-level rows get an empty species.
        silva_map = silva.lineage_by_taxid()

        for rank in ["species"] + list(SILVA_TAXMAP_RANKS):
            mapped = "division" if rank == "domain" else rank
            df[mapped] = df["taxid"].map(
                lambda t, r=rank: silva_map.get(str(t), {}).get(r, ""),
            )
        df["organism"] = df["species"]

    # ==============================================================
    # Read sorting by taxonomy
    # ==============================================================

    def _sort_reads_by_taxonomy(
        self,
        reads_gz: Path,
        annot_table: pd.DataFrame,
        reads_dir: Path,
        toptaxrank: str = "division",
    ) -> None:
        """Sort reads into taxonomy directory hierarchy.

        Port of ``sort_reads_bytaxa.pl``.
        Creates: reads_dir/{division}/{family}/{genus}/{species}.fq
        """
        from emtrace.io.fastq import iter_fastq

        # Build classification map: qseqid -> (div, fam, genus, species)
        required = [
            "qseqid", toptaxrank,
            "family", "genus", "species",
        ]
        missing = [
            c for c in required
            if c not in annot_table.columns
        ]
        if missing:
            self.logger.warning(
                "Cannot sort reads — missing columns: %s",
                missing,
            )
            return

        # For each qseqid, classify by majority taxonomy
        class_map = self._classify_reads(
            annot_table, toptaxrank,
        )
        if not class_map:
            self.logger.info(
                "No reads to sort by taxonomy"
            )
            return

        # Read all FASTQ records into memory
        seqs: Dict[str, tuple] = {}
        for header, seq, plus, qual in iter_fastq(
            reads_gz,
        ):
            rid = header.lstrip("@").split()[0]
            seqs[rid] = (header, seq, qual)

        # Write sorted reads. Each file is truncated on its first write in this
        # run, so running the report again does not duplicate its reads.
        sorted_count = 0
        started: Set[Path] = set()
        for rid, (div, fam, gen, sp) in class_map.items():
            if rid not in seqs:
                continue

            header, seq, qual = seqs[rid]

            # Sanitize directory names
            div_s = self._sanitize_dirname(div)
            fam_s = self._sanitize_dirname(fam)
            gen_s = self._sanitize_dirname(gen)
            sp_s = self._sanitize_dirname(sp)

            # Write at species level
            sp_dir = (
                reads_dir / div_s / fam_s / gen_s
            )
            ensure_dir(sp_dir)
            sp_fq = sp_dir / f"{sp_s}.fq"
            with open(sp_fq, "a" if sp_fq in started else "w", newline="\n") as fh:
                fh.write(
                    f"@{header.lstrip('@')}\n"
                    f"{seq}\n+\n{qual}\n"
                )
            started.add(sp_fq)
            sorted_count += 1

        self.logger.info(
            "Sorted %d reads into taxonomy dirs",
            sorted_count,
        )

    def _classify_reads(
        self,
        annot: pd.DataFrame,
        toptaxrank: str,
    ) -> Dict[str, tuple]:
        """Classify each qseqid by majority taxonomy.

        Port of classify_read() from sort_reads_bytaxa.pl.
        For each rank, if 1 unique value → use it,
        2 → join with '+', 3+ → 'mixture', 0 → 'unknown'.
        """
        result: Dict[str, tuple] = {}
        ranks = [toptaxrank, "family", "genus", "species"]

        grouped = annot.groupby("qseqid")
        for qseqid, group in grouped:
            labels: List[str] = []
            for rank in ranks:
                vals = set(
                    v for v in group[rank].dropna()
                    if str(v).strip()
                )
                if len(vals) == 0:
                    labels.append("unknown")
                elif len(vals) == 1:
                    labels.append(vals.pop())
                elif len(vals) == 2:
                    labels.append(
                        "+".join(sorted(vals))
                    )
                else:
                    labels.append("mixture")
            result[str(qseqid)] = tuple(labels)

        return result

    @staticmethod
    def _sanitize_dirname(name: str) -> str:
        """Replace chars unsafe for directory names."""
        if not name:
            return "unknown"
        return re.sub(r'[\s/\\()\[\]]+', '_', name)

    # ==============================================================
    # Helpers
    # ==============================================================

    def _find_reads_for_sorting(
        self, sample: PipelineJob,
    ) -> Path | None:
        """Find the reads file to sort by taxonomy."""
        # Prefer host-filtered, fall back through preprocessing outputs
        for key in (
            "host_filtered_reads",
            "chopped_reads",
            "porechop_reads",
            "trimmed_reads",
            "current_reads",
        ):
            f = sample.get_file(key)
            if f and f.exists():
                return f
        return sample.input_reads

    def _find_phage_filter(self) -> Path | None:
        """Locate the ICTV phage filter TSV."""
        # Check in package data or config
        import importlib.resources

        try:
            pkg = importlib.resources.files(
                "emtrace",
            )
            candidate = pkg / "data" / (
                "ICTV.bphage.filter.tsv"
            )
            with importlib.resources.as_file(
                candidate,
            ) as p:
                if p.exists():
                    return p
        except (TypeError, FileNotFoundError):
            pass
        return None

    def _to_excel(
        self, df: pd.DataFrame, path: Path,
    ) -> None:
        """Write DataFrame to Excel (.xlsx) via openpyxl."""
        try:
            df.to_excel(
                path, index=False, engine="openpyxl",
            )
            self.logger.info("Excel export: %s", path)
        except ImportError:
            self.logger.warning(
                "openpyxl not installed — "
                "skipping Excel export"
            )
        except Exception as exc:
            self.logger.warning(
                "Excel export failed: %s", exc,
            )

    def _write_section_csv(
        self,
        rows: list,
        headers: list,
        path: Path,
    ) -> None:
        """Write a flat CSV section; no-op when rows is empty."""
        import csv
        if not rows:
            return
        with open(path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(headers)
            for row in rows:
                writer.writerow([row.get(h, "") for h in headers])
        self.logger.info("Stats -> %s", path)

    def _collect_fastq_stages(self, sample: "PipelineJob") -> List:
        """Collect (label, path) for each existing FASTQ pipeline stage."""
        stages: List = []
        if sample.input_reads and sample.input_reads.exists():
            stages.append(("input", sample.input_reads))
        porechop = sample.reads_dir / "reads.porechop.fq.gz"
        if porechop.exists():
            stages.append(("porechop", porechop))
        chopped = sample.reads_dir / "reads.chopped.fq.gz"
        if chopped.exists():
            stages.append(("chopped", chopped))
        hflt = sample.reads_dir / "reads.trim.hflt.fq.gz"
        if hflt.exists():
            stages.append(("host_filtered", hflt))
        for f in sorted(sample.reads_dir.glob("hflt.*.reads.fq.gz")):
            label = f.name.replace("hflt.", "host_").replace(".reads.fq.gz", "")
            stages.append((label, f))
        for rnd in ("ann1", "ann2"):
            for tag in ("vi", "ph", "ab", "euk", "un"):
                fgz = sample.reads_dir / f"reads.{rnd}.{tag}.fq.gz"
                fq  = sample.reads_dir / f"reads.{rnd}.{tag}.fq"
                p = fgz if fgz.exists() else fq
                if p.exists():
                    stages.append((f"{rnd}_{tag}", p))
        return stages

    def _collect_assembly_stats(self, sample: "PipelineJob") -> list:
        """Collect BAM stats dicts for all deduplicated assemblies."""
        rows = []
        asm_dir = sample.assembly_dir
        if not asm_dir.exists():
            return rows
        for bam in sorted(asm_dir.glob("*.md.bam")):
            try:
                self.logger.info("BAM stats: %s", bam.name)
                rows.append(_bam_stats(bam, num_threads=sample.num_threads))
            except Exception as exc:
                self.logger.warning("Failed BAM stats for %s: %s", bam.name, exc)
        return rows
