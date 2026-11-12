"""Pydantic configuration models for EMtrace.

Ports the Perl config.yaml (YAML::Tiny) to typed, validated Python models
with environment variable expansion ($data, $databases, $hostgenomes).
"""

from __future__ import annotations

import importlib.resources
import logging
import os
import re
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import yaml
from pydantic import BaseModel, Field, field_validator, model_validator

logger = logging.getLogger("emtrace.core.config")


class TaxonomyType(str, Enum):
    """Supported taxonomy systems."""
    NCBI = "NCBI"
    SILVA = "SILVA"


class RuntimeMode(str, Enum):
    """Pipeline operating mode controlling fallback behaviour for assembly references.

    online:  try local assembly DB first; fall back to NCBI Entrez if the accession
             is absent from the DB.
    offline: try local assembly DB only; warn and abort if the accession is missing
             (no network calls are made).
    """
    ONLINE = "online"
    OFFLINE = "offline"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _deep_merge(base: Dict[str, Any], overlay: Dict[str, Any]) -> Dict[str, Any]:
    """Overlay `overlay` onto `base`, recursing into mappings only.

    Scalars and lists replace wholesale: a profile naming ``parasite_genera``
    or ``anns.options`` sets that list rather than appending to it, which is
    the only rule that stays predictable as the config grows. `base` is not
    modified.
    """
    out = dict(base)
    for key, value in overlay.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _deep_merge(out[key], value)
        else:
            out[key] = value
    return out


def _expand(value: str) -> str:
    """Expand $ENV_VARS and ~ in a string."""
    return os.path.expandvars(os.path.expanduser(value))


def _expand_path(value: str | None) -> Path | None:
    """Expand env vars and return a Path, or None."""
    if value is None or value == "" or value == "0":
        return None
    return Path(_expand(str(value)))


# ---------------------------------------------------------------------------
# Sub-models
# ---------------------------------------------------------------------------

class ChopperConfig(BaseModel):
    """Chopper quality/length filtering settings."""

    quality: int = 10                       # min Phred average quality score
    minlength: int = 400                    # min read length (bp)
    maxlength: Union[int, float] = float("inf")  # max read length, inf = no limit
    contam: Optional[str] = None            # path to contaminant FASTA
    trim_approach: str = "trim-by-quality"  # fixed-crop|trim-by-quality|best-read-segment|split-by-low-quality
    cutoff: int = 10                        # quality threshold for end trimming
    headcrop: int = 0                       # trim N bases from read start
    tailcrop: int = 0                       # trim N bases from read end

    @field_validator("maxlength", mode="before")
    @classmethod
    def parse_inf(cls, v: Any) -> Any:
        if isinstance(v, str) and v.upper() == "INF":
            return float("inf")
        # 0 in YAML means "no limit" (consistent with GUI spinbox convention)
        if v == 0:
            return float("inf")
        return v

    def to_cli_args(self) -> list[str]:
        """Build Chopper CLI arguments from config values."""
        args = [
            "--quality", str(self.quality),
            "--minlength", str(self.minlength),
            "--headcrop", str(self.headcrop),
            "--tailcrop", str(self.tailcrop),
        ]
        if self.maxlength != float("inf"):
            args.extend(["--maxlength", str(int(self.maxlength))])
        if self.contam:
            args.extend(["--contam", self.contam])
        if self.trim_approach:
            args.extend(["--trim-approach", self.trim_approach])
        if self.cutoff:
            args.extend(["--cutoff", str(self.cutoff)])
        return args


class PorechopABIConfig(BaseModel):
    """Porechop_ABI adapter removal settings."""

    # general
    verbosity: int = 1                      # output detail: 0=none, 1=some, 2=lots, 3=full
    # ab initio options
    ab_initio: bool = False                 # infer adapters from reads (recommended)
    custom_adapters: Optional[str] = None   # path to custom adapter file
    discard_database: bool = False          # ignore built-in adapter database
    window_size: int = 3                    # smoothing window for drop-cut algorithm
    # consensus mode options
    number_of_run: int = 10                 # sampling iterations for initial consensus
    consensus_run: int = 20                 # additional runs if no stable consensus found
    all_above_x: int = 10                   # min % of total adapters for consensus selection
    best_of_x: int = 0                      # keep only best N consensus seqs, 0 = all
    # adapter search
    adapter_threshold: float = 90.0         # min % identity to label adapter as present
    check_reads: int = 10000                # reads sampled for adapter detection
    scoring_scheme: str = "3,-6,-5,-2"      # alignment scores: match,mismatch,gap-open,gap-extend
    # end adapter trimming
    end_size: int = 150                     # bp searched at each read end for adapters
    min_trim_size: int = 4                  # ignore adapter alignments smaller than this
    extra_end_trim: int = 2                 # extra bases removed next to end adapters
    end_threshold: float = 75.0             # min % identity for end adapter removal
    # middle adapter settings
    no_split: bool = True                   # skip splitting reads at middle adapters
    discard_middle: bool = False            # discard reads with middle adapters
    middle_threshold: float = 90.0          # min % identity for middle adapter detection
    extra_middle_trim_good_side: int = 10   # extra bases removed on good side of middle adapter
    extra_middle_trim_bad_side: int = 100   # extra bases removed on bad side of middle adapter
    min_split_read_size: int = 1000         # min bp for post-split read pieces

    def to_cli_args(self) -> list[str]:
        """Build Porechop_ABI CLI arguments from config values."""
        args: list[str] = []
        if self.ab_initio:
            args.append("--ab_initio")
        if self.custom_adapters:
            # absolute: porechop_abi runs in the output dir (wrappers/porechop.py)
            args.extend(["--custom_adapters",
                         str(Path(self.custom_adapters).expanduser().resolve())])
        if self.discard_database:
            args.append("--discard_database")
        args.extend([
            "--verbosity", str(self.verbosity),
            "--window_size", str(self.window_size),
            "--number_of_run", str(self.number_of_run),
            "--consensus_run", str(self.consensus_run),
            "--all_above_x", str(self.all_above_x),
            "--adapter_threshold", str(self.adapter_threshold),
            "--check_reads", str(self.check_reads),
            "--scoring_scheme", self.scoring_scheme,
            "--end_size", str(self.end_size),
            "--min_trim_size", str(self.min_trim_size),
            "--extra_end_trim", str(self.extra_end_trim),
            "--end_threshold", str(self.end_threshold),
            "--middle_threshold", str(self.middle_threshold),
            "--extra_middle_trim_good_side", str(self.extra_middle_trim_good_side),
            "--extra_middle_trim_bad_side", str(self.extra_middle_trim_bad_side),
            "--min_split_read_size", str(self.min_split_read_size),
        ])
        if self.best_of_x > 0:
            args.extend(["--best_of_x", str(self.best_of_x)])
        if self.no_split:
            args.append("--no_split")
        if self.discard_middle:
            args.append("--discard_middle")
        return args


    def to_cli_args_porechop(self) -> list[str]:
        """Build Porechop CLI arguments from subset of Porechop_ABI config values."""
        args: list[str] = []
        args.extend([
            "--adapter_threshold", str(self.adapter_threshold),
            "--check_reads", str(self.check_reads),
            "--scoring_scheme", self.scoring_scheme,
            "--end_size", str(self.end_size),
            "--min_trim_size", str(self.min_trim_size),
            "--extra_end_trim", str(self.extra_end_trim),
            "--end_threshold", str(self.end_threshold),
            "--middle_threshold", str(self.middle_threshold),
            "--extra_middle_trim_good_side", str(self.extra_middle_trim_good_side),
            "--extra_middle_trim_bad_side", str(self.extra_middle_trim_bad_side),
            "--min_split_read_size", str(self.min_split_read_size),
            "--verbosity", str(self.verbosity)
        ])
        if self.no_split:
            args.append("--no_split")
        if self.discard_middle:
            args.append("--discard_middle")

        return args


class BgFilterMinimapConfig(BaseModel):
    """minimap2 alignment parameters for host filtering."""
    c: bool = False              # Generate DP aligment and CIGAR
    m: int = 40                  # -m Chain score threshold passed to minimap2
    s: int = 80                  # -s DP alignment score threshold passed to minimap2
    k: int = 0                   # -k k-mer size; 0 = minimap/index default (15)
    w: int = 0                   # -w minimizer window; 0 = minimap/index default (10)
    save_mmi_index: bool = False # build and reuse *.mmi index alongside the FASTA
    params: str = ""             # additional CLI flags (appended verbatim)


class BgFilterConfig(BaseModel):
    """Background (host) filtering settings."""

    bgfilter: str = "Homo_sapiens"      # CSV of host genome keys to filter against
    minimap: BgFilterMinimapConfig = Field(default_factory=BgFilterMinimapConfig)


class SubsampleConfig(BaseModel):
    """Read subsampling settings — fraction (0 < size ≤ 1) or absolute count (size > 1)."""
    active: bool = False
    size: float = 0.2


class PreprocessConfig(BaseModel):
    """Preprocessing pipeline settings (Porechop_ABI + Chopper + host filtering)."""

    subsample: SubsampleConfig = Field(default_factory=SubsampleConfig)
    porechop_abi: PorechopABIConfig = Field(default_factory=PorechopABIConfig)
    chopper: ChopperConfig = Field(default_factory=ChopperConfig)
    bgfilter: BgFilterConfig = Field(default_factory=BgFilterConfig)


class BlastnConfig(BaseModel):
    """BLASTN annotation search settings."""

    bitscore: int = 200                     # min bitscore to keep a hit
    evalue: float = 10                      # max e-value
    max_target_seqs: int = 5                # max target sequences per query
    save_hits: bool = False                 # write raw BLAST tabular to disk alongside annot1/annot2.tsv


class MinimapConfig(BaseModel):
    """Minimap2 annotation search settings."""

    s: int = 200                            # min DP alignment score (AS:i:, only effective when c=True)
    m: int = 100                            # -m chain score threshold (always applied)
    c: bool = False                         # -c: full DP alignment; enables AS:i: score and accurate qcov/scov
    f: float = 0.0002                       # -f ignore minimizers occurring in top FLOAT fraction (0.0002–1.0)
    k: int = 0                              # -k k-mer size (0 = use preset default, e.g. 15 for map-ont; try 11–13 for divergent viruses)
    w: int = 0                              # -w minimizer window size (0 = use preset default; try 5–7 when reducing k)
    evalue: float = 10                      # max e-value (not used by minimap2 directly)
    max_target_seqs: int = 5                # max target sequences per query
    preset: str = "map-ont"
    secondary: str = "yes"
    K: str = "500M"                         # -K minibatch size for mapping (bases); lower (e.g. "100M") to cut peak RAM on large-index searches
    save_hits: bool = False                 # write raw PAF to disk alongside annot1/annot2.tsv


class DiamondBlastpConfig(BaseModel):
    """Diamond BLASTP/BLASTX annotation search settings."""

    bitscore: int = 80                      # min bitscore to keep a hit (benchmarked chain2)
    evalue: float = 10                      # max e-value
    mode: str = "very-sensitive"            # sensitivity mode
    max_target_seqs: int = 5               # max target sequences per query
    save_hits: bool = False                 # write raw Diamond tabular to disk alongside annot1/annot2.tsv


class OrfipyConfig(BaseModel):
    """Orfipy ORF detection settings."""

    min_len: int = 90                       # min ORF length (nt; benchmarked chain2)
    max_len: Union[int, float] = float("inf")  # max ORF length (nt), INF = no limit
    table: int = 1                          # genetic code table (NCBI)
    partial_5: bool = False                 # allow ORFs missing start codon
    partial_3: bool = True                  # allow ORFs missing stop codon
    strand: str = "b"                       # strand: f=forward, r=reverse, b=both

    @field_validator("max_len", mode="before")
    @classmethod
    def parse_inf(cls, v: Any) -> Any:
        if isinstance(v, str) and v.upper() == "INF":
            return float("inf")
        return v


class AnnotationConfig(BaseModel):
    """Read annotation settings (search tools + annotation rounds)."""

    blastn: BlastnConfig = Field(default_factory=BlastnConfig)
    minimap: MinimapConfig = Field(default_factory=MinimapConfig)
    diamond_blastp: DiamondBlastpConfig = Field(default_factory=DiamondBlastpConfig)
    orfipy: OrfipyConfig = Field(default_factory=OrfipyConfig)
    ann1: str = "minimap.nt.abv"            # round 1 database key(s)
    ann2: str = ""                          # round 2 database key(s)
    anns: Optional[str] = None              # named annotation strategy
    anns_options: Dict[str, str] = Field(default_factory=dict)  # strategy presets
    # Keep only each read's top-scoring hits (ties retained) after every search.
    # False (default) keeps all hits so EM can weigh near-best alternatives;
    # read partitions and the annotation/abundance tables use top hits either way.
    filter_tophits: bool = False

    @field_validator("anns", mode="before")
    @classmethod
    def parse_none(cls, v: Any) -> Any:
        if v is None or (isinstance(v, str) and v.lower() == "none"):
            return None
        return v


class RefAutoMinimapConfig(BaseModel):
    """minimap2 seeding for the assembly-DB search of ``ref_auto``."""
    k: int = 0                   # -k k-mer size; 0 = index default. Non-zero bypasses the .mmi (slower)
    w: int = 0                   # -w minimizer window; 0 = index default


class RefAutoConfig(BaseModel):
    """Automatic reference selection (see ``analysis/ref_select.py``)."""

    enabled: bool = False                # pick references when refseqs/refaccs are empty
    min_reads: int = 10                  # min distinct reads for a species / an extra accession (segment)
    min_read_frac: float = 0.01          # min share of the source's viral reads for a species
    max_species: int = 3                 # co-infections: at most this many species
    max_refs_per_species: int = 8        # segments per species (greedy read cover)
    quality_metric: str = "qcov_identity"  # qcov_identity | pident_qcov | bitscore | bits_per_base | reads_won
    secondary_margin: float = 1.0        # how much better a secondary pick must align to replace a primary pick
    use_assembly_db: bool = True         # primary source: search assembly.database
    use_annotation: bool = True          # secondary source: annotation hits (annot2/annot1), nucleotide searches only
    minimap: RefAutoMinimapConfig = Field(default_factory=RefAutoMinimapConfig)

    @field_validator("quality_metric")
    @classmethod
    def _check_metric(cls, v: str) -> str:
        from emtrace.analysis.ref_select import QUALITY_METRICS
        if v not in QUALITY_METRICS:
            raise ValueError(f"ref_auto.quality_metric must be one of {QUALITY_METRICS}, got {v!r}")
        return v

    @model_validator(mode="after")
    def _check_sources(self) -> "RefAutoConfig":
        if not (self.use_assembly_db or self.use_annotation):
            raise ValueError("ref_auto: use_assembly_db and use_annotation cannot both be false")
        return self


class AssemblyMinimapConfig(BaseModel):
    """minimap2 settings for mapping reads to the assembly reference(s).

    Tuned for divergent samples: small seeds and minimap2's own score
    thresholds, since the reference is small (so -k/-w cost nothing) and
    reads 15-20% away from it fail default k=15 seeding and stricter -m/-s.
    """
    c: bool = True               # -c base-level alignment (needed for calling)
    m: int = 0                   # -m chain score threshold; 0 = minimap2 default (40 for map-ont)
    s: int = 0                   # -s DP alignment score threshold; 0 = minimap2 default (40)
    k: int = 11                  # -k k-mer size; 0 = preset default (15 for map-ont)
    w: int = 5                   # -w minimizer window; 0 = preset default (10 for map-ont)
    secondary: bool = False      # --secondary=yes; mpileup/depth skip secondary alignments anyway
    n: int = 5                   # -N max secondary alignments (with secondary: true)
    params: str = ""             # unused; kept so older YAML files still load


class AssemblyMedakaConfig(BaseModel):
    """Optional medaka polishing of the bcftools consensus (ONT only)."""

    enabled: bool = True        # polish the consensus with medaka; CLI --medaka / --no-medaka
    model: str = "auto"         # medaka model name/path; "auto" = from the read headers
    fallback_model: str = ""    # used by "auto" when the reads carry no basecaller model; "" = medaka default
    path: str = ""              # medaka executable (or its bin/ dir); "" = on PATH
    batch_size: int = 100       # medaka inference --batch_size
    cpu: bool = False           # run inference on the CPU (no GPU/MPS)

    @field_validator("path", mode="before")
    @classmethod
    def expand_env(cls, v: Any) -> Any:
        if isinstance(v, str):
            return _expand(v)
        return v


class AssemblyConfig(BaseModel):
    """Reference-based assembly settings."""

    minimap: AssemblyMinimapConfig = Field(default_factory=AssemblyMinimapConfig)
    refseqs: str = ""                   # CSV of FASTA paths for reference-based assembly
    refaccs: str = ""                   # CSV of Entrez accessions for reference-based assembly
    database: str = ""                  # databases.<key> used as ref source in offline mode (e.g. 'minimap.assembly')
    min_coverage: int = 10                   # min depth to include a position in consensus
    min_af: float = 0.5                      # min AF for an allele to enter the CONSENSUS (strict >)
    report_min_af: float = 0.2               # min AF for an allele to appear in the VARIANT TABLE (>=)
    min_ac: int = 5                          # min reads supporting an allele (both outputs)
    vcf_flt: str = ""                        # optional extra bcftools expression, ANDed into both filters
    # Expert setting, deliberately not in default.yaml: bcftools mpileup -X
    # profile. "auto" = ont-sup (ONT indel calling) for nanopore, none for
    # Illumina; "" = none. Override only for e.g. old R9 / fast-basecalled data.
    mpileup_preset: str = "auto"
    # Expert setting, not in default.yaml: which reads are mapped to the
    # reference(s). "bgfiltered" = all host-filtered reads (annotation-
    # independent: divergent reads the annotation misses still count);
    # "viral" = the annotation's viral partition (reads.ann{2,1}.vi.fq).
    input_reads: str = "bgfiltered"

    @field_validator("input_reads")
    @classmethod
    def _check_input_reads(cls, v: str) -> str:
        if v not in ("bgfiltered", "viral"):
            raise ValueError(f"assembly.input_reads must be 'bgfiltered' or 'viral', got {v!r}")
        return v
    ref_auto: RefAutoConfig = Field(default_factory=RefAutoConfig)  # automatic reference selection
    medaka: AssemblyMedakaConfig = Field(default_factory=AssemblyMedakaConfig)  # consensus polishing

    @model_validator(mode="after")
    def _check_af_tiers(self) -> "AssemblyConfig":
        """The consensus is built from the reported VCF, so the reporting
        threshold must be the looser of the two — otherwise a dominant allele
        could be filtered out of the variant table before it ever reaches
        ``bcftools consensus``, silently leaving the reference base in place."""
        if self.report_min_af > self.min_af:
            raise ValueError(
                f"assembly.report_min_af ({self.report_min_af}) must be <= "
                f"assembly.min_af ({self.min_af}): the consensus is derived "
                f"from the reported variants."
            )
        return self


class TaxonomySourceConfig(BaseModel):
    """Single taxonomy source (NCBI or SILVA)."""

    db: str = ""
    toptaxrank: str = "division"
    taxgroups: str = "all"
    name: str = ""
    url: str = ""
    # SILVA-specific
    taxmap_slv: Optional[str] = None
    tax_slv: Optional[str] = None
    # URL of tax_slv_ssu_<release>.txt.gz, SILVA's per-path rank declarations.
    # Required from release 144 on: without it the lineage parser falls back to
    # positional rank assignment, which 144 broke by inserting a kingdom level.
    tax_slv_url: Optional[str] = None

    model_config = {"extra": "ignore"}

    @field_validator("db", "taxmap_slv", "tax_slv", mode="before")
    @classmethod
    def expand_env(cls, v: Any) -> Any:
        if isinstance(v, str):
            return _expand(v)
        return v


class TaxonomyConfig(BaseModel):
    """Taxonomy sources (NCBI / SILVA).

    There is no global taxonomy switch. Each annotation round reports through
    the taxonomy of the database it searched (see
    ``EMtraceConfig.taxonomy_for_database``): SILVA databases are keyed by
    SILVA accession and only resolve in SILVA taxonomy, every other database
    resolves in NCBI taxonomy.
    """

    # Factories, not instances: an instance would expand $EMTRACE_DB_DIR once,
    # at import, before the user config has set it.
    NCBI: TaxonomySourceConfig = Field(default_factory=lambda: TaxonomySourceConfig(
        db="$EMTRACE_DB_DIR/taxonomy/ncbi/taxonomy.db",
        toptaxrank="division",
        name="NCBI taxonomy",
        url="ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
    ))
    SILVA: TaxonomySourceConfig = Field(default_factory=lambda: TaxonomySourceConfig(
        db="$EMTRACE_DB_DIR/taxonomy/silva/taxonomy.db",
        toptaxrank="domain",
        name="SILVA taxonomy",
        url="https://ftp.arb-silva.de/release_144/Exports/taxonomy/taxmap_slv_ssu_ref144.txt.gz",
        tax_slv_url="https://ftp.arb-silva.de/release_144/Exports/taxonomy/tax_slv_ssu_144.txt.gz",
    ))

    def source(self, tax_type: TaxonomyType) -> TaxonomySourceConfig:
        """Return the taxonomy source for *tax_type*."""
        return getattr(self, tax_type.value)


class DatabaseEntry(BaseModel):
    """Single annotation/host database definition."""

    accession: str = "NA"
    db: str
    name: str
    search: str  # Minimap2 | Blastn | Blastp | Diamond
    url: str = "NA"
    urls: List[str] = Field(default_factory=list)  # multiple source URLs (minimap2 multi-file)
    group: str = "general"  # display group: "general" | "16S" | user-defined
    # Taxonomy the database's accessions resolve in. Unset means inferred from
    # the database name ("silva" -> SILVA, else NCBI); set it for a
    # SILVA-derived database whose name does not say so.
    taxonomy: Optional[TaxonomyType] = None

    @field_validator("taxonomy", mode="before")
    @classmethod
    def normalise_taxonomy(cls, v: Any) -> Any:
        if isinstance(v, str):
            return v.strip().upper() or None
        return v

    @field_validator("db", mode="before")
    @classmethod
    def expand_db_path(cls, v: Any) -> Any:
        if isinstance(v, str):
            return _expand(v)
        return v

    @property
    def url_list(self) -> List[str]:
        """Return effective URL list: ``urls`` if set, else ``[url]`` if valid, else ``[]``."""
        if self.urls:
            return list(self.urls)
        if self.url and self.url not in ("NA", "0", ""):
            return [self.url]
        return []



class EMConfig(BaseModel):
    """EM abundance estimation parameters."""

    temp: float = 1.0
    max_iter: int = 1000
    log_L_threshold: float = 0.01
    exclude_no_species: bool = False
    # Two-pass EM (as Emu): >0 drops taxa with fewer expected reads after the
    # first EM and re-runs it, re-assigning their reads. 0 = one pass.
    prune_min_reads: float = 0.0


class ReportConfig(BaseModel):
    """Reporting step parameters."""

    # min query coverage (fraction 0-1) for annotation/abundance tables; the only
    # coverage filter in the pipeline. Rows with unknown qcov (e.g. minimap2
    # without -c) are kept, not treated as 0.
    min_qcov: float = 0.0
    min_bitscore: float = 0.0       # min bitscore for annotation/abundance tables
    # Cumulative-share tail trim, as a FRACTION: 0.01 drops the smallest taxa
    # that together hold the last 1% of reads. Ranks on readn_em when EM ran,
    # else readn. Computed on the abundance table and propagated to
    # annot_table, sorted reads, stats and plots, so every report agrees on
    # what was detected; reads of trimmed taxa count as unmapped in batch
    # roll-ups. The raw annotation/annot*.tsv tables are never filtered.
    # 0 disables. Unrelated to preprocess.chopper.tailcrop (read-end trimming).
    tail: float = 0.0
    # False keeps the viral division out of the trim -- what virus discovery
    # needs, since the cumulative tail is where rare true positives live.
    # Viral reads still count toward the share denominator either way, so on a
    # virus-dominated library a small `tail` can still remove much of the
    # non-viral background. No effect when tail is 0.
    tail_viruses: bool = False
    # drop taxa with no resolved binomial species from annot/abund tables
    # (independent of em.exclude_no_species, which filters the EM input)
    exclude_no_species: bool = False
    export_excel: bool = True       # export annot/abund tables as .xlsx (requires openpyxl)
    sort_reads: bool = True         # sort reads into taxonomy directory hierarchy
    toptaxrank: str = "auto"        # top taxonomy rank for read sorting; "auto" detects from annotation db column
    sankey: bool = True             # generate Sankey diagram in the plot step
    sunburst: bool = True           # generate sunburst chart in the plot step
    stacked_bar: bool = True        # generate stacked bar chart in the plot step
    heatmap_genus: bool = False     # generate genus-level heatmap in the plot step
    heatmap_species: bool = False   # generate species-level heatmap in the plot step
    # invertebrate genera classified as parasites/pathogens in `batch-abundance`
    # (splits the NCBI "Invertebrates" division). Empty = use the built-in default
    # in emtrace.batch.abundance.PARASITE_GENERA.
    parasite_genera: List[str] = Field(default_factory=list)


class StatsConfig(BaseModel):
    """Statistics step parameters."""

    output_format: str = "csv"          # output format: "csv" or "yaml"
    plot_histograms: bool = True        # generate read-length histogram plots (requires matplotlib)


class RuntimeConfig(BaseModel):
    """Runtime parameters."""

    num_threads: int = 8                    # threads for BLAST, Diamond
    output_root: str = "results"            # root directory for all output (-o/--outdir default)
    temp_dir: str = "$TMPDIR"               # temp file directory
    entrez_email: str = ""                  # email for NCBI Entrez API (required for accession fetch)
    entrez_api_key: str = ""                # optional NCBI API key (raises rate limit to 10 req/s)
    mode: RuntimeMode = RuntimeMode.ONLINE  # online: Entrez fallback; offline: local DB only

    @field_validator("temp_dir", mode="before")
    @classmethod
    def expand_temp_dir(cls, v: Any) -> Any:
        if isinstance(v, str):
            return _expand(v)
        return v

    @property
    def temp_dir_path(self) -> Path:
        return Path(self.temp_dir)


class OutputConfig(BaseModel):
    """Output/pack format settings."""

    fasta_wrap: bool = False                # false = single-line FASTA; true = wrapped
    fasta_line_width: int = 80              # characters per line when fasta_wrap is true
    pack_reads: bool = False                # include all reads (FASTQ) in packed *.tar.gz
    pack_reads_vi: bool = True              # include viral reads (FASTQ) in packed *.tar.gz
    pack_assembly_bam: bool = False         # include assembly BAM + index in packed *.tar.gz
    realtime_abund_table: bool = True       # write abund_table.tsv after each realtime batch



# ---------------------------------------------------------------------------
# Annotation database defaults and YAML parsers (module-level helpers)
# ---------------------------------------------------------------------------

#: Default annotation database entries keyed by database name.
#: ``db`` holds the default filesystem subpath, resolved at runtime under
#: ``$EMTRACE_DB_DIR``.  URLs are read from the ``url`` / ``urls`` fields
#: directly on each entry in the ``databases:`` YAML section.
_DEFAULT_DB_ENTRIES: Dict[str, DatabaseEntry] = {
    "minimap.refseq.abv":    DatabaseEntry(name="minimap.refseq.abv",    search="Minimap2", db="minimap.refseq.abv.fa",         url=""),
    "minimap.refseq.abv.1g": DatabaseEntry(name="minimap.refseq.abv.1g", search="Minimap2", db="minimap.refseq.abv.1g",         url=""),
    "minimap.refseq.vi":     DatabaseEntry(name="minimap.refseq.vi",     search="Minimap2", db="minimap.refseq.vi.fa",          url=""),
    "blastn.refseq.abv":     DatabaseEntry(name="blastn.refseq.abv",     search="Blastn",   db="blastn.refseq.abv",             url=""),
    "blastn.refseq.abv.1g":  DatabaseEntry(name="blastn.refseq.abv.1g",  search="Blastn",   db="blastn.refseq.abv.1g",          url=""),
    "blastn.refseq.vi":      DatabaseEntry(name="blastn.refseq.vi",      search="Blastn",   db="blastn.refseq.vi",              url=""),
    "diamond.nrclust.vi":    DatabaseEntry(name="diamond.nrclust.vi",    search="Diamond",  db="nr_cluster_seq.vi.dmnd",        url=""),
    "diamond.uniref.vi":     DatabaseEntry(name="diamond.uniref.vi",     search="Diamond",  db="diamond.uniref100.vi.dmnd",     url=""),
    "minimap.silva":         DatabaseEntry(name="minimap.silva",         search="Minimap2", db="minimap.SILVA_SSURef.fasta",    url="", group="16S", taxonomy=TaxonomyType.SILVA),
    "blastn.silva":          DatabaseEntry(name="blastn.silva",          search="Blastn",   db="blastn.SILVA_SSURef",           url="", group="16S", taxonomy=TaxonomyType.SILVA),
    "minimap.assembly":      DatabaseEntry(name="minimap.assembly",      search="Minimap2", db="minimap.refseq.vi.234.fa",      url="https://a3s.fi/swift/v1/lazypipe-db/2026/minimap.refseq.vi.234.tar.gz", group="assembly"),
}


def _make_default_databases() -> Dict[str, DatabaseEntry]:
    """Return a fresh copy of the default annotation database registry."""
    return dict(_DEFAULT_DB_ENTRIES)


def _parse_databases_section(section: Dict[str, Any]) -> Dict[str, DatabaseEntry]:
    """Parse the ``databases:`` YAML section into a ``Dict[str, DatabaseEntry]``.

    When the YAML section is present and non-empty, it is **authoritative**:
    only the keys listed in YAML appear in the result.  Defaults fill in
    missing fields (url, search, db) for entries that exist in both places.
    When the YAML section is absent or empty, the full defaults are returned.
    """
    defaults = _make_default_databases()

    if not section:
        return defaults

    databases: Dict[str, DatabaseEntry] = {}
    for key, val in section.items():
        if not isinstance(val, dict):
            continue
        default = defaults.get(key)
        url = str(val.get("url", default.url if default else ""))
        urls_raw = val.get("urls", [])
        urls: List[str] = list(urls_raw) if isinstance(urls_raw, list) else []
        databases[key] = DatabaseEntry(
            name=key,
            search=str(val.get("search", default.search if default else key.split(".")[0])),
            db=str(val.get("db", default.db if default else key)),
            url=url,
            urls=urls,
            group=str(val.get("group", default.group if default else "general")),
            taxonomy=val.get("taxonomy", default.taxonomy if default else None),
        )

    return databases


_ACCESSION_RE = re.compile(r"(GC[AF]_\d+\.\d+)")


def _accession_from_db_path(db: str) -> str:
    """Extract a GenBank/RefSeq accession from a database path filename.

    Matches ``GCF_XXXXXXXXX.V`` or ``GCA_XXXXXXXXX.V`` in the filename.
    Returns ``"NA"`` when no accession is found.
    """
    m = _ACCESSION_RE.search(Path(db).name)
    return m.group(1) if m else "NA"


# Settings profiles: "default" is default.yaml itself; every other profile is a
# bundled overlay (config/<name>.yaml). Aliases are accepted wherever a profile
# name is (CLI --profile, the config's `profile` key).
PROFILE_ALIASES = {"virome": "default"}


def available_profiles() -> list[str]:
    """``"default"`` followed by the bundled overlay profiles (e.g. ``"16S"``)."""
    config_dir = importlib.resources.files("emtrace") / "config"
    overlays = sorted(
        f.name[: -len(".yaml")]
        for f in config_dir.iterdir()
        if f.name.endswith(".yaml") and not f.name.startswith("default")
    )
    return ["default", *overlays]


def canonical_profile(name: str) -> str:
    """Resolve a profile name or alias, case-insensitively: ``virome`` → ``default``."""
    from emtrace.core.exceptions import ConfigurationError

    key = PROFILE_ALIASES.get(name.lower(), name).lower()
    profiles = available_profiles()
    for profile in profiles:
        if profile.lower() == key:
            return profile
    aliases = ", ".join(f"{a} = {p}" for a, p in PROFILE_ALIASES.items())
    raise ConfigurationError(
        f"Unknown profile '{name}'. Available: {', '.join(profiles)} (alias: {aliases})"
    )


def _preset_profile(anns: str) -> str:
    """The profile an annotation preset belongs to: ``16S.*`` → 16S, else default."""
    return "16S" if anns.startswith("16S.") else "default"


def _leaf_paths(d: Dict[str, Any], prefix: tuple = ()) -> list[tuple]:
    """Key paths of the non-mapping values in a nested dict."""
    paths = []
    for k, v in d.items():
        if isinstance(v, dict):
            paths.extend(_leaf_paths(v, prefix + (k,)))
        else:
            paths.append(prefix + (k,))
    return paths


def _get_path(obj: Any, path: tuple) -> Any:
    for k in path:
        obj = getattr(obj, k)
    return obj


def _set_path(obj: Any, path: tuple, value: Any) -> None:
    for k in path[:-1]:
        obj = getattr(obj, k)
    setattr(obj, path[-1], value)


def _drop_bwa(section: Dict[str, Any], where: str) -> Dict[str, Any]:
    """Strip the retired bwa settings (``aligner``, ``bwa_mem``) from a config section.

    Configs written before EMtrace 1.0.1 still carry them; they load, with a
    warning, and minimap2 is used.
    """
    section = dict(section)
    aligner = section.pop("aligner", None)
    bwa_mem = section.pop("bwa_mem", None)
    if aligner is not None and str(aligner).lower().startswith("bwa"):
        logger.warning(
            "%s.aligner: bwa is not supported any more; using minimap2", where,
        )
    elif bwa_mem is not None:
        logger.warning(
            "%s.bwa_mem: bwa is not supported any more; settings ignored", where,
        )
    return section


def _parse_host_databases_section(section: Dict[str, Any]) -> Dict[str, DatabaseEntry]:
    """Parse the ``host.databases:`` YAML section into a ``Dict[str, DatabaseEntry]``.

    Host entries only require ``db`` (and optionally ``url``).  ``name``
    defaults to the YAML key (e.g. ``"Homo_sapiens"``) and ``search``
    defaults to ``"Minimap2"`` — the aligner used for host read filtering.
    The ``accession`` field is auto-extracted from the ``db`` path filename
    when not explicitly supplied (looks for ``GCF_``/``GCA_`` patterns).
    """
    host_databases: Dict[str, DatabaseEntry] = {}
    for key, val in section.items():
        if key == "sources" or not isinstance(val, dict):
            continue
        val = dict(val)
        if str(val.get("search", "")).lower() == "bwa":
            logger.warning(
                "host.databases.%s: bwa is not supported any more; entry skipped", key,
            )
            continue
        val.setdefault("name", key)
        val.setdefault("search", "Minimap2")
        entry = DatabaseEntry(**val)
        if entry.accession == "NA" and entry.db:
            entry = entry.model_copy(update={"accession": _accession_from_db_path(entry.db)})
        host_databases[key] = entry
    return host_databases


def _resolve_db_path(entry: DatabaseEntry) -> Optional[Path]:
    """Return the absolute path for a ``DatabaseEntry``.

    If ``entry.db`` is relative, it is resolved under ``$EMTRACE_DB_DIR``.
    Returns ``None`` when the path is relative and the env var is unset.
    """
    db = _expand(entry.db)
    p = Path(db)
    if p.is_absolute():
        return p
    db_dir = os.environ.get("EMTRACE_DB_DIR", "")
    if db_dir:
        return Path(db_dir) / db
    return None


# ---------------------------------------------------------------------------
# Top-level config
# ---------------------------------------------------------------------------

class EMtraceConfig(BaseModel):
    """Top-level configuration for EMtrace.

    Mirrors the Perl config.yaml structure with Pydantic validation.
    """

    # Preprocessing (Porechop_ABI + Chopper + host filtering)
    preprocess: PreprocessConfig = Field(default_factory=PreprocessConfig)

    # Annotation (search tool settings + round config)
    annotation: AnnotationConfig = Field(default_factory=AnnotationConfig)

    # Reference-based assembly
    assembly: AssemblyConfig = Field(default_factory=AssemblyConfig)

    # Default pipeline options
    logs: str = "logs"
    res: str = "results"
    pipe: str = "main"

    # Runtime
    runtime: RuntimeConfig = Field(default_factory=RuntimeConfig)

    # Output/pack format
    output: OutputConfig = Field(default_factory=OutputConfig)

    # Taxonomy
    taxonomy: TaxonomyConfig = Field(default_factory=TaxonomyConfig)

    # Annotation databases: key → DatabaseEntry (defaults + YAML overrides)
    databases: Dict[str, DatabaseEntry] = Field(
        default_factory=_make_default_databases,
    )

    # Host genome databases: key → DatabaseEntry (from host.databases YAML)
    host_databases: Dict[str, DatabaseEntry] = Field(
        default_factory=dict,
    )

    # Annotation strategies — keyed by strategy name
    ann_strategies: Dict[str, str] = Field(default_factory=dict)

    # EM abundance estimation
    em: EMConfig = Field(default_factory=EMConfig)

    # Reporting
    report: ReportConfig = Field(default_factory=ReportConfig)

    # Statistics
    stats: StatsConfig = Field(default_factory=StatsConfig)

    # Settings profile the config was built from: "default" or a bundled
    # overlay such as "16S" (CLI --profile, GUI Profile switch). None = not
    # recorded (configs written before 1.0.2); active_profile then infers it.
    profile: Optional[str] = None

    # Database root directory — persists EMTRACE_DB_DIR across sessions.
    # When set, applied to the env var at config load time (shell export wins).
    db_dir: Optional[str] = None

    # ------------------------------------------------------------------ #
    # Convenience accessor — bgfilter lives inside preprocess
    # ------------------------------------------------------------------ #

    @property
    def bgfilter(self) -> BgFilterConfig:
        """Shortcut: ``config.bgfilter`` → ``config.preprocess.bgfilter``."""
        return self.preprocess.bgfilter

    # ------------------------------------------------------------------ #
    # Validators
    # ------------------------------------------------------------------ #

    @field_validator("profile")
    @classmethod
    def _check_profile(cls, v: Optional[str]) -> Optional[str]:
        return canonical_profile(v) if v else None

    @model_validator(mode="after")
    def _apply_db_dir(self) -> "EMtraceConfig":
        """Set EMTRACE_DB_DIR from db_dir if not already in the environment."""
        if self.db_dir:
            db_dir = _expand(self.db_dir)
            os.environ.setdefault("EMTRACE_DB_DIR", db_dir)
        return self

    # ------------------------------------------------------------------ #
    # Constructors
    # ------------------------------------------------------------------ #

    @classmethod
    def from_yaml(cls, path: Path | str) -> EMtraceConfig:
        """Load configuration from a YAML file."""
        path = Path(path)
        with open(path) as fh:
            raw = yaml.safe_load(fh)
        return cls._from_raw(raw)

    @classmethod
    def create_default(cls) -> EMtraceConfig:
        """Load the bundled default config (the ``default`` profile)."""
        cfg = cls._from_raw(cls._default_raw())
        cfg.profile = "default"
        return cfg

    @classmethod
    def from_profile(cls, name: str) -> EMtraceConfig:
        """Load a bundled profile, overlaid on ``default.yaml``.

        A profile is a small YAML file under ``emtrace/config/`` holding
        only the keys that differ — see ``16S.yaml``. It is deep-merged onto
        the defaults, so a profile never has to restate the database,
        taxonomy or host blocks.

        This is deliberately NOT how ``--config`` behaves. ``from_yaml``
        replaces the defaults outright, and anything a config file omits falls
        back to the *model* defaults; changing that would silently alter what
        every existing full config inherits. Profiles are the additive path,
        config files are the exact one.

        Merge rule: mappings merge key by key; scalars and lists replace
        wholesale, so a profile that names ``parasite_genera`` or
        ``anns.options`` sets the whole list rather than appending to it.

        *name* may be an alias (``virome`` = ``default``); the result records
        the canonical name in ``profile``.
        """
        name = canonical_profile(name)
        cfg = cls._from_raw(_deep_merge(cls._default_raw(), cls._profile_raw(name)))
        cfg.profile = name
        return cfg

    @classmethod
    def _profile_raw(cls, name: str) -> Dict[str, Any]:
        """A bundled profile's overlay, parsed but not merged (``{}`` for default)."""
        name = canonical_profile(name)
        if name == "default":
            return {}
        profile_yaml = importlib.resources.files("emtrace") / "config" / f"{name}.yaml"
        with importlib.resources.as_file(profile_yaml) as path:
            with open(path) as fh:
                return yaml.safe_load(fh) or {}

    # ------------------------------------------------------------------
    # Settings profile: default (default.yaml) or a bundled overlay (16S)
    # ------------------------------------------------------------------

    @property
    def active_profile(self) -> str:
        """The recorded ``profile``, else inferred from the annotation.

        Configs without a ``profile`` key (written before 1.0.2) count as
        ``"16S"`` when annotation targets a 16S database (a ``16S.*`` preset,
        or round-1 databases of group 16S), else as ``"default"``.
        """
        if self.profile:
            return self.profile
        if self.annotation.anns:
            return _preset_profile(self.annotation.anns)
        keys = [k.strip() for k in (self.annotation.ann1 or "").split(",") if k.strip()]
        groups = {self.databases[k].group for k in keys if k in self.databases}
        return "16S" if groups == {"16S"} else "default"

    @classmethod
    def _profile_paths(cls) -> list[tuple]:
        """Key paths any bundled overlay profile sets, in file order."""
        paths: list[tuple] = []
        for name in available_profiles()[1:]:
            for path in _leaf_paths(cls._profile_raw(name)):
                if path not in paths:
                    paths.append(path)
        return paths

    def with_profile(self, name: str) -> tuple[EMtraceConfig, list[tuple[str, Any, Any]]]:
        """Return a copy switched to profile *name*, and the settings that changed.

        The settings involved are exactly those the overlay profiles name;
        they take the values of *name* (default.yaml for ``"default"``) and
        every other setting is kept. The annotation preset is kept when it
        already belongs to *name* (e.g. ``16S.emu``), otherwise it becomes
        the profile's preset. The copy records ``profile = name``.

        Changes are ``(dotted key, old, new)``.
        """
        name = canonical_profile(name)
        source = EMtraceConfig.from_profile(name)
        new = self.model_copy(deep=True)
        changes: list[tuple[str, Any, Any]] = []
        for path in self._profile_paths():
            if path == ("annotation", "anns"):
                if self.annotation.anns and _preset_profile(self.annotation.anns) == name:
                    continue
                new.annotation.ann1 = ""
                new.annotation.ann2 = ""
            old, value = _get_path(self, path), _get_path(source, path)
            if old != value:
                _set_path(new, path, value)
                changes.append((".".join(path), old, value))
        new.profile = name
        return new, changes

    def profile_deviations(self) -> list[tuple[str, Any, Any]]:
        """Profile settings that differ from ``active_profile``'s values.

        Covers the keys the overlay profiles name. Returns
        ``(dotted key, value here, profile value)``, empty when the config
        matches its profile.
        """
        source = EMtraceConfig.from_profile(self.active_profile)
        out = []
        for path in self._profile_paths():
            value, expected = _get_path(self, path), _get_path(source, path)
            if value != expected:
                out.append((".".join(path), value, expected))
        return out

    @classmethod
    def _default_raw(cls) -> Dict[str, Any]:
        """The bundled default.yaml, parsed but not yet built."""
        config_dir = importlib.resources.files("emtrace") / "config"
        default_yaml = config_dir / "default.yaml"
        with importlib.resources.as_file(default_yaml) as p:
            with open(p) as fh:
                return yaml.safe_load(fh) or {}

    @classmethod
    def _from_raw(cls, raw: Dict[str, Any]) -> EMtraceConfig:
        """Parse a raw YAML dict into an EMtraceConfig."""
        from emtrace.utils.user_config import load_user_config

        general_raw = raw.get("general.parameters", {})
        taxonomy_raw = raw.get("taxonomy", {})
        host_db_raw = raw.get("host.databases") or raw.get("host_databases", {})
        # db_dir: project config overrides user config; both override shell env var
        # (shell env var precedence is enforced in _apply_db_dir via setdefault).
        user_cfg = load_user_config()
        db_dir: Optional[str] = raw.get("db_dir") or user_cfg.get("db_dir") or None
        # Export it before any sub-model is built: taxonomy paths expand
        # $EMTRACE_DB_DIR when constructed, so leaving this to _apply_db_dir
        # kept a literal "$EMTRACE_DB_DIR" on the first load of a fresh process.
        if db_dir:
            os.environ.setdefault("EMTRACE_DB_DIR", _expand(db_dir))

        # Build preprocess config (bgfilter nested under preprocess; falls back
        # to top-level "bgfilter:" for backwards-compatible YAML files)
        pre_raw = raw.get("preprocess", {})
        chopper_raw = pre_raw.get("chopper", {})
        pabi_raw = pre_raw.get("porechop_abi", {})
        bgf_raw = _drop_bwa(pre_raw.get("bgfilter", raw.get("bgfilter", {})), "bgfilter")
        pre_sub_raw = pre_raw.get("subsample", {})
        preprocess = PreprocessConfig(
            subsample=SubsampleConfig(**pre_sub_raw) if pre_sub_raw else SubsampleConfig(),
            porechop_abi=PorechopABIConfig(**pabi_raw) if pabi_raw else PorechopABIConfig(),
            chopper=ChopperConfig(**chopper_raw) if chopper_raw else ChopperConfig(),
            bgfilter=BgFilterConfig(**bgf_raw) if bgf_raw else BgFilterConfig(),
        )

        # Build annotation config
        ann_raw = raw.get("annotation", {})
        # "anns.options" is the original dotted YAML key; "anns_options" is what
        # model_dump/to_yaml writes back (Pydantic field name); "ann_strategies"
        # is the top-level EMtraceConfig field also written by to_yaml.
        strategies_raw = (
            ann_raw.get("anns.options")
            or ann_raw.get("anns_options")
            or raw.get("ann_strategies")
            or {}
        )
        if ann_raw:
            # Search-level qcov keys were never applied (no search tool filters
            # on coverage); only report.min_qcov filters. Say so, then drop them.
            stale_qcov = [
                f"annotation.{tool}.qcov" for tool in ("blastn", "minimap", "diamond_blastp")
                if isinstance(ann_raw.get(tool), dict) and "qcov" in ann_raw[tool]
            ]
            if stale_qcov:
                logger.warning(
                    "Ignoring %s: search tools do not filter on query coverage. "
                    "Use report.min_qcov (or --min_qcov) to filter the "
                    "annotation/abundance tables.",
                    ", ".join(stale_qcov),
                )
                ann_raw = {
                    k: ({kk: vv for kk, vv in v.items() if kk != "qcov"}
                        if k in ("blastn", "minimap", "diamond_blastp") and isinstance(v, dict)
                        else v)
                    for k, v in ann_raw.items()
                }
            annotation = AnnotationConfig(
                blastn=BlastnConfig(**ann_raw.get("blastn", {})),
                minimap=MinimapConfig(**ann_raw.get("minimap", {})),
                diamond_blastp=DiamondBlastpConfig(**ann_raw.get("diamond_blastp", {})),
                orfipy=OrfipyConfig(**ann_raw.get("orfipy", {})),
                ann1=str(ann_raw.get("ann1", "minimap.nt.abv")),
                ann2=str(ann_raw.get("ann2", "")),
                anns=ann_raw.get("anns"),
                anns_options=strategies_raw,
                filter_tophits=bool(ann_raw.get("filter_tophits", False)),
            )
        else:
            annotation = AnnotationConfig()

        # Build assembly config
        asm_raw = _drop_bwa(raw.get("assembly", {}), "assembly")
        assembly = AssemblyConfig(**asm_raw) if asm_raw else AssemblyConfig()

        # Build taxonomy config
        # YAML layout: taxonomy.ncbi / taxonomy.silva (lowercase keys; uppercase
        # accepted for older files). There is no taxonomy switch: each round
        # reports through its database's taxonomy. Older files carry a
        # "selected" key (boolean per source, or a top-level string); it no
        # longer does anything, so say so rather than ignore it silently.
        if taxonomy_raw:
            ncbi_raw = dict(taxonomy_raw.get("ncbi", taxonomy_raw.get("NCBI", {})))
            silva_raw = dict(taxonomy_raw.get("silva", taxonomy_raw.get("SILVA", {})))
            stale = [
                key for key, present in (
                    ("taxonomy.selected", "selected" in taxonomy_raw),
                    ("taxonomy.ncbi.selected", "selected" in ncbi_raw),
                    ("taxonomy.silva.selected", "selected" in silva_raw),
                ) if present
            ]
            if stale:
                logger.warning(
                    "Ignoring %s: taxonomy is no longer selectable. SILVA "
                    "databases report through SILVA taxonomy, all others "
                    "through NCBI taxonomy. Remove the key from your config.",
                    ", ".join(stale),
                )
            ncbi_raw.pop("selected", None)
            silva_raw.pop("selected", None)
            # Merge YAML values onto class defaults so fields absent from YAML
            # (e.g. "db") retain their hardcoded defaults rather than becoming "".
            _defaults = TaxonomyConfig()
            taxonomy = TaxonomyConfig(
                NCBI=TaxonomySourceConfig(**{**_defaults.NCBI.model_dump(), **ncbi_raw}),
                SILVA=TaxonomySourceConfig(**{**_defaults.SILVA.model_dump(), **silva_raw}),
            )
        else:
            taxonomy = TaxonomyConfig()

        # Build annotation databases dict — YAML is authoritative when present
        databases = _parse_databases_section(raw.get("databases", {}))

        # Build host databases dict
        host_databases = _parse_host_databases_section(host_db_raw)

        # Build runtime config — prefer nested 'runtime:' section, fall back to legacy flat keys
        rt_raw = raw.get("runtime", {})
        runtime = RuntimeConfig(
            num_threads=rt_raw.get("num_threads", raw.get("numth", 8)),
            output_root=rt_raw.get("output_root", raw.get("res", "results")),
            temp_dir=rt_raw.get("temp_dir", raw.get("tmpdir", "$TMPDIR")),
            entrez_email=rt_raw.get("entrez_email", ""),
            entrez_api_key=rt_raw.get("entrez_api_key", ""),
            mode=rt_raw.get("mode", RuntimeMode.ONLINE),
        )

        # Build output config
        out_raw = raw.get("output", {})
        output = OutputConfig(
            fasta_wrap=bool(out_raw.get("fasta_wrap", False)),
            fasta_line_width=int(out_raw.get("fasta_line_width", 80)),
            pack_reads=bool(out_raw.get("pack_reads", False)),
            pack_reads_vi=bool(out_raw.get("pack_reads_vi", True)),
            pack_assembly_bam=bool(out_raw.get("pack_assembly_bam", False)),
            realtime_abund_table=bool(out_raw.get("realtime_abund_table", True)),
        )

        # Build EM config
        em_raw = raw.get("em", raw.get("emabund", {}))
        em = EMConfig(**em_raw) if em_raw else EMConfig()

        # Build Report config
        rep_raw = raw.get("report", {})
        report = ReportConfig(**rep_raw) if rep_raw else ReportConfig()

        # Build Stats config
        sta_raw = raw.get("stats", {})
        stats = StatsConfig(**sta_raw) if sta_raw else StatsConfig()

        return cls(
            preprocess=preprocess,
            annotation=annotation,
            databases=databases,
            assembly=assembly,
            output=output,
            logs=str(raw.get("logs", "logs")),
            res=str(raw.get("res", "results")),
            pipe=str(raw.get("pipe", "main")),
            runtime=runtime,
            taxonomy=taxonomy,
            host_databases=host_databases,
            ann_strategies=strategies_raw or {},
            em=em,
            report=report,
            stats=stats,
            db_dir=db_dir,
            profile=raw.get("profile"),
        )

    # ------------------------------------------------------------------ #
    # Serialization
    # ------------------------------------------------------------------ #

    def to_yaml(self, path: Path | str) -> None:
        """Write configuration to a YAML file."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as fh:
            yaml.dump(
                self.model_dump(mode="json"),
                fh,
                default_flow_style=False,
                sort_keys=False,
            )

    # ------------------------------------------------------------------ #
    # Lookup helpers
    # ------------------------------------------------------------------ #

    def get_host_database(self, key: str) -> DatabaseEntry:
        """Look up a host database by key; raise if missing."""
        from emtrace.core.exceptions import DatabaseNotConfiguredError

        if key not in self.host_databases:
            raise DatabaseNotConfiguredError(key)
        return self.host_databases[key]

    def get_database(self,key: str) -> DatabaseEntry:
        """Look up annotation database by key; raise if missing.

        Falls back to a case-insensitive match so that registry keys carrying
        uppercase (``minimap.rrnDB``, ``minimap.refseq16S``) stay reachable —
        Pipeline._parse_ann1/_parse_ann2 lowercase every ``--ann1``/``--ann2``
        token before looking it up.
        """
        from emtrace.core.exceptions import DatabaseNotConfiguredError

        if key in self.databases:
            return self.databases[key]
        key_lc = key.lower()
        for name, entry in self.databases.items():
            if name.lower() == key_lc:
                return entry
        raise DatabaseNotConfiguredError(key)

    def taxonomy_for_database(self, name: str) -> TaxonomyType:
        """Return the taxonomy that database *name* reports through.

        *name* is a registry key (``minimap.silva``) or a database path, as
        found in an annotation table's ``db`` column. A registry entry's
        explicit ``taxonomy`` wins; otherwise any name or path containing
        ``silva`` is SILVA and everything else is NCBI.
        """
        from emtrace.core.exceptions import DatabaseNotConfiguredError

        name = str(name or "")
        try:
            entry = self.get_database(name)
        except DatabaseNotConfiguredError:
            entry = None
        if entry is not None and entry.taxonomy is not None:
            return entry.taxonomy
        if "silva" in name.lower():
            return TaxonomyType.SILVA
        return TaxonomyType.NCBI

    def resolve_db_path(self, key: str) -> Optional[Path]:
        """Return the absolute path for an annotation or host database entry."""
        entry = self.databases.get(key) or self.host_databases.get(key)
        return _resolve_db_path(entry) if entry else None

    def get_db_url(self, key: str) -> Optional[str]:
        """Return the download URL for an annotation database, or None."""
        entry = self.databases.get(key)
        return entry.url if entry and entry.url else None

    def is_db_installed(self, key: str) -> bool:
        """Return True if the annotation database path exists on disk."""
        path = self.resolve_db_path(key)
        if path is None:
            return False
        if path.exists():
            return True
        return any(path.parent.glob(path.name + "*"))

    def resolve_strategy(self, strategy_name: str) -> str:
        """Return the annotation CLI string for a named strategy."""
        from emtrace.core.exceptions import ConfigurationError

        if strategy_name not in self.ann_strategies:
            available = ", ".join(sorted(self.ann_strategies.keys()))
            raise ConfigurationError(
                f"Unknown annotation strategy: '{strategy_name}'. "
                f"Available: {available}"
            )
        return self.ann_strategies[strategy_name]

    def list_host_databases(self) -> list[str]:
        """Return sorted list of available host genome keys."""
        return sorted(self.host_databases.keys())

    def list_databases(self) -> list[str]:
        """Return sorted list of available annotation database keys."""
        return sorted(self.databases.keys())

    def list_strategies(self) -> list[str]:
        """Return sorted list of available annotation strategies."""
        return sorted(self.ann_strategies.keys())
