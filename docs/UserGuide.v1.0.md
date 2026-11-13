# EMtrace User Guide

*Running EMtrace on macOS or Linux.*

**Version:** 1.0 · **Last updated:** 2026-11-05

## Table of Contents

- [About EMtrace](#about-emtrace)
- [Installing EMtrace](#installing-emtrace)
  - [System requirements](#system-requirements)
  - [Setting up directories](#setting-up-directories)
  - [Using a prebuilt bundle (macOS)](#using-a-prebuilt-bundle-macos)
    - [Opening the app](#opening-the-app)
  - [Building from source](#building-from-source)
  - [Installing dependencies](#installing-dependencies)
    - [Tested versions](#tested-versions)
  - [Installing databases](#installing-databases)
  - [Installing host and background filters](#installing-host-and-background-filters)
  - [Creating custom databases and filters](#creating-custom-databases-and-filters)
- [Running EMtrace (GUI)](#running-emtrace-gui)
  - [Creating a job](#creating-a-job)
  - [Configuring settings](#configuring-settings)
- [Annotation strategies and databases](#annotation-strategies-and-databases)
  - [Virus surveillance studies](#virus-surveillance-studies)
  - [Virus metagenomics](#virus-metagenomics)
- [Example 1: A virus surveillance run](#example-1-a-virus-surveillance-run)
- [Running EMtrace (CLI)](#running-emtrace-cli)
- [Real-time monitoring](#real-time-monitoring)
- [Citing EMtrace](#citing-emtrace)
- [Contact](#contact)

---

<a id="about-emtrace"></a>
## About EMtrace

EMtrace is a bioinformatics pipeline for taxonomic binning and reference-based
consensus assembly of long-read metagenomic data, designed primarily for Oxford
Nanopore (MinION/Flongle) viral sequencing. It is implemented in Python and
provides a graphical interface and a real-time mode that processes basecalled
reads as they are produced during a sequencing run.

A typical analysis proceeds through the following stages:

1. **Quality preprocessing** — adapter trimming (*Porechop* / *Porechop_ABI*) and
   quality- and length-based filtering (*Chopper*).
2. **Host and background depletion** — alignment-based removal of host and
   background reads against one or more reference genomes (*Minimap2*).
3. **Taxonomic annotation** — read classification against curated reference
   databases using *Minimap2* and *DIAMOND*, in up to two
   chained rounds.
4. **Abundance estimation** — an expectation–maximization (EM) algorithm resolves
   multi-mapping reads into per-taxon abundance probabilities.
5. **Reference-based assembly** *(optional)* — consensus reconstruction against
   user-supplied references, or references selected automatically from the
   reads' hits (`--ref_auto`).
6. **Reporting** — interactive HTML reports, taxonomic abundance tables, and
   quality-control statistics.

EMtrace is intended for offline, self-contained use on a laptop or
workstation: reference databases are downloaded once and stored locally, so no
network connection is required at analysis time. This makes it suitable for
surveillance and field applications where reproducibility and data locality
matter.

The pipeline has three interfaces that share the same engine and
configuration files:

1. **command-line interface** (`emtrace run`)
2. **desktop GUI** (`emtrace gui`)
3. **real-time monitor** (`emtrace watch`)

---

<a id="installing-emtrace"></a>
## Installing EMtrace

<a id="system-requirements"></a>
### System requirements

EMtrace needs no special hardware. It was developed and benchmarked on a laptop
(Apple M2 Pro, 12 cores, 16 GB RAM), on which most EMtrace analyses of the
EMtrace paper ran, with a peak memory of 3–10 GB. The two largest datasets,
untargeted libraries in a human background, ran on a Linux HPC server and
peaked at 18–19 GB; for such libraries, plan for about 20 GB of RAM. Disk space depends on the reference databases installed:
from 0.24 GB for a 16S database to about 38 GB for
`minimap.genbank.vi+refseq.abv.1g` (sizes in
[Installing databases](#installing-databases)). The
[tutorial](../demo/DEMO.md) needs about 35 GB in all, almost all of it for its
databases. Medaka polishing (optional) can use a GPU but also runs on the CPU.

<a id="setting-up-directories"></a>
### Setting up directories

EMtrace separates *software*, *reference data*, and *results*. We recommend a
layout such as:

```
EMtrace/                 # the application (source or bundle)
EMtrace-databases/       # annotation databases and host/background filters
EMtrace-projects/        # one subdirectory per project (study), one folder per sample
```

The location of the reference data is provided to the pipeline through the
`EMTRACE_DB_DIR` environment variable (or directly in the configuration
file), which database paths in the configuration expand at load time:

```bash
export EMTRACE_DB_DIR=/path/to/EMtrace-databases
```

**NOTE:**
GUI users do not need to set this variable manually. In the desktop
application, open the **Databases** manager and click **Browse** next to the
`EMTRACE_DB_DIR` field to choose the database root. The selection is saved to
a per-user configuration file and applied automatically on all
subsequent runs.  
An `EMTRACE_DB_DIR` value set in the shell environment takes precedence over
the stored setting.

<a id="using-a-prebuilt-bundle-macos"></a>
### Using a prebuilt bundle (macOS)

Self-contained macOS bundles (`.app` in a `.dmg`) for Apple silicon (`arm64`) and
Intel (`x86_64`) Macs are attached to each
[GitHub release](https://github.com/weinstein-bioinfo/EMtrace/releases).
They do not need a separate Python installation.

Windows is not supported in this version; Windows support is planned for a future
release. The [installation from source](#building-from-source) may work on Windows
under WSL2 (Windows Subsystem for Linux), but this is untested.

<a id="opening-the-app"></a>
#### Opening the app

The bundles are not signed with an Apple certificate, so macOS warns the
first time they are opened. This is expected.

| Platform | What you see | What to do |
|:--|:--|:--|
| macOS 15 and later | "Apple could not verify 'EMtrace' is free of malware…" | Close the dialog. Open **System Settings → Privacy & Security**, scroll down to the message about EMtrace and click **Open Anyway**, then confirm. |
| macOS 14 and earlier | The same warning | Control-click (right-click) the app, choose **Open**, then **Open** again. |

<br>

Advanced macOS users can instead remove the download flag in a terminal:

```bash
xattr -dr com.apple.quarantine /Applications/EMtrace.app
```

Each release lists SHA-256 checksums of the bundles. To check a download:

```bash
shasum -a 256 EMtrace-macos-arm64-1.0.5.dmg
```

<a id="building-from-source"></a>
### Building from source

EMtrace is a standard Python package (Python ≥ 3.10). Clone the repository and
install it into a virtual environment:

```bash
git clone https://github.com/weinstein-bioinfo/EMtrace.git
cd EMtrace
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[gui,plots,excel]"   # editable install, with the optional extras
```

The extras are optional: `gui` installs PyQt6 for the desktop application,
`plots` the figure libraries (plotly, kaleido, pysam and others; needed for the
coverage, Sankey and sunburst plots), and `excel` writes `.xlsx` copies of the
tables. `pip install -e .` installs the command line only.

To use the Conda environment of the external tools (below) instead of a
virtual environment, activate it and run the same `pip install` in the cloned
folder.

Verify the installation:

```bash
emtrace --help
```

<a id="installing-dependencies"></a>
### Installing dependencies

EMtrace calls a small set of established command-line tools, which must be
available on the `PATH`.

The Python dependencies (NumPy, pandas, Biopython,
Pydantic, and others) are installed automatically by `pip`.

| Tool | Purpose | Required for | Website | Original article |
|------|---------|--------------|---------|------------------|
| *Chopper* | Quality/length filtering | Preprocessing | https://github.com/wdecoster/chopper | https://doi.org/10.1093/bioinformatics/btad311 |
| *DIAMOND* | Translated protein classification | Protein strategies | https://github.com/bbuchfink/diamond | https://doi.org/10.1038/s41592-021-01101-x |
| *Minimap2* | Alignment, host filtering, classification | All runs | https://github.com/lh3/minimap2 | https://doi.org/10.1093/bioinformatics/bty191 |
| *orfipy* | ORF detection for DIAMOND blastp | Protein strategies | https://github.com/urmi-21/orfipy | https://doi.org/10.1093/bioinformatics/btab090 |
| *Porechop* / *Porechop_ABI* | Adapter trimming | Preprocessing | https://github.com/bonsai-team/Porechop_ABI |  https://doi.org/10.1101/2022.07.07.499093 |
| *SAMtools*, *BCFtools* | Consensus assembly | Assembly step | http://www.htslib.org/ | https://doi.org/10.1093/gigascience/giab008 |
| *[BLAST+]* | Nucleotide/protein classification | BLAST-based strategies [Optional] | https://blast.ncbi.nlm.nih.gov/ | https://doi.org/10.1186/1471-2105-10-421 |
| *[Medaka]* | Consensus polishing (on by default; skipped with a warning if missing) | Assembly step [Optional] | https://github.com/nanoporetech/medaka | — |


#### Installing dependencies with Conda

The external tools are available from Bioconda. We recommend a dedicated
environment, created from `environment.yml` in the repository's top folder,
which pins the [tested versions](#tested-versions):

```bash
conda env create -f environment.yml
conda activate emtrace
```

Without a clone of the repository, download the file of your release first:

```bash
curl -LO https://raw.githubusercontent.com/weinstein-bioinfo/EMtrace/v1.0.5/environment.yml
```

**Installation time.** On the laptop of the
[System requirements](#system-requirements), with empty package caches, creating
the environment took about 30 s and `pip install ".[gui,plots,excel]"` about
1 min. Database downloads take longer and depend on their size: about 15 min for
the 16 GB download of `minimap.refseq.abv.1g` (download sizes and times in the
[tutorial](../demo/DEMO.md#before-you-start)).

Medaka has heavy dependencies and is best kept in its own environment; point
EMtrace to it with `assembly.medaka.path` in the configuration (for example
`~/miniconda3/envs/medaka/bin`):

```bash
conda create -n medaka -c conda-forge -c bioconda medaka
```

#### Installing dependencies manually

Each tool may also be installed from its upstream release and placed on the
`PATH`. EMtrace checks for required tools at startup and reports any that are
missing.

<a id="tested-versions"></a>
#### Tested versions

EMtrace 1.0.5 was tested on 13 November 2026 by running the
[tutorial](../demo/DEMO.md) demos with the external tools of `environment.yml`:
Demos 1–5 with the macOS app, and Demo 1 with the source install in the same
Conda environment. The results matched the tutorial's expected output.
EMtrace 1.0.4 was tested on Linux, on the CSC Roihu supercomputer.

| Component | macOS (EMtrace 1.0.5) | Linux (EMtrace 1.0.4) |
|:--|:--|:--|
| Operating system | macOS 26.6, Apple silicon (native `osx-arm64` builds) | Red Hat Enterprise Linux 9.8, x86_64 |
| Python | 3.12.16 | 3.12.16 |
| *Minimap2* | 2.31 | 2.31 |
| *Chopper* | 0.15.0 | 0.13.0 |
| *Porechop_ABI* | 0.5.1 | 0.5.1 |
| *DIAMOND* | 2.2.8 | 2.2.4 |
| *orfipy* | 0.0.4 | 0.0.4 |
| *SAMtools*, *BCFtools* | 1.25 | 1.23.1 |
| Python packages | NumPy 2.5.4, pandas 3.1.1, Biopython 1.88, Pydantic 2.14.3, PyQt6 6.12.0, matplotlib 3.11.2, plotly 7.1.0, kaleido 1.5.1, pysam 0.24.1 | NumPy 1.26.4, pandas 3.0.6, Biopython 1.85, Pydantic 2.14.1, matplotlib 3.11.2, plotly 5.24.1, pysam 0.24.1 |

<br>

`environment.yml` pins the macOS versions and resolves to the same tool
versions on Intel macOS (`osx-64`) and Linux (`linux-64`). The optional tools
(*Medaka*, *BLAST+*) are not used by the tutorial.

<a id="installing-databases"></a>
### Installing databases

Annotation databases are managed through the `databases` subcommand:

```bash
emtrace databases --list                 # show available / installed databases
emtrace databases --install minimap.refseq.abv.1g
```

Databases are downloaded from public object storage and stored under
`EMTRACE_DB_DIR`.

Database names encode their taxonomic scope with the
suffixes used throughout EMtrace:

- **vi** — viruses only;
- **ab** — archaea and bacteria;
- **abv** — archaea, bacteria, and viruses (RefSeq representative/reference
  genomes);
- **abv.1g** — as *abv*, but reduced to one genome per genus for bacteria and
  archaea (largest genome), while all available viral genomes are retained.

Supported reference databases include:

| Name | Search | Size (GB) | Description | Applications |
|------|--------|-----------|-------------|--------------|
| `minimap.refseq.abv.1g` | Minimap2 | ~23 | RefSeq archaeal, bacterial, and viral genomes, one genome per bacterial/archaeal genus (all viruses retained); distributed as a prebuilt index. | Surveillance / clinical identification of known viruses and bacterial genera. |
| `minimap.refseq.abv` | Minimap2 | ~28 (FASTA) | Full RefSeq archaeal, bacterial, and viral genomes (all representative species). Distributed as an un-indexed FASTA and indexed at runtime. | Adds bacterial species-level resolution, at the cost of a large Minimap2 index (≈70 GB, estimated, built at runtime). |
| `minimap.refseq.vi` | Minimap2 | ~2.2 | RefSeq viral genomes; distributed as a prebuilt index. | Fast, virus-only first-pass classification; the default assembly database for automatic reference selection. |
| `minimap.genbank.vi+refseq.abv.1g` | Minimap2 | ~38 | GenBank viral sequences combined with the `refseq.abv.1g` set; distributed as a prebuilt index. | Viral metagenomics. |
| `blastn.refseq.abv` | BLASTN | ~7 | BLASTN index of the `refseq.abv` set (see above). | As `minimap.refseq.abv`, for slower but more sensitive BLASTN search. |
| `blastn.refseq.abv.1g` | BLASTN | ~2.3 | BLASTN index of the `refseq.abv.1g` set (see above). | As `minimap.refseq.abv.1g`, for slower but more sensitive BLASTN search. |
| `blastn.refseq.vi` | BLASTN | ~0.14 | BLASTN index of RefSeq viral genomes. | Second-round confirmation of viral reads (`vi.2round`-type strategies). |
| `diamond.nrclust.vi` | DIAMOND | ~3.1 | Clustered NCBI nr viral proteins. | Protein-level detection of divergent viruses (`vi.chain1`). |
| `diamond.uniref.vi` | DIAMOND | ~1.6 | UniRef100 viral proteins. | Protein-level detection of divergent viruses (`vi.chain2`; used in the EMtrace manuscript). |
| **16S rRNA databases** | | | | |
| `minimap.silva` | Minimap2 | ~6.3 | SILVA SSU rRNA reference sequences (SSURef). | SILVA SSU rRNA ribosomal profiling. |
| `minimap.rrnDB.dedup` | Minimap2 | ~0.4 | rrnDB 5.10 16S rRNA genes, de-duplicated (one representative per named species, else genus or family, per cluster of identical sequences; 78,966 sequences). | 16S profiling with NCBI taxonomy (`16S.rrnDB.dedup`; the 16S profile default). |
| `minimap.emu` | Minimap2 | ~0.24 | Emu's own 16S database (rrnDB 5.6 + NCBI 16S RefSeq, one sequence per species; 49,243 sequences), unchanged. | 16S profiling directly comparable to Emu (`16S.emu`). |


<a id="installing-host-and-background-filters"></a>
### Installing host and background filters

Host and background genomes are used to deplete non-target reads
before annotation. They are provided either as configuration keys (e.g.
`Homo_sapiens`) or as direct paths to FASTA files. Multiple filters can be
chained — for example, the human genome together with one or more arthropod
vector genomes (e.g. `Homo_sapiens,Ixodes_scapularis` )

<a id="creating-custom-databases-and-filters"></a>
### Creating custom databases and filters

Any FASTA file can be used directly as a host filter or as a Minimap2 annotation
reference; EMtrace builds and caches the required index on first use. Guidance
on assembling custom curated databases will be added in a later revision.

---

<a id="running-emtrace-gui"></a>
## Running EMtrace GUI

Start the desktop application in one of two ways:

- **Installed app (prebuilt bundle):** double-click **EMtrace** in
  *Applications* (or find it with Spotlight, ⌘-Space). The first start shows a
  security warning;
  see [Opening the app](#opening-the-app).
- **From a terminal (installed from source):** activate the environment
  EMtrace was installed into and run:

```bash
emtrace gui
```

The window is organised into three panels: a **step list** (left), a **content
area** with dashboard, configuration, plots, and tables (centre), and a **live
log** (bottom).

<a id="creating-a-job"></a>
### Creating a job

A **job** analyses one sample: one FASTQ file, typically one barcode of a
sequencing run, or a directory of FASTQ files for real-time monitoring. From the
top menu select **New Job**, select the input, choose a results directory and give
the sample name. EMtrace creates `<results directory>/<sample name>/` with
subdirectories for reads, annotation, assembly and reports.

A **project** is simply a results directory shared by several jobs. Depending on
the study design, it can hold many samples from one sequencing run, or samples
from several runs. For example:

```
EMtrace-projects/
└── surveillance-2026/          # project
    ├── run01_barcode01/        # job (sample)
    ├── run01_barcode02/
    └── run02_barcode01/
```

Use sample names that are unique within the project, e.g. by including the run
name, because each job writes to a folder named after its sample. The
cross-sample summaries (`emtrace batch-summary`, `batch-taxa`,
`batch-abundance`; see [Example 1](#example-1-a-virus-surveillance-run)) take the
project directory and report all samples in it together, across runs.

**NOTE**: *Reads R2 (paired-end)* are only provided for reference-based assembly of Illumina reads. For Nanopore jobs this can be left empty.

<a id="configuring-settings"></a>
### Configuring settings

Pipeline parameters are edited in the configuration panel and saved with the
job. Every parameter has a sensible default; the groups below describe those
most often adjusted. The complete set is documented in the bundled
`default.yaml`.

**Profile.** The *Profile* drop-down at the top of the panel selects a set of
defaults: *default (virome)* for viral metagenomics, or *16S* for full-length
16S rRNA amplicons (the benchmarked 16S settings: `16S.rrnDB.dedup`, no host
filtering, 1,200–1,800 bp reads, 1 % tail trim). It is the GUI equivalent of
the CLI's `--profile`. Switching changes only the settings the profile defines
and keeps all others; picking a `16S.*` annotation preset switches to *16S*,
any other preset back to *default*. The profile is saved in `config.yaml`
(`profile:`) and recorded in each job's `job.yaml`. When profile settings
have been changed by hand, a *modified* note appears next to the drop-down;
hover over it to see which settings differ from the profile. The run log lists
the same differences.

**Runtime** — the number of threads and the run *Mode*. In *online* mode the
pipeline may fall back to NCBI Entrez to fetch a reference that is not present in
the local database (a valid contact e-mail is required); in *offline* mode it
never accesses the network and aborts if a reference is missing.

**Preprocessing.** Three sub-steps run in sequence — adapter trimming → quality
filtering → host/background depletion — and each is enabled by including its step
in the pipeline.

- *Subsample reads* — when enabled, a random subsample of the input (a fraction
  between 0 and 1, or an absolute read count) is drawn before any other step.
  This is useful for rapid pre-flight tests on a fraction of a run.
- *Porechop / Porechop_ABI* — adapter trimming. *Porechop_ABI* can infer
  adapters *ab initio* from the reads (`ab_initio`, disabled by default); a custom adapter file may be
  supplied instead. Reads containing internal adapters are split (disabled by default). The principal thresholds are the end- and middle-adapter
  percent-identity cutoffs (`end_threshold`, `middle_threshold`).
- *Chopper* — quality and length filtering. The main parameters are the minimum
  mean Phred quality (`quality`, default 10), the minimum read length
  (`minlength`), an optional maximum length, and head/tail
  cropping. The end-trimming strategy (`trim_approach`) determines how
  low-quality read ends are handled.
- *Host and background filtering* — depletion of non-target reads by alignment
  against one or more reference genomes, given as configuration keys or FASTA
  paths (comma-separated and applied in turn — e.g. a host genome together with a
  vector genome). Reads are aligned with *Minimap2*, and a read is removed once
  its alignment score reaches the configured minimum (`m` / `s`). A reusable
  `.mmi` index is built and cached alongside each genome.

**NOTE:**
The default minimum read length (`chopper.minlength: 400`) suits libraries of
mostly long reads. For libraries with shorter reads, e.g. degraded RNA, short
amplicons or field samples, it can discard a large share of the data. Check the
read lengths before a run with [SeqKit](https://bioinf.shenwei.me/seqkit/)
(`conda install -c bioconda seqkit`):

```bash
seqkit stats -a -T reads.fastq.gz
```

The `Q1`, `Q2` (median) and `N50` columns give the length distribution. If the
median read length is near or below 400 bp, lower `minlength` (e.g. to 150, or
100 for very short reads).

**Annotation.** Reads are classified against the databases configured for round 1
(`ann1`) and, optionally, round 2 (`ann2`); a named *strategy* (`anns`, e.g.
`abv.1g`) can set both at once (see [Annotation strategies and
databases](#annotation-strategies-and-databases)). Round 1 classifies all reads,
while round 2 re-searches a chosen partition (viral, phage, or bacterial) against
a more sensitive database. Per-tool settings control sensitivity and hit
filtering:

- *Minimap2* — `m` (chain-score threshold) and `s` (DP alignment-score threshold,
  effective only in full-DP mode, `c: true`); the minimizer parameters `k` / `w`
  (smaller values such as `k = 11, w = 5` increase sensitivity for divergent
  viruses); and `max_target_seqs`.
- *BLASTN* — `bitscore`, `evalue`, query coverage (`qcov`), and `max_target_seqs`.
- *DIAMOND* (translated protein search) — `bitscore`, `evalue`, and the
  sensitivity `mode` (e.g. `very-sensitive`). ORFs are first detected with
  *orfipy* (minimum ORF length `min_len`, genetic-code `table`, and strand).

Score and coverage thresholds can also be overridden globally from the command
line (`--min_minimap_score`, `--min_blastn_bits`, `--min_diamond_bits`,
`--min_qcov`).

**Assembly** *(optional)*. Reference-based consensus reconstruction. References
are supplied as FASTA paths (`refseqs`) or NCBI accessions (`refaccs`, fetched
through Entrez in online mode), or selected automatically (`ref_auto`, below).
Reads are mapped with *Minimap2*; a consensus is then called
and positions below the minimum read depth (`min_coverage`, default 10) are masked
as `N`. In the GUI, automatic selection is the *Auto-select references* checkbox;
a reference FASTA file or accessions, when given, always take precedence.

**Automatic reference selection** (`assembly.ref_auto`, `--ref_auto`) applies
when neither `refseqs` nor `refaccs` is set. It uses two sources of hits:

- **primary**: the reads are mapped to the assembly database (`assembly.database`,
  RefSeq viral genomes by default). The hits are saved to
  `annotation/assembly_candidates.tsv`, which the GUI's *Pick from viral hits*
  list also reads;
- **secondary**: the annotation hits (`annot2.tsv`, else `annot1.tsv`).

Only hits from nucleotide searches (Minimap2, Blastn) are candidates. Protein
searches (DIAMOND, Blastp) hit protein accessions, which cannot serve as
nucleotide references, so their hits are ignored, including those in a chained
round such as `vi.chain2`.

In each source, reference accessions are chosen by greedy read cover: the
accession that aligns best to the reads not yet covered is taken, as long as it
covers at least `min_reads` of them. The segments of one genome recruit
different reads, so each segment gets a reference. Strains of one segment share
reads, so only the best-aligning strain is taken. Up to `max_species` viruses
are kept for co-infections.

A secondary pick replaces a primary pick only if it aligns better, on the reads
both share, by more than `secondary_margin`; the alignment quality is set by
`quality_metric`, by default matched bases per read length. A virus that the
primary search misses is added from the secondary source. With
`--no-ref_auto_assembly_db` or `--no-ref_auto_annotation`, only one source is
used.

Every choice, and every rejected candidate with the reason, is written to
`assembly/ref_selection.tsv`.

If no sequence can be obtained for the chosen references (for example, an
accession missing from the local databases with no Entrez access), the
assembly step stops with an error instead of producing an all-`N` consensus.

Variants are then reported at two independent tiers, because the variant table
and the consensus sequence have different requirements. The **variant table**
(`vcf.gz`) lists every allele reaching `report_min_af` (default 0.2), so
sub-dominant alleles from co-infections or intra-host variation remain visible.
The **consensus FASTA** can hold only one base per position, so only an allele
strictly exceeding `min_af` (default 0.5 — the dominant allele) is written into
it. Both tiers additionally require `min_ac` supporting reads (default 5), and
`report_min_af` must not exceed `min_af`. An optional `vcf_flt` expression is
ANDed into both filters for advanced use; it is empty by default.

**Medaka polishing** (`assembly.medaka`, Nanopore single-end reads only) refines
the consensus with *Medaka* after variant calling. It is on by default
(`enabled`, CLI `--medaka` / `--no-medaka`) and is skipped with a warning if
Medaka is not installed. The model is read from the read headers (`model:
auto`, CLI `--medaka_model`); `fallback_model` is used when the headers carry
none (empty = Medaka's default). `path` points to the Medaka executable or its
`bin/` directory when Medaka lives in its own environment (empty = on `PATH`),
`batch_size` sets the inference batch size (lower it if GPU memory runs out),
and `cpu` forces inference on the CPU. In the GUI these settings are under
*Assembly › Medaka polishing (ONT)*.

**EM abundance.** Many reads align ambiguously to several taxa. An
expectation–maximization algorithm iteratively redistributes them to estimate
per-taxon relative abundances. The main controls are the maximum number of
iterations (`max_iter`), the convergence threshold on the change in
log-likelihood (`log_L_threshold`), and a temperature parameter (`temp`) that
smooths the conversion of alignment scores to probabilities. Reads with no
species-level assignment can optionally be excluded (`exclude_no_species`).

**Reporting.** Controls the report contents and the output tables.

- Annotation and abundance tables can be filtered by minimum query coverage
  (`min_qcov`) and minimum bitscore (`min_bitscore`), and exported to Excel
  (`export_excel`).
- Interactive figures are toggled individually: a Sankey diagram (`sankey`), a
  sunburst chart (`sunburst`), and a stacked-bar chart (`stacked_bar`); genus- and
  species-level heatmaps are available as PDFs (`heatmap_genus`,
  `heatmap_species`).
- Reads may optionally be sorted into a taxonomic directory hierarchy
  (`sort_reads`); the field used for this is often the "domain" but can be set to a proper taxonomic rank or `auto detect`
- A companion statistics step records per-stage read counts and length
  distributions, writing them as CSV (or YAML) with optional read-length
  histograms.

Click **Save Config** from the top menu to write your settings to a `config.yaml`
in the job's directory; when the job is opened again, its `config.yaml` takes
priority over the bundled defaults.

---

<a id="annotation-strategies-and-databases"></a>
## Annotation strategies and databases

The choice of reference database determines both sensitivity and runtime, and
should be matched to the study design and sample type. EMtrace supports two
broad use cases.

<a id="virus-surveillance-studies"></a>
### Virus surveillance studies

Targeted detection and characterization of known viruses, e.g. in vector
surveillance or outbreak response.

**Tissue and vector samples.** These samples have a well-defined host or carrier,
so background filtering is effective and a compact viral reference suffices:

- *Background filtering:* the host/background genome(s) — for example
  `Homo_sapiens` together with the relevant vector genome.
- *Annotation database:* `minimap.refseq.abv.1g` — a compact RefSeq-derived viral
  index giving fast, specific first-pass classification.

**Feces / environmental samples.**

> *[Planned]* These sample types carry a complex, poorly defined background and
> generally require large general-purpose databases (e.g. NCBI `nt`) for reliable
> classification. A curated database and recommended workflow for fecal and
> environmental samples are scheduled; this section is a placeholder.
>
> *Disclaimer:* with the compact viral databases alone, false negatives are
> likely due to host/environment eukaryotic mRNA in these samples.

<a id="virus-metagenomics"></a>
### Virus metagenomics

Discovery-oriented profiling, including detection of divergent or previously
uncharacterized viruses.

**Tissue and vector samples.**

- *Background filtering:* the host/carrier genome(s), as above.
- *Annotation strategy:* `--anns vi.chain3` (recommended). It chains a
  nucleotide search to a protein search: reads are first searched against
  `minimap.genbank.vi+refseq.abv.1g`, a broad index that adds GenBank viral
  sequences to the RefSeq set, and the reads it leaves unclassified are translated
  into ORFs (*orfipy*) and searched with *DIAMOND* against viral proteins
  (`diamond.uniref.vi`), which tolerates greater genetic distance.

Predefined chained strategies (all run in annotation round 1; the protein search
receives only the reads the nucleotide search left unclassified):

| Strategy (`--anns`) | Nucleotide search | Protein search | Use |
|:--|:--|:--|:--|
| `vi.chain3` | `minimap.genbank.vi+refseq.abv.1g` (~38 GB) | `diamond.uniref.vi` | Discovery: broadest nucleotide coverage of known viral diversity |
| `vi.chain2` | `minimap.refseq.abv.1g` (~23 GB) | `diamond.uniref.vi` | Divergent viruses with a smaller index; benchmarked in the EMtrace manuscript |
| `vi.chain1` | `minimap.refseq.abv` (full RefSeq, indexed at runtime) | `diamond.nrclust.vi` | Protein search against clustered NCBI nr viral proteins |

<br>

The protein search defaults, DIAMOND bitscore ≥ 80 and ORFs ≥ 90 nt
(`annotation.diamond_blastp.bitscore`, `annotation.orfipy.min_len`), are the
settings benchmarked with `vi.chain2` in the EMtrace manuscript. Nanopore read
errors break reads into short ORFs, whose hits score low, so stricter values
discard most genuine protein hits. The low threshold also admits more spurious
hits: treat a species supported only by protein hits as tentative, and confirm it
with a nucleotide hit, an assembly or a second method.

**Feces / environmental samples.**

> *[Planned]* As above, robust analysis of fecal and environmental samples
> requires large general-purpose databases. This section is a placeholder.
>
> *Disclaimer:* proper analysis most likely requires large databases such as
> NCBI `nt`.

---

<a id="example-1-a-virus-surveillance-run"></a>
## Example 1: A virus surveillance run

This example illustrates a complete surveillance analysis of a single barcode.

**Sample data.** A basecalled, demultiplexed FASTQ file from a MinION run, e.g.
`barcode01.fastq.gz`.

**Pipeline settings.**

- Steps: `main` (preprocessing → host filtering → annotation → EM abundance →
  reporting → packaging).
- Background filter: `Homo_sapiens`.
- Annotation database (round 1): `minimap.refseq.abv.1g`.

**Running (CLI).**

```bash
emtrace run \
  --reads   barcode01.fastq.gz \
  --sample  barcode01 \
  --res     ./EMtrace-projects/run01 \
  --pipe    main \
  --bgfilter Homo_sapiens \
  --ann1    minimap.refseq.abv.1g \
  --numth   8 \
  -v
```

**Result description.**

The run produces, under
`EMtrace-projects/run01/barcode01/`:

- `reads/` — preprocessed and partitioned read sets (viral, phage, bacterial,
  and unannotated bins).
- `annotation/` — raw per-read annotation tables (detailed below).
- `results/reports/` and `results/plots/` — tables, figures, and the HTML report
  (detailed below).

#### Annotation tables (`annotation/`)

The annotation step writes tab-delimited tables of per-read hits. `annot1.tsv`
holds all hits of the annotation; most analyses use this single annotation round.
Each search also writes its own table, `annot1.<database>.tsv` (for example
`annot1.minimap.refseq.abv.1g.tsv`). In a chained strategy such as `vi.chain2`
(Minimap2, then DIAMOND), this keeps the results of the first search even if a
later search fails. An optional second annotation round (`--ann2`) writes the
same files with the prefix `annot2`. All tables share the same column schema:

| Column | Description |
|--------|-------------|
| `search` | Search program (Minimap2, BLASTN, DIAMOND blastp/blastx). |
| `db` | Reference database (configuration key). |
| `qseqid` | Query (read) identifier. |
| `orf` | ORF identifier for protein searches; `NA` for nucleotide searches. |
| `qseqlen` | Full read length (bp). |
| `sseqid` | Subject (reference) sequence identifier. |
| `bitscore` | Alignment bitscore. |
| `evalue` | Expectation value. |
| `alen` | Alignment length (bp). |
| `pident` | Percent identity. |
| `qlen` | Aligned query-region length (bp). |
| `qcov` | Query coverage (fraction, 0–1). |
| `slen` | Subject sequence length (bp). |
| `scov` | Subject coverage (fraction, 0–1). |
| `staxid` | Subject taxonomy ID (NCBI or SILVA). |
| `sname` | Subject name / title. |
| `division` | NCBI division (e.g. *Viruses*, *Bacteria*) or, for SILVA, the domain. |
| `bphage` | Bacteriophage flag (`yes`/`no`); assigned for NCBI taxonomy only. |

The same schema is used by `assembly_candidates.tsv`, a viral-only hit list
written for the assembly step.

#### Reports and figures (`results/reports/`, `results/plots/`)

| File | Description |
|------|-------------|
| `reports/abund_table.tsv` | Per-taxon abundance table: read counts, EM probabilities, and taxonomic lineage. |
| `reports/annot_table.tsv` | Deduplicated, taxonomy-linked annotation table (nucleotide and protein hits merged) used for reporting. |
| `reports/stats_fastq.csv` | Per-stage read statistics: read count, min/max/median/mean length, total bp, and N50. |
| `reports/stats_div.csv` | Read counts per taxonomic division. |
| `reports/stats_report.html` | Interactive HTML summary report (embeds the figures below). |
| `EM.taxid_prob.tsv` | Final EM abundance probability per taxon (written in `results/`). |
| `EM.taxid_prob.byiter.tsv` | EM probability per taxon at each iteration (convergence trace). |
| `plots/stacked_bar.png`, `.html` | Taxonomic composition as a stacked bar chart. |
| `plots/sunburst.html` | Hierarchical taxonomic composition (interactive sunburst). |
| `plots/sankey.html` | Taxonomic read flow (interactive Sankey diagram). |
| `plots/qc/read_lengths.png` | Read-length distribution histogram. |

Annotation and abundance tables are additionally exported as `.xlsx` when
`export_excel` is enabled.

Genus- and species-level heatmaps
(`plots/heatmap_*.pdf`) are produced when enabled in the reporting settings.

**Cross-sample summaries.**

Summaries can be computed over one or several samples or runs:

```bash
emtrace batch-summary ./EMtrace-projects/run01 -o run01.summary.tsv
emtrace batch-taxa    ./EMtrace-projects/run01 --target vi -o run01.viral.tsv
```

<br>
`batch-summary` writes one row per sample, tracking read flow through the
pipeline:

| Column | Description |
|--------|-------------|
| `Sample` | Sample name. |
| `reads` | Input read count. |
| `reads.porechop` | Reads passing adapter trimming. |
| `reads.chopped` | Reads passing quality/length filtering. |
| `reads.hostflt` | Reads remaining after host/background depletion. |
| `n50` | Read-length N50 of the input. |
| `reads.{division}` | Reads assigned to each taxonomic division (one column per division observed, e.g. `reads.Viruses`, `reads.Bacteria`). |
| `reads.mapped` | Total annotated (classified) reads. |
| `reads.unmapped` | Reads with no annotation. |

<br>
`batch-taxa` writes one row per taxon per sample, filtered by `--target`
(`vi` viruses, `ba` bacteria/archaea, `ph` phages):

| Column | Description |
|--------|-------------|
| `Sample` | Sample name. |
| `taxid` | NCBI/SILVA taxonomy ID. |
| `readn` | Reads assigned to the taxon: EM-adjusted counts (`readn_em`, rounded to whole reads), as benchmarked; raw read counts where `readn_em` is not available. |
| `prob` | EM-estimated relative abundance. |
| `organism` | Scientific name. |
| `species` | Species-rank lineage name. |
| `genus` | Genus. |
| `family` | Family. |
| `order` | Order. |
| `class` | Class. |
| `division` | NCBI division (used by the `--target` filter). |

---

<a id="running-emtrace-cli"></a>
## Running EMtrace (CLI)

The full pipeline is run with `emtrace run`. The most frequently used options
are:

| Option | Description |
|--------|-------------|
| `-1, --reads` | Input FASTQ file (plain or gzipped). |
| `-s, --sample` | Sample name / label. |
| `-r, --res` | Root results directory. |
| `-p, --pipe` | Workflow name (`main`, `annotate`, `report`, `realtime`) or comma-separated steps (`pore,chop,flt,ann,em,ass,rep,plot,pack,clean`). |
| `--bgfilter` | Host/background genome key(s) or FASTA path(s). |
| `--ann1` / `--ann2` | Annotation database(s) for round 1 / `target:db` pairs for round 2. |
| `--anns` | Named annotation strategy (overrides `--ann1`/`--ann2`). |
| `--refaccs` / `--refseqs` | References for reference-based assembly. |
| `--ref_auto` | Select assembly references automatically (see *Assembly*). |
| `--medaka` / `--no-medaka` | Polish the consensus with *Medaka* (on by default); `--medaka_model` sets the model (`auto` = from the read headers). |
| `-t, --numth` | Number of threads. |
| `--config` | Configuration YAML; replaces the bundled defaults. |
| `--profile` | Settings profile: `default` (alias `virome`) or `16S`, overlaid on the bundled defaults (see *Configuring settings*). Cannot be combined with `--config`. |
| `-v` | Increase verbosity (`-v` info, `-vv` debug). |

<br>
Selected steps may be run in isolation — for example, to regenerate reports from
an existing annotation without repeating the search:

```bash
emtrace run --reads barcode01.fastq.gz --sample barcode01 \
  --res ./EMtrace-projects/run01 --pipe rep,pack
```

**Pre-building host-filter indices.** When background filtering with *Minimap2*
(`bgfilter.minimap.save_mmi_index: true`), the host genome is indexed to a `.mmi`
file on first use and reused thereafter. For a large collection run as many
parallel jobs (e.g. an HPC array job), build the indices once, before launching
the jobs:

```bash
emtrace databases --build-index hostdbs --config emtrace.config.yaml
```

This downloads each configured host genome (if not already present), concatenates
it, and builds the `.mmi` index at the exact path the runs expect. Pass a single
host key (e.g. `--build-index Homo_sapiens`) to build just one, and `--force` to
rebuild an existing index.

Without this step the first jobs build the index themselves. This is safe: an
index is written under a temporary name and only appears once complete, and
concurrent jobs wait for one build and then reuse it. But every waiting job holds
its slot (and memory) meanwhile, and on file systems without file locking (some
network file systems) each job builds its own copy.

Other subcommands:

- `emtrace databases --list` / `--install KEY` — list configured databases and
  their install status, and download and install a database.
- `emtrace databases --build-index KEY` — pre-build the *Minimap2*
  `.mmi` index for one host filter (`KEY`) or for all configured host filters
  (`--build-index hostdbs`); see *Pre-building host-filter indices* above.
- `emtrace batch-summary` / `batch-taxa` — cross-sample QC and taxonomic
  summaries.
- `emtrace batch-abundance RESDIR -o OUTDIR` — cross-sample OTU/feature
  tables (rows = species, columns = samples), one file per division group
  (bacteria, bacteriophages, viruses, invertebrates_parasites,
  invertebrates_other). Each table has `species` + `taxid` as the first two
  columns and a `taxpath` lineage as the last. Options: `--value
  readn|prob|both` (read counts and/or the EM assignment probability; default
  `both`). Read counts are the EM-adjusted counts (`readn_em`, rounded to whole
  reads), as benchmarked; raw `readn` is used where `readn_em` is not available.
  `--xlsx` to also write a single styled Excel
  workbook of the read-count tables (one sheet per division, header rotated,
  compact sample columns); `--config FILE` to override the invertebrate
  parasite genus allowlist via `report.parasite_genera`.
- `emtrace gui` — launch the desktop interface.

---

<a id="real-time-monitoring"></a>
## Real-time monitoring

> *[Experimental]* Real-time mode is implemented but **not yet validated on live
> sequencing runs**; treat results as preliminary and verify with a standard
> (post-run) analysis.

`emtrace watch` monitors a directory for new FASTQ files (e.g. MinKNOW
output) and processes each batch through preprocessing, host filtering,
annotation, and EM abundance estimation, updating a running abundance table.
Stopping (Ctrl+C) merges the results into the job and writes the report; see
[Stopping and next steps](#watch-stopping).

```bash
emtrace watch \
  --dir   /path/to/minknow/output \
  --sample run01 \
  --hostgen Homo_sapiens \
  --ann1   minimap.refseq.abv.1g \
  --report \
  -v
```

Settings come from `--config` or `--profile` (e.g. `--profile 16S` for 16S
amplicons), as for `emtrace run`. After each batch the console prints the top
taxa with their names.

Files may be written, copied, or moved/renamed into the watched directory
(e.g. a Finder drag, `mv`, or `rsync`, which writes a temporary file and
renames it). A file counts once it has stopped growing for `--debounce`
seconds. The directory is also rescanned every 30 s, which catches files whose
change notifications were missed (e.g. on network drives). A file that was
already processed is not processed again after a rename. FASTQ files already in
the directory when monitoring starts are ignored unless `--scan-existing` is
given, except when an earlier watch of the job is continued (see below).

**In the GUI**, the **Realtime** tab runs the same monitor. Set the
*Watch directory*, the *Batch size* (files per batch) and *Debounce* (seconds),
choose whether to *Process files already in the directory* (on by default), and
click **▶ Start Watching**. The tab shows the current batch, total reads, files
and taxa detected, and a *Top taxa* list by EM abundance, with species names
(hover for the taxid). **■ Stop** ends monitoring after the current batch.

Watching needs a job of its own, because its results are merged into the job.
If the open job is a watch job, the tab continues it. Otherwise, **▶ Start
Watching** asks for a name and output directory and creates one, without asking
for input reads. A job with its own input reads or results is never overwritten.
The same applies to `emtrace watch`, which uses the job directory
`<res>/<sample>` and stops with an error if that job has its own input or
results. The job directory and the watched directory may not be inside one
another, because the watcher also looks in subdirectories.

<a id="watch-stopping"></a>
### Stopping and next steps

While monitoring, batches and the live tables are in `<job>/watch`. Stopping
then does the following:

1. Finishes the current batch and any files waiting for a full batch. Other
   files still in the queue are left for the next watch of the job.
2. Merges all completed batches into the job's standard folders: `reads/`,
   including the input as `reads/reads.input.fq.gz`, and `annotation/`.
3. Runs EM and the report on the merged data, writing
   `results/reports/stats_report.html` and the other standard outputs. With no
   hits, the report holds QC statistics only.

EM runs from scratch here, as in `emtrace run`, so the final numbers can differ
slightly from the live ones.

The job can then continue like any other. In the GUI, select **assembly**,
**report**, **plot**, **pack** or **clean** in the step list and click **Run**.
Preprocessing, host filtering and annotation are disabled for a watch job,
because the next merge would replace their results. On the command line, point
`emtrace run` at the same job:

```bash
emtrace run \
  -1 <res>/run01/reads/reads.input.fq.gz \
  -s run01 -r <res> \
  -p ass,rep,plot,pack,clean
```

To add more data, start watching the same job again. Earlier batches are
reloaded, files already processed are skipped, and files that arrived in the
meantime are picked up, whatever the *Process files already in the directory*
setting or `--scan-existing` say. Stopping then updates the results from all
batches. If EMtrace was closed without stopping, for example after a crash, do
the same: start watching the job and stop it to merge its batches.

A batch that fails is reported and left out of the results, and its files are
processed again the next time the job is watched.

**clean** removes the batches' intermediate files, but keeps what continuing the
watch needs. Merged reads exist only in `reads/`, so clearing them from the
Dashboard cannot be undone.

---

<a id="citing-emtrace"></a>
## Citing EMtrace

Article in preparation.

---

<a id="contact"></a>
## Contact

Email: emtrace-support@helsinki.fi
