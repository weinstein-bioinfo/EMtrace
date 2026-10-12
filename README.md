# EMtrace

A bioinformatics pipeline for taxonomic binning and reference-based consensus
assembly of Oxford Nanopore metagenomic sequencing reads, with command-line,
graphical, and real-time interfaces.

## Documentation

Full documentation is available in [docs/UserGuide.v1.0.md](docs/UserGuide.v1.0.md).

## Features

- **Quality preprocessing** — adapter trimming with *Porechop* / *Porechop_ABI*
  and quality- and length-based filtering with *Chopper*, with per-stage read
  statistics (read counts, length distributions, N50).
- **Host and background depletion** — alignment-based removal of host and
  background reads with *Minimap2* or *BWA-MEM* against one or more reference
  genomes.
- **Taxonomic annotation** — read classification with *Minimap2*,
  *BLASTN* / *BLASTP*, and *DIAMOND* (BLASTX / BLASTP, with *orfipy* ORF
  detection), organised into up to two configurable, chained annotation rounds;
  per-database hit tables are retained.
- **Taxonomic resolution** — NCBI-division and SILVA-domain assignment with
  bacteriophage flagging, partitioning reads into viral, phage, bacterial, and
  unannotated bins.
- **Abundance estimation** — expectation–maximization (EM) estimation of
  per-taxon relative abundances from ambiguously assigned reads.
- **Reference-based assembly** — consensus assembly with references suggested
  from annotation hits (*Minimap2* / *BWA-MEM* alignment, *SAMtools*,
  *BCFtools* consensus).
- **Curated databases** — RefSeq- and GenBank-derived databases supporting
  offline metagenomic, clinical, and surveillance applications, with built-in
  download and management.
- **Real-time operation** — live monitoring and incremental binning of
  *MinION* / *MinKNOW* runs as basecalled FASTQ files are produced, with a final
  report generated on completion.
- **Interactive reporting** — HTML reports with stacked-bar, sunburst, and
  Sankey visualisations, alongside tabular (TSV / CSV) annotation, abundance, and
  quality-control outputs.
- **Multiple interfaces** — a command-line interface, a *PyQt6* desktop GUI, and
  cross-sample batch summaries.

## Citing EMtrace

Article in preparation.

## Installation

From source (Python ≥ 3.10), with the external tools from Bioconda:

```bash
conda create -n emtrace -c conda-forge -c bioconda \
    minimap2 chopper porechop_abi diamond orfipy samtools bcftools
conda activate emtrace
pip install "git+https://github.com/weinstein-bioinfo/EMtrace.git@v1.0.0"
emtrace --help
```

Self-contained macOS and Windows bundles are attached to each
[release](https://github.com/weinstein-bioinfo/EMtrace/releases). They are not
signed, so the operating system warns on first launch; see
[Opening the app](docs/UserGuide.v1.0.md#opening-the-app).

Databases, optional tools (BLAST+, BWA, Medaka) and the GUI extra are described in
the [User Guide](docs/UserGuide.v1.0.md#installing-emtrace).

## License

GNU GPLv3 — see LICENSE and COPYRIGHT.

## Contact

Email: emtrace-support@helsinki.fi
