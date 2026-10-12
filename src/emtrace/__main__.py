"""Command-line interface for EMtrace.

Usage examples::

    emtrace run -1 reads.fq.gz -s sample01 -p main \\
        --hostgen Homo_sapiens -t 16
    emtrace run -1 reads.fq.gz -s sample01 \\
        --anns vi.nt --hostgen Homo_sapiens
    emtrace databases --list
    emtrace databases --install minimap.nt.abv
    emtrace gui
"""

from __future__ import annotations

import logging
import os
import sys
from pathlib import Path
from typing import Optional

import click

from emtrace import __version__
from emtrace.core.config import EMtraceConfig
from emtrace.core.exceptions import ConfigurationError
from emtrace.core.pipeline import Pipeline
from emtrace.core.pipelinejob import PipelineJob
from emtrace.utils.db_download import (
    DatabaseAlreadyInstalledError,
    DatabaseInstallError,
    DatabaseNotConfiguredError,
    database_status,
    install_by_key,
    list_installed_files,
    prebuild_host_filters,
    remove_database,
)
from emtrace.utils.logging import (
    quiet_external_loggers,
    setup_logging,
)


# ------------------------------------------------------------------
# CLI group
# ------------------------------------------------------------------

@click.group()
@click.version_option(version=__version__, prog_name="emtrace")
def cli():
    """EMtrace: Nanopore read-binning and taxonomic classification."""
    # Strip any inherited LD_PRELOAD before we spawn external tools.
    # On CSC Puhti the `python-data` tykky container injects
    # LD_PRELOAD=/appl/soft/ai/lib/fake_tcp_ulp.so (a TCP shim built against a
    # newer glibc). None of emtrace's external tools (minimap2, samtools,
    # diamond, tar, plotting backends, …) need it, and passing it through to a
    # host /usr/bin/bash subprocess triggers a "GLIBC_2.34 not found" error.
    # Clearing it here keeps every child process clean on any node.
    preload = os.environ.pop("LD_PRELOAD", None)
    if preload:
        logging.getLogger("emtrace").debug(
            "Cleared inherited LD_PRELOAD for subprocesses: %s", preload,
        )


# ------------------------------------------------------------------
# run
# ------------------------------------------------------------------

def _load_config(config_file, profile):
    """Build a config from --config (exact) or --profile (overlay), else defaults.

    The two differ on purpose. `--config` replaces the bundled defaults
    outright, which is what a frozen, fully-specified config file needs.
    `--profile` deep-merges a small bundled overlay onto default.yaml, so a
    profile only states what it changes. They are mutually exclusive because
    combining them has no single sensible meaning.
    """
    if config_file and profile:
        raise click.UsageError(
            "--config and --profile are mutually exclusive: --config replaces "
            "the bundled defaults, --profile overlays them. Pick one."
        )
    if config_file:
        return EMtraceConfig.from_yaml(config_file)
    if profile:
        return EMtraceConfig.from_profile(profile)
    return EMtraceConfig.create_default()


@cli.command()
# --- Required ---
@click.option(
    "-1", "--reads", "reads_file", required=True,
    type=click.Path(exists=True),
    help="Input FASTQ file (plain or gzipped).",
)
@click.option(
    "-s", "--sample", "sample_name", required=True,
    help="PipelineJob name / label.",
)
# --- Output paths ---
@click.option(
    "-r", "--res", "res_dir", default=None,
    type=click.Path(),
    help="Root results directory [current dir].",
)
@click.option(
    "--logs", "logs_dir", default=None,
    type=click.Path(),
    help="Root logs directory [logs].",
)
# --- Pipeline control ---
@click.option(
    "-p", "--pipe", "pipeline", default="main",
    show_default=True,
    help=(
        "Workflow name or comma-separated steps "
        "(pore,chop,flt,ann,em,ass,rep,plot,pack,clean). "
        "Workflow names: main, annotate, report, realtime."
    ),
)
@click.option(
    "--append", "append", is_flag=True, default=False,
    help="Append annotation to existing annot1/annot2 TSV.",
)
# --- Assembly ---
@click.option(
    "--refseqs", default=None,
    help="Comma-separated FASTA paths for reference-based assembly (overrides config).",
)
@click.option(
    "--refaccs", default=None,
    help="Comma-separated NCBI accessions for reference-based assembly (overrides config).",
)
@click.option(
    "--ref_auto/--no-ref_auto", "ref_auto", default=None,
    help="Select assembly references automatically when --refseqs/--refaccs are not given "
         "(assembly.ref_auto.enabled).",
)
@click.option(
    "--ref_auto_min_reads", type=int, default=None,
    help="Min distinct reads for a species or an extra accession in automatic reference selection.",
)
@click.option(
    "--secondary_margin", type=float, default=None,
    help="How much better an annotation hit must align to replace an assembly-DB pick "
         "(assembly.ref_auto.secondary_margin).",
)
@click.option(
    "--ref_auto_assembly_db/--no-ref_auto_assembly_db", "ref_auto_assembly_db", default=None,
    help="Search assembly.database for reference candidates (assembly.ref_auto.use_assembly_db).",
)
@click.option(
    "--ref_auto_annotation/--no-ref_auto_annotation", "ref_auto_annotation", default=None,
    help="Use annotation hits as reference candidates (assembly.ref_auto.use_annotation).",
)
# --- Host filtering ---
@click.option(
    "--hostgen", default=None,
    help=(
        "Host genome key(s) or FASTA path(s), "
        "comma-separated (e.g. Homo_sapiens)."
    ),
)
@click.option(
    "--bgfilter", default=None,
    help=(
        "Alias for --hostgen: bgfilter database key(s) or FASTA path(s), "
        "comma-separated. Takes precedence over --hostgen if both are given."
    ),
)
# --- Annotation ---
@click.option(
    "--ann1", default=None,
    help=(
        "Comma-separated annotation database keys for "
        "round 1 (e.g. minimap.nt.abv)."
    ),
)
@click.option(
    "--ann2", default=None,
    help=(
        "Comma-separated target:db pairs for round 2 "
        "(e.g. vi:blastn.nt,ab:diamondx.nr)."
    ),
)
@click.option(
    "--anns", default=None,
    help=(
        "Named annotation strategy from config "
        "(e.g. vi.nt). Overrides --ann1/--ann2."
    ),
)
@click.option(
    "--tophits", is_flag=True, default=False,
    help=(
        "Keep only each read's top-scoring hits (ties retained) after every "
        "search, as annotation.filter_tophits: true. By default all hits are "
        "kept so EM can weigh near-best alternatives; read partitions and the "
        "annotation/abundance tables use top hits either way."
    ),
)
# --- Scoring thresholds ---
@click.option(
    "--min_blastn_bits", type=int, default=None,
    help="Minimum bitscore for blastn hits.",
)
@click.option(
    "--min_minimap_score", type=int, default=None,
    help="Minimum DP peak score for minimap2 hits.",
)
@click.option(
    "--min_diamond_bits", type=int, default=None,
    help="Minimum bitscore for diamond hits.",
)
@click.option(
    "--min_qcov", type=float, default=None,
    help=(
        "Minimum query coverage fraction [0-1] for the annotation and "
        "abundance tables (report.min_qcov). Hits with unknown coverage "
        "(minimap2 without -c) are not filtered."
    ),
)
@click.option(
    "--min_coverage", type=int, default=None,
    help="Minimum read coverage for assembly masking.",
)
@click.option(
    "--medaka/--no-medaka", "medaka", default=None,
    help="Polish the assembly consensus with medaka (assembly.medaka.enabled; ONT only).",
)
@click.option(
    "--medaka_model", default=None,
    help="Medaka model name or path, or 'auto' to read it from the read headers "
         "(assembly.medaka.model).",
)
@click.option(
    "--min_af", type=float, default=None,
    help="Minimum allele frequency for an allele to be written into the consensus.",
)
@click.option(
    "--report_min_af", type=float, default=None,
    help="Minimum allele frequency for a variant to be listed in the variant table.",
)
@click.option(
    "--min_read2hostgen_score", type=int, default=None,
    help="Post-alignment minimum AS score for host read filtering.",
)
# --- Reporting ---
@click.option(
    "--tail", "report_tail", type=float, default=None,
    help=(
        "Cumulative-share tail trim as a FRACTION [0-1] for the abundance "
        "table (report.tail): 0.01 drops the smallest taxa that together "
        "hold the last 1% of reads, 0 disables. The surviving taxa also "
        "define annot_table, the sorted reads, the division stats and the "
        "plots, so every report agrees; the raw annotation/annot*.tsv "
        "tables are unaffected. Not related to --tailcrop, which trims "
        "bases from read ends."
    ),
)
@click.option(
    "--tail-viruses/--no-tail-viruses", "report_tail_viruses", default=None,
    help=(
        "Also trim the viral division (report.tail_viruses). Off by "
        "default: the cumulative tail is where rare viral true positives "
        "live. Viral reads count toward the share either way."
    ),
)
@click.option(
    "--nospnorep", is_flag=True, default=False,
    help=(
        "Exclude hits with no species-level assignment (empty, or "
        "non-binomial such as 'Salmonella sp.') from the EM input and "
        "from the annotation/abundance reports. The raw "
        "annotation/annot*.tsv tables are unaffected."
    ),
)
# --- EM abundance ---
@click.option(
    "--score2prob_temp", type=float, default=None,
    help="Temperature for alignment score → probability conversion.",
)
# --- Packing ---
@click.option(
    "--pack_reads", is_flag=True, default=False,
    help="Include all reads in output tarball.",
)
@click.option(
    "--pack_reads_vi", is_flag=True, default=False,
    help="Include viral reads in output tarball.",
)
# --- Runtime ---
@click.option(
    "-t", "--numth", type=int, default=None,
    help="Number of threads.",
)
@click.option(
    "--config", "config_file", default=None,
    type=click.Path(exists=True),
    help="Config YAML file (default: bundled default.yaml).",
)
@click.option(
    "--profile", default=None,
    help=(
        "Bundled settings profile overlaid on the defaults, e.g. '16S'. "
        "Unlike --config, which replaces the defaults outright, a profile "
        "only changes the keys it names. Mutually exclusive with --config."
    ),
)
@click.option(
    "-v", "--verbose", "verbosity", count=True,
    help="Increase verbosity (-v info, -vv debug).",
)
def run(
    reads_file: str,
    sample_name: str,
    res_dir: Optional[str],
    logs_dir: Optional[str],
    pipeline: str,
    append: bool,
    refseqs: Optional[str],
    refaccs: Optional[str],
    ref_auto: Optional[bool],
    ref_auto_min_reads: Optional[int],
    secondary_margin: Optional[float],
    ref_auto_assembly_db: Optional[bool],
    ref_auto_annotation: Optional[bool],
    hostgen: Optional[str],
    bgfilter: Optional[str],
    ann1: Optional[str],
    ann2: Optional[str],
    anns: Optional[str],
    tophits: bool,
    min_blastn_bits: Optional[int],
    min_minimap_score: Optional[int],
    min_diamond_bits: Optional[int],
    min_qcov: Optional[float],
    min_coverage: Optional[int],
    medaka: Optional[bool],
    medaka_model: Optional[str],
    min_af: Optional[float],
    report_min_af: Optional[float],
    min_read2hostgen_score: Optional[int],
    report_tail: Optional[float],
    report_tail_viruses: Optional[bool],
    nospnorep: bool,
    score2prob_temp: Optional[float],
    pack_reads: bool,
    pack_reads_vi: bool,
    numth: Optional[int],
    config_file: Optional[str],
    profile: Optional[str],
    verbosity: int,
):
    """Run the EMtrace pipeline on a FASTQ file."""
    # ----------------------------------------------------------
    # Logging
    # ----------------------------------------------------------
    log_level = (
        logging.DEBUG if verbosity >= 2 else
        logging.INFO if verbosity >= 1 else
        logging.WARNING
    )
    logger = setup_logging(level=log_level)
    quiet_external_loggers()

    # ----------------------------------------------------------
    # Config
    # ----------------------------------------------------------
    try:
        cfg = _load_config(config_file, profile)
    except Exception as exc:
        logger.error("Failed to load config: %s", exc)
        sys.exit(1)

    # Apply CLI overrides to config
    if numth is not None:
        cfg.runtime.num_threads = numth
    if refseqs is not None:
        cfg.assembly.refseqs = refseqs
    if refaccs is not None:
        cfg.assembly.refaccs = refaccs
    ra = cfg.assembly.ref_auto
    if ref_auto is not None:
        ra.enabled = ref_auto
    if ref_auto_min_reads is not None:
        ra.min_reads = ref_auto_min_reads
    if secondary_margin is not None:
        ra.secondary_margin = secondary_margin
    if ref_auto_assembly_db is not None:
        ra.use_assembly_db = ref_auto_assembly_db
    if ref_auto_annotation is not None:
        ra.use_annotation = ref_auto_annotation
    if not (ra.use_assembly_db or ra.use_annotation):
        raise click.BadParameter(
            "--no-ref_auto_assembly_db and --no-ref_auto_annotation cannot be combined: "
            "automatic reference selection needs at least one hit source."
        )

    # Scoring thresholds
    if min_blastn_bits is not None:
        cfg.annotation.blastn.bitscore = min_blastn_bits
    if min_minimap_score is not None:
        cfg.annotation.minimap.s = min_minimap_score
    if min_diamond_bits is not None:
        cfg.annotation.diamond_blastp.bitscore = min_diamond_bits
    if min_qcov is not None:
        cfg.report.min_qcov = min_qcov
    if min_coverage is not None:
        cfg.assembly.min_coverage = min_coverage
    if medaka is not None:
        cfg.assembly.medaka.enabled = medaka
    if medaka_model is not None:
        cfg.assembly.medaka.model = medaka_model
    if min_af is not None:
        cfg.assembly.min_af = min_af
    if report_min_af is not None:
        cfg.assembly.report_min_af = report_min_af
    if min_read2hostgen_score is not None:
        cfg.bgfilter.min_score = min_read2hostgen_score

    # Packing
    if pack_reads:
        cfg.output.pack_reads = True
    if pack_reads_vi:
        cfg.output.pack_reads_vi = True

    # Cumulative-share tail trim. CLI beats the config file, which beats the
    # model default (off).
    if report_tail is not None:
        cfg.report.tail = report_tail
    if report_tail_viruses is not None:
        cfg.report.tail_viruses = report_tail_viruses

    # No-species filtering. One switch drives both consumers: the EM input
    # and the reported annotation/abundance tables. Set on cfg (not on
    # sample.metadata) because the steps read their settings from config.
    if nospnorep:
        cfg.em.exclude_no_species = True
        cfg.report.exclude_no_species = True

    # Top-hit filtering at annotation (opt-in; YAML annotation.filter_tophits).
    if tophits:
        cfg.annotation.filter_tophits = True

    # ----------------------------------------------------------
    # Output directories
    # ----------------------------------------------------------
    out_root = Path(res_dir) if res_dir else Path.cwd()
    sample_outdir = out_root / sample_name
    effective_logs = logs_dir or "logs"

    # ----------------------------------------------------------
    # File logging (to {logs}/{sample}/run.log)
    # ----------------------------------------------------------
    log_path = (
        Path(effective_logs) / sample_name / "run.log"
    )
    setup_logging(
        level=log_level,
        log_file=log_path,
        log_file_level=logging.DEBUG,
    )

    # ----------------------------------------------------------
    # Build PipelineJob
    # ----------------------------------------------------------
    sample = PipelineJob(
        name=sample_name,
        input_reads=Path(reads_file),
        output_dir=sample_outdir,
        logs_dir=effective_logs,
        results_dir="results",
        pipe=pipeline,
        append=append,
        num_threads=cfg.runtime.num_threads,
        verbose=verbosity,
    )
    if score2prob_temp is not None:
        sample.metadata["score2prob_temp"] = score2prob_temp

    sample.init_dirs()

    # ----------------------------------------------------------
    # Resolve annotation strategy + host databases
    # ----------------------------------------------------------
    pipe_obj = Pipeline(cfg)
    try:
        pipe_obj.setup_job(
            sample,
            anns=anns or "",
            ann1=ann1 or "",
            ann2=ann2 or "",
            bgfilter=bgfilter or hostgen or "",
        )
    except ConfigurationError as exc:
        logger.error("Configuration error: %s", exc)
        sys.exit(1)

    # ----------------------------------------------------------
    # Parse steps
    # ----------------------------------------------------------
    try:
        steps = Pipeline.parse_steps(pipeline)
    except Exception as exc:
        logger.error("Invalid pipeline steps: %s", exc)
        sys.exit(1)

    # ----------------------------------------------------------
    # Run
    # ----------------------------------------------------------
    logger.info(
        "emtrace %s — sample: %s", __version__, sample_name,
    )
    logger.info("Input:  %s", reads_file)
    logger.info("Output: %s", sample_outdir)
    logger.info(
        "Steps:  %s", " -> ".join(s.value for s in steps),
    )

    try:
        pipe_obj.run(sample, steps)
    except Exception as exc:
        logger.error("Pipeline failed: %s", exc)
        if verbosity >= 2:
            import traceback
            traceback.print_exc()
        sys.exit(1)

    logger.info("Done.")


# ------------------------------------------------------------------
# databases
# ------------------------------------------------------------------

@cli.command(name="databases")
@click.option(
    "--list", "list_dbs", is_flag=True,
    help="List all configured databases and installed status.",
)
@click.option(
    "--install", "install_key", default=None,
    metavar="KEY",
    help=(
        "Download and install a database by key "
        "(e.g. minimap.nt.abv, Homo_sapiens, "
        "taxonomy, hostdbs)."
    ),
)
@click.option(
    "--update", "update_key", default=None,
    metavar="KEY|all",
    help=(
        "Re-download a database (force overwrite). "
        "Use 'all' to update every installed database."
    ),
)
@click.option(
    "--remove", "remove_key", default=None,
    metavar="KEY",
    help="Remove an installed database by key.",
)
@click.option(
    "--build-index", "build_index_key", default=None,
    metavar="KEY|hostdbs",
    help=(
        "Pre-build the minimap2 .mmi index for a host filter (downloading the "
        "FASTA first if needed). Use 'hostdbs' for all host genomes. Run this "
        "serially before HPC array jobs so tasks reuse one index instead of "
        "racing to build it."
    ),
)
@click.option(
    "--force", is_flag=True, default=False,
    help="Overwrite existing files (--install / --build-index).",
)
@click.option(
    "--config", "config_file", default=None,
    type=click.Path(exists=True),
    help="Config YAML file.",
)
@click.option(
    "-v", "--verbose", "verbosity", count=True,
    help="Increase verbosity.",
)
def databases(
    list_dbs: bool,
    install_key: Optional[str],
    update_key: Optional[str],
    remove_key: Optional[str],
    build_index_key: Optional[str],
    force: bool,
    config_file: Optional[str],
    verbosity: int,
):
    """Manage reference databases."""
    log_level = (
        logging.DEBUG if verbosity >= 2 else
        logging.INFO if verbosity >= 1 else
        logging.WARNING
    )
    setup_logging(level=log_level)

    try:
        cfg = (
            EMtraceConfig.from_yaml(config_file)
            if config_file
            else EMtraceConfig.create_default()
        )
    except Exception as exc:
        click.echo(f"Error loading config: {exc}", err=True)
        sys.exit(1)

    if list_dbs:
        _list_databases(cfg)

    if install_key:
        _do_install(cfg, install_key, force=force)

    if update_key:
        # --update is --install --force
        if update_key.lower() == "all":
            ann_rows, host_rows = database_status(cfg)
            for row in ann_rows:
                if row["installed"]:
                    _do_install(
                        cfg, row["key"], force=True,
                    )
            for row in host_rows:
                if row["installed"]:
                    _do_install(
                        cfg, row["key"], force=True,
                    )
        else:
            _do_install(cfg, update_key, force=True)

    if remove_key:
        _do_remove(cfg, remove_key)

    if build_index_key:
        _do_build_index(cfg, build_index_key, force=force)

    if not any([list_dbs, install_key, update_key, remove_key, build_index_key]):
        click.echo(
            "No action specified. Use --list, --install, "
            "--update, --remove, or --build-index.",
        )


def _list_databases(cfg: EMtraceConfig) -> None:
    """Print installed/available database status to stdout."""
    ann_rows, host_rows = database_status(cfg)

    click.echo("\nAnnotation Databases\n" + "-" * 72)
    click.echo(
        f"  {'KEY':<28} {'SEARCH':<12} {'STATUS':<10} NAME",
    )
    click.echo("  " + "-" * 68)
    for row in ann_rows:
        status = "installed" if row["installed"] else "not installed"
        url_tag = "" if row["url"] not in ("NA", "", "0") else " [no URL]"
        click.echo(
            f"  {row['key']:<28} {row['search']:<12} "
            f"{status:<10} {row['name']}{url_tag}",
        )

    click.echo("\nHost Genome Databases\n" + "-" * 72)
    click.echo(
        f"  {'KEY':<28} {'STATUS':<14} SPECIES",
    )
    click.echo("  " + "-" * 68)
    for row in host_rows:
        status = "installed" if row["installed"] else "not installed"
        click.echo(
            f"  {row['key']:<28} {status:<14} {row['name']}",
        )

    click.echo("\nAnnotation Strategies\n" + "-" * 72)
    for key in cfg.list_strategies():
        click.echo(f"  {key:<20} {cfg.ann_strategies[key]}")
    click.echo()


def _do_install(
    cfg: EMtraceConfig,
    key: str,
    force: bool = False,
) -> None:
    """Run install_by_key with CLI progress output."""
    def _progress(done: int, total: int) -> None:
        if total > 0:
            pct = done / total * 100
            mb = done / 1e6
            click.echo(
                f"\r  {pct:5.1f}%  {mb:.1f} MB",
                nl=False,
            )
        elif total == -1:
            click.echo(
                f"\r  {done / 1e6:.1f} MB", nl=False,
            )

    click.echo(f"Installing '{key}'...")
    try:
        install_by_key(
            cfg, key, progress_cb=_progress, force=force,
        )
        click.echo(f"\nDone: '{key}' installed.")
    except DatabaseAlreadyInstalledError as exc:
        click.echo(
            f"\nSkipped: {exc}\nUse --force to overwrite.",
            err=True,
        )
        sys.exit(1)
    except DatabaseNotConfiguredError as exc:
        click.echo(f"\nError: {exc}", err=True)
        sys.exit(1)
    except DatabaseInstallError as exc:
        click.echo(f"\nInstall failed: {exc}", err=True)
        sys.exit(1)
    except Exception as exc:
        click.echo(f"\nUnexpected error: {exc}", err=True)
        sys.exit(1)


def _do_build_index(
    cfg: EMtraceConfig,
    key: str,
    force: bool = False,
) -> None:
    """Pre-build host-filter minimap2 indices with CLI progress output."""
    def _progress(done: int, total: int) -> None:
        if total > 0:
            click.echo(f"\r  {done / total * 100:5.1f}%  {done / 1e6:.1f} MB", nl=False)
        elif total == -1:
            click.echo(f"\r  {done / 1e6:.1f} MB", nl=False)

    keys = None if key.lower() == "hostdbs" else [key]
    click.echo(f"Building host-filter index for '{key}'...")
    try:
        prebuild_host_filters(cfg, keys, force=force, progress_cb=_progress)
        click.echo(f"\nDone: index build for '{key}' complete.")
    except DatabaseNotConfiguredError as exc:
        click.echo(f"\nError: {exc}", err=True)
        sys.exit(1)
    except DatabaseInstallError as exc:
        click.echo(f"\nIndex build failed: {exc}", err=True)
        sys.exit(1)


def _do_remove(cfg: EMtraceConfig, key: str) -> None:
    """Remove database files matching the config entry for *key*."""
    from emtrace.utils.db_download import (
        _find_host_key,
    )

    db_path = None

    if key in cfg.databases:
        db_path = Path(cfg.databases[key].db)
    else:
        host_key = _find_host_key(cfg, key)
        if host_key:
            db_path = Path(cfg.host_databases[host_key].db)

    if db_path is None:
        click.echo(
            f"Error: '{key}' not found in config.", err=True,
        )
        sys.exit(1)

    files = list_installed_files(db_path)
    if not files:
        click.echo(f"Nothing to remove for '{key}'.")
        return

    click.echo(f"Removing {len(files)} file(s) for '{key}':")
    for f in files:
        click.echo(f"  {f}")

    if not click.confirm("Proceed?"):
        click.echo("Aborted.")
        return

    n = remove_database(db_path)
    click.echo(f"Removed {n} file(s).")


# ------------------------------------------------------------------
# watch (real-time mode)
# ------------------------------------------------------------------

@cli.command()
@click.option(
    "--dir", "watch_dir", required=True,
    type=click.Path(exists=True),
    help="Directory to monitor for FASTQ files.",
)
@click.option(
    "-s", "--sample", "sample_name", default="realtime",
    show_default=True,
    help="PipelineJob name / label.",
)
@click.option(
    "-r", "--res", "res_dir", default=None,
    type=click.Path(),
    help="Root output directory [current dir].",
)
# --- Annotation ---
@click.option(
    "--ann1", default=None,
    help="Annotation database keys for round 1.",
)
@click.option(
    "--ann2", default=None,
    help="target:db pairs for round 2.",
)
@click.option(
    "--anns", default=None,
    help="Named annotation strategy.",
)
# --- Host filtering ---
@click.option(
    "--hostgen", default=None,
    help="Host genome key(s), comma-separated.",
)
# --- Batch ---
@click.option(
    "--batch-size", type=int, default=1,
    show_default=True,
    help="Number of files per processing batch.",
)
@click.option(
    "--debounce", type=float, default=2.0,
    show_default=True,
    help=(
        "Seconds to wait after a file stops changing "
        "before processing."
    ),
)
@click.option(
    "--scan-existing", is_flag=True, default=False,
    help="Process FASTQ files already in the directory.",
)
@click.option(
    "--report", "run_report", is_flag=True, default=False,
    help="Run report step after each batch.",
)
# --- Runtime ---
@click.option(
    "-t", "--numth", type=int, default=None,
    help="Number of threads.",
)
@click.option(
    "--config", "config_file", default=None,
    type=click.Path(exists=True),
    help="Config YAML file.",
)
@click.option(
    "--profile", default=None,
    help=(
        "Bundled settings profile overlaid on the defaults, e.g. '16S'. "
        "Unlike --config, which replaces the defaults outright, a profile "
        "only changes the keys it names. Mutually exclusive with --config."
    ),
)
@click.option(
    "-v", "--verbose", "verbosity", count=True,
    help="Increase verbosity.",
)
def watch(
    watch_dir: str,
    sample_name: str,
    res_dir: Optional[str],
    ann1: Optional[str],
    ann2: Optional[str],
    anns: Optional[str],
    hostgen: Optional[str],
    batch_size: int,
    debounce: float,
    scan_existing: bool,
    run_report: bool,
    numth: Optional[int],
    config_file: Optional[str],
    profile: Optional[str],
    verbosity: int,
):
    """Monitor a directory for FASTQ files and process in real-time.

    Watches the specified directory for new .fastq/.fastq.gz files
    (e.g. MinKNOW output), processes each batch through preprocessing,
    host filtering, annotation, and EM abundance estimation.

    Press Ctrl+C to stop — a final report will be generated from
    all accumulated data.

    \b
    Example:
        emtrace watch --dir /path/to/minknow/output \\
            --hostgen Homo_sapiens --anns vi.nt -v
    """
    import signal
    import threading

    log_level = (
        logging.DEBUG if verbosity >= 2 else
        logging.INFO if verbosity >= 1 else
        logging.WARNING
    )
    log = setup_logging(level=log_level)
    quiet_external_loggers()

    # Config
    try:
        cfg = _load_config(config_file, profile)
    except Exception as exc:
        log.error("Failed to load config: %s", exc)
        sys.exit(1)

    if numth is not None:
        cfg.runtime.num_threads = numth

    # Output directory
    out_root = Path(res_dir) if res_dir else Path.cwd()
    output_dir = out_root / sample_name

    # Stop event for graceful shutdown
    stop_event = threading.Event()

    def _on_signal(signum, frame):
        click.echo(
            "\nShutting down (generating final report)...",
        )
        stop_event.set()

    signal.signal(signal.SIGINT, _on_signal)
    signal.signal(signal.SIGTERM, _on_signal)

    # Console callback: print running abundance table
    def _on_batch(batch_num, state):
        snap = state.get_snapshot()
        click.echo(
            f"\n--- Batch {snap['batch_count']} ---"
        )
        click.echo(
            f"Files: {snap['total_files']}  "
            f"Reads: {snap['total_reads']}  "
            f"Taxa: {snap['n_taxa']}  "
            f"Annotations: {snap['n_annotations']}"
        )
        if snap["log_L"] is not None:
            click.echo(
                f"EM logL: {snap['log_L']:.4f}"
            )
        if snap["top_taxa"]:
            click.echo("Top taxa:")
            for taxid, prob in snap["top_taxa"][:10]:
                click.echo(
                    f"  {taxid:<12} {prob:.6f}"
                )

    # Run
    log.info(
        "emtrace watch — sample: %s", sample_name,
    )
    log.info("Watch dir: %s", watch_dir)
    log.info("Output:    %s", output_dir)

    from emtrace.realtime.processor import (
        run_realtime,
    )

    state = run_realtime(
        watch_dir=watch_dir,
        output_dir=output_dir,
        config=cfg,
        job_name=sample_name,
        ann1=ann1 or cfg.annotation.ann1,
        ann2=ann2 or cfg.annotation.ann2,
        anns=anns or "",
        hostgen=hostgen or cfg.bgfilter.bgfilter,
        batch_size=batch_size,
        debounce=debounce,
        scan_existing=scan_existing,
        on_batch_complete=_on_batch,
        stop_event=stop_event,
        run_report=run_report,
    )

    # Summary
    snap = state.get_snapshot()
    click.echo("\n=== FINAL SUMMARY ===")
    click.echo(
        f"Batches: {snap['batch_count']}  "
        f"Files: {snap['total_files']}  "
        f"Reads: {snap['total_reads']}"
    )
    click.echo(
        f"Taxa: {snap['n_taxa']}  "
        f"Annotations: {snap['n_annotations']}"
    )
    if snap["log_L"] is not None:
        click.echo(f"Final EM logL: {snap['log_L']:.4f}")
    click.echo(f"Results in: {output_dir}")


# ------------------------------------------------------------------
# batch-summary
# ------------------------------------------------------------------

@cli.command("batch-summary")
@click.argument("resdir", type=click.Path(exists=True, file_okay=False, path_type=Path))
@click.option(
    "--output", "-o",
    type=click.Path(path_type=Path),
    default=None,
    help="Output file (.tsv, .csv, .xlsx). Default: TSV to stdout.",
)
@click.option(
    "--max-samples",
    default=None,
    type=int,
    metavar="N",
    help="Limit to the first N samples (useful for testing).",
)
@click.option("-v", "--verbose", is_flag=True)
def batch_summary_cmd(
    resdir: Path,
    output: Optional[Path],
    max_samples: Optional[int],
    verbose: bool,
):
    """Cross-sample QC and taxonomic composition summary (one row per sample).

    Reads reports/stats_fastq.csv and reports/abund_table.tsv from each
    subdirectory of RESDIR and produces one row per sample with read-count
    QC columns (reads, reads.flt, reads.hgflt, n50) and per-division read
    counts (reads.Viruses, reads.Bacteria, …).
    """
    from emtrace.batch.summary import build_summary

    log_level = logging.DEBUG if verbose else logging.INFO
    setup_logging(level=log_level)

    df = build_summary(resdir, max_samples=max_samples, verbose=verbose)

    if df.empty:
        click.echo("No data found.", err=True)
        sys.exit(1)

    _write_batch_output(df, output)


# ------------------------------------------------------------------
# batch-taxa
# ------------------------------------------------------------------

@cli.command("batch-taxa")
@click.argument("resdir", type=click.Path(exists=True, file_okay=False, path_type=Path))
@click.option(
    "--target", "-t",
    default="vi",
    show_default=True,
    type=click.Choice(["vi", "ba", "ph", "inv"]),
    help="Target organism type: vi=Viruses, ba=Bacteria, ph=Phages, "
         "inv=Invertebrates (protozoa + helminths; eukaryotic parasites, "
         "excluding fungi).",
)
@click.option(
    "--output", "-o",
    type=click.Path(path_type=Path),
    default=None,
    help="Output file (.tsv, .csv, .xlsx). Default: TSV to stdout.",
)
@click.option(
    "--max-samples",
    default=None,
    type=int,
    metavar="N",
    help="Limit to the first N samples.",
)
@click.option("-v", "--verbose", is_flag=True)
def batch_taxa_cmd(
    resdir: Path,
    target: str,
    output: Optional[Path],
    max_samples: Optional[int],
    verbose: bool,
):
    """Cross-sample species abundance table, filtered by target organism type.

    Reads reports/abund_table.tsv from each subdirectory of RESDIR and
    produces one row per taxon per sample; readn is the EM-adjusted read
    count (readn_em, rounded; raw readn where missing).  Use --target to restrict to
    Viruses (vi), Bacteria (ba), Phages (ph), or Invertebrates (inv;
    eukaryotic parasites/pathogens such as protozoa and helminths, excluding
    fungi).
    """
    from emtrace.batch.taxa import build_taxa_summary

    log_level = logging.DEBUG if verbose else logging.INFO
    setup_logging(level=log_level)

    df = build_taxa_summary(resdir, target=target, max_samples=max_samples, verbose=verbose)

    if df.empty:
        click.echo(
            f"No data found for target '{target}'.", err=True,
        )
        sys.exit(1)

    _write_batch_output(df, output)


# ------------------------------------------------------------------
# batch-abundance
# ------------------------------------------------------------------

@cli.command("batch-abundance")
@click.argument("resdir", type=click.Path(exists=True, file_okay=False, path_type=Path))
@click.option(
    "--outdir", "-o",
    type=click.Path(file_okay=False, path_type=Path),
    required=True,
    help="Output directory for the per-division abundance tables (CSV).",
)
@click.option(
    "--value",
    type=click.Choice(["readn", "prob", "both"]),
    default="both",
    show_default=True,
    help="Which matrices to write: readn=EM-adjusted read counts (raw readn "
         "where readn_em is missing), prob=EM probability, both.",
)
@click.option(
    "--xlsx",
    is_flag=True,
    default=False,
    help="Also write a single styled Excel workbook (read counts) with one "
         "sheet per division (requires openpyxl).",
)
@click.option(
    "--config", "config_file", default=None,
    type=click.Path(exists=True, path_type=Path),
    help="Config YAML; report.parasite_genera overrides the built-in "
         "invertebrate-parasite genus allowlist.",
)
@click.option(
    "--max-samples",
    default=None,
    type=int,
    metavar="N",
    help="Limit to the first N samples.",
)
@click.option("-v", "--verbose", is_flag=True)
def batch_abundance_cmd(
    resdir: Path,
    outdir: Path,
    value: str,
    xlsx: bool,
    config_file: Optional[Path],
    max_samples: Optional[int],
    verbose: bool,
):
    """Cross-sample OTU/feature tables: rows = species, columns = samples.

    Reads reports/abund_table.tsv from each subdirectory of RESDIR and writes
    one table per division group (bacteria, bacteriophages, viruses,
    invertebrates_parasites, invertebrates_other) under OUTDIR. Samples are
    ordered by batch; species rows by total read count. Two matrices per group:
    EM-adjusted read counts (.counts.csv; readn_em, raw readn where missing)
    and the EM assignment probability (.prob.csv; NOT a statistical p-value).
    """
    from emtrace.batch.abundance import build_abundance, write_abundance_xlsx

    log_level = logging.DEBUG if verbose else logging.INFO
    setup_logging(level=log_level)

    parasite_genera = None
    if config_file is not None:
        cfg = EMtraceConfig.from_yaml(config_file)
        parasite_genera = cfg.report.parasite_genera or None

    tables = build_abundance(
        resdir, max_samples=max_samples, verbose=verbose,
        parasite_genera=parasite_genera,
    )
    if not tables:
        click.echo("No abundance data found.", err=True)
        sys.exit(1)

    outdir.mkdir(parents=True, exist_ok=True)
    kinds = ["counts", "prob"] if value == "both" else \
        (["counts"] if value == "readn" else ["prob"])

    n_written = 0
    for group, mats in tables.items():
        for kind in kinds:
            df = mats[kind]
            path = outdir / f"abundance.{group}.{kind}.csv"
            df.to_csv(path, index=False, na_rep="")
            n_samples = df.shape[1] - 3  # minus species, taxid, taxpath
            click.echo(f"Wrote {len(df)} taxa x {n_samples} samples -> {path}")
            n_written += 1

    if xlsx:
        xlsx_path = outdir / "abundance.all_divisions.counts.xlsx"
        try:
            n_sheets = write_abundance_xlsx(tables, xlsx_path)
            click.echo(f"Wrote workbook ({n_sheets} sheets) -> {xlsx_path}")
        except ImportError:
            click.echo(
                "xlsx export needs openpyxl (pip install openpyxl).", err=True,
            )

    click.echo(f"Done: {n_written} tables in {outdir}")


def _write_batch_output(df: "pd.DataFrame", output: Optional[Path]) -> None:  # noqa: F821
    """Write DataFrame to file or stdout based on output path extension."""
    if output is None:
        df.to_csv(sys.stdout, sep="\t", index=False)
        return

    suffix = output.suffix.lower()
    if suffix == ".xlsx":
        df.to_excel(output, index=False)
    elif suffix == ".csv":
        df.to_csv(output, sep=",", index=False)
    else:
        df.to_csv(output, sep="\t", index=False)
    click.echo(f"Wrote {len(df)} rows to {output}")


# ------------------------------------------------------------------
# gui
# ------------------------------------------------------------------

@cli.command()
@click.option(
    "--config", "config_file", default=None,
    type=click.Path(exists=True),
    help="Config YAML file.",
)
def gui(config_file: Optional[str]):
    """Launch the desktop GUI (Phase 3)."""
    try:
        from emtrace.gui.main_window import launch_gui  # noqa: F401
        launch_gui(config_file=config_file)
    except ImportError:
        click.echo(
            "GUI not available. Install with: "
            "pip install emtrace[gui]",
            err=True,
        )
        sys.exit(1)


# ------------------------------------------------------------------
# Entry point
# ------------------------------------------------------------------

def main():
    """Package entry point."""
    cli()


if __name__ == "__main__":
    main()
