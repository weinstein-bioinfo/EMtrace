"""Cross-sample QC and taxonomic composition summary.

Ports LazypipeX ``batch_sample_summary.pl`` to Python.
Reads ``reports/stats_fastq.csv`` and ``reports/abund_table.tsv``
from each sample subdirectory and produces one row per sample.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

from emtrace.models.abundance import read_abundance_tsv

logger = logging.getLogger("emtrace.batch.summary")

# Stages written to stats_fastq.csv by the report step
_STAGE_INPUT = "input"
_STAGE_PORECHOP = "porechop"
_STAGE_CHOPPED = "chopped"
_STAGE_HOST_FLT = "host_filtered"

# Rows with an unassigned taxonomy path come back with an empty division
# string; relabel them so they appear as a clear "reads.unassigned" column
# rather than a nameless "reads." column.
_DIVISION_UNASSIGNED = "unassigned"

# The Viruses division is always reported (zero-filled when no sample has viral
# reads) so the batch summary reliably carries a "reads.Viruses" column.
_DIVISION_VIRUSES = "Viruses"


def _natural_key(s: str) -> list:
    """Split a string into text/integer runs for natural sort."""
    return [int(c) if c.isdigit() else c.lower() for c in re.split(r"(\d+)", s)]


def find_sample_dirs(
    resdir: Path,
    max_samples: Optional[int] = None,
) -> List[Tuple[str, Path]]:
    """Return naturally-sorted [(sample_name, sample_path)] for all valid subdirs.

    A subdirectory is valid if it contains an abundance table
    (``reports/abund_table.tsv`` or a top-level ``abund_table.tsv`` legacy
    layout) **or** FASTQ QC stats (``reports/stats_fastq.csv``). The latter lets
    samples with no annotations still appear in the summary — the report step
    skips writing ``abund_table.tsv`` when the abundance table is empty, but it
    still emits ``stats_fastq.csv``. Hidden directories (names starting with
    '.') are skipped.
    """
    candidates = []
    for p in sorted(resdir.iterdir(), key=lambda x: _natural_key(x.name)):
        if not p.is_dir() or p.name.startswith("."):
            continue
        if _abund_path(p) is not None or _fastq_stats_path(p) is not None:
            candidates.append((p.name, p))

    if max_samples is not None:
        candidates = candidates[:max_samples]
    return candidates


def _abund_path(sample_dir: Path) -> Optional[Path]:
    """Return path to abund_table.tsv, checking known layout variants."""
    for candidate in (
        sample_dir / "results" / "reports" / "abund_table.tsv",  # standard run output
        sample_dir / "reports" / "abund_table.tsv",              # flat layout
        sample_dir / "abund_table.tsv",                          # legacy
    ):
        if candidate.exists():
            return candidate
    return None


def _fastq_stats_path(sample_dir: Path) -> Optional[Path]:
    for candidate in (
        sample_dir / "results" / "reports" / "stats_fastq.csv",
        sample_dir / "reports" / "stats_fastq.csv",
    ):
        if candidate.exists():
            return candidate
    return None


def load_fastq_stats(sample_dir: Path) -> Dict[str, object]:
    """Read ``reports/stats_fastq.csv`` and extract per-stage read counts.

    Returns dict with keys: reads, reads.flt, reads.hgflt, n50.
    Returns empty dict (with warning) if the file is missing.
    """
    path = _fastq_stats_path(sample_dir)
    if path is None:
        logger.warning("%s: reports/stats_fastq.csv not found", sample_dir.name)
        return {}

    try:
        df = pd.read_csv(path, dtype=str, keep_default_na=False)
    except Exception as exc:
        logger.warning("%s: could not read stats_fastq.csv: %s", sample_dir.name, exc)
        return {}

    if "stage" not in df.columns or "read_count" not in df.columns:
        logger.warning("%s: stats_fastq.csv missing expected columns", sample_dir.name)
        return {}

    df["read_count"] = pd.to_numeric(df["read_count"], errors="coerce")
    stage_map = df.set_index("stage")["read_count"].to_dict()

    def _get(stage: str) -> Optional[int]:
        v = stage_map.get(stage)
        return int(v) if pd.notna(v) else None

    result: Dict[str, object] = {
        "reads": _get(_STAGE_INPUT),
        "reads.porechop": _get(_STAGE_PORECHOP),
        "reads.chopped": _get(_STAGE_CHOPPED),
        "reads.hostflt": _get(_STAGE_HOST_FLT),
    }

    # n50 from the input stage row
    if "n50" in df.columns:
        df["n50"] = pd.to_numeric(df["n50"], errors="coerce")
        n50_map = df.set_index("stage")["n50"].to_dict()
        result["n50"] = int(n50_map[_STAGE_INPUT]) if pd.notna(n50_map.get(_STAGE_INPUT)) else None

    return result


def load_abund_table(sample_dir: Path) -> pd.DataFrame:
    """Read the abundance table for a sample.

    Returns empty DataFrame (with warning) if the file is missing.
    """
    path = _abund_path(sample_dir)
    if path is None:
        logger.warning("%s: abund_table.tsv not found", sample_dir.name)
        return pd.DataFrame()

    try:
        return read_abundance_tsv(path)
    except Exception as exc:
        logger.warning("%s: could not read abund_table.tsv: %s", sample_dir.name, exc)
        return pd.DataFrame()


def aggregate_divisions(df: pd.DataFrame) -> Dict[str, float]:
    """Sum readn per division. Returns {division_name: total_readn}.

    Rows whose division is blank (taxids that didn't resolve to an NCBI
    division) are folded into the "Unassigned" division rather than producing a
    nameless "reads." column.
    """
    if df.empty or "division" not in df.columns or "readn" not in df.columns:
        return {}
    numeric = pd.to_numeric(df["readn"], errors="coerce").fillna(0)
    division = df["division"].fillna("").astype(str).str.strip()
    division = division.where(division != "", _DIVISION_UNASSIGNED)
    return (
        df.assign(_readn=numeric, _division=division)
        .groupby("_division")["_readn"]
        .sum()
        .to_dict()
    )


def _safe_int(v) -> Optional[int]:
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def build_summary(
    resdir: Path,
    max_samples: Optional[int] = None,
    verbose: bool = False,
) -> pd.DataFrame:
    """Build the cross-sample summary DataFrame.

    Two-pass approach:
      Pass 1 — load all samples, collect per-sample data and global division totals.
      Pass 2 — assemble output rows with dynamic division columns sorted by global read count.
    """
    sample_dirs = find_sample_dirs(resdir, max_samples)
    if not sample_dirs:
        logger.warning("No valid sample directories found in %s", resdir)
        return pd.DataFrame()

    # Pass 1
    sample_data: List[Dict] = []
    global_div_totals: Dict[str, float] = {}

    try:
        from rich.progress import track as _track
        iter_samples = _track(sample_dirs, description="Reading samples…")
    except ImportError:
        iter_samples = sample_dirs

    for sample_name, sample_dir in iter_samples:
        if verbose:
            logger.info("Loading %s", sample_name)

        fq = load_fastq_stats(sample_dir)
        abund = load_abund_table(sample_dir)
        divs = aggregate_divisions(abund)

        for div, cnt in divs.items():
            global_div_totals[div] = global_div_totals.get(div, 0.0) + cnt

        sample_data.append({
            "sample": sample_name,
            "fastq": fq,
            "divs": divs,
        })

    # For an NCBI run (detected by any division data — SILVA runs report
    # 'domain' instead and yield none here), always report the Viruses division
    # so a batch with no viral reads still carries a zero-filled reads.Viruses
    # column.
    if global_div_totals:
        global_div_totals.setdefault(_DIVISION_VIRUSES, 0.0)

    # Sort divisions by ascending global read count (matching Perl behaviour)
    ordered_divs = sorted(global_div_totals, key=lambda d: global_div_totals[d])

    # Build column header names: spaces → underscores
    div_col = {d: "reads." + d.replace(" ", "_") for d in ordered_divs}

    # Pass 2 — assemble rows
    rows = []
    for sd in sample_data:
        fq = sd["fastq"]
        divs = sd["divs"]

        mapped = round(sum(divs.values()))
        hgflt = _safe_int(fq.get("reads.hostflt"))
        unmapped = (hgflt - mapped) if hgflt is not None else None

        row: Dict[str, object] = {"Sample": sd["sample"]}
        row["reads"] = _safe_int(fq.get("reads"))
        row["reads.porechop"] = _safe_int(fq.get("reads.porechop"))
        row["reads.chopped"] = _safe_int(fq.get("reads.chopped"))
        row["reads.hostflt"] = hgflt
        row["n50"] = _safe_int(fq.get("n50"))

        for div in ordered_divs:
            row[div_col[div]] = round(divs.get(div, 0.0))

        row["reads.mapped"] = mapped
        row["reads.unmapped"] = unmapped

        rows.append(row)

    return pd.DataFrame(rows)
