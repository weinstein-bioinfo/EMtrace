"""Cross-sample species-level abundance summary, filtered by organism type.

Ports LazypipeX ``batch_abund_summary.pl`` to Python.
Reads ``reports/abund_table.tsv`` from each sample subdirectory and
produces one row per taxon per sample, filtered by ``--target vi/ba/ph/inv``.
``readn`` is the EM-adjusted read count (``readn_em``, rounded); raw ``readn``
where ``readn_em`` is missing (see :func:`em_read_counts`).
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

from emtrace.models.abundance import read_abundance_tsv

logger = logging.getLogger("emtrace.batch.taxa")

# Division-level filters for each target type.
# Uses the NCBI division column value from abund_table.tsv.
_TARGET_FILTERS: Dict[str, object] = {
    "vi": lambda div: str(div).lower() == "viruses",
    "ba": lambda div: str(div).lower() == "bacteria",
    "ph": lambda div: "phage" in str(div).lower() or str(div).lower() == "phages",
    # The NCBI "Invertebrates" division collects invertebrate animals and,
    # by inheritance, the protists — so it captures most eukaryotic parasites
    # and pathogens (protozoa such as Plasmodium/Babesia/Trypanosoma and
    # helminths) for downstream surveillance summaries. Fungal pathogens are
    # excluded: they fall under the separate "Plants and Fungi" division.
    "inv": lambda div: str(div).lower() == "invertebrates",
}

# Output column order (contign / assembly.size dropped — absent from EMtrace)
_OUTPUT_COLS = [
    "Sample",
    "taxid",
    "readn",
    "prob",
    "organism",
    "species",
    "genus",
    "family",
    "order",
    "class",
    "division",
]


def _natural_key(s: str) -> list:
    return [int(c) if c.isdigit() else c.lower() for c in re.split(r"(\d+)", s)]


def find_sample_dirs(
    resdir: Path,
    max_samples: Optional[int] = None,
) -> List[Tuple[str, Path]]:
    """Return naturally-sorted [(sample_name, sample_path)] for all valid subdirs."""
    candidates = []
    for p in sorted(resdir.iterdir(), key=lambda x: _natural_key(x.name)):
        if not p.is_dir() or p.name.startswith("."):
            continue
        has_reports = (
            (p / "results" / "reports" / "abund_table.tsv").exists()
            or (p / "reports" / "abund_table.tsv").exists()
        )
        has_legacy = (p / "abund_table.tsv").exists()
        if has_reports or has_legacy:
            candidates.append((p.name, p))

    if max_samples is not None:
        candidates = candidates[:max_samples]
    return candidates


def _abund_path(sample_dir: Path) -> Optional[Path]:
    for candidate in (
        sample_dir / "results" / "reports" / "abund_table.tsv",
        sample_dir / "reports" / "abund_table.tsv",
        sample_dir / "abund_table.tsv",
    ):
        if candidate.exists():
            return candidate
    return None


def load_abund_table(sample_dir: Path) -> pd.DataFrame:
    path = _abund_path(sample_dir)
    if path is None:
        logger.warning("%s: abund_table.tsv not found", sample_dir.name)
        return pd.DataFrame()
    try:
        return read_abundance_tsv(path)
    except Exception as exc:
        logger.warning("%s: could not read abund_table.tsv: %s", sample_dir.name, exc)
        return pd.DataFrame()


def em_read_counts(df: pd.DataFrame, sample_name: str = "") -> pd.Series:
    """Per-taxon read counts: EM-adjusted ``readn_em`` rounded to whole reads.

    Falls back to raw ``readn`` for rows with no ``readn_em`` value, and for
    whole tables written before the ``readn_em`` column existed.
    """
    raw = pd.to_numeric(df["readn"], errors="coerce")
    if "readn_em" in df.columns:
        em = pd.to_numeric(df["readn_em"], errors="coerce")
    else:
        em = pd.Series(float("nan"), index=df.index)
    if em.isna().all():
        logger.info("%s: no readn_em values, using raw readn", sample_name)
    return em.round().where(em.notna(), raw).fillna(0).astype(int)


def filter_target(df: pd.DataFrame, target: str) -> pd.DataFrame:
    """Filter rows by target organism type using the division column."""
    if df.empty or "division" not in df.columns:
        return df
    fn = _TARGET_FILTERS.get(target)
    if fn is None:
        raise ValueError(f"Unknown target '{target}'. Use vi, ba, ph, or inv.")
    mask = df["division"].apply(fn)
    return df[mask].copy()


def build_taxa_summary(
    resdir: Path,
    target: str = "vi",
    max_samples: Optional[int] = None,
    verbose: bool = False,
) -> pd.DataFrame:
    """Build the cross-sample species-level summary DataFrame.

    For each sample: load abund_table → filter by target → tag with Sample name.
    Returns a long-format DataFrame concatenated across all samples.
    """
    sample_dirs = find_sample_dirs(resdir, max_samples)
    if not sample_dirs:
        logger.warning("No valid sample directories found in %s", resdir)
        return pd.DataFrame(columns=_OUTPUT_COLS)

    try:
        from rich.progress import track as _track
        iter_samples = _track(sample_dirs, description="Reading samples…")
    except ImportError:
        iter_samples = sample_dirs

    parts: List[pd.DataFrame] = []
    for sample_name, sample_dir in iter_samples:
        if verbose:
            logger.info("Loading %s", sample_name)

        df = load_abund_table(sample_dir)
        if df.empty:
            continue

        df = filter_target(df, target)
        if df.empty:
            logger.info("%s: no rows after filtering for target '%s'", sample_name, target)
            continue

        df = df.copy()
        df["readn"] = em_read_counts(df, sample_name)
        df.insert(0, "Sample", sample_name)
        parts.append(df)

    if not parts:
        return pd.DataFrame(columns=_OUTPUT_COLS)

    combined = pd.concat(parts, ignore_index=True)

    # Keep only the defined output columns (add missing ones as empty string)
    for col in _OUTPUT_COLS:
        if col not in combined.columns:
            combined[col] = ""
    return combined[_OUTPUT_COLS]
