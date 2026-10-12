"""
TSV/CSV I/O utilities — pandas-based port of Lazypipe::Utils TSV functions.

Perl originals: read_tsv2hash, read_tsv2hashtable, read_tsv2arraytable,
nlines, colind, mcolind (all in Lazypipe/Utils.pm).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def read_tsv(
    path: Path | str,
    sep: str = "\t",
    **kwargs,
) -> pd.DataFrame:
    """Read a TSV into a DataFrame (thin wrapper around pd.read_csv)."""
    return pd.read_csv(
        path, sep=sep, dtype=str, keep_default_na=False, **kwargs,
    )


def write_tsv(
    df: pd.DataFrame,
    path: Path | str,
    sep: str = "\t",
    index: bool = False,
    **kwargs,
) -> None:
    """Write a DataFrame to a TSV file."""
    df.to_csv(path, sep=sep, index=index, **kwargs)


def filter_tophits(
    df: pd.DataFrame,
    query_col: str = "qseqid",
    score_col: str = "bitscore",
    retain_ties: bool = False,
) -> pd.DataFrame:
    """Keep only the top-scoring hit per query (port of filter_tophits()).

    Input order is not assumed — the DataFrame is sorted internally.

    Parameters
    ----------
    df : pd.DataFrame
        Hit table (any row order accepted).
    query_col : str
        Column containing query sequence IDs.
    score_col : str
        Column containing alignment scores (higher = better).
    retain_ties : bool
        If True, keep all rows that share the top score for a query.
        If False (default), keep only the single top hit.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame, sorted by *query_col* then *score_col* desc.
    """
    df = df.copy()
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    df = df.sort_values(
        [query_col, score_col], ascending=[True, False],
    )

    if retain_ties:
        max_scores = df.groupby(query_col)[score_col].transform("max")
        return df[df[score_col] == max_scores].reset_index(drop=True)
    else:
        return (
            df.drop_duplicates(subset=[query_col], keep="first")
            .reset_index(drop=True)
        )


def top_hits_per_search(
    df: pd.DataFrame,
    query_col: str = "qseqid",
    score_col: str = "bitscore",
    search_col: str = "db",
) -> pd.DataFrame:
    """Keep each read's top-scoring hits, ties retained, separately per search.

    Equivalent to running ``filter_tophits(retain_ties=True)`` on each
    search's hits (grouped by *search_col*, normally the annotation ``db``),
    which is what the annotate step does when ``annotation.filter_tophits``
    is on. Consumers that need per-read top classifications (read
    partitioning, the annotation and abundance tables) call this so they
    behave identically whether or not sub-best hits were kept at annotation.
    Idempotent; input row order is preserved. Rows with a non-numeric score
    are dropped, as in ``filter_tophits``.
    """
    if df.empty or query_col not in df.columns or score_col not in df.columns:
        return df
    keys = [c for c in (search_col, query_col) if c in df.columns]
    score = pd.to_numeric(df[score_col], errors="coerce")
    best = score.groupby([df[k] for k in keys]).transform("max")
    return df[score == best]
