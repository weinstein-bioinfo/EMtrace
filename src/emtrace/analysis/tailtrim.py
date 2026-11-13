"""Cumulative-share tail trim — the shared definition.

This file is duplicated **verbatim** in two packages, which cannot import each
other:

    EMtrace/src/emtrace/analysis/tailtrim.py
    EMtrace-benchmark/ngs_benchmark/tailtrim.py

`EMtrace-benchmark/tests/test_tailtrim_equivalence.py` asserts the two are
byte-identical, so a change to one is a red test until it is copied to the
other. Keep them that way: the whole point is that "1 % tail trim" means one
thing in the pipeline and in the scorers.

Semantics
---------
``tail`` is a FRACTION: 0.01 drops the smallest taxa that together hold at most
the last 1 % of reads, and nothing else.

Sort taxa by abundance descending. A row is kept when the taxa ranked ABOVE it
hold less than ``1 - tail`` of the total (an exclusive running share). This is
the same as walking up from the smallest taxon and dropping rows while their
combined share stays within ``tail``. A taxon that straddles the line is kept,
so the top taxon, a single-taxon table, or two taxa of 50 % each are never
trimmed.

(Before EMtrace v1.0.5 the share was inclusive — a row's own count counted —
and rows were kept at ``<= 1 - tail``. That dropped the straddling taxon
whatever its size, and always the last row: two taxa of 50 % each lost the
second. Benchmark scores with ``tail > 0`` computed before the fix used that
rule.)

The denominator is **every row in the frame**, always. Grouping (``group_col`` /
``group_tails``) changes only each row's *threshold*, never the denominator —
so "protect viruses" is expressed as a tail of 0.0 for the viral group, whose
rows are then never trimmed, while the remaining taxa are still measured
against the full library rather than against each other.

This is a filter over **species-granularity rows** — one row per taxon, as the
abundance table has. Applying it to an already-collapsed genus or family table
gives a different and wrong answer; collapse the survivors instead.

Not to be confused with ``preprocess.chopper.tailcrop``, which trims bases from
the ends of reads.
"""

from __future__ import annotations

from typing import Mapping, Optional

import pandas as pd

__all__ = ["resolve_count_col", "tail_trim"]


def resolve_count_col(
    df: pd.DataFrame,
    prefer: str = "readn_em",
    fallback: str = "readn",
) -> str:
    """Return the count column to filter on: `prefer` if usable, else `fallback`.

    "Usable" is: present, and holding at least one value that survives numeric
    coercion. That is what detects the literal ``"NA"`` column the pipeline
    writes for ``readn_em`` when EM did not run — every value coerces to NaN,
    so the column is rejected and raw read counts are used instead.

    The predicate is ``.any()``, not ``.all()``: a column with a few
    unparseable rows is still the right basis, and those rows count as 0.
    """
    if prefer in df.columns:
        if pd.to_numeric(df[prefer], errors="coerce").notna().any():
            return prefer
    return fallback


def tail_trim(
    df: pd.DataFrame,
    tail: float,
    count_col: Optional[str] = None,
    *,
    group_col: Optional[str] = None,
    group_tails: Optional[Mapping[str, float]] = None,
    preserve_order: bool = True,
) -> pd.DataFrame:
    """Drop the rows holding the smallest `tail` share of `count_col`.

    Parameters
    ----------
    df
        One row per taxon (species granularity). See the module docstring.
    tail
        Fraction to drop, 0.01 = 1 %. The threshold for any row not covered by
        `group_tails`. ``<= 0`` disables the trim for those rows.
    count_col
        Column to rank and accumulate on. ``None`` resolves it with
        :func:`resolve_count_col` (``readn_em``, falling back to ``readn``).
    group_col, group_tails
        Per-group thresholds: a row in group *g* is kept when the taxa ranked
        above it hold less than ``1 - group_tails[g]``, defaulting to `tail`
        for groups the mapping does not name; a tail of ``<= 0`` keeps the row. The denominator stays global either way. Pass
        ``group_tails={"Viruses": 0.0}`` to exempt viruses from the trim.
    preserve_order
        ``True`` returns the surviving rows in the input's own order, with
        every value and dtype untouched — use this when the frame is about to
        be written out. ``False`` returns them sorted by `count_col` with that
        column coerced to numeric, reproducing `filter_cell`'s return value
        exactly.

    Returns
    -------
    A new frame, or `df` itself when nothing can be trimmed.
    """
    if len(df) == 0:
        return df

    thresholds = _row_tails(df, tail, group_col, group_tails)
    # Every row exempt: return the input untouched. This is not an
    # optimisation. Without it, a tail of 0 would still coerce the count
    # column and re-sort the frame, and floating-point cumsum can put the
    # final row's share at 1.0 + 1e-16 and drop it.
    if not (thresholds > 0).any():
        return df

    if count_col is None:
        count_col = resolve_count_col(df)
    if count_col not in df.columns:
        return df

    counts = pd.to_numeric(df[count_col], errors="coerce").fillna(0)
    total = counts.sum()
    if total <= 0:
        return df

    out = df.copy()
    out[count_col] = counts
    out = out.sort_values(count_col, ascending=False)
    # Share held by the taxa ranked above each row (exclusive of the row).
    above = (out[count_col].cumsum() - out[count_col]) / total
    row_tail = thresholds.reindex(out.index)
    # Exempt rows (tail <= 0, e.g. protected viruses) are always kept.
    keep = (above < (1.0 - row_tail)) | (row_tail <= 0)

    if preserve_order:
        return df.loc[df.index.isin(out.index[keep])]
    return out[keep]


def _row_tails(
    df: pd.DataFrame,
    tail: float,
    group_col: Optional[str],
    group_tails: Optional[Mapping[str, float]],
) -> pd.Series:
    """Per-row tail threshold, `tail` except where `group_tails` overrides it."""
    thresholds = pd.Series(float(tail), index=df.index, dtype=float)
    if group_col and group_tails and group_col in df.columns:
        groups = df[group_col].astype("string")
        for group, group_tail in group_tails.items():
            thresholds = thresholds.mask(groups == group, float(group_tail))
    return thresholds
