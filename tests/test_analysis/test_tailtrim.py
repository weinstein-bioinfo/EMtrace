"""Shared unit tests for the cumulative-share tail trim.

This file is duplicated verbatim beside the other copy of the module it tests:

    EMtrace/tests/test_analysis/test_tailtrim.py
    EMtrace-benchmark/tests/test_tailtrim.py

Only the import line differs. `test_tailtrim_equivalence.py` in the benchmark
repo asserts the two implementations agree; these tests pin what they agree on.
"""

from __future__ import annotations

import pandas as pd
import pytest

from emtrace.analysis.tailtrim import resolve_count_col, tail_trim


def frame(counts, divisions=None, col="readn"):
    """One row per taxon, named t0..tN, in the given count order."""
    data = {"taxid": [f"t{i}" for i in range(len(counts))], col: counts}
    if divisions is not None:
        data["division"] = divisions
    return pd.DataFrame(data)


# ------------------------------------------------------------------ #
# resolve_count_col
# ------------------------------------------------------------------ #

def test_prefers_readn_em_when_numeric():
    df = pd.DataFrame({"readn": [5, 3], "readn_em": [4.5, 3.5]})
    assert resolve_count_col(df) == "readn_em"


def test_falls_back_when_readn_em_is_the_literal_NA_column():
    # What the pipeline writes when EM did not run: the whole column is "NA".
    df = pd.DataFrame({"readn": [5, 3], "readn_em": ["NA", "NA"]})
    assert resolve_count_col(df) == "readn"


def test_falls_back_when_readn_em_absent():
    assert resolve_count_col(pd.DataFrame({"readn": [1]})) == "readn"


def test_one_numeric_value_is_enough():
    # The predicate is .any(), not .all(): a partly-unparseable column is
    # still the right basis, and its bad rows count as 0.
    df = pd.DataFrame({"readn": [5, 3], "readn_em": ["NA", 3.5]})
    assert resolve_count_col(df) == "readn_em"


# ------------------------------------------------------------------ #
# Boundary: drop only taxa that TOGETHER hold at most `tail`
# ------------------------------------------------------------------ #

def test_drops_the_smallest_taxa_holding_at_most_the_tail():
    # counts 1000,100,10,1 over 1111. From the bottom: 1 (0.09 %), 1+10
    # (0.99 %) fit in a 1 % tail; adding 100 (10.9 %) does not.
    out = tail_trim(frame([1000, 100, 10, 1]), 0.01)
    assert list(out["taxid"]) == ["t0", "t1"]


def test_the_straddling_taxon_is_kept():
    # A 10 % tail drops t1..t3 (111 reads, 9.99 %) and keeps the top taxon,
    # whatever its own share.
    assert list(tail_trim(frame([1000, 100, 10, 1]), 0.10)["taxid"]) == ["t0"]


def test_tiny_tail_keeps_everything_larger_than_it():
    # The smallest taxon alone is 0.09 % > 0.05 %, so nothing is dropped.
    out = tail_trim(frame([1000, 100, 10, 1]), 0.0005)
    assert list(out["taxid"]) == ["t0", "t1", "t2", "t3"]


def test_two_equal_halves_are_both_kept():
    # The pre-v1.0.5 rule dropped the second 50 % taxon.
    assert len(tail_trim(frame([500, 500]), 0.01)) == 2


def test_a_large_last_taxon_is_never_dropped():
    assert len(tail_trim(frame([600, 300, 100]), 0.01)) == 3


def test_demo5_tail():
    # Demo 5 (16S gut standard): the 1.2 % genus is kept; the 0.93 % genus and
    # three single reads (67 reads, 0.97 % together) are dropped.
    out = tail_trim(frame([6726, 81, 64, 1, 1, 1]), 0.01)
    assert list(out["taxid"]) == ["t0", "t1"]


@pytest.mark.parametrize("tail", [0.001, 0.01, 0.05, 0.2])
def test_dropped_mass_never_exceeds_the_tail(tail):
    counts = [1000, 500, 100, 50, 10, 5, 1]
    out = tail_trim(frame(counts), tail, "readn")
    dropped = sum(counts) - out["readn"].sum()
    assert dropped <= tail * sum(counts)
    # ... and the smallest survivor would not have fitted.
    assert dropped + out["readn"].min() > tail * sum(counts)


def test_uniform_counts_drop_exactly_the_last_rows():
    # 100 taxa of 1 read each: shares are 0.01 .. 1.00, tail 0.05 keeps 95.
    out = tail_trim(frame([1] * 100), 0.05)
    assert len(out) == 95


# ------------------------------------------------------------------ #
# tail <= 0 is identity, including dtypes and row order
# ------------------------------------------------------------------ #

def test_zero_tail_returns_the_input_unchanged():
    df = frame([3, 1000, 7])            # deliberately unsorted
    out = tail_trim(df, 0.0)
    pd.testing.assert_frame_equal(out, df)


def test_negative_tail_is_identity():
    df = frame([3, 1000, 7])
    pd.testing.assert_frame_equal(tail_trim(df, -1.0), df)


def test_zero_tail_does_not_coerce_the_count_column():
    # Guards the short-circuit: a tail of 0 must not rewrite "NA" as 0.0.
    df = pd.DataFrame({"taxid": ["a"], "readn": [5], "readn_em": ["NA"]})
    out = tail_trim(df, 0.0, "readn_em")
    assert out["readn_em"].tolist() == ["NA"]


def test_empty_frame_is_identity():
    df = frame([])
    pd.testing.assert_frame_equal(tail_trim(df, 0.01), df)


def test_all_zero_counts_is_identity():
    df = frame([0, 0, 0])
    pd.testing.assert_frame_equal(tail_trim(df, 0.01), df)


def test_missing_count_column_is_identity():
    df = frame([5, 3])
    pd.testing.assert_frame_equal(tail_trim(df, 0.01, "nosuchcol"), df)


def test_single_row_survives_any_tail_below_one():
    # A one-taxon table has no tail to trim.
    assert len(tail_trim(frame([42]), 0.5)) == 1


# ------------------------------------------------------------------ #
# preserve_order
# ------------------------------------------------------------------ #

def test_preserve_order_keeps_input_order_and_values():
    df = pd.DataFrame({
        "taxid": ["a", "b", "c", "d"],
        "readn": [10, 1000, 1, 100],        # unsorted
        "readn_em": ["10.0", "1000.0", "1.0", "100.0"],  # strings on purpose
    })
    out = tail_trim(df, 0.01, "readn", preserve_order=True)
    assert list(out["taxid"]) == ["b", "d"]            # input order preserved
    assert out["readn_em"].tolist() == ["1000.0", "100.0"]  # strings untouched


def test_preserve_order_false_sorts_and_coerces():
    df = pd.DataFrame({"taxid": ["a", "b", "c"], "readn": ["10", "1000", "1"]})
    out = tail_trim(df, 0.002, "readn", preserve_order=False)
    assert list(out["taxid"]) == ["b", "a"]            # descending by count
    assert out["readn"].tolist() == [1000.0, 10.0]     # coerced to numeric


def test_both_orders_keep_the_same_rows():
    df = frame([7, 1000, 3, 250, 1])
    a = tail_trim(df, 0.02, "readn", preserve_order=True)
    b = tail_trim(df, 0.02, "readn", preserve_order=False)
    assert set(a["taxid"]) == set(b["taxid"])


# ------------------------------------------------------------------ #
# Grouping: per-group THRESHOLD, global DENOMINATOR
# ------------------------------------------------------------------ #

VI, BA = "Viruses", "Bacteria"


def test_group_tail_zero_protects_every_row_in_that_group():
    # The two viral taxa are the smallest in the library and would be the
    # first thing a 1 % trim removed. Exempting them keeps both.
    df = frame([1000, 100, 2, 1], [BA, BA, VI, VI])
    out = tail_trim(df, 0.01, "readn", group_col="division",
                    group_tails={VI: 0.0})
    assert set(out[out["division"] == VI]["taxid"]) == {"t2", "t3"}


def test_group_tail_applies_when_not_exempted():
    df = frame([1000, 100, 2, 1], [BA, BA, VI, VI])
    out = tail_trim(df, 0.01, "readn", group_col="division",
                    group_tails={VI: 0.01})
    assert VI not in set(out["division"])


def test_denominator_is_global_not_per_group():
    # THE design property, and its sharp edge. Bacterial counts are identical
    # in both frames; the second adds a large PROTECTED viral mass.
    #
    # Because the share is measured against the whole library, the viral mass
    # is consumed first and pushes every bacterial row's cumulative share up,
    # so MORE of them fall past the threshold -- here, all of them. A
    # per-group denominator would have kept the same three either way.
    #
    # In other words: in a virus-dominated sample, a global tail can remove
    # the entire bacterial background. That is the literal meaning of "drop
    # the taxa holding the last N % of the library", and it is why `tail`
    # ships at 0 and is switched on per assay.
    bact = [100, 10, 5, 3, 2, 1]
    alone = tail_trim(frame(bact, [BA] * 6), 0.05, "readn",
                      group_col="division", group_tails={VI: 0.0})
    withvi = tail_trim(frame(bact + [10000], [BA] * 6 + [VI]), 0.05, "readn",
                       group_col="division", group_tails={VI: 0.0})

    # 100,10,5,3,2,1 over 121: from the bottom 1+2+3 = 6 (4.96 %) fit in a
    # 5 % tail; adding 5 (9.1 %) does not.
    assert set(alone["taxid"]) == {"t0", "t1", "t2"}
    assert set(withvi["taxid"]) == {"t6"}, "only the exempt viral row survives"


def test_unnamed_groups_use_the_default_tail():
    df = frame([1000, 100, 10, 1], [BA, BA, "Eukaryota", VI])
    out = tail_trim(df, 0.01, "readn", group_col="division",
                    group_tails={VI: 0.0})
    # Eukaryota is not named, so it gets tail=0.01 like Bacteria.
    assert "Eukaryota" not in set(out["division"])
    assert VI in set(out["division"])


def test_all_groups_exempt_is_identity():
    df = frame([1000, 100, 10, 1], [BA, BA, VI, VI])
    out = tail_trim(df, 0.0, "readn", group_col="division",
                    group_tails={VI: 0.0})
    pd.testing.assert_frame_equal(out, df)


def test_group_col_absent_falls_back_to_the_flat_tail():
    # No division column: every row gets `tail`. The caller is responsible for
    # refusing to trim when it needs the groups -- see ReportStep.
    df = frame([1000, 100, 10, 1])
    out = tail_trim(df, 0.01, "readn", group_col="division",
                    group_tails={VI: 0.0})
    assert list(out["taxid"]) == ["t0", "t1"]


# ------------------------------------------------------------------ #
# readn_em basis
# ------------------------------------------------------------------ #

def test_default_count_col_is_resolved_not_assumed():
    # readn and readn_em rank the taxa in opposite orders; the trim must
    # follow readn_em, which resolve_count_col prefers.
    df = pd.DataFrame({
        "taxid": ["a", "b", "c"],
        "readn": [1, 500, 1000],
        "readn_em": [1000.0, 500.0, 1.0],
    })
    assert list(tail_trim(df, 0.01)["taxid"]) == ["a", "b"]
    assert list(tail_trim(df, 0.01, "readn")["taxid"]) == ["b", "c"]


def test_na_readn_em_falls_back_to_readn():
    df = pd.DataFrame({
        "taxid": ["a", "b", "c"],
        "readn": [1000, 500, 1],
        "readn_em": ["NA", "NA", "NA"],
    })
    assert list(tail_trim(df, 0.01)["taxid"]) == ["a", "b"]


@pytest.mark.parametrize("tail", [0.0, 0.001, 0.01, 0.05, 0.5])
def test_survivors_are_always_a_prefix_of_the_ranking(tail):
    df = frame([1000, 500, 100, 50, 10, 5, 1])
    out = tail_trim(df, tail, "readn")
    ranked = list(df.sort_values("readn", ascending=False)["taxid"])
    assert list(out["taxid"]) == ranked[:len(out)]
