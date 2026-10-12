"""`report.tail` — the cumulative-share tail trim on the abundance table.

The trim is computed on the abundance table and then propagated, so the
annotation table, the sorted reads and every plot describe the same taxa. Two
pieces carry that: `_apply_tail_trim` decides what survives, and
`_restrict_annot_to_taxa` pushes the decision into the annotation table.

The share itself is `emtrace.analysis.tailtrim`, covered by
tests/test_analysis/test_tailtrim.py; what is pinned here is the policy —
off by default, viruses protected, and a refusal to trim blind.
"""

import logging

import pandas as pd
import pytest

from emtrace.core.config import EMtraceConfig, ReportConfig
from emtrace.steps.report import ReportStep

VI, BA = "Viruses", "Bacteria"


@pytest.fixture
def step():
    return ReportStep(config=EMtraceConfig())


def abund(counts, divisions=None, readn_em=None):
    """Abundance table: one row per taxon, taxids t0..tN."""
    data = {"taxid": [f"t{i}" for i in range(len(counts))], "readn": counts}
    data["division"] = divisions if divisions is not None else [BA] * len(counts)
    if readn_em is not None:
        data["readn_em"] = readn_em
    return pd.DataFrame(data)


# A long tail: t0 holds 90 % and the rest trail off.
LONG_TAIL = abund([900, 50, 30, 10, 5, 3, 2])


# ------------------------------------------------------------------ #
# Off by default
# ------------------------------------------------------------------ #

def test_default_config_does_not_trim(step):
    cfg = ReportConfig()
    assert cfg.tail == 0.0
    pd.testing.assert_frame_equal(step._apply_tail_trim(LONG_TAIL, cfg), LONG_TAIL)


def test_shipped_default_yaml_does_not_trim():
    # The model default and default.yaml must agree, or a --config run would
    # behave differently from a bare one.
    cfg = EMtraceConfig.create_default()
    assert cfg.report.tail == 0.0
    assert cfg.report.tail_viruses is False


def test_explicit_zero_is_identity(step):
    cfg = ReportConfig(tail=0.0, tail_viruses=True)
    pd.testing.assert_frame_equal(step._apply_tail_trim(LONG_TAIL, cfg), LONG_TAIL)


def test_empty_table_is_identity(step):
    empty = abund([])
    cfg = ReportConfig(tail=0.01)
    pd.testing.assert_frame_equal(step._apply_tail_trim(empty, cfg), empty)


# ------------------------------------------------------------------ #
# Trimming, and virus protection
# ------------------------------------------------------------------ #

def test_trim_drops_the_smallest_taxa(step):
    out = step._apply_tail_trim(LONG_TAIL, ReportConfig(tail=0.01))
    assert len(out) < len(LONG_TAIL)
    # Survivors are a prefix of the abundance ranking.
    assert list(out["taxid"]) == list(LONG_TAIL["taxid"])[:len(out)]


def test_viruses_are_protected_by_default(step):
    # The two viral taxa are the smallest in the library, so an unprotected
    # 1 % trim would take them first.
    df = abund([900, 50, 30, 10, 2, 1], [BA, BA, BA, BA, VI, VI])
    out = step._apply_tail_trim(df, ReportConfig(tail=0.01))
    assert set(out[out["division"] == VI]["taxid"]) == {"t4", "t5"}


def test_tail_viruses_true_trims_them(step):
    df = abund([900, 50, 30, 10, 2, 1], [BA, BA, BA, BA, VI, VI])
    out = step._apply_tail_trim(df, ReportConfig(tail=0.01, tail_viruses=True))
    assert VI not in set(out["division"])


def test_protection_does_not_spare_the_background(step):
    # Viral reads still count toward the share, so protecting viruses makes
    # the trim harsher on everything else, not gentler. Documented because it
    # surprises people.
    df = abund([9000, 50, 30, 10], [VI, BA, BA, BA])
    out = step._apply_tail_trim(df, ReportConfig(tail=0.01))
    assert set(out["division"]) == {VI}


def test_readn_em_is_preferred_when_present(step):
    # readn and readn_em rank these in opposite orders; the trim must follow
    # readn_em.
    df = abund([1, 50, 900], readn_em=[900.0, 50.0, 1.0])
    out = step._apply_tail_trim(df, ReportConfig(tail=0.01, tail_viruses=True))
    assert list(out["taxid"]) == ["t0"]


def test_na_readn_em_falls_back_to_readn(step):
    df = abund([900, 50, 1], readn_em=["NA", "NA", "NA"])
    out = step._apply_tail_trim(df, ReportConfig(tail=0.01, tail_viruses=True))
    assert list(out["taxid"]) == ["t0"]


def test_row_values_are_not_rewritten(step):
    # preserve_order=True: the surviving rows must be byte-for-byte what came
    # in, so the written TSV is the input minus rows.
    df = abund([900, 50, 1], readn_em=["900.0", "50.0", "1.0"])
    out = step._apply_tail_trim(df, ReportConfig(tail=0.01, tail_viruses=True))
    assert out["readn_em"].tolist() == ["900.0"]


# ------------------------------------------------------------------ #
# Refuse to trim blind
# ------------------------------------------------------------------ #

def test_no_division_column_refuses_to_trim(step, caplog):
    df = LONG_TAIL.drop(columns=["division"])
    with caplog.at_level(logging.WARNING):
        out = step._apply_tail_trim(df, ReportConfig(tail=0.01))
    pd.testing.assert_frame_equal(out, df)
    assert "NOT trimming" in caplog.text


def test_empty_division_column_refuses_to_trim(step, caplog):
    # What a run whose taxonomy linkage failed looks like.
    df = LONG_TAIL.assign(division="")
    with caplog.at_level(logging.WARNING):
        out = step._apply_tail_trim(df, ReportConfig(tail=0.01))
    pd.testing.assert_frame_equal(out, df)
    assert "NOT trimming" in caplog.text


def test_refusal_holds_even_when_viruses_are_not_protected(step, caplog):
    # tail_viruses=true would not need `division`, but trimming a table with
    # no taxonomy at all is not something to do silently.
    df = LONG_TAIL.drop(columns=["division"])
    with caplog.at_level(logging.WARNING):
        out = step._apply_tail_trim(df, ReportConfig(tail=0.01, tail_viruses=True))
    pd.testing.assert_frame_equal(out, df)


def test_dropping_rows_is_logged_at_warning(step, caplog):
    with caplog.at_level(logging.WARNING):
        step._apply_tail_trim(LONG_TAIL, ReportConfig(tail=0.01))
    assert "report.tail" in caplog.text
    assert "viral taxa exempt" in caplog.text


# ------------------------------------------------------------------ #
# Propagation into the annotation table
# ------------------------------------------------------------------ #

def annot(staxids):
    return pd.DataFrame({
        "qseqid": [f"r{i}" for i in range(len(staxids))],
        "staxid": staxids,
        "bitscore": [100] * len(staxids),
    })


def test_annot_is_cut_to_the_surviving_taxa():
    kept = ReportStep._restrict_annot_to_taxa(
        annot(["t0", "t1", "t2", "t2"]), abund([9, 9]),  # taxids t0, t1
    )
    assert set(kept["staxid"]) == {"t0", "t1"}
    assert len(kept) == 2


def test_propagation_compares_as_strings():
    # The abundance table carries taxid as text; annotation tables have been
    # seen with it as int. A dtype mismatch must not silently drop everything.
    a = pd.DataFrame({"qseqid": ["r0", "r1"], "staxid": [1234, 5678]})
    b = pd.DataFrame({"taxid": ["1234"], "readn": [10]})
    kept = ReportStep._restrict_annot_to_taxa(a, b)
    assert kept["staxid"].tolist() == [1234]


def test_propagation_is_a_noop_without_the_join_columns():
    a = annot(["t0", "t1"])
    pd.testing.assert_frame_equal(
        ReportStep._restrict_annot_to_taxa(a, pd.DataFrame({"readn": [1]})), a,
    )


def test_trim_and_propagation_agree_end_to_end(step):
    """The invariant the reordering exists to guarantee."""
    ab = abund([900, 50, 30, 10, 5, 3, 2])
    an = annot([f"t{i}" for i in range(7)] + ["t6", "t5"])

    trimmed = step._apply_tail_trim(ab, ReportConfig(tail=0.01, tail_viruses=True))
    kept = ReportStep._restrict_annot_to_taxa(an, trimmed)

    assert set(kept["staxid"]) <= set(trimmed["taxid"].astype(str))
    assert set(kept["staxid"]) == set(trimmed["taxid"].astype(str))
