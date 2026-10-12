"""Tests for the no-species filter in analysis/taxonomy_utils.py.

Regression cover for the benchmarking bug where ``em.exclude_no_species``
filtered the ``sname`` column with ``r"\\w+ sp\\."``.  For minimap2 hits
``sname`` is the subject accession (``paf_to_annotation`` sets
``sname = sseqid``), so the filter matched nothing on every minimap2-based
run and ``silva.silva.em.nosp`` came out a duplicate of ``silva.silva``.

The filter now tests the ``species`` rank of the assigned taxpath.
"""

import pandas as pd
import pytest

from emtrace.analysis.em import em_loop
from emtrace.analysis.taxonomy_utils import (
    filter_no_species,
    is_no_species,
)


# ------------------------------------------------------------------
# Predicate
# ------------------------------------------------------------------

@pytest.mark.parametrize(
    "species",
    [
        "Escherichia coli",
        "Staphylococcus aureus",
        # Viral names are non-Linnaean and must survive — EMtrace is
        # primarily a virus pipeline.
        "Tulip top breaking virus",
        "Bacteroides phage PhiCrAssBcn18",
        "Invertebrate iridescent virus 1",
        "Malignant rabbit fibroma virus",
    ],
)
def test_resolved_species_kept(species):
    assert not is_no_species(pd.Series([species])).iloc[0]


@pytest.mark.parametrize(
    "species",
    [
        # Non-binomial NCBI names — 49.4% of the installed taxonomy DB
        "Vitreoscilla sp.",
        "Caulobacter sp.",
        "Methylobacterium sp. E-123",
        "Eulophidae sp. BIOUG26785-A12",
        "Rhizobium cf. leguminosarum",
        "Drosophila aff. melanogaster",
        # Unresolved markers
        "uncultured bacterium",
        "unidentified adenovirus",
        "environmental samples",
        "unclassified Bacteria",
    ],
)
def test_non_binomial_dropped(species):
    assert is_no_species(pd.Series([species])).iloc[0]


@pytest.mark.parametrize(
    "species", ["", "   ", None, "NA", "nan", "none", "-", "null"],
)
def test_empty_and_placeholders_dropped(species):
    """species is NULL for every taxon above species rank (7.6% of NCBI)."""
    assert is_no_species(pd.Series([species])).iloc[0]


def test_species_pattern_needs_word_boundary():
    """'sp.' must not match inside a longer word."""
    assert not is_no_species(pd.Series(["Genus crisp.ata"])).iloc[0]


def test_bacterium_suffix_is_not_a_marker():
    """A broader pattern would discard 111k validly-named taxa."""
    assert not is_no_species(
        pd.Series(["Candidatus Carsonella ruddii"])
    ).iloc[0]


# ------------------------------------------------------------------
# filter_no_species
# ------------------------------------------------------------------

def test_filter_drops_only_unresolved_rows():
    df = pd.DataFrame({
        "qseqid": ["r1", "r2", "r3", "r4"],
        "species": ["Escherichia coli", "Vitreoscilla sp.", "", "Tulip top breaking virus"],
    })
    out = filter_no_species(df)
    assert list(out["qseqid"]) == ["r1", "r4"]


def test_filter_missing_column_is_a_noop(caplog):
    """An old annot table must degrade to 'keep everything', not drop all."""
    df = pd.DataFrame({"qseqid": ["r1", "r2"], "sname": ["ACC1", "ACC2"]})
    out = filter_no_species(df)
    assert len(out) == 2
    assert "not found" in caplog.text.lower()


def test_filter_does_not_mutate_input():
    df = pd.DataFrame({"species": ["Escherichia coli", "Vitreoscilla sp."]})
    before = df.copy()
    filter_no_species(df)
    assert df.equals(before)


# ------------------------------------------------------------------
# The regression itself
# ------------------------------------------------------------------

def _minimap_frame():
    """A minimap2-style annotation frame.

    ``sname`` holds SILVA accessions, exactly as ``paf_to_annotation``
    writes them — this is what the old filter looked at.
    """
    return pd.DataFrame({
        "qseqid": ["r1", "r2", "r3"],
        "staxid": ["101", "102", "103"],
        "bitscore": [100.0, 90.0, 80.0],
        "sname": [
            "HG345195.864.2657",
            "AB000001.1.1450",
            "AY123456.1.1500",
        ],
        "species": ["Escherichia coli", "Vitreoscilla sp.", ""],
    })


def test_minimap_accessions_do_not_defeat_the_filter():
    """The exact regression: accession-valued sname, binomial species.

    Fails on the old sname-based implementation, which dropped nothing.
    """
    out = filter_no_species(_minimap_frame())
    assert list(out["qseqid"]) == ["r1"]


def test_em_loop_excludes_by_species_not_sname():
    df = _minimap_frame()
    kept = em_loop(df, exclude_no_species=True)
    allrows = em_loop(df, exclude_no_species=False)

    # taxid 101 is the only row with a resolved binomial species
    assert set(kept.F) == {"101"}
    # The filter must actually change the outcome — the bug was that it
    # produced an identical result either way.
    assert set(allrows.F) != set(kept.F)


def test_em_loop_leaves_caller_frame_intact():
    """EM filtering must not reach the annotation tables on disk."""
    df = _minimap_frame()
    before = df.copy()
    em_loop(df, exclude_no_species=True)
    assert df.equals(before)


def test_em_loop_without_species_column_keeps_all_taxa():
    df = _minimap_frame().drop(columns=["species"])
    res = em_loop(df, exclude_no_species=True)
    assert set(res.F) == {"101", "102", "103"}
