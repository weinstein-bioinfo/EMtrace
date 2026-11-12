"""Automatic reference selection (assembly.ref_auto, analysis/ref_select.py).

Primary hits come from the assembly-DB search, secondary hits from annotation.
A secondary pick replaces a primary pick only if it aligns better, on the reads
both share, by more than secondary_margin; secondary picks the primary search
missed are added.
"""

import os
import time
from types import SimpleNamespace

import pandas as pd
import pytest
from pydantic import ValidationError

from emtrace.analysis.ref_select import candidates_are_current, select_references
from emtrace.core.config import EMtraceConfig, RefAutoConfig
from emtrace.core.exceptions import StepValidationError
from emtrace.steps.assembly import AssemblyStep


def _hits(acc, reads, species="Lassa virus", pident=95.0, alen=1000, qseqlen=1000,
          bitscore=None, db="minimap.assembly", taxid="11620", search="Minimap2"):
    """One hit of each read in *reads* to *acc*."""
    return pd.DataFrame([{
        "qseqid": r, "sseqid": acc, "species": species, "staxid": taxid, "db": db,
        "pident": pident, "alen": alen, "qseqlen": qseqlen, "qlen": alen,
        "qcov": alen / qseqlen, "bitscore": bitscore if bitscore is not None else alen * pident / 100,
        "division": "Viruses", "search": search,
    } for r in reads])


def _reads(prefix, n):
    return [f"{prefix}{i}" for i in range(n)]


def _cfg(**kw):
    return RefAutoConfig(enabled=True, **kw)


def _selected(table):
    return set(table.loc[table["selected"] == True, "accession"])  # noqa: E712


R = _reads("r", 30)


def test_secondary_within_margin_keeps_primary():
    picks, table = select_references(_hits("NC_1", R, pident=95.0),
                                     _hits("GB_1", R, pident=95.5), _cfg())
    assert [p.accession for p in picks] == ["NC_1"]
    assert "kept" in table.loc[table["accession"] == "NC_1", "reason"].iloc[0]


def test_secondary_beyond_margin_replaces_primary():
    picks, table = select_references(_hits("NC_1", R, pident=90.0),
                                     _hits("GB_1", R, pident=98.0), _cfg())
    assert [p.accession for p in picks] == ["GB_1"]
    reason = table.loc[table["accession"] == "GB_1", "reason"].iloc[0]
    assert reason.startswith("replaces primary NC_1")


def test_short_high_pident_alignment_does_not_win():
    # 100% identity over 30% of each read vs 90% over the whole read
    picks, _ = select_references(_hits("NC_1", R, pident=90.0, alen=1000),
                                 _hits("GB_1", R, pident=100.0, alen=300), _cfg())
    assert [p.accession for p in picks] == ["NC_1"]


def test_comparison_uses_shared_reads_only():
    # the secondary has more reads, but aligns worse on the reads both share
    prim = _hits("NC_1", R, pident=97.0)
    sec = pd.concat([_hits("GB_1", R, pident=93.0),
                     _hits("GB_1", _reads("x", 40), pident=99.0)])
    picks, _ = select_references(prim, sec, _cfg())
    assert "NC_1" in [p.accession for p in picks]
    assert "GB_1" not in [p.accession for p in picks]


def test_segments_with_disjoint_reads_are_both_picked():
    prim = pd.concat([_hits("NC_L", _reads("l", 40)), _hits("NC_S", _reads("s", 20))])
    picks, _ = select_references(prim, None, _cfg())
    assert {p.accession for p in picks} == {"NC_L", "NC_S"}


def test_strains_sharing_reads_pick_one():
    prim = pd.concat([_hits("NC_A", R, pident=95.0), _hits("NC_B", R, pident=90.0)])
    picks, table = select_references(prim, None, _cfg())
    assert [p.accession for p in picks] == ["NC_A"]
    assert "already covered" in table.loc[table["accession"] == "NC_B", "reason"].iloc[0]


def test_coinfection_picks_both_species():
    prim = pd.concat([_hits("NC_1", _reads("a", 30), species="Virus A", taxid="1"),
                      _hits("NC_2", _reads("b", 20), species="Virus B", taxid="2")])
    picks, _ = select_references(prim, None, _cfg())
    assert {p.accession for p in picks} == {"NC_1", "NC_2"}


def test_species_missed_by_primary_added_from_secondary_closest_genbank():
    prim = _hits("NC_1", _reads("a", 30), species="Virus A", taxid="1")
    sec = pd.concat([_hits("NC_1", _reads("a", 30), species="Virus A", taxid="1", db="ann"),
                     _hits("GB_far", _reads("b", 20), species="Virus B", taxid="2", pident=85.0, db="ann"),
                     _hits("GB_near", _reads("b", 20), species="Virus B", taxid="2", pident=97.0, db="ann")])
    picks, table = select_references(prim, sec, _cfg())
    assert {p.accession for p in picks} == {"NC_1", "GB_near"}
    added = table.loc[table["accession"] == "GB_near", "reason"].iloc[0]
    assert "not found by the primary search" in added


def test_without_primary_uses_annotation_hits_only():
    picks, table = select_references(None, _hits("GB_1", R, db="ann"), _cfg())
    assert [p.accession for p in picks] == ["GB_1"]
    assert (table["source"] == "secondary").all()


def test_species_below_min_reads_dropped():
    prim = pd.concat([_hits("NC_1", _reads("a", 30), species="Virus A", taxid="1"),
                      _hits("NC_2", _reads("b", 5), species="Virus B", taxid="2")])
    picks, table = select_references(prim, None, _cfg(min_reads=10))
    assert [p.accession for p in picks] == ["NC_1"]
    assert table["reason"].str.contains("below min_reads").any()


def test_both_sources_off_is_a_config_error():
    with pytest.raises(ValidationError):
        RefAutoConfig(use_assembly_db=False, use_annotation=False)


def test_candidates_reuse_requires_current_base_level_file(tmp_path):
    reads = tmp_path / "reads.fq"
    reads.write_text("@r\nA\n+\nI\n")
    cand = tmp_path / "assembly_candidates.tsv"
    _hits("NC_1", R).to_csv(cand, sep="\t", index=False)
    assert candidates_are_current(cand, reads)

    chain_only = _hits("NC_1", R).assign(qcov="")          # written without -c
    chain_only.to_csv(cand, sep="\t", index=False)
    assert not candidates_are_current(cand, reads)

    _hits("NC_1", R).to_csv(cand, sep="\t", index=False)
    old = time.time() - 3600
    os.utime(cand, (old, old))                             # older than the reads
    assert not candidates_are_current(cand, reads)


def test_explicit_refaccs_bypass_ref_auto():
    cfg = EMtraceConfig()
    cfg.assembly.refaccs = "NC_1"
    cfg.assembly.ref_auto.enabled = True
    step = AssemblyStep(cfg)

    def _fail(*a, **k):
        raise AssertionError("ref_auto must not run when refaccs are set")

    step._select_refs = _fail
    step._resolve_assembly_reads = lambda sample: (None, None)   # stop right after the check
    sample = SimpleNamespace()
    assert step.run(sample) is sample


# --- protein-search hits are not nucleotide references ---------------------

def test_protein_search_hits_are_not_candidates():
    primary = _hits("NC_1", R, pident=75.0, alen=330)
    secondary = _hits("UniRef100_X", R, pident=95.0, alen=1000,
                      db="diamond.uniref.vi", search="Diamond blastx")
    _, table = select_references(primary, secondary, _cfg())
    assert _selected(table) == {"NC_1"}
    prot = table[table["reason"].str.contains("protein search", na=False)]
    assert len(prot) == 1
    row = prot.iloc[0]
    assert row["db"] == "diamond.uniref.vi" and row["source"] == "secondary"
    assert row["reads"] == len(R) and not row["selected"]


def test_chained_round_keeps_nucleotide_rows():
    secondary = pd.concat([
        _hits("NC_2", R, pident=90.0, db="minimap.refseq.abv.1g"),
        _hits("UniRef100_Y", R, pident=99.0, db="diamond.uniref.vi", search="Diamond blastx"),
    ])
    picks, table = select_references(None, secondary, _cfg())
    assert [p.accession for p in picks] == ["NC_2"]
    assert "UniRef100_Y" not in set(table["accession"])


def test_hits_without_search_column_still_used():
    secondary = _hits("NC_3", R).drop(columns="search")
    picks, _ = select_references(None, secondary, _cfg())
    assert [p.accession for p in picks] == ["NC_3"]


def test_empty_combined_reference_raises(tmp_path):
    step = AssemblyStep(EMtraceConfig())
    step._fetch_accessions = lambda refaccs, asm_dir, sample: None
    with pytest.raises(StepValidationError, match="UniRef100_X"):
        step._build_combined_ref([], ["UniRef100_X"], tmp_path, SimpleNamespace())
