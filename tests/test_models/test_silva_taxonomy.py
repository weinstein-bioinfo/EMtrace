"""SILVA lineage rank assignment.

Regression cover for a silent, severe bug: SILVA 144 inserted a ``kingdom``
level into every bacterial lineage, so the positional parser shifted each rank
down one and wrote the *family* name into the *genus* column. A database built
that way reported zero rows for genus *Bacillus* and half a million for genus
*Enterobacteriaceae*, and nothing failed — genus is the headline rank of the
16S benchmark, so every SILVA number would have been wrong but plausible.
"""

from pathlib import Path

import pytest

from emtrace.models.taxonomy import SILVA_TAXMAP_RANKS, SILVATaxonomy


# Real SILVA 144 lineages, seven levels: domain;kingdom;phylum;class;order;family;genus
ECOLI_144 = (
    "Bacteria;Pseudomonadati;Pseudomonadota;Gammaproteobacteria;"
    "Enterobacterales;Enterobacteriaceae;Escherichia-Shigella;"
)
# The same organism in SILVA 138.2, six levels — no kingdom.
ECOLI_138 = (
    "Bacteria;Pseudomonadota;Gammaproteobacteria;"
    "Enterobacterales;Enterobacteriaceae;Escherichia-Shigella;"
)


def _rank_of(path: str, ranks: list[str]) -> dict[str, str]:
    """Build a prefix -> rank map the way tax_slv_ssu_*.txt declares one."""
    parts = [p for p in path.rstrip(";").split(";") if p]
    out, prefix = {}, ""
    for part, rank in zip(parts, ranks):
        prefix = f"{prefix}{part};"
        out[prefix] = rank
    return out


RANKS_144 = ["domain", "kingdom", "phylum", "class", "order", "family", "genus"]
RANKS_138 = ["domain", "phylum", "class", "order", "family", "genus"]


@pytest.fixture
def tax():
    return SILVATaxonomy(Path("unused.db"))


def test_ranked_parse_skips_the_kingdom_level(tax):
    """SILVA 144: each element lands in the rank SILVA declares, kingdom dropped."""
    got = tax._parse_taxpath_ranked(ECOLI_144, _rank_of(ECOLI_144, RANKS_144))
    assert got == {
        "domain": "Bacteria",
        "phylum": "Pseudomonadota",
        "class": "Gammaproteobacteria",
        "order": "Enterobacterales",
        "family": "Enterobacteriaceae",
        "genus": "Escherichia-Shigella",
    }
    # kingdom is not a column in this schema and must not displace anything.
    assert "Pseudomonadati" not in got.values()


def test_positional_parse_is_wrong_for_silva_144(tax):
    """The documented failure mode, pinned so the fallback's limits stay visible."""
    got = tax._parse_taxpath(ECOLI_144)
    assert got["genus"] == "Enterobacteriaceae"   # the family name
    assert got["phylum"] == "Pseudomonadati"      # the kingdom name


def test_both_parsers_agree_on_silva_138(tax):
    """Pre-144 lineages have no kingdom, so positional and ranked coincide."""
    ranked = tax._parse_taxpath_ranked(ECOLI_138, _rank_of(ECOLI_138, RANKS_138))
    assert ranked == tax._parse_taxpath(ECOLI_138)
    assert ranked["genus"] == "Escherichia-Shigella"


def test_truncated_lineage_leaves_lower_ranks_unset(tax):
    """An unclassified path must not have its tail guessed into a genus."""
    path = "Bacteria;environmental samples;"
    rank_of = {"Bacteria;": "domain"}  # SILVA declares no rank for the second element
    got = tax._parse_taxpath_ranked(path, rank_of)
    assert got == {"domain": "Bacteria"}
    assert "genus" not in got


def test_every_assigned_rank_is_in_the_schema(tax):
    """Ranks SILVA declares but this schema lacks are skipped, never stored."""
    got = tax._parse_taxpath_ranked(ECOLI_144, _rank_of(ECOLI_144, RANKS_144))
    assert set(got).issubset(set(SILVA_TAXMAP_RANKS))


# --- Species keys -----------------------------------------------------------
# A SILVA taxid identifies a lineage, not an organism: the Escherichia-Shigella
# taxid alone carries 3,851 organism names in SILVA 144. The report used to label
# each taxid with one arbitrary organism_name ("Solea senegalensis" for E. coli).

from emtrace.models.taxonomy import (  # noqa: E402
    assign_silva_species_taxids,
    clean_silva_species,
)


@pytest.mark.parametrize("organism, genus, expected", [
    ("Escherichia coli O157:H7 str. Sakai", "Escherichia-Shigella", "Escherichia coli"),
    ("Shigella flexneri 2a", "Escherichia-Shigella", "Shigella flexneri"),
    ("Solea senegalensis (Senegalese sole)", "Escherichia-Shigella", ""),
    ("Rhizoctonia solani", "Pseudomonas", ""),
    ("Bacillus sp. KL-2", "Bacillus", ""),
    ("uncultured Bacillus sp.", "Bacillus", ""),
    ("gut metagenome", "Bacillus", ""),
    ("Lactobacillus fermentum", "Limosilactobacillus", ""),
    ("Candidatus Pelagibacter ubique HTCC1062", "Candidatus Pelagibacter", "Candidatus Pelagibacter ubique"),
    ("[Eubacterium] rectale", "[Eubacterium] eligens group", "Eubacterium rectale"),
    ("Salmonella phage", "Salmonella", ""),
    ("", "Salmonella", ""),
])
def test_clean_silva_species(organism, genus, expected):
    assert clean_silva_species(organism, genus) == expected


def test_species_keys_are_deterministic_and_keep_lineage_taxid():
    pairs = {("100", "Shigella flexneri"), ("100", "Escherichia coli"), ("100", ""), ("200", "")}
    got = assign_silva_species_taxids(pairs)
    assert got == {
        ("100", "Escherichia coli"): "100.1",
        ("100", "Shigella flexneri"): "100.2",
        ("100", ""): "100",
        ("200", ""): "200",
    }


def test_build_and_link_species(tmp_path):
    """End to end: build from a taxmap, then staxid and lineage resolve to species."""
    pd = pytest.importorskip("pandas")
    taxmap = tmp_path / "taxmap.txt"
    taxmap.write_text(
        "primaryAccession\tstart\tstop\tpath\torganism_name\ttaxid\n"
        f"A1\t1\t1500\t{ECOLI_144}\tEscherichia coli K-12\t100\n"
        f"A2\t1\t1500\t{ECOLI_144}\tSolea senegalensis (Senegalese sole)\t100\n"
        f"A3\t1\t1500\t{ECOLI_144}\tShigella sonnei\t100\n"
    )
    tax_slv = tmp_path / "tax_slv.txt"
    prefix, lines = "", []
    for part, rank in zip(ECOLI_144.rstrip(";").split(";"), RANKS_144):
        prefix += f"{part};"
        lines.append(f"{prefix}\t1\t{rank}\t\t144")
    tax_slv.write_text("\n".join(lines) + "\n")

    silva = SILVATaxonomy(tmp_path / "taxonomy.db")
    silva.build_database(taxmap, tax_path=tax_slv)
    assert silva.has_species()

    df = pd.DataFrame({"sseqid": ["A1.1.1500", "A2.1.1500", "A3.1.1500"]})
    silva.add_lineage_to_df(df)
    assert list(df["staxid"]) == ["100.1", "100", "100.2"]
    assert list(df["species"]) == ["Escherichia coli", "", "Shigella sonnei"]

    lineage = silva.lineage_by_taxid()
    assert lineage["100.1"]["species"] == "Escherichia coli"
    assert lineage["100.1"]["genus"] == "Escherichia-Shigella"
    assert lineage["100"]["species"] == ""   # genus-only row, never a fish
