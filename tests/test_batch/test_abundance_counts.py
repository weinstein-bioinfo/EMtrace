"""batch-abundance counts use EM-adjusted reads (readn_em), raw readn as fallback."""

from pathlib import Path

from emtrace.batch.abundance import build_abundance

HEADER = "taxid\treadn\tbitscore\tprob\treadn_em\tspecies\tgenus\tfamily\torder\tclass\tdivision\torganism\n"
HEADER_NO_EM = "taxid\treadn\tbitscore\tprob\tspecies\tgenus\tfamily\torder\tclass\tdivision\torganism\n"


def _write(resdir: Path, sample: str, text: str) -> None:
    d = resdir / sample / "results" / "reports"
    d.mkdir(parents=True)
    (d / "abund_table.tsv").write_text(text)


def test_counts_use_readn_em_with_readn_fallback(tmp_path: Path):
    _write(tmp_path, "s1", HEADER
           + "11\t10\t1.0\t0.7\t12.6\tVirus alpha\tAlphavirus\tF\tO\tC\tViruses\tVirus alpha\n"
           + "12\t5\t1.0\t0.3\t\tVirus beta\tBetavirus\tF\tO\tC\tViruses\tVirus beta\n")
    _write(tmp_path, "s2", HEADER_NO_EM
           + "11\t7\t1.0\t1.0\tVirus alpha\tAlphavirus\tF\tO\tC\tViruses\tVirus alpha\n")

    counts = build_abundance(tmp_path)["viruses"]["counts"].set_index("species")

    assert counts.loc["Virus alpha", "s1"] == 13   # readn_em 12.6, rounded
    assert counts.loc["Virus beta", "s1"] == 5     # readn_em missing: raw readn
    assert counts.loc["Virus alpha", "s2"] == 7    # no readn_em column: raw readn


def test_batch_taxa_readn_uses_readn_em(tmp_path: Path):
    from emtrace.batch.taxa import build_taxa_summary

    _write(tmp_path, "s1", HEADER
           + "11\t10\t1.0\t0.7\t12.6\tVirus alpha\tAlphavirus\tF\tO\tC\tViruses\tVirus alpha\n"
           + "12\t5\t1.0\t0.3\t\tVirus beta\tBetavirus\tF\tO\tC\tViruses\tVirus beta\n")
    _write(tmp_path, "s2", HEADER_NO_EM
           + "11\t7\t1.0\t1.0\tVirus alpha\tAlphavirus\tF\tO\tC\tViruses\tVirus alpha\n")

    df = build_taxa_summary(tmp_path, target="vi")
    got = {(r.Sample, str(r.taxid)): int(r.readn) for r in df.itertuples()}
    assert got == {("s1", "11"): 13, ("s1", "12"): 5, ("s2", "11"): 7}
