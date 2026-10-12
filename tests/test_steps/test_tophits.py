"""Top-hit filtering is opt-in (annotation.filter_tophits / --tophits).

All hits are kept by default so EM can weigh alternatives; read partitions and
the annotation/abundance tables select top hits themselves, so they come out
the same whether or not sub-best hits were kept at annotation.
"""

from pathlib import Path
from types import SimpleNamespace

import pandas as pd
from click.testing import CliRunner

from emtrace.__main__ import cli
from emtrace.analysis.em import em_loop, log_hit_structure
from emtrace.core.config import EMtraceConfig
from emtrace.io.tsv import top_hits_per_search
from emtrace.steps.annotate import AnnotateStep
from emtrace.steps.report import ReportStep


def _hits(rows):
    """Annotation rows: (db, qseqid, staxid, bitscore[, qcov[, search]])."""
    out = []
    for r in rows:
        db, q, t, b = r[:4]
        qcov = r[4] if len(r) > 4 else 1.0
        search = r[5] if len(r) > 5 else "Minimap2"
        out.append({"search": search, "db": db, "qseqid": q, "orf": "NA",
                    "staxid": t, "bitscore": b, "qcov": qcov})
    return pd.DataFrame(out)


# All hits for two reads: r1 ties 1|2 at the top with a sub-best 3; r2 has a
# top hit 4 and sub-best 5.
ALL = _hits([
    ("db", "r1", "1", 100), ("db", "r1", "2", 100), ("db", "r1", "3", 90),
    ("db", "r2", "4", 80), ("db", "r2", "5", 70),
])
TOP = ALL[ALL["staxid"].isin(["1", "2", "4"])].reset_index(drop=True)


def _report():
    return ReportStep(config=EMtraceConfig())


def _em_file(path: Path, probs: dict, readn_em: bool = True) -> Path:
    df = pd.DataFrame({"taxid": list(probs), "prob": list(probs.values())})
    if readn_em:
        df["readn_em"] = df["prob"] * 2
    df.to_csv(path, sep="\t", index=False)
    return path


class TestTopHitsPerSearch:
    def test_keeps_ties_drops_sub_best(self):
        assert set(top_hits_per_search(ALL)["staxid"]) == {"1", "2", "4"}

    def test_per_search(self):
        df = _hits([("a", "r1", "1", 100), ("b", "r1", "2", 5)])
        assert len(top_hits_per_search(df)) == 2

    def test_idempotent(self):
        once = top_hits_per_search(ALL)
        assert top_hits_per_search(once).equals(once)


class TestConfigAndCli:
    def test_default_keeps_all_hits(self):
        assert EMtraceConfig().annotation.filter_tophits is False
        assert EMtraceConfig.create_default().annotation.filter_tophits is False

    def test_yaml_opt_in(self, tmp_dir: Path):
        p = tmp_dir / "cfg.yaml"
        p.write_text("annotation:\n  filter_tophits: true\n")
        assert EMtraceConfig.from_yaml(p).annotation.filter_tophits is True

    def test_cli_flag(self):
        out = CliRunner().invoke(cli, ["run", "--help"]).output
        assert "--tophits" in out
        assert "--allhits" not in out


class TestPartition:
    def test_sub_best_hit_does_not_move_read(self, tmp_dir: Path, sample_fastq: Path):
        # read1's top hit is bacterial; a weaker viral hit must not put it in vi.
        annot = _hits([("db", "read1", "1", 100), ("db", "read1", "2", 50)])
        annot["division"] = ["Bacteria", "Viruses"]
        annot["bphage"] = "no"
        files = {}
        sample = SimpleNamespace(reads_dir=tmp_dir, set_file=files.__setitem__)
        AnnotateStep(EMtraceConfig())._partition_reads(sample_fastq, annot, "ann1", sample)
        assert "read1" not in (tmp_dir / "reads.ann1.vi.fq").read_text()
        assert "read1" in (tmp_dir / "reads.ann1.ab.fq").read_text()


class TestReportTables:
    def test_annotation_table_same_for_top_and_all_hits(self):
        rep = _report()
        a = rep._build_annotation_table(ALL, 0, 0)
        b = rep._build_annotation_table(TOP, 0, 0)
        assert a[["qseqid", "staxid"]].values.tolist() == b[["qseqid", "staxid"]].values.tolist()
        assert set(a["staxid"]) == {"1", "2", "4"}

    def test_abundance_same_for_top_and_all_hits(self, tmp_dir: Path):
        em = _em_file(tmp_dir / "em.tsv", {"1": 0.2, "2": 0.5, "4": 0.3})
        rep = _report()
        a = rep._build_abundance_table(ALL, 0, 0, em_prob_path=em)
        b = rep._build_abundance_table(TOP, 0, 0, em_prob_path=em)
        cols = ["taxid", "readn", "prob", "readn_em"]
        assert a[cols].sort_values("taxid").values.tolist() == b[cols].sort_values("taxid").values.tolist()

    def test_top_hit_failing_qcov_is_not_promoted(self):
        df = _hits([("db", "r1", "1", 100, 0.1), ("db", "r1", "2", 90, 0.9)])
        rep = _report()
        assert rep._build_annotation_table(df, 0.5, 0).empty
        assert rep._build_abundance_table(df, 0.5, 0).empty

    def test_tie_goes_to_most_abundant_taxon(self, tmp_dir: Path):
        em = _em_file(tmp_dir / "em.tsv", {"1": 0.2, "2": 0.5, "4": 0.3})
        abund = _report()._build_abundance_table(ALL, 0, 0, em_prob_path=em)
        readn = dict(zip(abund["taxid"], abund["readn"]))
        assert readn == {"2": 1, "4": 1}
        assert set(abund.columns) >= {"prob", "readn_em"}

    def test_tie_without_em_goes_to_lowest_taxid(self):
        abund = _report()._build_abundance_table(ALL.sample(frac=1, random_state=3), 0, 0)
        assert dict(zip(abund["taxid"], abund["readn"])) == {"1": 1, "4": 1}
        assert (abund["readn_em"] == "NA").all()

    def test_prot_tied_read_counted_once(self):
        df = _hits([("dp", "r1", "7", 50, 1.0, "Diamond blastx"),
                    ("dp", "r1", "8", 50, 1.0, "Diamond blastx")])
        abund = _report()._build_abundance_table(df, 0, 0)
        assert abund["readn"].sum() == 1


class TestEmHitStructure:
    def test_counts_and_readn_em(self):
        df = _hits([
            ("db", "q1", "A", 100),                          # single
            ("db", "q2", "A", 100), ("db", "q2", "B", 100),  # tied
            ("db", "q3", "A", 100), ("db", "q3", "B", 60),   # ranked
        ])
        res = em_loop(df, temp=10.0)
        assert (res.n_queries, res.n_single, res.n_tied, res.n_ranked) == (3, 1, 1, 1)
        assert abs(sum(res.readn_em().values()) - 3) < 1e-9

    def test_warns_when_only_ties(self, caplog):
        df = _hits([("db", "q1", "A", 100), ("db", "q1", "B", 100), ("db", "q2", "A", 50)])
        import logging
        log = logging.getLogger("test.em")
        assert log_hit_structure(em_loop(df), log, filter_tophits=True)
        assert "--tophits" in caplog.text
        caplog.clear()
        assert log_hit_structure(em_loop(df), log, filter_tophits=False)
        assert "No read has ranked alternative hits" in caplog.text

    def test_no_warning_with_ranked_hits(self, caplog):
        import logging
        res = em_loop(_hits([("db", "q1", "A", 100), ("db", "q1", "B", 60)]))
        assert not log_hit_structure(res, logging.getLogger("test.em"), filter_tophits=False)
        assert "WARNING" not in caplog.text
