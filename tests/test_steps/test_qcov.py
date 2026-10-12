"""Query coverage is filtered only at reporting (report.min_qcov).

Search tools never filtered on coverage, so the per-search qcov keys are gone;
hits with unknown coverage (minimap2 without -c) skip the filter.
"""

from pathlib import Path

import pandas as pd
from click.testing import CliRunner

from emtrace.__main__ import cli
from emtrace.core.config import EMtraceConfig
from emtrace.steps.report import ReportStep


def _annot(qcovs):
    return pd.DataFrame({
        "search": "Minimap2", "db": "db", "orf": "NA", "bitscore": 100,
        "qseqid": [f"r{i}" for i in range(len(qcovs))],
        "staxid": [str(10 + i) for i in range(len(qcovs))],
        "qcov": qcovs,
    })


def test_search_configs_have_no_qcov():
    ann = EMtraceConfig().annotation
    for tool in (ann.blastn, ann.minimap, ann.diamond_blastp):
        assert not hasattr(tool, "qcov")


def test_legacy_search_qcov_keys_warn(tmp_dir: Path, caplog):
    p = tmp_dir / "cfg.yaml"
    p.write_text("annotation:\n  minimap:\n    qcov: 0.4\n    s: 150\n")
    cfg = EMtraceConfig.from_yaml(p)
    assert cfg.annotation.minimap.s == 150
    assert "annotation.minimap.qcov" in caplog.text


def test_min_qcov_help_points_at_report():
    out = CliRunner().invoke(cli, ["run", "--help"]).output
    assert "report.min_qcov" in out


def test_unknown_qcov_is_not_filtered():
    # r0 known and too low, r1 known and passing, r2/r3 unknown (NA, empty).
    annot = _annot([0.1, 0.9, pd.NA, ""])
    rep = ReportStep(config=EMtraceConfig())
    table = rep._build_annotation_table(annot, 0.4, 0)
    abund = rep._build_abundance_table(annot, 0.4, 0)
    assert set(table["qseqid"]) == {"r1", "r2", "r3"}
    assert set(abund["taxid"]) == {"11", "12", "13"}
