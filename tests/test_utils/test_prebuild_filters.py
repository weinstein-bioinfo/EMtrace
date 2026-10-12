"""Tests for host-filter index pre-building (databases --build-index)."""

from __future__ import annotations

from pathlib import Path

import pytest

from emtrace.core.config import DatabaseEntry, EMtraceConfig
from emtrace.utils import db_download
from emtrace.utils.db_download import _host_mmi_path, prebuild_host_index
from emtrace.wrappers.minimap2 import Minimap2


@pytest.mark.parametrize(
    "k, w, expected_suffix",
    [
        (0, 0, ".mmi"),
        (19, 0, ".k19w0.mmi"),
        (0, 11, ".k0w11.mmi"),
        (19, 11, ".k19w11.mmi"),
    ],
)
def test_host_mmi_path_matches_bgfilter_suffix(k, w, expected_suffix):
    """The prebuilt index path must equal the suffix BgFilterStep computes.

    This mirrors the inline logic in ``BgFilterStep._align_minimap2`` — if the
    two ever diverge the run would never find the prebuilt index.
    """
    db = Path("/data/hostgenomes/host.fna.gz")
    assert _host_mmi_path(db, k, w) == Path(str(db) + expected_suffix)


def test_prebuild_host_index_builds_at_runtime_path(tmp_path, monkeypatch):
    """prebuild builds the .mmi at the exact path the filter step will read."""
    fasta = tmp_path / "test.fna.gz"
    fasta.write_text(">x\nACGT\n")  # present -> download is skipped

    cfg = EMtraceConfig.create_default()
    cfg.host_databases["TestHost"] = DatabaseEntry(
        db=str(fasta), name="TestHost", search="Minimap2",
    )

    install_called = []
    monkeypatch.setattr(
        db_download, "install_host_minimap2",
        lambda *a, **k: install_called.append(a),
    )
    monkeypatch.setattr(Minimap2, "_check_tool", lambda self: "minimap2")

    built = []

    def fake_build(self, db, mmi_path, k=0, w=0):
        Path(mmi_path).write_text("INDEX")
        built.append(Path(mmi_path))

    monkeypatch.setattr(Minimap2, "build_index", fake_build)

    prebuild_host_index(cfg, "TestHost")

    expected = Path(str(fasta) + ".mmi")  # k=w=0 default
    assert built == [expected]
    assert expected.exists()
    assert not install_called  # FASTA already present


def test_prebuild_host_index_skips_existing_index(tmp_path, monkeypatch):
    """An existing index is left untouched unless force=True."""
    fasta = tmp_path / "test.fna.gz"
    fasta.write_text(">x\nACGT\n")
    (Path(str(fasta) + ".mmi")).write_text("OLD")

    cfg = EMtraceConfig.create_default()
    cfg.host_databases["TestHost"] = DatabaseEntry(
        db=str(fasta), name="TestHost", search="Minimap2",
    )
    monkeypatch.setattr(Minimap2, "_check_tool", lambda self: "minimap2")

    def fail_build(self, *a, **k):
        raise AssertionError("build_index should not be called")

    monkeypatch.setattr(Minimap2, "build_index", fail_build)

    prebuild_host_index(cfg, "TestHost")  # no exception => skipped
    assert (Path(str(fasta) + ".mmi")).read_text() == "OLD"
