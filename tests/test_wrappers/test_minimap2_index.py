"""Tests for safe (array-job-friendly) minimap2 host-index building."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from emtrace.wrappers import minimap2 as mm2_mod
from emtrace.wrappers.minimap2 import Minimap2


@pytest.fixture
def mm2(monkeypatch):
    """A Minimap2 wrapper that doesn't require the binary to be installed."""
    monkeypatch.setattr(Minimap2, "_check_tool", lambda self: "minimap2")
    return Minimap2(num_threads=1)


def _index_target(cmd: list) -> Path:
    """Return the path following ``-d`` in a minimap2 command list."""
    return Path(cmd[cmd.index("-d") + 1])


def test_build_index_atomic_and_unique_temp(mm2, monkeypatch, tmp_path):
    """build_index writes a per-process temp then atomically renames it."""
    seen_targets = []

    def fake_run(cmd, **kwargs):
        target = _index_target(cmd)
        seen_targets.append(target)
        # The temp must be unique (carry the pid) and not the final path.
        assert target.name != "host.fna.gz.mmi"
        assert ".tmp." in target.name
        target.write_text("INDEX")  # simulate minimap2 -d output
        return SimpleNamespace(stdout="", stderr="")

    monkeypatch.setattr(mm2_mod, "shell_run", fake_run)

    db = tmp_path / "host.fna.gz"
    db.write_text(">x\nACGT\n")
    mmi = tmp_path / "host.fna.gz.mmi"

    mm2.build_index(db, mmi)
    mm2.build_index(db, mmi)

    assert mmi.exists() and mmi.read_text() == "INDEX"
    # No temp files left behind, and the two builds used distinct temp names.
    assert not list(tmp_path.glob("*.tmp.*"))
    assert seen_targets[0] != seen_targets[1]


def test_build_index_cleans_temp_on_failure(mm2, monkeypatch, tmp_path):
    """A failed build leaves no partial index at the final path or temp."""
    def boom(cmd, **kwargs):
        _index_target(cmd).write_text("PARTIAL")
        raise RuntimeError("minimap2 died")

    monkeypatch.setattr(mm2_mod, "shell_run", boom)

    db = tmp_path / "host.fna.gz"
    db.write_text(">x\nACGT\n")
    mmi = tmp_path / "host.fna.gz.mmi"

    with pytest.raises(RuntimeError):
        mm2.build_index(db, mmi)

    assert not mmi.exists()
    assert not list(tmp_path.glob("*.tmp.*"))


def test_host_filter_builds_only_when_index_missing(mm2, monkeypatch, tmp_path):
    """host_filter builds the index once (under lock) and reuses it after."""
    build_calls = []

    def fake_build(self, db, mmi_path, k=0, w=0):
        build_calls.append(Path(mmi_path))
        Path(mmi_path).write_text("INDEX")

    # The mapping run itself just needs to produce empty PAF output.
    monkeypatch.setattr(Minimap2, "build_index", fake_build)
    monkeypatch.setattr(
        mm2_mod, "shell_run",
        lambda cmd, **kw: SimpleNamespace(stdout="", stderr=""),
    )

    db = tmp_path / "host.fna.gz"
    db.write_text(">x\nACGT\n")
    mmi = tmp_path / "host.fna.gz.mmi"
    reads = tmp_path / "reads.fq"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    out = tmp_path / "out.paf"

    # First call: index absent -> built once.
    mm2.host_filter(db, reads, out, mmi_path=mmi, save_mmi=True)
    assert build_calls == [mmi]
    assert mmi.exists()

    # Second call: index present -> no rebuild.
    mm2.host_filter(db, reads, out, mmi_path=mmi, save_mmi=True)
    assert build_calls == [mmi]  # unchanged
