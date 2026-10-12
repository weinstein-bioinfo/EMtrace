"""minimap2 command construction: -N is emitted for any non-negative n."""

import pytest

from emtrace.wrappers.minimap2 import Minimap2


@pytest.fixture
def mm2(monkeypatch):
    monkeypatch.setattr(Minimap2, "_check_tool", lambda self: "minimap2")
    return Minimap2(num_threads=1)


@pytest.mark.parametrize("n, expected", [(0, "0"), (5, "5"), (50, "50")])
def test_secondary_alignments_flag_emitted(mm2, n, expected):
    cmd = mm2._build_command(db="db.mmi", query="reads.fq", n=n)
    assert cmd[cmd.index("-N") + 1] == expected


def test_negative_n_omits_flag(mm2):
    cmd = mm2._build_command(db="db.mmi", query="reads.fq", n=-1)
    assert "-N" not in cmd
