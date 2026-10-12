"""The taxonomy CLI flags are gone: taxonomy follows the annotation database."""

import pytest
from click.testing import CliRunner

from emtrace.__main__ import cli


@pytest.mark.parametrize("flag", ["--NCBItax", "--SILVAtax"])
def test_taxonomy_flags_rejected(tmp_path, flag):
    reads = tmp_path / "reads.fq"
    reads.write_text("@r1\nACGT\n+\nIIII\n")
    result = CliRunner().invoke(
        cli,
        ["run", "-1", str(reads), "-s", "s1", "-r", str(tmp_path / "out"),
         "--anns", "16S.silva", flag],
    )
    assert result.exit_code == 2
    assert "No such option" in result.output
    assert flag in result.output
    # Rejected at argument parsing, before any output is created.
    assert not (tmp_path / "out").exists()


def test_help_does_not_mention_taxonomy_flags():
    result = CliRunner().invoke(cli, ["run", "--help"])
    assert result.exit_code == 0
    assert "NCBItax" not in result.output
    assert "SILVAtax" not in result.output
