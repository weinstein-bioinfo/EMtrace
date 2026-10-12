"""Tests for wrappers/samtools.py — command construction."""

from unittest.mock import MagicMock, patch

import pytest

from emtrace.wrappers.samtools import Samtools


@pytest.fixture
def mock_samtools():
    """Create Samtools wrapper with mocked tool check."""
    with patch("shutil.which", return_value="/usr/bin/samtools"):
        return Samtools(num_threads=4)


class TestSamtoolsCommandConstruction:
    """Test that samtools subcommands build correct args."""

    def test_sort_command(self, mock_samtools):
        cmd = mock_samtools._build_command(
            subcmd="sort",
            args=["-O", "bam", "-o", "out.bam", "in.sam"],
        )
        assert cmd[0] == "/usr/bin/samtools"
        assert cmd[1] == "sort"
        assert "-O" in cmd
        assert "bam" in cmd

    def test_index_command(self, mock_samtools):
        cmd = mock_samtools._build_command(
            subcmd="index",
            args=["-@", "4", "test.bam"],
        )
        assert cmd[1] == "index"
        assert "test.bam" in cmd

    def test_view_command(self, mock_samtools):
        cmd = mock_samtools._build_command(
            subcmd="view",
            args=["-O", "bam", "-o", "out.bam", "in.sam"],
        )
        assert cmd[1] == "view"
