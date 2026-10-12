"""Tests for wrappers/base.py — ToolWrapper base class."""

from unittest.mock import patch

import pytest

from emtrace.core.exceptions import DependencyError
from emtrace.wrappers.base import ToolWrapper


class _DummyWrapper(ToolWrapper):
    """Minimal concrete ToolWrapper for testing."""

    tool_name = "echo"
    install_hint = "should be available everywhere"

    def _build_command(self, **kwargs):
        msg = kwargs.get("msg", "hello")
        return [self._tool_path, msg]


class TestToolWrapperInit:
    """Test wrapper initialization."""

    def test_finds_tool_on_path(self):
        w = _DummyWrapper(num_threads=4)
        assert w._tool_path is not None
        assert w.num_threads == 4

    def test_missing_tool_raises(self):
        class _Missing(ToolWrapper):
            tool_name = "nonexistent_tool_xyz_999"
            install_hint = "install it"

            def _build_command(self, **kwargs):
                return [self._tool_path]

        with pytest.raises(DependencyError):
            _Missing()


class TestToolWrapperRun:
    """Test command execution."""

    def test_run_returns_stdout(self):
        w = _DummyWrapper()
        out = w.run(msg="world")
        assert "world" in out

    def test_build_command(self):
        w = _DummyWrapper()
        cmd = w._build_command(msg="test")
        assert cmd[-1] == "test"
