"""Tests for core/exceptions.py — Exception hierarchy."""

from emtrace.core.exceptions import (
    ConfigurationError,
    DatabaseNotConfiguredError,
    DependencyError,
    ExternalToolError,
    EMtraceError,
    SampleFileNotFoundError,
    StepNotFoundError,
    StepValidationError,
)


class TestExceptions:
    """Test exception construction and attributes."""

    def test_base_exception(self):
        e = EMtraceError("test")
        assert str(e) == "test"

    def test_external_tool_error(self):
        e = ExternalToolError(
            "failed", tool="minimap2",
            returncode=1, stderr="error output",
        )
        assert e.tool == "minimap2"
        assert e.returncode == 1
        assert "error output" in e.stderr

    def test_dependency_error(self):
        e = DependencyError(
            "samtools", install_hint="conda install samtools",
        )
        assert e.tool == "samtools"
        assert "conda" in str(e)

    def test_database_not_configured(self):
        e = DatabaseNotConfiguredError("test_db")
        assert e.db_key == "test_db"

    def test_step_not_found(self):
        e = StepNotFoundError("xyz", available=["a", "b"])
        assert e.step_name == "xyz"
        assert "a" in str(e)

    def test_step_validation_error(self):
        e = StepValidationError("pre", reason="no reads")
        assert e.step_name == "pre"
        assert "no reads" in str(e)

    def test_sample_file_not_found(self):
        e = SampleFileNotFoundError("/a/b.fq", step="ann")
        assert e.path == "/a/b.fq"
        assert e.step == "ann"

    def test_configuration_error(self):
        e = ConfigurationError("bad config")
        assert isinstance(e, EMtraceError)
