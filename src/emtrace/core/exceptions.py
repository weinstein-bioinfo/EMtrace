"""EMtrace exception hierarchy."""


class EMtraceError(Exception):
    """Base exception for all EMtrace errors."""


class ConfigurationError(EMtraceError):
    """Invalid or missing configuration."""


class SampleError(EMtraceError):
    """Error related to sample state or files."""


class SampleNotInitializedError(SampleError):
    """Sample has not been initialized."""


class SampleFileNotFoundError(SampleError):
    """Expected sample file is missing."""

    def __init__(self, path: str, step: str = ""):
        msg = f"File not found: {path}"
        if step:
            msg += f" (required by step '{step}')"
        super().__init__(msg)
        self.path = path
        self.step = step


class ExternalToolError(EMtraceError):
    """An external tool (minimap2, blastn, etc.) failed."""

    def __init__(
        self,
        message: str,
        tool: str = "",
        returncode: int | None = None,
        stderr: str = "",
    ):
        super().__init__(message)
        self.tool = tool
        self.returncode = returncode
        self.stderr = stderr


class DependencyError(EMtraceError):
    """A required external tool is not installed or not on PATH."""

    def __init__(self, tool: str, install_hint: str = ""):
        msg = f"Required tool not found: {tool}"
        if install_hint:
            msg += f"\n  Install: {install_hint}"
        super().__init__(msg)
        self.tool = tool
        self.install_hint = install_hint


class DatabaseNotConfiguredError(EMtraceError):
    """A required database is not configured."""

    def __init__(self, db_key: str):
        super().__init__(f"Database not configured: {db_key}")
        self.db_key = db_key


class DatabaseNotFoundError(EMtraceError):
    """A configured database file does not exist on disk."""

    def __init__(self, db_key: str, path: str):
        super().__init__(f"Database '{db_key}' not found at: {path}")
        self.db_key = db_key
        self.path = path


class StepNotFoundError(EMtraceError):
    """Unknown pipeline step name."""

    def __init__(self, step_name: str, available: list[str] | None = None):
        msg = f"Unknown step: '{step_name}'"
        if available:
            msg += f"\n  Available: {', '.join(available)}"
        super().__init__(msg)
        self.step_name = step_name
        self.available = available


class StepValidationError(EMtraceError):
    """Step input validation failed."""

    def __init__(self, step_name: str, reason: str = ""):
        msg = f"Validation failed for step '{step_name}'"
        if reason:
            msg += f": {reason}"
        super().__init__(msg)
        self.step_name = step_name
        self.reason = reason
