"""Base class for external tool wrappers."""

import gzip
import logging
import shutil
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Optional, Sequence

from emtrace.core.exceptions import DependencyError
from emtrace.utils.shell import run as shell_run


def count_fastq_reads(fq_path: Path) -> int:
    """Count reads in a FASTQ file (plain or gzipped) by counting header lines."""
    fq_path = Path(fq_path)
    count = 0
    open_fn = gzip.open if str(fq_path).endswith(".gz") else open
    with open_fn(fq_path, "rt", errors="replace") as fh:
        for i, _ in enumerate(fh):
            if i % 4 == 0:
                count += 1
    return count


class ToolWrapper(ABC):
    """Abstract base for bioinformatics tool wrappers.

    Subclasses must set ``tool_name`` and implement ``_build_command``.
    The wrapper checks tool availability at init and provides a
    ``run()`` method that logs commands, captures output, and raises
    ``ExternalToolError`` on failure.
    """

    tool_name: str = ""
    install_hint: str = ""

    def __init__(
        self,
        num_threads: int = 1,
        logger: Optional[logging.Logger] = None,
    ):
        self.num_threads = num_threads
        self.logger = logger or logging.getLogger(
            f"emtrace.wrappers.{self.tool_name}"
        )
        self._tool_path = self._check_tool()
        self.last_stderr: str = ""

    def _check_tool(self) -> str:
        """Verify the tool is on PATH; return its absolute path."""
        path = shutil.which(self.tool_name)
        if path is None:
            raise DependencyError(self.tool_name, self.install_hint)
        self.logger.debug("Found %s at %s", self.tool_name, path)
        return path

    @abstractmethod
    def _build_command(self, **kwargs) -> list[str]:
        """Build the command-line argument list."""

    def run(
        self,
        *,
        extra_args: Optional[Sequence[str]] = None,
        cwd: Optional[Path] = None,
        stdin_data: Optional[str] = None,
        timeout: Optional[int] = None,
        check: bool = True,
        log_file: Optional[Path] = None,
        **kwargs,
    ) -> str:
        """Build and execute the tool command.

        Parameters
        ----------
        log_file : Path, optional
            If provided, stderr from the tool is written to this file after
            the run completes (the file is created/overwritten).

        Returns
        -------
        str
            stdout of the completed process.
        """
        cmd = self._build_command(**kwargs)
        if extra_args:
            cmd.extend(extra_args)

        result = shell_run(
            cmd,
            cwd=cwd,
            check=check,
            capture=True,
            stdin_data=stdin_data,
            timeout=timeout,
            tool_name=self.tool_name,
        )

        self.last_stderr = result.stderr or ""

        if self.last_stderr:
            for line in self.last_stderr.splitlines():
                line = line.strip()
                if not line:
                    continue
                if "[WARNING]" in line or "warning" in line.lower():
                    self.logger.warning("%s: %s", self.tool_name, line)
                else:
                    self.logger.debug("%s: %s", self.tool_name, line)

        if log_file is not None and self.last_stderr:
            Path(log_file).write_text(self.last_stderr, encoding="utf-8")

        return result.stdout or ""
