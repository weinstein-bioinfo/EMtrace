"""Porechop / Porechop_ABI wrappers — Nanopore adapter removal.

Both tools write progress and adapter information to stderr.  Pass
``log_file`` to have that output saved to disk.  ``PorechopABI`` also
accepts ``adapters_file`` (written by the step layer after parsing).

Both tools write a temp file (``TEMP_<pid>.fastq``) to the working
directory, so they run in the output directory with absolute paths:
the caller's working directory may be read-only (``/`` for a GUI app
launched from Finder) or shared with another run.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Optional

from emtrace.wrappers.base import ToolWrapper, count_fastq_reads

if TYPE_CHECKING:
    from emtrace.core.config import PorechopABIConfig


class Porechop(ToolWrapper):
    """Original Porechop adapter removal (fallback when porechop_abi unavailable)."""

    tool_name = "porechop"
    install_hint = "conda install -c bioconda porechop"

    def _build_command(
        self,
        *,
        input_fq: str | Path,
        output_fq: str | Path,
        extra_params: str = "",
    ) -> list[str]:
        cmd = [
            self._tool_path,
            "-t", str(self.num_threads),
            "-i", str(input_fq),
            "-o", str(output_fq),
            "--format", "fastq.gz",
        ]
        if extra_params:
            cmd.extend(extra_params.split())
        return cmd

    def trim(
        self,
        input_fq: str | Path,
        output_fq: str | Path,
        extra_args: list[str] | None = None,
        log_file: Optional[Path] = None,
    ) -> Path:
        """Remove adapters using the built-in adapter database."""
        input_fq = Path(input_fq).resolve()
        output_fq = Path(output_fq).resolve()
        n_in = count_fastq_reads(input_fq)
        extra_params = " ".join(extra_args) if extra_args else ""
        self.run(
            input_fq=input_fq,
            output_fq=output_fq,
            extra_params=extra_params,
            log_file=log_file,
            cwd=output_fq.parent,
        )
        n_out = count_fastq_reads(output_fq)
        self.logger.info("porechop: kept %d / %d reads", n_out, n_in)
        return output_fq

    def trim_from_config(
        self,
        input_fq: str | Path,
        output_fq: str | Path,
        cfg: "PorechopABIConfig",
        log_file: Optional[Path] = None,
    ) -> Path:
        """Trim using a subset of PorechopABIConfig options supported by original Porechop.

        Ab-initio flags are ignored (not supported by original Porechop).
        """
        return self.trim(
            input_fq=input_fq,
            output_fq=output_fq,
            extra_args=cfg.to_cli_args_porechop(),
            log_file=log_file,
        )


class PorechopABI(ToolWrapper):
    """Porechop_ABI — ab initio adapter discovery + trimming."""

    tool_name = "porechop_abi"
    install_hint = "conda install -c conda-forge -c bioconda porechop_abi"

    def _build_command(
        self,
        *,
        input_fq: str | Path,
        output_fq: str | Path,
        extra_args: list[str] | None = None,
    ) -> list[str]:
        cmd = [
            self._tool_path,
            "-t", str(self.num_threads),
            "-i", str(input_fq),
            "-o", str(output_fq),
            "--format", "fastq.gz",
        ]
        if extra_args:
            cmd.extend(extra_args)
        return cmd

    def trim(
        self,
        input_fq: str | Path,
        output_fq: str | Path,
        extra_args: list[str] | None = None,
        log_file: Optional[Path] = None,
    ) -> Path:
        """Remove adapters using ab initio discovery.

        Parameters
        ----------
        log_file : Path, optional
            If provided, stderr from the tool is written to this file.
        """
        input_fq = Path(input_fq).resolve()
        output_fq = Path(output_fq).resolve()
        n_in = count_fastq_reads(input_fq)
        self.run(
            input_fq=input_fq,
            output_fq=output_fq,
            extra_args=extra_args,
            log_file=log_file,
            cwd=output_fq.parent,
        )
        n_out = count_fastq_reads(output_fq)
        self.logger.info("porechop_abi: kept %d / %d reads", n_out, n_in)
        return output_fq

    def trim_from_config(
        self,
        input_fq: str | Path,
        output_fq: str | Path,
        cfg: "PorechopABIConfig",
        log_file: Optional[Path] = None,
        adapters_file: Optional[Path] = None,  # accepted for API consistency; written by step layer
    ) -> Path:
        """Trim using settings from a :class:`PorechopABIConfig` instance."""
        return self.trim(
            input_fq=input_fq,
            output_fq=output_fq,
            extra_args=cfg.to_cli_args(),
            log_file=log_file,
        )
