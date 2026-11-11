"""Chopper wrapper — Nanopore read quality/length filtering.

Chopper reads from stdin and writes to stdout, so this wrapper
pipes data through the tool using the shell.

Log file: chopper writes progress/stats to stderr.  Pass ``log_file``
to redirect stderr to a file (``2>log_file`` in the shell command).
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Optional

from emtrace.wrappers.base import ToolWrapper, count_fastq_reads
from emtrace.utils.shell import run_shell, shell_path

if TYPE_CHECKING:
    from emtrace.core.config import ChopperConfig


class Chopper(ToolWrapper):
    """Chopper quality/length filtering wrapper."""

    tool_name = "chopper"
    install_hint = "conda install -c bioconda chopper"

    def _build_command(self, **kwargs) -> list[str]:
        # Not used directly — filter() uses a shell pipeline
        raise NotImplementedError("Use filter() instead")

    def filter(
        self,
        input_fq: str | Path,
        output_fq: str | Path,
        quality: int | None = None,
        min_length: int | None = None,
        max_length: int | float | None = None,
        trim: int | None = None,
        extra_params: str = "",
        gzip_cmd: str = "pigz",
        log_file: Optional[Path] = None,
    ) -> Path:
        """Filter reads with Chopper.

        Handles gzipped input transparently.  Output is always gzipped.

        Parameters
        ----------
        input_fq : Path
            Input FASTQ file (plain or .gz).
        output_fq : Path
            Output FASTQ file (will be gzipped).
        quality : int, optional
            Minimum quality score.
        min_length : int, optional
            Minimum read length (bp).
        max_length : int or float, optional
            Maximum read length (bp).  ``0`` and ``float("inf")`` both mean
            no limit (``--maxlength`` flag is omitted from the command).
        trim : int, optional
            Trim N bases from start/end.
        extra_params : str
            Additional chopper parameters string.
        gzip_cmd : str
            Compression command (pigz or gzip).
        log_file : Path, optional
            If given, chopper stderr is redirected to this file.

        Returns
        -------
        Path
            Path to the output file.
        """
        input_fq = Path(input_fq)
        output_fq = Path(output_fq)

        params = extra_params
        if quality is not None:
            params += f" --quality {quality}"
        if min_length is not None:
            params += f" --minlength {min_length}"
        # 0 and inf both mean "no limit" — omit the flag
        if max_length is not None and max_length != 0 and max_length != float("inf"):
            params += f" --maxlength {int(max_length)}"
        if trim is not None:
            params += f" --trim {trim}"

        if str(input_fq).endswith(".gz"):
            cat_cmd = f"{gzip_cmd} -dc {shell_path(input_fq)}"
        else:
            cat_cmd = f"cat {shell_path(input_fq)}"

        stderr_redirect = f" 2>{shell_path(log_file)}" if log_file else ""

        cmd = (
            f"{cat_cmd} | "
            f"{shell_path(self._tool_path)} {params} -t {self.num_threads}{stderr_redirect} | "
            f"{gzip_cmd} -c > {shell_path(output_fq)}"
        )

        n_in = count_fastq_reads(input_fq)
        run_shell(cmd, tool_name="chopper")
        n_out = count_fastq_reads(output_fq)
        self.logger.info("chopper: kept %d / %d reads", n_out, n_in)
        return output_fq

    def filter_from_config(
        self,
        input_fq: str | Path,
        output_fq: str | Path,
        cfg: "ChopperConfig",
        gzip_cmd: str = "pigz",
        log_file: Optional[Path] = None,
    ) -> Path:
        """Filter reads using settings from a :class:`ChopperConfig` instance."""
        args = cfg.to_cli_args()
        extra_params = " ".join(args)
        return self.filter(
            input_fq=input_fq,
            output_fq=output_fq,
            extra_params=extra_params,
            gzip_cmd=gzip_cmd,
            log_file=log_file,
        )
