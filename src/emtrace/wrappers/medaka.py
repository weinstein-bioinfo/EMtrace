"""Medaka wrapper — neural-network polishing of an ONT consensus.

Used in the ref-based assembly step (``assembly.medaka``). The step supplies
its own read-to-draft alignment, so only medaka's two back-end commands run:

  medaka inference  BAM consensus.hdf --model M   (per-position probabilities)
  medaka sequence   consensus.hdf draft.fa polished.fa

``medaka_consensus`` is deliberately not used: it re-aligns the reads with
minimap2 defaults (k15), losing the divergent reads the assembly step maps
with small seeds (``assembly.minimap`` k/w).

Medaka pins its own torch/python stack, so it is usually installed in a
separate environment; ``path`` points at that environment's executable.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from emtrace.core.exceptions import DependencyError, ExternalToolError
from emtrace.utils import wsl
from emtrace.wrappers.base import ToolWrapper


class Medaka(ToolWrapper):
    """Medaka consensus polishing wrapper."""

    tool_name = "medaka"
    install_hint = (
        "conda create -n medaka -c conda-forge -c bioconda medaka, then set "
        "assembly.medaka.path to <env>/bin/medaka"
    )

    def __init__(self, num_threads: int = 1, path: str = "", logger=None):
        self._path_hint = path
        super().__init__(num_threads=num_threads, logger=logger)

    def _check_tool(self) -> str:
        """Use the configured executable if given, else look on PATH.

        On Windows a configured path starting with ``/`` is a Linux path
        inside WSL (an executable, or the environment's ``bin`` directory).
        """
        if self._path_hint:
            if wsl.IS_WINDOWS and self._path_hint.startswith("/"):
                for cand in (self._path_hint, self._path_hint.rstrip("/") + "/medaka"):
                    if wsl.is_executable(cand):
                        self.via_wsl = True
                        return cand
                raise DependencyError(
                    self.tool_name, f"assembly.medaka.path not found in WSL: {self._path_hint}",
                )
            p = Path(self._path_hint).expanduser()
            if p.is_dir():
                p = p / "medaka"
            if p.is_file():
                return str(p)
            raise DependencyError(self.tool_name, f"assembly.medaka.path not found: {p}")
        return super()._check_tool()

    def _build_command(self, *, subcmd: list[str], args: list[str]) -> list[str]:
        return [self._tool_path, *subcmd, *args]

    def resolve_model(self, model: str, reads: Optional[str | Path] = None,
                      fallback: str = "") -> str:
        """Resolve *model* to a medaka model (name or file path).

        ``auto`` reads the basecaller model from *reads* (Dorado header /
        ``basecall_model_version_id``); when the reads carry none (simulated,
        older or re-written FASTQ), *fallback* is used, and when that is empty,
        ``""`` — medaka's own default model (no ``--model`` is passed).
        """
        if model != "auto":
            return model
        if reads is not None:
            try:
                out = self.run(subcmd=["tools", "resolve_model"],
                               args=["--auto_model", "consensus", str(reads)])
                resolved = out.strip().splitlines()[-1].strip()
                if resolved:
                    self.logger.info("medaka model from read headers: %s", Path(resolved).name)
                    return resolved
            except (ExternalToolError, IndexError):
                pass
        if fallback:
            self.logger.warning(
                "medaka: no basecaller model in the read headers — using %s "
                "(set assembly.medaka.model to match the basecaller)", fallback,
            )
            return fallback
        self.logger.warning(
            "medaka: no basecaller model in the read headers — using medaka's "
            "default model (set assembly.medaka.model to match the basecaller)"
        )
        return ""

    def inference(self, bam: str | Path, output_hdf: str | Path, model: str,
                  regions: Optional[str | Path] = None,
                  batch_size: int = 100, cpu: bool = False,
                  log_file: Optional[Path] = None) -> Path:
        """``medaka inference``: consensus probabilities for the aligned reads.

        *regions* (a BED file) restricts inference to those draft regions.
        """
        output_hdf = Path(output_hdf)
        output_hdf.unlink(missing_ok=True)
        args = [str(bam), str(output_hdf), "--batch_size", str(batch_size),
                "--bam_workers", str(max(1, min(self.num_threads, 4)))]
        if model:
            args += ["--model", model]
        if regions is not None:
            args += ["--regions", str(regions)]
        if cpu:
            args.append("--cpu")
        self.run(subcmd=["inference"], args=args, log_file=log_file)
        return output_hdf

    def sequence(self, hdf: str | Path, draft: str | Path, output_fa: str | Path,
                 fill_char: str = "", log_file: Optional[Path] = None) -> Path:
        """``medaka sequence``: apply the probabilities to *draft*.

        Draft regions without inference output are filled with the draft
        sequence (medaka's default), or with *fill_char* when given, so the
        contigs keep their full span either way.
        """
        output_fa = Path(output_fa)
        args = [str(hdf), str(draft), str(output_fa),
                "--threads", str(max(1, self.num_threads))]
        if fill_char:
            args += ["--fill_char", fill_char]
        self.run(subcmd=["sequence"], args=args, log_file=log_file)
        self.logger.info("Polished consensus -> %s", output_fa)
        return output_fa
