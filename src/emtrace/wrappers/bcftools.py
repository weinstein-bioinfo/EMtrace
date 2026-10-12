"""Bcftools wrapper — variant calling and consensus generation.

Used in the ref-based assembly step:
  bcftools mpileup | bcftools call -> VCF
  bcftools filter -> filtered VCF
  bcftools consensus -> consensus FASTA
"""

from __future__ import annotations

import subprocess
from functools import lru_cache
from pathlib import Path
from typing import Optional

from emtrace.wrappers.base import ToolWrapper
from emtrace.utils.shell import run_shell


class Bcftools(ToolWrapper):
    """Bcftools variant calling and consensus wrapper."""

    tool_name = "bcftools"
    install_hint = "conda install -c bioconda bcftools"

    def _build_command(self, *, subcmd: str, args: list[str]) -> list[str]:
        return [self._tool_path, subcmd] + args

    def has_mpileup_preset(self, preset: str) -> bool:
        """True if this bcftools knows the ``mpileup -X`` profile *preset*."""
        return preset in _mpileup_presets(self._tool_path)

    def version(self) -> str:
        """The bcftools version string, e.g. ``"1.21"`` (empty if unknown)."""
        return _bcftools_version(self._tool_path)

    def mpileup_call(
        self,
        bam: str | Path,
        ref: str | Path,
        output_vcf: str | Path,
        max_depth: int = 1000,
        ploidy: int = 1,
        no_baq: bool = True,
        split_multiallelic: bool = True,
        preset: str = "",
    ) -> Path:
        """Run mpileup | call | norm to produce a VCF of candidate alleles.

          bcftools mpileup [-X preset] -B -d1000 -a INFO/AD -Ov -f ref bam |
          bcftools call --ploidy 1 -mv -Ov |
          bcftools norm -m -any -f ref -Ov -o output.vcf

        ``no_baq`` passes ``-B`` to mpileup. BAQ is a per-base realignment model
        built for Illumina indel artefacts; on ONT data it downgrades base
        qualities around homopolymers — precisely where the real signal is — and
        suppresses genuine calls.

        ``split_multiallelic`` pipes the calls through ``bcftools norm -m -any``
        so every ALT allele gets its own record. This is required before any
        ``INFO/AD[1]`` filtering: on a packed multiallelic record ``AD[1]`` is
        only the *first* ALT, so a minority allele can be mistaken for the
        dominant one (see docs/refbased.assembly.md).

        ``preset`` is passed as ``mpileup -X`` (a platform profile, see
        ``bcftools mpileup -X list``).  ``ont-sup`` enables indel calling with
        an ONT error model; without a profile mpileup calls almost no indels on
        ONT data.  Explicit options after it (``-B``, ``-d``) still apply.
        """
        output_vcf = Path(output_vcf)
        t = str(self.num_threads)
        baq = "-B " if no_baq else ""
        xopt = f"-X {preset} " if preset else ""
        cmd = (
            f"{self._tool_path} mpileup "
            f"--threads {t} {xopt}{baq}-d{max_depth} -a INFO/AD -Ov "
            f"-f {ref} {bam} | "
            f"{self._tool_path} call "
            f"--threads {t} --ploidy {ploidy} -mv -Ov"
        )
        if split_multiallelic:
            cmd += (
                f" | {self._tool_path} norm "
                f"--threads {t} -m -any -f {ref} -Ov"
            )
        cmd += f" -o {output_vcf}"
        run_shell(cmd, tool_name="bcftools")
        self.logger.info("Variant calling -> %s", output_vcf)
        return output_vcf

    def filter(
        self,
        input_vcf: str | Path,
        output_vcf: str | Path,
        expression: str = "",
        output_type: str = "z",
    ) -> Path:
        """Filter VCF variants by expression.

        Parameters
        ----------
        expression : str
            Filter expression (e.g. ``"INFO/AD[1]>=5"``).
        output_type : str
            ``"z"`` for compressed VCF, ``"v"`` for uncompressed.
        """
        output_vcf = Path(output_vcf)
        self.run(
            subcmd="filter",
            args=[
                "-i", expression,
                f"-O{output_type}",
                str(input_vcf),
                "-o", str(output_vcf),
            ],
        )
        self.logger.info("Filtered VCF -> %s", output_vcf)
        return output_vcf

    def consensus(
        self,
        vcf: str | Path,
        ref: str | Path,
        output_fa: str | Path,
        include: Optional[str] = None,
    ) -> Path:
        """Generate consensus FASTA from VCF + reference.

        ``include`` is a bcftools expression (``-i``) selecting which records are
        applied to the reference. The consensus can carry only one base per
        position, so this must be a *dominant allele* rule (AF strictly above
        0.5); the VCF itself may list additional sub-dominant alleles.
        """
        output_fa = Path(output_fa)
        args = ["-f", str(ref)]
        if include:
            args += ["-i", include]
        args += [str(vcf), "-o", str(output_fa)]
        self.run(subcmd="consensus", args=args)
        self.logger.info("Consensus -> %s", output_fa)
        return output_fa

    def index_vcf(self, vcf: str | Path) -> Path:
        """Index a VCF/BCF file (creates .tbi or .csi)."""
        vcf = Path(vcf)
        self.run(subcmd="index", args=[str(vcf)])
        return vcf


@lru_cache(maxsize=None)
def _mpileup_presets(tool: str) -> frozenset:
    """Profile names listed by ``bcftools mpileup -X list`` (empty on failure)."""
    try:
        res = subprocess.run([tool, "mpileup", "-X", "list"], capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return frozenset()
    names = set()
    for line in (res.stdout + res.stderr).splitlines():
        if line and not line[0].isspace() and not line.startswith("Configuration"):
            names.update(n.strip() for n in line.split(",") if n.strip())
    return frozenset(names)


@lru_cache(maxsize=None)
def _bcftools_version(tool: str) -> str:
    try:
        out = subprocess.run([tool, "--version"], capture_output=True, text=True, timeout=30).stdout
    except (OSError, subprocess.SubprocessError):
        return ""
    first = out.splitlines()[0] if out else ""
    return first.split()[-1] if first else ""
