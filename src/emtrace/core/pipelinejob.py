"""PipelineJob state object — replaces the Perl %opt hash.

Tracks input reads, output paths, run parameters, and intermediate
file locations as the job progresses through pipeline steps.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml

from emtrace import __version__
from emtrace.utils.paths import ensure_dir

logger = logging.getLogger(__name__)

JOB_MANIFEST_NAME = "job.yaml"


@dataclass
class AnnotationRound:
    """Configuration for a single annotation round (ann1 or ann2)."""

    databases: List[Dict[str, str]] = field(default_factory=list)
    # Each entry: {"db": "...", "name": "...", "search": "minimap|blastn|...", "target": "vi|ab|un|..."}


@dataclass
class PipelineJob:
    """Mutable state object that accompanies a job through the pipeline.

    Replaces the Perl ``%opt`` hash.  Each pipeline step reads from and
    writes to this object so downstream steps can find intermediate files.
    """

    # -- Identity --
    name: str = ""
    input_reads: Optional[Path] = None
    # Optional second mate for paired-end Illumina input. When set, the job
    # is treated as `platform == "illumina_pe"` regardless of the explicit
    # platform field. Nanopore jobs leave this as None.
    input_reads_r2: Optional[Path] = None
    # Sequencing platform — controls aligner preset selection in the
    # assembly step and which preprocessing steps are applicable.
    # Currently supported: "nanopore" (default), "illumina_pe".
    platform: str = "nanopore"

    # -- Directory layout --
    output_dir: Path = Path(".")
    logs_dir: str = "logs"
    results_dir: str = "results"

    # -- Pipeline control --
    pipe: str = "main"

    # -- Background filtering --
    bgfilter_databases: List[Dict[str, str]] = field(default_factory=list)

    # -- Annotation --
    ann1: AnnotationRound = field(default_factory=AnnotationRound)
    ann2: AnnotationRound = field(default_factory=AnnotationRound)
    anns: str = "0"
    append: bool = False

    # -- Runtime --
    num_threads: int = 16
    verbose: int = 0

    # -- Intermediate file tracking --
    # Populated by each step as it produces output
    files: Dict[str, Path] = field(default_factory=dict)

    # -- Arbitrary metadata --
    metadata: Dict[str, Any] = field(default_factory=dict)

    # ------------------------------------------------------------------ #
    # Directory helpers
    # ------------------------------------------------------------------ #

    @property
    def reads_dir(self) -> Path:
        return self.output_dir / "reads"

    @property
    def annotation_dir(self) -> Path:
        return self.output_dir / "annotation"

    @property
    def assembly_dir(self) -> Path:
        return self.output_dir / "assembly"

    @property
    def reports_dir(self) -> Path:
        return self.output_dir / self.results_dir

    @property
    def log_path(self) -> Path:
        return self.output_dir / self.logs_dir

    def init_dirs(self) -> None:
        """Create the standard output directory tree and write the job manifest."""
        for d in (self.reads_dir, self.annotation_dir, self.assembly_dir,
                  self.reports_dir, self.log_path):
            ensure_dir(d)
        self._write_job_manifest()

    # ------------------------------------------------------------------ #
    # Job manifest (input-reads ↔ output-dir link)
    # ------------------------------------------------------------------ #

    @property
    def manifest_path(self) -> Path:
        return self.output_dir / JOB_MANIFEST_NAME

    def _write_job_manifest(self) -> None:
        """Write ``<output_dir>/job.yaml`` so Open Job can recover input_reads.

        Preserves the original ``created_at`` if a manifest already exists,
        and only updates ``input_reads`` and ``emtrace_version``.
        """
        existing: Dict[str, Any] = {}
        if self.manifest_path.exists():
            existing = self.read_job_manifest(self.output_dir) or {}
        payload = {
            "name": self.name,
            "input_reads": str(self.input_reads) if self.input_reads else "",
            "input_reads_r2": str(self.input_reads_r2) if self.input_reads_r2 else "",
            "platform": self.platform,
            "created_at": existing.get("created_at") or datetime.now().isoformat(timespec="seconds"),
            "emtrace_version": __version__,
        }
        try:
            with self.manifest_path.open("w") as f:
                yaml.safe_dump(payload, f, sort_keys=False)
        except OSError as exc:
            logger.warning("Could not write job manifest %s: %s", self.manifest_path, exc)

    @classmethod
    def read_job_manifest(cls, outdir: Path) -> Optional[Dict[str, Any]]:
        """Read ``<outdir>/job.yaml`` if present; return ``None`` on missing/invalid."""
        path = outdir / JOB_MANIFEST_NAME
        if not path.exists():
            return None
        try:
            with path.open("r") as f:
                data = yaml.safe_load(f) or {}
        except (OSError, yaml.YAMLError) as exc:
            logger.warning("Could not read job manifest %s: %s", path, exc)
            return None
        if not isinstance(data, dict):
            return None
        return data

    # ------------------------------------------------------------------ #
    # File registry helpers
    # ------------------------------------------------------------------ #

    def set_file(self, key: str, path: Path) -> None:
        """Register an output file produced by a step."""
        self.files[key] = path

    def get_file(self, key: str) -> Optional[Path]:
        """Get a registered file path, or None."""
        return self.files.get(key)

    def has_file(self, key: str) -> bool:
        """Check whether a file has been registered *and* exists on disk."""
        p = self.files.get(key)
        return p is not None and p.exists()

    # ------------------------------------------------------------------ #
    # Convenience: current reads path
    # ------------------------------------------------------------------ #

    @property
    def current_reads(self) -> Optional[Path]:
        """The most-recently-produced reads file.

        Steps should update this via ``set_file("current_reads", path)``.
        Falls back to ``input_reads`` if no steps have run yet.
        """
        return self.files.get("current_reads", self.input_reads)

    @property
    def is_paired(self) -> bool:
        """True when an R2 mate is registered (paired-end Illumina input)."""
        return self.input_reads_r2 is not None
