"""Abstract base class for pipeline steps."""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import TYPE_CHECKING, Optional

from emtrace.core.exceptions import StepValidationError

if TYPE_CHECKING:
    from emtrace.core.config import EMtraceConfig
    from emtrace.core.pipelinejob import PipelineJob


def apply_preprocess_subsample(
    sample: "PipelineJob",
    config: "EMtraceConfig",
    logger: logging.Logger,
) -> Path:
    """Subsample input reads before annotation.

    Idempotent: if ``reads/input.subsample.fq`` already exists on disk the
    function re-registers it as ``current_reads`` and returns without
    re-running subsampling.

    Used by the annotate step only.  Preprocessing steps (PORE/CHOP/FLT)
    use :func:`resolve_preprocess_input` instead, which resolves input from
    ``input_reads`` rather than ``current_reads``.

    Returns the path to use as ``reads_in`` (subsampled file or
    ``current_reads`` when subsampling is disabled).
    """
    sub = config.preprocess.subsample
    if not sub.active or sub.size <= 0:
        return Path(sample.current_reads)

    sub_path = sample.reads_dir / "input.subsample.fq"

    if sub_path.exists():
        sample.set_file("current_reads", sub_path)
        logger.debug("Subsampled reads already present: %s", sub_path)
        return sub_path

    from emtrace.io.fastq import subsample_fastq
    from emtrace.utils.paths import ensure_dir
    ensure_dir(sample.reads_dir)
    n = subsample_fastq(Path(sample.current_reads), sub_path, size=sub.size)
    logger.info("Subsampled %d reads (size=%s) → %s", n, sub.size, sub_path)
    sample.set_file("current_reads", sub_path)
    return sub_path


def resolve_preprocess_input(
    sample: "PipelineJob",
    config: "EMtraceConfig",
    logger: logging.Logger,
    predecessors: tuple = (),
) -> Path:
    """Resolve the correct input reads for a preprocessing step.

    Always derives from ``input_reads`` (not ``current_reads``), so that
    re-running a step after a previous full run still uses the correct source.

    Resolution order:
      1. First *predecessors* filename found in ``sample.reads_dir`` (non-empty).
      2. ``reads/input.subsample.fq`` when subsampling is active and the file
         exists (created by a prior step in the same run).
      3. Create ``reads/input.subsample.fq`` from ``input_reads`` when
         subsampling is active but the file does not yet exist.
      4. ``input_reads`` when subsampling is disabled.

    Parameters
    ----------
    predecessors:
        Filenames (not full paths) to check in ``sample.reads_dir`` in order.
        The first existing non-empty file is returned without further logic.
        Pass ``()`` (the default) for Porechop, which has no predecessor.
    """
    for fname in predecessors:
        p = sample.reads_dir / fname
        if p.exists() and p.stat().st_size > 0:
            logger.debug("Preprocessing input resolved from predecessor: %s", p)
            return p

    # No predecessor on disk — resolve from raw input_reads
    sub = config.preprocess.subsample
    if not sub.active or sub.size <= 0:
        return Path(sample.input_reads)

    sub_path = sample.reads_dir / "input.subsample.fq"
    if sub_path.exists():
        logger.debug("Subsampled reads already present: %s", sub_path)
        return sub_path

    from emtrace.io.fastq import subsample_fastq
    from emtrace.utils.paths import ensure_dir
    ensure_dir(sample.reads_dir)
    n = subsample_fastq(Path(sample.input_reads), sub_path, size=sub.size)
    logger.info("Subsampled %d reads (size=%s) → %s", n, sub.size, sub_path)
    sample.set_file("current_reads", sub_path)
    return sub_path


class PipelineStep(ABC):
    """Abstract base class for all pipeline steps.

    Each step must implement ``run()`` and ``validate_inputs()``.
    The ``execute()`` method provides the full lifecycle:
    validate -> pre_run -> run -> post_run.
    """

    step_name: str = "base"
    step_description: str = "Base pipeline step"

    def __init__(
        self,
        config: EMtraceConfig,
        logger: Optional[logging.Logger] = None,
    ):
        self.config = config
        self.logger = logger or logging.getLogger(
            f"emtrace.steps.{self.step_name}"
        )

    @abstractmethod
    def run(self, sample: PipelineJob, **kwargs) -> PipelineJob:
        """Execute the step's core logic.

        Parameters
        ----------
        sample : PipelineJob
            Current sample state.

        Returns
        -------
        PipelineJob
            Updated sample state.
        """

    @abstractmethod
    def validate_inputs(self, sample: PipelineJob, **kwargs) -> bool:
        """Check that required inputs exist before running.

        Returns
        -------
        bool
            True if inputs are valid.
        """

    def pre_run(self, sample: PipelineJob, **kwargs) -> None:
        """Hook called before ``run()``. Override for setup."""
        self.logger.info("Starting step: %s", self.step_name)

    def post_run(self, sample: PipelineJob, **kwargs) -> None:
        """Hook called after ``run()``. Override for cleanup."""
        self.logger.info("Completed step: %s", self.step_name)

    def execute(self, sample: PipelineJob, **kwargs) -> PipelineJob:
        """Full step lifecycle: validate -> pre_run -> run -> post_run."""
        if not self.validate_inputs(sample, **kwargs):
            raise StepValidationError(self.step_name)
        self.pre_run(sample, **kwargs)
        sample = self.run(sample, **kwargs)
        self.post_run(sample, **kwargs)
        return sample
