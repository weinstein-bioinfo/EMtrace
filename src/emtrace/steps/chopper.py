"""Chopper quality/length filtering step.

Registered as :data:`StepName.CHOP`.

Behaviour:
  - ``maxlength == 0`` (from YAML ``maxlength: 0``) is treated as no limit,
    consistent with the Config editor spinbox convention (handled by the
    ChopperConfig validator in ``config.py``).
  - Tool log (chopper stderr) is saved to ``<logs>/chopper.log``.
  - Warns when at least ``SHORT_READ_WARN_FRAC`` of the input reads are
    shorter than ``minlength``: the default (400 bp) discards much of a
    short-read library.

Output file key: ``"chopped_reads"`` / ``"current_reads"``
"""

from __future__ import annotations

import statistics
from array import array
from pathlib import Path
from typing import TYPE_CHECKING

from emtrace.core.exceptions import StepValidationError
from emtrace.io.fastq import iter_fastq
from emtrace.steps import StepName, register_step
from emtrace.steps.base import PipelineStep, resolve_preprocess_input
from emtrace.utils.paths import ensure_dir
from emtrace.wrappers.chopper import Chopper

if TYPE_CHECKING:
    from emtrace.core.pipelinejob import PipelineJob

SHORT_READ_WARN_FRAC = 0.25   # warn when this share of reads is below minlength


@register_step(StepName.CHOP)
class ChopperStep(PipelineStep):
    """Quality/length filtering with Chopper."""

    step_name = "chopper"
    step_description = "Quality/length filtering (Chopper)"

    def validate_inputs(self, sample: PipelineJob, **kwargs) -> bool:
        if sample.input_reads is None or not Path(sample.input_reads).exists():
            raise StepValidationError(self.step_name, "No input reads available")
        return True

    def run(self, sample: PipelineJob, **kwargs) -> PipelineJob:
        cfg = self.config.preprocess.chopper
        ensure_dir(sample.reads_dir)
        ensure_dir(sample.log_path)
        reads_in = resolve_preprocess_input(
            sample, self.config, self.logger,
            predecessors=("reads.porechop.fq.gz",),
        )

        reads_out = sample.reads_dir / "reads.chopped.fq.gz"
        log_file = sample.log_path / "chopper.log"

        self._warn_short_reads(reads_in, cfg.minlength)
        self.logger.info("Running Chopper: %s -> %s", reads_in, reads_out)
        wrapper = Chopper(num_threads=sample.num_threads)
        wrapper.filter_from_config(
            input_fq=reads_in,
            output_fq=reads_out,
            cfg=cfg,
            log_file=log_file,
        )
        sample.set_file("chopped_reads", reads_out)
        sample.set_file("current_reads", reads_out)
        return sample

    def _warn_short_reads(self, reads_in: Path, minlength: int) -> None:
        """Warn when many input reads are shorter than *minlength*.

        Checked on chopper's input, so the share counts only what the length
        cut-off removes (chopper's own in/out counts mix in the quality filter).
        """
        if minlength <= 0:
            return
        lengths = array("L", (len(seq) for _, seq, _, _ in iter_fastq(reads_in)))
        if not lengths:
            return
        frac = sum(1 for n in lengths if n < minlength) / len(lengths)
        if frac >= SHORT_READ_WARN_FRAC:
            self.logger.warning(
                "%.0f%% of reads < minlength %d bp (median %d bp): "
                "consider lowering chopper.minlength",
                100 * frac, minlength, statistics.median(lengths),
            )
