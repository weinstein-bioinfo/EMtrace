"""Chopper quality/length filtering step.

Registered as :data:`StepName.CHOP`.

Behaviour:
  - ``maxlength == 0`` (from YAML ``maxlength: 0``) is treated as no limit,
    consistent with the Config editor spinbox convention (handled by the
    ChopperConfig validator in ``config.py``).
  - Tool log (chopper stderr) is saved to ``<logs>/chopper.log``.

Output file key: ``"chopped_reads"`` / ``"current_reads"``
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from emtrace.core.exceptions import StepValidationError
from emtrace.steps import StepName, register_step
from emtrace.steps.base import PipelineStep, resolve_preprocess_input
from emtrace.utils.paths import ensure_dir
from emtrace.wrappers.chopper import Chopper

if TYPE_CHECKING:
    from emtrace.core.pipelinejob import PipelineJob


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
