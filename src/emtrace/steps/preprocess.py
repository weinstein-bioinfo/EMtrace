"""Preprocessing step — superseded by PorechopStep (PORE) + ChopperStep (CHOP).

This module is kept for reference but is no longer registered in the pipeline.
New code should use:
  - ``emtrace.steps.porechop.PorechopStep``  (StepName.PORE)
  - ``emtrace.steps.chopper.ChopperStep``     (StepName.CHOP)
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from emtrace.core.config import ChopperConfig, PorechopABIConfig
from emtrace.core.exceptions import DependencyError, StepValidationError
from emtrace.io.fastq import copy_fastq_gz
from emtrace.steps.base import PipelineStep
from emtrace.utils.paths import ensure_dir
from emtrace.wrappers.chopper import Chopper
from emtrace.wrappers.porechop import Porechop, PorechopABI

if TYPE_CHECKING:
    from emtrace.core.pipelinejob import PipelineJob


class PreprocessStep(PipelineStep):
    """Preprocess raw Nanopore reads with Porechop_ABI and/or Chopper.
    """

    step_name = "preprocess"
    step_description = (
        "Adapter removal and quality filtering (Porechop_ABI + Chopper)"
    )

    def validate_inputs(self, sample: PipelineJob, **kwargs) -> bool:
        if sample.input_reads is None:
            raise StepValidationError(
                self.step_name, "sample.input_reads is not set"
            )
        if not Path(sample.input_reads).exists():
            raise StepValidationError(
                self.step_name,
                f"Input reads not found: {sample.input_reads}",
            )
        return True

    def run(self, sample: PipelineJob, **kwargs) -> PipelineJob:
        reads_in = Path(sample.input_reads or "")
        ensure_dir(sample.reads_dir)
        ensure_dir(sample.log_path)

        out_final = sample.reads_dir / "reads.trim.fq.gz"
        pre_cfg = self.config.preprocess
        num_threads = sample.num_threads

        self.logger.info("Preprocessing input: %s", reads_in)

        # --- Step 1: adapter trimming ---
        reads_in = self._run_adapter_trimming(
            reads_in=reads_in,
            reads_out=sample.reads_dir / "reads.adaptertrim.tmp.fq.gz",
            cfg=pre_cfg.porechop_abi,
            num_threads=num_threads,
        )

        # --- Step 2: quality/length filtering ---
        reads_in = self._run_chopper(
            reads_in=reads_in,
            reads_out=out_final,
            cfg=pre_cfg.chopper,
            num_threads=num_threads,
        )
        if reads_in != out_final:
            # neither chopper ran, just move/copy to final destination
            if reads_in.suffix == out_final.suffix:
                reads_in.rename(out_final)
            else:
                copy_fastq_gz(reads_in, out_final)

        # clean up intermediate file if it still exists
        tmp = sample.reads_dir / "reads.adaptertrim.tmp.fq.gz"
        if tmp.exists() and tmp != out_final:
            tmp.unlink()

        self.logger.info("Trimmed reads -> %s", out_final)
        sample.set_file("trimmed_reads", out_final)
        sample.set_file("current_reads", out_final)
        return sample

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _run_adapter_trimming(
        self,
        reads_in: Path,
        reads_out: Path,
        cfg: PorechopABIConfig,
        num_threads: int,
    ) -> Path:
        """Run adapter trimming; try porechop_abi, fall back to porechop."""
        # Try porechop_abi first
        try:
            wrapper = PorechopABI(num_threads=num_threads)
            self.logger.info(
                "Running Porechop_ABI: %s -> %s", reads_in, reads_out
            )
            wrapper.trim_from_config(
                input_fq=reads_in,
                output_fq=reads_out,
                cfg=cfg,
            )
            return reads_out
        except DependencyError:
            self.logger.warning(
                "porechop_abi not found — falling back to porechop "
                "(ab initio options will be skipped)"
            )

        # Fall back to original porechop with compatible options
        try:
            wrapper_legacy = Porechop(num_threads=num_threads)
            self.logger.info(
                "Running Porechop (fallback): %s -> %s", reads_in, reads_out
            )
            wrapper_legacy.trim_from_config(
                input_fq=reads_in,
                output_fq=reads_out,
                cfg=cfg,
            )
            return reads_out
        except DependencyError:
            self.logger.warning(
                "porechop not found either — skipping adapter trimming"
            )
            return reads_in  # pass through unchanged

    def _run_chopper(
        self,
        reads_in: Path,
        reads_out: Path,
        cfg: ChopperConfig,
        num_threads: int,
    ) -> Path:
        self.logger.info("Running Chopper: %s -> %s", reads_in, reads_out)
        wrapper = Chopper(num_threads=num_threads)
        wrapper.filter_from_config(
            input_fq=reads_in,
            output_fq=reads_out,
            cfg=cfg,
        )
        return reads_out
