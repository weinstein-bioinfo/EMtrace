"""Porechop / Porechop_ABI adapter-removal step.

Registered as :data:`StepName.PORE`.

Behaviour:
  - Tries ``porechop_abi`` first; falls back to ``porechop``; if neither is
    on PATH, logs a warning and forwards reads unchanged.
  - Log (stdout/stderr from the tool) is saved to ``<logs>/porechop.log``.
  - Adapter trim counts (start/end/middle) are emitted at INFO.

Output file key: ``"porechop_reads"`` / ``"current_reads"``
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import TYPE_CHECKING, Optional

from emtrace.core.config import PorechopABIConfig
from emtrace.core.exceptions import DependencyError, StepValidationError
from emtrace.steps import StepName, register_step
from emtrace.steps.base import PipelineStep, resolve_preprocess_input
from emtrace.utils.paths import ensure_dir
from emtrace.wrappers.porechop import Porechop, PorechopABI

if TYPE_CHECKING:
    from emtrace.core.pipelinejob import PipelineJob


@register_step(StepName.PORE)
class PorechopStep(PipelineStep):
    """Adapter removal using Porechop_ABI (fallback: Porechop)."""

    step_name = "porechop"
    step_description = "Adapter removal (Porechop_ABI with fallback to Porechop)"

    def validate_inputs(self, sample: PipelineJob, **kwargs) -> bool:
        if sample.input_reads is None or not Path(sample.input_reads).exists():
            raise StepValidationError(self.step_name, "No input reads available")
        return True

    def run(self, sample: PipelineJob, **kwargs) -> PipelineJob:
        ensure_dir(sample.reads_dir)
        ensure_dir(sample.log_path)

        reads_in = resolve_preprocess_input(sample, self.config, self.logger)

        cfg = self.config.preprocess.porechop_abi

        reads_out = sample.reads_dir / "reads.porechop.fq.gz"
        log_file = sample.log_path / "porechop.log"
        adapters_file = sample.log_path / "porechop_adapters.txt"

        reads_out = self._run_adapter_trimming(
            reads_in=reads_in,
            reads_out=reads_out,
            cfg=cfg,
            num_threads=sample.num_threads,
            log_file=log_file,
            adapters_file=adapters_file,
        )

        self.logger.info("Adapter-trimmed reads -> %s", reads_out)
        sample.set_file("porechop_reads", reads_out)
        sample.set_file("current_reads", reads_out)
        return sample

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _run_adapter_trimming(
        self,
        reads_in: Path,
        reads_out: Path,
        cfg: PorechopABIConfig,
        num_threads: int,
        log_file: Path,
        adapters_file: Path,
    ) -> Path:
        # Try Porechop_ABI first
        try:
            wrapper = PorechopABI(num_threads=num_threads)
            self.logger.info(
                "Running Porechop_ABI: %s -> %s", reads_in, reads_out
            )
            wrapper.trim_from_config(
                input_fq=reads_in,
                output_fq=reads_out,
                cfg=cfg,
                log_file=log_file,
                adapters_file=adapters_file,
            )
            self._log_adapter_summary(log_file, adapters_file, tool="porechop_abi")
            return reads_out
        except DependencyError:
            self.logger.warning(
                "porechop_abi not found — falling back to porechop "
                "(ab initio options will be skipped)"
            )

        # Fall back to original Porechop
        try:
            wrapper_legacy = Porechop(num_threads=num_threads)
            self.logger.info(
                "Running Porechop (fallback): %s -> %s", reads_in, reads_out
            )
            wrapper_legacy.trim_from_config(
                input_fq=reads_in,
                output_fq=reads_out,
                cfg=cfg,
                log_file=log_file,
            )
            self._log_adapter_summary(log_file, adapters_file=None, tool="porechop")
            return reads_out
        except DependencyError:
            self.logger.warning(
                "porechop not found either — skipping adapter removal"
            )
            return reads_in  # pass through unchanged

    def _log_adapter_summary(
        self,
        log_file: Path,
        adapters_file: Optional[Path],
        tool: str,
    ) -> None:
        """Parse the tool log and emit trim statistics at INFO level."""
        if not log_file.exists():
            return
        log_text = log_file.read_text(errors="replace")
        stats = _parse_trim_stats(log_text)

        if stats:
            start = stats.get("start", "")
            end = stats.get("end", "")
            middle = stats.get("middle", "")
            if start:
                self.logger.info("%s: reads with start adapters trimmed: %s", tool, start)
            if end:
                self.logger.info("%s: reads with end adapters trimmed:   %s", tool, end)
            if middle:
                self.logger.info("%s: reads split on middle adapters:     %s", tool, middle)
        else:
            self.logger.info(
                "%s: no trim statistics found in log (see %s for details)",
                tool, log_file,
            )


# ---------------------------------------------------------------------------
# Log parser
# ---------------------------------------------------------------------------

def _parse_trim_stats(log_text: str) -> dict:
    """Parse adapter-trimming counts from Porechop / Porechop_ABI output.

    Looks for lines of the form (on stdout):
      651 / 933 reads had adapters trimmed from their start (36,193 bp removed)
      516 / 933 reads had adapters trimmed from their end (43,737 bp removed)
       21 / 933 reads were split based on middle adapters

    Returns a dict with keys "start", "end", "middle" whose values are
    human-readable strings like "651 / 933 (36,193 bp removed)".
    """
    stats: dict = {}

    # Start / end adapter trimming
    for m in re.finditer(
        r"([\d,]+)\s*/\s*([\d,]+)\s+reads had adapters trimmed from their\s+(start|end)"
        r"(?:\s+\(([\d,]+\s+bp\s+removed)\))?",
        log_text,
    ):
        side = m.group(3)
        trimmed, total = m.group(1), m.group(2)
        bp_info = m.group(4)  # may be None for older porechop without bp count
        if bp_info:
            stats[side] = f"{trimmed} / {total} ({bp_info})"
        else:
            stats[side] = f"{trimmed} / {total}"

    # Middle adapter splits
    m = re.search(
        r"([\d,]+)\s*/\s*([\d,]+)\s+reads were split based on middle adapters",
        log_text,
    )
    if m:
        stats["middle"] = f"{m.group(1)} / {m.group(2)}"

    return stats
