"""QThread wrapper for pipeline execution with Qt signal feedback.

Runs the chosen sequence of pipeline steps in a background thread so the
GUI stays responsive.  Each step boundary emits a signal the MainWindow
wires to :class:`StepPanel` and :class:`DashboardWidget` for live status
updates.

Cancellation has two modes:

* :meth:`PipelineWorker.request_cancel` — cooperative.  Sets a flag
  that is checked between steps; the active step finishes naturally.
* :meth:`PipelineWorker.force_quit` — immediate.  Sets the cancel flag
  *and* signals the running subprocess group via
  :func:`emtrace.utils.shell.kill_active`, so even a long-running
  external tool (annotation, assembly) is aborted within seconds.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional

from PyQt6.QtCore import QObject, QThread, pyqtSignal

from emtrace.core.config import EMtraceConfig
from emtrace.core.pipeline import Pipeline
from emtrace.core.pipelinejob import PipelineJob
from emtrace.steps import STEP_REGISTRY, StepName
from emtrace.utils import shell
from emtrace.utils.paths import working_dir

logger = logging.getLogger(__name__)


class WorkerSignals(QObject):
    """Signals emitted by :class:`PipelineWorker`."""

    step_started = pyqtSignal(str)            # step value (e.g. "pore")
    step_finished = pyqtSignal(str, bool)     # step value, success
    pipeline_finished = pyqtSignal(bool)      # overall success
    progress = pyqtSignal(int, int)           # current_step (1-based), total
    error = pyqtSignal(str)                   # human-readable error message
    job_updated = pyqtSignal(object)          # final PipelineJob after run


class PipelineWorker(QThread):
    """Run a list of pipeline steps in a background thread.

    Usage::

        worker = PipelineWorker(config, job, [StepName.PORE, StepName.CHOP, StepName.FLT])
        worker.signals.step_started.connect(panel.mark_running)
        worker.signals.step_finished.connect(panel.mark_done)
        worker.signals.pipeline_finished.connect(self._on_done)
        worker.signals.error.connect(self._on_error)
        worker.start()
    """

    def __init__(
        self,
        config: EMtraceConfig,
        job: PipelineJob,
        steps: List[StepName],
        kwargs: Optional[Dict] = None,
        parent: Optional[QObject] = None,
    ) -> None:
        super().__init__(parent)
        self.signals = WorkerSignals()
        self._config = config
        self._job = job
        self._steps = list(steps)
        self._kwargs = kwargs or {}
        self._cancel_requested = False

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    @property
    def job(self) -> PipelineJob:
        """Return the (possibly mutated) job after execution."""
        return self._job

    def request_cancel(self) -> None:
        """Request cancellation; the worker stops after the current step."""
        self._cancel_requested = True
        logger.info("Pipeline cancellation requested")

    def force_quit(self) -> bool:
        """Set the cancel flag and signal the running subprocess group.

        Returns ``True`` if a live subprocess was signalled, ``False`` if
        no subprocess was active (e.g. between steps).  The worker's
        ``run()`` loop will exit naturally once the killed subprocess
        raises and the next-step check sees ``_cancel_requested``.
        """
        self._cancel_requested = True
        logger.warning("Pipeline force-quit requested")
        return shell.kill_active()

    # ------------------------------------------------------------------
    # QThread entry point
    # ------------------------------------------------------------------

    def run(self) -> None:
        """Execute the configured step sequence (called by ``QThread.start``).

        Runs inside the job directory so external tools never inherit the
        read-only ``/`` working directory of a Finder-launched app.
        """
        with working_dir(self._job.output_dir):
            self._run()

    def _run(self) -> None:
        pipeline = Pipeline(self._config)
        pipeline.setup_job(self._job)
        total = len(self._steps)
        success = True

        try:
            for index, step_name in enumerate(self._steps, start=1):
                if self._cancel_requested:
                    self.signals.error.emit("Pipeline cancelled by user.")
                    success = False
                    break

                if step_name not in STEP_REGISTRY:
                    msg = f"Step '{step_name.value}' is not registered"
                    logger.warning(msg)
                    self.signals.error.emit(msg)
                    self.signals.step_finished.emit(step_name.value, False)
                    success = False
                    break

                self.signals.progress.emit(index, total)
                self.signals.step_started.emit(step_name.value)

                try:
                    step = pipeline._steps[step_name]
                    self._job = step.execute(self._job, **self._kwargs)
                    self.signals.step_finished.emit(step_name.value, True)
                except Exception as exc:  # noqa: BLE001
                    if self._cancel_requested:
                        # Killed intentionally — no traceback spam.
                        logger.info("Step %s aborted: %s", step_name.value, exc)
                    else:
                        logger.exception("Step %s failed", step_name.value)
                    self.signals.step_finished.emit(step_name.value, False)
                    self.signals.error.emit(f"{step_name.value}: {exc}")
                    success = False
                    break
        finally:
            self.signals.job_updated.emit(self._job)
            self.signals.pipeline_finished.emit(success)
