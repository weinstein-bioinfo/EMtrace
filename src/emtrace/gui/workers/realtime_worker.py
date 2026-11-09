"""QThread wrapper for real-time FASTQ directory watching and batch processing.

Drives :class:`~emtrace.realtime.watcher.FastqWatcher` and
:class:`~emtrace.realtime.processor.BatchProcessor` in a background thread,
emitting Qt signals at each batch boundary so the GUI stays responsive.

The main loop mirrors :func:`~emtrace.realtime.processor.run_realtime` but
is inlined here to allow ``batch_started`` to fire *before* processing and
``batch_finished`` + ``abundance_updated`` to fire *after*.
"""

from __future__ import annotations

import logging
import threading
from pathlib import Path
from typing import Optional

from PyQt6.QtCore import QObject, QThread, pyqtSignal

from emtrace.core.config import EMtraceConfig
from emtrace.core.pipelinejob import PipelineJob
from emtrace.realtime.processor import BatchProcessor
from emtrace.realtime.watcher import FastqWatcher
from emtrace.utils.paths import working_dir

logger = logging.getLogger(__name__)


class RealtimeWorkerSignals(QObject):
    """Qt signals emitted by :class:`RealtimeWorker`."""

    batch_started = pyqtSignal(int)           # batch_num (1-based)
    batch_finished = pyqtSignal(int, int, int) # batch_num, reads_in_batch, total_reads
    abundance_updated = pyqtSignal(dict)       # RealtimeState.get_snapshot()
    error = pyqtSignal(str)
    # "watching" | "stopping" | "finalizing" | "stopped" | NO_REPORT_STATUS
    watch_status_changed = pyqtSignal(str)
    job_updated = pyqtSignal(object)           # PipelineJob after the final merge


#: Final status when no batch completed, so there is nothing to report.
NO_REPORT_STATUS = "stopped: no batch completed, so no report"


class RealtimeWorker(QThread):
    """Background thread that watches a directory and processes FASTQ batches.

    Results go to the job: batches to ``<job>/watch``, merged into the job's
    standard layout when stopped. A job watched before resumes where it ended.

    Usage::

        worker = RealtimeWorker(config, watch_dir="/data/nanopore", job=job)
        worker.signals.batch_finished.connect(panel.on_batch_finished)
        worker.signals.abundance_updated.connect(panel.on_abundance_updated)
        worker.signals.watch_status_changed.connect(panel.on_watch_status_changed)
        worker.signals.job_updated.connect(self._on_job_updated)
        worker.finished.connect(self._on_watch_done)
        worker.start()
        ...
        worker.request_stop()
    """

    def __init__(
        self,
        config: EMtraceConfig,
        watch_dir: Path,
        job: PipelineJob,
        *,
        ann1: str = "",
        ann2: str = "",
        anns: str = "",
        hostgen: str = "",
        batch_size: int = 1,
        debounce: float = 2.0,
        scan_existing: bool = False,
        run_report: bool = False,
        parent: Optional[QObject] = None,
    ) -> None:
        super().__init__(parent)
        self.signals = RealtimeWorkerSignals()
        self._config = config
        self._watch_dir = Path(watch_dir)
        self._job = job
        self._ann1 = ann1
        self._ann2 = ann2
        self._anns = anns
        self._hostgen = hostgen
        self._batch_size = batch_size
        self._debounce = debounce
        self._scan_existing = scan_existing
        self._run_report = run_report
        self._stop_event = threading.Event()

    def request_stop(self) -> None:
        """Signal the watch loop to exit after the current batch completes."""
        self._stop_event.set()
        self.signals.watch_status_changed.emit("stopping")
        logger.info("Realtime watch stop requested")

    def run(self) -> None:
        """QThread entry point — blocks until stop_event is set.

        Runs inside the job directory so external tools never inherit the
        read-only ``/`` working directory of a Finder-launched app.
        """
        with working_dir(self._job.output_dir):
            self._run()

    def _run(self) -> None:
        try:
            processor = BatchProcessor(
                config=self._config,
                job_dir=self._job.output_dir,
                job_name=self._job.name,
                ann1=self._ann1,
                ann2=self._ann2,
                anns=self._anns,
                hostgen=self._hostgen,
                batch_size=self._batch_size,
                run_report=self._run_report,
            )
            processed = processor.resume()
            # A resumed watch catches up on files that arrived while stopped
            watcher = FastqWatcher(
                self._watch_dir,
                debounce_seconds=self._debounce,
                scan_existing=self._scan_existing or bool(processed),
                already_processed=processed,
            )
            watcher.start()
        except Exception as exc:
            logger.error("Realtime setup failed: %s", exc)
            self.signals.error.emit(f"Setup failed: {exc}")
            self.signals.watch_status_changed.emit("stopped")
            return

        self.signals.watch_status_changed.emit("watching")
        snap = processor.state.get_snapshot()
        prev_total = snap["total_reads"]
        if processed:
            self.signals.abundance_updated.emit(snap)
        batch_buffer = []

        try:
            while not self._stop_event.is_set():
                path = watcher.get(timeout=1.0)
                if path is None:
                    continue

                batch_buffer.append(path)

                if len(batch_buffer) >= self._batch_size:
                    batch_num = processor.next_batch
                    self.signals.batch_started.emit(batch_num)
                    try:
                        processor.process_files(batch_buffer)
                    except Exception as exc:
                        logger.error("Batch %d failed: %s", batch_num, exc)
                        self.signals.error.emit(f"Batch {batch_num}: {exc}")
                    batch_buffer.clear()

                    snap = processor.state.get_snapshot()
                    reads_this = snap["total_reads"] - prev_total
                    prev_total = snap["total_reads"]
                    self.signals.batch_finished.emit(batch_num, reads_this, snap["total_reads"])
                    self.signals.abundance_updated.emit(snap)

        except Exception as exc:
            logger.exception("Realtime worker unexpected error")
            self.signals.error.emit(str(exc))
        finally:
            watcher.stop()
            status = self._finish(processor, watcher, batch_buffer)
            self.signals.job_updated.emit(processor.job)
            self.signals.watch_status_changed.emit(status)

    def _finish(
        self, processor: BatchProcessor, watcher: FastqWatcher, batch_buffer: list,
    ) -> str:
        """Process buffered files and finalize the job; returns the final status."""
        if batch_buffer:
            logger.info("Processing %d remaining file(s) before shutdown", len(batch_buffer))
            try:
                processor.process_files(batch_buffer)
            except Exception as exc:
                logger.error("Final batch failed: %s", exc)
                self.signals.error.emit(f"Final batch: {exc}")
            self.signals.abundance_updated.emit(processor.state.get_snapshot())
        if watcher.queued_count:
            logger.info(
                "%d queued file(s) not processed; starting this job again "
                "picks them up", watcher.queued_count,
            )

        self.signals.watch_status_changed.emit("finalizing")
        try:
            return "stopped" if processor.finalize() else NO_REPORT_STATUS
        except Exception as exc:
            logger.exception("Final report failed")
            self.signals.error.emit(f"Final report failed: {exc}")
            return "stopped"
