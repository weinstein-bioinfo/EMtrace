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
from emtrace.realtime.processor import BatchProcessor
from emtrace.realtime.watcher import FastqWatcher

logger = logging.getLogger(__name__)


class RealtimeWorkerSignals(QObject):
    """Qt signals emitted by :class:`RealtimeWorker`."""

    batch_started = pyqtSignal(int)           # batch_num (1-based)
    batch_finished = pyqtSignal(int, int, int) # batch_num, reads_in_batch, total_reads
    abundance_updated = pyqtSignal(dict)       # RealtimeState.get_snapshot()
    error = pyqtSignal(str)
    watch_status_changed = pyqtSignal(str)     # "watching" | "stopped"
    worker_finished = pyqtSignal()


class RealtimeWorker(QThread):
    """Background thread that watches a directory and processes FASTQ batches.

    Usage::

        worker = RealtimeWorker(config, watch_dir="/data/nanopore", output_dir="/out")
        worker.signals.batch_finished.connect(panel.on_batch_finished)
        worker.signals.abundance_updated.connect(panel.on_abundance_updated)
        worker.signals.watch_status_changed.connect(panel.on_watch_status_changed)
        worker.signals.worker_finished.connect(self._on_watch_done)
        worker.start()
        ...
        worker.request_stop()
    """

    def __init__(
        self,
        config: EMtraceConfig,
        watch_dir: Path,
        output_dir: Path,
        job_name: str = "realtime",
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
        self._output_dir = Path(output_dir)
        self._job_name = job_name
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
        logger.info("Realtime watch stop requested")

    def run(self) -> None:
        """QThread entry point — blocks until stop_event is set."""
        try:
            watcher = FastqWatcher(
                self._watch_dir,
                debounce_seconds=self._debounce,
                scan_existing=self._scan_existing,
            )
            processor = BatchProcessor(
                config=self._config,
                output_dir=self._output_dir,
                job_name=self._job_name,
                ann1=self._ann1,
                ann2=self._ann2,
                anns=self._anns,
                hostgen=self._hostgen,
                batch_size=self._batch_size,
                run_report=self._run_report,
            )
        except Exception as exc:
            logger.error("Realtime setup failed: %s", exc)
            self.signals.error.emit(f"Setup failed: {exc}")
            self.signals.worker_finished.emit()
            return

        watcher.start()
        self.signals.watch_status_changed.emit("watching")
        batch_num = 0
        prev_total = 0
        batch_buffer = []

        try:
            while not self._stop_event.is_set():
                path = watcher.get(timeout=1.0)
                if path is None:
                    continue

                batch_buffer.append(path)

                if len(batch_buffer) >= self._batch_size:
                    batch_num += 1
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
            if batch_buffer:
                logger.info("Processing %d remaining file(s) before shutdown", len(batch_buffer))
                try:
                    processor.process_files(batch_buffer)
                except Exception as exc:
                    logger.error("Final batch failed: %s", exc)
            try:
                processor.generate_final_report()
            except Exception as exc:
                logger.error("Final report failed: %s", exc)

            watcher.stop()
            self.signals.watch_status_changed.emit("stopped")
            self.signals.worker_finished.emit()
