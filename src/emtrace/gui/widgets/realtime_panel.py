"""Real-time monitoring panel for the EMtrace desktop GUI.

Displays watch status, batch progress, throughput metrics, and a live
top-taxa list driven by :class:`~emtrace.gui.workers.realtime_worker.RealtimeWorker`
signals.  The panel owns its Start / Stop controls and emits
:pyqtSignal:`start_requested` / :pyqtSignal:`stop_requested` for the
:class:`~emtrace.gui.main_window.MainWindow` to act on.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import (
    QCheckBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)


class RealtimePanelWidget(QWidget):
    """Panel combining watch controls, live status, and top-taxa display."""

    #: Emitted when the user clicks "Start Watching".
    #: Arguments: watch_dir, batch_size, debounce_seconds, scan_existing
    start_requested = pyqtSignal(Path, int, float, bool)

    #: Emitted when the user clicks "Stop".
    stop_requested = pyqtSignal()

    _STATUS_STYLE: dict[str, str] = {
        "watching":   "color: green;  font-size: 16px;",
        "starting":   "color: orange; font-size: 16px;",
        "stopping":   "color: orange; font-size: 16px;",
        "finalizing": "color: orange; font-size: 16px;",
        "stopped":    "color: grey;   font-size: 16px;",
        "idle":       "color: grey;   font-size: 16px;",
        "error":      "color: red;    font-size: 16px;",
    }

    _STATUS_TEXT: dict[str, str] = {
        "starting":   "Starting…",
        "stopping":   "Stopping: finishing the current batch…",
        "finalizing": "Finalizing: merging batches and writing the report…",
    }

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._job_dir: Optional[Path] = None
        self._build_ui()
        self.set_job(None)
        self._set_controls("idle")

    # ==================================================================
    # UI construction
    # ==================================================================

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(8)

        # ── Control form ─────────────────────────────────────────────
        ctrl_form = QFormLayout()
        ctrl_form.setContentsMargins(0, 0, 0, 0)

        self._watch_edit = QLineEdit()
        self._watch_edit.setPlaceholderText("Directory to monitor for FASTQ files…")
        watch_browse = QPushButton("Browse…")
        watch_browse.clicked.connect(self._pick_watch_dir)
        ctrl_form.addRow("Watch directory:", self._make_browse_row(self._watch_edit, watch_browse))

        self._results_label = QLabel()
        self._results_label.setWordWrap(True)
        ctrl_form.addRow("Results:", self._results_label)

        self._scan_check = QCheckBox("Process files already in the directory")
        self._scan_check.setChecked(True)
        self._scan_check.setToolTip(
            "If unchecked, only files that arrive after Start are processed. "
            "A job watched before always catches up on files in the directory."
        )
        ctrl_form.addRow("", self._scan_check)

        self._batch_spin = QSpinBox()
        self._batch_spin.setRange(1, 1000)
        self._batch_spin.setValue(1)
        self._batch_spin.setToolTip("Accumulate this many FASTQ files before processing a batch")

        self._debounce_spin = QDoubleSpinBox()
        self._debounce_spin.setRange(0.5, 60.0)
        self._debounce_spin.setSingleStep(0.5)
        self._debounce_spin.setValue(2.0)
        self._debounce_spin.setSuffix(" s")
        self._debounce_spin.setToolTip(
            "Wait this many seconds after a file stops growing before treating it as complete"
        )

        params_widget = QWidget()
        params_layout = QHBoxLayout(params_widget)
        params_layout.setContentsMargins(0, 0, 0, 0)
        params_layout.addWidget(QLabel("Batch size:"))
        params_layout.addWidget(self._batch_spin)
        params_layout.addSpacing(20)
        params_layout.addWidget(QLabel("Debounce:"))
        params_layout.addWidget(self._debounce_spin)
        params_layout.addStretch()
        ctrl_form.addRow("", params_widget)

        layout.addLayout(ctrl_form)

        # ── Start / Stop buttons ──────────────────────────────────────
        btn_layout = QHBoxLayout()
        btn_layout.setContentsMargins(0, 4, 0, 4)
        self._start_btn = QPushButton("▶  Start Watching")
        self._stop_btn = QPushButton("■  Stop")
        self._start_btn.clicked.connect(self._on_start_clicked)
        self._stop_btn.clicked.connect(self._on_stop_clicked)
        btn_layout.addWidget(self._start_btn)
        btn_layout.addWidget(self._stop_btn)
        btn_layout.addStretch()
        layout.addLayout(btn_layout)

        layout.addWidget(self._separator())

        # ── Status row ───────────────────────────────────────────────
        status_layout = QHBoxLayout()
        status_layout.setContentsMargins(0, 0, 0, 0)
        self._indicator = QLabel("●")
        self._indicator.setStyleSheet(self._STATUS_STYLE["idle"])
        self._status_text = QLabel("Idle")
        status_layout.addWidget(self._indicator)
        status_layout.addWidget(self._status_text)
        status_layout.addStretch()
        layout.addLayout(status_layout)

        # ── Stats form ───────────────────────────────────────────────
        stats_form = QFormLayout()
        stats_form.setContentsMargins(0, 0, 0, 0)
        self._batch_label = QLabel("—")
        self._reads_label = QLabel("—")
        self._files_label = QLabel("—")
        self._taxa_count_label = QLabel("—")
        stats_form.addRow("Batch:", self._batch_label)
        stats_form.addRow("Total reads:", self._reads_label)
        stats_form.addRow("Total files:", self._files_label)
        stats_form.addRow("Taxa detected:", self._taxa_count_label)
        layout.addLayout(stats_form)

        layout.addWidget(self._separator())

        # ── Top taxa list ────────────────────────────────────────────
        layout.addWidget(QLabel("Top taxa (by EM abundance):"))
        self._taxa_list = QListWidget()
        self._taxa_list.setMaximumHeight(220)
        layout.addWidget(self._taxa_list)

        layout.addStretch()

    @staticmethod
    def _make_browse_row(line_edit: QLineEdit, browse_btn: QPushButton) -> QWidget:
        w = QWidget()
        h = QHBoxLayout(w)
        h.setContentsMargins(0, 0, 0, 0)
        h.addWidget(line_edit, stretch=1)
        h.addWidget(browse_btn)
        return w

    @staticmethod
    def _separator() -> QFrame:
        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        return sep

    # ==================================================================
    # Slots — connected from RealtimeWorker signals
    # ==================================================================

    def on_batch_started(self, batch_num: int) -> None:
        """Fired before a batch is processed."""
        self._batch_label.setText(f"{batch_num}  (processing…)")

    def on_batch_finished(self, batch_num: int, reads_in_batch: int, total_reads: int) -> None:
        """Fired after a batch completes."""
        self._batch_label.setText(f"{batch_num}  (+{reads_in_batch:,} reads this batch)")
        self._reads_label.setText(f"{total_reads:,}")

    def on_abundance_updated(self, snapshot: dict) -> None:
        """Refresh stats and top-taxa list from a RealtimeState snapshot."""
        self._reads_label.setText(f"{snapshot.get('total_reads', 0):,}")
        self._files_label.setText(str(snapshot.get("total_files", "—")))
        self._taxa_count_label.setText(str(snapshot.get("n_taxa", "—")))

        top = snapshot.get("top_taxa", [])
        names = snapshot.get("taxa_names", {})
        self._taxa_list.clear()
        for taxid, prob in top[:10]:
            pct = prob * 100
            item = QListWidgetItem(f"{names.get(taxid) or taxid}   {pct:.1f} %")
            item.setToolTip(f"taxid {taxid}")
            self._taxa_list.addItem(item)

    def on_watch_status_changed(self, status: str) -> None:
        """Update the status indicator and toggle control availability."""
        style = self._STATUS_STYLE.get(status, self._STATUS_STYLE["idle"])
        self._indicator.setStyleSheet(style)
        self._status_text.setText(self._STATUS_TEXT.get(status, status.capitalize()))
        self._set_controls(status)

    def on_error(self, msg: str) -> None:
        QMessageBox.warning(self, "Realtime error", msg)

    # ==================================================================
    # Called by MainWindow
    # ==================================================================

    def set_job(self, job) -> None:
        """Show where results go; clear the live numbers when the job changes."""
        job_dir = job.output_dir if job is not None else None
        if job_dir == self._job_dir and job_dir is not None:
            return
        self._job_dir = job_dir
        self._results_label.setText(
            str(job.watch_dir) if job is not None
            else "<job directory>/watch (a job is created when you start watching)"
        )
        for label in (self._batch_label, self._reads_label,
                      self._files_label, self._taxa_count_label):
            label.setText("—")
        self._taxa_list.clear()

    # ==================================================================
    # Internal helpers
    # ==================================================================

    def _set_controls(self, status: str) -> None:
        """Enable controls for *status*: only Stop while watching, none in transitions."""
        busy = status in ("starting", "watching", "stopping", "finalizing")
        self._start_btn.setEnabled(not busy)
        self._stop_btn.setEnabled(status == "watching")
        for w in (self._watch_edit, self._scan_check,
                  self._batch_spin, self._debounce_spin):
            w.setEnabled(not busy)

    def _on_stop_clicked(self) -> None:
        self._stop_btn.setEnabled(False)
        self.stop_requested.emit()

    def _pick_watch_dir(self) -> None:
        path = QFileDialog.getExistingDirectory(
            self, "Select Watch Directory", self._watch_edit.text() or "",
        )
        if path:
            self._watch_edit.setText(path)

    def _on_start_clicked(self) -> None:
        watch_str = self._watch_edit.text().strip()
        if not watch_str:
            QMessageBox.warning(self, "Missing directory", "Select a watch directory first.")
            return
        watch_dir = Path(watch_str)
        if not watch_dir.is_dir():
            QMessageBox.warning(self, "Invalid directory", f"Not a directory:\n{watch_dir}")
            return

        self.start_requested.emit(
            watch_dir,
            self._batch_spin.value(),
            self._debounce_spin.value(),
            self._scan_check.isChecked(),
        )
