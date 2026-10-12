"""QApplication + :class:`MainWindow` — the root container of the EMtrace GUI.

Layout
------
    +-------------------------------------------------------------+
    |  toolbar (New / Open / Save / Run / Stop)                   |
    +--------------------+----------------------------------------+
    |  StepPanel         |  MainContentStack                      |
    |  (left, fixed)     |  Dashboard | Config | Results | Files  |
    |                    |                                        |
    +--------------------+----------------------------------------+
    |  LogPanel (bottom, resizable)                               |
    +-------------------------------------------------------------+
    |  status bar                                                 |
    +-------------------------------------------------------------+

Adapted from the ViRefSet ``app.py`` pattern.  EMtrace has no
``PipelineStatus`` tracker yet, so step state is reflected directly via
the :class:`PipelineWorker` signals (``step_started`` / ``step_finished``)
and the dashboard's in-memory failure set.
"""

from __future__ import annotations

import logging
import re
import sys
from pathlib import Path
from typing import Optional

from PyQt6.QtCore import Qt
from PyQt6.QtGui import (
    QBrush,
    QColor,
    QFont,
    QKeySequence,
    QShortcut,
    QTextCharFormat,
)
from PyQt6.QtWidgets import (
    QApplication,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMenu,
    QMessageBox,
    QPushButton,
    QSplitter,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from emtrace import __version__
from emtrace.core.config import EMtraceConfig
from emtrace.core.pipelinejob import PipelineJob
from emtrace.gui.models.log_handler import QtLogHandler
from emtrace.gui.platform import (
    add_recent,
    bundle_root,
    configure_platform,
    load_recent,
)
from emtrace.gui.widgets.log_panel import LogPanel
from emtrace.gui.widgets.main_content import MainContentStack
from emtrace.gui.widgets.run_toolbar import RunToolbar
from emtrace.gui.widgets.step_panel import StepPanel
from emtrace.gui.workers.pipeline_worker import PipelineWorker
from emtrace.gui.workers.realtime_worker import RealtimeWorker
from emtrace.steps import STEP_REGISTRY, StepName


# ---------------------------------------------------------------------------
# ANSI escape-code renderer for log viewer
# ---------------------------------------------------------------------------

# Standard 8 + bright-8 foreground colour palette (codes 30-37, 90-97).
# Chosen to be readable on both light and dark backgrounds.
_ANSI_FG: dict[int, QColor] = {
    30: QColor("#3c3c3c"),   # black  → dark gray (readable on light bg)
    31: QColor("#cc0000"),   # red
    32: QColor("#4e9a06"),   # green
    33: QColor("#9a6e00"),   # yellow → darkened for readability
    34: QColor("#3465a4"),   # blue
    35: QColor("#75507b"),   # magenta
    36: QColor("#06989a"),   # cyan
    37: QColor("#888888"),   # white  → mid-gray
    90: QColor("#555555"),   # bright black
    91: QColor("#ef2929"),   # bright red
    92: QColor("#8ae234"),   # bright green
    93: QColor("#c4a000"),   # bright yellow
    94: QColor("#729fcf"),   # bright blue
    95: QColor("#ad7fa8"),   # bright magenta
    96: QColor("#34e2e2"),   # bright cyan
    97: QColor("#d3d7cf"),   # bright white
}

_ANSI_RE = re.compile(r"\x1b\[([0-9;]*)m")


def _render_ansi(text_edit: QTextEdit, text: str) -> None:
    """Insert *text* into *text_edit* with ANSI colour/style codes rendered.

    Unknown codes are silently ignored; the current format carries forward
    unchanged, so malformed sequences degrade gracefully.
    """
    cursor = text_edit.textCursor()
    cursor.movePosition(cursor.MoveOperation.End)

    fmt = QTextCharFormat()
    pos = 0

    for m in _ANSI_RE.finditer(text):
        # Plain text segment before this escape
        if m.start() > pos:
            cursor.insertText(text[pos : m.start()], fmt)

        # Parse semicolon-separated code list; bare ESC[m means reset (code 0)
        raw = m.group(1)
        codes = [int(c) for c in raw.split(";") if c.isdigit()] if raw else [0]

        for code in codes:
            if code == 0:
                fmt = QTextCharFormat()                        # reset all
            elif code == 1:
                fmt.setFontWeight(QFont.Weight.Bold)
            elif code == 22:
                fmt.setFontWeight(QFont.Weight.Normal)
            elif code == 4:
                fmt.setFontUnderline(True)
            elif code == 24:
                fmt.setFontUnderline(False)
            elif code in _ANSI_FG:
                fmt.setForeground(QBrush(_ANSI_FG[code]))
            elif code == 39:                                   # default fg
                fmt.setForeground(QBrush())

        pos = m.end()

    # Remaining text after the last escape sequence
    if pos < len(text):
        cursor.insertText(text[pos:], fmt)


# ---------------------------------------------------------------------------
# Log file viewer dialog
# ---------------------------------------------------------------------------

class _LogViewerDialog(QDialog):
    """Non-modal dialog that displays a log file with ANSI colour rendering."""

    def __init__(self, path: Path, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle(path.name)
        self.resize(820, 560)
        self.setWindowModality(Qt.WindowModality.NonModal)

        from emtrace.gui.platform import monospace_font

        self._text = QTextEdit()
        self._text.setReadOnly(True)
        self._text.setFont(monospace_font(10))

        try:
            content = path.read_text(encoding="utf-8", errors="replace")
        except OSError as exc:
            content = f"Could not read file:\n{exc}"

        _render_ansi(self._text, content)
        # Scroll to end — most recent log lines are at the bottom
        self._text.moveCursor(self._text.textCursor().MoveOperation.End)

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.accept)

        layout = QVBoxLayout(self)
        layout.addWidget(QLabel(str(path)))
        layout.addWidget(self._text, stretch=1)
        btn_row = QHBoxLayout()
        btn_row.addStretch()
        btn_row.addWidget(close_btn)
        layout.addLayout(btn_row)


# ---------------------------------------------------------------------------
# New-sample dialog (small inline helper)
# ---------------------------------------------------------------------------

class _NewSampleDialog(QDialog):
    """Modal dialog collecting a sample name, input FASTQ, and output dir."""

    def __init__(
        self,
        default_output_dir: Optional[Path] = None,
        parent: Optional[QWidget] = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("New Job")
        self.setModal(True)
        self.resize(520, 0)

        self._name_edit = QLineEdit()
        self._name_edit.setPlaceholderText("e.g. sample01")

        self._reads_edit = QLineEdit()
        self._reads_edit.setPlaceholderText("Path to FASTQ (.fastq, .fastq.gz)")
        reads_browse = QPushButton("Browse\u2026")
        reads_browse.clicked.connect(self._pick_reads)
        reads_row = QHBoxLayout()
        reads_row.setContentsMargins(0, 0, 0, 0)
        reads_row.addWidget(self._reads_edit, stretch=1)
        reads_row.addWidget(reads_browse)
        reads_widget = QWidget()
        reads_widget.setLayout(reads_row)

        # Optional R2 mate (paired-end Illumina). When filled, the job is
        # treated as illumina_pe and Nanopore-only preprocessing steps are
        # auto-disabled.
        self._reads_r2_edit = QLineEdit()
        self._reads_r2_edit.setPlaceholderText(
            "Optional \u2014 R2 FASTQ for paired-end Illumina input"
        )
        reads_r2_browse = QPushButton("Browse\u2026")
        reads_r2_browse.clicked.connect(self._pick_reads_r2)
        reads_r2_row = QHBoxLayout()
        reads_r2_row.setContentsMargins(0, 0, 0, 0)
        reads_r2_row.addWidget(self._reads_r2_edit, stretch=1)
        reads_r2_row.addWidget(reads_r2_browse)
        reads_r2_widget = QWidget()
        reads_r2_widget.setLayout(reads_r2_row)

        self._out_edit = QLineEdit()
        if default_output_dir is not None:
            self._out_edit.setText(str(default_output_dir))
        out_browse = QPushButton("Browse\u2026")
        out_browse.clicked.connect(self._pick_outdir)
        out_row = QHBoxLayout()
        out_row.setContentsMargins(0, 0, 0, 0)
        out_row.addWidget(self._out_edit, stretch=1)
        out_row.addWidget(out_browse)
        out_widget = QWidget()
        out_widget.setLayout(out_row)

        form = QFormLayout()
        form.addRow("Job name:", self._name_edit)
        form.addRow("Input reads:", reads_widget)
        form.addRow("Reads R2 (paired-end):", reads_r2_widget)
        form.addRow("Output base directory:", out_widget)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok
            | QDialogButtonBox.StandardButton.Cancel,
            self,
        )
        buttons.accepted.connect(self._on_accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(buttons)

    def _pick_reads(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Input FASTQ",
            self._reads_edit.text() or "",
            "FASTQ files (*.fastq *.fq *.fastq.gz *.fq.gz);;All files (*)",
        )
        if path:
            self._reads_edit.setText(path)

    def _pick_reads_r2(self) -> None:
        start = self._reads_r2_edit.text() or self._reads_edit.text() or ""
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select R2 FASTQ (paired-end Illumina)",
            start,
            "FASTQ files (*.fastq *.fq *.fastq.gz *.fq.gz);;All files (*)",
        )
        if path:
            self._reads_r2_edit.setText(path)

    def _pick_outdir(self) -> None:
        path = QFileDialog.getExistingDirectory(
            self,
            "Select Output Directory",
            self._out_edit.text() or "",
        )
        if path:
            self._out_edit.setText(path)

    def _on_accept(self) -> None:
        name = self._name_edit.text().strip()
        reads = self._reads_edit.text().strip()
        reads_r2 = self._reads_r2_edit.text().strip()
        outdir = self._out_edit.text().strip()
        if not name:
            QMessageBox.warning(self, "Missing name", "Job name is required.")
            return
        if not reads or not Path(reads).is_file():
            QMessageBox.warning(self, "Invalid reads", "Pick an existing FASTQ file.")
            return
        if reads_r2 and not Path(reads_r2).is_file():
            QMessageBox.warning(self, "Invalid R2 reads", "Pick an existing R2 FASTQ or leave it blank.")
            return
        if not outdir:
            QMessageBox.warning(self, "Missing output", "Select an output directory.")
            return
        Path(outdir).mkdir(parents=True, exist_ok=True)
        self.accept()

    def result_values(self) -> tuple[str, Path, Path, Optional[Path]]:
        r2_text = self._reads_r2_edit.text().strip()
        return (
            self._name_edit.text().strip(),
            Path(self._reads_edit.text().strip()),
            Path(self._out_edit.text().strip()),
            Path(r2_text) if r2_text else None,
        )


# ---------------------------------------------------------------------------
# Main window
# ---------------------------------------------------------------------------

class MainWindow(QMainWindow):
    """Top-level EMtrace GUI window."""

    def __init__(
        self,
        config: EMtraceConfig,
        config_source: Optional[Path] = None,
    ) -> None:
        super().__init__()

        self._config = config
        self._config_source = config_source
        self._sample: Optional[PipelineJob] = None
        self._worker: Optional[PipelineWorker] = None
        self._realtime_worker: Optional[RealtimeWorker] = None
        self._closing: bool = False

        self.setWindowTitle("EMtrace")
        self.resize(1200, 800)

        self._build_menus()

        # --- Toolbar ---
        self._toolbar = RunToolbar(self)
        self.addToolBar(self._toolbar)

        # --- Centre: horizontal splitter (StepPanel | MainContentStack) ---
        self._step_panel = StepPanel(self)
        self._main_content = MainContentStack(self)

        hsplitter = QSplitter(Qt.Orientation.Horizontal)
        hsplitter.addWidget(self._step_panel)
        hsplitter.addWidget(self._main_content)
        hsplitter.setStretchFactor(0, 0)
        hsplitter.setStretchFactor(1, 1)
        hsplitter.setSizes([220, 980])

        # --- Bottom: log panel ---
        self._log_panel = LogPanel(self)

        vsplitter = QSplitter(Qt.Orientation.Vertical)
        vsplitter.addWidget(hsplitter)
        vsplitter.addWidget(self._log_panel)
        vsplitter.setStretchFactor(0, 1)
        vsplitter.setStretchFactor(1, 0)
        vsplitter.setSizes([580, 220])

        self.setCentralWidget(vsplitter)

        # --- Status bar ---
        self._status_label = QLabel("Ready")
        self.statusBar().addPermanentWidget(self._status_label)

        # --- Logging bridge (Python logging → LogPanel) ---
        self._log_handler = QtLogHandler()
        self._log_handler.signals.log_record.connect(self._log_panel.append_record)
        rt_logger = logging.getLogger("emtrace")
        rt_logger.setLevel(logging.DEBUG)
        rt_logger.addHandler(self._log_handler)

        # --- Toolbar wiring ---
        self._toolbar.new_clicked.connect(self._on_new)
        self._toolbar.open_clicked.connect(self._on_open)
        self._toolbar.save_clicked.connect(self._on_save)
        self._toolbar.run_clicked.connect(self._on_run)
        self._toolbar.stop_clicked.connect(self._on_stop)
        self._toolbar.force_quit_clicked.connect(self._on_force_quit)

        # Dashboard's empty-state "New Sample" button → same handler.
        self._main_content.dashboard.new_sample_requested.connect(self._on_new)

        # Dashboard "Clear files…" buttons per phase group.
        db = self._main_content.dashboard
        db.clear_preprocessing_requested.connect(self._on_clear_preprocessing)
        db.clear_annotation_requested.connect(self._on_clear_annotation)
        db.clear_em_requested.connect(self._on_clear_em)
        db.clear_assembly_requested.connect(self._on_clear_assembly)
        db.clear_reporting_requested.connect(self._on_clear_reporting)

        # Config editor Apply → update in-memory config.
        self._main_content.config_page.config_changed.connect(self._on_config_changed)

        # Step-panel double-click on status column → switch to Files tab.
        self._step_panel.view_requested.connect(self._on_view_step)

        # File browser double-click → route to appropriate viewer tab.
        self._main_content.file_browser.file_activated.connect(self._on_file_activated)

        # Realtime panel Start / Stop → worker lifecycle.
        rp = self._main_content.realtime_panel
        rp.start_requested.connect(self._on_watch_start)
        rp.stop_requested.connect(self._on_watch_stop)

        # --- Keyboard shortcuts: Ctrl+1..8 jump between tabs ---
        for i in range(8):
            QShortcut(
                QKeySequence(f"Ctrl+{i + 1}"), self,
                lambda idx=i: self._main_content.switch_to(idx),
            )
        # F1 opens the user guide
        QShortcut(QKeySequence("F1"), self, self._on_open_manual)

        # --- Disable steps that are declared in StepName but not registered ---
        # (e.g. ``chi`` has no implementation yet — checking it would only
        # produce a "step not registered" error when Run is clicked.)
        for step in StepName:
            if step not in STEP_REGISTRY:
                self._step_panel.set_step_enabled(step.value, False)

        # Pre-populate the Config editor with the startup config.
        self._main_content.set_config(self._config)

    # ==================================================================
    # Menu bar
    # ==================================================================

    def _build_menus(self) -> None:
        bar = self.menuBar()

        file_menu = bar.addMenu("&File")
        file_menu.addAction("&New Job\u2026", self._on_new)
        file_menu.addAction("&Open Job\u2026", self._on_open)
        self._recent_menu: QMenu = file_menu.addMenu("Open &Recent")
        self._recent_menu.aboutToShow.connect(self._populate_recent_menu)
        file_menu.addSeparator()
        file_menu.addAction("&Save Config\u2026", self._on_save)
        file_menu.addSeparator()
        file_menu.addAction("&Quit", self.close)

        view_menu = bar.addMenu("&View")
        M = MainContentStack
        view_menu.addAction("&Dashboard", lambda: self._main_content.switch_to(M.TAB_DASHBOARD))
        view_menu.addAction("&Config",    lambda: self._main_content.switch_to(M.TAB_CONFIG))
        view_menu.addAction("&Abundance", lambda: self._main_content.switch_to(M.TAB_ABUNDANCE))
        view_menu.addAction("&Assembly",  lambda: self._main_content.switch_to(M.TAB_ASSEMBLY))
        view_menu.addAction("Fi&gures",   lambda: self._main_content.switch_to(M.TAB_FIGURES))
        view_menu.addAction("Fi&les",     lambda: self._main_content.switch_to(M.TAB_FILES))
        view_menu.addAction("&QC",        lambda: self._main_content.switch_to(M.TAB_QC))
        view_menu.addAction("&Tables",    lambda: self._main_content.switch_to(M.TAB_TABLES))

        run_menu = bar.addMenu("&Run")
        run_menu.addAction("&Run Pipeline", self._on_run)
        run_menu.addAction("&Stop Pipeline", self._on_stop)
        run_menu.addAction("&Force Stop Pipeline", self._on_force_quit)

        tools_menu = bar.addMenu("&Tools")
        tools_menu.addAction("&Database Manager\u2026", self._on_database_manager)

        help_menu = bar.addMenu("&Help")
        help_menu.addAction("&User Guide\tF1", self._on_open_manual)
        help_menu.addSeparator()
        help_menu.addAction("&About", self._on_about)

    # ==================================================================
    # Sample / config slots
    # ==================================================================

    def _on_new(self) -> None:
        """File ▶ New Job — collect identity + paths and create a job."""
        recent = load_recent()
        default_out = (
            Path(recent[0]["path"]).parent if recent else Path.cwd()
        )

        dlg = _NewSampleDialog(default_output_dir=default_out, parent=self)
        if dlg.exec() != QDialog.DialogCode.Accepted:
            return

        name, reads_path, out_root, reads_r2 = dlg.result_values()
        sample_dir = out_root / name

        try:
            self._sample = PipelineJob(
                name=name,
                input_reads=reads_path,
                input_reads_r2=reads_r2,
                platform="illumina_pe" if reads_r2 is not None else "nanopore",
                output_dir=sample_dir,
                num_threads=self._config.runtime.num_threads,
            )
            self._sample.init_dirs()
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(
                self, "Error",
                f"Failed to create job:\n{exc}",
            )
            return

        self._activate_sample(self._sample)
        self._main_content.set_config(self._config)
        add_recent(str(sample_dir), name, input_reads=str(reads_path))
        logging.getLogger("emtrace").info("Created job: %s", sample_dir)

    def _on_open(self) -> None:
        """File ▶ Open Job — pick a results directory; auto-resolve input reads."""
        recent = load_recent()
        default_out = recent[0].get("path", "") if recent else ""

        outdir = QFileDialog.getExistingDirectory(
            self, "Select Job Directory", default_out,
        )
        if not outdir:
            return
        sample_dir = Path(outdir)

        reads_path = self._resolve_input_reads(sample_dir)
        if reads_path is None:
            return  # user cancelled the fallback prompt

        # R2 + platform are read directly from the manifest; if absent the job
        # is treated as single-end Nanopore (the historical default).
        manifest = PipelineJob.read_job_manifest(sample_dir) or {}
        r2_recorded = manifest.get("input_reads_r2") or ""
        reads_r2 = Path(r2_recorded) if r2_recorded and Path(r2_recorded).is_file() else None
        platform = manifest.get("platform") or ("illumina_pe" if reads_r2 else "nanopore")

        self._activate_job(reads_path, sample_dir, reads_r2=reads_r2, platform=platform)

    def _resolve_input_reads(self, sample_dir: Path) -> Optional[Path]:
        """Resolve the input FASTQ for *sample_dir* via manifest → recent → prompt.

        Returns the resolved path, or ``None`` if the user cancels the
        fallback file dialog.  Always re-prompts (and returns the user's
        pick) when the recorded path no longer exists on disk.
        """
        # 1. Job manifest written into the output dir
        manifest = PipelineJob.read_job_manifest(sample_dir)
        recorded = manifest.get("input_reads") if manifest else ""

        # 2. Fallback to the recent-jobs index
        if not recorded:
            sample_dir_s = str(sample_dir)
            for entry in load_recent():
                if entry.get("path") == sample_dir_s:
                    recorded = entry.get("input_reads", "")
                    break

        recorded_path = Path(recorded) if recorded else None
        if recorded_path and recorded_path.is_file():
            return recorded_path

        # 3. Re-prompt: nothing recorded, or recorded path is stale.
        title = "Select Input FASTQ"
        if recorded_path is not None:
            QMessageBox.warning(
                self, "Input FASTQ not found",
                f"The recorded input reads file no longer exists:\n\n{recorded_path}\n\n"
                "Pick the new location to relink it to this job.",
            )
        start_dir = str(recorded_path.parent) if recorded_path else ""
        path, _ = QFileDialog.getOpenFileName(
            self, title, start_dir,
            "FASTQ files (*.fastq *.fq *.fastq.gz *.fq.gz);;All files (*)",
        )
        if not path:
            return None
        return Path(path)

    def _activate_job(
        self,
        reads_path: Path,
        sample_dir: Path,
        reads_r2: Optional[Path] = None,
        platform: str = "nanopore",
    ) -> None:
        """Create a PipelineJob from *reads_path* + *sample_dir* and activate it.

        Loads an adjacent ``config.yaml`` if present, updates the UI,
        appends to the recent list, and (re)writes ``job.yaml`` so the
        manifest stays in sync with the resolved path.
        """
        try:
            self._sample = PipelineJob(
                name=sample_dir.name,
                input_reads=reads_path,
                input_reads_r2=reads_r2,
                platform=platform,
                output_dir=sample_dir,
                num_threads=self._config.runtime.num_threads,
            )
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, "Error", f"Failed to open job:\n{exc}")
            return

        sample_config = sample_dir / "config.yaml"
        if sample_config.exists():
            try:
                self._config = EMtraceConfig.from_yaml(sample_config)
                self._config_source = sample_config
            except Exception as exc:  # noqa: BLE001
                logging.getLogger("emtrace").warning(
                    "Could not load %s: %s", sample_config, exc,
                )

        # Refresh the manifest so legacy / re-linked / cold-open jobs gain one.
        self._sample._write_job_manifest()

        self._scan_existing_files(self._sample)
        self._activate_sample(self._sample)
        self._main_content.set_config(self._config)
        add_recent(str(sample_dir), sample_dir.name, input_reads=str(reads_path))
        logging.getLogger("emtrace").info("Opened job: %s", sample_dir)

    def _scan_existing_files(self, job: PipelineJob) -> None:
        """Populate job.files from files already on disk in the output directory.

        Called when opening an existing job so the QC panel and other
        viewers can find intermediate FASTQ files without a pipeline run.
        """
        rd = job.reads_dir

        # Input reads — expose in files dict so QC dropdown can find it
        if job.input_reads and job.input_reads.exists():
            job.set_file("input_reads", job.input_reads)
        if job.input_reads_r2 and job.input_reads_r2.exists():
            job.set_file("input_reads_r2", job.input_reads_r2)

        # Preprocessing outputs — register each known filename if present
        _reads_candidates = [
            ("reads.porechop.fq.gz", "porechop_reads"),
            ("reads.chopped.fq.gz",  "chopped_reads"),
        ]
        for fname, key in _reads_candidates:
            p = rd / fname
            if p.exists():
                job.set_file(key, p)

        # Host-filtered reads
        hflt = rd / "reads.trim.hflt.fq.gz"
        if hflt.exists():
            job.set_file("host_filtered_reads", hflt)

        # Per-host background filter outputs  hflt.{label}.reads.fq.gz
        for p in sorted(rd.glob("hflt.*.reads.fq.gz")):
            label = p.name[len("hflt."):-len(".reads.fq.gz")]
            job.set_file(f"hflt.{label}", p)

        # Annotation partition reads (plain .fq or .fq.gz)
        for rnd in ("ann1", "ann2"):
            for tag in ("vi", "ph", "ab", "euk", "un"):
                for suffix in (".fq.gz", ".fq"):
                    p = rd / f"reads.{rnd}.{tag}{suffix}"
                    if p.exists():
                        job.set_file(f"{rnd}_{tag}_reads", p)
                        break  # prefer .gz

        # current_reads: point to the most-processed reads file available,
        # mirroring the order steps would have set it during a live run.
        for key in ("host_filtered_reads", "chopped_reads", "porechop_reads"):
            if job.has_file(key):
                job.set_file("current_reads", job.files[key])
                break

        # Annotation outputs
        ann_dir = job.annotation_dir
        for key, fname in [("annot1", "annot1.tsv"), ("annot2", "annot2.tsv")]:
            p = ann_dir / fname
            if p.exists():
                job.set_file(key, p)

        # EM abundance outputs
        res_dir = job.reports_dir
        for key, fname in [
            ("em_taxid_prob",        "EM.taxid_prob.tsv"),
            ("em_taxid_prob_byiter", "EM.taxid_prob.byiter.tsv"),
        ]:
            p = res_dir / fname
            if p.exists():
                job.set_file(key, p)

        # Assembly outputs
        asm_dir = job.assembly_dir
        for key, fname in [
            ("assembly_consensus", "consensus.fa"),
            ("assembly_bam",       "reads.refgen.md.bam"),
        ]:
            p = asm_dir / fname
            if p.exists():
                job.set_file(key, p)

        # Plots output
        sunburst = res_dir / "plots" / "sunburst.html"
        if sunburst.exists():
            job.set_file("plots_sunburst", sunburst)

        # Reporting outputs
        reports_subdir = res_dir / "reports"
        for key, fname in [
            ("annot_table", "annot_table.tsv"),
            ("abund_table", "abund_table.tsv"),
            ("stats_html",  "stats_report.html"),
        ]:
            p = reports_subdir / fname
            if p.exists():
                job.set_file(key, p)

        # Pack output — archive lives directly in output_dir
        archive = job.output_dir / f"{job.name or job.reports_dir.name}.tar.gz"
        if archive.exists():
            job.set_file("pack_archive", archive)

    def _populate_recent_menu(self) -> None:
        """Rebuild the Open Recent submenu from the persisted recent list."""
        self._recent_menu.clear()
        entries = load_recent()
        if not entries:
            action = self._recent_menu.addAction("(no recent jobs)")
            action.setEnabled(False)
            return
        for entry in entries:
            name = entry.get("name", "Unknown")
            path = entry.get("path", "")
            input_reads = entry.get("input_reads", "")
            date = entry.get("last_opened", "")[:10]
            label = f"{name}  \u2014  {date}" if date else name
            action = self._recent_menu.addAction(label)
            action.setToolTip(path)
            action.triggered.connect(
                lambda _checked=False, p=path, ir=input_reads, n=name:
                    self._on_open_recent(p, ir, n)
            )
        self._recent_menu.addSeparator()
        self._recent_menu.addAction("Clear Recent", self._on_clear_recent)

    def _on_open_recent(self, path: str, input_reads: str, name: str) -> None:
        """Reopen a job from the recent list."""
        sample_dir = Path(path)
        if not sample_dir.exists():
            QMessageBox.warning(
                self, "Directory not found",
                f"Results directory no longer exists:\n{path}",
            )
            return

        reads_path = Path(input_reads) if input_reads else None
        if reads_path is None or not reads_path.is_file():
            # Reads path missing or stale — ask the user to locate it.
            reads_str, _ = QFileDialog.getOpenFileName(
                self,
                f"Locate Input Reads for \u2018{name}\u2019",
                str(sample_dir),
                "FASTQ files (*.fastq *.fq *.fastq.gz *.fq.gz);;All files (*)",
            )
            if not reads_str:
                return
            reads_path = Path(reads_str)

        manifest = PipelineJob.read_job_manifest(sample_dir) or {}
        r2_recorded = manifest.get("input_reads_r2") or ""
        reads_r2 = Path(r2_recorded) if r2_recorded and Path(r2_recorded).is_file() else None
        platform = manifest.get("platform") or ("illumina_pe" if reads_r2 else "nanopore")

        self._activate_job(reads_path, sample_dir, reads_r2=reads_r2, platform=platform)

    def _on_clear_recent(self) -> None:
        """Erase the persisted recent-jobs list."""
        from emtrace.gui.platform import recent_file
        try:
            recent_file().write_text("[]", encoding="utf-8")
        except OSError:
            pass

    def _activate_sample(self, sample: PipelineJob) -> None:
        """Push *sample* into every UI panel and switch to the Dashboard tab."""
        self._main_content.set_job(sample)
        self._step_panel.reset_status()
        self._step_panel.apply_platform(sample.platform)
        self._step_panel.update_quick_info(
            name=sample.name,
            base_dir=str(sample.output_dir),
            reads=None,
        )
        self._main_content.switch_to(MainContentStack.TAB_DASHBOARD)
        self._status_label.setText(f"Job: {sample.name}")

    def _on_save(self) -> None:
        """File ▶ Save Config — write the active config to a YAML file."""
        default_dir = (
            str(self._sample.output_dir) if self._sample else str(Path.cwd())
        )
        save_path, _ = QFileDialog.getSaveFileName(
            self, "Save Configuration",
            str(Path(default_dir) / "config.yaml"),
            "YAML files (*.yaml *.yml)",
        )
        if not save_path:
            return
        try:
            self._config.to_yaml(Path(save_path))
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(
                self, "Save failed", f"Could not write config:\n{exc}",
            )
            return
        self._config_source = Path(save_path)
        self._status_label.setText(f"Config saved: {Path(save_path).name}")
        logging.getLogger("emtrace").info("Config saved to %s", save_path)

    def _on_clear_preprocessing(self) -> None:
        """Delete preprocessing intermediate files after user confirmation."""
        if self._sample is None:
            return

        reads_dir = self._sample.reads_dir
        candidates = [
            reads_dir / "input.subsample.fq",
            reads_dir / "reads.porechop.fq.gz",
            reads_dir / "reads.chopped.fq.gz",
            reads_dir / "reads.trim.hflt.fq.gz",
        ]
        candidates += list(reads_dir.glob("hflt.*.reads.fq.gz")) if reads_dir.exists() else []
        to_delete = [p for p in candidates if p.exists()]

        if not to_delete:
            QMessageBox.information(
                self, "Nothing to clear",
                "No preprocessing intermediate files found.",
            )
            return

        file_list = "\n".join(f"  • {p.name}" for p in to_delete)
        reply = QMessageBox.question(
            self,
            "Clear preprocessing files?",
            f"The following files will be deleted:\n\n{file_list}\n\n"
            "Preprocessing steps will re-run from the original input reads.",
            QMessageBox.StandardButton.Ok | QMessageBox.StandardButton.Cancel,
        )
        if reply != QMessageBox.StandardButton.Ok:
            return

        removed, failed = 0, []
        for p in to_delete:
            try:
                p.unlink()
                removed += 1
            except OSError as exc:
                failed.append(f"{p.name}: {exc}")

        # Remove stale registry keys so current_reads falls back to input_reads
        for key in ("porechop_reads", "chopped_reads", "host_filtered_reads", "current_reads"):
            self._sample.files.pop(key, None)

        self._main_content.dashboard.refresh()

        if failed:
            QMessageBox.warning(
                self, "Partial clear",
                f"Deleted {removed} file(s). Could not delete:\n" + "\n".join(failed),
            )
        else:
            logging.getLogger("emtrace").info(
                "Cleared %d preprocessing file(s)", removed
            )

    # ------------------------------------------------------------------
    # Shared clear helper
    # ------------------------------------------------------------------

    def _clear_section_files(
        self,
        title: str,
        to_delete: list,
        post_message: str = "",
    ) -> None:
        """Confirm with the user, then delete *to_delete* and refresh.

        *to_delete* is a list of :class:`Path` objects that exist on disk.
        Registry keys whose paths are in *to_delete* are removed automatically.
        """
        if not to_delete:
            QMessageBox.information(self, "Nothing to clear",
                                    f"No {title.lower()} output files found.")
            return

        file_list = "\n".join(f"  • {p.name}" for p in to_delete)
        body = f"The following files will be deleted:\n\n{file_list}"
        if post_message:
            body += f"\n\n{post_message}"
        reply = QMessageBox.question(
            self, f"Clear {title} files?", body,
            QMessageBox.StandardButton.Ok | QMessageBox.StandardButton.Cancel,
        )
        if reply != QMessageBox.StandardButton.Ok:
            return

        delete_set = set(to_delete)
        removed, failed = 0, []
        for p in to_delete:
            try:
                p.unlink()
                removed += 1
            except OSError as exc:
                failed.append(f"{p.name}: {exc}")

        # Clear any registry keys pointing at deleted files.
        if self._sample is not None:
            for key in list(self._sample.files.keys()):
                if self._sample.files.get(key) in delete_set:
                    self._sample.files.pop(key, None)

        self._main_content.dashboard.refresh()

        if failed:
            QMessageBox.warning(
                self, "Partial clear",
                f"Deleted {removed} file(s). Could not delete:\n" + "\n".join(failed),
            )
        else:
            logging.getLogger("emtrace").info(
                "Cleared %d %s file(s)", removed, title.lower()
            )

    def _on_clear_annotation(self) -> None:
        """Delete all files in the annotation output directory."""
        if self._sample is None:
            return
        ann_dir = self._sample.annotation_dir
        to_delete = sorted(ann_dir.rglob("*") if ann_dir.exists() else [],
                           key=lambda p: p.name)
        to_delete = [p for p in to_delete if p.is_file()]
        self._clear_section_files(
            "Annotation", to_delete,
            "Annotation steps will need to re-run to regenerate these outputs.",
        )

    def _on_clear_em(self) -> None:
        """Delete EM abundance output files."""
        if self._sample is None:
            return
        rep_dir = self._sample.reports_dir
        candidates = [
            rep_dir / "EM.taxid_prob.tsv",
            rep_dir / "EM.taxid_prob.byiter.tsv",
        ]
        to_delete = [p for p in candidates if p.exists()]
        self._clear_section_files(
            "EM Abundance", to_delete,
            "The EM step will need to re-run to regenerate these outputs.",
        )

    def _on_clear_assembly(self) -> None:
        """Delete all files in the assembly output directory."""
        if self._sample is None:
            return
        asm_dir = self._sample.assembly_dir
        to_delete = sorted(asm_dir.rglob("*") if asm_dir.exists() else [],
                           key=lambda p: p.name)
        to_delete = [p for p in to_delete if p.is_file()]
        self._clear_section_files(
            "Assembly", to_delete,
            "The Assembly step will need to re-run to regenerate these outputs.",
        )

    def _on_clear_reporting(self) -> None:
        """Delete all reporting and plot output files."""
        if self._sample is None:
            return
        rep_dir = self._sample.reports_dir
        candidates: list = []
        for subdir in ("reports", "plots"):
            d = rep_dir / subdir
            if d.exists():
                candidates.extend(p for p in d.rglob("*") if p.is_file())
        stats_yaml = rep_dir / "stats.yaml"
        if stats_yaml.exists():
            candidates.append(stats_yaml)
        to_delete = sorted(set(candidates), key=lambda p: p.name)
        self._clear_section_files(
            "Reporting", to_delete,
            "The Reporting and Plots steps will need to re-run to regenerate these outputs.",
        )

    # ==================================================================
    # Pipeline run / stop
    # ==================================================================

    def _on_run(self) -> None:
        """Launch a :class:`PipelineWorker` for the checked steps."""
        if self._worker and self._worker.isRunning():
            return

        if self._sample is None:
            QMessageBox.warning(
                self, "No Job",
                "Create or open a job first (File \u25b6 New Job).",
            )
            return

        step_values = self._step_panel.checked_steps()
        if not step_values:
            QMessageBox.warning(
                self, "No Steps Selected",
                "Check at least one step on the left to run.",
            )
            return

        steps = [StepName(v) for v in step_values]

        # Flush any pending config-editor changes before starting the run so
        # the user never has to remember to click "Apply Changes" separately.
        self._main_content.config_page._on_apply()

        self._main_content.dashboard.clear_failures()
        self._step_panel.reset_status()

        self._worker = PipelineWorker(
            self._config, self._sample, steps, parent=self,
        )
        self._worker.signals.step_started.connect(self._on_step_started)
        self._worker.signals.step_finished.connect(self._on_step_finished)
        self._worker.signals.pipeline_finished.connect(self._on_pipeline_done)
        self._worker.signals.error.connect(self._on_pipeline_error)
        self._worker.signals.job_updated.connect(self._on_job_updated)
        self._worker.signals.progress.connect(self._on_progress)

        self._toolbar.set_running(True)
        self._status_label.setText("Running\u2026")
        self._worker.start()

    def _on_stop(self) -> None:
        """Cooperative cancel — finishes the current step, then stops."""
        if self._worker and self._worker.isRunning():
            self._worker.request_cancel()
            self._status_label.setText("Cancelling\u2026")

    def _on_force_quit(self) -> None:
        """Immediately kill the running step (with confirmation)."""
        if not (self._worker and self._worker.isRunning()):
            return
        reply = QMessageBox.question(
            self,
            "Force Stop Pipeline?",
            "This will immediately kill the running external process.\n"
            "Partial output files may be left on disk.\n\nProceed?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.Yes,
        )
        if reply != QMessageBox.StandardButton.Yes:
            return
        signalled = self._worker.force_quit()
        self._status_label.setText(
            "Force-killing\u2026" if signalled
            else "No active subprocess; cancelling\u2026"
        )

    # ------------------------------------------------------------------
    # Worker signal handlers
    # ------------------------------------------------------------------

    def _on_step_started(self, step_value: str) -> None:
        self._step_panel.mark_running(step_value)
        self._status_label.setText(f"Running: {step_value}")

    def _on_step_finished(self, step_value: str, success: bool) -> None:
        self._step_panel.mark_done(step_value, success)
        if not success:
            self._main_content.dashboard.mark_step_failed(step_value)

    def _on_progress(self, current: int, total: int) -> None:
        self._status_label.setText(f"Step {current}/{total}")

    def _on_job_updated(self, sample: PipelineJob) -> None:
        self._sample = sample
        self._main_content.set_job(sample)

    def _on_pipeline_done(self, success: bool) -> None:
        self._toolbar.set_running(False)
        self._status_label.setText(
            "Pipeline completed" if success else "Pipeline stopped",
        )
        # Refresh result viewers to pick up newly generated outputs
        if self._sample is not None:
            self._main_content.refresh_results(self._sample)

    def _on_pipeline_error(self, message: str) -> None:
        self._toolbar.set_running(False)
        if not self._closing:
            QMessageBox.critical(self, "Pipeline Error", message)

    # ==================================================================
    # Misc slots
    # ==================================================================

    def _on_config_changed(self, config: EMtraceConfig) -> None:
        """Config editor Apply → update active config in memory."""
        self._config = config
        self._status_label.setText("Config updated (not saved to disk).")
        logging.getLogger("emtrace").info("Config updated via editor.")

    def _on_view_step(self, step_value: str) -> None:
        """Step-panel double-click on the status column → jump to Files tab."""
        logging.getLogger("emtrace").debug(
            "View requested for step: %s", step_value,
        )
        self._main_content.switch_to(MainContentStack.TAB_FILES)

    def _on_file_activated(self, path: str, kind: str) -> None:
        """Route a double-clicked file from the file browser to the correct viewer."""
        p = Path(path)
        if kind == "table":
            self._main_content.table_viewer.load_file(p)
            self._main_content.switch_to(MainContentStack.TAB_TABLES)
        elif kind == "figure":
            self._main_content.figure_viewer.load_figure(p)
            self._main_content.switch_to(MainContentStack.TAB_FIGURES)
        elif kind == "fastq":
            self._main_content.qc_plot.load_fastq(p)
            self._main_content.switch_to(MainContentStack.TAB_QC)
        elif kind in ("log", "text"):
            dlg = _LogViewerDialog(p, parent=self)
            if kind == "text":
                dlg._text.moveCursor(dlg._text.textCursor().MoveOperation.Start)
            dlg.show()
        else:
            # For unknown types, try opening as a table (TSV/CSV are common)
            self._main_content.table_viewer.load_file(p)
            self._main_content.switch_to(MainContentStack.TAB_TABLES)

    def _on_database_manager(self) -> None:
        """Tools ▶ Database Manager…"""
        from emtrace.gui.widgets.db_manager import DatabaseManagerDialog

        dlg = DatabaseManagerDialog(
            self._config,
            config_source=self._config_source,
            parent=self,
        )
        dlg.config_updated.connect(lambda cfg: setattr(self, "_config", cfg))
        dlg.exec()

    def _on_open_manual(self) -> None:
        """Help > User Guide — open docs/UserGuide.v1.0.md with the system default app."""
        import os
        import subprocess
        import sys
        from emtrace.gui.platform import bundle_root
        manual = bundle_root() / "docs" / "UserGuide.v1.0.md"
        if not manual.exists():
            QMessageBox.warning(
                self, "User guide not found",
                f"Could not find the user guide at:\n{manual}",
            )
            return
        try:
            if sys.platform == "darwin":
                subprocess.Popen(["open", str(manual)])
            elif sys.platform.startswith("linux"):
                subprocess.Popen(["xdg-open", str(manual)])
            else:
                os.startfile(str(manual))  # type: ignore[attr-defined]
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(
                self, "Cannot open user guide",
                f"Failed to open the user guide:\n{exc}",
            )

    def _on_about(self) -> None:
        QMessageBox.about(
            self, "About EMtrace",
            f"<b>EMtrace</b> v{__version__}<br><br>"
            "Nanopore viral metagenomics pipeline.<br>"
            "Desktop GUI powered by PyQt6.",
        )

    # ==================================================================
    # Lifecycle
    # ==================================================================

    # ==================================================================
    # Realtime watch handlers
    # ==================================================================

    def _on_watch_start(self, watch_dir: Path, output_dir: Path,
                        batch_size: int, debounce: float) -> None:
        if self._realtime_worker and self._realtime_worker.isRunning():
            return
        cfg = self._config
        ann1 = cfg.annotation.ann1
        ann2 = cfg.annotation.ann2
        anns = cfg.annotation.anns
        hostgen = cfg.preprocess.bgfilter.bgfilter

        self._realtime_worker = RealtimeWorker(
            config=cfg,
            watch_dir=watch_dir,
            output_dir=output_dir,
            ann1=ann1,
            ann2=ann2,
            anns=anns,
            hostgen=hostgen,
            batch_size=batch_size,
            debounce=debounce,
        )
        rp = self._main_content.realtime_panel
        w = self._realtime_worker
        w.signals.batch_started.connect(rp.on_batch_started)
        w.signals.batch_finished.connect(rp.on_batch_finished)
        w.signals.abundance_updated.connect(rp.on_abundance_updated)
        w.signals.watch_status_changed.connect(rp.on_watch_status_changed)
        w.signals.error.connect(rp.on_error)
        w.signals.worker_finished.connect(self._on_watch_done)
        w.start()
        self._main_content.switch_to(MainContentStack.TAB_REALTIME)

    def _on_watch_stop(self) -> None:
        if self._realtime_worker:
            self._realtime_worker.request_stop()

    def _on_watch_done(self) -> None:
        self._realtime_worker = None

    def closeEvent(self, event) -> None:  # noqa: N802 (Qt naming)
        """Detach the log handler before exit so the root logger is clean."""
        try:
            logging.getLogger("emtrace").removeHandler(self._log_handler)
        except Exception:  # noqa: BLE001
            pass

        if self._worker and self._worker.isRunning():
            self._worker.request_cancel()
            if not self._worker.wait(2000):
                # Cooperative cancel didn't release within 2 s — there's
                # an active subprocess. Force-kill it so the GUI can exit
                # cleanly instead of orphaning the child tree.
                self._worker.force_quit()
                self._worker.wait(3000)

        super().closeEvent(event)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def launch_gui(config_file: Optional[str] = None) -> int:
    """Construct the QApplication and run the EMtrace GUI event loop.

    Returns the application's exit code, mirroring ``QApplication.exec()``.
    """
    config_source: Optional[Path] = None
    if config_file:
        config_source = Path(config_file)
        config = EMtraceConfig.from_yaml(config_source)
    else:
        config = EMtraceConfig.create_default()

    app = QApplication.instance() or QApplication(sys.argv)
    app.setApplicationName("EMtrace")
    app.setOrganizationName("NGS")
    configure_platform(app)

    window = MainWindow(config, config_source=config_source)
    window.show()

    # Help users find the bundled docs directory in the log if anything
    # below the GUI tries to load resources from it.
    logging.getLogger("emtrace").debug(
        "Bundle root: %s", bundle_root(),
    )

    return app.exec()
