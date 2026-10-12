"""Dashboard widget — sample state at a glance.

Two states:

  - **Empty**: shown before any sample is loaded; offers a "New Sample"
    call-to-action that re-emits the same signal as the toolbar's New
    button.
  - **Loaded**: shows sample identity + per-step status cards grouped
    by phase.

Status per step:
  - ``\u2713``  green   — output file present in ``job.files`` and on disk
  - ``\u2717``  red     — step ran but failed in this session
  - ``\u25cb``  grey    — step has not run yet (or output not produced)
"""

from __future__ import annotations

from typing import Dict, List, Optional, Set, Tuple

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (
    QFrame,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from emtrace.core.pipelinejob import PipelineJob
from emtrace.steps import StepName

# ---------------------------------------------------------------------------
# Glyphs
# ---------------------------------------------------------------------------
_OK = "\u2713"        # ✓
_FAIL = "\u2717"      # ✗
_PENDING = "\u25cb"   # ○

_COLOR_OK = "#16a34a"
_COLOR_FAIL = "#dc2626"
_COLOR_PENDING = "#9ca3af"


# Map StepName → (display label, list of file keys whose presence proves
# the step ran).  A step shows ✓ when *any* of the listed keys is on disk.
_STEP_OUTPUTS: Dict[StepName, Tuple[str, List[str]]] = {
    StepName.PORE:  ("Adapter removal (Porechop/Porechop_ABI)", ["porechop_reads"]),
    StepName.CHOP:  ("Quality filtering (Chopper)", ["chopped_reads"]),
    StepName.FLT:   ("Host filtering", ["host_filtered_reads"]),
    StepName.ANN:   ("Annotation (rounds 1 & 2)", ["annot1", "annot2", "annotation_tsv"]),
    StepName.EM:    ("EM abundance", ["em_taxid_prob"]),
    StepName.ASS:   ("Reference assembly", ["assembly_consensus", "assembly_bam"]),
    StepName.REP:   ("Reports + Statistics", ["abund_table", "annot_table", "stats_html"]),
    StepName.PLOT:  ("Plots", ["plots_sunburst"]),
    StepName.PACK:  ("Pack output", ["pack_archive"]),
    StepName.CLEAN: ("Clean intermediates", []),
}

# Phase grouping mirrors the StepPanel groups for visual consistency.
_PHASE_GROUPS: List[Tuple[str, List[StepName]]] = [
    ("Preprocessing",   [StepName.PORE, StepName.CHOP, StepName.FLT]),
    ("Annotation",      [StepName.ANN]),
    ("EM Abundance",    [StepName.EM]),
    ("Assembly",        [StepName.ASS]),
    ("Reporting",       [StepName.REP, StepName.PLOT]),
    ("Pack & Clean",    [StepName.PACK, StepName.CLEAN]),
]


class _StatusCard(QWidget):
    """Single icon + label status row inside a phase group."""

    def __init__(self, label: str, state: str, parent=None) -> None:
        super().__init__(parent)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(4, 2, 4, 2)
        layout.setSpacing(6)

        symbol, color = {
            "ok":      (_OK, _COLOR_OK),
            "fail":    (_FAIL, _COLOR_FAIL),
            "pending": (_PENDING, _COLOR_PENDING),
        }[state]

        icon = QLabel(symbol)
        icon.setFixedWidth(16)
        icon.setAlignment(Qt.AlignmentFlag.AlignCenter)
        icon.setStyleSheet(f"color: {color}; font-weight: bold;")

        text = QLabel(label)
        text.setStyleSheet("color: #333; font-size: 12px;")

        layout.addWidget(icon)
        layout.addWidget(text, stretch=1)


class DashboardWidget(QWidget):
    """Empty/loaded dashboard for the active :class:`PipelineJob`.

    Signals
    -------
    new_sample_requested
        Emitted when the user clicks the "New Sample" button in the
        empty-state placeholder.
    clear_preprocessing_requested
        Emitted when the user clicks "Clear" in the Preprocessing group.
        MainWindow handles the confirmation dialog and file deletion.
    """

    new_sample_requested = pyqtSignal()
    clear_preprocessing_requested = pyqtSignal()
    clear_annotation_requested = pyqtSignal()
    clear_em_requested = pyqtSignal()
    clear_assembly_requested = pyqtSignal()
    clear_reporting_requested = pyqtSignal()

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._job: Optional[PipelineJob] = None
        self._failed_steps: Set[str] = set()

        self._root_layout = QVBoxLayout(self)
        self._root_layout.setContentsMargins(12, 12, 12, 12)
        self._root_layout.setSpacing(8)

        # --- Empty state ---
        self._empty_widget = QWidget()
        empty_layout = QVBoxLayout(self._empty_widget)
        empty_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        empty_layout.setSpacing(12)
        prompt = QLabel("No pipeline job loaded")
        prompt.setAlignment(Qt.AlignmentFlag.AlignCenter)
        prompt.setStyleSheet("color: #999; font-size: 16px;")
        hint = QLabel("Use File \u203a New Job to start.")
        hint.setAlignment(Qt.AlignmentFlag.AlignCenter)
        hint.setStyleSheet("color: #aaa; font-size: 12px;")
        new_btn = QPushButton("New Job\u2026")
        new_btn.setFixedWidth(160)
        new_btn.clicked.connect(self.new_sample_requested.emit)
        empty_layout.addStretch(1)
        empty_layout.addWidget(prompt)
        empty_layout.addWidget(hint)
        empty_layout.addWidget(new_btn, alignment=Qt.AlignmentFlag.AlignCenter)
        empty_layout.addStretch(2)
        self._root_layout.addWidget(self._empty_widget, stretch=1)

        # --- Loaded-state container (hidden initially) ---
        self._content = QWidget()
        self._content_layout = QVBoxLayout(self._content)
        self._content_layout.setContentsMargins(0, 0, 0, 0)
        self._content_layout.setSpacing(8)
        self._content.setVisible(False)
        self._root_layout.addWidget(self._content, stretch=1)

    # ==================================================================
    # Public slots
    # ==================================================================

    def set_job(self, job: Optional[PipelineJob]) -> None:
        """Switch to *job* (or back to the empty state when ``None``).

        Failure marks accumulated during the current run are preserved
        when *job* is the same instance as the current one — the worker
        re-emits the (mutated) sample on every ``job_updated`` and we
        don't want that to wipe the ✗ icons from steps that just failed.
        """
        if self._job is not job:
            self._failed_steps.clear()
        self._job = job
        self.refresh()

    def mark_step_failed(self, step_value: str) -> None:
        """Record a failure so the matching card switches to ✗ on next refresh."""
        self._failed_steps.add(step_value)
        self.refresh()

    def clear_failures(self) -> None:
        """Reset session-level failure flags (e.g. before a fresh run)."""
        self._failed_steps.clear()
        self.refresh()

    def refresh(self) -> None:
        """Rebuild the dashboard from the current job state."""
        if self._job is None:
            self._empty_widget.setVisible(True)
            self._content.setVisible(False)
            return

        self._clear_layout(self._content_layout)

        # --- Header: name + paths ---
        title = QLabel(f"<b>{self._job.name or 'Sample'}</b>")
        title.setStyleSheet("font-size: 16px; color: #222;")
        self._content_layout.addWidget(title)

        meta_lines = []
        if self._job.input_reads:
            r1_label = "R1:" if self._job.is_paired else "Input:"
            meta_lines.append(f"{r1_label} {self._job.input_reads}")
        if self._job.is_paired and self._job.input_reads_r2:
            meta_lines.append(f"R2: {self._job.input_reads_r2}")
        platform_badge = (
            "Platform: Illumina paired-end"
            if self._job.platform == "illumina_pe"
            else "Platform: Nanopore"
        )
        meta_lines.append(platform_badge)
        meta_lines.append(f"Output: {self._job.output_dir}")
        meta_label = QLabel("<br>".join(meta_lines))
        meta_label.setStyleSheet("font-size: 11px; color: #888;")
        meta_label.setTextInteractionFlags(
            Qt.TextInteractionFlag.TextSelectableByMouse,
        )
        self._content_layout.addWidget(meta_label)

        # --- Separator ---
        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        sep.setStyleSheet("color: #ddd;")
        self._content_layout.addWidget(sep)

        # --- Status grid ---
        self._content_layout.addWidget(self._build_status_grid())
        self._content_layout.addStretch(1)

        self._empty_widget.setVisible(False)
        self._content.setVisible(True)

    # ==================================================================
    # Builders
    # ==================================================================

    def _build_status_grid(self) -> QWidget:
        container = QWidget()
        grid = QGridLayout(container)
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(16)
        grid.setVerticalSpacing(8)

        for index, (group_label, steps) in enumerate(_PHASE_GROUPS):
            box = QGroupBox(group_label)
            box_layout = QVBoxLayout(box)
            box_layout.setContentsMargins(6, 6, 6, 6)
            box_layout.setSpacing(2)

            for step in steps:
                label, _ = _STEP_OUTPUTS.get(step, (step.value, []))
                state = self._step_state(step)
                box_layout.addWidget(_StatusCard(label, state))

            signal_map = {
                "Preprocessing": (self.clear_preprocessing_requested,
                                  [StepName.PORE, StepName.CHOP, StepName.FLT]),
                "Annotation":    (self.clear_annotation_requested, [StepName.ANN]),
                "EM Abundance":  (self.clear_em_requested,         [StepName.EM]),
                "Assembly":      (self.clear_assembly_requested,   [StepName.ASS]),
                "Reporting":     (self.clear_reporting_requested,  [StepName.REP, StepName.PLOT]),
            }
            if group_label in signal_map:
                signal, check_steps = signal_map[group_label]
                clear_btn = QPushButton("Clear files…")
                clear_btn.setEnabled(self._any_steps_done(check_steps))
                clear_btn.clicked.connect(signal.emit)
                box_layout.addWidget(clear_btn)

            box_layout.addStretch(1)
            row, col = divmod(index, 3)
            grid.addWidget(box, row, col)

        return container

    def _any_preprocessing_done(self) -> bool:
        return self._any_steps_done([StepName.PORE, StepName.CHOP, StepName.FLT])

    def _any_steps_done(self, steps: List[StepName]) -> bool:
        """Return True when at least one output file for *steps* exists on disk."""
        if self._job is None:
            return False
        for step in steps:
            _, keys = _STEP_OUTPUTS.get(step, ("", []))
            for key in keys:
                if self._job.has_file(key):
                    return True
        return False

    def _step_state(self, step: StepName) -> str:
        """Return ``"ok" | "fail" | "pending"`` for a step."""
        if step.value in self._failed_steps:
            return "fail"
        if self._job is None:
            return "pending"
        _, file_keys = _STEP_OUTPUTS.get(step, ("", []))
        for key in file_keys:
            if self._job.has_file(key):
                return "ok"
        return "pending"

    @staticmethod
    def _clear_layout(layout) -> None:
        """Remove and delete every child widget from *layout*."""
        while layout.count():
            item = layout.takeAt(0)
            w = item.widget()
            if w is not None:
                w.deleteLater()
