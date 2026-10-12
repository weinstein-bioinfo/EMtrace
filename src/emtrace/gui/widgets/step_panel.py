"""Left panel — EMtrace step list with checkboxes grouped by category.

Status indicators (right-aligned column 1):

    idle     — (empty)
    running  — ●  blue
    success  — ✓  green
    error    — ✗  red

Adapted from the ViRefSet StepPanel pattern; group structure mirrors
the actual ``StepName`` enum in :mod:`emtrace.steps`.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Set

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QBrush, QColor
from PyQt6.QtWidgets import (
    QHeaderView,
    QLabel,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from emtrace.steps import STEP_REGISTRY, StepName

# ---------------------------------------------------------------------------
# Glyphs
# ---------------------------------------------------------------------------
_COLLAPSED = "\u25b6"        # ▶
_EXPANDED = "\u25bc"         # ▼

_STATUS_IDLE = ""
_STATUS_RUNNING = "\u25cf"   # ●
_STATUS_SUCCESS = "\u2713"   # ✓
_STATUS_ERROR = "\u2717"     # ✗

_COLOR_RUNNING = QColor("#2563eb")   # blue
_COLOR_SUCCESS = QColor("#16a34a")   # green
_COLOR_ERROR = QColor("#dc2626")     # red
_COLOR_DEFAULT = QColor("#222222")
_COLOR_DISABLED = QColor("#aaaaaa")  # platform-disabled steps


class StepPanel(QWidget):
    """Step tree with checkboxes, status icons, and a quick-info summary.

    Signals
    -------
    view_requested(str)
        Emitted when the user clicks the status column for a step
        (intended to switch the main view to that step's results).
    """

    view_requested = pyqtSignal(str)

    # Group → ordered list of StepName values.  Keep ``CHI`` collapsed
    # under preprocessing and ``ASS`` separate; ``PLOT`` lives with
    # reporting since both consume the abundance table.
    _STEP_GROUPS: Dict[str, List[StepName]] = {
        "Preprocessing": [StepName.PORE, StepName.CHOP, StepName.FLT],
        "Annotation": [StepName.ANN],
        "EM Abundance": [StepName.EM],
        "Assembly": [StepName.ASS],
        "Reporting": [StepName.REP, StepName.PLOT],
        "Pack & Clean": [StepName.PACK, StepName.CLEAN],
    }

    _STEP_DATA_ROLE = Qt.ItemDataRole.UserRole
    _GROUP_NAME_ROLE = Qt.ItemDataRole.UserRole + 1

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setMinimumWidth(180)
        self.setMaximumWidth(280)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        # Two-column tree: [name] [status]
        self._tree = QTreeWidget()
        self._tree.setColumnCount(2)
        self._tree.setHeaderHidden(True)
        self._tree.setRootIsDecorated(False)
        self._tree.setIndentation(16)

        header = self._tree.header()
        header.setStretchLastSection(False)
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.Fixed)
        header.resizeSection(1, 28)

        self._tree.itemClicked.connect(self._on_item_clicked)
        self._tree.itemChanged.connect(self._on_item_changed)
        self._tree.itemExpanded.connect(self._on_expanded)
        self._tree.itemCollapsed.connect(self._on_collapsed)

        # itemChanged fires before itemClicked for direct checkbox clicks;
        # this flag prevents double-toggling when the user clicks the box.
        self._native_check_toggle = False

        # Steps disabled because they are incompatible with the current job's
        # platform (e.g. Nanopore-only steps for illumina_pe jobs). Tracked
        # separately from the Qt item flags so clicks can be reliably blocked
        # regardless of how the platform renders disabled tree items.
        self._platform_disabled: Set[str] = set()

        self._step_items: Dict[str, QTreeWidgetItem] = {}
        self._build_tree()

        layout.addWidget(self._tree, stretch=1)

        # Bottom quick-info label
        self._info_label = QLabel("No sample loaded")
        self._info_label.setWordWrap(True)
        self._info_label.setStyleSheet("padding: 6px; color: #666;")
        layout.addWidget(self._info_label)

    # ==================================================================
    # Tree construction
    # ==================================================================

    def _build_tree(self) -> None:
        for group_name, steps in self._STEP_GROUPS.items():
            label = f"{_EXPANDED}  {group_name}"
            group_item = QTreeWidgetItem(self._tree, [label, ""])
            group_item.setData(0, self._GROUP_NAME_ROLE, group_name)
            group_item.setFlags(
                group_item.flags() & ~Qt.ItemFlag.ItemIsSelectable,
            )
            font = group_item.font(0)
            font.setBold(True)
            group_item.setFont(0, font)

            for step in steps:
                cls = STEP_REGISTRY.get(step)
                display = cls.step_name if cls else step.value
                tooltip = cls.step_description if cls else ""
                item = QTreeWidgetItem(group_item, [display, ""])
                item.setFlags(
                    item.flags() | Qt.ItemFlag.ItemIsUserCheckable,
                )
                item.setCheckState(0, Qt.CheckState.Unchecked)
                item.setData(0, self._STEP_DATA_ROLE, step.value)
                item.setTextAlignment(1, Qt.AlignmentFlag.AlignCenter)
                if tooltip:
                    item.setToolTip(0, tooltip)
                self._step_items[step.value] = item

            group_item.setExpanded(True)

    def _on_expanded(self, item: QTreeWidgetItem) -> None:
        name = item.data(0, self._GROUP_NAME_ROLE)
        if name:
            item.setText(0, f"{_EXPANDED}  {name}")

    def _on_collapsed(self, item: QTreeWidgetItem) -> None:
        name = item.data(0, self._GROUP_NAME_ROLE)
        if name:
            item.setText(0, f"{_COLLAPSED}  {name}")

    # ==================================================================
    # Public API
    # ==================================================================

    def checked_steps(self) -> List[str]:
        """Return checked step names in canonical ``StepName`` order."""
        checked = {
            step_val
            for step_val, item in self._step_items.items()
            if item.checkState(0) == Qt.CheckState.Checked
        }
        return [s.value for s in StepName if s.value in checked]

    def set_checked(self, step_values: List[str], checked: bool = True) -> None:
        """Set the check state for a list of step values."""
        state = Qt.CheckState.Checked if checked else Qt.CheckState.Unchecked
        for step_val in step_values:
            item = self._step_items.get(step_val)
            if item is not None:
                item.setCheckState(0, state)

    def set_step_enabled(self, step_value: str, enabled: bool) -> None:
        """Enable or disable a single step (e.g. when its config is missing).

        Disabled steps are unchecked and greyed out.
        """
        item = self._step_items.get(step_value)
        if item is None:
            return
        if enabled:
            item.setFlags(
                item.flags()
                | Qt.ItemFlag.ItemIsEnabled
                | Qt.ItemFlag.ItemIsUserCheckable
            )
            item.setForeground(0, QBrush(_COLOR_DEFAULT))
        else:
            item.setCheckState(0, Qt.CheckState.Unchecked)
            item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEnabled)
            item.setForeground(0, QBrush(_COLOR_DISABLED))

    def mark_running(self, step_value: str) -> None:
        """Show the running indicator (blue ●) and reveal the parent group."""
        self._set_status(step_value, _STATUS_RUNNING, _COLOR_RUNNING)
        item = self._step_items.get(step_value)
        if item and item.parent():
            item.parent().setExpanded(True)

    def mark_done(self, step_value: str, success: bool) -> None:
        """Show success (✓ green) or error (✗ red).  Uncheck on success."""
        if success:
            self._set_status(step_value, _STATUS_SUCCESS, _COLOR_SUCCESS)
            item = self._step_items.get(step_value)
            if item:
                item.setCheckState(0, Qt.CheckState.Unchecked)
        else:
            self._set_status(step_value, _STATUS_ERROR, _COLOR_ERROR)
            item = self._step_items.get(step_value)
            if item and item.parent():
                item.parent().setExpanded(True)

    def reset_all(self) -> None:
        """Clear status indicators and uncheck every step."""
        for item in self._step_items.values():
            item.setText(1, _STATUS_IDLE)
            item.setForeground(0, QBrush(_COLOR_DEFAULT))
            item.setCheckState(0, Qt.CheckState.Unchecked)

    def reset_status(self) -> None:
        """Clear status indicators only (preserve check states)."""
        for step_val, item in self._step_items.items():
            item.setText(1, _STATUS_IDLE)
            color = _COLOR_DISABLED if step_val in self._platform_disabled else _COLOR_DEFAULT
            item.setForeground(0, QBrush(color))

    def update_quick_info(
        self,
        name: str,
        base_dir: str,
        reads: Optional[int] = None,
    ) -> None:
        """Update the bottom info label with sample-level summary."""
        parts = [f"<b>{name}</b>", f"Dir: {base_dir}"]
        if reads is not None:
            parts.append(f"Reads: {reads:,}")
        self._info_label.setText("<br>".join(parts))

    # Steps that currently have no paired-end / Illumina implementation —
    # disabled when the active job is `platform == "illumina_pe"` so the
    # user cannot tick a step that would crash on PE input.
    _NANOPORE_ONLY_STEPS = (
        StepName.PORE,
        StepName.CHOP,
        StepName.FLT,
        StepName.ANN,
        StepName.EM,
    )

    def apply_platform(self, platform: str) -> None:
        """Enable/disable steps based on the job's sequencing platform.

        For ``illumina_pe`` the Nanopore-only preprocessing/annotation
        steps are disabled (and unchecked); for ``nanopore`` they are
        re-enabled provided the step has a registered implementation.
        """
        nanopore_only = {s.value for s in self._NANOPORE_ONLY_STEPS}
        is_pe = platform == "illumina_pe"

        # Block signals so programmatic setCheckState calls don't trip
        # _native_check_toggle and swallow the next user click.
        self._tree.blockSignals(True)
        try:
            for step_val in nanopore_only:
                if is_pe:
                    self._platform_disabled.add(step_val)
                    self.set_step_enabled(step_val, False)
                else:
                    self._platform_disabled.discard(step_val)
                    # Re-enable only when the step is actually registered.
                    step = next((s for s in StepName if s.value == step_val), None)
                    if step is not None and step in STEP_REGISTRY:
                        self.set_step_enabled(step_val, True)
        finally:
            self._tree.blockSignals(False)

    # ==================================================================
    # Event handlers
    # ==================================================================

    def _set_status(self, step_value: str, symbol: str, color: QColor) -> None:
        item = self._step_items.get(step_value)
        if item is None:
            return
        item.setText(1, symbol)
        item.setForeground(1, QBrush(color))
        font = item.font(1)
        font.setBold(True)
        item.setFont(1, font)

    def _on_item_changed(self, item: QTreeWidgetItem, column: int) -> None:
        if column == 0 and item.data(0, self._STEP_DATA_ROLE):
            self._native_check_toggle = True

    def _on_item_clicked(self, item: QTreeWidgetItem, column: int) -> None:
        # Group header → toggle expand/collapse
        if item.data(0, self._GROUP_NAME_ROLE):
            item.setExpanded(not item.isExpanded())
            return

        step_val = item.data(0, self._STEP_DATA_ROLE)
        if not step_val:
            return

        # Platform-disabled steps are non-interactive — undo any Qt-side
        # checkbox change that may have slipped through and clear the toggle flag.
        if step_val in self._platform_disabled:
            self._tree.blockSignals(True)
            item.setCheckState(0, Qt.CheckState.Unchecked)
            self._tree.blockSignals(False)
            self._native_check_toggle = False
            return

        # Status column → "view this step's output"
        if column == 1:
            self.view_requested.emit(step_val)
            return

        # Column 0: native checkbox click already toggled state via Qt.
        # If the click landed on the row text, toggle it ourselves.
        if self._native_check_toggle:
            self._native_check_toggle = False
        else:
            current = item.checkState(0)
            item.setCheckState(
                0,
                Qt.CheckState.Unchecked
                if current == Qt.CheckState.Checked
                else Qt.CheckState.Checked,
            )
            self._native_check_toggle = False
