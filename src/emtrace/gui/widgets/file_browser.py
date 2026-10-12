"""File Browser widget — QTreeView over the sample output directory.

Adapted from ViRefSet ``file_browser.py``.  Adds a virtual *Kind*
column to the standard :class:`QFileSystemModel` columns so the user
can see at a glance whether a file is a Table, Figure, FASTQ, etc.

Double-clicking a file emits :pyqtSignal:`file_activated` with the
absolute path and a kind hint (``'table'``, ``'figure'``, ``'fastq'``,
or ``'other'``).  The main window uses the hint to jump to the
appropriate viewer tab.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from PyQt6.QtCore import QModelIndex, QSortFilterProxyModel, QUrl, Qt, pyqtSignal
from PyQt6.QtGui import QDesktopServices, QFileSystemModel
from PyQt6.QtWidgets import (
    QLabel,
    QTreeView,
    QVBoxLayout,
    QWidget,
)


# ── Kind classification ──────────────────────────────────────────────

_KIND_MAP: dict[str, str] = {
    # Tables
    ".csv": "Table",
    ".tsv": "Table",
    ".txt": "Table",
    ".xlsx": "Table",
    # Figures
    ".png": "Figure",
    ".svg": "Figure",
    ".pdf": "Figure",
    # Sequences
    ".fasta": "FASTA",
    ".fa": "FASTA",
    ".fna": "FASTA",
    ".fq": "FASTQ",
    ".fastq": "FASTQ",
    # Alignment / BAM
    ".bam": "BAM",
    ".bai": "Index",
    ".sam": "SAM",
    ".paf": "PAF",
    ".vcf": "VCF",
    # Config / structured data
    ".yaml": "Config",
    ".yml": "Config",
    ".json": "Config",
    # Interactive figures
    ".html": "Figure",
    ".htm": "Figure",
    # Logs
    ".log": "Log",
    # Compressed
    ".gz": "Archive",
    ".tar": "Archive",
    ".bed": "BED",
}


def _kind_for_path(p: Path) -> str:
    """Classify a file by its suffix (handles double-suffixes like .fq.gz)."""
    name = p.name.lower()
    # Two-part extensions
    if name.endswith(".fq.gz") or name.endswith(".fastq.gz"):
        return "FASTQ"
    if name.endswith(".fa.gz") or name.endswith(".fasta.gz") or name.endswith(".fna.gz"):
        return "FASTA"
    if name.endswith(".vcf.gz"):
        return "VCF"
    return _KIND_MAP.get(p.suffix.lower(), "Other")


# ── QFileSystemModel subclass with virtual Kind column ────────────────


class _KindFileSystemModel(QFileSystemModel):
    """QFileSystemModel extended with a virtual Kind column at index 4."""

    _KIND_COL = 4  # after Name(0), Size(1), Type(2), Date Modified(3)

    def columnCount(self, parent: QModelIndex = QModelIndex()) -> int:  # noqa: N802
        return self._KIND_COL + 1

    def headerData(self, section: int, orientation, role=Qt.ItemDataRole.DisplayRole):  # noqa: N802
        if (
            orientation == Qt.Orientation.Horizontal
            and role == Qt.ItemDataRole.DisplayRole
            and section == self._KIND_COL
        ):
            return "Kind"
        return super().headerData(section, orientation, role)

    def data(self, index: QModelIndex, role=Qt.ItemDataRole.DisplayRole):
        if not index.isValid():
            return None
        if index.column() == self._KIND_COL:
            if role in (Qt.ItemDataRole.DisplayRole, Qt.ItemDataRole.EditRole):
                fp = self.filePath(index.siblingAtColumn(0))
                if fp:
                    p = Path(fp)
                    return "" if p.is_dir() else _kind_for_path(p)
            return None
        return super().data(index, role)

    def flags(self, index: QModelIndex) -> Qt.ItemFlag:
        if index.column() == self._KIND_COL:
            return Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable
        return super().flags(index)


# ── Sort proxy: correct ordering for the Kind column ─────────────────


class _KindSortProxy(QSortFilterProxyModel):
    """Sorts the Kind column lexicographically; delegates others to source."""

    def lessThan(self, left: QModelIndex, right: QModelIndex) -> bool:  # noqa: N802
        if left.column() == _KindFileSystemModel._KIND_COL:
            l_val = self.sourceModel().data(left) or ""
            r_val = self.sourceModel().data(right) or ""
            return l_val < r_val
        return super().lessThan(left, right)


# ── File browser widget ─────────────────────────────────────────────


class FileBrowserWidget(QWidget):
    """QTreeView showing the sample output directory structure.

    Double-clicking a file emits ``file_activated`` with the absolute
    path and a kind hint string (``'table'``, ``'figure'``, ``'fastq'``,
    or ``'other'``).
    """

    # (absolute_path, kind_hint)
    file_activated = pyqtSignal(str, str)

    _TABLE_SUFFIXES = {".csv", ".tsv", ".txt"}
    _EXTERNAL_SUFFIXES = {".xlsx", ".fa", ".fasta", ".fna", ".bam", ".sam"}
    _FIGURE_SUFFIXES = {".png", ".svg", ".pdf", ".html", ".htm"}
    _FASTQ_SUFFIXES = {".fq", ".fastq"}
    _LOG_SUFFIXES = {".log"}
    _TEXT_SUFFIXES = {".yaml", ".yml", ".json"}

    def __init__(self, parent=None):
        super().__init__(parent)

        self._model = _KindFileSystemModel()
        self._model.setRootPath("")
        self._model.setReadOnly(True)

        self._proxy = _KindSortProxy(self)
        self._proxy.setSourceModel(self._model)

        self._tree = QTreeView()
        self._tree.setModel(self._proxy)
        self._tree.setColumnWidth(0, 280)
        self._tree.setColumnWidth(self._model._KIND_COL, 80)
        self._tree.hideColumn(2)  # hide generic 'Type' column
        self._tree.setAnimated(True)
        self._tree.setSortingEnabled(True)
        self._tree.sortByColumn(0, Qt.SortOrder.AscendingOrder)

        # Empty-state label
        self._empty_label = QLabel(
            "No sample loaded.\nOpen or create a sample to browse output files."
        )
        self._empty_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._empty_label.setStyleSheet("color: #999; font-size: 14px;")

        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addWidget(self._tree, stretch=1)
        layout.addWidget(self._empty_label, stretch=1)

        self._tree.setVisible(False)
        self._empty_label.setVisible(True)

        self._tree.doubleClicked.connect(self._on_double_click)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def set_root(self, path: Optional[Path | str]) -> None:
        """Set the root directory shown in the tree.

        Pass ``None`` to revert to the empty state.
        """
        if path is None:
            self._tree.setVisible(False)
            self._empty_label.setVisible(True)
            return

        src_idx = self._model.setRootPath(str(path))
        self._tree.setRootIndex(self._proxy.mapFromSource(src_idx))
        self._tree.setVisible(True)
        self._empty_label.setVisible(False)

    # ------------------------------------------------------------------
    # Private
    # ------------------------------------------------------------------

    def _on_double_click(self, index: QModelIndex) -> None:
        src_index = self._proxy.mapToSource(index)
        file_path = self._model.filePath(src_index)
        if not file_path:
            return
        p = Path(file_path)
        if p.is_dir():
            return

        name = p.name.lower()
        if (
            p.suffix.lower() in self._EXTERNAL_SUFFIXES
            or name.endswith((".fa.gz", ".fasta.gz", ".fna.gz"))
        ):
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(p)))
            return

        kind = self._classify(p)
        self.file_activated.emit(str(p), kind)

    def _classify(self, p: Path) -> str:
        """Map a file path to a kind hint for viewer routing."""
        name = p.name.lower()
        # Two-part suffixes
        if name.endswith((".fq.gz", ".fastq.gz")):
            return "fastq"
        suffix = p.suffix.lower()
        if suffix in self._TABLE_SUFFIXES:
            return "table"
        if suffix in self._FIGURE_SUFFIXES:
            return "figure"
        if suffix in self._FASTQ_SUFFIXES:
            return "fastq"
        if suffix in self._LOG_SUFFIXES:
            return "log"
        if suffix in self._TEXT_SUFFIXES:
            return "text"
        return "other"
