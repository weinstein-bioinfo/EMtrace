"""Table Viewer widget — CSV/TSV file selector + searchable QTableView.

Adapted from ViRefSet ``table_viewer.py``.  Scans the sample output
directory tree for ``*.csv`` and ``*.tsv`` files, displays them in a
sortable + filterable table, and supports export.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd
from PyQt6.QtCore import QSortFilterProxyModel, Qt
from PyQt6.QtWidgets import (
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QPushButton,
    QTableView,
    QVBoxLayout,
    QWidget,
)

from emtrace.gui.models.pandas_model import PandasTableModel

# File patterns to discover (relative to sample output_dir)
_TABLE_PATTERNS = ("*.csv", "*.tsv")


class _TableSortProxy(QSortFilterProxyModel):
    """Proxy that sorts numeric columns numerically.

    Column 0 is the virtual '#' row-index column (always integer).
    All other columns attempt float conversion first; only if that
    fails for either value does the comparison fall back to the
    default case-insensitive string sort.
    """

    def lessThan(self, left, right):  # noqa: N802
        lv, rv = left.data(), right.data()
        try:
            return float(lv) < float(rv)
        except (ValueError, TypeError):
            pass
        return super().lessThan(left, right)


class TableViewerWidget(QWidget):
    """Interactive CSV / TSV table viewer.

    Top row: file selector combo + filter text field + row count + export.
    Body: sortable, filterable QTableView backed by PandasTableModel.
    """

    def __init__(self, parent=None):
        super().__init__(parent)

        # ── File selector ────────────────────────────────────────────
        self._file_combo = QComboBox()
        self._file_combo.setSizeAdjustPolicy(
            QComboBox.SizeAdjustPolicy.AdjustToContents
        )
        self._file_combo.setMinimumWidth(200)

        # ── Search / filter ──────────────────────────────────────────
        self._search = QLineEdit()
        self._search.setPlaceholderText("Filter rows\u2026")
        self._search.setClearButtonEnabled(True)

        # ── Row count label ──────────────────────────────────────────
        self._row_label = QLabel()

        # ── Export button ────────────────────────────────────────────
        self._export_btn = QPushButton("Export\u2026")
        self._export_btn.setFixedWidth(80)

        # ── Table view ───────────────────────────────────────────────
        self._model = PandasTableModel()
        self._proxy = _TableSortProxy()
        self._proxy.setSourceModel(self._model)
        self._proxy.setFilterCaseSensitivity(
            Qt.CaseSensitivity.CaseInsensitive
        )
        self._proxy.setFilterKeyColumn(-1)  # search all columns

        self._table = QTableView()
        self._table.setModel(self._proxy)
        self._table.setSortingEnabled(True)
        self._table.setAlternatingRowColors(True)
        self._table.setSelectionBehavior(
            QTableView.SelectionBehavior.SelectRows
        )
        self._table.horizontalHeader().setStretchLastSection(True)
        self._table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.Interactive
        )
        self._table.verticalHeader().setDefaultSectionSize(22)

        # ── Empty-state label ────────────────────────────────────────
        self._empty_label = QLabel(
            "No result tables found.\n"
            "Run pipeline steps to generate output."
        )
        self._empty_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._empty_label.setStyleSheet("color: #999; font-size: 14px;")

        # ── Layout ───────────────────────────────────────────────────
        top_row = QHBoxLayout()
        top_row.setContentsMargins(4, 2, 4, 2)
        top_row.addWidget(QLabel("File:"))
        top_row.addWidget(self._file_combo, stretch=1)
        top_row.addWidget(self._search, stretch=1)
        top_row.addWidget(self._row_label)
        top_row.addWidget(self._export_btn)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addLayout(top_row)
        layout.addWidget(self._table, stretch=1)
        layout.addWidget(self._empty_label, stretch=1)

        self._set_has_data(False)

        # ── Connections ──────────────────────────────────────────────
        self._file_combo.currentIndexChanged.connect(self._on_file_changed)
        self._search.textChanged.connect(self._on_filter_changed)
        self._export_btn.clicked.connect(self._on_export)
        self._proxy.rowsInserted.connect(self._update_row_count)
        self._proxy.rowsRemoved.connect(self._update_row_count)
        self._proxy.modelReset.connect(self._update_row_count)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def load_sample_tables(self, output_dir: Path) -> None:
        """Scan sample output directory for CSV/TSV files and populate combo."""
        self._file_combo.blockSignals(True)
        self._file_combo.clear()

        paths: list[Path] = []
        for pat in _TABLE_PATTERNS:
            paths.extend(sorted(output_dir.rglob(pat)))

        for path in paths:
            try:
                if path.stat().st_size == 0:
                    continue
            except OSError:
                continue
            rel = path.relative_to(output_dir)
            self._file_combo.addItem(str(rel), str(path))

        self._file_combo.blockSignals(False)

        has_files = self._file_combo.count() > 0
        self._set_has_data(has_files)
        if has_files:
            self._file_combo.setCurrentIndex(0)
            self._on_file_changed(0)

    def load_file(self, path: Path) -> None:
        """Load a specific CSV/TSV file into the viewer."""
        idx = self._file_combo.findData(str(path))
        if idx >= 0:
            self._file_combo.setCurrentIndex(idx)
        else:
            self._load_csv(path)

    def refresh(self) -> None:
        """Re-scan the current output_dir for new result files."""
        if self._file_combo.count() == 0:
            return
        # Recover output_dir from first item's absolute path
        first_abs = self._file_combo.itemData(0)
        first_rel = self._file_combo.itemText(0)
        if first_abs and first_rel:
            output_dir = Path(first_abs[: -len(first_rel)])
            self.load_sample_tables(output_dir)

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _set_has_data(self, has: bool) -> None:
        self._file_combo.setVisible(has)
        self._search.setVisible(has)
        self._row_label.setVisible(has)
        self._export_btn.setVisible(has)
        self._table.setVisible(has)
        self._empty_label.setVisible(not has)

    def _on_file_changed(self, index: int) -> None:
        if index < 0:
            return
        path = Path(self._file_combo.itemData(index))
        self._load_csv(path)

    def _load_csv(self, path: Path) -> None:
        """Read a CSV/TSV into the model."""
        if not path.exists():
            return
        sep = "\t" if path.suffix == ".tsv" else ","
        try:
            df = pd.read_csv(
                path, sep=sep, keep_default_na=False, na_values=[""], comment="#"
            )
        except (pd.errors.EmptyDataError, pd.errors.ParserError, ValueError):
            df = pd.DataFrame()

        self._table.horizontalHeader().setSortIndicator(
            -1, Qt.SortOrder.AscendingOrder
        )
        self._model.set_dataframe(df)
        self._table.resizeColumnsToContents()
        self._search.clear()
        self._update_row_count()

    def _on_filter_changed(self, text: str) -> None:
        self._proxy.setFilterFixedString(text)
        self._update_row_count()

    def _update_row_count(self) -> None:
        shown = self._proxy.rowCount()
        total = self._model.rowCount()
        if shown == total:
            self._row_label.setText(f"{total} rows")
        else:
            self._row_label.setText(f"{shown}/{total} rows")

    def _on_export(self) -> None:
        """Export visible (filtered) rows to a new CSV file."""
        if self._proxy.rowCount() == 0:
            return

        save_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Table",
            "",
            "CSV Files (*.csv);;TSV Files (*.tsv)",
        )
        if not save_path:
            return

        rows: list[int] = []
        for proxy_row in range(self._proxy.rowCount()):
            source_idx = self._proxy.mapToSource(
                self._proxy.index(proxy_row, 0)
            )
            rows.append(source_idx.row())

        df = self._model.dataframe.iloc[rows].copy()
        sep = "\t" if save_path.endswith(".tsv") else ","
        df.to_csv(save_path, sep=sep, index=False)
