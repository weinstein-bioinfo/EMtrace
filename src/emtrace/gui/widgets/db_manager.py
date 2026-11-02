"""Database Manager dialog — install and update EMtrace databases.

Three tabs:
- **Annotation databases** — minimap2/BLASTN/Diamond targets
- **Host genomes** — multi-file FASTA concatenated for minimap2
- **Taxonomy** — NCBI taxdump → SQLite; SILVA taxmap → SQLite

Each tab has its own table (Key / Status / Path / Source URL / Action),
Install All Missing button, Update All button, and progress bar.
Downloads run in :class:`DbInstallWorker` (background QThread) with
progress routed to the tab that initiated the download.
The DB-directory picker and the Close button are shared at dialog level.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QColor
from PyQt6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from emtrace.core.config import EMtraceConfig

# ---------------------------------------------------------------------------
# Column indices  (shared by all tables)
# ---------------------------------------------------------------------------
_COL_KEY    = 0
_COL_STATUS = 1
_COL_PATH   = 2
_COL_SOURCE = 3
_COL_ACTION = 4

_GREEN   = QColor("#16a34a")
_RED     = QColor("#dc2626")
_GREY    = QColor("#9ca3af")
_ENV_VAR = "EMTRACE_DB_DIR"


class DatabaseManagerDialog(QDialog):
    """Modal dialog for viewing and installing EMtrace databases.

    Signals
    -------
    config_updated
        Emitted after the DB directory is changed and config reloaded,
        so MainWindow can update its own ``_config`` reference.
    """

    config_updated = pyqtSignal(object)  # EMtraceConfig

    def __init__(
        self,
        config: EMtraceConfig,
        config_source: Optional[Path] = None,
        parent: Optional[QWidget] = None,
    ) -> None:
        super().__init__(parent)
        self._config = config
        self._config_source = config_source
        self._worker = None

        # key → QTableWidget; row looked up dynamically for sort safety
        self._key_table: Dict[str, QTableWidget] = {}
        # key → action button
        self._action_buttons: Dict[str, QPushButton] = {}
        # table → (last sorted column, last sort order)
        self._sort_state: Dict[QTableWidget, Tuple[int, Qt.SortOrder]] = {}

        # Progress widgets for the currently active download — set at launch
        self._active_progress: Optional[QProgressBar] = None
        self._active_status_label: Optional[QLabel] = None

        self.setWindowTitle("Database Manager")
        self.setMinimumWidth(1050)
        self.setMinimumHeight(600)

        self._build_ui()
        self._populate_tables()

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setSpacing(8)

        title = QLabel("<b>Database Manager</b>")
        title.setStyleSheet("font-size: 14px;")
        layout.addWidget(title)

        # Shared DB directory picker
        dir_row = QHBoxLayout()
        dir_row.addWidget(QLabel(f"{_ENV_VAR}:"))
        self._db_dir_edit = QLineEdit(os.environ.get(_ENV_VAR, ""))
        self._db_dir_edit.setReadOnly(True)
        self._db_dir_edit.setPlaceholderText("(not set — click Browse to choose)")
        self._style_db_dir()
        dir_row.addWidget(self._db_dir_edit, stretch=1)
        browse_btn = QPushButton("Browse\u2026")
        browse_btn.clicked.connect(self._on_browse)
        dir_row.addWidget(browse_btn)
        apply_btn = QPushButton("Apply")
        apply_btn.clicked.connect(self._on_apply)
        dir_row.addWidget(apply_btn)
        layout.addLayout(dir_row)

        # Tab widget — three tabs
        self._tabs = QTabWidget()
        self._ann_table, ann_tab = self._make_tab("Annotation databases")
        self._mm2_table, mm2_tab = self._make_tab("Host genomes")
        self._tax_table, tax_tab = self._make_tab("Taxonomy")

        # Per-tab controls
        (self._install_all_ann_btn, self._update_all_ann_btn,
         self._progress_ann, self._status_label_ann) = self._add_tab_controls(
            ann_tab, self._on_install_all_missing_ann, self._on_update_all_ann)

        (self._install_all_mm2_btn, self._update_all_mm2_btn,
         self._progress_mm2, self._status_label_mm2) = self._add_tab_controls(
            mm2_tab, self._on_install_all_missing_mm2, self._on_update_all_mm2)

        (self._install_all_tax_btn, self._update_all_tax_btn,
         self._progress_tax, self._status_label_tax) = self._add_tab_controls(
            tax_tab, self._on_install_all_missing_tax, self._on_update_all_tax)

        layout.addWidget(self._tabs, stretch=1)

        # Shared Close button
        close_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        close_box.rejected.connect(self.reject)
        layout.addWidget(close_box)

    def _make_tab(self, label: str) -> Tuple[QTableWidget, QWidget]:
        """Create a tab widget with a 5-column table; return (table, tab_widget)."""
        tab = QWidget()
        tab_layout = QVBoxLayout(tab)
        tab_layout.setContentsMargins(0, 4, 0, 0)

        table = QTableWidget(0, 5)
        table.setHorizontalHeaderLabels(["Key", "Status", "Path", "Source URL", "Action"])
        hdr = table.horizontalHeader()
        hdr.setSectionResizeMode(_COL_PATH,   QHeaderView.ResizeMode.Interactive)
        hdr.setSectionResizeMode(_COL_SOURCE, QHeaderView.ResizeMode.Interactive)
        hdr.setStretchLastSection(False)
        table.setSortingEnabled(False)
        hdr.setSectionsClickable(True)
        hdr.sectionClicked.connect(
            lambda col, t=table: self._on_header_clicked(t, col)
        )
        table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        table.setSelectionMode(QTableWidget.SelectionMode.NoSelection)
        table.verticalHeader().setVisible(False)
        table.setAlternatingRowColors(True)
        tab_layout.addWidget(table, stretch=1)

        self._tabs.addTab(tab, label)
        return table, tab

    def _add_tab_controls(
        self,
        tab: QWidget,
        install_slot,
        update_slot,
    ) -> Tuple[QPushButton, QPushButton, QProgressBar, QLabel]:
        """Add Install All Missing + Update All + progress row to *tab*."""
        tab_layout = tab.layout()

        btn_row = QHBoxLayout()
        install_btn = QPushButton("Install All Missing")
        install_btn.clicked.connect(install_slot)
        btn_row.addWidget(install_btn)
        update_btn = QPushButton("Update All")
        update_btn.clicked.connect(update_slot)
        btn_row.addWidget(update_btn)
        btn_row.addStretch()
        tab_layout.addLayout(btn_row)

        prog_row = QHBoxLayout()
        progress = QProgressBar()
        progress.setRange(0, 100)
        progress.setValue(0)
        progress.setVisible(False)
        status_label = QLabel("")
        prog_row.addWidget(progress, stretch=1)
        prog_row.addWidget(status_label)
        tab_layout.addLayout(prog_row)

        return install_btn, update_btn, progress, status_label

    # ------------------------------------------------------------------
    # DB directory
    # ------------------------------------------------------------------

    def _style_db_dir(self) -> None:
        color = "#16a34a" if self._db_dir_edit.text() else "#9ca3af"
        self._db_dir_edit.setStyleSheet(f"color: {color};")
        self._db_dir_edit.setToolTip(self._db_dir_edit.text())

    def _on_browse(self) -> None:
        current = self._db_dir_edit.text() or str(Path.home())
        chosen = QFileDialog.getExistingDirectory(
            self, "Select Database Directory", current
        )
        if chosen:
            self._db_dir_edit.setText(chosen)
            self._style_db_dir()

    def _on_apply(self) -> None:
        path = self._db_dir_edit.text().strip()
        if not path:
            QMessageBox.warning(self, "No Directory", "Browse to a directory first.")
            return
        if not Path(path).is_dir():
            QMessageBox.warning(
                self, "Invalid Directory",
                f"'{path}' does not exist or is not a directory.",
            )
            return

        os.environ[_ENV_VAR] = path

        from emtrace.utils.user_config import load_user_config, save_user_config
        cfg = load_user_config()
        cfg["db_dir"] = path
        save_user_config(cfg)

        if self._config_source and self._config_source.exists():
            self._config = EMtraceConfig.from_yaml(self._config_source)
        else:
            self._config = EMtraceConfig.create_default()

        self.config_updated.emit(self._config)
        self._style_db_dir()
        self._clear_tables()
        self._populate_tables()

    def _fmt_path(self, path: str) -> str:
        """Abbreviate the DB dir prefix to $EMTRACE_DB_DIR for display."""
        db_dir = os.environ.get(_ENV_VAR, "")
        if db_dir and path.startswith(db_dir):
            rel = path[len(db_dir):].lstrip("/")
            return f"${_ENV_VAR}/{rel}"
        return path or "(not configured)"

    # ------------------------------------------------------------------
    # Table population
    # ------------------------------------------------------------------

    def _clear_tables(self) -> None:
        for table in (self._ann_table, self._mm2_table, self._tax_table):
            table.clearContents()
            table.setRowCount(0)
        self._key_table.clear()
        self._action_buttons.clear()

    def _populate_tables(self) -> None:
        from emtrace.utils.db_download import database_status, taxonomy_status

        ann_rows, host_rows = database_status(self._config)
        tax_rows = taxonomy_status(self._config)

        self._fill_table(self._ann_table, ann_rows)
        self._fill_table(self._mm2_table, host_rows)
        self._fill_table(self._tax_table, tax_rows)

    def _fill_table(self, table: QTableWidget, rows: List[dict]) -> None:
        table.setRowCount(len(rows))

        for row_idx, row in enumerate(rows):
            key = row["key"]
            self._key_table[key] = table

            # Key
            table.setItem(row_idx, _COL_KEY, QTableWidgetItem(key))

            # Status
            self._set_status_cell(key, row["installed"])

            # Path — abbreviated display, full path in tooltip
            raw_path = row.get("path", "")
            path_display = self._fmt_path(raw_path)
            path_item = QTableWidgetItem(path_display)
            if not raw_path:
                path_item.setForeground(_GREY)
            else:
                path_item.setToolTip(raw_path)
            table.setItem(row_idx, _COL_PATH, path_item)

            # Source URL — handle single URL, multi-URL (minimap2), and none
            urls: List[str] = row.get("urls", [])
            url: str = row.get("url", "") or ""
            if len(urls) > 1:
                src_item = QTableWidgetItem(f"{len(urls)} files")
                src_item.setToolTip("\n".join(urls))
            elif url:
                src_item = QTableWidgetItem(url)
                src_item.setToolTip(url)
            else:
                src_item = QTableWidgetItem("(manual install only)")
                src_item.setForeground(_GREY)
            table.setItem(row_idx, _COL_SOURCE, src_item)

            # Action button — shown when any URL is configured
            has_url = bool(urls or url)
            if has_url:
                btn_url = urls[0] if urls else url   # worker uses install_by_key; URL is display only
                label = "Update" if row["installed"] else "Install"
                btn = QPushButton(label)
                btn.clicked.connect(
                    lambda checked, k=key, u=btn_url, p=row.get("path", ""):
                        self._on_install_single(k, u, p)
                )
                self._action_buttons[key] = btn
                table.setCellWidget(row_idx, _COL_ACTION, btn)

        table.resizeColumnToContents(_COL_KEY)
        table.resizeColumnToContents(_COL_STATUS)
        table.resizeColumnToContents(_COL_ACTION)
        table.setColumnWidth(_COL_PATH,   300)
        table.setColumnWidth(_COL_SOURCE, 300)

    def _find_row(self, key: str) -> Optional[int]:
        """Return the current visual row index for *key*, or None."""
        table = self._key_table.get(key)
        if table is None:
            return None
        matches = table.findItems(key, Qt.MatchFlag.MatchExactly)
        for item in matches:
            if item.column() == _COL_KEY:
                return item.row()
        return None

    def _on_header_clicked(self, table: QTableWidget, col: int) -> None:
        """Sort by *col* if it is Key or Status; ignore other columns."""
        if col not in (_COL_KEY, _COL_STATUS):
            return
        last_col, last_order = self._sort_state.get(
            table, (-1, Qt.SortOrder.AscendingOrder)
        )
        if last_col == col:
            order = (
                Qt.SortOrder.DescendingOrder
                if last_order == Qt.SortOrder.AscendingOrder
                else Qt.SortOrder.AscendingOrder
            )
        else:
            order = Qt.SortOrder.AscendingOrder
        self._sort_state[table] = (col, order)
        table.sortItems(col, order)
        hdr = table.horizontalHeader()
        hdr.setSortIndicator(col, order)
        hdr.setSortIndicatorShown(True)

    def _set_status_cell(self, key: str, installed: bool) -> None:
        table   = self._key_table.get(key)
        row_idx = self._find_row(key)
        if table is None or row_idx is None:
            return
        item = QTableWidgetItem("\u2713 Installed" if installed else "\u2717 Missing")
        item.setForeground(_GREEN if installed else _RED)
        table.setItem(row_idx, _COL_STATUS, item)

    # ------------------------------------------------------------------
    # Helpers for install / update tab actions
    # ------------------------------------------------------------------

    def _get_rows_for_tab(self, tab: str) -> List[dict]:
        """Return status rows for *tab* ('ann', 'mm2', or 'tax')."""
        from emtrace.utils.db_download import database_status, taxonomy_status
        ann_rows, host_rows = database_status(self._config)
        tax_rows = taxonomy_status(self._config)
        return {"ann": ann_rows, "mm2": host_rows, "tax": tax_rows}.get(tab, [])

    def _tab_label(self, tab: str) -> str:
        return {
            "ann": "annotation",
            "mm2": "host genome",
            "tax": "taxonomy",
        }.get(tab, tab)

    def _tab_progress(self, tab: str) -> Tuple[QProgressBar, QLabel]:
        return {
            "ann": (self._progress_ann, self._status_label_ann),
            "mm2": (self._progress_mm2, self._status_label_mm2),
            "tax": (self._progress_tax, self._status_label_tax),
        }[tab]

    # ------------------------------------------------------------------
    # Install actions — per tab
    # ------------------------------------------------------------------

    def _on_install_all_missing_tax(self) -> None:
        self._install_all_missing_for_tab("tax")

    def _on_update_all_tax(self) -> None:
        self._update_all_for_tab("tax")

    def _on_install_all_missing_ann(self) -> None:
        self._install_all_missing_for_tab("ann")

    def _on_install_all_missing_mm2(self) -> None:
        self._install_all_missing_for_tab("mm2")

    def _on_update_all_ann(self) -> None:
        self._update_all_for_tab("ann")

    def _on_update_all_mm2(self) -> None:
        self._update_all_for_tab("mm2")

    def _update_all_for_tab(self, tab: str) -> None:
        if self._worker and self._worker.isRunning():
            return
        if not self._check_db_dir():
            return

        source_rows = self._get_rows_for_tab(tab)
        installs = [
            (r["key"], r.get("url", ""), r.get("path", ""))
            for r in source_rows
            if r.get("url") or r.get("urls")
        ]

        if not installs:
            QMessageBox.information(
                self, "Nothing to update",
                "No databases with a source URL are configured.",
            )
            return

        label = self._tab_label(tab)
        reply = QMessageBox.question(
            self, "Confirm Update All",
            f"Re-download and overwrite {len(installs)} {label} database(s)?\n"
            "This may download several gigabytes of data.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
        )
        if reply != QMessageBox.StandardButton.Yes:
            return

        progress, status_label = self._tab_progress(tab)
        self._launch_worker(installs, force=True,
                            progress=progress, status_label=status_label)

    def _install_all_missing_for_tab(self, tab: str) -> None:
        if self._worker and self._worker.isRunning():
            return
        if not self._check_db_dir():
            return

        source_rows = self._get_rows_for_tab(tab)
        installs = [
            (r["key"], r.get("url", ""), r.get("path", ""))
            for r in source_rows
            if not r["installed"] and (r.get("url") or r.get("urls"))
        ]

        if not installs:
            QMessageBox.information(
                self, "Nothing to install",
                "All databases with a source URL are already installed.",
            )
            return

        label = self._tab_label(tab)
        reply = QMessageBox.question(
            self, "Confirm Install",
            f"Install {len(installs)} missing {label} database(s)?\n"
            "This may download several gigabytes of data.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
        )
        if reply != QMessageBox.StandardButton.Yes:
            return

        progress, status_label = self._tab_progress(tab)
        self._launch_worker(installs, force=False,
                            progress=progress, status_label=status_label)

    def _on_install_single(self, key: str, url: str, path: str) -> None:
        if self._worker and self._worker.isRunning():
            return
        if not self._check_db_dir():
            return

        table   = self._key_table.get(key)
        row_idx = self._find_row(key)
        installed = False
        if table is not None and row_idx is not None:
            status_item = table.item(row_idx, _COL_STATUS)
            installed = status_item is not None and "\u2713" in status_item.text()

        force = False
        if installed:
            reply = QMessageBox.question(
                self, "Re-download?",
                f"'{key}' is already installed.\nRe-download and overwrite?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            )
            if reply != QMessageBox.StandardButton.Yes:
                return
            force = True

        if table is self._ann_table:
            progress, status_label = self._progress_ann, self._status_label_ann
        elif table is self._mm2_table:
            progress, status_label = self._progress_mm2, self._status_label_mm2
        else:
            progress, status_label = self._progress_tax, self._status_label_tax

        self._launch_worker([(key, url, path)], force=force,
                            progress=progress, status_label=status_label)

    def _check_db_dir(self) -> bool:
        """Return True if EMTRACE_DB_DIR is set and is a directory."""
        db_dir = os.environ.get(_ENV_VAR, "").strip()
        if not db_dir or not Path(db_dir).is_dir():
            QMessageBox.warning(
                self, "DB Directory Not Set",
                f"Please set {_ENV_VAR} to a valid directory using the "
                "Browse + Apply controls above before installing databases.",
            )
            return False
        return True

    def _launch_worker(
        self,
        installs: list,
        force: bool,
        progress: QProgressBar,
        status_label: QLabel,
    ) -> None:
        from emtrace.gui.workers.db_worker import DbInstallWorker

        if self._worker and self._worker.isRunning():
            return

        self._active_progress     = progress
        self._active_status_label = status_label

        self._set_buttons_enabled(False)
        progress.setVisible(True)
        progress.setValue(0)
        status_label.setText("Starting\u2026")

        self._worker = DbInstallWorker(self._config, installs, force=force, parent=self)
        self._worker.signals.progress.connect(self._on_progress)
        self._worker.signals.db_finished.connect(self._on_db_finished)
        self._worker.signals.all_finished.connect(self._on_all_finished)
        self._worker.signals.error.connect(self._on_error)
        self._worker.start()

    # ------------------------------------------------------------------
    # Worker signal handlers
    # ------------------------------------------------------------------

    def _on_progress(self, key: str, done: int, total: int) -> None:
        progress     = self._active_progress
        status_label = self._active_status_label
        if progress is None or status_label is None:
            return
        if total > 0:
            progress.setRange(0, 100)
            progress.setValue(int(done * 100 / total))
            status_label.setText(f"{key}: {done / 1e6:.1f} / {total / 1e6:.1f} MB")
        else:
            progress.setRange(0, 0)
            status_label.setText(
                f"{key}: {done / 1e6:.1f} MB" if done else f"{key}: connecting\u2026"
            )

    def _on_db_finished(self, key: str, success: bool, message: str) -> None:
        if key in self._key_table:
            from emtrace.utils.db_download import database_status, taxonomy_status
            ann_rows, host_rows = database_status(self._config)
            tax_rows = taxonomy_status(self._config)
            all_rows = {r["key"]: r for r in ann_rows + host_rows + tax_rows}
            installed = all_rows.get(key, {}).get("installed", False)
            self._set_status_cell(key, installed)
            btn = self._action_buttons.get(key)
            if btn:
                btn.setText("Update" if installed else "Install")
                btn.setEnabled(True)

        if not success:
            QMessageBox.warning(self, "Install Failed", f"{key}: {message}")

    def _on_all_finished(self, installed: int, skipped: int, failed: int) -> None:
        progress     = self._active_progress
        status_label = self._active_status_label
        if progress is not None:
            progress.setRange(0, 100)
            progress.setValue(100 if failed == 0 else 0)
            progress.setVisible(False)
        if status_label is not None:
            status_label.setText(
                f"Done: {installed} installed, {skipped} skipped, {failed} failed"
            )
        self._set_buttons_enabled(True)

    def _on_error(self, message: str) -> None:
        QMessageBox.warning(self, "Error", message)
        self._set_buttons_enabled(True)

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _set_buttons_enabled(self, enabled: bool) -> None:
        for btn in (
            self._install_all_ann_btn, self._update_all_ann_btn,
            self._install_all_mm2_btn, self._update_all_mm2_btn,
            self._install_all_tax_btn, self._update_all_tax_btn,
        ):
            btn.setEnabled(enabled)
        for btn in self._action_buttons.values():
            btn.setEnabled(enabled)
