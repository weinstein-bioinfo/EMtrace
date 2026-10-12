"""Bottom panel — color-coded log viewer with level filter and auto-scroll.

Connects to :class:`emtrace.gui.models.log_handler.QtLogHandler` so the
pipeline (running on a worker thread) can stream records here without
direct cross-thread mutation.
"""

from __future__ import annotations

import logging
from datetime import datetime
from html import escape
from pathlib import Path
from typing import List

from PyQt6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QPlainTextEdit,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from emtrace.gui.platform import monospace_font


class LogPanel(QWidget):
    """Color-coded log viewer with level filter, auto-scroll, and Save."""

    _LEVEL_COLORS = {
        "DEBUG": "#888888",
        "INFO": "#1a7f37",
        "WARNING": "#9a6700",
        "ERROR": "#cf222e",
        "CRITICAL": "#cf222e",
    }

    _LEVEL_VALUES = {
        "DEBUG": logging.DEBUG,
        "INFO": logging.INFO,
        "WARNING": logging.WARNING,
        "ERROR": logging.ERROR,
        "CRITICAL": logging.CRITICAL,
    }

    _MAX_RECORDS = 10_000

    def __init__(self, parent=None) -> None:
        super().__init__(parent)

        # Ring buffer of received records (used to re-render on filter change)
        self._records: List[logging.LogRecord] = []

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        # Read-only text area
        self._text = QPlainTextEdit()
        self._text.setReadOnly(True)
        self._text.setMaximumBlockCount(self._MAX_RECORDS)
        self._text.setFont(monospace_font(10))
        layout.addWidget(self._text, stretch=1)

        # Controls
        controls = QHBoxLayout()

        self._clear_btn = QPushButton("Clear")
        self._clear_btn.clicked.connect(self._on_clear)
        controls.addWidget(self._clear_btn)

        self._save_btn = QPushButton("Save Log")
        self._save_btn.clicked.connect(self._on_save)
        controls.addWidget(self._save_btn)

        controls.addStretch()

        self._level_combo = QComboBox()
        self._level_combo.addItems(["DEBUG", "INFO", "WARNING", "ERROR"])
        self._level_combo.setCurrentText("INFO")
        self._level_combo.currentTextChanged.connect(self._refilter)
        controls.addWidget(self._level_combo)

        self._auto_scroll = QCheckBox("Auto-scroll")
        self._auto_scroll.setChecked(True)
        controls.addWidget(self._auto_scroll)

        layout.addLayout(controls)

    # ==================================================================
    # Public slot — connect to QtLogHandler.signals.log_record
    # ==================================================================

    def append_record(self, record: logging.LogRecord) -> None:
        """Store *record* in the ring buffer and render it if it passes the filter."""
        self._records.append(record)
        if len(self._records) > self._MAX_RECORDS:
            # Drop oldest in chunks of 1k to keep amortised cost low
            self._records = self._records[-self._MAX_RECORDS:]
        if self._should_show(record.levelname):
            self._render(record)

    # ==================================================================
    # Rendering
    # ==================================================================

    def _render(self, record: logging.LogRecord) -> None:
        level = record.levelname
        timestamp = datetime.fromtimestamp(record.created).strftime("%H:%M:%S")
        # Strip the ``emtrace.`` prefix to keep the logger column compact
        name = record.name
        if name.startswith("emtrace."):
            name = name[len("emtrace."):]
        color = self._LEVEL_COLORS.get(level, "#000000")
        message = escape(record.getMessage())
        html = (
            f'<span style="color:{color}">'
            f"{timestamp} | {level:8s} | {escape(name)} | {message}"
            f"</span>"
        )
        self._text.appendHtml(html)
        if self._auto_scroll.isChecked():
            self._text.ensureCursorVisible()

    def _refilter(self) -> None:
        """Re-render the buffered records using the current level filter."""
        self._text.clear()
        for record in self._records:
            if self._should_show(record.levelname):
                self._render(record)

    def _should_show(self, level_name: str) -> bool:
        threshold = self._LEVEL_VALUES.get(
            self._level_combo.currentText(), logging.INFO,
        )
        record_level = self._LEVEL_VALUES.get(level_name, logging.DEBUG)
        return record_level >= threshold

    # ==================================================================
    # Buttons
    # ==================================================================

    def _on_clear(self) -> None:
        """Clear both the on-screen text and the ring buffer."""
        self._records.clear()
        self._text.clear()

    def _on_save(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Log", "emtrace.log", "Log Files (*.log *.txt)",
        )
        if path:
            Path(path).write_text(
                self._text.toPlainText(), encoding="utf-8", newline="\n",
            )
