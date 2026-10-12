"""Top toolbar: New, Open, Save Config, Run, Stop, Force Stop.

Pure signal-emitting widget — the MainWindow wires each signal to the
appropriate handler.  The ``set_running`` method toggles between the
two operational modes (idle: Run enabled, Stop/Force Stop disabled /
running: vice versa).

Stop is cooperative (waits for the current step to finish).  Force Stop
signals the running subprocess group immediately and is gated by a
confirmation dialog in the MainWindow handler.
"""

from __future__ import annotations

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtGui import QKeySequence
from PyQt6.QtWidgets import QToolBar


class RunToolbar(QToolBar):
    """Main application toolbar with file and pipeline actions."""

    new_clicked = pyqtSignal()
    open_clicked = pyqtSignal()
    save_clicked = pyqtSignal()
    run_clicked = pyqtSignal()
    stop_clicked = pyqtSignal()
    force_quit_clicked = pyqtSignal()

    def __init__(self, parent=None) -> None:
        super().__init__("Run", parent)
        self.setMovable(False)

        self._new_action = self.addAction("New Job")
        self._new_action.setShortcut(QKeySequence.StandardKey.New)
        self._new_action.triggered.connect(self.new_clicked.emit)

        self._open_action = self.addAction("Open\u2026")
        self._open_action.setShortcut(QKeySequence.StandardKey.Open)
        self._open_action.triggered.connect(self.open_clicked.emit)

        self._save_action = self.addAction("Save Config")
        self._save_action.setShortcut(QKeySequence.StandardKey.Save)
        self._save_action.triggered.connect(self.save_clicked.emit)

        self.addSeparator()

        self._run_action = self.addAction("Run \u25b6")
        self._run_action.setShortcut(QKeySequence("Ctrl+R"))
        self._run_action.triggered.connect(self.run_clicked.emit)

        self._stop_action = self.addAction("Stop \u25a0")
        self._stop_action.setEnabled(False)
        self._stop_action.setToolTip(
            "Cancel after the current step finishes (preserves partial outputs)"
        )
        self._stop_action.triggered.connect(self.stop_clicked.emit)

        self._force_quit_action = self.addAction("Force Stop \u26a1")
        self._force_quit_action.setEnabled(False)
        self._force_quit_action.setToolTip(
            "Immediately kill the running step (with confirmation)"
        )
        self._force_quit_action.triggered.connect(self.force_quit_clicked.emit)

    def set_running(self, running: bool) -> None:
        """Toggle action enabled states between running and idle."""
        self._run_action.setEnabled(not running)
        self._stop_action.setEnabled(running)
        self._force_quit_action.setEnabled(running)
        # Disable file actions while a pipeline is in flight to prevent
        # the user from swapping the sample mid-run.
        self._new_action.setEnabled(not running)
        self._open_action.setEnabled(not running)
