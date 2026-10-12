"""Bridges Python logging → Qt signal for thread-safe GUI updates.

The pipeline and worker threads use the standard ``logging`` module.
``QtLogHandler`` converts each ``LogRecord`` into a Qt signal so it can
be consumed by a widget on the GUI main thread without cross-thread
state mutation.
"""

from __future__ import annotations

import logging

from PyQt6.QtCore import QObject, pyqtSignal


class QtLogHandlerSignals(QObject):
    """Signals emitted by :class:`QtLogHandler`."""

    log_record = pyqtSignal(object)  # logging.LogRecord


class QtLogHandler(logging.Handler):
    """Logging handler that re-emits records as a Qt signal.

    Usage::

        handler = QtLogHandler()
        logging.getLogger("emtrace").addHandler(handler)
        handler.signals.log_record.connect(log_panel.append_record)

    The signal queue is thread-safe: worker threads may log freely and
    the slot is invoked on the receiver's thread (typically the GUI
    main thread).
    """

    def __init__(self, level: int = logging.DEBUG) -> None:
        super().__init__(level)
        self.signals = QtLogHandlerSignals()

    def emit(self, record: logging.LogRecord) -> None:
        self.signals.log_record.emit(record)
