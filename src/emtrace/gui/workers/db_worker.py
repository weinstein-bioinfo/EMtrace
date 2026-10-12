"""QThread worker for database download and installation."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional, Tuple

from PyQt6.QtCore import QObject, QThread, pyqtSignal

from emtrace.core.config import EMtraceConfig

logger = logging.getLogger(__name__)

# (key, url, dest_path_str)
InstallRecord = Tuple[str, str, str]


class DbWorkerSignals(QObject):
    """Signals emitted by :class:`DbInstallWorker`."""

    # (key, bytes_done, bytes_total)  — bytes_total=-1 when unknown
    # Use 'object' for byte counts: C++ int is 32-bit, overflows above ~2 GB.
    progress = pyqtSignal(str, object, object)

    # (key, success, message)
    db_finished = pyqtSignal(str, bool, str)

    # (installed_count, skipped_count, failed_count)
    all_finished = pyqtSignal(int, int, int)

    error = pyqtSignal(str)


class DbInstallWorker(QThread):
    """Install one or more databases in a background thread.

    Usage::

        worker = DbInstallWorker(config, installs, force=False)
        worker.signals.progress.connect(on_progress)
        worker.signals.db_finished.connect(on_db_done)
        worker.signals.all_finished.connect(on_all_done)
        worker.start()
    """

    def __init__(
        self,
        config: EMtraceConfig,
        installs: List[InstallRecord],
        force: bool = False,
        parent: Optional[QObject] = None,
    ) -> None:
        super().__init__(parent)
        self.signals = DbWorkerSignals()
        self._config = config
        self._installs = installs
        self._force = force

    def run(self) -> None:
        from emtrace.utils.db_download import (
            DatabaseAlreadyInstalledError,
            DatabaseInstallError,
            DatabaseNotConfiguredError,
            install_by_key,
        )

        installed = skipped = failed = 0

        for key, _url, _path in self._installs:
            def _cb(done: int, total: int, _k: str = key) -> None:
                self.signals.progress.emit(_k, done, total)

            try:
                install_by_key(
                    self._config,
                    key,
                    progress_cb=_cb,
                    force=self._force,
                )
                self.signals.db_finished.emit(key, True, "installed")
                installed += 1

            except DatabaseAlreadyInstalledError:
                self.signals.db_finished.emit(key, True, "skipped")
                skipped += 1

            except (DatabaseNotConfiguredError, DatabaseInstallError, Exception) as exc:  # noqa: BLE001
                logger.error("Failed to install %s: %s", key, exc)
                self.signals.db_finished.emit(key, False, str(exc))
                failed += 1

        self.signals.all_finished.emit(installed, skipped, failed)
