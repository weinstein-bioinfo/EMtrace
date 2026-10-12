"""Abundance tab container — sub-tabs for stacked bar and sunburst views."""

from __future__ import annotations

from pathlib import Path

from PyQt6.QtWidgets import (
    QStackedWidget,
    QTabBar,
    QVBoxLayout,
    QWidget,
)

from emtrace.gui.widgets.abundance_plot import AbundancePlotWidget
from emtrace.gui.widgets.sunburst_plot import SunburstPlotWidget


class AbundanceTabWidget(QWidget):
    """Wraps the bar-chart and sunburst widgets behind a small sub-tab bar."""

    def __init__(self, parent=None):
        super().__init__(parent)

        self.bar = AbundancePlotWidget()
        self.sunburst = SunburstPlotWidget()

        self._tab_bar = QTabBar()
        self._tab_bar.addTab("Stacked bar")
        self._tab_bar.addTab("Sunburst")

        self._stack = QStackedWidget()
        self._stack.addWidget(self.bar)
        self._stack.addWidget(self.sunburst)

        self._tab_bar.currentChanged.connect(self._stack.setCurrentIndex)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self._tab_bar)
        layout.addWidget(self._stack, stretch=1)

    # Forward the public surface used by MainContentStack
    def load_sample(self, output_dir: Path) -> None:
        self.bar.load_sample(output_dir)
        self.sunburst.load_sample(output_dir)

    def refresh(self) -> None:
        self.bar.refresh()
        self.sunburst.refresh()
