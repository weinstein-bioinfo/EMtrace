"""Centre panel — :class:`QStackedWidget` with a top tab bar.

Pages:
  0 Dashboard  — sample overview + step status
  1 Config     — YAML config editor
  2 Tables     — CSV/TSV result table viewer
  3 Figures    — static PNG/SVG/PDF figure viewer
  4 Abundance  — embedded matplotlib abundance bar chart
  5 QC         — read length + quality histograms
  6 Files      — output directory tree browser
  7 Assembly   — coverage strip, variant table, segment summary (disabled until
                 assembly step produces consensus.fa)
"""

from __future__ import annotations

from typing import Optional

from PyQt6.QtWidgets import (
    QStackedWidget,
    QTabBar,
    QVBoxLayout,
    QWidget,
)

from emtrace.core.config import EMtraceConfig
from emtrace.core.pipelinejob import PipelineJob
from emtrace.gui.widgets.abundance_tab import AbundanceTabWidget
from emtrace.gui.widgets.assembly_view import AssemblyViewWidget
from emtrace.gui.widgets.config_editor import ConfigEditorWidget
from emtrace.gui.widgets.dashboard import DashboardWidget
from emtrace.gui.widgets.figure_viewer import FigureViewerWidget
from emtrace.gui.widgets.file_browser import FileBrowserWidget
from emtrace.gui.widgets.qc_plot import QCPlotWidget
from emtrace.gui.widgets.realtime_panel import RealtimePanelWidget
from emtrace.gui.widgets.table_viewer import TableViewerWidget


class MainContentStack(QWidget):
    """Tab bar + stacked pages for all content views."""

    TAB_DASHBOARD = 0
    TAB_CONFIG = 1
    TAB_ABUNDANCE = 2
    TAB_ASSEMBLY = 3
    TAB_FIGURES = 4
    TAB_FILES = 5
    TAB_QC = 6
    TAB_REALTIME = 7
    TAB_TABLES = 8

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self._tab_bar = QTabBar()
        self._stack = QStackedWidget()

        # Pages — keep names public so MainWindow can wire signals.
        self.dashboard = DashboardWidget()
        self.config_page = ConfigEditorWidget()
        self.table_viewer = TableViewerWidget()
        self.figure_viewer = FigureViewerWidget()
        self.abundance_plot = AbundanceTabWidget()
        self.qc_plot = QCPlotWidget()
        self.file_browser = FileBrowserWidget()
        self.assembly_view = AssemblyViewWidget()
        self.realtime_panel = RealtimePanelWidget()

        for label, widget in (
            ("Dashboard", self.dashboard),
            ("Config", self.config_page),
            ("Abundance", self.abundance_plot),
            ("Assembly", self.assembly_view),
            ("Figures", self.figure_viewer),
            ("Files", self.file_browser),
            ("QC", self.qc_plot),
            ("Realtime", self.realtime_panel),
            ("Tables", self.table_viewer),
        ):
            self._tab_bar.addTab(label)
            self._stack.addWidget(widget)

        self._tab_bar.currentChanged.connect(self._stack.setCurrentIndex)

        layout.addWidget(self._tab_bar)
        layout.addWidget(self._stack, stretch=1)

    # ==================================================================
    # Public helpers
    # ==================================================================

    def switch_to(self, index: int) -> None:
        """Programmatically switch tabs (use ``TAB_*`` constants)."""
        if 0 <= index < self._tab_bar.count():
            self._tab_bar.setCurrentIndex(index)

    def set_job(self, job: Optional[PipelineJob]) -> None:
        """Propagate a new (or cleared) job to all sub-pages."""
        self.dashboard.set_job(job)
        self.assembly_view.set_job(job)
        self.realtime_panel.set_job(job)

        if job is not None:
            out = job.output_dir
            self.file_browser.set_root(out)
            self.table_viewer.load_sample_tables(out)
            self.figure_viewer.load_sample_figures(out)
            self.abundance_plot.load_sample(out)
            self.qc_plot.set_job_files(job.files, output_dir=job.output_dir)
            self.config_page.set_output_dir(out)
        else:
            self.file_browser.set_root(None)
            self.config_page.set_output_dir(None)

        self._update_assembly_tab_state(job)

    def set_config(self, config: EMtraceConfig) -> None:
        """Load *config* into the Config editor tab."""
        self.config_page.load_config(config)

    def refresh_results(self, job: Optional[PipelineJob] = None) -> None:
        """Re-scan output directory for new result files and plots.

        Called after pipeline completion to pick up newly generated outputs.
        """
        if job is not None:
            self.set_job(job)
        else:
            self.table_viewer.refresh()
            self.abundance_plot.refresh()

    def _update_assembly_tab_state(self, job: Optional[PipelineJob]) -> None:
        """Enable the Assembly tab only when consensus.fa exists on disk."""
        has_assembly = (
            job is not None
            and (job.assembly_dir / "consensus.fa").exists()
        )
        self._tab_bar.setTabEnabled(self.TAB_ASSEMBLY, has_assembly)
