"""QC Plot widget — read length and quality histograms + PDF report.

Reads a FASTQ file (plain or gzipped) and renders:
  - (a) Read length distribution histogram
  - (b) Per-read mean quality score distribution histogram
  - (c) 5′ head per-position quality + coverage
  - (d) 3′ tail per-position quality + coverage

PDF report: multi-page A4, two figures per page, scientific captions.
Falls back gracefully if matplotlib is not available.
"""

from __future__ import annotations

import gzip
import logging
import statistics
import textwrap
from datetime import datetime
from pathlib import Path
from typing import Optional

from PyQt6.QtCore import QEvent, QObject, QThread, Qt, pyqtSignal
from PyQt6.QtWidgets import (
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

logger = logging.getLogger(__name__)

try:
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.figure import Figure
    from matplotlib.ticker import FuncFormatter as _FuncFormatter

    _HAS_MPL = True
    # Y-axis formatter: raw counts → thousands (e.g. 25000 → "25")
    _KFMT = _FuncFormatter(lambda x, _: f"{x/1000:g}")
except ImportError:  # pragma: no cover
    _HAS_MPL = False
    _KFMT = None

try:
    import numpy as _np
    _HAS_NP = True
except ImportError:
    _HAS_NP = False


def _trim_path(path_str: str) -> str:
    """Shorten a long path to …/parent/parent/filename (last 3 components)."""
    if not path_str or path_str == "—":
        return path_str
    parts = Path(path_str).parts
    if len(parts) <= 3:
        return path_str
    return "…/" + "/".join(parts[-3:])


# Number of positions captured from the 5′ head and 3′ tail of each read
_POSITION_LEN = 500

# A4 dimensions in inches
_A4_W, _A4_H = 8.27, 11.69

# (path, mtime) → scan-data dict
_ScanCache = dict[tuple[Path, float], dict]

# ── FASTQ stage keys (displayed in the file combo) ──────────────────

_FASTQ_STAGES = [
    ("Input reads",             "input_reads"),
    ("Current reads",           "current_reads"),
    ("Porechop-filtered reads", "porechop_reads"),
    ("Chopper-filtered reads",  "chopped_reads"),
    ("Host-filtered reads",     "host_filtered_reads"),
    ("Annotated vi reads",      "ann1_vi_reads"),
    ("Annotated ba reads",      "ann1_ab_reads"),
    ("Annotated ph reads",      "ann1_ph_reads"),
    ("Unannotated reads",       "ann1_un_reads"),
]

# "Current reads" is a convenience alias — excluded from the PDF report
# to avoid duplicating whichever stage it points to.
_REPORT_STAGES = [s for s in _FASTQ_STAGES if s[1] != "current_reads"]


# ── Shared scan kernel ────────────────────────────────────────────────


def _scan_fastq_data(path: Path, max_reads: int = 250_000) -> dict:
    """Scan *path* synchronously; call from a worker thread."""
    P = _POSITION_LEN
    lengths: list[int] = []
    qualities: list[float] = []

    if _HAS_NP:
        head_sum = _np.zeros(P, dtype=_np.float64)
        head_cnt = _np.zeros(P, dtype=_np.int64)
        tail_sum = _np.zeros(P, dtype=_np.float64)
        tail_cnt = _np.zeros(P, dtype=_np.int64)
    else:
        import array as _array
        head_sum = _array.array("d", [0.0] * P)
        head_cnt = _array.array("l", [0] * P)
        tail_sum = _array.array("d", [0.0] * P)
        tail_cnt = _array.array("l", [0] * P)

    try:
        fh = (
            gzip.open(path, "rt")
            if str(path).endswith(".gz")
            else open(path, "r")
        )
        with fh:
            n = 0
            while n < max_reads:
                header = fh.readline()
                if not header:
                    break
                seq = fh.readline().rstrip("\n")
                _plus = fh.readline()
                qual = fh.readline().rstrip("\n")
                if not (header and seq and _plus and qual):
                    break
                rlen = len(seq)
                lengths.append(rlen)
                if qual:
                    if _HAS_NP:
                        q = _np.frombuffer(qual.encode(), dtype=_np.uint8).astype(_np.float32) - 33
                        qualities.append(float(q.mean()))
                        head_end = min(rlen, P)
                        head_sum[:head_end] += q[:head_end]
                        head_cnt[:head_end] += 1
                        tail_start = max(0, rlen - P)
                        tail_len = rlen - tail_start
                        tail_sum[:tail_len] += q[tail_start:]
                        tail_cnt[:tail_len] += 1
                    else:
                        q_scores = [ord(c) - 33 for c in qual]
                        qualities.append(sum(q_scores) / len(q_scores))
                        head_end = min(rlen, P)
                        for i in range(head_end):
                            head_sum[i] += q_scores[i]
                            head_cnt[i] += 1
                        tail_start = max(0, rlen - P)
                        for i, pos in enumerate(range(tail_start, rlen)):
                            tail_sum[i] += q_scores[pos]
                            tail_cnt[i] += 1
                n += 1
    except Exception as exc:
        logger.warning("FASTQ scan failed (%s): %s", path.name, exc)

    return {
        "lengths": lengths,
        "qualities": qualities,
        "head_sum": list(head_sum),
        "head_cnt": list(head_cnt),
        "tail_sum": list(tail_sum),
        "tail_cnt": list(tail_cnt),
    }


# ── Shared plotting helpers ───────────────────────────────────────────


def _draw_positional(ax, pos_sum, pos_cnt, title: str, color: str) -> None:
    """Mean per-position quality on *ax* with a coverage fill on a twin axis."""
    from matplotlib.patches import Patch

    last = next((i for i in range(len(pos_cnt) - 1, -1, -1) if pos_cnt[i] > 0), -1)
    ax.set_title(title, fontsize=8)
    ax.set_xlabel("Position (bp)", fontsize=7)
    ax.set_ylabel("Mean Phred quality", color=color, fontsize=7)
    ax.tick_params(axis="both", labelsize=7)
    if last < 0:
        return

    xs = list(range(last + 1))
    mean_q = [pos_sum[i] / pos_cnt[i] if pos_cnt[i] else 0.0 for i in xs]
    coverage = [pos_cnt[i] for i in xs]

    ax.plot(xs, mean_q, color=color, linewidth=1.0, label="Mean Q")
    ax.tick_params(axis="y", labelcolor=color)
    ax.set_ylim(bottom=0)

    ax2 = ax.twinx()
    ax2.fill_between(xs, coverage, alpha=0.35, color="#aaaaaa")
    ax2.set_ylabel("Reads (thousands)", color="#777777", fontsize=7)
    ax2.tick_params(axis="y", labelcolor="#777777", labelsize=7)
    ax2.set_ylim(bottom=0)
    if _KFMT is not None:
        ax2.yaxis.set_major_formatter(_KFMT)

    handles = ax.get_lines() + [Patch(facecolor="#aaaaaa", alpha=0.5, label="Coverage")]
    ax.legend(handles=handles, fontsize=6)


def _fill_figure(fig: Figure, data: dict, stage_label: str = "", path_str: str = "") -> dict:
    """Populate *fig* with the 2×2 QC grid.  Returns a stats dict for captions.

    Panel labels (a)–(d) are drawn inside each subplot.
    """
    lengths = data.get("lengths", [])

    if not lengths:
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, "No reads in this file",
                ha="center", va="center", fontsize=12, color="#999999",
                transform=ax.transAxes)
        ax.axis("off")
        return {"empty": True, "n_reads": 0}

    qualities = data.get("qualities", [])
    head_sum  = data.get("head_sum", [])
    head_cnt  = data.get("head_cnt", [])
    tail_sum  = data.get("tail_sum", [])
    tail_cnt  = data.get("tail_cnt", [])

    median_len = statistics.median(lengths)
    median_q   = statistics.median(qualities) if qualities else None

    _PLK = dict(fontsize=9, fontweight="bold", va="top", ha="left")

    # hspace/wspace only — never pass left/right/top/bottom here.
    # When fig is a SubFigure, add_gridspec propagates any explicit margin
    # params to the parent Figure via subplots_adjust, overriding page margins.
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.28)

    # (a) Length distribution
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.hist(lengths, bins=80, color="#1f77b4", edgecolor="white", linewidth=0.3)
    ax1.set_xlabel("Read length (bp)", fontsize=7)
    ax1.set_ylabel("Count (thousands)", fontsize=7)
    ax1.set_title("Read length distribution", fontsize=8)
    ax1.tick_params(labelsize=7)
    if _KFMT is not None:
        ax1.yaxis.set_major_formatter(_KFMT)
    ax1.axvline(median_len, color="#d62728", linestyle="--", linewidth=1,
                label=f"Median: {median_len:,.0f} bp\nn = {len(lengths):,}")
    ax1.legend(fontsize=6)
    ax1.text(0.02, 0.97, "(a)", transform=ax1.transAxes, **_PLK)

    # (b) Quality distribution
    ax2 = fig.add_subplot(gs[0, 1])
    if qualities:
        ax2.hist(qualities, bins=60, color="#2ca02c", edgecolor="white", linewidth=0.3)
        ax2.axvline(median_q, color="#d62728", linestyle="--", linewidth=1,
                    label=f"Median: {median_q:.1f}")
        ax2.legend(fontsize=6)
    ax2.set_xlabel("Mean Phred quality", fontsize=7)
    ax2.set_ylabel("Count (thousands)", fontsize=7)
    ax2.set_title("Quality score distribution", fontsize=8)
    ax2.tick_params(labelsize=7)
    if _KFMT is not None:
        ax2.yaxis.set_major_formatter(_KFMT)
    ax2.text(0.02, 0.97, "(b)", transform=ax2.transAxes, **_PLK)

    # (c) Head quality + coverage
    ax3 = fig.add_subplot(gs[1, 0])
    _draw_positional(ax3, head_sum, head_cnt,
                     title=f"5′ head quality (first {_POSITION_LEN} bp)",
                     color="#1f77b4")
    ax3.text(0.02, 0.97, "(c)", transform=ax3.transAxes, **_PLK)

    # (d) Tail quality + coverage
    ax4 = fig.add_subplot(gs[1, 1])
    _draw_positional(ax4, tail_sum, tail_cnt,
                     title=f"3′ tail quality (last {_POSITION_LEN} bp)",
                     color="#9467bd")
    ax4.text(0.02, 0.97, "(d)", transform=ax4.transAxes, **_PLK)

    return {
        "empty": False,
        "n_reads": len(lengths),
        "median_len": median_len,
        "median_q": median_q,
    }


def _draw_caption(cap_ax, fig_num: int, stage_label: str, stats: dict) -> None:
    """Render a scientific figure caption inside *cap_ax*."""
    cap_ax.axis("off")
    empty = stats.get("empty", False)

    if empty:
        desc = "No reads found in this file."
    else:
        desc = (
            "(a) Read length distribution."
            " (b) Per-read mean Phred quality score distribution."
            " (c) Mean quality per position from the 5′ end."
            " (d) As (c) from the 3′ end. Grey fill shows read coverage."
        )

    full_text = f"Figure {fig_num}: {stage_label}. " + desc
    wrapped   = textwrap.fill(full_text, width=115)

    cap_ax.text(0.0, 0.98, wrapped,
                fontsize=8, va="top", ha="left", linespacing=1.3,
                transform=cap_ax.transAxes, clip_on=False)


# ── Background workers ────────────────────────────────────────────────


class _FastqScanWorker(QThread):
    """Scan a single FASTQ in a background thread."""

    finished = pyqtSignal(object)  # emits a dict

    def __init__(self, path: Path, max_reads: int = 250_000, parent=None):
        super().__init__(parent)
        self._path = path
        self._max_reads = max_reads

    def run(self) -> None:
        self.finished.emit(_scan_fastq_data(self._path, self._max_reads))


class _PdfReportWorker(QThread):
    """Scan all tracked FASTQ stages and write a multi-page A4 PDF report.

    Layout: cover page with metadata, then two figures per A4 page.
    Each figure = 2×2 subplot grid with a scientific caption below.
    """

    progress    = pyqtSignal(str)           # status message for the UI
    report_done = pyqtSignal(str, object)   # (pdf_path, {key: (path, data)})
    report_error = pyqtSignal(str)

    def __init__(
        self,
        stages: list[tuple[str, str]],
        job_files: dict[str, Path],
        cache_snapshot: dict,
        output_path: str,
        job_name: str = "",
        input_str: str = "",
        output_dir_str: str = "",
        max_reads: int = 250_000,
        parent=None,
    ):
        super().__init__(parent)
        self._stages         = stages
        self._job_files      = dict(job_files)
        self._cache          = dict(cache_snapshot)
        self._output_path    = output_path
        self._job_name       = job_name
        self._input_str      = input_str
        self._output_dir_str = output_dir_str
        self._max_reads      = max_reads

    # ------------------------------------------------------------------

    def run(self) -> None:
        from matplotlib.backends.backend_pdf import PdfPages

        # ── 1. Collect + scan all available stages ─────────────────
        scanned: list[tuple[str, str, dict]] = []  # (label, path_str, data)
        new_scans: dict[str, tuple[Path, dict]] = {}

        for label, key in self._stages:
            raw = self._job_files.get(key)
            if raw is None:
                continue
            path = Path(raw)

            self.progress.emit(f"Scanning {label}…")

            data: dict | None = None
            try:
                mtime = path.stat().st_mtime
                data = self._cache.get((path, mtime))
            except OSError:
                pass

            if data is None:
                try:
                    is_empty = not path.exists() or path.stat().st_size == 0
                except OSError:
                    is_empty = True
                if is_empty:
                    data = {}
                else:
                    data = _scan_fastq_data(path, self._max_reads)
                    if data.get("lengths"):
                        new_scans[key] = (path, data)

            scanned.append((label, str(path), data))

        if not scanned:
            self.report_error.emit("No FASTQ files found for this job.")
            return

        # ── 2. Render PDF ──────────────────────────────────────────
        self.progress.emit("Rendering PDF…")
        try:
            with PdfPages(self._output_path) as pdf:
                # Cover page
                pdf.savefig(self._make_cover(len(scanned)))

                # Content pages: 2 figures per A4 page
                fig_num = 1
                for i in range(0, len(scanned), 2):
                    batch = scanned[i : i + 2]
                    pdf.savefig(self._make_content_page(batch, fig_num))
                    fig_num += len(batch)

            self.report_done.emit(self._output_path, new_scans)
        except Exception as exc:
            self.report_error.emit(str(exc))

    # ------------------------------------------------------------------

    def _make_cover(self, n_figures: int) -> Figure:
        fig = Figure(figsize=(_A4_W, _A4_H))
        fig.subplots_adjust(left=0.121, right=0.879, top=0.915, bottom=0.086)
        ax = fig.add_subplot(111)
        ax.axis("off")

        # Title
        ax.text(0.5, 0.82, "QC Report",
                ha="center", va="center", fontsize=32, fontweight="bold",
                transform=ax.transAxes)

        # Horizontal rule
        ax.plot([0.05, 0.95], [0.75, 0.75], color="#cccccc", linewidth=1,
                transform=ax.transAxes, clip_on=False)

        # Metadata block
        now = datetime.now().strftime("%Y-%m-%d  %H:%M")
        meta = [
            ("Job name",   self._job_name                   or "—"),
            ("Input data", _trim_path(self._input_str)      or "—"),
            ("Output dir", _trim_path(self._output_dir_str) or "—"),
            ("Date",       now),
            ("Figures",    str(n_figures)),
        ]
        label_x, value_x, y0, dy = 0.08, 0.32, 0.70, 0.038
        for i, (lbl, val) in enumerate(meta):
            y = y0 - i * dy
            ax.text(label_x, y, lbl + ":",
                    ha="left", va="top", fontsize=11, fontweight="bold",
                    transform=ax.transAxes)
            ax.text(value_x, y, val,
                    ha="left", va="top", fontsize=11,
                    transform=ax.transAxes)

        fig.tight_layout()
        return fig

    def _make_content_page(
        self,
        batch: list[tuple[str, str, dict]],
        fig_num_start: int,
    ) -> Figure:
        """Render 1 or 2 figures on a single A4 page using subfigures."""
        # Each figure slot: plot area + caption. Two slots per A4 page.
        # 1-inch margins: left/right = 1/8.27 ≈ 0.121, top/bottom = 1/11.69 ≈ 0.086
        ratios = [10, 1.2, 10, 1.2]
        page = Figure(figsize=(_A4_W, _A4_H))
        subfigs = page.subfigures(nrows=4, ncols=1,
                                  height_ratios=ratios, hspace=0.0)

        for slot, (label, path_str, data) in enumerate(batch):
            plot_sf = subfigs[slot * 2]
            cap_sf  = subfigs[slot * 2 + 1]

            stats = _fill_figure(plot_sf, data, label, path_str)

            # Caption axis spans full width of the subfigure
            cap_ax = cap_sf.add_subplot(111)
            _draw_caption(cap_ax, fig_num_start + slot, label, stats)

        # Leave unused slots blank (single-figure last page)
        for slot in range(len(batch), 2):
            subfigs[slot * 2].set_visible(False)
            subfigs[slot * 2 + 1].set_visible(False)

        # Set page margins LAST — _fill_figure's add_gridspec on the SubFigure
        # propagates subplot params to the parent figure, overriding any margins
        # set earlier.  Calling subplots_adjust here ensures it wins.
        page.subplots_adjust(left=0.121, right=0.879, top=0.915, bottom=0.086)

        return page


# ── Widget ────────────────────────────────────────────────────────────


class QCPlotWidget(QWidget):
    """Read length and quality distribution histograms with PDF export."""

    def __init__(self, parent=None):
        super().__init__(parent)

        self._worker: Optional[_FastqScanWorker] = None
        self._pdf_worker: Optional[_PdfReportWorker] = None
        self._job_files: dict[str, Path] = {}
        self._output_dir: Optional[Path] = None
        self._current_path: Optional[Path] = None
        self._scan_cache: _ScanCache = {}

        # ── Controls ─────────────────────────────────────────────────
        self._stage_combo = QComboBox()
        self._stage_combo.setMinimumWidth(180)
        for label, key in _FASTQ_STAGES:
            self._stage_combo.addItem(label, key)

        self._refresh_btn = QPushButton("Refresh")
        self._refresh_btn.setFixedWidth(80)

        self._save_btn = QPushButton("Save figure")
        self._save_btn.setFixedWidth(90)
        self._save_btn.setEnabled(False)

        self._report_btn = QPushButton("PDF report")
        self._report_btn.setFixedWidth(90)
        self._report_btn.setEnabled(False)

        self._status_label = QLabel()
        self._status_label.setStyleSheet("color: #888; font-size: 11px;")

        self._file_label = QLabel()
        self._file_label.setStyleSheet(
            "color: #555; font-size: 11px;"
            " font-family: 'Menlo', 'Consolas', 'Courier New', 'DejaVu Sans Mono';"
            " padding: 2px 4px; background: #f5f5f5; border-radius: 3px;"
        )
        self._file_label.setWordWrap(False)

        ctrl_row = QHBoxLayout()
        ctrl_row.setContentsMargins(4, 2, 4, 2)
        ctrl_row.addWidget(QLabel("FASTQ stage:"))
        ctrl_row.addWidget(self._stage_combo)
        ctrl_row.addWidget(self._refresh_btn)
        ctrl_row.addWidget(self._save_btn)
        ctrl_row.addWidget(self._report_btn)
        ctrl_row.addWidget(self._status_label, stretch=1)

        # ── Canvas / empty state ─────────────────────────────────────
        self._empty_label = QLabel(
            "No FASTQ data available.\n"
            "Load a sample and run preprocessing steps."
            if _HAS_MPL
            else "matplotlib is not installed.\n"
            "Install matplotlib to view QC plots."
        )
        self._empty_label_default = (
            "No FASTQ data available.\n"
            "Load a sample and run preprocessing steps."
        )
        self._empty_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._empty_label.setStyleSheet("color: #999; font-size: 14px;")

        self._canvas: Optional[FigureCanvas] = None
        self._fig: Optional[Figure] = None
        self._scroll: Optional[QScrollArea] = None
        if _HAS_MPL:
            self._fig = Figure(figsize=(10, 8), dpi=100, layout="constrained")
            self._canvas = FigureCanvas(self._fig)
            self._canvas.setMinimumHeight(750)
            self._scroll = QScrollArea()
            self._scroll.setWidgetResizable(True)
            self._scroll.setWidget(self._canvas)
            # Forward wheel events from the canvas to the scroll area
            self._canvas.installEventFilter(self)

        # ── Layout ───────────────────────────────────────────────────
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addLayout(ctrl_row)
        layout.addWidget(self._file_label)
        if self._scroll is not None:
            layout.addWidget(self._scroll, stretch=1)
        layout.addWidget(self._empty_label, stretch=1)

        self._set_has_data(False)

        # ── Connections ──────────────────────────────────────────────
        self._stage_combo.currentIndexChanged.connect(self._on_stage_changed)
        self._refresh_btn.clicked.connect(self._on_refresh)
        self._save_btn.clicked.connect(self._on_save_figure)
        self._report_btn.clicked.connect(self._on_export_pdf_report)

    # ------------------------------------------------------------------
    # Event filter — forward canvas wheel events to the scroll area
    # ------------------------------------------------------------------

    def eventFilter(self, obj: QObject, event: QEvent) -> bool:
        if obj is self._canvas and event.type() == QEvent.Type.Wheel:
            if self._scroll is not None:
                self._scroll.wheelEvent(event)
            return True
        return super().eventFilter(obj, event)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def set_job_files(
        self,
        files: dict[str, Path],
        output_dir: Optional[Path] = None,
    ) -> None:
        """Update the map of available FASTQ files and rebuild the combo."""
        self._job_files = dict(files)
        self._output_dir = output_dir

        self._stage_combo.blockSignals(True)
        prev_key = self._stage_combo.currentData()
        self._stage_combo.clear()
        for label, key in _FASTQ_STAGES:
            p = files.get(key)
            if p and Path(p).exists():
                self._stage_combo.addItem(label, key)

        restored = False
        if prev_key:
            for i in range(self._stage_combo.count()):
                if self._stage_combo.itemData(i) == prev_key:
                    self._stage_combo.setCurrentIndex(i)
                    restored = True
                    break
        if not restored and self._stage_combo.count() > 0:
            self._stage_combo.setCurrentIndex(0)
        self._stage_combo.blockSignals(False)

        self._report_btn.setEnabled(self._stage_combo.count() > 0)
        self._on_stage_changed()

    def load_fastq(self, path: Path) -> None:
        """Directly load and plot a specific FASTQ file."""
        self._stage_combo.setCurrentIndex(-1)
        self._scan_fastq(path)

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _set_has_data(self, has: bool) -> None:
        if self._scroll is not None:
            self._scroll.setVisible(has)
        self._empty_label.setVisible(not has)
        self._save_btn.setEnabled(has)

    def _get_cached(self, path: Path) -> dict | None:
        try:
            mtime = path.stat().st_mtime
        except OSError:
            return None
        return self._scan_cache.get((path, mtime))

    def _store_cache(self, path: Path, data: dict) -> None:
        try:
            mtime = path.stat().st_mtime
            self._scan_cache[(path, mtime)] = data
        except OSError:
            pass

    def _format_path(self, path: Path) -> str:
        if self._output_dir is not None:
            try:
                return str(path.relative_to(self._output_dir))
            except ValueError:
                pass
        return str(path)

    def _set_file_label(self, path: Optional[Path]) -> None:
        if path is None:
            self._file_label.setText("")
            self._file_label.setToolTip("")
        else:
            self._file_label.setText(self._format_path(path))
            self._file_label.setToolTip(str(path))

    # ------------------------------------------------------------------
    # Slots — interactive viewer
    # ------------------------------------------------------------------

    def _on_refresh(self) -> None:
        if self._current_path:
            stale = [k for k in self._scan_cache if k[0] == self._current_path]
            for k in stale:
                del self._scan_cache[k]
        self._on_stage_changed()

    def _on_save_figure(self) -> None:
        if self._fig is None:
            return
        stem = self._current_path.stem if self._current_path else "qc_plot"
        default_name = f"{stem}.qc.png"
        if self._output_dir is not None:
            figures_dir = self._output_dir / "results" / "plots" / "qc"
            figures_dir.mkdir(parents=True, exist_ok=True)
            default_name = str(figures_dir / default_name)
        save_path, _ = QFileDialog.getSaveFileName(
            self, "Save Figure", default_name,
            "PNG image (*.png);;PDF document (*.pdf);;SVG image (*.svg)",
        )
        if not save_path:
            return
        try:
            self._fig.savefig(save_path, dpi=150, bbox_inches="tight")
            self._status_label.setText(f"Saved → {Path(save_path).name}")
        except Exception as exc:
            self._status_label.setText(f"Save failed: {exc}")

    def _on_stage_changed(self, _index: int = 0) -> None:
        key = self._stage_combo.currentData()
        if not key:
            return
        path = self._job_files.get(key)
        if path and Path(path).exists():
            self._scan_fastq(Path(path))
        else:
            self._set_has_data(False)
            self._set_file_label(None)
            self._status_label.setText(f"File not found for stage: {key}")

    def _scan_fastq(self, path: Path) -> None:
        if not _HAS_MPL:
            return
        if self._worker is not None and self._worker.isRunning():
            return

        self._current_path = path
        self._set_file_label(path)

        try:
            empty = path.stat().st_size == 0
        except OSError:
            empty = False

        if empty:
            self._empty_label.setText(f"FASTQ file is empty.\n{self._format_path(path)}")
            self._set_has_data(False)
            self._status_label.setText("File is empty.")
            return

        cached = self._get_cached(path)
        if cached is not None:
            self._status_label.setText(f"{path.name} (cached)")
            self._empty_label.setText(self._empty_label_default)
            self._on_scan_done(cached)
            return

        self._empty_label.setText(self._empty_label_default)
        self._status_label.setText(f"Scanning {path.name}…")
        self._worker = _FastqScanWorker(path, parent=self)
        self._worker.finished.connect(self._on_scan_done)
        self._worker.start()

    def _on_scan_done(self, data: dict) -> None:
        lengths = data.get("lengths", [])
        if not lengths:
            self._empty_label.setText(
                f"FASTQ file is empty.\n"
                f"{self._format_path(self._current_path) if self._current_path else ''}"
            )
            self._set_has_data(False)
            self._status_label.setText("No reads found.")
            return

        if self._current_path:
            self._store_cache(self._current_path, data)
        self._set_has_data(True)
        self._status_label.setText(f"{len(lengths):,} reads scanned.")
        self._plot(data)

    def _plot(self, data: dict) -> None:
        if self._fig is None or self._canvas is None:
            return
        self._fig.clear()
        _fill_figure(self._fig, data,
                     stage_label=self._stage_combo.currentText() or "",
                     path_str=str(self._current_path or ""))
        self._canvas.draw()

    # ------------------------------------------------------------------
    # Slots — PDF report
    # ------------------------------------------------------------------

    def _on_export_pdf_report(self) -> None:
        if self._pdf_worker is not None and self._pdf_worker.isRunning():
            self._status_label.setText("Report already in progress…")
            return

        sample_name = self._output_dir.name if self._output_dir else "qc_report"
        default_name = f"{sample_name}.qc_report.pdf"
        if self._output_dir is not None:
            reports_dir = self._output_dir / "results" / "reports"
            reports_dir.mkdir(parents=True, exist_ok=True)
            default_name = str(reports_dir / default_name)

        save_path, _ = QFileDialog.getSaveFileName(
            self, "Export PDF Report", default_name, "PDF document (*.pdf)",
        )
        if not save_path:
            return

        # Derive job metadata from available context
        job_name      = self._output_dir.name if self._output_dir else ""
        input_reads   = self._job_files.get("input_reads")
        input_str     = str(input_reads) if input_reads else ""
        output_dir_str = str(self._output_dir) if self._output_dir else ""

        self._report_btn.setEnabled(False)
        self._status_label.setText("Building PDF report…")

        self._pdf_worker = _PdfReportWorker(
            stages=_REPORT_STAGES,
            job_files=self._job_files,
            cache_snapshot=self._scan_cache,
            output_path=save_path,
            job_name=job_name,
            input_str=input_str,
            output_dir_str=output_dir_str,
            parent=self,
        )
        self._pdf_worker.progress.connect(self._status_label.setText)
        self._pdf_worker.report_done.connect(self._on_pdf_done)
        self._pdf_worker.report_error.connect(self._on_pdf_error)
        self._pdf_worker.start()

    def _on_pdf_done(self, pdf_path: str, new_scans: object) -> None:
        self._report_btn.setEnabled(True)
        self._status_label.setText(f"Report saved → {Path(pdf_path).name}")
        for path, data in (new_scans or {}).values():
            self._store_cache(path, data)

    def _on_pdf_error(self, msg: str) -> None:
        self._report_btn.setEnabled(True)
        self._status_label.setText(f"Report failed: {msg}")
        logger.error("PDF report failed: %s", msg)
