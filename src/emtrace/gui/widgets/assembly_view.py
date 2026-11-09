"""Assembly tab widget — coverage strip, variant table, segment summary.

Concepts implemented
--------------------
A  Coverage strip: matplotlib canvas showing depth fill, low-coverage BED
   overlay (red shade), min-coverage threshold line, and VCF variant ticks
   in a separate sub-track.
C  Variant table: sortable QTableView of VCF records; clicking a row scrolls
   the coverage strip to that locus.
D  Per-segment summary: compact stats table (length, mean depth, low-cov %,
   variant count, completeness) with row click to switch segments.
"""

from __future__ import annotations

import gzip
import logging
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from PyQt6.QtCore import (
    QAbstractTableModel,
    QModelIndex,
    QSortFilterProxyModel,
    QThread,
    Qt,
    pyqtSignal,
)
from PyQt6.QtWidgets import (
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QMessageBox,
    QPushButton,
    QSizePolicy,
    QSplitter,
    QTableView,
    QVBoxLayout,
    QWidget,
)

from emtrace.gui.models.pandas_model import PandasTableModel

logger = logging.getLogger(__name__)

# ── Optional dependencies ────────────────────────────────────────────────────

try:
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as _FigureCanvas
    from matplotlib.figure import Figure as _Figure
    import matplotlib.patches as _mpatches
    import numpy as _np
    _HAS_MPL = True
except ImportError:
    _HAS_MPL = False

# ── Palette ──────────────────────────────────────────────────────────────────

_COL_DEPTH  = "#4C9BE8"
_COL_LOWCOV = "#E84C4C"
_COL_MINCOV = "#888888"
_COL_SNP    = "#F59D1E"
_COL_DEL    = "#E84C4C"
_COL_INS    = "#4CAF50"
_COL_MNP    = "#9C27B0"

_TYPE_COLORS = {"SNP": _COL_SNP, "DEL": _COL_DEL, "INS": _COL_INS, "MNP": _COL_MNP}


# ── Parsing helpers ──────────────────────────────────────────────────────────

def _parse_fasta_lengths(fa_path: Path) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    name: Optional[str] = None
    length = 0
    opener = gzip.open if str(fa_path).endswith(".gz") else open
    with opener(fa_path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = length
                name = line[1:].split()[0]
                length = 0
            else:
                length += len(line.strip())
    if name is not None:
        lengths[name] = length
    return lengths


def _parse_bed(bed_path: Path) -> Dict[str, List[Tuple[int, int]]]:
    """Return {chrom: [(start0, end0), ...]} (0-based half-open BED coords)."""
    intervals: Dict[str, List[Tuple[int, int]]] = {}
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            intervals.setdefault(chrom, []).append((start, end))
    return intervals


def _parse_vcf(vcf_path: Path) -> pd.DataFrame:
    """Parse VCF (plain or gzipped) into a DataFrame of variant records."""
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    rows: List[dict] = []
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom = parts[0]
            pos   = int(parts[1])
            ref   = parts[3]
            alt   = parts[4].split(",")[0]   # first ALT allele
            qual  = parts[5] if len(parts) > 5 and parts[5] != "." else ""
            if len(ref) == 1 and len(alt) == 1:
                vtype = "SNP"
            elif len(ref) > len(alt):
                vtype = "DEL"
            elif len(ref) < len(alt):
                vtype = "INS"
            else:
                vtype = "MNP"
            rows.append({
                "CHROM": chrom, "POS": pos,
                "REF": ref, "ALT": alt,
                "QUAL": qual, "TYPE": vtype,
            })
    if rows:
        return pd.DataFrame(rows)
    return pd.DataFrame(columns=["CHROM", "POS", "REF", "ALT", "QUAL", "TYPE"])


# ── Background depth loader ──────────────────────────────────────────────────

class _DepthLoader(QThread):
    """Run ``samtools depth -a -r <chrom> <bam>`` in a background thread."""

    finished = pyqtSignal(str, object, object)  # chrom, pos_array, depth_array
    error    = pyqtSignal(str, str)             # chrom, message

    def __init__(self, bam: Path, chrom: str, parent=None):
        super().__init__(parent)
        self._bam   = bam
        self._chrom = chrom

    def run(self) -> None:
        if not _HAS_MPL:
            self.error.emit(self._chrom, "numpy/matplotlib not available")
            return
        samtools = shutil.which("samtools")
        if not samtools:
            self.error.emit(self._chrom, "samtools not found on PATH")
            return
        try:
            out = subprocess.check_output(
                [samtools, "depth", "-a", "-r", self._chrom, str(self._bam)],
                stderr=subprocess.DEVNULL,
                timeout=300,
            )
            positions: List[int] = []
            depths: List[int] = []
            for line in out.decode().splitlines():
                parts = line.split("\t")
                if len(parts) >= 3:
                    positions.append(int(parts[1]))
                    depths.append(int(parts[2]))
            self.finished.emit(
                self._chrom,
                _np.array(positions, dtype=_np.int32),
                _np.array(depths,    dtype=_np.int32),
            )
        except Exception as exc:
            self.error.emit(self._chrom, str(exc))


# ── Segment summary model (Concept D) ───────────────────────────────────────

_SEG_HEADERS = ["Segment", "Length", "Mean depth", "Low-cov %", "Variants", "Completeness"]


class _SegmentModel(QAbstractTableModel):
    """Per-reference segment statistics (Concept D)."""

    def __init__(self, rows: List[dict], parent=None):
        super().__init__(parent)
        self._rows = rows

    def rowCount(self, parent: QModelIndex = QModelIndex()) -> int:
        return 0 if parent.isValid() else len(self._rows)

    def columnCount(self, parent: QModelIndex = QModelIndex()) -> int:
        return 0 if parent.isValid() else len(_SEG_HEADERS)

    def headerData(self, section: int, orientation, role=Qt.ItemDataRole.DisplayRole):
        if (
            orientation == Qt.Orientation.Horizontal
            and role == Qt.ItemDataRole.DisplayRole
        ):
            return _SEG_HEADERS[section]
        return None

    def data(self, index: QModelIndex, role=Qt.ItemDataRole.DisplayRole):
        if not index.isValid() or role != Qt.ItemDataRole.DisplayRole:
            return None
        row = self._rows[index.row()]
        col = index.column()
        if col == 0:
            return row.get("segment", "")
        if col == 1:
            return f"{row.get('length', 0):,}"
        if col == 2:
            d = row.get("mean_depth", 0.0)
            return f"{d:.1f}×" if d else "—"
        if col == 3:
            return f"{row.get('lowcov_pct', 0.0):.1f}%"
        if col == 4:
            return str(row.get("variants", 0))
        if col == 5:
            c = row.get("completeness", 0.0)
            return f"{c:.1f}%"
        return None

    def segment_name(self, row: int) -> str:
        return self._rows[row]["segment"] if 0 <= row < len(self._rows) else ""

    def update_mean_depth(self, chrom: str, mean_depth: float) -> None:
        for i, row in enumerate(self._rows):
            if row["segment"] == chrom:
                row["mean_depth"] = mean_depth
                tl = self.index(i, 2)
                self.dataChanged.emit(tl, tl)
                break


# ── Coverage canvas (Concept A) ─────────────────────────────────────────────

class _CoverageCanvas(QWidget):
    """Matplotlib coverage strip embedded in a QWidget.

    Two axes share the x-axis:
      - ax_depth: filled depth area, low-cov BED overlay, min-cov threshold.
      - ax_vars : vertical tick per variant coloured by type.
    """

    locus_clicked = pyqtSignal(int)  # emits 1-based genomic position

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        self._ax_depth = None
        self._ax_vars  = None
        self._chrom    = ""
        self._pos      = None
        self._dep      = None
        self._bed:     Dict[str, List[Tuple[int, int]]] = {}
        self._vcf_df:  Optional[pd.DataFrame] = None
        self._min_cov  = 10

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        if _HAS_MPL:
            self._fig    = _Figure(figsize=(8, 3))
            self._canvas = _FigureCanvas(self._fig)
            self._canvas.setSizePolicy(
                QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
            )
            layout.addWidget(self._canvas)
            self._canvas.mpl_connect("button_press_event", self._on_mpl_click)
            self._draw_placeholder("Select a segment to view coverage.")
        else:
            lbl = QLabel("matplotlib / numpy not available — install them to see coverage plots.")
            lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
            lbl.setStyleSheet("color:#999; font-size:12px;")
            layout.addWidget(lbl)

    # ------------------------------------------------------------------ #
    # Public API
    # ------------------------------------------------------------------ #

    def set_data(
        self,
        chrom: str,
        positions,            # numpy int32 array | None
        depths,               # numpy int32 array | None
        bed: Dict[str, List[Tuple[int, int]]],
        vcf_df: Optional[pd.DataFrame],
        min_cov: int = 10,
    ) -> None:
        self._chrom   = chrom
        self._pos     = positions
        self._dep     = depths
        self._bed     = bed
        self._vcf_df  = vcf_df
        self._min_cov = min_cov
        self._render()

    def show_loading(self, chrom: str) -> None:
        self._draw_placeholder(f"Loading depth for {chrom}…")

    def scroll_to(self, pos: int, window: int = 2000) -> None:
        """Zoom x-axis to ±window bp around *pos* (1-based)."""
        if not _HAS_MPL or self._ax_depth is None:
            return
        half = window // 2
        xlo = max(1, pos - half)
        xhi = pos + half
        self._ax_depth.set_xlim(xlo, xhi)
        self._ax_vars.set_xlim(xlo, xhi)
        self._canvas.draw_idle()

    def reset_zoom(self) -> None:
        if not _HAS_MPL or self._ax_depth is None or self._pos is None or not len(self._pos):
            return
        self._ax_depth.set_xlim(int(self._pos[0]), int(self._pos[-1]))
        self._ax_vars.set_xlim(int(self._pos[0]), int(self._pos[-1]))
        self._canvas.draw_idle()

    def save_figure(self, path: str) -> None:
        """Save the current coverage figure to *path* (PNG/PDF/SVG).

        Raises ``RuntimeError`` when matplotlib is unavailable.
        """
        if not _HAS_MPL:
            raise RuntimeError("matplotlib is not available")
        self._fig.savefig(path, bbox_inches="tight", dpi=150)

    # ------------------------------------------------------------------ #
    # Private
    # ------------------------------------------------------------------ #

    def _draw_placeholder(self, msg: str) -> None:
        if not _HAS_MPL:
            return
        self._ax_depth = None
        self._ax_vars  = None
        self._fig.clear()
        ax = self._fig.add_subplot(111)
        ax.text(0.5, 0.5, msg, ha="center", va="center",
                color="#aaa", fontsize=11, transform=ax.transAxes)
        ax.axis("off")
        self._canvas.draw_idle()

    def _render(self) -> None:
        if not _HAS_MPL:
            return
        self._fig.clear()

        gs = self._fig.add_gridspec(
            2, 1, height_ratios=[3, 1], hspace=0.05,
            left=0.07, right=0.97, top=0.93, bottom=0.12,
        )
        ax_d = self._fig.add_subplot(gs[0])
        ax_v = self._fig.add_subplot(gs[1], sharex=ax_d)
        self._ax_depth = ax_d
        self._ax_vars  = ax_v

        pos = self._pos
        dep = self._dep
        has_data = pos is not None and len(pos) > 0

        if has_data:
            # Depth fill
            ax_d.fill_between(pos, dep, color=_COL_DEPTH, alpha=0.65, linewidth=0)
            ax_d.plot(pos, dep, color=_COL_DEPTH, linewidth=0.5)

            # Low-coverage threshold line
            ax_d.axhline(
                self._min_cov, color=_COL_MINCOV, linewidth=0.9,
                linestyle="--", label=f"min cov = {self._min_cov}×",
            )

            # Low-cov BED overlay
            for start, end in self._bed.get(self._chrom, []):
                ax_d.axvspan(start, end, color=_COL_LOWCOV, alpha=0.18, linewidth=0)

            ax_d.set_xlim(int(pos[0]), int(pos[-1]))
            ax_d.set_ylim(bottom=0)
        else:
            ax_d.text(
                0.5, 0.5, "Depth data unavailable (BAM not found or samtools error)",
                ha="center", va="center", transform=ax_d.transAxes,
                color="#aaa", fontsize=10,
            )
            # Still show BED and VCF on the axes with arbitrary extent
            ref_len = 1
            for start, end in self._bed.get(self._chrom, []):
                ax_d.axvspan(start, end, color=_COL_LOWCOV, alpha=0.18, linewidth=0)
                ref_len = max(ref_len, end)
            ax_d.set_xlim(1, ref_len)

        ax_d.set_ylabel("Depth (×)", fontsize=8)
        ax_d.tick_params(labelsize=7)
        ax_d.legend(fontsize=7, loc="upper right", framealpha=0.6)
        ax_d.xaxis.set_visible(False)

        # ── Variant ticks ──────────────────────────────────────────────
        if self._vcf_df is not None and not self._vcf_df.empty:
            seg_vcf = self._vcf_df[self._vcf_df["CHROM"] == self._chrom]
            for _, vrow in seg_vcf.iterrows():
                col = _TYPE_COLORS.get(str(vrow.get("TYPE", "")), "#666")
                ax_v.vlines(int(vrow["POS"]), 0, 1, colors=col, linewidth=1.5, alpha=0.85)
            handles = [
                _mpatches.Patch(color=c, label=t)
                for t, c in _TYPE_COLORS.items()
            ]
            ax_v.legend(handles=handles, fontsize=6, loc="upper right",
                        ncol=4, framealpha=0.6)

        ax_v.set_yticks([])
        ax_v.set_ylim(0, 1.05)
        ax_v.set_xlabel(f"Position in {self._chrom} (bp)", fontsize=8)
        ax_v.tick_params(labelsize=7)
        for spine in ("left", "top", "right"):
            ax_v.spines[spine].set_visible(False)

        self._canvas.draw_idle()

    def _on_mpl_click(self, event) -> None:
        if event.inaxes in (self._ax_depth, self._ax_vars) and event.xdata is not None:
            self.locus_clicked.emit(int(round(event.xdata)))


# ── Main assembly view widget ────────────────────────────────────────────────

class AssemblyViewWidget(QWidget):
    """Assembly tab combining Concepts A (coverage strip), C (variant table),
    and D (per-segment summary)."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._job:          object = None
        self._asm_dir:      Optional[Path] = None
        self._bam:          Optional[Path] = None
        self._bed_data:     Dict[str, List[Tuple[int, int]]] = {}
        self._vcf_df:       Optional[pd.DataFrame] = None
        self._ref_lens:     Dict[str, int] = {}
        self._depth_cache:  Dict[str, Tuple] = {}
        self._loader:       Optional[_DepthLoader] = None
        self._min_cov:      int = 10

        self._build_ui()

    # ------------------------------------------------------------------ #
    # UI construction
    # ------------------------------------------------------------------ #

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(6, 6, 6, 6)
        root.setSpacing(4)

        # ── Toolbar ───────────────────────────────────────────────────
        tb = QHBoxLayout()
        tb.setSpacing(6)

        self._dir_label = QLabel("No assembly loaded")
        self._dir_label.setStyleSheet("color: #777; font-size: 11px;")
        self._dir_label.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        tb.addWidget(self._dir_label, stretch=1)

        tb.addWidget(QLabel("Segment:"))
        self._seg_combo = QComboBox()
        self._seg_combo.setMinimumWidth(220)
        self._seg_combo.currentTextChanged.connect(self._on_segment_changed)
        tb.addWidget(self._seg_combo)

        self._refresh_btn = QPushButton("Refresh")
        self._refresh_btn.setFixedWidth(90)
        self._refresh_btn.clicked.connect(self._load_assembly)
        tb.addWidget(self._refresh_btn)

        self._save_fig_btn = QPushButton("Save figure…")
        self._save_fig_btn.setFixedWidth(100)
        self._save_fig_btn.setEnabled(False)
        self._save_fig_btn.clicked.connect(self._on_save_figure)
        tb.addWidget(self._save_fig_btn)

        root.addLayout(tb)

        # ── Empty-state label ─────────────────────────────────────────
        self._empty_label = QLabel(
            "No assembly output found.\n"
            "Run the Assembly step to generate consensus and coverage data."
        )
        self._empty_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._empty_label.setStyleSheet("color: #999; font-size: 13px;")
        root.addWidget(self._empty_label, stretch=1)

        # ── Main splitter ─────────────────────────────────────────────
        self._splitter = QSplitter(Qt.Orientation.Vertical)
        self._splitter.setChildrenCollapsible(False)

        # D — segment summary table (fixed height)
        seg_container = QWidget()
        seg_layout = QVBoxLayout(seg_container)
        seg_layout.setContentsMargins(0, 0, 0, 0)
        seg_layout.setSpacing(2)
        seg_layout.addWidget(QLabel("Segments"))

        self._seg_table = QTableView()
        self._seg_table.setSelectionBehavior(QTableView.SelectionBehavior.SelectRows)
        self._seg_table.setSelectionMode(QTableView.SelectionMode.SingleSelection)
        self._seg_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.ResizeToContents
        )
        self._seg_table.verticalHeader().setVisible(False)
        self._seg_table.setAlternatingRowColors(True)
        self._seg_table.setMaximumHeight(150)
        self._seg_table.clicked.connect(self._on_seg_row_clicked)
        seg_layout.addWidget(self._seg_table)
        self._splitter.addWidget(seg_container)

        # A — coverage canvas
        self._cov_canvas = _CoverageCanvas()
        self._cov_canvas.setMinimumHeight(200)
        self._cov_canvas.locus_clicked.connect(self._on_locus_clicked)
        self._splitter.addWidget(self._cov_canvas)

        # C — variant table
        var_container = QWidget()
        var_layout = QVBoxLayout(var_container)
        var_layout.setContentsMargins(0, 0, 0, 0)
        var_layout.setSpacing(2)

        var_hdr = QHBoxLayout()
        self._var_hdr_label = QLabel("Variants")
        var_hdr.addWidget(self._var_hdr_label)
        self._var_count_label = QLabel("")
        self._var_count_label.setStyleSheet("color: #777; font-size: 11px;")
        var_hdr.addWidget(self._var_count_label, stretch=1)
        reset_btn = QPushButton("Reset zoom")
        reset_btn.setFixedWidth(90)
        reset_btn.clicked.connect(self._cov_canvas.reset_zoom)
        var_hdr.addWidget(reset_btn)
        var_layout.addLayout(var_hdr)

        self._var_table = QTableView()
        self._var_table.setSortingEnabled(True)
        self._var_table.setAlternatingRowColors(True)
        self._var_table.verticalHeader().setVisible(False)
        self._var_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.ResizeToContents
        )
        self._var_table.setSelectionBehavior(QTableView.SelectionBehavior.SelectRows)

        self._var_model = PandasTableModel()
        self._var_proxy = QSortFilterProxyModel()
        self._var_proxy.setSourceModel(self._var_model)
        self._var_table.setModel(self._var_proxy)
        self._var_table.clicked.connect(self._on_variant_clicked)
        var_layout.addWidget(self._var_table, stretch=1)
        self._splitter.addWidget(var_container)

        self._splitter.setSizes([130, 280, 180])
        root.addWidget(self._splitter, stretch=1)

        self._set_content_visible(False)

    def _set_content_visible(self, visible: bool) -> None:
        self._splitter.setVisible(visible)
        self._empty_label.setVisible(not visible)
        self._save_fig_btn.setEnabled(visible and _HAS_MPL)

    # ------------------------------------------------------------------ #
    # Public API
    # ------------------------------------------------------------------ #

    def set_job(self, job) -> None:
        self._job = job
        if job is None:
            self._clear()
        else:
            self._load_assembly()

    # ------------------------------------------------------------------ #
    # Loading
    # ------------------------------------------------------------------ #

    def _clear(self) -> None:
        self._asm_dir   = None
        self._bam       = None
        self._bed_data  = {}
        self._vcf_df    = None
        self._ref_lens  = {}
        self._depth_cache.clear()
        self._dir_label.setText("No assembly loaded")
        self._seg_combo.blockSignals(True)
        self._seg_combo.clear()
        self._seg_combo.blockSignals(False)
        self._var_model.set_dataframe(pd.DataFrame())
        self._var_count_label.setText("")
        self._set_content_visible(False)
        self._cov_canvas._draw_placeholder("No assembly loaded.")

    def _load_assembly(self) -> None:
        if self._job is not None:
            self._asm_dir = self._job.assembly_dir

            bam_from_job = self._job.get_file("assembly_bam")
            if bam_from_job and bam_from_job.exists():
                self._bam = bam_from_job

        if self._asm_dir is None or not self._asm_dir.exists():
            self._clear()
            return

        consensus = self._asm_dir / "consensus.fa"
        if not consensus.exists():
            self._set_content_visible(False)
            return

        # Locate BAM if not already resolved from job registry
        if self._bam is None or not self._bam.exists():
            candidate = self._asm_dir / "reads.refgen.md.bam"
            if candidate.exists():
                self._bam = candidate

        # Reference lengths
        ref_fa = self._asm_dir / "combined_ref.fa"
        try:
            self._ref_lens = _parse_fasta_lengths(
                ref_fa if ref_fa.exists() else consensus
            )
        except Exception:
            logger.exception("Could not parse reference FASTA lengths")
            self._ref_lens = {}

        # Low-coverage BED
        lowcov = self._asm_dir / "reads.refgen.lowcovmask.bed"
        try:
            self._bed_data = _parse_bed(lowcov) if lowcov.exists() else {}
        except Exception:
            logger.exception("Could not parse BED file")
            self._bed_data = {}

        # VCF
        vcf_path = self._asm_dir / "vcf.gz"
        try:
            self._vcf_df = _parse_vcf(vcf_path) if vcf_path.exists() else None
        except Exception:
            logger.exception("Could not parse VCF")
            self._vcf_df = None

        # Invalidate old depth cache (new assembly data on disk)
        self._depth_cache.clear()

        self._dir_label.setText(str(self._asm_dir))

        # Populate segment combo (suppress signal during rebuild)
        current = self._seg_combo.currentText()
        self._seg_combo.blockSignals(True)
        self._seg_combo.clear()
        for name in self._ref_lens:
            self._seg_combo.addItem(name)
        self._seg_combo.blockSignals(False)

        self._rebuild_segment_summary()

        if not self._ref_lens:
            self._set_content_visible(False)
            return

        self._set_content_visible(True)

        # Restore previous segment or pick first
        target = current if current in self._ref_lens else next(iter(self._ref_lens))
        self._seg_combo.setCurrentText(target)
        self._load_segment(target)

    def _rebuild_segment_summary(self) -> None:
        rows: List[dict] = []
        for seg, length in self._ref_lens.items():
            bed_intervals = self._bed_data.get(seg, [])
            lowcov_bp     = sum(e - s for s, e in bed_intervals)
            lowcov_pct    = (lowcov_bp / length * 100) if length else 0.0
            n_variants    = 0
            if self._vcf_df is not None and not self._vcf_df.empty:
                n_variants = int((self._vcf_df["CHROM"] == seg).sum())
            rows.append({
                "segment":      seg,
                "length":       length,
                "mean_depth":   0.0,
                "lowcov_pct":   lowcov_pct,
                "variants":     n_variants,
                "completeness": max(0.0, 100.0 - lowcov_pct),
            })
        model = _SegmentModel(rows)
        self._seg_table.setModel(model)
        self._seg_table.resizeColumnsToContents()
        self._seg_table.horizontalHeader().setSectionResizeMode(
            0, QHeaderView.ResizeMode.Stretch
        )

    def _load_segment(self, chrom: str) -> None:
        if not chrom:
            return

        self._update_variant_table(chrom)

        if chrom in self._depth_cache:
            pos, dep = self._depth_cache[chrom]
            self._render_coverage(chrom, pos, dep)
            return

        self._cov_canvas.show_loading(chrom)

        if self._bam and self._bam.exists():
            # Cancel any running loader for a previous segment
            if self._loader and self._loader.isRunning():
                self._loader.finished.disconnect()
                self._loader.error.disconnect()
                # Don't wait — let it finish silently; its result goes to cache
                self._loader = None

            self._loader = _DepthLoader(self._bam, chrom, parent=self)
            self._loader.finished.connect(self._on_depth_loaded)
            self._loader.error.connect(self._on_depth_error)
            self._loader.start()
        else:
            # No BAM available — render BED + VCF without depth fill
            self._render_coverage(chrom, None, None)

    def _update_variant_table(self, chrom: str) -> None:
        if self._vcf_df is not None and not self._vcf_df.empty:
            seg_df = self._vcf_df[self._vcf_df["CHROM"] == chrom].copy()
        else:
            seg_df = pd.DataFrame(
                columns=["CHROM", "POS", "REF", "ALT", "QUAL", "TYPE"]
            )
        self._var_model.set_dataframe(seg_df)
        n = len(seg_df)
        self._var_count_label.setText(
            f"{n} variant{'s' if n != 1 else ''}" if n else "No variants"
        )

    def _render_coverage(self, chrom: str, positions, depths) -> None:
        self._cov_canvas.set_data(
            chrom, positions, depths,
            self._bed_data, self._vcf_df, self._min_cov,
        )

        # Update mean depth in segment summary
        if _HAS_MPL and positions is not None and len(positions):
            try:
                mean_d = float(_np.mean(depths))
                model  = self._seg_table.model()
                if isinstance(model, _SegmentModel):
                    model.update_mean_depth(chrom, mean_d)
            except Exception:
                pass

    # ------------------------------------------------------------------ #
    # Signal handlers
    # ------------------------------------------------------------------ #

    def _on_segment_changed(self, chrom: str) -> None:
        self._load_segment(chrom)
        # Keep segment table selection in sync
        model = self._seg_table.model()
        if isinstance(model, _SegmentModel):
            for i, row in enumerate(model._rows):
                if row["segment"] == chrom:
                    self._seg_table.selectRow(i)
                    break

    def _on_seg_row_clicked(self, index: QModelIndex) -> None:
        model = self._seg_table.model()
        if isinstance(model, _SegmentModel):
            chrom = model.segment_name(index.row())
            if chrom and chrom != self._seg_combo.currentText():
                self._seg_combo.setCurrentText(chrom)

    def _on_variant_clicked(self, proxy_index: QModelIndex) -> None:
        src = self._var_proxy.mapToSource(proxy_index)
        df  = self._var_model.dataframe
        if df is None or df.empty:
            return
        # PandasTableModel has virtual '#' column at index 0; data cols start at 1
        df_col = src.column() - 1
        if df_col < 0:
            return
        try:
            pos = int(df.iloc[src.row()]["POS"])
        except (KeyError, ValueError, IndexError):
            return
        self._cov_canvas.scroll_to(pos)

    def _on_locus_clicked(self, pos: int) -> None:
        logger.debug("Coverage strip: position %d clicked", pos)

    def _on_depth_loaded(self, chrom: str, positions, depths) -> None:
        self._depth_cache[chrom] = (positions, depths)
        if self._seg_combo.currentText() == chrom:
            self._render_coverage(chrom, positions, depths)

    def _on_depth_error(self, chrom: str, msg: str) -> None:
        logger.warning("Depth load failed for %s: %s", chrom, msg)
        if self._seg_combo.currentText() == chrom:
            self._render_coverage(chrom, None, None)

    def _on_save_figure(self) -> None:
        """Open a save dialog and export the coverage figure."""
        chrom = self._seg_combo.currentText() or "coverage"
        default_name = f"coverage_{chrom}.pdf"
        if self._asm_dir is not None:
            default_path = str(self._asm_dir / default_name)
        else:
            default_path = default_name

        path, selected_filter = QFileDialog.getSaveFileName(
            self,
            "Save Coverage Figure",
            default_path,
            "PDF (*.pdf);;PNG (*.png);;SVG (*.svg)",
        )
        if not path:
            return

        # Append extension when the user didn't type one
        ext_map = {"PDF (*.pdf)": ".pdf", "PNG (*.png)": ".png", "SVG (*.svg)": ".svg"}
        ext = ext_map.get(selected_filter, ".pdf")
        if not any(path.lower().endswith(e) for e in (".pdf", ".png", ".svg")):
            path += ext

        try:
            self._cov_canvas.save_figure(path)
            logger.info("Coverage figure saved to %s", path)
        except Exception as exc:
            QMessageBox.critical(self, "Save failed", f"Could not save figure:\n{exc}")
