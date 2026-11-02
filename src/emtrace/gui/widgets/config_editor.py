"""Config editor widget — Phase 3.2.

Scrollable form with collapsible sections, auto-generated from EMtraceConfig.
Emits ``config_changed(EMtraceConfig)`` when the user clicks Apply.

Sections exposed:
  - Runtime          — threads, output dir, Entrez credentials
  - Preprocessing    — Chopper + Porechop_ABI (basic + Advanced sub-section)
  - Background Filtering   — host genome checkboxes + alignment thresholds
  - Annotation       — strategy selector, ann1/ann2, tool settings
  - Assembly         — minimap2 mapping, coverage/variant thresholds
  - EM Abundance     — EM parameters
  - Reporting        — thresholds, Excel export, read sorting
  - Output           — pack / FASTA / realtime table flags
"""

from __future__ import annotations

import logging
import re
import types as _types
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, get_args, get_origin

import pandas as pd
from pydantic import BaseModel
from emtrace.analysis.taxonomy_utils import is_silva_annotation
from pydantic.fields import FieldInfo
from PyQt6.QtCore import Qt, QTimer, pyqtSignal
from PyQt6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from emtrace.core.config import (
    AnnotationConfig,
    AssemblyConfig,
    BgFilterConfig,
    AssemblyMinimapConfig,
    BgFilterMinimapConfig,
    ChopperConfig,
    DiamondBlastpConfig,
    EMConfig,
    EMtraceConfig,
    MinimapConfig,
    BlastnConfig,
    OrfipyConfig,
    OutputConfig,
    PorechopABIConfig,
    PreprocessConfig,
    ReportConfig,
    RuntimeConfig,
    SubsampleConfig,
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Triangle collapse glyphs
# ---------------------------------------------------------------------------
_COLLAPSED = "\u25b6"   # ▶
_EXPANDED = "\u25bc"    # ▼

# ---------------------------------------------------------------------------
# Enumerated string field options (not Literal in Pydantic, but fixed sets)
# ---------------------------------------------------------------------------
_CHOPPER_TRIM_APPROACHES = [
    "trim-by-quality", "fixed-crop", "best-read-segment", "split-by-low-quality",
]
_MINIMAP_PRESETS = ["map-ont", "map-pb", "ava-ont", "ava-pb", "sr", "map-hifi"]
_MINIMAP_SECONDARY = ["yes", "no"]
_DIAMOND_MODES = [
    "ultra-sensitive", "very-sensitive", "more-sensitive",
    "sensitive", "mid-sensitive", "fast",
]
_ORFIPY_STRANDS = ["b", "f", "r"]
_TOPTAXRANK_OPTIONS = ["auto", "division", "domain", "kingdom", "order", "family"]

# ---------------------------------------------------------------------------
# Field skip lists (replaced by custom widgets or not user-facing)
# ---------------------------------------------------------------------------
_SKIP_FIELDS: Dict[str, set] = {
    "runtime":        {"temp_dir"},           # shown as dir-picker, add manually
    "chopper":        {"contam"},             # rarely used
    "porechop_abi":   set(),                  # split into basic / advanced
    "bgfilter":       {"bgfilter"},           # replaced by HostGenomePicker
    "annotation":     {"anns", "anns_options", "ann1", "ann2"},  # handled specially
    "blastn":         set(),
    "minimap":        {"evalue"},             # minimap2 doesn't use e-value
    "diamond_blastp": set(),
    "orfipy":         set(),
    "assembly":       {"refseqs", "refaccs", "vcf_flt"},  # job-level / advanced
    "em":             set(),
    "output":         set(),
}

# Porechop fields shown only in the "Advanced" sub-section
_PABI_ADVANCED = frozenset({
    "number_of_run", "consensus_run", "all_above_x", "best_of_x",
    "discard_database", "custom_adapters",
    "end_size", "min_trim_size", "extra_end_trim",
    "extra_middle_trim_good_side", "extra_middle_trim_bad_side",
    "min_split_read_size", "no_split", "discard_middle",
    "scoring_scheme", "window_size",
})

# Fields where Path should open a directory picker (not a file picker)
_DIR_PATH_FIELDS = frozenset({"runtime.output_root", "runtime.temp_dir"})

# Display label overrides: "section.field" → human label
_FIELD_LABELS: Dict[str, str] = {
    "runtime.num_threads":        "Threads",
    "runtime.output_root":        "Default output directory",
    "runtime.temp_dir":           "Temp directory",
    "runtime.entrez_email":       "NCBI Entrez email",
    "runtime.entrez_api_key":     "NCBI API key (optional)",
    "runtime.mode":               "Mode (online/offline)",
    "assembly.database":          "Assembly reference DB",
    "chopper.quality":            "Min quality (Phred)",
    "chopper.minlength":          "Min length (bp)",
    "chopper.maxlength":          "Max length (bp, 0 = no limit)",
    "chopper.headcrop":           "Head crop (bp)",
    "chopper.tailcrop":           "Tail crop (bp)",
    "chopper.trim_approach":      "Trim approach",
    "chopper.cutoff":             "Quality cutoff",
    "porechop_abi.verbosity":     "Verbosity (0-3)",
    "porechop_abi.ab_initio":     "Ab initio adapter detection",
    "porechop_abi.adapter_threshold": "Adapter identity threshold (%)",
    "porechop_abi.check_reads":   "Reads sampled for detection",
    "porechop_abi.end_threshold": "End adapter identity threshold (%)",
    "porechop_abi.middle_threshold": "Middle adapter threshold (%)",
    "porechop_abi.window_size":   "Smoothing window",
    "porechop_abi.number_of_run": "Initial consensus iterations",
    "porechop_abi.consensus_run": "Extra consensus runs",
    "porechop_abi.all_above_x":  "Min adapter % for consensus",
    "porechop_abi.best_of_x":    "Keep best N seqs (0=all)",
    "porechop_abi.discard_database": "Ignore built-in adapter DB",
    "porechop_abi.custom_adapters": "Custom adapter file",
    "porechop_abi.end_size":     "End search size (bp)",
    "porechop_abi.min_trim_size": "Min adapter alignment (bp)",
    "porechop_abi.extra_end_trim": "Extra bases trimmed at ends",
    "porechop_abi.extra_middle_trim_good_side": "Extra trim — good side",
    "porechop_abi.extra_middle_trim_bad_side":  "Extra trim — bad side",
    "porechop_abi.min_split_read_size": "Min post-split read size",
    "porechop_abi.no_split":     "Skip middle-adapter splitting",
    "porechop_abi.discard_middle": "Discard reads with middle adapters",
    "porechop_abi.scoring_scheme": "Scoring scheme (match,mm,go,ge)",
    "bgfilter.minimap.c":              "Full DP alignment + CIGAR (-c, ~2-5× slower)",
    "bgfilter.minimap.m":              "Chaining score threshold (-m)",
    "bgfilter.minimap.s":              "DP alignment score threshold (-s, requires -c)",
    "bgfilter.minimap.save_mmi_index": "Save/reuse minimap2 index (.mmi)",
    "annotation.ann1":           "Round 1 databases (CSV of keys)",
    "annotation.ann2":           "Round 2 databases (CSV of keys)",
    "blastn.bitscore":           "Min bitscore",
    "blastn.evalue":             "Max e-value",
    "blastn.max_target_seqs":    "Max target sequences",
    "blastn.save_hits":          "Save raw BLAST output to disk (for debugging)",
    "minimap.s":                 "Min DP alignment score (-s, only effective with -c)",
    "minimap.m":                 "Chain score threshold (-m, always applied)",
    "minimap.f":                 "Ignore top fraction of most frequent minimizers (-f, 0.0002–0.0010)",
    "minimap.c":                 "Full DP alignment (-c, enables AS:i: score + accurate coverage, ~2-5× slower)",
    "minimap.k":                 "K-mer size (-k, 0 = preset default ~15; use 11–13 for divergent viruses)",
    "minimap.w":                 "Minimizer window size (-w, 0 = preset default; use 5–7 when lowering k)",
    "minimap.preset":            "Preset",
    "minimap.secondary":         "Report secondary alignments",
    "minimap.max_target_seqs":   "Max target sequences",
    "minimap.save_hits":         "Save raw PAF to disk (for debugging)",
    "diamond_blastp.bitscore":   "Min bitscore",
    "diamond_blastp.evalue":     "Max e-value",
    "diamond_blastp.mode":       "Sensitivity mode",
    "diamond_blastp.max_target_seqs": "Max target sequences",
    "diamond_blastp.save_hits":  "Save raw Diamond output to disk (for debugging)",
    "orfipy.min_len":            "Min ORF length (nt)",
    "orfipy.max_len":            "Max ORF length (nt, 0 = no limit)",
    "orfipy.table":              "Genetic code table",
    "orfipy.partial_5":          "Allow partial 5′ ORFs",
    "orfipy.partial_3":          "Allow partial 3′ ORFs",
    "orfipy.strand":             "Strand",
    "assembly.min_coverage":           "Min depth for consensus",
    "assembly.min_af":                 "Min AF for the consensus base (dominant allele)",
    "assembly.report_min_af":          "Min AF to list a variant in the report",
    "assembly.min_ac":                 "Min allele count",
    "assembly.minimap.c":              "Full DP alignment + CIGAR (-c, ~2-5× slower)",
    "assembly.minimap.m":              "Chaining score threshold (-m, 0 = minimap2 default)",
    "assembly.minimap.s":              "DP alignment score threshold (-s, requires -c, 0 = default)",
    "assembly.minimap.k":              "k-mer size (-k, 0 = default 15)",
    "assembly.minimap.w":              "Minimizer window (-w, 0 = default 10)",
    "em.temp":                   "Temperature",
    "em.max_iter":               "Max iterations",
    "em.log_L_threshold":        "Log-likelihood threshold",
    "em.exclude_no_species":     "Exclude reads with no species hit",
    "report.exclude_no_species": "Exclude taxa with no species-level assignment",
    "output.pack_reads":         "Pack all reads in archive",
    "output.pack_reads_vi":      "Pack viral reads in archive",
    "output.fasta_wrap":         "Wrap FASTA sequences",
    "output.fasta_line_width":   "FASTA line width",
    "output.realtime_abund_table": "Write abundance table after each batch",
    "report.min_qcov":           "Min query coverage (%)",
    "report.min_bitscore":       "Min bitscore",
    "report.export_excel":       "Export tables as Excel (.xlsx)",
    "report.sort_reads":         "Sort reads by taxonomy",
    "report.toptaxrank":         "Top taxonomy rank for read sorting",
}


# ===========================================================================
# Collapsible section widget (same pattern as ViRefSet)
# ===========================================================================

class _CollapsibleSection(QWidget):
    """Clickable triangle header + collapsible content area."""

    def __init__(self, label: str, parent=None) -> None:
        super().__init__(parent)
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 2, 0, 2)
        outer.setSpacing(0)

        self._label = label
        self._expanded = False

        self._toggle_btn = QPushButton(f"{_COLLAPSED}  {label}")
        self._toggle_btn.setFlat(True)
        self._toggle_btn.setStyleSheet(
            "text-align: left; font-weight: bold; padding: 4px 6px;"
        )
        self._toggle_btn.clicked.connect(self._on_toggle)

        self._content = QWidget()
        self._content.setVisible(False)

        outer.addWidget(self._toggle_btn)
        outer.addWidget(self._content)

    def set_content_layout(self, layout) -> None:
        self._content.setLayout(layout)

    def set_content_widget(self, widget: QWidget) -> None:
        vbox = QVBoxLayout(self._content)
        vbox.setContentsMargins(0, 0, 0, 0)
        vbox.setSpacing(0)
        vbox.addWidget(widget)

    def expand(self) -> None:
        if not self._expanded:
            self._on_toggle()

    def _on_toggle(self) -> None:
        self._expanded = not self._expanded
        symbol = _EXPANDED if self._expanded else _COLLAPSED
        self._toggle_btn.setText(f"{symbol}  {self._label}")
        self._content.setVisible(self._expanded)


# ===========================================================================
# Special-purpose sub-widgets
# ===========================================================================

class _HostGenomePicker(QWidget):
    """Checkable list of host genome keys (bgfilter.bgfilter CSV)."""

    def __init__(
        self,
        host_db_keys: List[str],
        selected_csv: str,
        parent=None,
    ) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        selected = {k.strip() for k in selected_csv.split(",") if k.strip()}

        self._list = QListWidget()
        self._list.setMaximumHeight(120)
        self._list.setAlternatingRowColors(True)

        if not host_db_keys:
            placeholder = QListWidgetItem("(no host genomes configured)")
            placeholder.setFlags(Qt.ItemFlag.NoItemFlags)
            self._list.addItem(placeholder)
        else:
            for key in sorted(host_db_keys):
                item = QListWidgetItem(key)
                item.setCheckState(
                    Qt.CheckState.Checked if key in selected else Qt.CheckState.Unchecked
                )
                self._list.addItem(item)

        layout.addWidget(self._list)
        if not host_db_keys:
            note = QLabel("Install host genomes via Tools \u203a Database Manager.")
            note.setStyleSheet("color: #999; font-size: 11px;")
            layout.addWidget(note)

    def get_value(self) -> str:
        """Return comma-separated string of checked host genome keys."""
        checked = []
        for i in range(self._list.count()):
            item = self._list.item(i)
            if item.checkState() == Qt.CheckState.Checked:
                checked.append(item.text())
        return ",".join(checked)


def _make_section_label(text: str) -> QLabel:
    """Return a bold QLabel used as a visual sub-section header in a QFormLayout."""
    lbl = QLabel(text)
    font = lbl.font()
    font.setBold(True)
    lbl.setFont(font)
    return lbl


class _RefseqPathWidget(QWidget):
    """Single-file path selector for a reference FASTA file.

    Rendered as a QLineEdit + Browse button, matching the standard
    path-selector pattern used elsewhere in the form.
    """

    def __init__(self, value: str = "", parent=None) -> None:
        super().__init__(parent)
        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)

        self._edit = QLineEdit(value)
        self._edit.setPlaceholderText("Path to reference FASTA (.fa / .fna / .faa / .fasta)")
        self._edit.setToolTip("Reference FASTA file for reference-based assembly")

        browse_btn = QPushButton("Browse…")
        browse_btn.setFixedWidth(72)
        browse_btn.clicked.connect(self._on_browse)

        row.addWidget(self._edit, stretch=1)
        row.addWidget(browse_btn)

    def _on_browse(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Select Reference FASTA",
            self._edit.text() or "",
            "FASTA files (*.fa *.faa *.fna *.fasta *.fa.gz *.faa.gz *.fna.gz *.fasta.gz);;All files (*)",
        )
        if path:
            self._edit.setText(path)

    def get_value(self) -> str:
        return self._edit.text().strip()


class _ViralHitsCombo(QComboBox):
    """QComboBox pre-populated with viral sseqid values from annotation results.

    Reads ``annotation/annot2.tsv`` then ``annotation/annot1.tsv`` (first found),
    filters rows where ``division == 'Viruses'``, and lists unique ``sseqid``
    values as options.  Disabled when no annotation file exists or when the
    annotation is SILVA (detected via presence of ``domain`` / ``organism_name``
    columns, absence of ``division``).

    Selecting an item copies the accession into the *refaccs* QLineEdit partner.
    """

    def __init__(self, annot_dir: Optional[Path] = None,
                 refaccs_edit: Optional[QLineEdit] = None,
                 parent=None) -> None:
        super().__init__(parent)
        self._refaccs_edit = refaccs_edit
        self.setToolTip("Viral sequence IDs from annotation results (NCBI only)")
        self.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
        self.setMinimumWidth(180)

        self.currentTextChanged.connect(self._on_selected)
        self.load_from_annot_dir(annot_dir)

    def load_from_annot_dir(
        self,
        annot_dir: Optional[Path],
    ) -> None:
        """Reload combo entries from *annot_dir*.  Pass ``None`` to disable.

        Prefers ``assembly_candidates.tsv`` (written by the assembly-DB search
        during annotation). Falls back to ``annot2.tsv`` / ``annot1.tsv``.
        """
        self.blockSignals(True)
        self.clear()
        self.setEnabled(False)

        if annot_dir is None:
            self.addItem("(no annotation)")
            self.blockSignals(False)
            return

        tsv: Optional[Path] = None
        for name in ("assembly_candidates.tsv", "annot2.tsv", "annot1.tsv"):
            p = annot_dir / name
            if p.exists() and p.stat().st_size > 0:
                tsv = p
                break

        if tsv is None:
            self.addItem("(no annotation)")
            self.blockSignals(False)
            return

        try:
            df = pd.read_csv(tsv, sep="\t", dtype=str, keep_default_na=False)
        except Exception:
            self.addItem("(read error)")
            self.blockSignals(False)
            return

        if is_silva_annotation(df):
            self.addItem("(SILVA — N/A)")
            self.blockSignals(False)
            return

        if "division" not in df.columns or "sseqid" not in df.columns:
            self.addItem("(no viral data)")
            self.blockSignals(False)
            return

        viral = df[df["division"].str.strip().str.lower() == "viruses"]
        viral = viral[viral["sseqid"].str.strip() != ""]
        # Most-supported accessions first: the candidate table keeps secondary
        # hits, so it lists more accessions than reads' top hits alone.
        if "qseqid" in viral.columns:
            counts = viral.groupby(viral["sseqid"].str.strip())["qseqid"].nunique()
            ids = counts.sort_values(ascending=False, kind="stable").index.tolist()
        else:
            ids = sorted({s.strip() for s in viral["sseqid"]})

        if not ids:
            self.addItem("(no viral hits)")
            self.blockSignals(False)
            return

        self.addItem("")  # blank = no selection
        for sid in ids:
            self.addItem(sid)
        self.setEnabled(True)
        self.blockSignals(False)

    def _on_selected(self, text: str) -> None:
        if self._refaccs_edit is not None and text and not text.startswith("("):
            self._refaccs_edit.setText(text)


def _make_maxlen_spinbox(value: Any) -> QSpinBox:
    """Return a QSpinBox for a max-length field where 0 means "no limit".

    ``float("inf")`` and ``None`` are both shown as 0.
    On apply, 0 is converted back to ``float("inf")``.
    """
    sb = QSpinBox()
    sb.setRange(0, 999_999_999)
    sb.setSpecialValueText("0 (no limit)")
    display = 0 if (value is None or value == float("inf")) else int(value)
    sb.setValue(display)
    return sb


_ANN2_TARGETS = [
    ("vi",  "Viral"),
    ("ph",  "Phage"),
    ("ab",  "Bacterial"),
    ("un",  "Unannotated"),
]


class _Ann1Picker(QWidget):
    """Ordered list of annotation database keys for Round 1.

    Each database may appear at most once; order matters for pipeline execution.
    Emits ``changed`` when the list is modified.
    """

    changed = pyqtSignal()

    def __init__(
        self,
        databases: Dict,   # Dict[str, DatabaseEntry]
        installed: set,
        ann1_csv: str,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self._databases = databases
        self._installed = installed

        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(4)

        # --- Add row: combo of remaining DBs + Add button ---
        add_row = QHBoxLayout()
        self._add_combo = QComboBox()
        self._add_combo.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        add_btn = QPushButton("Add")
        add_btn.setFixedWidth(48)
        add_btn.clicked.connect(self._add_item)
        add_row.addWidget(self._add_combo)
        add_row.addWidget(add_btn)
        root.addLayout(add_row)

        # --- List widget + reorder/remove buttons on the right ---
        list_row = QHBoxLayout()
        self._list = QListWidget()
        self._list.setMaximumHeight(100)
        list_row.addWidget(self._list, stretch=1)

        btn_col = QVBoxLayout()
        self._up_btn   = QPushButton("↑")
        self._down_btn = QPushButton("↓")
        self._del_btn  = QPushButton("✕")
        for b in (self._up_btn, self._down_btn, self._del_btn):
            b.setFixedWidth(28)
            btn_col.addWidget(b)
        btn_col.addStretch()
        list_row.addLayout(btn_col)
        root.addLayout(list_row)

        self._up_btn.clicked.connect(self._move_up)
        self._down_btn.clicked.connect(self._move_down)
        self._del_btn.clicked.connect(self._remove_item)

        self.set_value(ann1_csv)

    def _refresh_add_combo(self) -> None:
        """Rebuild add combo with DBs not yet in the list."""
        current = set(self._get_items())
        self._add_combo.clear()
        mid = self._add_combo.palette().color(
            self._add_combo.palette().ColorRole.Mid
        )
        for key in sorted(self._databases.keys()):
            if key in current:
                continue
            self._add_combo.addItem(key, key)
            if key not in self._installed:
                item = self._add_combo.model().item(self._add_combo.count() - 1)
                if item is not None:
                    item.setEnabled(False)
                    item.setForeground(mid)

    def _get_items(self) -> list[str]:
        return [self._list.item(i).text() for i in range(self._list.count())]

    def _add_item(self) -> None:
        key = self._add_combo.currentData()
        if not key or key in self._get_items():
            return
        self._list.addItem(key)
        self._refresh_add_combo()
        self.changed.emit()

    def _move_up(self) -> None:
        row = self._list.currentRow()
        if row <= 0:
            return
        item = self._list.takeItem(row)
        self._list.insertItem(row - 1, item)
        self._list.setCurrentRow(row - 1)
        self.changed.emit()

    def _move_down(self) -> None:
        row = self._list.currentRow()
        if row < 0 or row >= self._list.count() - 1:
            return
        item = self._list.takeItem(row)
        self._list.insertItem(row + 1, item)
        self._list.setCurrentRow(row + 1)
        self.changed.emit()

    def _remove_item(self) -> None:
        row = self._list.currentRow()
        if row < 0:
            return
        self._list.takeItem(row)
        self._refresh_add_combo()
        self.changed.emit()

    def get_value(self) -> str:
        return ",".join(self._get_items())

    def set_value(self, csv: str) -> None:
        self._list.clear()
        for key in csv.split(","):
            key = key.strip()
            if key and key in self._databases:
                self._list.addItem(key)
        self._refresh_add_combo()


class _Ann2Picker(QWidget):
    """List of ``target:database`` pairs for Round 2 annotation.

    The same target and/or database may appear multiple times.
    Emits ``changed`` when the list is modified.
    """

    changed = pyqtSignal()

    def __init__(
        self,
        databases: Dict,   # Dict[str, DatabaseEntry]
        installed: set,
        ann2_csv: str,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self._databases = databases
        self._installed = installed

        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(4)

        # --- Add row: target combo + DB combo + Add button ---
        add_row = QHBoxLayout()

        self._target_combo = QComboBox()
        for code, label in _ANN2_TARGETS:
            self._target_combo.addItem(f"{label} ({code})", code)
        self._target_combo.setFixedWidth(140)

        self._db_combo = QComboBox()
        self._db_combo.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        mid = self._db_combo.palette().color(
            self._db_combo.palette().ColorRole.Mid
        )
        for key in sorted(databases.keys()):
            self._db_combo.addItem(key, key)
            if key not in installed:
                item = self._db_combo.model().item(self._db_combo.count() - 1)
                if item is not None:
                    item.setEnabled(False)
                    item.setForeground(mid)

        add_btn = QPushButton("Add")
        add_btn.setFixedWidth(48)
        add_btn.clicked.connect(self._add_item)

        add_row.addWidget(self._target_combo)
        add_row.addWidget(self._db_combo)
        add_row.addWidget(add_btn)
        root.addLayout(add_row)

        # --- List widget + Remove button on the right ---
        list_row = QHBoxLayout()
        self._list = QListWidget()
        self._list.setMaximumHeight(100)
        list_row.addWidget(self._list, stretch=1)

        del_btn = QPushButton("✕")
        del_btn.setFixedWidth(28)
        del_btn.clicked.connect(self._remove_item)
        btn_col = QVBoxLayout()
        btn_col.addWidget(del_btn)
        btn_col.addStretch()
        list_row.addLayout(btn_col)
        root.addLayout(list_row)

        self.set_value(ann2_csv)

    def _add_item(self) -> None:
        target = self._target_combo.currentData()
        db = self._db_combo.currentData()
        if not target or not db:
            return
        item = QListWidgetItem(f"{target} : {db}")
        item.setData(Qt.ItemDataRole.UserRole, (target, db))
        self._list.addItem(item)
        self.changed.emit()

    def _remove_item(self) -> None:
        row = self._list.currentRow()
        if row >= 0:
            self._list.takeItem(row)
            self.changed.emit()

    def get_value(self) -> str:
        parts = []
        for i in range(self._list.count()):
            data = self._list.item(i).data(Qt.ItemDataRole.UserRole)
            if data:
                target, db = data
                parts.append(f"{target}:{db}")
        return ",".join(parts)

    def set_value(self, csv: str) -> None:
        self._list.clear()
        for token in csv.split(","):
            token = token.strip()
            if not token:
                continue
            if ":" in token:
                target, db = token.split(":", 1)
            else:
                target, db = "vi", token   # legacy bare key → vi
            if db not in self._databases:
                continue
            item = QListWidgetItem(f"{target} : {db}")
            item.setData(Qt.ItemDataRole.UserRole, (target, db))
            self._list.addItem(item)


class _AnnotationStrategyBar(QWidget):
    """Preset selector + Round-1 ordered DB list + Round-2 target:db pair list.

    * Selecting a preset populates both Round-1 and Round-2 pickers.
    * Any manual change to Round-1 or Round-2 resets the preset to "(custom)".
    * Presets that reference uninstalled databases are shown greyed and
      non-selectable.

    Signals
    -------
    preset_selected
        The preset name chosen in the drop-down (``None`` for "(custom)").
    """

    preset_selected = pyqtSignal(object)  # Optional[str]

    def __init__(
        self,
        strategies: Dict[str, str],
        databases: Dict,           # Dict[str, DatabaseEntry]
        installed_dbs: set,
        current_anns: Optional[str],
        ann1: str,
        ann2: str,
        parent=None,
    ) -> None:
        super().__init__(parent)
        layout = QFormLayout(self)
        layout.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)

        self._strategies = strategies
        self._installed = installed_dbs
        self._blocking = False   # prevents re-entrant signal loops

        # ------------------------------------------------------------------
        # Preset dropdown
        # ------------------------------------------------------------------
        self._strategy_combo = QComboBox()
        self._strategy_combo.addItem("(custom)", None)
        for name in sorted(strategies.keys()):
            self._strategy_combo.addItem(name, name)
            idx = self._strategy_combo.count() - 1
            if not self._strategy_all_installed(strategies[name]):
                item = self._strategy_combo.model().item(idx)
                if item is not None:
                    item.setEnabled(False)
                    item.setForeground(
                        self._strategy_combo.palette().color(
                            self._strategy_combo.palette().ColorRole.Mid
                        )
                    )
        layout.addRow("Annotation preset:", self._strategy_combo)

        # ------------------------------------------------------------------
        # Round 1 — ordered database list
        # ------------------------------------------------------------------
        self._ann1_picker = _Ann1Picker(databases, installed_dbs, ann1)
        layout.addRow("Round 1 databases:", self._ann1_picker)
        self._ann1_picker.changed.connect(self._on_db_changed)

        # ------------------------------------------------------------------
        # Round 2 — list of target:db pairs
        # ------------------------------------------------------------------
        self._ann2_picker = _Ann2Picker(databases, installed_dbs, "")
        layout.addRow("Round 2 target:databases:", self._ann2_picker)
        self._ann2_picker.changed.connect(self._on_db_changed)

        # ------------------------------------------------------------------
        # Wire preset signal and load initial state
        # ------------------------------------------------------------------
        self._strategy_combo.currentIndexChanged.connect(self._on_strategy_changed)

        if current_anns and current_anns in strategies:
            i = self._strategy_combo.findData(current_anns)
            if i >= 0:
                self._strategy_combo.setCurrentIndex(i)
                return   # _on_strategy_changed populates pickers

        # No matching preset — stay at (custom), populate from raw values
        self._ann2_picker.set_value(ann2)

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _strategy_all_installed(self, option_str: str) -> bool:
        """Return True if every DB key referenced by *option_str* is installed."""
        for key in self._parse_option(option_str, "ann1").split(","):
            key = key.strip()
            if key and key not in self._installed:
                return False
        for token in self._parse_option(option_str, "ann2").split(","):
            token = token.strip()
            dbkey = token.split(":", 1)[-1] if ":" in token else token
            if dbkey and dbkey not in self._installed:
                return False
        return True

    @staticmethod
    def _parse_option(option_str: str, key: str) -> str:
        m = re.search(rf"--{key}\s+(\S+)", option_str)
        return m.group(1) if m else ""

    # ------------------------------------------------------------------
    # Slots
    # ------------------------------------------------------------------

    def _on_strategy_changed(self, _index: int) -> None:
        if self._blocking:
            return
        name = self._strategy_combo.currentData()
        self.preset_selected.emit(name)
        if name is None:
            return  # (custom) — leave pickers as-is

        option_str = self._strategies.get(name, "")
        ann1_str = self._parse_option(option_str, "ann1")
        ann2_str = self._parse_option(option_str, "ann2")

        self._blocking = True
        try:
            self._ann1_picker.set_value(ann1_str)
            self._ann2_picker.set_value(ann2_str)
        finally:
            self._blocking = False

    def _on_db_changed(self) -> None:
        if self._blocking:
            return
        self._blocking = True
        try:
            self._strategy_combo.setCurrentIndex(0)  # "(custom)"
        finally:
            self._blocking = False

    # ------------------------------------------------------------------
    # Read-out API (called by _build_config)
    # ------------------------------------------------------------------

    def get_strategy(self) -> Optional[str]:
        return self._strategy_combo.currentData()   # None for (custom)

    def get_ann1(self) -> str:
        return self._ann1_picker.get_value()

    def get_ann2(self) -> str:
        return self._ann2_picker.get_value()


# ===========================================================================
# Main config editor
# ===========================================================================

class ConfigEditorWidget(QScrollArea):
    """Scrollable, collapsible-section config editor.

    Signals
    -------
    config_changed
        Emitted with the new :class:`EMtraceConfig` after Apply is clicked,
        and after a switch of analysis mode.
    status_message
        A one-line message for the main window's status bar.
    """

    config_changed = pyqtSignal(object)  # EMtraceConfig
    status_message = pyqtSignal(str)

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setWidgetResizable(True)

        self._config: Optional[EMtraceConfig] = None
        self._field_widgets: Dict[str, QWidget] = {}
        self._host_picker: Optional[_HostGenomePicker] = None        # legacy ref kept for compat
        self._host_picker_mm2: Optional[_HostGenomePicker] = None  # minimap2 host picker
        self._strategy_bar: Optional[_AnnotationStrategyBar] = None
        self._annot_dir: Optional[Path] = None
        self._viral_hits_combo: Optional[_ViralHitsCombo] = None
        self._mode_combo: Optional[QComboBox] = None

        # Outer container
        self._container = QWidget()
        self._vbox = QVBoxLayout(self._container)
        self._vbox.setContentsMargins(10, 10, 10, 10)
        self._vbox.setSpacing(4)

        # Empty-state label
        self._empty_label = QLabel("No configuration loaded.")
        self._empty_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._empty_label.setStyleSheet("color: #999; font-size: 14px;")
        self._vbox.addWidget(self._empty_label)
        self._vbox.addStretch(1)

        self.setWidget(self._container)

    # ==================================================================
    # Public API
    # ==================================================================

    def load_config(self, config: EMtraceConfig) -> None:
        """Populate the form from *config*."""
        self._config = config
        self._clear()

        ann = config.annotation

        # --- Analysis mode: virome / 16S ---
        self._vbox.addWidget(self._build_mode_row(config))

        # --- Section: Runtime ---
        self._vbox.addWidget(self._build_section("runtime", "Runtime", config.runtime))

        # --- Section: Preprocessing (Porechop → Chopper → Background Filtering) ---
        self._vbox.addWidget(self._build_preprocess_section(config.preprocess, config))

        # --- Section: Annotation ---
        # Determine which annotation DBs are installed for the dropdown
        _installed_ann_dbs = {
            key for key in config.databases
            if config.is_db_installed(key)
        }
        self._vbox.addWidget(
            self._build_annotation_section(
                ann,
                config.ann_strategies,
                databases=config.databases,
                installed_dbs=_installed_ann_dbs,
            )
        )

        # --- Section: Assembly ---
        self._vbox.addWidget(self._build_assembly_section(config.assembly))

        # --- Section: EM Abundance ---
        self._vbox.addWidget(self._build_section("em", "EM Abundance", config.em))

        # --- Section: Report ---
        self._vbox.addWidget(self._build_report_section(config.report))

        # --- Section: Output ---
        self._vbox.addWidget(self._build_section("output", "Output Format", config.output))

        # --- Runtime temp_dir (special: dir-picker, added after auto sections) ---
        # Already handled inside _build_section via dir-mode detection

        # --- Apply button ---
        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        sep.setStyleSheet("color: #ddd;")
        self._vbox.addWidget(sep)

        btn_row = QHBoxLayout()
        self._apply_btn = QPushButton("Apply Changes")
        self._apply_btn.setFixedWidth(140)
        self._apply_btn.clicked.connect(self._on_apply)
        self._status_label = QLabel("")
        self._status_label.setStyleSheet("color: #888; font-size: 11px;")
        btn_row.addWidget(self._apply_btn)
        btn_row.addWidget(self._status_label)
        btn_row.addStretch(1)
        self._vbox.addLayout(btn_row)
        self._vbox.addStretch(1)

        self._empty_label.setVisible(False)

    def set_output_dir(self, output_dir: Optional[Path]) -> None:
        """Notify the editor of the current job's output directory.

        Updates the viral hits combo with the annotation directory so that
        viral sseqid values can be offered as reference accession shortcuts.
        """
        self._annot_dir = output_dir / "annotation" if output_dir is not None else None
        if self._viral_hits_combo is not None:
            self._viral_hits_combo.load_from_annot_dir(self._annot_dir)

    def _on_assembly_db_changed(self, combo: "QComboBox", new_idx: int) -> None:
        """Confirm before changing the assembly reference database."""
        prev_idx = combo.property("_prev_index") or 0
        if new_idx == prev_idx:
            return
        new_name = combo.itemText(new_idx) or "(none)"
        answer = QMessageBox.question(
            self,
            "Change assembly reference database?",
            f"The assembly reference database is used to look up viral reference\n"
            f"genomes for reference-based assembly. In most cases 'minimap.assembly'\n"
            f"is the correct choice.\n\n"
            f"Change to: {new_name}\n\n"
            f"Are you sure?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        if answer == QMessageBox.StandardButton.Yes:
            combo.setProperty("_prev_index", new_idx)
        else:
            combo.blockSignals(True)
            combo.setCurrentIndex(prev_idx)
            combo.blockSignals(False)

    # ==================================================================
    # Analysis mode (virome / 16S)
    # ==================================================================

    def _build_mode_row(self, config: EMtraceConfig) -> QWidget:
        """Mode selector at the top of the form: virome or 16S settings."""
        row = QWidget()
        form = QFormLayout(row)
        form.setContentsMargins(0, 0, 0, 6)
        self._mode_combo = QComboBox()
        for label, mode in _MODE_LABELS.items():
            self._mode_combo.addItem(label, mode)
        self._mode_combo.setToolTip(
            "Virome: default settings for viral metagenomics.\n"
            "16S rRNA: the benchmarked 16S profile (full-length amplicons;\n"
            "-N 25, EM temperature 10, no host filter, 1200-1800 bp reads,\n"
            "1 % tail trim). Switching changes only these settings."
        )
        self._mode_combo.setCurrentIndex(self._mode_combo.findData(config.mode))
        self._mode_combo.currentIndexChanged.connect(self._on_mode_selected)
        form.addRow("Analysis mode:", self._mode_combo)
        return row

    def _on_mode_selected(self, _index: int) -> None:
        mode = self._mode_combo.currentData()
        # Deferred: the switch rebuilds the form, including this combo.
        QTimer.singleShot(0, lambda: self._switch_mode(mode, "mode switch"))

    def _on_preset_selected(self, name: Optional[str]) -> None:
        """A 16S.* preset switches to 16S mode, any other preset back to virome."""
        if name is None or self._mode_combo is None:
            return
        mode = "16S" if name.startswith("16S.") else "virome"
        if mode != self._mode_combo.currentData():
            QTimer.singleShot(0, lambda: self._switch_mode(mode, f"preset {name}"))

    def _switch_mode(self, mode: str, trigger: str) -> None:
        """Apply *mode* to the current form values and reload the form."""
        if self._config is None:
            return
        try:
            new_cfg, changes = self._build_config().with_mode(mode)
        except Exception as exc:  # noqa: BLE001
            logger.error("Mode switch failed: %s", exc)
            self.status_message.emit(f"Could not switch to {mode} mode: {exc}")
            return
        self.load_config(new_cfg)
        self.config_changed.emit(new_cfg)
        label = next(k for k, v in _MODE_LABELS.items() if v == mode)
        if changes:
            details = ", ".join(
                f"{_MODE_KEY_LABELS.get(key, key)} {_fmt_mode_value(old, key)}→{_fmt_mode_value(new, key)}"
                for key, old, new in changes
            )
            msg = f"{label} mode ({trigger}): {details}"
        else:
            msg = f"{label} mode ({trigger}): settings already match"
        self._status_label.setText(f"{label} mode applied.")
        self._status_label.setStyleSheet("color: #16a34a; font-size: 11px;")
        self.status_message.emit(msg)

    # ==================================================================
    # Apply
    # ==================================================================

    def _on_apply(self) -> None:
        """Read all form fields → emit config_changed(new_config)."""
        if self._config is None:
            return
        try:
            new_cfg = self._build_config()
            self._config = new_cfg
            self.config_changed.emit(new_cfg)
            self._status_label.setText("Applied.")
            self._status_label.setStyleSheet("color: #16a34a; font-size: 11px;")
        except Exception as exc:  # noqa: BLE001
            logger.error("Config apply failed: %s", exc)
            self._status_label.setText(f"Error: {exc}")
            self._status_label.setStyleSheet("color: #dc2626; font-size: 11px;")

    def _build_config(self) -> EMtraceConfig:
        """Reconstruct EMtraceConfig from form widgets."""
        cfg = self._config

        # Runtime
        rt_data = cfg.runtime.model_dump()
        for k in ("num_threads", "output_root", "entrez_email", "entrez_api_key", "mode"):
            w = self._field_widgets.get(f"runtime.{k}")
            if w is not None:
                rt_data[k] = self._read_widget_value(w)
        # temp_dir special: dir-picker container
        temp_w = self._field_widgets.get("runtime.temp_dir")
        if temp_w is not None:
            rt_data["temp_dir"] = self._read_widget_value(temp_w) or "$TMPDIR"
        runtime = RuntimeConfig(**rt_data)

        # Chopper
        ch_data = cfg.preprocess.chopper.model_dump()
        for field in ("quality", "minlength", "headcrop", "tailcrop", "cutoff", "trim_approach"):
            w = self._field_widgets.get(f"chopper.{field}")
            if w is not None:
                ch_data[field] = self._read_widget_value(w)
        # maxlength: spinbox value 0 means no limit → float("inf")
        maxlen_w = self._field_widgets.get("chopper.maxlength")
        if maxlen_w is not None:
            raw = self._read_widget_value(maxlen_w)
            ch_data["maxlength"] = float("inf") if raw == 0 else raw
        chopper = ChopperConfig(**ch_data)

        # Porechop_ABI
        pabi_data = cfg.preprocess.porechop_abi.model_dump()
        for field_name in PorechopABIConfig.model_fields:
            w = self._field_widgets.get(f"porechop_abi.{field_name}")
            if w is not None:
                pabi_data[field_name] = self._read_widget_value(w)
        porechop_abi = PorechopABIConfig(**pabi_data)

        # BgFilter (nested under preprocess)
        bgf_data = cfg.preprocess.bgfilter.model_dump()
        _selected_hosts: list[str] = []
        if self._host_picker_mm2 is not None:
            _selected_hosts.extend(v for v in self._host_picker_mm2.get_value().split(",") if v)
        if _selected_hosts:
            bgf_data["bgfilter"] = ",".join(_selected_hosts)
        elif self._host_picker is not None:
            bgf_data["bgfilter"] = self._host_picker.get_value()
        # Nested minimap sub-fields
        for field in ("c", "m", "s", "save_mmi_index"):
            w = self._field_widgets.get(f"bgfilter.minimap.{field}")
            if w is not None:
                bgf_data["minimap"][field] = self._read_widget_value(w)
        bgfilter = BgFilterConfig(**bgf_data)

        # Preprocessing subsample
        pre_sub_data = cfg.preprocess.subsample.model_dump()
        pre_sub_check_w = self._field_widgets.get("pre.subsample.active")
        pre_sub_combo_w = self._field_widgets.get("pre.subsample.size")
        if pre_sub_check_w is not None:
            pre_sub_data["active"] = self._read_widget_value(pre_sub_check_w)
        if pre_sub_combo_w is not None:
            pre_sub_data["size"] = self._read_widget_value(pre_sub_combo_w)
        pre_subsample = SubsampleConfig(**pre_sub_data)

        preprocess = PreprocessConfig(
            subsample=pre_subsample,
            porechop_abi=porechop_abi,
            chopper=chopper,
            bgfilter=bgfilter,
        )

        # Annotation
        ann_data = cfg.annotation.model_dump()
        if self._strategy_bar is not None:
            ann_data["anns"] = self._strategy_bar.get_strategy()
            ann_data["ann1"] = self._strategy_bar.get_ann1()
            ann_data["ann2"] = self._strategy_bar.get_ann2()
        for sub, model in (("blastn", BlastnConfig), ("minimap", MinimapConfig),
                           ("diamond_blastp", DiamondBlastpConfig), ("orfipy", OrfipyConfig)):
            sub_data = ann_data[sub]
            for field_name in model.model_fields:
                w = self._field_widgets.get(f"{sub}.{field_name}")
                if w is not None:
                    sub_data[field_name] = self._read_widget_value(w)
        # orfipy max_len: spinbox value 0 means no limit → float("inf")
        orf_maxlen_w = self._field_widgets.get("orfipy.max_len")
        if orf_maxlen_w is not None:
            raw = self._read_widget_value(orf_maxlen_w)
            ann_data["orfipy"]["max_len"] = float("inf") if raw == 0 else raw
        annotation = AnnotationConfig(**ann_data)

        # Assembly
        asm_data = cfg.assembly.model_dump()
        for field_name in ("min_coverage", "min_ac", "refaccs", "database"):
            w = self._field_widgets.get(f"assembly.{field_name}")
            if w is not None:
                asm_data[field_name] = self._read_widget_value(w)
        refseqs_w = self._field_widgets.get("assembly.refseqs")
        if refseqs_w is not None and hasattr(refseqs_w, "get_value"):
            asm_data["refseqs"] = refseqs_w.get_value()
        # min_af from combo (typed float data)
        af_w = self._field_widgets.get("assembly.min_af")
        if af_w is not None:
            asm_data["min_af"] = self._read_widget_value(af_w)
        raf_w = self._field_widgets.get("assembly.report_min_af")
        if raf_w is not None:
            asm_data["report_min_af"] = self._read_widget_value(raf_w)
        # Nested minimap sub-fields
        for field_name in ("c", "m", "s", "k", "w"):
            w = self._field_widgets.get(f"assembly.minimap.{field_name}")
            if w is not None:
                asm_data["minimap"][field_name] = self._read_widget_value(w)
        assembly = AssemblyConfig(**asm_data)

        # EM
        em_data = cfg.em.model_dump()
        for field_name in EMConfig.model_fields:
            w = self._field_widgets.get(f"em.{field_name}")
            if w is not None:
                em_data[field_name] = self._read_widget_value(w)
        em = EMConfig(**em_data)

        # Report
        rep_data = cfg.report.model_dump()
        for field_name in ReportConfig.model_fields:
            w = self._field_widgets.get(f"report.{field_name}")
            if w is not None:
                rep_data[field_name] = self._read_widget_value(w)
        report = ReportConfig(**rep_data)

        # Output
        out_data = cfg.output.model_dump()
        for field_name in OutputConfig.model_fields:
            w = self._field_widgets.get(f"output.{field_name}")
            if w is not None:
                out_data[field_name] = self._read_widget_value(w)
        output = OutputConfig(**out_data)

        return EMtraceConfig(
            preprocess=preprocess,
            annotation=annotation,
            assembly=assembly,
            runtime=runtime,
            output=output,
            taxonomy=cfg.taxonomy,
            databases=cfg.databases,
            host_databases=cfg.host_databases,
            ann_strategies=cfg.ann_strategies,
            em=em,
            report=report,
            db_dir=cfg.db_dir,
        )

    # ==================================================================
    # Section builders
    # ==================================================================

    def _build_section(
        self,
        attr_name: str,
        label: str,
        model: BaseModel,
    ) -> _CollapsibleSection:
        """Generic collapsible section for a flat Pydantic model."""
        section = _CollapsibleSection(label)
        form = QFormLayout()
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)

        skip = _SKIP_FIELDS.get(attr_name, set())

        for field_name, field_info in type(model).model_fields.items():
            if field_name in skip:
                continue
            key = f"{attr_name}.{field_name}"
            current_value = getattr(model, field_name, None)
            widget = self._create_widget(field_info, field_name, key, current_value)
            self._field_widgets[key] = widget
            display_label = _FIELD_LABELS.get(key) or field_name.replace("_", " ").title()
            form.addRow(display_label + ":", widget)

        section.set_content_layout(form)
        return section

    def _build_preprocess_section(
        self,
        preprocess: PreprocessConfig,
        config: EMtraceConfig,
    ) -> _CollapsibleSection:
        """Preprocessing: Porechop_ABI → Chopper → Background Filtering (in step order)."""
        section = _CollapsibleSection("Preprocessing")
        container = QWidget()
        vbox = QVBoxLayout(container)
        vbox.setContentsMargins(12, 0, 0, 0)
        vbox.setSpacing(0)

        # --- Subsampling (applied before any preprocessing step) ---
        _SUBSAMPLE_OPTIONS = [
            ("10%",    0.1),
            ("20%",    0.2),
            ("30%",    0.3),
            ("40%",    0.4),
            ("50%",    0.5),
            ("5,000",  5000.0),
            ("10,000", 10000.0),
        ]
        pre_sub = preprocess.subsample
        pre_sub_check = QCheckBox()
        pre_sub_check.setChecked(pre_sub.active)
        self._field_widgets["pre.subsample.active"] = pre_sub_check

        pre_sub_combo = QComboBox()
        for _lbl, _val in _SUBSAMPLE_OPTIONS:
            pre_sub_combo.addItem(_lbl, _val)
        _cur_pre_size = float(pre_sub.size)
        for _i in range(pre_sub_combo.count()):
            if pre_sub_combo.itemData(_i) == _cur_pre_size:
                pre_sub_combo.setCurrentIndex(_i)
                break
        pre_sub_combo.setEnabled(pre_sub.active)
        pre_sub_combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
        pre_sub_combo.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self._field_widgets["pre.subsample.size"] = pre_sub_combo

        pre_sub_check.toggled.connect(pre_sub_combo.setEnabled)

        pre_sub_row = QWidget()
        pre_sub_row_layout = QHBoxLayout(pre_sub_row)
        pre_sub_row_layout.setContentsMargins(0, 0, 0, 0)
        pre_sub_row_layout.addWidget(pre_sub_check)
        pre_sub_row_layout.addWidget(pre_sub_combo, stretch=0)

        pre_sub_form = QFormLayout()
        pre_sub_form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        pre_sub_form.addRow("Subsample reads:", pre_sub_row)
        pre_sub_form_widget = QWidget()
        pre_sub_form_widget.setLayout(pre_sub_form)
        vbox.addWidget(pre_sub_form_widget)

        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        sep.setStyleSheet("color: #eee;")
        vbox.addWidget(sep)

        # --- Step 1: Porechop_ABI ---
        pabi_sec = _CollapsibleSection("Porechop_ABI")
        pabi_outer = QWidget()
        pabi_vbox = QVBoxLayout(pabi_outer)
        pabi_vbox.setContentsMargins(12, 0, 0, 0)
        pabi_vbox.setSpacing(0)

        pa = preprocess.porechop_abi

        basic_form = QFormLayout()
        basic_form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)

        for field_name in ("ab_initio", "verbosity", "adapter_threshold",
                           "check_reads", "end_threshold", "middle_threshold"):
            key = f"porechop_abi.{field_name}"
            w = self._create_widget(
                PorechopABIConfig.model_fields[field_name],
                field_name, key,
                getattr(pa, field_name),
            )
            self._field_widgets[key] = w
            lbl = _FIELD_LABELS.get(key) or field_name.replace("_", " ").title()
            basic_form.addRow(lbl + ":", w)

        basic_widget = QWidget()
        basic_widget.setLayout(basic_form)
        pabi_vbox.addWidget(basic_widget)

        adv_sec = _CollapsibleSection("Advanced Porechop options")
        adv_form = QFormLayout()
        adv_form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)

        for field_name in sorted(_PABI_ADVANCED):
            if field_name not in PorechopABIConfig.model_fields:
                continue
            key = f"porechop_abi.{field_name}"
            w = self._create_widget(
                PorechopABIConfig.model_fields[field_name],
                field_name, key,
                getattr(pa, field_name),
            )
            self._field_widgets[key] = w
            lbl = _FIELD_LABELS.get(key) or field_name.replace("_", " ").title()
            adv_form.addRow(lbl + ":", w)

        adv_sec.set_content_layout(adv_form)
        pabi_vbox.addWidget(adv_sec)
        pabi_sec.set_content_widget(pabi_outer)
        vbox.addWidget(pabi_sec)

        # --- Step 2: Chopper ---
        chopper_sec = _CollapsibleSection("Chopper")
        ch_form = QFormLayout()
        ch_form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        ch = preprocess.chopper

        trim_combo = QComboBox()
        for opt in _CHOPPER_TRIM_APPROACHES:
            trim_combo.addItem(opt)
        idx = trim_combo.findText(ch.trim_approach)
        if idx >= 0:
            trim_combo.setCurrentIndex(idx)
        self._field_widgets["chopper.trim_approach"] = trim_combo
        ch_form.addRow("Trim approach:", trim_combo)

        for field_name in ("quality", "minlength"):
            w = self._create_widget(
                ChopperConfig.model_fields[field_name],
                field_name, f"chopper.{field_name}",
                getattr(ch, field_name),
            )
            self._field_widgets[f"chopper.{field_name}"] = w
            label = _FIELD_LABELS.get(f"chopper.{field_name}") or field_name.replace("_", " ").title()
            ch_form.addRow(label + ":", w)

        maxlen_sb = _make_maxlen_spinbox(ch.maxlength)
        self._field_widgets["chopper.maxlength"] = maxlen_sb
        ch_form.addRow(_FIELD_LABELS["chopper.maxlength"] + ":", maxlen_sb)

        for field_name in ("cutoff", "headcrop", "tailcrop"):
            w = self._create_widget(
                ChopperConfig.model_fields[field_name],
                field_name, f"chopper.{field_name}",
                getattr(ch, field_name),
            )
            self._field_widgets[f"chopper.{field_name}"] = w
            label = _FIELD_LABELS.get(f"chopper.{field_name}") or field_name.replace("_", " ").title()
            ch_form.addRow(label + ":", w)

        chopper_sec.set_content_layout(ch_form)
        vbox.addWidget(chopper_sec)

        # --- Step 3: Background Filtering ---
        hflt_sec = _CollapsibleSection("Background Filtering")
        hf_form = QFormLayout()
        hf_form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)

        bgf = preprocess.bgfilter

        mm2_keys = [k for k, e in config.host_databases.items() if e.search.lower() in ("minimap2", "minimap")]
        selected_csv = bgf.bgfilter
        selected_set = {v.strip() for v in selected_csv.split(",") if v.strip()}

        self._host_picker_mm2 = _HostGenomePicker(
            host_db_keys=mm2_keys,
            selected_csv=",".join(k for k in mm2_keys if k in selected_set),
        )
        hf_form.addRow("Minimap2 bgfilters:", self._host_picker_mm2)


        # Minimap2 sub-section
        hf_form.addRow(_make_section_label("Minimap2 options:"), QWidget())
        for key, field_info, cur_val in [
            ("bgfilter.minimap.c",              BgFilterMinimapConfig.model_fields["c"],              bgf.minimap.c),
            ("bgfilter.minimap.m",              BgFilterMinimapConfig.model_fields["m"],              bgf.minimap.m),
            ("bgfilter.minimap.s",              BgFilterMinimapConfig.model_fields["s"],              bgf.minimap.s),
            ("bgfilter.minimap.save_mmi_index", BgFilterMinimapConfig.model_fields["save_mmi_index"], bgf.minimap.save_mmi_index),
        ]:
            field_name = key.split(".")[-1]
            w = self._create_widget(field_info, field_name, key, cur_val)
            self._field_widgets[key] = w
            lbl = _FIELD_LABELS.get(key) or field_name.replace("_", " ").title()
            hf_form.addRow("  " + lbl + ":", w)


        hflt_sec.set_content_layout(hf_form)
        vbox.addWidget(hflt_sec)

        section.set_content_widget(container)
        return section

    def _build_assembly_section(self, asm: AssemblyConfig) -> _CollapsibleSection:
        """Assembly: minimap2 mapping, consensus thresholds."""
        section = _CollapsibleSection("Assembly")
        form = QFormLayout()
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)

        # Reference FASTA file (single path selector)
        refseqs_w = _RefseqPathWidget(asm.refseqs)
        self._field_widgets["assembly.refseqs"] = refseqs_w
        form.addRow("Reference fasta file:", refseqs_w)

        # Reference accessions
        refaccs_w = QLineEdit(asm.refaccs)
        refaccs_w.setPlaceholderText("e.g. NC_001802,MN908947 (comma-separated)")
        refaccs_w.setToolTip("Entrez accessions for reference-based assembly (comma-separated)")
        self._field_widgets["assembly.refaccs"] = refaccs_w
        form.addRow("Reference accessions:", refaccs_w)

        # Viral hits picker — populated from annotation results
        viral_combo = _ViralHitsCombo(
            annot_dir=self._annot_dir,
            refaccs_edit=refaccs_w,
        )
        self._viral_hits_combo = viral_combo
        form.addRow("Pick from viral hits:", viral_combo)

        # Assembly reference database — source for local ref-genome extraction
        db_combo = QComboBox()
        db_combo.setToolTip(
            "Minimap2 database used to look up viral reference genomes for assembly.\n"
            "Sequences are extracted from its FASTA sidecar; the .mmi index is used\n"
            "for the annotation candidate search after each annotation run.\n"
            "Leave blank to always fetch references from NCBI Entrez."
        )
        db_combo.addItem("")  # blank = use Entrez only
        if self._config is not None:
            for name, entry in self._config.databases.items():
                if str(entry.search).lower() in ("minimap", "minimap2"):
                    db_combo.addItem(name)
        idx = db_combo.findText(asm.database)
        if idx >= 0:
            db_combo.setCurrentIndex(idx)
        db_combo.setProperty("_prev_index", db_combo.currentIndex())
        db_combo.currentIndexChanged.connect(
            lambda new_idx, c=db_combo: self._on_assembly_db_changed(c, new_idx)
        )
        self._field_widgets["assembly.database"] = db_combo
        form.addRow("Assembly reference DB:", db_combo)

        sep0 = QFrame()
        sep0.setFrameShape(QFrame.Shape.HLine)
        form.addRow(sep0)

        # Minimap2 sub-section (no save_mmi_index)
        form.addRow(_make_section_label("Minimap2 options:"), QWidget())
        for key, field_info, cur_val in [
            ("assembly.minimap.c", AssemblyMinimapConfig.model_fields["c"], asm.minimap.c),
            ("assembly.minimap.m", AssemblyMinimapConfig.model_fields["m"], asm.minimap.m),
            ("assembly.minimap.s", AssemblyMinimapConfig.model_fields["s"], asm.minimap.s),
            ("assembly.minimap.k", AssemblyMinimapConfig.model_fields["k"], asm.minimap.k),
            ("assembly.minimap.w", AssemblyMinimapConfig.model_fields["w"], asm.minimap.w),
        ]:
            field_name = key.split(".")[-1]
            w = self._create_widget(field_info, field_name, key, cur_val)
            self._field_widgets[key] = w
            lbl = _FIELD_LABELS.get(key) or field_name.replace("_", " ").title()
            form.addRow("  " + lbl + ":", w)

        # Separator before consensus thresholds
        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        form.addRow(sep)

        # Min depth
        min_cov_w = self._create_widget(
            AssemblyConfig.model_fields["min_coverage"], "min_coverage",
            "assembly.min_coverage", asm.min_coverage,
        )
        self._field_widgets["assembly.min_coverage"] = min_cov_w
        form.addRow("Min depth for consensus:", min_cov_w)

        # Min allele frequency — fixed dropdown 0.0–1.0 step 0.1, no stretch
        af_combo = QComboBox()
        af_combo.setSizePolicy(
            af_combo.sizePolicy().horizontalPolicy(),
            af_combo.sizePolicy().verticalPolicy(),
        )
        af_combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
        for i in range(11):
            v = round(i * 0.1, 1)
            af_combo.addItem(str(v), v)
        cur_af = round(float(asm.min_af), 1)
        for i in range(af_combo.count()):
            if af_combo.itemData(i) == cur_af:
                af_combo.setCurrentIndex(i)
                break
        self._field_widgets["assembly.min_af"] = af_combo
        af_row = QHBoxLayout()
        af_row.setContentsMargins(0, 0, 0, 0)
        af_row.addWidget(af_combo, 0)          # stretch=0
        af_container = QWidget()
        af_container.setLayout(af_row)
        form.addRow("Min AF (consensus):", af_container)

        # Reporting AF — may be lower than min_af; sub-dominant alleles are
        # listed in the VCF but never written into the consensus FASTA.
        raf_combo = QComboBox()
        raf_combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
        for i in range(11):
            v = round(i * 0.1, 1)
            raf_combo.addItem(str(v), v)
        cur_raf = round(float(asm.report_min_af), 1)
        for i in range(raf_combo.count()):
            if raf_combo.itemData(i) == cur_raf:
                raf_combo.setCurrentIndex(i)
                break
        self._field_widgets["assembly.report_min_af"] = raf_combo
        raf_row = QHBoxLayout()
        raf_row.setContentsMargins(0, 0, 0, 0)
        raf_row.addWidget(raf_combo, 0)
        raf_container = QWidget()
        raf_container.setLayout(raf_row)
        form.addRow("Min AF (variant report):", raf_container)

        # Min allele count
        min_ac_w = self._create_widget(
            AssemblyConfig.model_fields["min_ac"], "min_ac",
            "assembly.min_ac", asm.min_ac,
        )
        self._field_widgets["assembly.min_ac"] = min_ac_w
        form.addRow("Min allele count:", min_ac_w)

        section.set_content_layout(form)
        return section

    def _build_report_section(self, rep: ReportConfig) -> _CollapsibleSection:
        """Build the Reporting configuration section."""
        section = _CollapsibleSection("Reporting")
        form = QFormLayout()
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)

        # Min query coverage
        qcov_w = QDoubleSpinBox()
        qcov_w.setRange(0.0, 100.0)
        qcov_w.setDecimals(1)
        qcov_w.setSingleStep(5.0)
        qcov_w.setValue(rep.min_qcov)
        self._field_widgets["report.min_qcov"] = qcov_w
        form.addRow("Min query coverage (%):", qcov_w)

        # Min bitscore
        bits_w = QDoubleSpinBox()
        bits_w.setRange(0.0, 10000.0)
        bits_w.setDecimals(1)
        bits_w.setSingleStep(10.0)
        bits_w.setValue(rep.min_bitscore)
        self._field_widgets["report.min_bitscore"] = bits_w
        form.addRow("Min bitscore:", bits_w)

        # Separator
        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        sep.setStyleSheet("color: #ddd;")
        form.addRow(sep)

        # Export Excel
        excel_w = QCheckBox()
        excel_w.setChecked(rep.export_excel)
        self._field_widgets["report.export_excel"] = excel_w
        form.addRow("Export tables as Excel (.xlsx):", excel_w)

        # Sort reads
        sort_w = QCheckBox()
        sort_w.setChecked(rep.sort_reads)
        self._field_widgets["report.sort_reads"] = sort_w
        form.addRow("Sort reads by taxonomy:", sort_w)

        # Top taxonomy rank (dropdown, enabled only when sort_reads is checked)
        rank_combo = QComboBox()
        _RANK_LABELS = {
            "auto":     "Auto (detect from data)",
            "division": "Division (NCBI)",
            "domain":   "Domain (SILVA)",
            "kingdom":  "Kingdom",
            "order":    "Order",
            "family":   "Family",
        }
        for r in _TOPTAXRANK_OPTIONS:
            rank_combo.addItem(_RANK_LABELS.get(r, r.capitalize()), r)
        cur_rank = rep.toptaxrank if rep.toptaxrank in _TOPTAXRANK_OPTIONS else "auto"
        rank_combo.setCurrentIndex(rank_combo.findData(cur_rank))
        rank_combo.setEnabled(rep.sort_reads)
        self._field_widgets["report.toptaxrank"] = rank_combo
        form.addRow("Top taxonomy rank:", rank_combo)

        # Link sort_reads checkbox → enable/disable rank combo
        sort_w.toggled.connect(rank_combo.setEnabled)

        section.set_content_layout(form)
        return section

    def _build_hostfilter_section(
        self,
        bgfilter: BgFilterConfig,
        config: EMtraceConfig,
    ) -> _CollapsibleSection:
        """Host filtering: host genome list and minimap2 thresholds."""
        section = _CollapsibleSection("Background Filtering")
        form = QFormLayout()
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)

        mm2_keys = [k for k, e in config.host_databases.items() if e.search.lower() in ("minimap2", "minimap")]
        selected_set = {v.strip() for v in bgfilter.bgfilter.split(",") if v.strip()}

        self._host_picker_mm2 = _HostGenomePicker(
            host_db_keys=mm2_keys,
            selected_csv=",".join(k for k in mm2_keys if k in selected_set),
        )
        form.addRow("Minimap2 bgfilters:", self._host_picker_mm2)


        # Minimap2 sub-section
        form.addRow(_make_section_label("Minimap2 options:"), QWidget())
        for key, field_info, cur_val in [
            ("bgfilter.minimap.c",              BgFilterMinimapConfig.model_fields["c"],              bgfilter.minimap.c),
            ("bgfilter.minimap.m",              BgFilterMinimapConfig.model_fields["m"],              bgfilter.minimap.m),
            ("bgfilter.minimap.s",              BgFilterMinimapConfig.model_fields["s"],              bgfilter.minimap.s),
            ("bgfilter.minimap.save_mmi_index", BgFilterMinimapConfig.model_fields["save_mmi_index"], bgfilter.minimap.save_mmi_index),
        ]:
            field_name = key.split(".")[-1]
            w = self._create_widget(field_info, field_name, key, cur_val)
            self._field_widgets[key] = w
            lbl = _FIELD_LABELS.get(key) or field_name.replace("_", " ").title()
            form.addRow("  " + lbl + ":", w)


        section.set_content_layout(form)
        return section

    def _build_annotation_section(
        self,
        ann: AnnotationConfig,
        strategies: Dict[str, str],
        databases: Dict,           # Dict[str, DatabaseEntry]
        installed_dbs: set,
    ) -> _CollapsibleSection:
        """Annotation: strategy selector + search tool sub-sections."""
        section = _CollapsibleSection("Annotation")
        container = QWidget()
        vbox = QVBoxLayout(container)
        vbox.setContentsMargins(12, 0, 0, 0)
        vbox.setSpacing(4)

        # Strategy / Round1 / Round2 bar
        # TODO: re-enable when 2-round annotation strategy GUI is complete
        self._strategy_bar = _AnnotationStrategyBar(
            strategies=strategies,
            databases=databases,
            installed_dbs=installed_dbs,
            current_anns=ann.anns,
            ann1=ann.ann1,
            ann2=ann.ann2,
        )
        self._strategy_bar.preset_selected.connect(self._on_preset_selected)
        vbox.addWidget(self._strategy_bar)

        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        sep.setStyleSheet("color: #eee;")
        vbox.addWidget(sep)

        # Subsampling row: checkbox + dropdown
        # BLASTN
        blastn_sec = _CollapsibleSection("BLASTN")
        blastn_form = QFormLayout()
        blastn_form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        bl = ann.blastn
        for field_name in ("bitscore", "evalue", "max_target_seqs", "save_hits"):
            if field_name in _SKIP_FIELDS.get("blastn", set()):
                continue
            key = f"blastn.{field_name}"
            w = self._create_widget(
                BlastnConfig.model_fields[field_name], field_name, key, getattr(bl, field_name)
            )
            self._field_widgets[key] = w
            lbl = _FIELD_LABELS.get(key) or field_name.replace("_", " ").title()
            blastn_form.addRow(lbl + ":", w)
        blastn_sec.set_content_layout(blastn_form)
        vbox.addWidget(blastn_sec)

        # Minimap2
        minimap_sec = _CollapsibleSection("Minimap2")
        mm_form = QFormLayout()
        mm_form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        mm = ann.minimap

        # preset — dropdown
        preset_combo = QComboBox()
        for opt in _MINIMAP_PRESETS:
            preset_combo.addItem(opt)
        idx = preset_combo.findText(mm.preset)
        if idx >= 0:
            preset_combo.setCurrentIndex(idx)
        self._field_widgets["minimap.preset"] = preset_combo
        mm_form.addRow("Preset:", preset_combo)

        # k / w — sensitivity override (0 = use preset default)
        mm_form.addRow(_make_section_label("Sensitivity override (0 = preset default):"), QWidget())
        for field_name in ("k", "w"):
            key = f"minimap.{field_name}"
            w_widget = self._create_widget(
                MinimapConfig.model_fields[field_name], field_name, key, getattr(mm, field_name)
            )
            self._field_widgets[key] = w_widget
            mm_form.addRow(_FIELD_LABELS[key] + ":", w_widget)

        # secondary — dropdown yes/no
        sec_combo = QComboBox()
        for opt in _MINIMAP_SECONDARY:
            sec_combo.addItem(opt)
        idx = sec_combo.findText(mm.secondary)
        if idx >= 0:
            sec_combo.setCurrentIndex(idx)
        self._field_widgets["minimap.secondary"] = sec_combo
        mm_form.addRow("Report secondary alignments:", sec_combo)

        mm_form.addRow(_make_section_label("Scoring:"), QWidget())
        for field_name in ("c", "m", "s", "max_target_seqs", "save_hits"):
            if field_name in _SKIP_FIELDS.get("minimap", set()):
                continue
            key = f"minimap.{field_name}"
            w = self._create_widget(
                MinimapConfig.model_fields[field_name], field_name, key, getattr(mm, field_name)
            )
            self._field_widgets[key] = w
            lbl = _FIELD_LABELS.get(key) or field_name.replace("_", " ").title()
            mm_form.addRow(lbl + ":", w)

        minimap_sec.set_content_layout(mm_form)
        vbox.addWidget(minimap_sec)

        # Diamond BLASTP
        diamond_sec = _CollapsibleSection("Diamond BLASTP")
        dm_form = QFormLayout()
        dm_form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        dm = ann.diamond_blastp

        # mode — dropdown
        mode_combo = QComboBox()
        for opt in _DIAMOND_MODES:
            mode_combo.addItem(opt)
        idx = mode_combo.findText(dm.mode)
        if idx >= 0:
            mode_combo.setCurrentIndex(idx)
        self._field_widgets["diamond_blastp.mode"] = mode_combo
        dm_form.addRow("Sensitivity mode:", mode_combo)

        for field_name in ("bitscore", "evalue", "max_target_seqs", "save_hits"):
            if field_name in _SKIP_FIELDS.get("diamond_blastp", set()):
                continue
            key = f"diamond_blastp.{field_name}"
            w = self._create_widget(
                DiamondBlastpConfig.model_fields[field_name], field_name, key,
                getattr(dm, field_name),
            )
            self._field_widgets[key] = w
            lbl = _FIELD_LABELS.get(key) or field_name.replace("_", " ").title()
            dm_form.addRow(lbl + ":", w)
        diamond_sec.set_content_layout(dm_form)
        vbox.addWidget(diamond_sec)

        # Orfipy
        orfipy_sec = _CollapsibleSection("Orfipy (ORF detection)")
        orf_form = QFormLayout()
        orf_form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        orf = ann.orfipy

        # min_len, then max_len (0 = no limit) immediately after
        for field_name in ("min_len",):
            key = f"orfipy.{field_name}"
            w = self._create_widget(
                OrfipyConfig.model_fields[field_name], field_name, key, getattr(orf, field_name)
            )
            self._field_widgets[key] = w
            lbl = _FIELD_LABELS.get(key) or field_name.replace("_", " ").title()
            orf_form.addRow(lbl + ":", w)

        # max_len — QSpinBox, 0 means no limit
        orf_maxlen_sb = _make_maxlen_spinbox(orf.max_len)
        self._field_widgets["orfipy.max_len"] = orf_maxlen_sb
        orf_form.addRow(_FIELD_LABELS["orfipy.max_len"] + ":", orf_maxlen_sb)

        for field_name in ("table", "partial_5", "partial_3"):
            key = f"orfipy.{field_name}"
            w = self._create_widget(
                OrfipyConfig.model_fields[field_name], field_name, key, getattr(orf, field_name)
            )
            self._field_widgets[key] = w
            lbl = _FIELD_LABELS.get(key) or field_name.replace("_", " ").title()
            orf_form.addRow(lbl + ":", w)

        # strand — dropdown
        strand_combo = QComboBox()
        for opt in _ORFIPY_STRANDS:
            strand_combo.addItem(opt)
        idx = strand_combo.findText(orf.strand)
        if idx >= 0:
            strand_combo.setCurrentIndex(idx)
        self._field_widgets["orfipy.strand"] = strand_combo
        orf_form.addRow("Strand:", strand_combo)

        orfipy_sec.set_content_layout(orf_form)
        vbox.addWidget(orfipy_sec)

        section.set_content_widget(container)
        return section

    # ==================================================================
    # Widget factory
    # ==================================================================

    def _create_widget(
        self,
        field_info: FieldInfo,
        field_name: str,
        key: str,
        current_value: Any,
    ) -> QWidget:
        """Create the appropriate Qt widget for a Pydantic field."""
        annotation = field_info.annotation
        real_type = _unwrap_optional(annotation)

        # Enum → QComboBox
        if isinstance(real_type, type) and issubclass(real_type, Enum):
            combo = QComboBox()
            for member in real_type:
                combo.addItem(member.value, member.value)
            current_text = (
                current_value.value if isinstance(current_value, Enum)
                else str(current_value) if current_value is not None else ""
            )
            idx = combo.findText(current_text)
            if idx >= 0:
                combo.setCurrentIndex(idx)
            return combo

        # bool → QCheckBox
        if real_type is bool:
            cb = QCheckBox()
            cb.setChecked(bool(current_value))
            return cb

        # int → QSpinBox
        if real_type is int:
            sb = QSpinBox()
            sb.setRange(0, 999_999)
            sb.setValue(int(current_value) if current_value is not None else 0)
            return sb

        # float → QLineEdit (handles e-values like 1e-10 + infinity-like values)
        if real_type is float:
            le = QLineEdit()
            le.setText(str(current_value) if current_value is not None else "")
            le.setPlaceholderText("e.g. 1e-5")
            le.setProperty("plain_float", True)
            return le

        # Path → dir/file picker
        if real_type is Path or (isinstance(real_type, type) and issubclass(real_type, Path)):
            return self._create_path_picker(current_value, dir_mode=key in _DIR_PATH_FIELDS)

        # Default: str → QLineEdit
        le = QLineEdit()
        le.setText(str(current_value) if current_value is not None else "")
        return le

    def _create_path_picker(self, current_value: Any, dir_mode: bool = False) -> QWidget:
        container = QWidget()
        layout = QHBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)
        le = QLineEdit()
        le.setText(str(current_value) if current_value is not None else "")
        le.setObjectName("_path_edit")
        btn = QPushButton("\u2026")
        btn.setFixedWidth(28)
        btn.clicked.connect(lambda: self._browse_path(le, dir_mode=dir_mode))
        layout.addWidget(le, stretch=1)
        layout.addWidget(btn)
        return container

    def _browse_path(self, line_edit: QLineEdit, dir_mode: bool = False) -> None:
        if dir_mode:
            path = QFileDialog.getExistingDirectory(self, "Select Directory", line_edit.text())
        else:
            path, _ = QFileDialog.getOpenFileName(self, "Select File", line_edit.text())
        if path:
            line_edit.setText(path)

    # ==================================================================
    # Value readers
    # ==================================================================

    @staticmethod
    def _read_widget_value(widget: QWidget) -> Any:
        if isinstance(widget, QCheckBox):
            return widget.isChecked()
        if isinstance(widget, QSpinBox):
            return widget.value()
        if isinstance(widget, QDoubleSpinBox):
            return widget.value()
        if isinstance(widget, QComboBox):
            data = widget.currentData()
            return data if data is not None else widget.currentText()
        if isinstance(widget, QLineEdit):
            text = widget.text().strip()
            if widget.property("plain_float"):
                try:
                    return float(text)
                except ValueError:
                    return None
            return text
        # Path picker container
        le = widget.findChild(QLineEdit, "_path_edit")
        if le is not None:
            return le.text().strip() or None
        return None

    # ==================================================================
    # Helpers
    # ==================================================================

    def _clear(self) -> None:
        self._field_widgets.clear()
        self._host_picker = None
        self._host_picker_mm2 = None
        self._strategy_bar = None
        self._mode_combo = None
        while self._vbox.count():
            item = self._vbox.takeAt(0)
            w = item.widget()
            if w and w is not self._empty_label:
                w.deleteLater()
            elif item.layout():
                _clear_layout(item.layout())


# ===========================================================================
# Module-level helpers
# ===========================================================================

_MODE_LABELS = {"Virome": "virome", "16S rRNA": "16S"}

# Short names for the settings an analysis-mode switch changes.
_MODE_KEY_LABELS = {
    "report.tail": "tail trim",
    "annotation.anns": "preset",
    "annotation.minimap.max_target_seqs": "minimap2 -N",
    "annotation.minimap.s": "-s",
    "annotation.minimap.m": "-m",
    "em.temp": "EM temp",
    "preprocess.bgfilter.bgfilter": "host filter",
    "preprocess.porechop_abi.no_split": "read splitting",
    "preprocess.porechop_abi.min_split_read_size": "min split piece",
    "preprocess.chopper.minlength": "min length",
    "preprocess.chopper.maxlength": "max length",
}


def _fmt_mode_value(value: Any, key: str = "") -> str:
    if isinstance(value, bool):
        if key.endswith(".no_split"):
            value = not value      # shown as "read splitting"
        return "on" if value else "off"
    if value == "0":
        return "none"
    if value == float("inf"):
        return "no limit"
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def _clear_layout(layout) -> None:
    while layout.count():
        item = layout.takeAt(0)
        w = item.widget()
        if w:
            w.deleteLater()


def _unwrap_optional(annotation: Any) -> Any:
    """Unwrap Optional[X] / Union[X, None] → X."""
    # Python 3.10+ X | None
    if isinstance(annotation, _types.UnionType):
        args = [a for a in get_args(annotation) if a is not type(None)]
        if len(args) == 1:
            return args[0]
    origin = get_origin(annotation)
    if origin is not None and "Union" in str(origin):
        args = [a for a in get_args(annotation) if a is not type(None)]
        if len(args) == 1:
            return args[0]
    return annotation
