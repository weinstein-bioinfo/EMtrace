"""Abundance bar chart widget — embedded matplotlib stacked bar.

Reads ``abund_table.tsv`` from the sample reports directory and
renders a stacked bar chart showing taxonomic composition.  Falls
back gracefully if matplotlib is not available.

In real-time mode the chart can be refreshed as new EM abundance
estimates arrive.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import pandas as pd
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QComboBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

logger = logging.getLogger(__name__)

try:
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.figure import Figure

    _HAS_MPL = True
except ImportError:  # pragma: no cover
    _HAS_MPL = False


# Canonical rank order (coarsest → finest); first match wins for defaults.
_RANK_CHAIN = ["division", "domain", "phylum", "class", "order", "family", "genus", "species", "organism"]

# Auto-mode parameters
_AUTO_CASCADE = ["organism", "species", "genus", "family"]   # finest → coarsest probe order
_AUTO_THRESHOLD = 10                # max entries that still pick a finer rank
_AUTO_TOPN_STEPS = (5, 10, 15, 20, 25, 30)   # round actual count up to one of these
_AUTO_FALLBACK_RANK = "genus"
_AUTO_FALLBACK_TOPN = 20

# Long-label handling
_LABEL_TRUNCATE_LEN = 50
_LABEL_LONG_THRESHOLD = 30          # any label > this drops fontsize one notch

# Neutral colour for "Other" bucket
_OTHER_COLOR = "#bbbbbb"

# ── Default colour palette (categorical-20, colourblind-safe subset) ──

_COLORS = [
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
    "#aec7e8",
    "#ffbb78",
    "#98df8a",
    "#ff9896",
    "#c5b0d5",
    "#c49c94",
    "#f7b6d2",
    "#c7c7c7",
    "#dbdb8d",
    "#9edae5",
]


class AbundancePlotWidget(QWidget):
    """Embedded stacked bar chart of taxonomic abundance.

    Reads ``abund_table.tsv`` from the sample's reports directory.
    The chart groups taxa by the selected rank (defaults to species)
    and shows the ``top_n`` most abundant taxa, with the rest collapsed
    into an "Other" category.
    """

    def __init__(self, parent=None):
        super().__init__(parent)

        # ── Controls ─────────────────────────────────────────────────
        self._top_n_combo = QComboBox()
        for n in (5, 10, 15, 20, 25, 30):
            self._top_n_combo.addItem(str(n), n)
        self._top_n_combo.setCurrentIndex(1)  # default 10
        self._top_n_combo.setMinimumWidth(60)
        self._top_n_combo.setToolTip("Number of top taxa to display")

        self._rank_combo = QComboBox()
        self._rank_combo.setMinimumWidth(90)
        self._rank_combo.setToolTip("Taxonomic rank to display")

        self._unit_combo = QComboBox()
        self._unit_combo.addItems(["Reads", "Percent"])
        self._unit_combo.setMinimumWidth(80)
        self._unit_combo.setToolTip("Display read counts or percentage of total")

        self._refresh_btn = QPushButton("Refresh")
        self._refresh_btn.setFixedWidth(80)

        self._status = QLabel("")
        self._status.setStyleSheet("color: #444; font-size: 11px;")

        ctrl_row = QHBoxLayout()
        ctrl_row.setContentsMargins(4, 2, 4, 2)
        ctrl_row.addWidget(QLabel("Top:"))
        ctrl_row.addWidget(self._top_n_combo)
        ctrl_row.addWidget(QLabel("Rank:"))
        ctrl_row.addWidget(self._rank_combo)
        # self._unit_combo intentionally hidden; keep for future use
        ctrl_row.addWidget(self._status, stretch=1)
        ctrl_row.addWidget(self._refresh_btn)

        # ── Canvas / empty state ─────────────────────────────────────
        self._empty_label = QLabel(
            "No abundance data found.\n"
            "Run the EM abundance and report steps first."
            if _HAS_MPL
            else "matplotlib is not installed.\n"
            "Install matplotlib to view abundance plots."
        )
        self._empty_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._empty_label.setStyleSheet("color: #999; font-size: 14px;")

        self._canvas: Optional[FigureCanvas] = None
        self._fig: Optional[Figure] = None
        if _HAS_MPL:
            self._fig = Figure(figsize=(8, 5), dpi=100)
            self._canvas = FigureCanvas(self._fig)

        # ── Layout ───────────────────────────────────────────────────
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addLayout(ctrl_row)
        if self._canvas is not None:
            layout.addWidget(self._canvas, stretch=1)
        layout.addWidget(self._empty_label, stretch=1)

        self._set_has_data(False)

        self._output_dir: Optional[Path] = None
        self._abund: Optional[pd.DataFrame] = None
        self._bar_patches: list = []   # list of (Rectangle, metadata_dict)

        # ── Connections ──────────────────────────────────────────────
        self._refresh_btn.clicked.connect(self._on_refresh)
        self._top_n_combo.currentIndexChanged.connect(self._on_replot)
        self._rank_combo.currentIndexChanged.connect(self._on_replot)
        self._unit_combo.currentIndexChanged.connect(self._on_replot)
        if self._canvas is not None:
            self._canvas.mpl_connect("motion_notify_event", self._on_hover)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def load_sample(self, output_dir: Path) -> None:
        """Set the sample output directory and attempt to load abundance data."""
        self._output_dir = output_dir
        self._on_refresh()

    def refresh(self) -> None:
        """Re-read data from disk and replot."""
        self._on_refresh()

    # ------------------------------------------------------------------
    # Private
    # ------------------------------------------------------------------

    def _set_has_data(self, has: bool) -> None:
        if self._canvas is not None:
            self._canvas.setVisible(has)
        self._top_n_combo.setVisible(has)
        self._rank_combo.setVisible(has)
        self._status.setVisible(has)
        self._empty_label.setVisible(not has)

    def _find_abund_table(self) -> Optional[Path]:
        """Locate abund_table.tsv under the output directory."""
        if self._output_dir is None:
            return None
        # Standard location
        candidate = self._output_dir / "results" / "abund_table.tsv"
        if candidate.exists():
            return candidate
        # Fallback: search recursively
        for p in self._output_dir.rglob("abund_table.tsv"):
            return p
        return None

    def _populate_rank_combo(self, abund: pd.DataFrame) -> None:
        """Fill rank combo with all ranks present in the data.

        "Auto" is always the first item and the default selection on (re)load.
        """
        previous = self._rank_combo.currentData() if self._rank_combo.count() else None

        available = [r for r in _RANK_CHAIN if r in abund.columns]

        self._rank_combo.blockSignals(True)
        self._rank_combo.clear()
        self._rank_combo.addItem("Auto", "auto")
        for r in available:
            self._rank_combo.addItem(r, r)

        # Restore previous selection if still available, else default to Auto.
        idx = self._rank_combo.findData(previous) if previous is not None else -1
        self._rank_combo.setCurrentIndex(idx if idx >= 0 else 0)
        self._rank_combo.blockSignals(False)

    @staticmethod
    def _auto_pick(abund: pd.DataFrame) -> tuple[str, int]:
        """Return ``(rank, top_n)`` chosen via the cascade rule.

        Cascade: finest rank in ``_AUTO_CASCADE`` whose distinct count is
        within ``_AUTO_THRESHOLD``.  Falls back to genus / top-20.
        ``top_n`` is rounded up to the next nice step.
        """
        for rank in _AUTO_CASCADE:
            if rank not in abund.columns:
                continue
            n = abund[rank].astype(str).str.strip().replace("", pd.NA).nunique()
            if 0 < n <= _AUTO_THRESHOLD:
                top_n = next((s for s in _AUTO_TOPN_STEPS if s >= n), _AUTO_TOPN_STEPS[-1])
                return rank, top_n
        return _AUTO_FALLBACK_RANK, _AUTO_FALLBACK_TOPN

    def _on_refresh(self) -> None:
        """Re-read data from disk, repopulate rank combo, and replot."""
        if not _HAS_MPL or self._fig is None:
            return

        path = self._find_abund_table()
        if path is None or not path.exists():
            self._abund = None
            self._set_has_data(False)
            return

        try:
            abund = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        except Exception:
            self._abund = None
            self._set_has_data(False)
            return

        if abund.empty:
            self._abund = None
            self._set_has_data(False)
            return

        for col in ("readn", "prob"):
            if col in abund.columns:
                abund[col] = pd.to_numeric(abund[col], errors="coerce").fillna(0)

        self._abund = abund
        self._populate_rank_combo(abund)
        self._on_replot()

    def _on_replot(self) -> None:
        """Replot using cached data and current control values."""
        if self._abund is None or not _HAS_MPL or self._fig is None:
            return
        rank = self._rank_combo.currentData() or "auto"
        top_n = int(self._top_n_combo.currentData())
        rank_is_auto = (rank == "auto")

        if rank_is_auto:
            label_col, _ = self._auto_pick(self._abund)
        else:
            label_col = rank

        if label_col not in self._abund.columns:
            return
        value_col = "readn" if "readn" in self._abund.columns else self._abund.columns[0]
        self._plot_bar(self._abund, value_col, label_col, top_n, rank_is_auto)

    def _plot_bar(
        self,
        abund: pd.DataFrame,
        value_col: str,
        label_col: str,
        top_n: int,
        is_auto: bool = False,
    ) -> None:
        """Render the stacked bar chart."""
        show_pct = self._unit_combo.currentText() == "Percent"
        has_prob = "prob" in abund.columns and abund["prob"].gt(0).any()

        # Aggregate readn + weighted-mean prob per taxon label
        if has_prob:
            tmp = abund.copy()
            tmp["_pw"] = tmp[value_col] * tmp["prob"]
            _agg = tmp.groupby(label_col, dropna=False).agg(
                _val=(value_col, "sum"), _pw=("_pw", "sum")
            )
            _agg["_prob"] = _agg["_pw"] / _agg["_val"].replace(0, 1)
            grouped = _agg["_val"].sort_values(ascending=False)
            prob_by_label = _agg["_prob"]
        else:
            grouped = (
                abund.groupby(label_col, dropna=False)[value_col]
                .sum()
                .sort_values(ascending=False)
            )
            prob_by_label = None

        # Top N + Other — but skip "Other" if it would only collapse 1 entry.
        if len(grouped) > top_n + 1:
            top = grouped.iloc[:top_n]
            other_val = grouped.iloc[top_n:].sum()
            top = pd.concat([top, pd.Series({"Other": other_val})])
            if prob_by_label is not None:
                rest = grouped.iloc[top_n:]
                rest_p = prob_by_label.reindex(rest.index, fill_value=0)
                other_prob = (rest * rest_p).sum() / rest.sum() if rest.sum() > 0 else 0.0
                prob_by_label = pd.concat([
                    prob_by_label.reindex(grouped.iloc[:top_n].index, fill_value=0),
                    pd.Series({"Other": other_prob}),
                ])
        else:
            top = grouped

        # Filter zeros
        top = top[top > 0]
        if top.empty:
            self._set_has_data(False)
            return

        # Convert to percentage of the whole table (not just top N)
        if show_pct:
            total = abund[value_col].sum()
            top = top / total * 100.0 if total > 0 else top * 0.0

        self._set_has_data(True)
        self._bar_patches = []
        self._status.setText("")

        self._fig.clear()
        ax = self._fig.add_subplot(111)

        full_labels = [str(l) if str(l) else "Unknown" for l in top.index]
        values = top.values

        # Neutral grey for "Other"; palette colour for everything else.
        colors = [
            _OTHER_COLOR if lbl == "Other" else _COLORS[i % len(_COLORS)]
            for i, lbl in enumerate(full_labels)
        ]

        # Truncate display labels but keep full names for hover tooltips.
        display_labels = [
            l if len(l) <= _LABEL_TRUNCATE_LEN else l[: _LABEL_TRUNCATE_LEN - 1] + "…"
            for l in full_labels
        ]
        fontsize = 7 if any(len(l) > _LABEL_LONG_THRESHOLD for l in full_labels) else 8

        # Keep raw counts for tooltip regardless of display unit.
        raw_grouped = abund.groupby(label_col, dropna=False)[value_col].sum()
        grand_total = abund[value_col].sum() or 1.0

        containers = ax.barh(range(len(full_labels)), values, color=colors, edgecolor="white", linewidth=0.5)
        x_max = float(values.max()) if len(values) > 0 else 1.0
        for rect, label in zip(containers, full_labels):
            raw_count = int(round(raw_grouped.get(label, 0)))
            pct = 100.0 * raw_count / grand_total
            prob_val = float(prob_by_label.get(label, 0)) if prob_by_label is not None else None
            self._bar_patches.append((rect, {
                "rank": label_col,
                "name": label,
                "count": raw_count,
                "pct": pct,
                "prob": prob_val,
            }))
            # Inline p= label for bars wide enough to hold text (>7% of x range)
            if prob_val and x_max > 0 and rect.get_width() / x_max >= 0.07:
                ax.text(
                    rect.get_width() / 2, rect.get_y() + rect.get_height() / 2,
                    f"p={prob_val:.2f}",
                    ha="center", va="center",
                    fontsize=6, color="white", fontweight="bold",
                )

        ax.set_yticks(range(len(full_labels)))
        ax.set_yticklabels(display_labels, fontsize=fontsize)
        ax.invert_yaxis()
        # Reserve y-space for `top_n` rows so a single bar doesn't fill the plot.
        n_slots = max(top_n, len(full_labels))
        ax.set_ylim(n_slots - 0.5, -0.5)
        ax.set_xlabel("%" if show_pct else "Read count")
        if is_auto:
            ax.set_title(f"Top {top_n} taxa by {label_col} (auto)")
        else:
            ax.set_title(f"Top {top_n} taxa by {label_col}")

        self._fig.tight_layout()
        self._canvas.draw()

    def _on_hover(self, event) -> None:
        if event.inaxes is None:
            self._status.setText("")
            return
        for rect, md in self._bar_patches:
            if rect.contains(event)[0]:
                prob_txt = (
                    f", P(EM)={md['prob']:.3f}"
                    if md.get("prob") is not None and md["prob"] > 0
                    else ""
                )
                self._status.setText(
                    f"{md['rank']}: {md['name']} — "
                    f"{md['count']:,} reads ({md['pct']:.1f}%){prob_txt}"
                )
                return
        self._status.setText("")
