"""Interactive sunburst widget — embedded matplotlib hierarchical chart.

Reads ``abund_table.tsv`` and draws concentric rings of taxonomic
abundance (division → class → order → family → genus → species).
Hover a wedge for read count + percentage; click to drill into a branch.

Lightweight in-app companion to the Plotly ``sunburst.html`` export —
no QWebEngineView (Chromium) dependency required.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from emtrace.steps.plots import _RANK_CHAIN, _RANK_COLORS

logger = logging.getLogger(__name__)

try:
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.figure import Figure
    from matplotlib.patches import Wedge
    import matplotlib.colors as mcolors

    _HAS_MPL = True
except ImportError:  # pragma: no cover
    _HAS_MPL = False

_RING_WIDTH = 1.0


class SunburstPlotWidget(QWidget):
    """Embedded interactive sunburst of taxonomic abundance.

    Hover wedges for tooltips; left-click to drill down into a branch.
    """

    def __init__(self, parent=None):
        super().__init__(parent)

        # ── Controls ─────────────────────────────────────────────────
        self._back_btn = QPushButton("← Back")
        self._reset_btn = QPushButton("Reset")
        for b in (self._back_btn, self._reset_btn):
            b.setFixedWidth(80)

        self._status = QLabel("")
        self._status.setStyleSheet("color: #444; font-size: 11px;")

        ctrl_row = QHBoxLayout()
        ctrl_row.setContentsMargins(4, 2, 4, 2)
        ctrl_row.addWidget(QLabel("Sunburst:"))
        ctrl_row.addWidget(self._status, stretch=1)
        ctrl_row.addWidget(self._back_btn)
        ctrl_row.addWidget(self._reset_btn)

        # ── Canvas / empty state ─────────────────────────────────────
        self._empty_label = QLabel(
            "No abundance data found.\n"
            "Run the report step first."
            if _HAS_MPL
            else "matplotlib is not installed.\n"
            "Install matplotlib to view the sunburst."
        )
        self._empty_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._empty_label.setStyleSheet("color: #999; font-size: 14px;")

        self._canvas: Optional[FigureCanvas] = None
        self._fig: Optional[Figure] = None
        if _HAS_MPL:
            self._fig = Figure(figsize=(7, 7), dpi=100)
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
        self._drill_path: List[Tuple[str, str]] = []
        self._wedges: List[Wedge] = []

        # ── Connections ──────────────────────────────────────────────
        self._back_btn.clicked.connect(self._on_back)
        self._reset_btn.clicked.connect(self._on_reset)
        if self._canvas is not None:
            self._canvas.mpl_connect("motion_notify_event", self._on_hover)
            self._canvas.mpl_connect("button_press_event", self._on_click)

    # ------------------------------------------------------------------
    # Public API (mirrors AbundancePlotWidget)
    # ------------------------------------------------------------------

    def load_sample(self, output_dir: Path) -> None:
        self._output_dir = output_dir
        self._drill_path = []
        self._abund = self._read_abund()
        self._redraw()

    def refresh(self) -> None:
        self._abund = self._read_abund()
        self._redraw()

    # ------------------------------------------------------------------
    # Private — data
    # ------------------------------------------------------------------

    def _set_has_data(self, has: bool) -> None:
        if self._canvas is not None:
            self._canvas.setVisible(has)
        self._back_btn.setVisible(has)
        self._reset_btn.setVisible(has)
        self._status.setVisible(has)
        self._empty_label.setVisible(not has)

    def _find_abund_table(self) -> Optional[Path]:
        if self._output_dir is None:
            return None
        candidate = self._output_dir / "results" / "abund_table.tsv"
        if candidate.exists():
            return candidate
        for p in self._output_dir.rglob("abund_table.tsv"):
            return p
        return None

    def _read_abund(self) -> Optional[pd.DataFrame]:
        path = self._find_abund_table()
        if path is None or not path.exists():
            return None
        try:
            df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        except Exception:
            return None
        if df.empty or "readn" not in df.columns:
            return None
        df["readn"] = pd.to_numeric(df["readn"], errors="coerce").fillna(0)
        return df

    # ------------------------------------------------------------------
    # Private — drawing
    # ------------------------------------------------------------------

    def _redraw(self) -> None:
        if not _HAS_MPL or self._fig is None or self._abund is None:
            self._set_has_data(False)
            return

        df = self._filter_by_drill_path(self._abund)
        ranks = [r for r in _RANK_CHAIN if r in df.columns][len(self._drill_path):]
        if not ranks or df.empty or float(df["readn"].sum()) <= 0:
            self._set_has_data(False)
            return

        self._set_has_data(True)
        self._fig.clear()
        ax = self._fig.add_subplot(111)
        self._wedges = []

        tree = self._build_hierarchy(df, ranks)
        # Color each top-level child differently; descendants inherit & lighten
        top_children = sorted(
            tree["children"].items(), key=lambda kv: -kv[1]["count"]
        )
        for i, (name, child) in enumerate(top_children):
            child["_color_seed"] = _RANK_COLORS[i % len(_RANK_COLORS)]

        grand_total = float(tree["count"]) or 1.0
        self._draw_node(ax, tree, ring=0, a0=0.0, extent=360.0, color_seed=None, grand_total=grand_total)

        n_rings = len(ranks)
        ax.set_xlim(-(n_rings + 0.2), n_rings + 0.2)
        ax.set_ylim(-(n_rings + 0.2), n_rings + 0.2)
        ax.set_aspect("equal")
        ax.axis("off")

        title = "Taxonomic Abundance — Sunburst"
        if self._drill_path:
            crumbs = " › ".join(name for _, name in self._drill_path)
            title += f"  [{crumbs}]"
        ax.set_title(title, fontsize=11)
        self._fig.tight_layout()
        self._canvas.draw_idle()
        self._status.setText("")

    def _filter_by_drill_path(self, df: pd.DataFrame) -> pd.DataFrame:
        out = df
        for rank, name in self._drill_path:
            if rank in out.columns:
                out = out[out[rank].astype(str).str.strip() == name]
        return out

    @staticmethod
    def _build_hierarchy(
        df: pd.DataFrame, ranks: List[str],
    ) -> Dict[str, Any]:
        has_organism = "organism" in df.columns
        root: Dict[str, Any] = {"count": 0.0, "children": {}, "rank": None, "organisms": set()}
        for _, row in df.iterrows():
            count = float(row["readn"])
            if count <= 0:
                continue
            organism = str(row.get("organism", "")).strip() if has_organism else ""
            node = root
            node["count"] += count
            if organism:
                node["organisms"].add(organism)
            for rank in ranks:
                key = str(row.get(rank, "")).strip() or "unknown"
                children = node["children"]
                if key not in children:
                    children[key] = {"count": 0.0, "children": {}, "rank": rank, "organisms": set()}
                node = children[key]
                node["count"] += count
                if organism:
                    node["organisms"].add(organism)
        return root

    def _draw_node(
        self,
        ax,
        node: Dict[str, Any],
        ring: int,
        a0: float,
        extent: float,
        color_seed: Optional[str],
        grand_total: float,
    ) -> None:
        children = node.get("children") or {}
        if not children:
            return
        items = sorted(children.items(), key=lambda kv: -kv[1]["count"])
        sibling_total = sum(c["count"] for _, c in items) or 1.0
        a = a0
        for name, child in items:
            child_extent = extent * child["count"] / sibling_total
            seed = child.get("_color_seed") or color_seed or _RANK_COLORS[0]
            face = _shade(seed, ring)
            w = Wedge(
                (0.0, 0.0),
                r=(ring + 1) * _RING_WIDTH,
                width=_RING_WIDTH,
                theta1=a,
                theta2=a + child_extent,
                facecolor=face,
                edgecolor="white",
                linewidth=0.5,
            )
            organisms = child.get("organisms") or set()
            w.metadata = {  # type: ignore[attr-defined]
                "rank": child["rank"],
                "name": name,
                "count": int(round(child["count"])),
                "pct": 100.0 * child["count"] / grand_total,
                "organisms": organisms,
            }
            ax.add_patch(w)
            self._wedges.append(w)
            self._add_label(ax, name, a, a + child_extent, ring, face)
            self._draw_node(ax, child, ring + 1, a, child_extent, seed, grand_total)
            a += child_extent

    @staticmethod
    def _add_label(
        ax,
        name: str,
        theta1: float,
        theta2: float,
        ring: int,
        face_rgb: Tuple[float, float, float],
    ) -> None:
        """Draw a taxonomic label centered in the wedge — Plotly-style.

        Tangent to the arc (text follows the ring), flipped to stay upright,
        auto-contrast color, truncated to fit, hidden when too narrow.
        """
        extent = theta2 - theta1
        if extent < 4.0 or not name:
            return  # too narrow to fit any text

        mid = (theta1 + theta2) / 2.0
        rad = math.radians(mid)
        r_mid = (ring + 0.5) * _RING_WIDTH
        x = r_mid * math.cos(rad)
        y = r_mid * math.sin(rad)

        # Wide wedges (≥180°) or near-center: render horizontally.
        if extent >= 180.0:
            rotation = 0.0
        else:
            # Tangent rotation: text follows the arc.
            # Upper-right (0–90°) and lower-left (180–270°) are already upright;
            # upper-left (90–180°) and lower-right (270–360°) need an extra 180°.
            rotation = mid - 90.0
            if mid > 180.0:
                rotation += 180.0

        # Truncate to roughly the available arc length (in plot units).
        # arc_len ≈ extent_rad * r_mid; ~6 chars per unit at default font.
        arc_len = math.radians(extent) * r_mid
        max_chars = max(3, int(arc_len * 5.5))
        label = name if len(name) <= max_chars else name[: max_chars - 1] + "…"

        # Luminance-based contrast (sRGB approximation).
        r, g, b = face_rgb[:3]
        text_color = "white" if (0.299 * r + 0.587 * g + 0.114 * b) < 0.55 else "black"

        fontsize = 9 if ring == 0 else max(6, 9 - ring)
        ax.text(
            x, y, label,
            ha="center", va="center",
            rotation=rotation, rotation_mode="anchor",
            fontsize=fontsize, color=text_color,
            clip_on=True,
        )

    # ------------------------------------------------------------------
    # Private — events
    # ------------------------------------------------------------------

    def _on_hover(self, event) -> None:
        if event.inaxes is None or self._fig is None:
            self._status.setText("")
            return
        for w in self._wedges:
            contains, _ = w.contains(event)
            if contains:
                md = w.metadata  # type: ignore[attr-defined]
                text = (
                    f"{md['rank']}: {md['name']} — "
                    f"{md['count']:,} reads ({md['pct']:.1f}%)"
                )
                organisms = md.get("organisms") or set()
                # Skip the trivial case where the only organism equals the wedge name.
                meaningful = {o for o in organisms if o and o != md["name"]}
                if len(meaningful) == 1:
                    text += f"  [organism: {next(iter(meaningful))}]"
                elif len(meaningful) > 1:
                    text += f"  [{len(meaningful)} organisms]"
                self._status.setText(text)
                return
        self._status.setText("")

    def _on_click(self, event) -> None:
        if event.inaxes is None or event.button != 1:
            return
        # Click in the centre circle → go back (mirrors Plotly sunburst behaviour).
        if math.hypot(event.xdata, event.ydata) < _RING_WIDTH:
            self._on_back()
            return
        for w in self._wedges:
            contains, _ = w.contains(event)
            if contains:
                md = w.metadata  # type: ignore[attr-defined]
                self._drill_path.append((md["rank"], md["name"]))
                self._redraw()
                return

    def _on_back(self) -> None:
        if self._drill_path:
            self._drill_path.pop()
            self._redraw()

    def _on_reset(self) -> None:
        if self._drill_path:
            self._drill_path = []
            self._redraw()


def _shade(hex_color: str, ring: int) -> Tuple[float, float, float]:
    """Mix with white based on ring depth — outer rings lighter."""
    if not _HAS_MPL:
        return (0.5, 0.5, 0.5)
    rgb = mcolors.to_rgb(hex_color)
    factor = min(0.45, 0.10 * ring)
    return tuple(c + (1.0 - c) * factor for c in rgb)
