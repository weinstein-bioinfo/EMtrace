"""Figure Viewer widget — file selector dropdown + zoomable QGraphicsView preview.

Scans the sample output directory for PNG, SVG, and PDF files, populates a
dropdown for file selection, and shows a zoomable preview below.
"""

from __future__ import annotations

from pathlib import Path

from PyQt6.QtCore import QRectF, QSize, QTimer, QUrl, Qt
from PyQt6.QtGui import QDesktopServices, QPainter, QPixmap
from PyQt6.QtWidgets import (
    QComboBox,
    QFileDialog,
    QGraphicsScene,
    QGraphicsView,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QStackedWidget,
    QVBoxLayout,
    QWidget,
)

try:
    from PyQt6.QtWebEngineWidgets import QWebEngineView
    _WEBENGINE_AVAILABLE = True
except ImportError:
    _WEBENGINE_AVAILABLE = False

_FIGURE_PATTERNS = ("*.png", "*.svg", "*.pdf", "*.html")


class FigureViewerWidget(QWidget):
    """Interactive figure viewer with file selector dropdown and zoomable preview.

    Supports PNG, SVG, and PDF formats.  Multi-page PDFs are rendered page by
    page stacked vertically; Fit uses width-fit for those and contain-fit for
    single figures.  PDF pages are rasterised via ``PyQt6.QtPdf`` when
    available; otherwise the file is opened in the system default viewer.
    """

    def __init__(self, parent=None):
        super().__init__(parent)

        # ── File selector ────────────────────────────────────────────
        self._file_combo = QComboBox()
        self._file_combo.setSizeAdjustPolicy(
            QComboBox.SizeAdjustPolicy.AdjustToContents
        )
        self._file_combo.setMinimumWidth(200)

        # ── Controls ─────────────────────────────────────────────────
        self._zoom_in_btn  = QPushButton("+")
        self._zoom_out_btn = QPushButton("−")
        self._orig_btn     = QPushButton("1:1")
        self._fit_btn      = QPushButton("Fit")
        self._export_btn   = QPushButton("Save PNG…")

        for btn, w in [
            (self._zoom_in_btn,  36),
            (self._zoom_out_btn, 36),
            (self._orig_btn,     40),
            (self._fit_btn,      40),
            (self._export_btn,   96),
        ]:
            btn.setFixedWidth(w)

        _zf = self._zoom_in_btn.font()
        _zf.setPointSize(_zf.pointSize() + 2)
        for _b in (self._zoom_in_btn, self._zoom_out_btn):
            _b.setFont(_zf)

        top_row = QHBoxLayout()
        top_row.setContentsMargins(4, 2, 4, 2)
        top_row.addWidget(QLabel("File:"))
        top_row.addWidget(self._file_combo, stretch=1)
        top_row.addWidget(self._zoom_in_btn)
        top_row.addWidget(self._zoom_out_btn)
        top_row.addWidget(self._orig_btn)
        top_row.addWidget(self._fit_btn)
        top_row.addWidget(self._export_btn)

        # ── Raster/PDF preview ───────────────────────────────────────
        self._scene = QGraphicsScene()
        self._view  = QGraphicsView(self._scene)
        self._view.setDragMode(QGraphicsView.DragMode.ScrollHandDrag)
        self._view.setRenderHints(
            QPainter.RenderHint.Antialiasing
            | QPainter.RenderHint.SmoothPixmapTransform
            | QPainter.RenderHint.LosslessImageRendering
        )

        # ── HTML / interactive viewer ────────────────────────────────
        if _WEBENGINE_AVAILABLE:
            self._web_view: QWebEngineView | None = QWebEngineView()
        else:
            self._web_view = None

        # ── Stacked widget — page 0: raster/PDF, page 1: HTML ────────
        self._stack = QStackedWidget()
        self._stack.addWidget(self._view)
        if self._web_view is not None:
            self._stack.addWidget(self._web_view)

        # ── Empty-state label ────────────────────────────────────────
        self._empty_label = QLabel(
            "No figures found.\n"
            "Run pipeline steps to generate plots."
        )
        self._empty_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._empty_label.setStyleSheet("color: #999; font-size: 14px;")

        # ── Layout ───────────────────────────────────────────────────
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addLayout(top_row)
        layout.addWidget(self._stack, stretch=1)
        layout.addWidget(self._empty_label, stretch=1)

        self._is_multipage: bool = False  # True for PDFs with >1 page

        self._set_has_figures(False)

        # ── Connections ──────────────────────────────────────────────
        self._file_combo.currentIndexChanged.connect(self._on_file_changed)
        self._zoom_in_btn.clicked.connect(lambda: self._view.scale(1.25, 1.25))
        self._zoom_out_btn.clicked.connect(lambda: self._view.scale(0.8, 0.8))
        self._orig_btn.clicked.connect(self._view.resetTransform)
        self._fit_btn.clicked.connect(self._fit_view)
        self._export_btn.clicked.connect(self._on_export)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def load_sample_figures(self, output_dir: Path) -> None:
        """Scan sample output directory for figure files and populate the combo."""
        self._file_combo.blockSignals(True)
        self._file_combo.clear()
        self._scene.clear()

        paths: list[Path] = []
        for pat in _FIGURE_PATTERNS:
            paths.extend(sorted(output_dir.rglob(pat)))

        # Deduplicate while preserving order
        seen: set[Path] = set()
        for p in paths:
            if p not in seen:
                seen.add(p)
                rel = p.relative_to(output_dir)
                self._file_combo.addItem(str(rel), str(p))

        self._file_combo.blockSignals(False)

        has_files = self._file_combo.count() > 0
        self._set_has_figures(has_files)
        if has_files:
            self._file_combo.setCurrentIndex(0)
            self._on_file_changed(0)

    def load_figure(self, path: Path) -> None:
        """Select and display a specific figure file.

        If the path is in the combo, selects it; otherwise loads it directly.
        """
        idx = self._file_combo.findData(str(path))
        if idx >= 0:
            self._file_combo.setCurrentIndex(idx)
        else:
            self._show_figure(path)

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _set_has_figures(self, has: bool) -> None:
        self._file_combo.setVisible(has)
        self._stack.setVisible(has)
        self._empty_label.setVisible(not has)
        self._set_raster_controls_visible(has)

    def _set_raster_controls_visible(self, visible: bool) -> None:
        """Show/hide controls that only make sense for raster/PDF files."""
        for btn in (
            self._zoom_in_btn, self._zoom_out_btn,
            self._orig_btn, self._fit_btn, self._export_btn,
        ):
            btn.setVisible(visible)

    def _on_file_changed(self, index: int) -> None:
        if index < 0:
            return
        path = Path(self._file_combo.itemData(index))
        self._show_figure(path)

    def _show_figure(self, path: Path) -> None:
        suffix = path.suffix.lower()
        if suffix == ".html":
            self._load_html(path)
            return

        self._is_multipage = False
        self._scene.clear()
        self._stack.setCurrentIndex(0)
        self._set_raster_controls_visible(True)

        if suffix == ".svg":
            self._load_svg(path)
        elif suffix == ".pdf":
            self._load_pdf(path)
        else:
            px = QPixmap(str(path))
            if not px.isNull():
                item = self._scene.addPixmap(px)
                item.setTransformationMode(
                    Qt.TransformationMode.SmoothTransformation
                )
            else:
                self._scene.addText(f"Cannot load: {path.name}")

        self._fit_view()

    def _load_html(self, path: Path) -> None:
        """Display an interactive HTML file (e.g. Plotly) in the web viewer."""
        if self._web_view is not None:
            self._web_view.load(QUrl.fromLocalFile(str(path)))
            self._stack.setCurrentIndex(1)
        else:
            # QWebEngineView unavailable — open in system browser
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(path)))
            self._scene.clear()
            self._scene.addText(
                f"Opened in browser: {path.name}\n"
                "(install PyQt6-WebEngine to view inside the app)"
            )
            self._stack.setCurrentIndex(0)
        self._set_raster_controls_visible(False)

    def _load_svg(self, path: Path) -> None:
        try:
            from PyQt6.QtSvgWidgets import QGraphicsSvgItem

            svg_item = QGraphicsSvgItem(str(path))
            self._scene.addItem(svg_item)
        except Exception as exc:  # noqa: BLE001
            self._scene.addText(f"SVG load failed: {exc}")

    def _load_pdf(self, path: Path) -> None:
        """Render all PDF pages via QPdfDocument stacked vertically.

        Each page is rasterised at ~150 DPI.  Falls back to the OS default
        viewer if QtPdf is unavailable.
        """
        _TARGET_W_PX = 1200  # raster width (~150 dpi at A4)
        _PAGE_GAP    = 12    # vertical gap between pages (px)

        rendered = False
        try:
            from PyQt6.QtPdf import QPdfDocument  # type: ignore[import]

            doc = QPdfDocument(self)
            doc.load(str(path))
            n_pages  = doc.pageCount()
            y_offset = 0.0

            for page_idx in range(n_pages):
                pts = doc.pagePointSize(page_idx)
                if pts.isEmpty():
                    continue
                scale    = _TARGET_W_PX / pts.width()
                img_size = QSize(int(pts.width() * scale), int(pts.height() * scale))
                img      = doc.render(page_idx, img_size)
                if img.isNull():
                    continue
                px   = QPixmap.fromImage(img)
                item = self._scene.addPixmap(px)
                item.setPos(0.0, y_offset)
                item.setTransformationMode(Qt.TransformationMode.SmoothTransformation)
                y_offset += px.height() + _PAGE_GAP
                rendered = True

            self._is_multipage = n_pages > 1
            doc.load("")   # release file handle
        except Exception:  # noqa: BLE001
            pass

        if not rendered:
            from PyQt6.QtCore import QUrl
            from PyQt6.QtGui import QDesktopServices

            QDesktopServices.openUrl(QUrl.fromLocalFile(str(path)))
            self._scene.addText(f"Opened externally: {path.name}")

    def showEvent(self, event) -> None:
        super().showEvent(event)
        if self._stack.currentIndex() == 0:
            QTimer.singleShot(0, self._fit_view)

    def _fit_view(self) -> None:
        if self._stack.currentIndex() != 0:
            return
        # itemsBoundingRect() recomputes from actual items; sceneRect() can be
        # stale after clear()+add, so it sometimes returns the wrong rect.
        if not self._scene.items():
            return
        rect = self._scene.itemsBoundingRect()
        if rect.isEmpty():
            return
        self._view.setSceneRect(rect)
        if self._is_multipage:
            # Multi-page PDF: fit to width so pages are readable; user scrolls down.
            vp_w = self._view.viewport().width()
            if rect.width() > 0 and vp_w > 0:
                self._view.resetTransform()
                self._view.scale(vp_w / rect.width(), vp_w / rect.width())
                self._view.horizontalScrollBar().setValue(0)
                self._view.verticalScrollBar().setValue(0)
        else:
            self._view.fitInView(rect, Qt.AspectRatioMode.KeepAspectRatio)

    def _on_export(self) -> None:
        """Export the current scene to a PNG file chosen by the user."""
        rect = self._scene.sceneRect()
        if rect.isEmpty():
            return
        save_path, _ = QFileDialog.getSaveFileName(
            self, "Export Figure", "", "PNG Images (*.png)"
        )
        if not save_path:
            return
        pixmap = QPixmap(int(rect.width()), int(rect.height()))
        pixmap.fill(Qt.GlobalColor.white)
        painter = QPainter(pixmap)
        self._scene.render(painter, source=QRectF(rect))
        painter.end()
        pixmap.save(save_path, "PNG")
