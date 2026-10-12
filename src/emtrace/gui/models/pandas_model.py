"""QAbstractTableModel backed by a pandas DataFrame.

Adapted from ViRefSet ``pandas_model.py``.

Column 0 is a virtual ``#`` column showing the original 1-based row
index from the input file.  Clicking that column header restores the
original file order.  All other DataFrame columns are shifted by +1.
"""

from __future__ import annotations

from typing import Any, Optional

import pandas as pd
from PyQt6.QtCore import QAbstractTableModel, QModelIndex, Qt


class PandasTableModel(QAbstractTableModel):
    """Read-only table model wrapping a pandas DataFrame.

    Supports column sorting via ``sort()`` (called by QTableView when
    ``setSortingEnabled(True)``).  All cell values are presented as
    strings through ``Qt.DisplayRole``.
    """

    def __init__(self, df: Optional[pd.DataFrame] = None, parent=None):
        super().__init__(parent)
        self._df: pd.DataFrame = df if df is not None else pd.DataFrame()

    # ------------------------------------------------------------------
    # Required overrides
    # ------------------------------------------------------------------

    def rowCount(self, parent: QModelIndex = QModelIndex()) -> int:  # noqa: N802
        if parent.isValid():
            return 0
        return len(self._df)

    def columnCount(self, parent: QModelIndex = QModelIndex()) -> int:  # noqa: N802
        if parent.isValid():
            return 0
        if self._df.empty and len(self._df.columns) == 0:
            return 0
        return len(self._df.columns) + 1  # +1 for virtual '#' column

    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:
        if not index.isValid():
            return None
        if role == Qt.ItemDataRole.DisplayRole:
            col = index.column()
            if col == 0:
                return str(self._df.index[index.row()] + 1)
            val = self._df.iloc[index.row(), col - 1]
            if pd.isna(val):
                return ""
            return str(val)
        if role == Qt.ItemDataRole.TextAlignmentRole:
            if index.column() == 0:
                return int(
                    Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
                )
            return int(
                Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter
            )
        return None

    def headerData(  # noqa: N802
        self,
        section: int,
        orientation: Qt.Orientation,
        role: int = Qt.ItemDataRole.DisplayRole,
    ) -> Any:
        if role != Qt.ItemDataRole.DisplayRole:
            return None
        if orientation == Qt.Orientation.Horizontal:
            if section == 0:
                return "#"
            return str(self._df.columns[section - 1])
        return str(section + 1)

    # ------------------------------------------------------------------
    # Sorting
    # ------------------------------------------------------------------

    def sort(self, column: int, order: Qt.SortOrder = Qt.SortOrder.AscendingOrder) -> None:
        if self._df.empty:
            return
        self.layoutAboutToBeChanged.emit()
        ascending = order == Qt.SortOrder.AscendingOrder
        if column == 0:
            self._df = self._df.sort_index(ascending=ascending)
        else:
            col = self._df.columns[column - 1]
            self._df = self._df.sort_values(
                col, ascending=ascending, key=_sort_key
            )
        self.layoutChanged.emit()

    # ------------------------------------------------------------------
    # Public helpers
    # ------------------------------------------------------------------

    def set_dataframe(self, df: pd.DataFrame) -> None:
        """Replace the backing DataFrame and refresh the view."""
        self.beginResetModel()
        self._df = df.reset_index(drop=True)
        self.endResetModel()

    @property
    def dataframe(self) -> pd.DataFrame:
        """Access the underlying DataFrame (e.g. for export)."""
        return self._df


def _sort_key(series: pd.Series) -> pd.Series:
    """Try numeric conversion so '100' sorts after '9'."""
    numeric = pd.to_numeric(series, errors="coerce")
    if numeric.notna().sum() > len(series) * 0.5:
        return numeric
    return series.str.lower() if series.dtype == object else series
