from Qt.QtCore import QAbstractTableModel, Qt, QModelIndex

class TableModel(QAbstractTableModel):
    """A model to interface a Qt view with pandas dataframe """

    def __init__(self, data, parent=None):
        QAbstractTableModel.__init__(self, parent)
        self._data = data;

    def rowCount(self, parent=QModelIndex()) -> int:
        """ Override method from QAbstractTableModel

        Return row count of the pandas DataFrame
        """
        if parent == QModelIndex():
            return len(self._data) - 1                

        return 0

    def columnCount(self, parent=QModelIndex()) -> int:
        """Override method from QAbstractTableModel

        Return column count of the pandas DataFrame
        """
        if parent == QModelIndex():
            if len(self._data) == 0:
               return 0
            else:
               return len(self._data[0])                

    def data(self, index: QModelIndex, role=Qt.ItemDataRole):
        """Override method from QAbstractTableModel

        Return data cell from the pandas DataFrame
        """
        if not index.isValid():
            return None

        if role == Qt.DisplayRole:
            return str(self._data[index.row() + 1][index.column()])
        # str(self._dataframe.iloc[index.row(), index.column()])

        return None

    def headerData(
        self, section: int, orientation: Qt.Orientation, role: Qt.ItemDataRole
    ):
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return str(self._data[0][section])
        # str(self._dataframe.columns[section])

            if orientation == Qt.Vertical:
                return ""
        # str(self._dataframe.index[section])

        return None