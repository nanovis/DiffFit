from Qt.QtCore import QAbstractTableModel, Qt, QModelIndex

class TableModel(QAbstractTableModel):
    """A model to interface a Qt view with pandas dataframe """

    def __init__(self, data, parent=None):
        QAbstractTableModel.__init__(self, parent)
        self._data = data;
        self._header = ["Mol Id", "Quat Id", "Shift Id",
                        "x", "y", "z",
                        "rw", "rx", "ry", "rz",
                        "Density", "Overlap", "Correlation", "Cam"]

    def rowCount(self, parent=QModelIndex()) -> int:
        """ Override method from QAbstractTableModel

        Return row count of the pandas DataFrame
        """
        if parent == QModelIndex():
            return len(self._data)                

        return 0

    def columnCount(self, parent=QModelIndex()) -> int:
        """Override method from QAbstractTableModel

        Return column count of the pandas DataFrame
        """
        if parent == QModelIndex():
            if len(self._data) == 0 or len(self._data[0]) == 0:
               return 0
            else:
               return len(self._data[0][0])                

    def data(self, index: QModelIndex, role=Qt.ItemDataRole):
        """Override method from QAbstractTableModel

        Return data cell from the pandas DataFrame
        """
        if not index.isValid():
            return None

        if role == Qt.DisplayRole:
            column = index.column()
            
            if column < 3:
                return int(self._data[index.row()][0][index.column()])
            else:
                return float(self._data[index.row()][0][index.column()])
                    
        return None

    def headerData(
        self, section: int, orientation: Qt.Orientation, role: Qt.ItemDataRole
    ):
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                if section <= len(self._header):
                    return self._header[section]
                else:
                    return "<unknown>"
                    
            if orientation == Qt.Vertical:
                return ""
        
        return None