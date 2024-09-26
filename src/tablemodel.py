from Qt.QtCore import QAbstractTableModel, Qt, QModelIndex

class TableModel(QAbstractTableModel):
    """A model to interface a Qt view with pandas dataframe """

    def __init__(self, sqd_cluster_data, sqd_data, parent=None):
        QAbstractTableModel.__init__(self, parent)
        self._sqd_data = sqd_data
        self._sqd_cluster_data = sqd_cluster_data

        self._header = ["Id", "Mol Id", "Hits",
                        "Density", "Overlap", "Correlation", "Cam", "Inside"]

        # mapping of columns (from view to data)
        # self._mapping = [-1, -1, 10, 11, 12, 13]

    def rowCount(self, parent=QModelIndex()) -> int:
        """ Override method from QAbstractTableModel

        Return row count of the pandas DataFrame
        """
        if parent == QModelIndex():
            return len(self._sqd_cluster_data)

        return 0

    def columnCount(self, parent=QModelIndex()) -> int:
        """Override method from QAbstractTableModel

        Return column count of the pandas DataFrame
        """
        if parent == QModelIndex():
            if len(self._sqd_cluster_data) == 0:
               return 0
            else:
               return len(self._header)

    def data(self, index: QModelIndex, role=Qt.ItemDataRole):
        """Override method from QAbstractTableModel

        Return data cell from the pandas DataFrame
        """
        if not index.isValid():
            return None

        if role == Qt.DisplayRole:
            column = index.column()

            mol_idx = int(self._sqd_cluster_data[index.row(), 0])
            record_idx = int(self._sqd_cluster_data[index.row(), 1])
            iter_idx = int(self._sqd_cluster_data[index.row(), 2])
            
            if column == 0:
                return int(index.row() + 1)
            elif column == 1:
                return mol_idx
            elif column == 2:
                return int(self._sqd_cluster_data[index.row(), 3])
            elif 3 <= column <= 7:
                record_row = self._sqd_data[mol_idx, record_idx, iter_idx]
                return float(round(float(record_row[index.column() + 4]) * 10000)) / 10000.0  # for 4 decimals
                # return float(record_row[index.column() + 4])

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