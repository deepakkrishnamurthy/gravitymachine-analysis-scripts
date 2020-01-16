import sys
from PyQt5.QtWidgets import QApplication, QWidget, QInputDialog, QLineEdit, QFileDialog

fileDialog = QFileDialog.getExistingDirectory(None, "Open dataset folder")

print(fileDialog)