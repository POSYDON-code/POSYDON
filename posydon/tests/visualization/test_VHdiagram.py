import unittest
from unittest.mock import patch
from posydon.utils.common_functions import PATH_TO_POSYDON
import os

from posydon.visualization.VH_diagram.Presenter import Presenter, PresenterMode
from PyQt5.QtWidgets import QApplication
from PyQt5.QtCore import QTimer

PATH_TO_DATASET = os.path.join(
    PATH_TO_POSYDON,
    "posydon",
    "tests",
    "data",
    "POSYDON-UNIT-TESTS",
    "visualization",
    "20000_binaries.csv.gz"
)

# https://stackoverflow.com/questions/60692711/cant-create-python-qapplication-in-github-action

# if not os.path.exists(PATH_TO_DATASET):
#     raise ValueError("Dataset for unit testing (VH diagram) was not found!")
#
#
# class TestVHdiagram(unittest.TestCase):
#     def test_termination_detailled_view(self):
#         app = QApplication.instance()
#         if not app:
#             app = QApplication([])
#
#         presenter = Presenter(PATH_TO_DATASET)
#
#         with patch('PyQt5.QtWidgets.QMainWindow.show') as show_patch:
#             presenter.present(19628, PresenterMode.DETAILED)
#
#             QTimer.singleShot(0, lambda : presenter.close() )
#
#             app.exec_()
#
#             assert show_patch.called
#
#     def test_termination_reduced_view(self):
#         app = QApplication.instance()
#         if not app:
#             app = QApplication([])
#
#         presenter = Presenter(PATH_TO_DATASET)
#
#         with patch('PyQt5.QtWidgets.QMainWindow.show') as show_patch:
#             presenter.present(19628, PresenterMode.REDUCED)
#
#             QTimer.singleShot(0, lambda : presenter.close() )
#
#             app.exec_()
#
#             assert show_patch.called
#
#     def test_termination_simplified_view(self):
#         app = QApplication.instance()
#         if not app:
#             app = QApplication([])
#
#         presenter = Presenter(PATH_TO_DATASET)
#
#         with patch('PyQt5.QtWidgets.QMainWindow.show') as show_patch:
#             presenter.present(19628, PresenterMode.SIMPLIFIED)
#
#             QTimer.singleShot(0, lambda : presenter.close() )
#
#             app.exec_()
#
#             assert show_patch.called
#
#     def test_termination_diagram_view(self):
#         app = QApplication.instance()
#         if not app:
#             app = QApplication([])
#
#         presenter = Presenter(PATH_TO_DATASET)
#
#         with patch('PyQt5.QtWidgets.QMainWindow.show') as show_patch:
#             presenter.present(19628, PresenterMode.DIAGRAM)
#
#             QTimer.singleShot(0, lambda : presenter.close() )
#
#             app.exec_()
#
#             assert show_patch.called
#
# if __name__ == "__main__":
#     unittest.main()
