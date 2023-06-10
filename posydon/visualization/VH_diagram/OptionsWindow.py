"""Options  Window handler for the VH diagram."""


__authors__ = [
    "Maxime Rambosson <Maxime.Rambosson@etu.unige.ch>",
]


from PyQt5.QtWidgets import (QApplication, QWidget, QVBoxLayout,
                             QCheckBox, QComboBox)
from PyQt5.QtCore import Qt, pyqtSignal

from .PresenterMode import PresenterMode


class OptionsWindow(QWidget):
    """Window to display available options, emit signal for selected option."""

    detail_requiered = pyqtSignal()
    reduce_required = pyqtSignal()
    simplify_required = pyqtSignal()
    diagram_required = pyqtSignal()

    distance_representation_required = pyqtSignal(bool)

    _mode_to_text = {
        PresenterMode.DETAILED: "Detailled view",
        PresenterMode.REDUCED: "Reduced view",
        PresenterMode.SIMPLIFIED: "Simplified view",
        PresenterMode.DIAGRAM: "Diagram view",
    }

    def __init__(self):
        """Initialize a OptionsWindow instance."""
        super(OptionsWindow, self).__init__()

        self.setWindowTitle("Visualizer Options")

        screen_geometry = QApplication.desktop().screenGeometry()
        self.resize(int(screen_geometry.width() * 0.20),
                    int(screen_geometry.height() * 0.20))

        self._dropdow_options_signals = {
            self._mode_to_text[PresenterMode.DETAILED]:
                self._detailled_callback,
            self._mode_to_text[PresenterMode.REDUCED]:
                self._reduced_callback,
            self._mode_to_text[PresenterMode.SIMPLIFIED]:
                self._simplified_callback,
            self._mode_to_text[PresenterMode.DIAGRAM]:
                self._diagram_callback,
        }

        self._main_layout = QVBoxLayout()
        self.setLayout(self._main_layout)

        self._setup_options()

    def _prepare_additional_option(self):
        """Prepare option showed only in special situation."""
        self._distance_checkbox = QCheckBox("Distance representation")
        self._distance_checkbox.stateChanged.connect(
            lambda state: self.distance_representation_required.emit(state)
        )

        self._main_layout.addWidget(self._distance_checkbox)

    def _show_additional_option(self, state: bool):
        self._distance_checkbox.setVisible(state)

    def _setup_options(self):
        """Prepare the different options there callback."""
        self._dropdow_options = QComboBox()
        self._main_layout.addWidget(self._dropdow_options)

        for option in self._dropdow_options_signals.keys():
            self._dropdow_options.addItem(option)

        self._dropdow_options.textActivated.connect(
            lambda text: self._dropdow_options_signals[text]()
        )

        self._prepare_additional_option()
        self._show_additional_option(False)

    def _diagram_callback(self, *, emit=True):
        self._show_additional_option(True)

        if emit:
            self.diagram_required.emit()

    def _simplified_callback(self, *, emit=True):
        self._show_additional_option(False)

        if emit:
            self.simplify_required.emit()

    def _reduced_callback(self, *, emit=True):
        self._show_additional_option(False)

        if emit:
            self.reduce_required.emit()

    def _detailled_callback(self, *, emit=True):
        self._show_additional_option(False)

        if emit:
            self.detail_requiered.emit()

    def set_distance_representation(self, state: bool):
        """Set the distance representation."""
        if state:
            self._distance_checkbox.setCheckState(Qt.Checked)
        else:
            self._distance_checkbox.setCheckState(Qt.Unchecked)

    def set_showed_mode(self, mode: PresenterMode):
        """Set the showed mode."""
        self._dropdow_options.setCurrentText(self._mode_to_text[mode])

        # Active callback without emitting
        self._dropdow_options_signals[self._mode_to_text[mode]](emit=False)
