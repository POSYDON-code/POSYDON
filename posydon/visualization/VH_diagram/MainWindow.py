"""Main Window handler for VH diagram."""


__authors__ = [
    "Maxime Rambosson <Maxime.Rambosson@etu.unige.ch>",
]

try:
    from PyQt5.QtCore import pyqtSignal
    from PyQt5.QtWidgets import QApplication, QMainWindow, QScrollArea
except ImportError:
    raise ImportError('PyQt5 is not installed. Please run `pip install .[vis]` in the POSYDON base directory')

from .GraphVisualizer import GraphVisualizer
from .OptionsWindow import OptionsWindow


# Indirection layer to combine options and widget
# Options's signals maybe replaced by Obvsevable pattern
class VisualizerInterface:
    """Handle for all callback triggered by user."""

    def __init__(self, main_win, visualizer, options):
        """Initialize a VisualizerInterface instance.

        Parameters
        ----------
        main_win : MainWindow
            The main window used.
        visualizer : GraphVisualizer
            Custom widget to display graph.
        options : OptionsWindow
            Window used to display options.

        """
        self._main_win = main_win
        self._visualizer = visualizer
        self._options = options

    def __call__(self):
        """Return the visualizer."""
        return self._visualizer

    def detail_required(self):
        return self._options.detail_requiered

    def reduce_required(self):
        return self._options.reduce_required

    def simplify_required(self):
        return self._options.simplify_required

    def diagram_required(self):
        return self._options.diagram_required

    def save_required(self):
        return self._main_win.save_requiered

    def distance_representation_required(self):
        return self._options.distance_representation_required

    def options(self):
        return self._options


class MainWindow(QMainWindow):
    """Custom window to display the application."""

    save_requiered = pyqtSignal()

    def __init__(self):
        """Initialize a MainWindow instance."""
        super(MainWindow, self).__init__()

        self.setWindowTitle("Simulation visualisation")

        screen_geometry = QApplication.desktop().screenGeometry()
        self.resize(int(screen_geometry.width() * 0.66),
                    int(screen_geometry.height() * 0.66))

        self._visualizer = GraphVisualizer()
        self._option_window = OptionsWindow()

    def start_visualization(self):
        """Star the visualization.

        Returns
        -------
        VisualizerInterface
            Return the interface with current window, binded visualizer widget
            and option window. Expose callbacks.

        """
        self.scroll_area = QScrollArea()
        self.scroll_area.setWidget(self._visualizer)
        self.scroll_area.setWidgetResizable(True)
        self.setCentralWidget(self.scroll_area)
        self._setup_visualizer_option()
        return VisualizerInterface(self, self._visualizer, self._option_window)

    def _setup_visualizer_option(self):
        """Prepare all displayed option on MenuBar."""
        bar = self.menuBar()
        bar.setNativeMenuBar(False)
        option_action = bar.addAction(" &Options")
        option_action.triggered.connect(self._show_option_window)
        option_action = bar.addAction(" &Save")
        option_action.triggered.connect(lambda: self.save_requiered.emit())
        self.setMenuBar(bar)

    def _show_option_window(self, _):
        self._option_window.showNormal()

    def closeEvent(self, event):    # Override closeEvent of QWidget
        """Close the option window."""
        self._option_window.close()
