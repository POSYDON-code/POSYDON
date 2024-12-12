"""Provide LaTeX-interpreted labels in the VH diagram."""


__authors__ = [
    "Maxime Rambosson <Maxime.Rambosson@etu.unige.ch>",
]

try:
    from PyQt5 import QtGui, QtWidgets
    from PyQt5.QtWidgets import QVBoxLayout
    from PyQt5.QtCore import Qt
except ImportError:
    raise ImportError('PyQt5 is not installed. Please run `pip install .[vis]` in the POSYDON base directory')

import matplotlib as mpl
from matplotlib.backends.backend_agg import FigureCanvasAgg


class MathTextLabel(QtWidgets.QWidget):
    """Custom label to display LaTeX formula."""

    def __init__(self, parent=None, **kwargs):
        """Initialize a MathTextLabel instance."""
        super(QtWidgets.QWidget, self).__init__(parent, **kwargs)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self._label = QtWidgets.QLabel()
        layout.addWidget(self._label)

    def setText(self, text):
        """Set the diplayed text.

        Parameters
        ----------
        text : str
            Text to display, interpreted as LaTeX formula.

        """
        if text == "":
            return

        fig = mpl.figure.Figure()
        fig.patch.set_facecolor('none')
        fig.set_canvas(FigureCanvasAgg(fig))
        renderer = fig.canvas.get_renderer()

        # plot the mathTex expression

        ax = fig.add_axes([0, 0, 1, 1])
        ax.axis('off')
        ax.patch.set_facecolor('none')
        t = ax.text(0, 0, text, ha='left', va='bottom', fontsize=11)

        # fit figure size to text artist

        fwidth, fheight = fig.get_size_inches()
        fig_bbox = fig.get_window_extent(renderer)

        text_bbox = t.get_window_extent(renderer)

        tight_fwidth = text_bbox.width * fwidth / fig_bbox.width
        tight_fheight = text_bbox.height * fheight / fig_bbox.height

        fig.set_size_inches(tight_fwidth, tight_fheight)

        # convert mpl figure to QPixmap

        buf, size = fig.canvas.print_to_buffer()
        qimage = QtGui.QImage.rgbSwapped(QtGui.QImage(
            buf, size[0], size[1], QtGui.QImage.Format_ARGB32))
        qpixmap = QtGui.QPixmap(qimage)

        self._label.setPixmap(qpixmap)
        self._label.setAlignment(Qt.AlignCenter)
