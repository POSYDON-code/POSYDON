"""The definition of the VHdiagram class."""

__authors__ = [
    "Maxime Rambosson <Maxime.Rambosson@etu.unige.ch>",
]


__credits__ = [
    "Simone Bavera <Simone.Bavera@unige.ch>"
]


from enum import Enum, auto

import matplotlib.pyplot as plt
from IPython.display import Image, display
from PyQt5.QtCore import QTimer
from PyQt5.QtWidgets import QApplication

from .VH_diagram.Presenter import Presenter, PresenterMode


class DisplayMode(Enum):
    """The different view for the Presenter."""

    WINDOW = auto()
    INLINE_S = auto()
    INLINE_B = auto()


class VHdiagram:
    """Handle a VH diagram."""

    def __init__(
        self,
        filename,
        path="./",
        index=0,
        *,
        presentMode=PresenterMode.DETAILED,
        displayMode=DisplayMode.WINDOW,
        figsize=(10, 8)
    ):
        """Initialize a VHdiagram instance."""
        self._app = (
            QApplication.instance()
        )  # Check if there is instance of QApplication
        if not self._app:  # if not, create it
            self._app = QApplication([])

        self._presenter = Presenter(filename=filename, path=path)

        self._presenter.present(index, presentMode)

        if displayMode == DisplayMode.INLINE_B:
            QTimer.singleShot(0, lambda: self._display_inline_b())
        elif displayMode == DisplayMode.INLINE_S:
            QTimer.singleShot(0, lambda: self._display_inline_s(figsize))

        self._app.exec_()

    def _display_inline_b(self):
        filepath = self._presenter.screen()
        self._presenter.close()
        display(Image(filename=filepath))

    def _display_inline_s(self, figsize):
        filepath = self._presenter.screen()
        self._presenter.close()

        image = plt.imread(filepath)

        plt.figure(figsize=figsize)
        plt.imshow(image, aspect="auto")
        plt.axis("off")
