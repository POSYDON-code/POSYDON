"""The definition of the PresenterMode class in VH diagram."""


__authors__ = [
    "Maxime Rambosson <Maxime.Rambosson@etu.unige.ch>",
]


from enum import Enum, auto


class PresenterMode(Enum):
    """The different view for the Presenter."""

    DIAGRAM = auto()
    REDUCED = auto()
    SIMPLIFIED = auto()
    DETAILED = auto()
