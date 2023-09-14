"""Graph Visualizer for the VH diagram."""


__authors__ = [
    "Maxime Rambosson <Maxime.Rambosson@etu.unige.ch>",
]


from PyQt5.QtWidgets import (QWidget, QFrame, QHBoxLayout, QVBoxLayout,
                             QLabel, QGridLayout)
from PyQt5.QtCore import Qt, QPoint, QSize
from PyQt5.QtGui import QPainter, QPixmap, QFont
from .MathTextLabel import MathTextLabel

from dataclasses import dataclass
from enum import Enum, auto


class INFOSTYPE(Enum):
    """Enumeration of different type of infos."""

    CASE = auto()
    POINT = auto()
    STATE = auto()


@dataclass
class Infos:
    """Mother-class with common informations uselful for all widget."""

    column_id: int
    infos_type: INFOSTYPE
    connected: bool

    def __init__(self, column_id, infos_type, connected):
        """Initialize an Infos instance."""
        self.column_id = column_id
        self.infos_type = infos_type
        self.connected = connected


@dataclass
class CaseInfos(Infos):
    """Informations to create a case widget."""

    centered_text: str
    bot_right_text: str
    bot_left_text: str
    top_right_text: str
    top_left_text: str

    border_width: int

    def __init__(
        self,
        column_id,
        *,
        centered_txt: str = "",
        bot_right_txt: str = "",
        bot_left_txt: str = "",
        top_right_txt: str = "",
        tot_left_txt: str = "",
        border_width: int = 0,
    ):
        """Initialize a CaseInfos instance."""
        super(CaseInfos, self).__init__(column_id, INFOSTYPE.CASE, True)

        self.centered_text = centered_txt
        self.bot_right_text = bot_right_txt
        self.bot_left_text = bot_left_txt
        self.top_right_text = top_right_txt
        self.top_left_text = tot_left_txt
        self.border_width = border_width


@dataclass
class PointInfos(Infos):
    """Informations to create a widget with a point drew."""

    text: str

    def __init__(self, column_id, text: str = ""):
        """Initialize a PointInfos instance."""
        super(PointInfos, self).__init__(column_id, INFOSTYPE.POINT, True)

        self.text = text


class GraphVisualizerItem(QWidget):
    """Define the mother-class for widget in GraphVisualizer.

    Attributes
    ----------
    connected : bool
        Indicate if this widget need to be connected with the previous one in
        the same column.

    """

    def __init__(self):
        """Initialize a GraphVisualizerItem instrance."""
        super(GraphVisualizerItem, self).__init__()
        self.connected = False

    def get_attach_point_top(self):
        """Get coordinate to connect this widget with the previous one.

        Returns
        -------
        int
            Coordinate to link the widget with the previous one.

        """
        return self.mapToGlobal(self.rect().center())  # Default behavior

    def get_attach_point_bot(self):
        """Get coordinate to connect this widget with the next one.

        Returns
        -------
        int
            Coordinate to link the widget with the next one.

        """
        return self.mapToGlobal(self.rect().center())  # Default behavior


class GraphVisualizerCase(GraphVisualizerItem):
    """Case widget in GraphVisualizer.

    Provides the ability to display 5 texts: 2 on the top 2 on the bottom,
    and 1 in center. Can have a border.
    """

    def __init__(self):
        """Initialize a GraphVisualizerCase instance."""
        super(GraphVisualizerCase, self).__init__()

        self._v_layout = QVBoxLayout()
        self.setLayout(self._v_layout)

        top_label_layout = QHBoxLayout()
        self._v_layout.addLayout(top_label_layout)

        self._top_left_label = MathTextLabel()
        top_label_layout.addWidget(self._top_left_label)
        self._top_right_label = MathTextLabel()
        top_label_layout.addWidget(self._top_right_label)

        self._center_widget = QFrame()
        self._center_text_label = QLabel()
        self._center_text_label.setAlignment(Qt.AlignCenter)
        layout = QHBoxLayout()
        layout.addWidget(self._center_text_label)
        self._center_widget.setLayout(layout)

        self._v_layout.addWidget(self._center_widget)

        bot_label_layout = QHBoxLayout()
        self._v_layout.addLayout(bot_label_layout)

        self._bot_left_label = MathTextLabel()
        bot_label_layout.addWidget(self._bot_left_label)
        self._bot_right_label = MathTextLabel()
        bot_label_layout.addWidget(self._bot_right_label)

    def set_central_text(self, text):
        """Set the text at the center."""
        self._center_text_label.setText(text)

    def set_bottom_left_text(self, text):
        """Set the text at the bottom left."""
        self._bot_left_label.setText(text)

    def set_bottom_rigth_text(self, text):
        """Set the text at the bottom right."""
        self._bot_right_label.setText(text)

    def set_top_left_text(self, text):
        """Set the text at the top left."""
        self._top_left_label.setText(text)

    def set_top_right_text(self, text):
        """Set the text at the top right."""
        self._top_right_label.setText(text)

    def set_central_border(self, width):
        """Set the text at the central border."""
        self._center_widget.setFrameShape(QFrame.Box)
        self._center_widget.setLineWidth(width)

    def get_attach_point_top(self):
        """Get coordinate to connect this widget with the previous one."""
        return self._center_widget.mapToGlobal(
            self._center_widget.rect().topLeft()
        ) + QPoint(int(self._center_widget.width() / 2), 0)

    def get_attach_point_bot(self):
        """Get coordinate to connect this widget with the next one."""
        return self._center_widget.mapToGlobal(
            self._center_widget.rect().bottomLeft()
        ) + QPoint(int(self._center_widget.width() / 2), 0)


def prepare_case(infos: CaseInfos):
    """Help to create GraphVisualizerCase from CaseInfos.

    Parameters
    ----------
    infos : CaseInfos
        Infos to create the widget.

    Returns
    -------
    GraphVisualizerCase
        Created widget with given infos.

    """
    case = GraphVisualizerCase()

    case.connected = infos.connected

    case.set_central_text(infos.centered_text)

    case.set_bottom_left_text(infos.bot_left_text)
    case.set_bottom_rigth_text(infos.bot_right_text)
    case.set_top_left_text(infos.top_left_text)
    case.set_top_right_text(infos.top_right_text)
    if infos.border_width != 0:
        case.set_central_border(infos.border_width)

    return case


class GraphVisualizerPointDraw(QWidget):
    """Define an empty widget with a point drew."""

    def __init__(self):
        """Initialize a GraphVisualizerPointDraw instance."""
        super(GraphVisualizerPointDraw, self).__init__()

        self.setMinimumSize(QSize(13, 13))
        self.setMaximumSize(QSize(13, 13))

    def paintEvent(self, event):  # override paintEvent of QWidget
        """Paint an event."""
        painter = QPainter(self)
        painter.drawEllipse(self.rect().center(), 6, 6)
        painter.setBrush(Qt.black)
        painter.drawEllipse(self.rect().center(), 2, 2)


class GraphVisualizerPoint(GraphVisualizerItem):
    """Widget containing GraphVisualizerPointDraw.

    Provides the ability to display 2 texts, one at each side.
    """

    def __init__(self):
        """Initialize a GraphVisualizerPoint instance."""
        super(GraphVisualizerPoint, self).__init__()
        self._layout = QGridLayout()
        self.setLayout(self._layout)

        self._point_draw = GraphVisualizerPointDraw()
        self._layout.addWidget(self._point_draw, 0, 0)

        self._label = MathTextLabel()
        self._layout.addWidget(self._label, 0, 1)

    def set_text(self, text):
        """Se the text of the label."""
        self._label.setText(text)

    def get_attach_point_top(self):
        """Get coordinate to connect this widget with the previous one."""
        return self._point_draw.mapToGlobal(self._point_draw.rect().center())

    def get_attach_point_bot(self):
        """Get coordinate to connect this widget with the next one."""
        return self._point_draw.mapToGlobal(self._point_draw.rect().center())


def prepare_point(infos: PointInfos):
    """Help to create GraphVisualizerPoint from PointInfos.

    Parameters
    ----------
    infos : PointInfos
        Infos to create the widget.

    Returns
    -------
    GraphVisualizerPoint
        Created widget with given infos.

    """
    point = GraphVisualizerPoint()
    point.connected = infos.connected
    point.set_text(infos.text)
    return point


@dataclass
class StateInfos(Infos):
    """Information to create a widget with a line of diagram inside."""

    S1_filename: str
    S2_filename: str
    event_filename : str
    distance: int
    top_texts: list
    bot_texts: list

    def __init__(
        self,
        column_id,
        *,
        S1_filename=None,
        S2_filename=None,
        event_filename=None,
        distance=1,
    ):
        """Initialize a StateInfos instance."""
        super(StateInfos, self).__init__(column_id, INFOSTYPE.STATE, False)

        self.distance = distance
        self.S1_filename = S1_filename
        self.S2_filename = S2_filename
        self.event_filename = event_filename
        self.top_texts = []
        self.bot_texts = []


class GraphVisualizerState(GraphVisualizerItem):
    """Widget containing drawings.

    Provides the ability to display 4 texts on top & bottom.
    """

    _max_height = 128
    _offset = 10

    def __init__(self):
        """Initialize a GraphVisualizerState instance."""
        super(GraphVisualizerState, self).__init__()

        self.done = False

        layout = QGridLayout()
        layout.setContentsMargins(0, 10, 0, 10)
        self.setLayout(layout)

        self._bg_container = QLabel()
        self._bg_container.resize(self._max_height, self._max_height)
        layout.addWidget(self._bg_container, 1, 0, 1, 4)

        self._bg = QPixmap(self._max_height, self._max_height)
        self._bg_container.setPixmap(self._bg)

        self._S1_pixmap = QPixmap()
        self._S2_pixmap = QPixmap()
        self._event_pixmap = QPixmap()

        self._distance = 1

        self._top_labels = []
        self._bot_labels = []
        for i in range(4):
            self._top_labels.append(MathTextLabel(self))
            layout.addWidget(
                self._top_labels[i], 0, i, alignment=Qt.AlignCenter)

            self._bot_labels.append(MathTextLabel(self))
            layout.addWidget(
                self._bot_labels[i], 2, i, alignment=Qt.AlignCenter)

    def set_first_star(self, filename):
        """Load the picture 'filename' and resize it if needed.

        Parameters
        ----------
        filename : str
            Name (+ path) to the picture.

        Returns
        -------
        bool
            Indicate if loading & resize success.

        """
        if not self._S1_pixmap.load(filename):
            return False

        if self._S1_pixmap.height() > self._max_height:
            self._S1_pixmap = self._S1_pixmap.scaledToHeight(self._max_height)

        return True

    def set_second_star(self, filename):
        """Load the picture 'filename' and resize it if needed.

        Parameters
        ----------
        filename : str
            Name (+ path) to the picture.

        Returns
        -------
        bool
            Indicate if loading & resize success.

        """
        if not self._S2_pixmap.load(filename):
            return False

        if self._S2_pixmap.height() > self._max_height:
            self._S2_pixmap = self._S2_pixmap.scaledToHeight(self._max_height)

        return True

    def set_event(self, filename):
        """Load the picture 'filename' and resize it if needed.

        Parameters
        ----------
        filename : str
            Name (+ path) to the picture.

        Returns
        -------
        bool
            Indicate if loading & resize success.

        """
        if not self._event_pixmap.load(filename):
            return False

        if self._event_pixmap.height() > self._max_height:
            self._event_pixmap = self._event_pixmap.scaledToHeight(
                self._max_height)

        return True

    def set_distance(self, dist):
        """Set the distance."""
        if dist > 1 or dist < 0:
            print("`dist` needs to be normalised.")
        self._distance = dist

    def set_top_text(self, text, index):
        """Set the text of the top label."""
        if index >= len(self._top_labels):
            print(f"Can't set text to top label {index}")
            return

        self._top_labels[index].setText(text)

    def set_bot_text(self, text, index):
        """Se the text of the bottom label."""
        if index >= len(self._bot_labels):
            print(f"Can't set text to bot label {index}")
            return
        
        self._bot_labels[index].setText(text)

    def get_attach_point_top(self):
        """Get coordinate to connect this widget with the previous one."""
        return self._bg_container.mapToGlobal(
            self._bg_container.rect().topLeft()
        ) + QPoint(int(self._bg_container.width() / 2), 0)

    def get_attach_point_bot(self):
        """Get coordinate to connect this widget with the next one."""
        return self._bg_container.mapToGlobal(
            self._bg_container.rect().bottomLeft()
        ) + QPoint(int(self._bg_container.width() / 2), 0)

    def resizeEvent(self, event):
        """Resize the event."""
        self._bg_container.resize(event.size().width(),
                                  self._bg_container.height())
        # Offset is used here to keep the possibility to scale down
        self._bg = QPixmap(self._bg_container.width() - self._offset,
                           self._max_height)

        self._bg_container.setPixmap(self._bg)

    def paintEvent(self, event):  # override paintEvent of QWidget
        """Paint the event."""
        self._bg.fill()
        painter = QPainter(self._bg)

        x_S1 = ((self._bg.width() / 2 - self._S1_pixmap.width())
                * (1 - self._distance))
        y_S1 = self._bg.height() / 2 - self._S1_pixmap.height() / 2

        painter.drawPixmap(int(x_S1), int(y_S1), self._S1_pixmap)

        x_S2 = (
            self._bg.width() / 2
            + (self._bg.width() / 2 - self._S2_pixmap.width()) * self._distance
        )
        y_S2 = self._bg.height() / 2 - self._S2_pixmap.height() / 2

        painter.drawPixmap(int(x_S2), int(y_S2), self._S2_pixmap)

        if not self._event_pixmap.isNull():
            x_event = self._bg.width() / 2 - self._event_pixmap.width() / 2

            painter.drawPixmap(int(x_event), 0, self._event_pixmap)

        self._bg_container.setPixmap(self._bg)


def prepare_state(infos: StateInfos):
    """Help to create GraphVisualizerState from StateInfos.

    Parameters
    ----------
    infos : StateInfos
        Infos to create the widget.

    Returns
    -------
    GraphVisualizerState
        Created widget with given infos.

    """
    state = GraphVisualizerState()

    state.connected = infos.connected

    if infos.S1_filename:
        if not state.set_first_star(infos.S1_filename):
            print(f"Can't load {infos.S1_filename} as S1")

    if infos.S2_filename:
        if not state.set_second_star(infos.S2_filename):
            print(f"Can't load {infos.S2_filename} as S2")

    if infos.event_filename:
        if not state.set_event(infos.event_filename):
            print(f"Can't load {infos.event_filename} as event")

    state.set_distance(infos.distance)

    for i in range(len(infos.top_texts)):
        state.set_top_text(infos.top_texts[i], i)

    for i in range(len(infos.bot_texts)):
        state.set_bot_text(infos.bot_texts[i], i)

    return state


class columnTYPE(Enum):
    """Enumeration of different column type."""

    TIMELINE = auto()
    CONNECTED = auto()


@dataclass
class ConnectedItem:
    """Represent a visual link between 2 widgets."""

    from_item: GraphVisualizerItem
    to_item: GraphVisualizerItem

    def __init__(self, from_item: GraphVisualizerItem,
                 to_item: GraphVisualizerItem):
        """Initialize a ConnectedItem instance."""
        self.from_item = from_item
        self.to_item = to_item


class GraphVisualizercolumn:
    """Mother-class of visual column in GraphVisualizer.

    Manage one visual column in the QGridLayout
    (can take several logical columns).
    """

    def __init__(self, grid, column_id, column_span):
        """Initialize a GraphVisualizercolumn instance.

        Parameters
        ----------
        grid : QGridLayout
            Grid where column take place.
        column_id : int
            Unique id of this column.
        column_span : int
            Number of logical column used.

        """
        self._column_id = column_id
        self._row_index = 0
        self._column_span = column_span

        self._grid = grid
        self._items = []

        self._create_title_label()

        self._connected_items = []
        self._last_item = None

    def set_title(self, title):
        """Set the title."""
        font = QFont()
        font.setPointSize(25)
        self._title_label.setFont(font)
        self._title_label.setText(title)

    def add_item(self, item):
        """Add one item in the column, connect it to previous if needed.

        Parameters
        ----------
        item : GraphVisualizerItem
            Widget to add in column.

        """
        self._add_widget(item)

        self._items.append(item)

        if (
            item.connected and self._last_item is not None
        ):  # It's not the first item added ::
            self._connected_items.append(ConnectedItem(self._last_item, item))

        self._last_item = item

    def clear(self):
        """Delete all widget in the column, keep the title."""
        self._last_item = None
        self._connected_items = []

        for item in self._items:
            item.deleteLater()

        self._items = []
        self._row_index = 1

    def reset(self):
        """Delete everythings (widget + title), reset logical column used."""
        self._last_item = None
        self._connected_items = []

        self._title_label.deleteLater()

        for item in self._items:
            item.deleteLater()

        for i in range(self._column_span):
            self._grid.setColumnStretch(self._column_id + i, 0)

        self._items = []
        self._row_index = 0

    def skip(self):
        """Skip one row."""
        self._row_index += 1

    def get_id(self):
        """Return the column id."""
        return self._column_id

    def _create_title_label(self):
        self._title_label = QLabel()
        self._title_label.setAlignment(Qt.AlignCenter)
        self._add_widget(self._title_label)

    def _add_widget(self, widget):
        """Add widget to the logical column used.

        Parameters
        ----------
        widget : QWidget
            Widget to add.

        """
        self._grid.addWidget(
            widget, self._row_index, self._column_id, 1, self._column_span
        )
        self._row_index += 1


class GraphVisualizerConnectedcolumn(GraphVisualizercolumn):
    """Simple visual column with arrow between connected widget."""

    def __init__(self, grid, column_id, column_span=1):
        """Initialize a GraphVisualizerConnectedcolumn instance."""
        super(GraphVisualizerConnectedcolumn, self).__init__(
            grid, column_id, column_span
        )

        for i in range(column_span):
            self._grid.setColumnStretch(self._column_id + i, 1)

    def draw(self, surface):
        """Draw the surface."""
        painter = QPainter(surface)

        for connection in self._connected_items:
            start = surface.mapFromGlobal(
                connection.from_item.get_attach_point_bot())
            end = surface.mapFromGlobal(
                connection.to_item.get_attach_point_top())

            painter.drawLine(start + QPoint(0, 36)
                             ,end + QPoint(0, -8))

            left = end + QPoint(-4, -8)
            rigth = end + QPoint(4, -8)

            painter.setBrush(Qt.black)
            painter.drawEllipse(start + QPoint(0, 36), 2, 2)
            painter.drawPolygon(end, left, rigth)


class GraphVisualizerTimeline(GraphVisualizercolumn):
    """Draw a visual column, compressed by another other column."""

    def __init__(self, grid, column_id, column_span=1):
        """Initialize a GraphVisualizerTimeline instance."""
        super(GraphVisualizerTimeline, self).__init__(
            grid, column_id, column_span)
        for i in range(column_span):
            self._grid.setColumnStretch(self._column_id + i, 0)

    def draw(self, surface):
        """Draw the surface."""
        painter = QPainter(surface)

        for connection in self._connected_items:
            start = surface.mapFromGlobal(
                connection.from_item.get_attach_point_bot())
            end = surface.mapFromGlobal(
                connection.to_item.get_attach_point_top())

            painter.drawLine(start, end)


class GraphVisualizer(QWidget):
    """Widget used to display the different columns and add widget in them."""

    def __init__(self):
        """Initialize a GraphVisualizer instance."""
        super(GraphVisualizer, self).__init__()

        self._layout = QGridLayout()
        self._layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(self._layout)

        self._next_column = 0

        self._columns = []

    def add_column(self, column_type, column_span=1):
        """Create a column according to the column_type.

        It takes account for the column_span of the logical column.

        Parameters
        ----------
        column_type : columnTYPE
            Type of the column to add.
        column_span : int
            Nb of logical column used.

        Returns
        -------
        int
            ID of created column.

        """
        if column_type == columnTYPE.TIMELINE:
            self._columns.append(
                GraphVisualizerTimeline(
                    self._layout, self._next_column, column_span)
            )
        elif column_type == columnTYPE.CONNECTED:
            self._columns.append(
                GraphVisualizerConnectedcolumn(
                    self._layout, self._next_column, column_span)
            )

        self._next_column += column_span

        return len(self._columns) - 1

    def get_column(self, column_id):
        """Return visual column with the given id.

        Parameters
        ----------
        column_id : int
            ID of searched column.

        Returns
        -------
        GraphVisualizercolumn
            column with the given ID.

        """
        if column_id < len(self._columns):
            return self._columns[column_id]
        return None

    def add_line(self, infos):
        """Add all lines based on all elements in `infos`.

        For each info in infos, create the corresponding widget and add it to
        the column with corresponding id (given in each info), only 1 widget
        per column by call, if one column haven't any widget associated,
        column skip this row.

        Parameters
        ----------
        infos : Array of Infos
            Array with derivated struct of Infos to create the different
            widget for this line.

        """
        for column in self._columns:
            info = next((x for x in infos if (
                lambda info: info.column_id == column.get_id())(x)), None)

            if info is not None:
                item = None
                if info.infos_type == INFOSTYPE.CASE:
                    item = prepare_case(info)
                elif info.infos_type == INFOSTYPE.POINT:
                    item = prepare_point(info)
                elif info.infos_type == INFOSTYPE.STATE:
                    item = prepare_state(info)

                if item is not None:
                    column.add_item(item)
            else:
                column.skip()

    def clear(self):
        """Clear each column."""
        for col in self._columns:
            col.clear()

    def reset(self):
        """Clear all delete each column."""
        for col in self._columns:
            col.reset()

        self._next_column = 0
        self._columns = []

    def paintEvent(self, event):  # override paintEvent of QWidget
        """Paint the event."""
        for col in self._columns:
            col.draw(self)

    def saveAsPicture(self, filename):
        """Render this widget in a picture and save it a filename.

        Parameters
        ----------
        filename : str
            Destination file of the picture.
        """
        pixmap = QPixmap(self.size())
        self.render(pixmap)

        if not pixmap.save(filename):
            print("Can't save as " + filename)
