"""The POSYDON exception class and subclasses for more specific errors."""


from posydon.binary_evol.binarystar import BinaryStar
from posydon.binary_evol.singlestar import SingleStar


class POSYDONError(Exception):
    """General POSYDON exception class."""

    def __init__(self, message="", objects=None):
        """General POSYDON exception.

        Parameters
        ----------
        message : str
            POSYDONError message.
        objects : None or list of objects.
            A list of accompanied objects (or None if not set), that will be
            handled differently when `str` method is used.

        """
        self.message = message
        self.objects = objects
        super().__init__(self.message)

    def __str__(self):
        """Create the text that accompanies this exception."""
        result = self.message
        if self.objects is not None:
            for i, obj in enumerate(self.objects):
                if isinstance(obj, (BinaryStar, SingleStar)):
                    result += f"\nOBJECT #{i+1} ({type(obj)}):\n{str(obj)}"
        return self.message


class GridError(POSYDONError):
    """POSYDON error specific for PSyGrid operations."""


class FlowError(POSYDONError):
    """POSYDON error specific for binary evolution flow errors."""
