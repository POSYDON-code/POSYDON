from posydon.binary_evol.binarystar import BinaryStar
from posydon.binary_evol.singlestar import SingleStar

class POSYDONError(Exception):
    """
    POSYDON Exception class that includes subclasses for POSYDON-specific errors.

    """
    def __init__(self, message="", objects=None):
        """General POSYDON exception.

        Parameters
        ----------
        message : str
            POSYDONError message.
        objects : None or list of objects.
            If not None, a list of accompanied objects.

        """
        self.message = message
        self.objects = objects
        super().__init__(self.message)

    def __str__(self):
        result = self.message
        if self.objects is not None:
            for i, object in enumerate(self.objects):
                if isinstance(object, (BinaryStar, SingleStar)):
                    result += (f"\nOBJECT #{i+1} ({type(object)})\n"
                               + str(object))
        return self.message

class GridError(POSYDONError):
    def __init__(self, err_message):
        super().__init__(message=err_message)

class FlowError(POSYDONError):
    def __init__(self, err_message):
        super().__init__(message=err_message)
