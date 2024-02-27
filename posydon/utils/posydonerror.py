class POSYDONError(Exception):
    """
    POSYDON Exception class that includes subclasses for POSYDON-specific errors.

    """
    def __init__(self,message=""):
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return self.message

class GridError(POSYDONError):
    def __init__(self, err_message):
        super().__init__(message=err_message)

class FlowError(POSYDONError):
    def __init__(self, err_message):
        super().__init__(message=err_message)
