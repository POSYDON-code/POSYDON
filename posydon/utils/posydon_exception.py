class POSYDON_Exception(Exception):
    """
    POSYDON Exception class that includes subclasses for POSYDON-specific errors.
    
    """
    def __init__(self, error_checking=False, message=""):
        self.error_checking = error_checking
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return self.message
    
class GridError(POSYDON_Exception):    
    def __init__(self, err_message):
        super().__init__(message=err_message)

class FlowError(POSYDON_Exception):    
    def __init__(self, err_message):
        super().__init__(message=err_message)

