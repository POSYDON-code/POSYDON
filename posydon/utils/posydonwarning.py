"""The POSYDON warnings class and subclasses for more detailed warning
handling.

Mainly the function `Pwarn` should be called to submit a warning, which behaves
like warnings.warn.
All the defined POSYDON warnings are collected in `_POSYDONWarning_subclasses`
for those will be checked automatically, when using `Pwarn`.
The functions `SetPOSYDONWarnings`, `AllPOSYDONWarnings`, and
`NoPOSYDONWarnings` allow to change the filter settings.
It should be noted, that we register POSYDON warnings by category, filename,
and lineno. We do not include the warning text in the identifier of a new
warning for the default filter, hence it behaves different to python warnings.
But we record how often each warning is issued with your identifier. Those
records can be requested via `get_stats` and printed via `print_stats`.
To catch warnings we have the context manager class `Catch_POSYDON_Warnings`.
It will record the warnings and issue them at the end of the context. To not
issue the recorded warnings, simply use the `reset_cache` function of
`Catch_POSYDON_Warnings` before the end of the context. You can get the list of
the catched warnings via `get_cache`.
Here an example, where the argument True of get_cache, will empty the list of
catched warnings after it is copy to be returned:
    with Catch_POSYDON_Warnings() as cpw:
        ...
        Pwarn("Test", "POSYDONWarning")
        ...
        print(cpw.get_cache(True))
"""


__authors__ = [
    "Monica Gallegos-Garcia <monicagallegosgarcia@u.northwestern.edu>",
    "Camille Liotine <cliotine@u.northwestern.edu>",
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",   
]


import copy
import sys
import warnings


class POSYDONWarning(Warning):
    """General POSYDON warning class."""
    def __init__(self, message=''):
        self.message = message
    
    def __str__(self):
        return repr(self.message)

class ApproximationWarning(POSYDONWarning):
    """Warning that a physical approximation was used during binary evolution."""
    def __init__(self, message=''):
        super().__init__(message)

class BinaryParsingWarning(POSYDONWarning):
    """Warnings related to parsing of binaries."""
    def __init__(self, message=''):
        super().__init__(message)

class ClassificationWarning(POSYDONWarning):
    """Warnings related to classification during binary evolution."""
    def __init__(self, message=''):
        super().__init__(message)

class EvolutionWarning(POSYDONWarning):
    """Warning that something unexpeted occurred during the binary evolution, but
    the evolution is able to continue (binary did not fail)."""
    def __init__(self, message=''):
        super().__init__(message)
        
class InappropriateValueWarning(POSYDONWarning):
    """Warnings that a strange value is used."""
    def __init__(self, message=''):
        super().__init__(message)

class IncompletenessWarning(POSYDONWarning):
    """Warnings when not all tasks could be done."""
    def __init__(self, message=''):
        super().__init__(message)

class InterpolationWarning(POSYDONWarning):
    """Warnings related to interpolation during binary evolution."""
    def __init__(self, message=''):
        super().__init__(message)

class MissingFilesWarning(POSYDONWarning):
    """Warnings related to missing files."""
    def __init__(self, message=''):
        super().__init__(message)

class OverwriteWarning(POSYDONWarning):
    """Warning that a data will get overwritten."""
    def __init__(self, message=''):
        super().__init__(message)

class ReplaceValueWarning(POSYDONWarning):
    """Warning that a value got replaced."""
    def __init__(self, message=''):
        super().__init__(message)

class UnsupportedModelWarning(POSYDONWarning):
    """Warnings related to selecting a model that is not supported."""
    def __init__(self, message=''):
        super().__init__(message)


# all POSYDON warnings subclasses should be defined beforehand
_POSYDONWarning_subclasses = {cls.__name__: cls for cls in\
                              POSYDONWarning.__subclasses__()}

def _get_POSYDONWarning_class(category):
    """Inferring the POSYDONWarning class.

        Parameters
        ----------
        category : str or a POSYDON warnings category class
            A string, which will be converted into a POSYDONWarning class.

        Returns
        -------
        POSYDONWarning or None
            A POSYDONWarning class or None if the input is neither a
            POSYDONWarning class nor convertable to one (all unknown strings
            are converted to the base POSYDONWarning).
    """
    global _POSYDONWarning_subclasses
    if isinstance(category, str):
        if category in _POSYDONWarning_subclasses.keys():
            category = _POSYDONWarning_subclasses[category]
        else:
            category = POSYDONWarning
    if isinstance(category, type) and issubclass(category, POSYDONWarning):
        return category
    else:
        return None


_POSYDON_WARNINGS_REGISTRY = {}

def get_stats():
    """Return the statistics of the warnings."""
    global _POSYDON_WARNINGS_REGISTRY
    return _POSYDON_WARNINGS_REGISTRY

def print_stats():
    """Prints the statistics of the warnings."""
    global _POSYDON_WARNINGS_REGISTRY
    if len(_POSYDON_WARNINGS_REGISTRY)==0:
        print("No warnings occured.")
    else:
        print("There have been warnings:\n", _POSYDON_WARNINGS_REGISTRY)

def _apply_POSYDON_filter(warning=dict(message="No warning")):
    """Filter a warning.

        Parameters
        ----------
        warning : dict (default: {message="No warning"})
            Dictionary containing all options passed to warnings.warn.

        Returns
        -------
        warning or None in case it got filtered out.
    """
    global _POSYDON_WARNINGS_REGISTRY
    if not isinstance(warning, dict):
        raise TypeError("warning must be a dictionary.")
    # get stack level
    stacklevel = warning.get('stacklevel', 0)
    if not isinstance(stacklevel, int):
        raise TypeError("stacklevel must be an integer.")
    # get category
    category = warning.get('category', None)
    if not(isinstance(category, type) and issubclass(category, Warning)):
        category = _get_POSYDONWarning_class(category)
    if not isinstance(category, type):
        category = UserWarning
    # get filename and lineno from frame at stacklevel
    try:
        frame = sys._getframe(stacklevel)
    except:
        g = sys.__dict__
        filename = "sys"
        lineno = 1
    else:
        g = frame.f_globals
        filename = frame.f_code.co_filename
        lineno = frame.f_lineno
    # get module
    if isinstance(g, dict):
        module = g.get('__name__', "posydonwarnings")
    else:
        module = "posydonwarnings"
    # check registry
    if not isinstance(_POSYDON_WARNINGS_REGISTRY, dict):
        print("Reset _POSYDON_WARNINGS_REGISTRY, old was:",
              _POSYDON_WARNINGS_REGISTRY)
        _POSYDON_WARNINGS_REGISTRY = {}
    # set key for registry
    key = (category, filename, lineno)
    # get message text
    text = warning.get('message', "")
    if not isinstance(text, str):
        raise TypeError("message must be an integer.")
    # search the filters
    for item in warnings.filters:
        action, msg, cat, mod, ln = item
        if ((msg is None or msg.match(text)) and
            issubclass(category, cat) and
            (mod is None or mod.match(module)) and
            (ln == 0 or lineno == ln)):
            break
    else:
        action = "default"
    # apply action
    if action == "ignore":
        return None
    elif action in ["always", "default", "module", "once"]:
        if key in _POSYDON_WARNINGS_REGISTRY:
            _POSYDON_WARNINGS_REGISTRY[key] += 1
            if action == "default":
                return None
        else:
            _POSYDON_WARNINGS_REGISTRY[key] = 1
    return warning

def _issue_warn(warning=dict(message="No warning")):
    """Issue a warning.

        Parameters
        ----------
        warning : dict (default: {message="No warning"})
            Dictionary containing all options passed to warnings.warn.
    """
    filtered_warning = _apply_POSYDON_filter(warning)
    if filtered_warning is None:
        return
    else:
        warnings.warn(**filtered_warning)


class _Catched_POSYDON_Warnings:
    """Class which stores catched warnings."""
    def __init__(self, catch_warnings=False, record=True, filter_first=True):
        """Constructor of the object.

        Parameters
        ----------
        catch_warnings : bool (default: False)
            Determines, whether warnings are catched.
        record : bool (default: True)
            Determines, whether warnings are recorded.
        filter_first : bool (default: True)
            Determines, whether warnings are filtered before recorded or
            discarded.
        """
        self.catch_warnings = catch_warnings
        self.catched_warnings = []
        self.record = record
        self.filter_first = filter_first
        self._got_called = False
        if not isinstance(catch_warnings, bool):
            raise TypeError("catch_warnings must be a boolean.")
        if not isinstance(record, bool):
            raise TypeError("record must be a boolean.")
        if not isinstance(filter_first, bool):
            raise TypeError("filter_first must be a boolean.")
    
    def __str__(self):
        """Return the status of the object as a string."""
        if self.catch_warnings:
            if self.record:
                ret = "POSYDON warnings will be catched and recorded."
            else:
                ret = "POSYDON warnings will be catched and discarded."
        else:
            ret = "POSYDON warnings are shown."
        ncatched = len(self.catched_warnings)
        if ncatched>0:
            if ncatched==1:
                ret += "\nThere is 1 warning recorded."
            else:
                ret += "\nThere are {} warnings recorded.".format(ncatched)
        return ret

    def __call__(self, new_warning=None, empty_cache=False,
                 change_settings=None):
        """Deal with changes and new warnings

        Parameters
        ----------
        new_warning : dict or None (default: None)
            Dictionary containing all options passed to warnings.warn.
        empty_cache : bool (default: False)
            If True, the catched warnings will be reset first.
        change_settings : dict or None (default: None)
            The dictionary can contain any attribute this class has. The
            attributes get the new values. (Changing the filter_first, while
            there are recorded warnings, the statistics may count warnings
            twice or not at all.)
        """
        self._got_called = True
        # Check if there is anything to do
        if not empty_cache and new_warning is None and change_settings is None:
            raise ValueError("Nothing to do: either empty_cache has to be True"
                            " or new_warning/change_settings needs to be set.")
        if empty_cache:
            # Clear the cache
            self.reset_cache()
        if change_settings is not None:
            if not isinstance(change_settings, dict):
                raise TypeError("change_settings has to be a dict or None.")
            # Change attributes
            for attr,val in change_settings.items():
                if attr=="catched_warnings":
                    # Protect the list of warnings from changes
                    continue
                if hasattr(self, attr):
                    if isinstance(val,type(getattr(self, attr))):
                        setattr(self, attr, val)
                    else:
                        raise TypeError(f"{attr} has to be a "
                                        f"{type(getattr(self, attr))}.")
                else:
                    raise AttributeError(f"{attr} unknown to "
                                         "_Catched_POSYDON_Warnings.")
            if ((self.catch_warnings==False) and
                (len(self.catched_warnings)>0)):
                # If there are recorded warnings issue them and empty the list
                for w in self.catched_warnings:
                    w["stacklevel"] += 2
                    if self.filter_first:
                        warnings.warn(**w)
                    else:
                        _issue_warn(w)
                self.catched_warnings = []
        if new_warning is not None:
            # Process new warning
            if self.catch_warnings:
                # Catch warning
                if self.filter_first:
                    # Apply POSYDON filtering/stats
                    new_warning = _apply_POSYDON_filter(new_warning)
                if self.record and new_warning is not None:
                    # Add warning to catched ones
                    self.catched_warnings.append(new_warning)
            else:
                # Issue warning
                new_warning["stacklevel"] += 2
                _issue_warn(new_warning)
    
    def __del__(self):
        """Destructor of the object. It will issue still recorded warnings."""
        if len(self.catched_warnings)>0:
            # If there are recorded warnings issue them.
            self.catch_warnings = False
            print("There are still recorded warnings:")
            for w in self.catched_warnings:
                w["stacklevel"] = 2
                if self.filter_first:
                    warnings.warn(**w)
                else:
                    _issue_warn(w)
    
    def got_called(self):
        """Returns, whether the object got called."""
        return self._got_called
    
    def has_records(self):
        """Checks whether there are records of catched warnings.

        Returns
        -------
        True if there are recorded warnings otherwise False.
        """
        return len(self.catched_warnings)>0
    
    def get_cache(self, empty_cache=False):
        """Get catched warnings.

        Parameters
        ----------
        empty_cache : bool (default: False)
            If True, the catched warnings will be reset.

        Returns
        -------
        List of recorded warnings.
        """
        cache = copy.copy(self.catched_warnings)
        if empty_cache:
            self.reset_cache()
        return cache
    
    def reset_cache(self):
        """Resets the catched warnings."""
        self.catched_warnings = []

_CATCHED_POSYDON_WARNINGS = _Catched_POSYDON_Warnings()


class Catch_POSYDON_Warnings:
    """Context manager class to catch POSYDON warnings."""
    def __init__(self, catch_warnings=True, record=True, filter_first=True):
        """Constructor of the object.

        Parameters
        ----------
        catch_warnings : bool (default: False)
            Determines, whether warnings are catched.
        record : bool (default: True)
            Determines, whether warnings are recorded.
        filter_first : bool (default: True)
            Determines, whether warnings are filtered before recorded or
            discarded.
        """
        self.catch_warnings = catch_warnings
        self.record = record
        self.filter_first = filter_first

    def __enter__(self):
        """Enable catching."""
        global _CATCHED_POSYDON_WARNINGS
        _CATCHED_POSYDON_WARNINGS(change_settings={
                                      'catch_warnings': self.catch_warnings,
                                      'record': self.record,
                                      'filter_first': self.filter_first,
                                      '_got_called': False})
        return _CATCHED_POSYDON_WARNINGS
    
    def __exit__(self, exc_type, exc_value, exc_traceback):
        """Disable catching."""
        global _CATCHED_POSYDON_WARNINGS
        _CATCHED_POSYDON_WARNINGS(change_settings={'catch_warnings': False})
        return False


def Pwarn(message, category=None, stacklevel=2, **kwargs):
    """Issueing a warning via warnings.warn.

        Parameters
        ----------
        message : str
            The message printed in the warning.
        category : str or a warnings category class (default: None)
            A string, which will be converted into a POSYDONWarning class.
        stacklevel : int (default: 2)
            The stack level passed to warnings.warn, defaults to 2.
        **kwargs : dict (optional)
            Dictionary containing extra options passed to warnings.warn.
    """
    global _CATCHED_POSYDON_WARNINGS
    if not isinstance(message, str):
        raise TypeError("message must be a string.")
    if not isinstance(stacklevel, int):
        raise TypeError("stacklevel must be an integer.")
    if not(isinstance(category, type) and issubclass(category, Warning)):
        category = _get_POSYDONWarning_class(category)
    if (isinstance(category, type) and issubclass(category, POSYDONWarning)):
        # deal with POSYDON warnings
        _CATCHED_POSYDON_WARNINGS(new_warning=dict({"message": message,
                                                   "category": category,
                                                   "stacklevel": stacklevel},
                                                  **kwargs))
    else:
        # deal with python warnings
        warnings.warn(message=message, category=category,
                      stacklevel=stacklevel, **kwargs)

def SetPOSYDONWarnings(action="default", category=POSYDONWarning, **kwargs):
    """Add the warnings filter for POSYDON warnings.

        Parameters
        ----------
        action : str (default: "default")
            The behaviour for those warnings. Should be one out of
            "default", "error", "ignore", "always", "module", or "once"
        category : str or a warnings category class (default: POSYDONWarning)
            A string, which will be converted into a POSYDONWarning class.
        **kwargs : dict (optional)
            Dictionary containing extra options passed to
            warnings.filterwarnings.
    """
    category = _get_POSYDONWarning_class(category)
    if isinstance(category, type) and issubclass(category, POSYDONWarning):
        warnings.filterwarnings(action=action, category=category, **kwargs)

def NoPOSYDONWarnings(category=POSYDONWarning):
    """Switch the warnings filter to ignore for POSYDON warnings.

        Parameters
        ----------
        category : str or a warnings category class (default: POSYDONWarning)
            A string, which will be converted into a POSYDONWarning class.
    """
    SetPOSYDONWarnings(action="ignore", category=category)

def AllPOSYDONWarnings(category=POSYDONWarning):
    """Switch the warnings filter to always for POSYDON warnings.

        Parameters
        ----------
        category : str or a warnings category class (default: POSYDONWarning)
            A string, which will be converted into a POSYDONWarning class.
    """
    SetPOSYDONWarnings(action="always", category=category)


# The base class of all warnings is "Warning", which is derived from Exception
#
# List of Python's default warnings (sub)classes in alphabetic order
# BytesWarning: Base category for warnings related to bytes and bytearray.
# DeprecationWarning: Base category for warnings about deprecated features when
#     those warnings are intended for other Python developers (ignored by
#     default, unless triggered by code in __main__).
# EncodingWarning: Base class for warnings related to encodings.
# FutureWarning: Base category for warnings about deprecated features when
#     those warnings are intended for end users of applications that are
#     written in Python.
# ImportWarning: Base category for warnings triggered during the process of
#     importing a module (ignored by default).
# PendingDeprecationWarning: Base category for warnings about features that
#     will be deprecated in the future (ignored by default).
# ResourceWarning: Base category for warnings related to resource usage
#     (ignored by default).
# RuntimeWarning: Base category for warnings about dubious runtime features.
# SyntaxWarning: Base category for warnings about dubious syntactic features.
# UnicodeWarning: Base category for warnings related to Unicode.
# UserWarning: The default category for warn().
# 
# for more details see https://docs.python.org/3/library/exceptions.html#warnings
# and https://docs.python.org/3/library/warnings.html
# https://github.com/python/cpython/blob/3.11/Lib/warnings.py
