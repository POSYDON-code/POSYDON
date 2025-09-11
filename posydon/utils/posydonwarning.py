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
the caught warnings via `get_cache`.
Here an example, where the argument True of get_cache, will empty the list of
caught warnings after it is copy to be returned:
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
    "Seth Gossage <seth.gossage@northwestern.edu>"   
]


import copy
import sys
import warnings


def nosrc_code_format(message, category, filename, lineno, line=None):
    """
    This sets the warning format to not include the source code line.
    """
    return f"{filename}:{lineno}: {category.__name__}: {message}\n"

# Setting the warning format to use the above format function
warnings.formatwarning = nosrc_code_format


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
        
class SFHModelWarning(POSYDONWarning):
    """Warnings related to the SFH model."""
    def __init__(self, message=''):
        super().__init__(message)
        
class ValueWarning(POSYDONWarning):
    """Warnings related to a ValueError."""
    def __init__(self, message=''):
        super().__init__(message)


# All POSYDON warnings subclasses should be defined beforehand
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


# Defining an own resistry for the warnings, which will additionally save how
# often a warning occured. This is simply a dictionary with the warning
# characteristics as key pointing to an integer given the count.
_POSYDON_WARNINGS_REGISTRY = {}

def get_stats():
    """Return the statistics of the POSYDON warnings."""
    global _POSYDON_WARNINGS_REGISTRY
    return _POSYDON_WARNINGS_REGISTRY

def print_stats():
    """Prints the statistics of the POSYDON warnings."""
    global _POSYDON_WARNINGS_REGISTRY
    if len(_POSYDON_WARNINGS_REGISTRY)==0:
        print("No POSYDON warnings occured.")
    else:
        print("There have been POSYDON warnings in the global registry:\n",
              _POSYDON_WARNINGS_REGISTRY)

def _apply_POSYDON_filter(warning=dict(message="No warning"), registry=None):
    # In python warnings, this functionality is spead over two functions. It
    # compares the characteristics of a warning with the ones in the registry.
    # We use the warnings class and the code position (filename + line number)
    # as characteristics, while python's standard omits the code filename but
    # therefore adds the full warnings text, which causes different keys for
    # "i is negative, i=-1" and "i is negative, i=-3" to be different warnings
    # in python, but considered to have the same warning characteristics for
    # POSYDON.
    """Filter a warning.

        Parameters
        ----------
        warning : dict (default: {message="No warning"})
            Dictionary containing all options passed to warnings.warn.
        registry : dict or None (default: None)
            Warnings registry. If None, use the global one.

        Returns
        -------
        warning or None in case it got filtered out.
    """
    # Check registry and warning
    if registry is None:
        global _POSYDON_WARNINGS_REGISTRY
        if not isinstance(_POSYDON_WARNINGS_REGISTRY, dict): # pragma: no cover
            raise TypeError("_POSYDON_WARNINGS_REGISTRY is corrupt and can't "
                            "be used for the registry, hence registry can't be"
                            " None.")
        registry = _POSYDON_WARNINGS_REGISTRY
    elif not isinstance(registry, dict):
        raise TypeError("registry must be a dictionary or None.")
    if not isinstance(warning, dict):
        raise TypeError("warning must be a dictionary.")
    # Get stack level
    stacklevel = warning.get('stacklevel', 0)
    if not isinstance(stacklevel, int):
        raise TypeError("stacklevel must be an integer.")
    # Get category
    category = warning.get('category', None)
    if not(isinstance(category, type) and issubclass(category, Warning)):
        category = _get_POSYDONWarning_class(category)
    if not isinstance(category, type):
        # A non valid category will default to a Userwarning
        category = UserWarning
    # Get filename and lineno from frame at stacklevel
    try:
        frame = sys._getframe(stacklevel)
    except: # pragma: no cover
        g = sys.__dict__
        filename = "sys"
        lineno = 1
    else:
        g = frame.f_globals
        filename = frame.f_code.co_filename
        lineno = frame.f_lineno
    # Get module
    if isinstance(g, dict):
        module = g.get('__name__', "posydonwarnings")
    else: # pragma: no cover
        module = "posydonwarnings"
    # Set key for registry:
    # We do not use the warnings text, to allow it to contain detailed
    # information, while still identifying warnings with same origin
    key = (category, filename, lineno)
    # Get message text
    text = warning.get('message', "")
    if not isinstance(text, str):
        raise TypeError("message must be a string.")
    # Search the filters:
    # Here we still use the python filters
    for item in warnings.filters:
        action, msg, cat, mod, ln = item
        if ((msg is None or msg.match(text)) and
            issubclass(category, cat) and
            (mod is None or mod.match(module)) and
            (ln == 0 or lineno == ln)):
            break
    else:
        action = "default"
    # Apply action
    if action == "ignore":
        return None
    elif action in ["always", "default", "module", "once"]:
        # Compared to standard python not only save a occurence but count them
        if key in registry:
            registry[key] += 1
            if action in ["default", "module", "once"]:
                return None
        else:
            registry[key] = 1
    return warning

def _issue_warn(warning=dict(message="No warning"), registry=None):
    """Issue a warning.

        Parameters
        ----------
        warning : dict (default: {message="No warning"})
            Dictionary containing all options passed to warnings.warn.
        registry : dict or None (default: None)
            Warnings registry. If None, use the global one.
    """
    filtered_warning = _apply_POSYDON_filter(warning, registry)
    if filtered_warning is None:
        return
    else:
        warnings.warn(**filtered_warning)


class _Caught_POSYDON_Warnings:
    """Class which stores caught warnings."""
    def __init__(self, catch_warnings=False, record=True, filter_first=True,
                 registry=None):
        """Constructor of the object.

        Parameters
        ----------
        catch_warnings : bool (default: False)
            Determines, whether warnings are caught.
        record : bool (default: True)
            Determines, whether warnings are recorded.
        filter_first : bool (default: True)
            Determines, whether warnings are filtered before recorded or
            discarded.
        registry : dict or None (default: None)
            Warnings registry. If None, use the global one.
        """
        self.catch_warnings = catch_warnings
        self.caught_warnings = []
        self.record = record
        self.filter_first = filter_first
        self._got_called = False
        if registry is None:
            global _POSYDON_WARNINGS_REGISTRY
            self.registry = _POSYDON_WARNINGS_REGISTRY
        else:
            self.registry = registry
        if not isinstance(catch_warnings, bool):
            raise TypeError("catch_warnings must be a boolean.")
        if not isinstance(record, bool):
            raise TypeError("record must be a boolean.")
        if not isinstance(filter_first, bool):
            raise TypeError("filter_first must be a boolean.")
        if not isinstance(self.registry, dict):
            raise TypeError("registry must be a dictionary.")
    
    def __str__(self):
        """Return the status of the object as a string."""
        if self.catch_warnings:
            if self.record:
                ret = "POSYDON warnings will be caught and recorded."
                if self.filter_first:
                    ret += " Filters are applied before recording."
            else:
                ret = "POSYDON warnings will be caught and discarded."
        else:
            ret = "POSYDON warnings are shown."
        ncaught = len(self.caught_warnings)
        if ncaught>0:
            if ncaught==1:
                ret += "\nThere is 1 warning recorded."
            else:
                ret += "\nThere are {} warnings recorded.".format(ncaught)
        global _POSYDON_WARNINGS_REGISTRY
        if self.registry!=_POSYDON_WARNINGS_REGISTRY:
            ret += " Currently a private registry is used, it contains:\n"
            ret += "{}".format(self.registry)
        return ret

    def __call__(self, new_warning=None, empty_cache=False,
                 change_settings=None):
        """Deal with changes and new warnings

        Parameters
        ----------
        new_warning : dict or None (default: None)
            Dictionary containing all options passed to warnings.warn.
        empty_cache : bool (default: False)
            If True, the caught warnings will be reset first.
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
                if attr=="caught_warnings":
                    # Protect the list of warnings from changes
                    continue
                if attr=="registry":
                    # special treatment of registry
                    if val is None:
                        global _POSYDON_WARNINGS_REGISTRY
                        self.registry = _POSYDON_WARNINGS_REGISTRY
                    else:
                        self.registry = change_settings[attr]
                    continue
                if hasattr(self, attr):
                    if isinstance(val,type(getattr(self, attr))):
                        setattr(self, attr, val)
                    else:
                        raise TypeError(f"{attr} has to be a "
                                        f"{type(getattr(self, attr))}.")
                else:
                    raise AttributeError(f"{attr} unknown to "
                                         "_Caught_POSYDON_Warnings.")
            if ((self.catch_warnings==False) and
                (len(self.caught_warnings)>0)):
                # If there are recorded warnings issue them and empty the list
                for w in self.caught_warnings:
                    if "stacklevel" in w:
                        w["stacklevel"] += 2
                    else:
                        w["stacklevel"] = 2
                    if self.filter_first:
                        warnings.warn(**w)
                    else:
                        _issue_warn(w, self.registry)
                self.caught_warnings = []
        if new_warning is not None:
            # Process new warning
            if self.catch_warnings:
                # Catch warning
                if self.filter_first:
                    # Apply POSYDON filtering/stats
                    new_warning = _apply_POSYDON_filter(new_warning,
                                                        self.registry)
                if self.record and new_warning is not None:
                    # Add warning to caught ones
                    self.caught_warnings.append(new_warning)
            else:
                # Issue warning
                new_warning["stacklevel"] += 2
                _issue_warn(new_warning, self.registry)
    
    def __del__(self):
        """Destructor of the object. It will issue still recorded warnings."""
        if len(self.caught_warnings)>0:
            # If there are recorded warnings issue them.
            self.catch_warnings = False
            print("There are still recorded warnings:", file=sys.stderr)
            for w in self.caught_warnings:
                w["stacklevel"] = 2
                if self.filter_first:
                    warnings.warn(**w)
                else:
                    _issue_warn(w, self.registry)
    
    def got_called(self):
        """Returns, whether the object got called."""
        return self._got_called
    
    def has_records(self):
        """Checks whether there are records of caught warnings.

        Returns
        -------
        True if there are recorded warnings otherwise False.
        """
        return len(self.caught_warnings)>0
    
    def get_cache(self, empty_cache=False):
        """Get caught warnings.

        Parameters
        ----------
        empty_cache : bool (default: False)
            If True, the caught warnings will be reset.

        Returns
        -------
        List of recorded warnings.
        """
        cache = copy.copy(self.caught_warnings)
        if empty_cache:
            self.reset_cache()
        return cache
    
    def reset_cache(self):
        """Resets the caught warnings."""
        self.caught_warnings = []

# Here we store all our caught POSYDON warnings
_CAUGHT_POSYDON_WARNINGS = _Caught_POSYDON_Warnings()


class Catch_POSYDON_Warnings:
    # We use our own context manager, which does not overwrite functions.
    # Instead it can only catch POSYDON warnings issued via the Pwarn function.
    """Context manager class to catch POSYDON warnings."""
    def __init__(self, catch_warnings=True, record=True, filter_first=True,
                 own_registry=False, use_python_catch=False):
        """Constructor of the object.

        Parameters
        ----------
        catch_warnings : bool (default: False)
            Determines, whether warnings are caught.
        record : bool (default: True)
            Determines, whether warnings are recorded.
        filter_first : bool (default: True)
            Determines, whether warnings are filtered before recorded or
            discarded.
        own_registry : bool (default: False)
            Determines, whether the global POSYDON warnings registry should be
            used or an own one only valid within the context.
        use_python_catch : bool (default: False)
            If enabled, it put the python catch_warnings on top of the POSYDON
            catches. This brings all the drawbacks of the standard catching,
            hence it is strongly recommended to be not used (properly written
            code will not contain any use case for this option, hence it is
            only for backward compatibility.).
        """
        self.catch_warnings = catch_warnings
        self.record = record
        self.filter_first = filter_first
        if own_registry:
            self.context_registry = {}
        else:
            # no own registry will use the global _POSYDON_WARNINGS_REGISTRY
            self.context_registry = None
        if use_python_catch:
            self.python_catch = warnings.catch_warnings(record=self.record)
        else:
            self.python_catch = None

    def __enter__(self):
        """Enable catching."""
        global _CAUGHT_POSYDON_WARNINGS
        _CAUGHT_POSYDON_WARNINGS(change_settings={
                                     'catch_warnings': self.catch_warnings,
                                     'record': self.record,
                                     'filter_first': self.filter_first,
                                     '_got_called': False,
                                     'registry': self.context_registry})
        if isinstance(self.python_catch, warnings.catch_warnings):
            # enter catch of python as well (should be done last)
            self.python_catch.__enter__()
        return _CAUGHT_POSYDON_WARNINGS
    
    def __exit__(self, exc_type, exc_value, exc_traceback):
        """Disable catching."""
        if isinstance(self.python_catch, warnings.catch_warnings):
            # exit catch of python as well (should be done first)
            self.python_catch.__exit__()
            # reset it (needed because it cannot enter the context for the same
            # object twice)
            self.python_catch = None
        global _CAUGHT_POSYDON_WARNINGS
        # If the cache is not cleared before, it will issue all recorded
        # warnings.
        _CAUGHT_POSYDON_WARNINGS(change_settings={'catch_warnings': False,
                                                  'registry': None})
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
    global _CAUGHT_POSYDON_WARNINGS
    if not isinstance(message, str):
        raise TypeError("message must be a string.")
    if not isinstance(stacklevel, int):
        raise TypeError("stacklevel must be an integer.")
    if not(isinstance(category, type) and issubclass(category, Warning)):
        category = _get_POSYDONWarning_class(category)
    if (isinstance(category, type) and issubclass(category, POSYDONWarning)):
        # deal with POSYDON warnings
        _CAUGHT_POSYDON_WARNINGS(new_warning=dict({"message": message,
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
