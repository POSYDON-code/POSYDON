"""Class for handling ignore reasons when creating new grids."""


IGNORE_REASONS_PRIORITY = [
    # missing data files
    'ignored_no_history1',              # history1 is always needed
    'ignored_no_binary_history',        # needed for all binaries
    # corrupted data files
    'corrupted_history1',
    'corrupted_binary_history',
    'corrupted_history2',
    # issues with specific data files
    'ignored_scrubbed_history',         # issues when scrubbing history data
    'ignored_no_final_profile',         # issues with final profile data
    # binaries to be ignored for specialized grids
    'ignored_no_RLO'
]


class IgnoreReason():
    """Class for handling ignore reasons when creating new grids."""

    def __init__(self):
        self.reason = None
        self.order = None

    def __bool__(self):
        """Return True if an ignore reason has been set."""
        return self.reason is not None

    def __setattr__(self, name, value):
        """Only allow to set an ignore reason, and only from a defined list."""
        if name != "reason":
            return

        # allow resetting the ignore reason
        if value is None:
            self.__dict__["reason"] = None
            self.__dict__["order"] = None
            return

        # check if new value is valid
        if value not in IGNORE_REASONS_PRIORITY:
            raise ValueError(f"Ignore reason `{value}` not recognized.")

        new_order = IGNORE_REASONS_PRIORITY.index(value)

        # if previous reason was None, or higher priority... update
        if ((self.__dict__["reason"] is None)
                or (new_order < self.__dict__["order"])):
            self.__dict__["reason"] = value
            self.__dict__["order"] = new_order
            return
