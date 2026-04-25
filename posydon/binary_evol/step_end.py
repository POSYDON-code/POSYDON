"""Class defining the default step at the end of evolution of a binary."""


__authors__ = [
    "Simone Bavera <Simone.Bavera@unige.ch>",
]


class step_end:
    """Default end step."""

    def __call__(self, binary):
        """Change the event of the binary to 'end'."""
        if binary.state == "disrupted":
            binary.update_star_states()
        binary.event = 'END'
