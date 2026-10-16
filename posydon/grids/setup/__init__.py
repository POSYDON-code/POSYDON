"""Building blocks for setting up MESA grids.

This package is the rewrite of ``posydon-setup-grid``.

Modules
-------
inlist
    Read MESA inlists into dictionaries and write them back.
layers
    Stack several inlists into one, in a fixed and named order.
config
    The settings of a grid (GridConfig), read from a new-format ini file.
legacy
    Read today's grid ini files into the same GridConfig.
grid_types
    The registry of grid types and what makes each one different.
scenario
    Get the inlist repository at the scenario's commit and list its files.
build
    Stack the layers, add the generated keys and write a grid's inlists.
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]
