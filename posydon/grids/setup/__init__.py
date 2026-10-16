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
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]
