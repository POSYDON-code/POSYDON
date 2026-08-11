.. _grid-types:

################
Grid Types
################

A POSYDON *grid* is a collection of MESA simulations of stellar or binary
systems that share a common set of physical assumptions and differ in
their initial conditions. Understanding how grids are organized is key to
knowing which data you are using and how you can extend it with your own runs.

What grids are there?
----------------------

POSYDON ships with different grid types:

- ``HMS-HMS``_ grids: hydrogen-rich main-sequence primaries with hydrogen-rich main-sequence secondaries
- ``HMS-HMS_RLO``_ grids: hydrogen-rich main-sequence primaries with hydrogen-rich main-sequence secondaries, starting from Roche-lobe overflow (RLO)
- ``HMS-CO_RLO``_ grids: hydrogen-rich main-sequence primaries with a compact object companion (CO), starting from Roche-lobe overflow (RLO)
- ``HeMS-CO``_ grids: helium-rich main-sequence primaries with a compact object companion (CO)
- ``HeMS-CO_RLO``_ grids: helium-rich main-sequence primaries with a compact object companion (CO), starting from Roche-lobe overflow (RLO)

- ``single_HMS``_ grids: single hydrogen-rich main-sequence stars
- ``single_HeMS``_ grids: single helium-rich main-sequence stars

What a grid contains
--------------------

Each grid contains a set of MESA simulations that are down-sampled to save disk space.
For each run, the grid contains the ``initial_values``, ``final_values``, and ``history``
of the simulation. The ``history`` is split between the ``binary_history`` and
the histories of the individual stars (``history1`` and ``history2``).

The grid therefore spans a defined region of this parameter space, sampled on a
rectangular (grid) of points. Each run type has its own file, so within a given
metallicity you will find one dataset per type.


Fixed grids
~~~~~~~~~~~

A **fixed** (or rectilinearly sampled) grid is one in which the initial
conditions are chosen on a fixed rectangular lattice of values ahead of time.
Every point in the lattice is simulated, giving a uniform and predictable
coverage of the parameter space. The published data releases (DR1, DR2) are
fixed grids. You can run your own fixed grids with the POSYDON MESA grid API &mdash;
see :ref:`fixed_grid` and :ref:`MESA-grids <MESA-grids>`.


Example organization of downloaded data
-----------------------------------------

A typical data directory is organized by run type and metallicity, e.g.::

    POSYDON_data/
        HMS-HMS/
            1e+00_Zsun.h5
            ...
        HMS-HMS_RLO/
            ...
        single_HMS/
            ...

where each ``.h5`` file contains one ``PSyGrid`` for a given run type and
metallicity.

Summary
-------

.. seealso::
   * :ref:`data-releases` for the published datasets and their storage costs.
   * :ref:`MESA-grids <MESA-grids>` for the tutorials on running your own grids.
