.. _grid-types:

################
Grid Types
################

A POSYDON *grid* is a collection of MESA simulations of stellar or binary
systems that share a common set of physical assumptions and differ only in
their initial conditions. Understanding how grids are organized is key to
knowing which data you are using and how you can extend it with your own runs.

What a grid contains
--------------------

Each run inside a grid corresponds to one set of initial conditions. For a
binary grid those are, at a minimum:

* the initial mass of the primary star,
* the initial mass of the secondary star,
* the initial orbital period (or separation),
* the initial metallicity.

The grid therefore spans a defined region of this parameter space, sampled on a
rectangular (grid) of points. The terms "mass ratio" and "orbital period" that
appear in the names of some datasets (such as ``HMS-HMS``) refer to how that
rectilinear sampling is laid out.

For single-star grids the analogous parameters are just the initial mass and
metallicity of the star.

Fixed vs. dynamic grids
-----------------------

POSYDON distinguishes two ways of choosing which initial conditions to
simulate:

Fixed grids
~~~~~~~~~~~

A **fixed** (or rectilinearly sampled) grid is one in which the initial
conditions are chosen on a fixed rectangular lattice of values ahead of time.
Every point in the lattice is simulated, giving a uniform and predictable
coverage of the parameter space. The published data releases (DR1, DR2) are
fixed grids. You can run your own fixed grids with the POSYDON MESA grid API &mdash;
see :ref:`fixed_grid` and :ref:`MESA-grids <MESA-grids>`.

Dynamic grids
~~~~~~~~~~~~~~~~

A **dynamic** grid starts from a sparse set of simulations and then uses the
results to decide, adaptively, where extra simulations would be most valuable,
adding them on the fly. This improves the coverage of regions where the
behavior is complex without wasting effort in homogeneous parts of the space,
and consequently improves classification and interpolation accuracy. Dynamic
grids are **experimental** in POSYDON and not part of the published data
releases. More details are in :ref:`dynamic_grid`.

Grid run types
-----------------

The ``HMS-HMS`` in the dataset names refers to the type of binaries being
simulated. The run type specifies which kind of star / binary the grid models
(e.g. ``HMS-HMS`` for a hydrogen-rich main-sequence primaries with a hydrogen-rich
main-sequence secondary, or a single-star type for single stars). Each run type
has its own file, so within a given metallicity you will find one dataset per
type. The run type defines both the physics adopted and the range of initial
conditions that the grid spans.

Example organization of downloaded data
-----------------------------------------

A typical data directory is organized by run type and metallicity, e.g.::

    POSYDON_data/
        HMS-HMS/
            1e+00_Zsun.h5
            1e-01_Z.LTE_nearby.h5
            ...
        He-He/
            ...
        single-star/
            ...

where each ``.h5`` file contains one ``PSyGrid`` for a given run type and
metallicity.

Summary
-------

* A grid = many MESA simulations over a range of initial conditions.
* **Fixed** grids sample on a fixed lattice and are what the published data
  uses.
* **Dynamic** grids refine adaptively and are experimental.
* **Run types** (e.g. ``HMS-HMS``) name what kind of binaries a grid simulates.

.. seealso::
   * :ref:`data-releases` for the published datasets and their storage costs.
   * :ref:`MESA-grids <MESA-grids>` for the tutorials on running your own grids.
