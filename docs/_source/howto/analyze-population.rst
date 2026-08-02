.. _analyze-population:

Analyzing a Population
======================

After a population run finishes, the results are stored as HDF5 files describing
the evolved binaries. This page walks you through loading them, filtering and
selecting subsets, and getting them into a form you can work with (in-memory
DataFrames). It also points you to the *Van den Heuvel diagram* (VHD)
visualization.

Loading the output
------------------

Population outputs are read through the :ref:`Population
<synthetic-population>` interface, which wraps the HDF5 files produced by a
run. Load a file by passing its path to the constructor:

.. code-block:: python

   from posydon.popsyn import Population
   pop = Population('path/to/population.h5')

The object does not read everything at once; data is loaded from the file on
demand. Access the full evolutionary history of every system through
``pop.history``.

You can then access the underlying histories of every binary from the
population structure.

Selecting and filtering
-------------------------------------

To work with a subset of binaries, use the ``select`` method on the
population's history and choose the binaries you want, e.g. binaries whose
remnant masses pass a threshold:

.. code-block:: python

   sel = pop.history.select(columns=['S1_mass'], where='S1_mass > 20')

``select`` supports the columns that supported by the selection, plus the indices and
column names. If you try to select on a column that is not supported, you will
get an error like::

   ValueError: The passed where expression S1_mass > 10
   contains an invalid variable reference. All of the variable references
   must be a reference to an axis (e.g. 'index' or 'columns'), or a data_column.
   The currently defined references are: index, columns, state, event,
   step_names, S1_state, S2_state

Only specific columns are supported for selection (string columns, the indices,
and the column names). For everything else, do the selection in-memory instead.

Working in memory with DataFrames
----------------------------------

If you need to filter on columns that ``where`` does not support, or you want
more flexible data manipulation, convert the relevant history to a Pandas
DataFrame and filter/print with the usual tools:

.. code-block:: python

   df = pop.history[:]
   my_subset = df[df['S1_mass'] > 10]
   print(my_subset)

Once you have a DataFrame you can sort, aggregate, and export the subset with
the full comfort of Pandas.

Visualizing with Van-den-Heuvel diagrams
----------------------------------------

POSYDON provides the Van den Heuvel diagram (VHD) for visualizing the
evolution of single binaries or whole populations, including several
interactive modes and analysis tools. See the VHD tutorial under
``tutorials-examples/population-synthesis/`` for a guided tour.

Related
---------------

* :ref:`synthetic-population` reference for the full analysis API.
* The :ref:`binary population synthesis tutorial <binary-pop-syn>` for
  end-to-end examples including rate and observation computation.
* The Van den Heuvel diagram (VHD) tutorial under
  ``tutorials-examples/population-synthesis/``.
