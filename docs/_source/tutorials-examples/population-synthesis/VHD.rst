.. _VHD:

Van den Heuvel diagrams with POSYDON
====================================

.. warning::
   This is an **experimental feature** of POSYDON. Unexpected behaviors or results may occur.

.. note::
   To use the visualization functionalities, install the optional visualization modules `pip install .[vis]`


`VHdiagram` provides a visual representation of individual POSYDON binaries, offering a more intuitive sense of their properties. This tutorial uses the 'population.h5' dataset as an example. TODO: this need to be tested on a v2.0.0 dataset.

Visualizing a Specific Binary
-----------------------------

Simple Usage
~~~~~~~~~~~~

To visualize a specific binary within an interactive window, use the following code:

.. code-block:: python

    from posydon.visualization.VHdiagram import VHdiagram
    VHdiagram('population.h5', path='./dataset/', index=18976)

.. image:: pictures/detailled_window.png

.. tip::
   Use the 'path' named parameter if your dataset resides in a different directory.

The options window, accessible via the 'option' button, lets you customize the view:

.. image:: pictures/option_window.png

Clicking the 'save' button captures a screenshot of the current view and saves it to a newly created 'screens' folder in the current directory.

Predefined Views
~~~~~~~~~~~~~~~~

You can choose from four predefined views for faster visualization:

.. code-block:: python

    from posydon.visualization.VHdiagram import VHdiagram
    from posydon.visualization.VH_diagram.PresenterMode import PresenterMode

    VHdiagram('population.h5', index=19628, presentMode=PresenterMode.DIAGRAM)

.. image:: pictures/diagram_window.png

Display Modes
~~~~~~~~~~~~~

You have the flexibility to display the diagram either within an interactive window or inline within a Jupyter notebook. The available display modes are:

.. code-block:: python

    from posydon.visualization.VHdiagram import VHdiagram, DisplayMode
    from posydon.visualization.VH_diagram.Presenter import PresenterMode

    VHdiagram(
        "population.h5",
        index=19628,
        presentMode=PresenterMode.DIAGRAM,
        displayMode=DisplayMode.INLINE_B,
    )

.. image:: pictures/diagram_inline.png

Advanced Visualization Techniques
---------------------------------

Visualizing Multiple Binaries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The `VDdiagramm_m` module allows multiple binary visualizations to be arranged horizontally in a single plot:

.. code-block:: python

    from posydon.visualization.VH_diagram.PresenterMultiple import VHdiagramm_m
    from posydon.visualization.VHdiagram import DisplayMode
    from posydon.visualization.VH_diagram.PresenterMode import PresenterMode

    VHD = VHdiagramm_m('./data/population.h5',
                    index=cnt[:,0],
                    frequency=parse_df.get_frequencies(),
                    hierarchy=False,
                    presentMode=PresenterMode.DIAGRAM,
                    displayMode=DisplayMode.INLINE_B)

.. image:: pictures/diagram_multiple.png

Hierarchical Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~

This visualization style aggregates identical steps into a tree plot where nodes represent common steps. Each node is labeled with percentages relative to the parent node percentage.

.. code-block:: python

    from posydon.visualization.VH_diagram.PresenterMultiple import VHdiagramm_m
    from posydon.visualization.VHdiagram import DisplayMode
    from posydon.visualization.VH_diagram.PresenterMode import PresenterMode

    VHD = VHdiagramm_m('./data/population.h5',
                    index=cnt[:,0],
                    frequency=parse_df.get_frequencies(),
                    hierarchy=True,
                    presentMode=PresenterMode.DIAGRAM,
                    displayMode=DisplayMode.INLINE_B)

.. image:: pictures/diagram_hierarchy.png

Binary Analysis Tools
---------------------

Counting Binaries
~~~~~~~~~~~~~~~~~

Use the `ParseDataFrame` class to iterate through the binary file, counting identical binary simulations:

.. code-block:: python

    from posydon.visualization.VH_diagram.ParseDataFrame import ParseDataFrame

    parse_df = ParseDataFrame('./data/population.h5')
    parse_df.count_dict

>>> Counter({...})

Sorting Binaries
~~~~~~~~~~~~~~~~

You can sort binaries based on filenames of images representing their steps:

.. code-block:: python

    VHD = VHdiagramm_m('./data/population.h5',
                    index=VHD.get_sorted_index(),
                    frequency=parse_df.get_frequencies(),
                    hierarchy=False,
                    presentMode=PresenterMode.DIAGRAM,
                    displayMode=DisplayMode.INLINE_B)

.. image:: pictures/diagram_multiple_sort.png

