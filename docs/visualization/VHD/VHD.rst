.. _VHD:

#######################
Van den Heuvel diagrams
#######################

Visualize specific index
========================

`VHdiagrams` allow individual POSYDON binaries to be viewed in a more intuitive
sense. We use 'population.h5' as an example dataset.


Simple usage
------------

The basic usage is to simply visualize one index inside an independent and
reactive window, with the following :

.. code-block:: python

    from posydon.visualization.VHdiagram import VHdiagram

    VHdiagram('population.h5', path='./dataset/', index=18976)

.. image:: pngs/detailled_window.png

The named parameter 'path' is useful if the dataset is in another directory.

We can choose the view wanted with the option's window, displayed by
'option' button :

.. image:: pngs/option_window.png

The 'save' button take a screen of the view currently displayed, and
save it inside a 'screens' folder, created in the current directory.

Set view
--------

For faster use, it's possible to directly choose the view displayed.
There are 4 available modes:

.. code-block:: python

    PresenterMode.DIAGRAM
    PresenterMode.REDUCED
    PresenterMode.SIMPLIFIED
    PresenterMode.DETAILED

We can specify which view we want to display with the named parameter
'presentMode' :

.. code-block:: python

    from posydon.visualization.VHdiagram import VHdiagram
    from posydon.visualization.VH_diagram.PresenterMode import PresenterMode

    VHdiagram('population.h5', index=19628, presentMode=PresenterMode.DIAGRAM)

.. image:: pngs/diagram_window.png

Set display mode
----------------

There are two ways to display the diagram : inside a window or as a screen
inside a Jupyter notebook. There are three available modes:

.. code-block:: python

    DisplayMode.INLINE_S
    DisplayMode.INLINE_B
    DisplayMode.WINDOW

'INLINE_S' is to do an inline display with matplotlib while 'INLINE_B' use IPython.


We can specify the display mode wanted with the named parameter 'displayMode':

.. code-block:: python

    from posydon.visualization.VHdiagram import VHdiagram, DisplayMode
    from posydon.visualization.VH_diagram.Presenter import PresenterMode

    VHdiagram(
        "population.h5",
        index=19628,
        presentMode=PresenterMode.DIAGRAM,
        displayMode=DisplayMode.INLINE_B,
    )

.. image:: pngs/diagram_inline.png

Visualize multiple indexes
==========================

The `VDdiagramm_m` module allows you to print more than one index horizontally in the
same plot.

Counting binaries populations
-----------------------------

You can loop through the binary file to count identical binary simulations with the `ParseDataFrame`
class. The counts are accessible in the `count_dict` attribute which is a `Counter` python object 
and the frequencies with the `get_frequencies()` method. You can also get the n most frequent
binaries using the `get_most_numpy(k)` method.

.. code-block:: python

    from posydon.visualization.VH_diagram.ParseDataFrame import ParseDataFrame

    parse_df = ParseDataFrame('./data/population.h5')
    parse_df.count_dict

    >>> Counter({0: 14,
         1: 16,
         3: 1,
         4: 9,
         5: 1,
         6: 1,
         7: 6,
         8: 19,
         21: 1,
         24: 4,
         25: 1,
         26: 1,
         30: 1,
         31: 10,
         34: 1,
         36: 5,
         38: 1,
         55: 1,
         58: 1,
         59: 1,
         62: 1,
         65: 1,
         74: 1,
         93: 1,
         95: 1})


Side by side visualization
--------------------------

You can get multiple binary simulations printed side by side by providing a list of their index
in wanted order. The initialization call requires a dict of frequenties as parameter `frequency`
witch is by default uniform. `hierarchy` parameter must be set to false.

.. code-block:: python

    from posydon.visualization.VH_diagram.PresenterMultiple import VHdiagramm_m
    from posydon.visualization.VHdiagram import DisplayMode
    from posydon.visualization.VH_diagram.PresenterMode import PresenterMode

    VHD = VHdiagramm_m('./data/population.h5',
                    index=cnt[:,0],
                    frequency=parse_df.get_frequencies(),
                    hierarchy=False,
                    presentMode = PresenterMode.DIAGRAM,
                    displayMode = DisplayMode.INLINE_B)

.. image:: pngs/diagram_multiple.png

Hierarchical visualization
--------------------------

The hierarchical visualization aims to "factorize" identical steps resulting in a tree plot
where the nodes are the common steps. They are labeled by precentages relatively to the parent
node precentage (witch is also uniform by default).

.. code-block:: python

    from posydon.visualization.VH_diagram.PresenterMultiple import VHdiagramm_m
    from posydon.visualization.VHdiagram import DisplayMode
    from posydon.visualization.VH_diagram.PresenterMode import PresenterMode

    VHD = VHdiagramm_m('./data/population.h5',
                    index=cnt[:,0],
                    frequency=parse_df.get_frequencies(),
                    hierarchy=True,
                    presentMode = PresenterMode.DIAGRAM,
                    displayMode = DisplayMode.INLINE_B)

.. image:: pngs/diagram_hierarchy.png

Steps differentiation is based on picture's file names. So binaries like 30 and 36
are considered different since they refer to distinct files.

.. code-block:: python

    >>> print(VHD._presenter._infos[10][1].event_filename.split("\\")[-1])
    H-rich_oCE1_H-rich_.png
    >>> print(VHD._presenter._infos[14][1].event_filename.split("\\")[-1])
    H-rich_oCE2_H-rich_.png


The `get_sorted_index()` method of the `VHdiagramm_m` object give a list of indexes
sorted regarding the filenames of the pictures representing their steps

.. code-block:: python

    VHD = VHdiagramm_m('./data/population.h5',
                    index=VHD.get_sorted_index(),
                    frequency=parse_df.get_frequencies(),
                    hierarchy=False,
                    presentMode = PresenterMode.DIAGRAM,
                    displayMode = DisplayMode.INLINE_B)

.. image:: pngs/diagram_multiple_sort.png