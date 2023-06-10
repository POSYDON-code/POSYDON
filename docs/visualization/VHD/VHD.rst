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

The `VDdiagramm_m` module allow you to print more than one indexes horizontally in the
same plot.

Counting binaries populations
-----------------------------

You can browse throught binaries file to get the count of all the distinct binaries occurences
with the 'ParseDataFrame' class. Counts can be accessed in the 'count_dict' witch is a python 
`Counter` object, and frequencies with the 'get_frequencies()' method.(Parsing may exceed
10 s/M lines for fixed format Dataframes).You can also get the n most frequent binaries by
using the 'get_most_numpy(k)' method.

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

Here you can simply provide a list of index instead of a single one and everything
act as with the previous class. Binaries are plot by their order in the provided list.
The class initialization requires a dict of frequenties distribution witch is by 
default generated as an uniform distribution. (Open the picture in a new tab ofyour
browser for a better visualization)

.. code-block:: python

    from posydon.visualization.VH_diagram.PresenterMultiple import VHdiagramm_m
    from posydon.visualization.VHdiagram import DisplayMode
    from posydon.visualization.VH_diagram.PresenterMode import PresenterMode

    cnt = parse_df.get_most_numpy(len(parse_df.count_dict))

    VHD = VHdiagramm_m('./data/population.h5',
                    index=cnt[:,0],
                    count_dict=parse_df.get_frequencies(),
                    hierarchy=False,
                    presentMode = PresenterMode.DIAGRAM,
                    displayMode = DisplayMode.INLINE_B)

.. image:: pngs/diagram_multiple.png

Hierarchical visualization
--------------------------

The hierarchical visualization tends to "factorize" identical steps resulting in a tree 
where the nodes are the common steps. The nodes are labeled by precentages relatives to the
parent node precentage (witch is also uniform by default).

.. code-block:: python

    from posydon.visualization.VH_diagram.PresenterMultiple import VHdiagramm_m
    from posydon.visualization.VHdiagram import DisplayMode
    from posydon.visualization.VH_diagram.PresenterMode import PresenterMode

    VHD = VHdiagramm_m('./data/population.h5',
                    index=cnt[:,0],
                    count_dict=parse_df.get_frequencies(),
                    hierarchy=True,
                    presentMode = PresenterMode.DIAGRAM,
                    displayMode = DisplayMode.INLINE_B)

.. image:: pngs/diagram_hierarchy.png

Note that binaries 36 and 30 are considered distinct since they refer to differents simulation
representation image filenames. Distincts states or event may look the same for distinct png's
haven't been found yet for all possible pictures. This is set to change with future updates 
where every state / event png representation should be different.

.. code-block:: python

    >>> print(VHD._presenter._infos[10][1].event_filename.split("\\")[-1])
    H-rich_oCE1_H-rich_.png
    >>> print(VHD._presenter._infos[14][1].event_filename.split("\\")[-1])
    H-rich_oCE2_H-rich_.png


Also You can see the effect of that factorisation with the representaion before the elimination
form the index list of binaries sorted hierarchycally with the `get_sorted_index()` method of the
`VHdiagramm_m` object

.. code-block:: python

    VHD = VHdiagramm_m('./data/population.h5',
                   index=VHD.get_sorted_index(),
                   count_dict=parse_df.get_frequencies(),
                   hierarchy=False,
                   presentMode = PresenterMode.DIAGRAM,
                  displayMode = DisplayMode.INLINE_B)

.. image:: pngs/diagram_multiple_sort.png