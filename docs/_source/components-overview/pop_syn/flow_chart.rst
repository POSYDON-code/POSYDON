.. _flow-chart:

The Flow Chart Object
---------------------

Along with the ``BinaryStar`` object the ``flow_chart`` object is one of the most 
important building blocks of POSYDON, such that it is featured in the 
center of the POSYDON logo. The ``flow_chart`` object maps the evolution of a 
``BinaryStar`` object to its corresponding evolutionary step.

To use the ``flow_chart`` object, import it using:

.. code-block:: python

    from posydon.binary_evol.flow_chart import flow_chart

Components of a Flow Chart
~~~~~~~~~~~~~~~~~~~~~~~~~~
At its heart, the ``flow_chart`` object is simply a dictionary. Its keys are 
tuples that describe the state of a binary. Its values are the corresponding 
evolutionary steps that these states map to.

The tuples that define a binary state are expected to conform to a particular 
format. Schematically, they should look like this:

.. code-block:: python

    EXAMPLE_STATE = (S1_state, S2_state, binary_state, binary_event)

Here, ``S1_state`` or ``S2_state`` is a valid ``SingleStar.state`` property of either 
star 1 or star 2 (see the State Options of :ref:`single-star`), 
``binary_state`` should be a valid ``BinaryStar.state`` property, and 
``binary_event`` a valid ``BinaryStar.event`` property (see the State 
and Event Options of :ref:`binary-star`). A tuple conforming to the 
convention above fully describes the evolutionary state of a binary 
star system as far as POSYDON is concerned. This may then be mapped 
to a particular evolution step that dictates how that state is 
evolved.

The evolutionary steps that a given tuple can match to are described 
in detail in the :ref:`pop-params-guide`. In summary, there are 
twelve evolutionary steps that POSYDON utilizes. These are:

.. list-table:: Evolution steps
  :header-rows: 1
  :widths: 50 150

  * - Step name
    - Description
  * - ``step_HMS_HMS``
    - Evolve the binary using a MESA grid of HMS-HMS binaries.
  * - ``step_CO_HeMS``
    - Evolve the binary using a MESA grid of CO-HeMS binaries.
  * - ``step_CO_HMS_RLO``
    - Evolve the binary using a MESA grid of CO-HMS binaries, starting from RLO.
  * - ``step_CO_HeMS_RLO``
    - Evolve the binary using a MESA grid of CO-HeMS binaries, starting from RLO.
  * - ``step_detached``
    - Evolve the binary through detached evolution, using MESA single star grids.
  * - ``step_disrupted``
    - Evolve the binary after disruption via a supernova; evolving the remaining star with MESA single star grids.
  * - ``step_merged``
    - Evolve the binary in the event that its components merge, forming a single stellar object.
  * - ``step_initially_single``
    - Evolve an initially single star using the MESA single star grids.
  * - ``step_dco``
    - Evolve the binary through the case where both components are compact objects.
  * - ``step_SN``
    - Evolve the binary through a supernova event.
  * - ``step_CE``
    - Evolve the binary through a common envelope event.
  * - ``step_end``
    - This represents a state where evolution comes to an end.

POSYDON builds its own default flow chart to handle evolution, but you can modify it to suit your 
own needs too.

Modifying the Flow Chart
~~~~~~~~~~~~~~~~~~~~~~~~






