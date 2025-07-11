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

Once imported, you can create an instance of POSYDON's default 
flow chart by invoking, the ``flow_chart()`` function, e.g.,

.. code-block:: python

  default_flow = flow_chart()

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
evolved. For example, this ``tuple`` might look like this:

.. code-block:: python

  EXAMPLE_STATE = ('H-rich_Core_H_burning','H-rich_Core_H_burning','detached', 'ZAMS')

which by default maps to ``step_HMS_HMS``. The evolutionary steps that a given 
tuple can match to are described in detail in the :ref:`pop-params-guide`. In 
summary, there are twelve evolutionary steps that POSYDON utilizes. These are:

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

A typical item of a flow chart may then look something like this:

.. code-block:: python

  FLOW_ITEM = {('H-rich_Core_H_burning','H-rich_Core_H_burning','detached', 'ZAMS') : 'step_HMS_HMS'}

POSYDON builds its own default flow chart to handle evolution with 
a complete set of state to evolution step mappings. If you are curious 
to see the default construction, check out ``posydon/binary_evol/flow_chart.py``.

Modifying the Flow Chart
~~~~~~~~~~~~~~~~~~~~~~~~

You can also modify the flow chart to suit your own needs. As a 
flow chart needs to have a mapping for every possible state 
combination ``tuple`` to an evolution step, it can be convenient 
to start by inheriting POSYDON's default flow chart and only 
modify certain state to step mappings that are of interest. 
There are a few ways to do this, one of which is covered in the tutorial on
:ref:`/tutorials-examples/population-synthesis/custom_step_and_flow.ipynb`, 
which also shows how customize an evolution step. For an example of how you might 
write your own importable ``flow_chart`` via :ref:`/components-overview/pop_syn/user_modules.rst`, check out 
``posydon/user_modules/my_flow_chart_example.py``.

However, as a brief example and alternative way to adjust the flow chart, 
let's say we wanted to evolve a single binary with a modified flow. As in 
the tutorial that describes how to :ref:`/tutorials-examples/population-synthesis/evolve_single_binaries.ipynb`, 
we will load in the default simulation properties from ``posydon/pop_syn/population_params_default.ini``. 
First, let's copy the default simulation parameter ``.ini`` file to our current 
directory so we can play around with it:

.. code-block:: python

  import os
  import shutil
  from posydon.config import PATH_TO_POSYDON

  path_to_params = os.path.join(PATH_TO_POSYDON, "posydon/popsyn/population_params_default.ini")
  shutil.copyfile(path_to_params, './population_params.ini')

Next, we can load the evolution steps that our simulation will use, 
providing a desired metallicity to evolve our binary with:

.. code-block:: python

    # load the function to load the simulation properties from the ini file
    from posydon.popsyn.io import simprop_kwargs_from_ini
    from posydon.binary_evol.simulationproperties import SimulationProperties

    # Load the simulation properties from the default ini file.
    sim_kwargs = simprop_kwargs_from_ini('population_params.ini')
    # manually add the metallicity to each step that requires it
    metallicity = {'metallicity':1}

    sim_kwargs['step_HMS_HMS'][1].update(metallicity)
    sim_kwargs['step_CO_HeMS'][1].update(metallicity)
    sim_kwargs['step_CO_HMS_RLO'][1].update(metallicity)
    sim_kwargs['step_CO_HeMS_RLO'][1].update(metallicity)
    sim_kwargs['step_detached'][1].update(metallicity)
    sim_kwargs['step_disrupted'][1].update(metallicity)
    sim_kwargs['step_merged'][1].update(metallicity)
    sim_kwargs['step_initially_single'][1].update(metallicity)

    sim_prop = SimulationProperties(**sim_kwargs)
    # Load the steps and required data
    sim_prop.load_steps(verbose=True)

At this point, we would be ready to simulate our binary at solar metallicity, 
however, let's modify the flow first. The ``population_params.ini`` loads the 
default POSYDON flow chart by default. You can access it via 

.. code-block:: python

  # access the loaded flow_chart dictionary
  sim_prop.flow

If you print this, you will see the dictionary that is the default 
POSYDON flow chart. Let's replace this flow chart with our own.

.. code-block:: python
   
  # import POSYDON default flow chart (we're going to modify it)
  from posydon.binary_evol.flow_chart import flow_chart

  # create and store a copy of the default flow chart
  my_flow = flow_chart()

There are two optional arguments in the ``flow_chart`` function 
that facilitate modifications. 

One is the ``FLOW_CHART`` argument, which is ``None`` by default. 
With this argument, you may provide an entirely new flow chart 
dictionary that you have constructed, overriding POSYDON's default. 
If this argument is left as ``None``, then POSYDON simply loads its 
own flow chart.

The other argument is ``CHANGE_FLOW_CHART``, which also expects a 
dictionary as input. The difference is that this argument expects 
only a subset of the dictionary items that comprise a flow chart. 
Therefore, this argument is ideal for modifying only specific 
state to evolution step mappings within a pre-exisiting flow chart. 

We can use the ``CHANGE_FLOW_CHART`` argument like so to overwrite 
the HMS-HMS, detached ZAMS state (which is often the first step of 
evolution that a binary will go through). Normally, this state is 
mapped to ``step_HMS_HMS`` (as in the ``EXAMPLE_ITEM`` above), 
which will be evolved using the MESA grids. Let's set this to map 
to ``step_detached`` instead. You can do that like this:

.. code-block:: python

  custom_mappings = {('H-rich_Core_H_burning','H-rich_Core_H_burning','detached', 'ZAMS') : 'step_detached'}
  my_flow = flow_chart(CHANGE_FLOW_CHART=custom_mappings)

This state will now be mapped to ``step_detached``, instead of the default 
``step_HMS_HMS``. Likewise, you can replace any state-evolution step 
mapping of your choosing. As a last step in this example, let's replace 
the ``flow_chart`` that is stored in our simulation properties:

.. code-block:: python
  
  sim_prop.flow = my_flow

Now, if you proceed to simulate binary star evolution, your flow chart 
will be used in place of POSYDON's default.



