.. _custom_pop_synth:

###########################
Custom population synthesis
###########################

The modularity of the POSYDON code allows the user to easily change
the default evolutionary tree dictated by the flow chart tree. The flow
chart tree is a list of 4D points composed of ``binary.star_1.state,
binary.star_2.state, binary.state, binary.event`` which map to a corresponding
evolutionary step. Here we show how the user can replace part or the entire
POSYDON flow chart tree and customize the evolution with user specified
steps.

Custom step
===========

The user can develop a ``POSYDON`` step which is a python class with a
``__call__`` method which updates the current properties of the binary objects
(see reference binary object). As an example, let's create our own treatment
for a common envelope. For simplicity sake we will halve the orbital separation
and set the post-CE donor star mass to the star's helium core mass.

.. code-block:: python
  :name: custom_CE_step.py

  from posydon.utils import common_functions as cf

  class my_CE_step(object):
    """Compute a fake CE event."""

    def __init__(self, verbose):
      self.verbose = verbose

    def __call__(self, binary):

      if self.verbose:
        print('The orbital separation post CE is half the pre CE orbital separation!')

      # Determine which star is the donor and which is the companion
      if binary.event in ["oCE1", "oDoubleCE1"]:
          donor_star = binary.star_1
          comp_star = binary.star_2
      elif binary.event in ["oCE2", "oDoubleCE2"]:
          donor_star = binary.star_2
          comp_star = binary.star_1
      else:
          raise ValueError("CEE does not apply if `event` is not "
                           "`oCE1`, 'oDoubleCE1' or `oCE2`, 'oDoubleCE1'")

      binary.separation /= 2.
      m1 = donor_star.he_core_mass
      m2 = comp_star.mass
      binary.separation = cf.orbital_period_from_separation(binary.separation, m1, m2)

The new custom step can be used by importing the new function into the
simulation properties as follows (see ref. the POSYDON API)

.. code-block:: python

  from custom_CE_step import my_CE_step

  CE_STEP = dict(verbose=False)

  # pass the simulation properties to each step
  sim_kwargs = dict(
      flow = (flow_chart, {}),
      step_HMS_HMS = (MS_MS_step, MESA_STEP),
      step_CO_HeMS = (CO_HeMS_step, MESA_STEP),
      step_CO_HMS_RLO = (CO_HMS_RLO_step, MESA_STEP),
      step_detached = (detached_step, DETACHED_STEP),
      step_CE = (my_CE_step, CE_STEP), # <--- WE ONLY EDITED THIS LINE HERE
      step_SN = (StepSN, SN_STEP),
      step_dco = (DoubleCO, DCO_STEP),
      step_end = (step_end, END_STEP),
      extra_hooks = [(StepNamesHooks, {})]
  )


Custom flow charts
==================

Let us now see how one can alter the default flow chart links. Let us assume
that you have a custom step that maps the evolution of Roche lobe overflowing
HMS-HeMS binaries called ``HMS_HeMS_RLO_step``. With the ``CHANGE_FLOW_CHART``
option of the ``flow_chart`` method we can replace the the tree links which
are currently mapped to the ``step_end``
as follows.

.. code-block:: python

  from posydon.binary_evol.flow_chart import flow_chart, POSYDON_FLOW_CHART, STAR_STATES_H_RICH, STAR_STATES_HE_RICH

    NEW_FLOW_CHART_LINKS = {}

    for s1 in STAR_STATES_H_RICH:
        for s2 in STAR_STATES_HE_RICH:
            NEW_FLOW_CHART_LINKS[(s1, s2, 'RLO1', 'oRLO1')] = 'step_HMS_HeMS_RLO'
            NEW_FLOW_CHART_LINKS[(s2, s1, 'RLO2', 'oRLO2')] = 'step_HMS_HeMS_RLO'

    new_flow_chart = flow_chart(FLOW_CHART=POSYDON_FLOW_CHART, CHANGE_FLOW_CHART=NEW_FLOW_CHART_LINKS)

We can load the new flow chart into the simulation properties as.

.. code-block:: python

  from custom_HMS_HeMS_step import HMS_HeMS_RLO_step

  HMS_HeMS_RLO_STEP = dict(verbose=False)

  sim_kwargs = dict(
      flow = (new_flow_chart, {}), # <--- WE ONLY EDITED THIS LINE HERE
      step_HMS_HMS = (MS_MS_step, MESA_STEP),
      step_CO_HeMS = (CO_HeMS_step, MESA_STEP),
      step_CO_HMS_RLO = (CO_HMS_RLO_step, MESA_STEP),
      step_HMS_HeMS_RLO = (HMS_HeMS_RLO_step, HMS_HeMS_RLO_STEP), # ADD NEW STEP HERE
      step_detached = (detached_step, DETACHED_STEP),
      step_CE = (StepCEE, CE_STEP),
      step_SN = (StepSN, SN_STEP),
      step_dco = (DoubleCO, DCO_STEP),
      step_end = (step_end, END_STEP),
      extra_hooks = [(StepNamesHooks, {})]
  )
