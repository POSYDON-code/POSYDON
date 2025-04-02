.. _custom_hooks:

##################
Evolutionary Hooks
##################

During the evolution of a binary, one may be interested in having access to
the ``BinaryStar`` instance before, during, and after the evolutionary loop.
You may check the attributes, or add your own data to the binary during the
evolution.
We use *hooks* to pass functions that are called during the evolution
of all binaries.

The evolutionary loop is in ``BinaryStar.evolve``, and the pseudo code looks
like this:

.. code-block:: python

  def evolve(*args, **kwargs):
      pre_evolve(binary)
          while binary.event != 'END':
          pre_step(binary, step_name)
          ...
          binary.evolve()
          ...
          post_step(binary, step_name)
          ...
      post_evolve(binary)


All hooks functions and classes should be provided in the
``SimulationProperties`` with the name ``extra_hooks``.
Then ``extra_hooks`` points to a ``list`` of ``tuples``, containing either:

#. The corresponding extra step name and a callable or,

#. a class deriving from ``EvolveHooks`` and a dictionary passed during initialization.


1. Passing extra hooks functions
================================

.. code-block:: python
  :linenos:
  :emphasize-lines: 12,13,15,16,29,30

  # import evolutionary steps
  from posydon.popsyn.binarypopulation import BinaryPopulation
  from posydon.binary_evol.simulationproperties import SimulationProperties
  from posydon.binary_evol.flow_chart import flow_chart
  from posydon.binary_evol.CE.step_CEE import StepCEE
  from posydon.binary_evol.SN.step_SN import StepSN
  from posydon.binary_evol.step_end import step_end
  from posydon.binary_evol.MESA.step_mesa import CO_HeMS_step, MS_MS_step, CO_HMS_RLO_step
  from posydon.binary_evol.DT.step_detached import detached_step
  from posydon.binary_evol.DT.double_CO import DoubleCO

  def print_info_before(binary, step_name):
      print( f'Step {step_name} for \nbefore:\t{binary}' )

  def print_info_after(binary, step_name):
      print( f'after:\t{binary}\n' )

  sim_kwargs = dict(
        flow = (flow_chart, {}),
        step_HMS_HMS = (MS_MS_step, {}),
        step_CO_HeMS = (CO_HeMS_step, {}),
        step_CO_HMS_RLO = (CO_HMS_RLO_step, {}),
        step_detached = (detached_step, {}),
        step_CE = (StepCEE, {}),
        step_SN = (StepSN, {}),
        step_dco = (DoubleCO, {}),
        step_end = (step_end, {}),
        # put hooks here
        extra_hooks = [('extra_pre_step', print_info_before),
                      ('extra_post_step', print_info_after)]
    )

    sim_prop = SimulationProperties(**sim_kwargs)

With the corresponding output from evolving a single binary:

.. code-block:: text

    Step step_HMS_HMS for
    before:	BinaryStar(detached, ZAMS, p=3.13, S1=(H-rich_Core_H_burning,M=58.14), S2=(H-rich_Core_H_burning,M=50.83))
    after:	BinaryStar(contact, oCE1, p=4.05, S1=(H-rich_Core_H_burning,M=42.97), S2=(H-rich_Core_H_burning,M=44.23))

    Step step_CE for
    before:	BinaryStar(contact, oCE1, p=4.05, S1=(H-rich_Core_H_burning,M=42.97), S2=(H-rich_Core_H_burning,M=44.23))
    after:	BinaryStar(merged, None, p=nan, S1=(H-rich_Core_H_burning,M=87.20), S2=(H-rich_Core_H_burning,M=nan))

    Step step_end for
    before:	BinaryStar(merged, None, p=nan, S1=(H-rich_Core_H_burning,M=87.20), S2=(H-rich_Core_H_burning,M=nan))
    after:	BinaryStar(merged, END, p=nan, S1=(H-rich_Core_H_burning,M=87.20), S2=(H-rich_Core_H_burning,M=nan))

The options for step names and their arguments include:

* ``'extra_pre_evolve'``, functional form: ``f(binary)``
* ``'extra_pre_step'``, functional form: ``f(binary, step_name)``
* ``'extra_post_step'``, functional form: ``f(binary, step_name)``
* ``'extra_post_evolve'``, functional form: ``f(binary)``


2. Passing extra hooks classes
==============================

For more complex tasks, you can provide a class deriving from ``EvolveHooks``.
We have implemented two examples of hooks classes which add step timing and
step names to the history of a binary in ``simulationproperties.py``:

.. code-block:: python
  :linenos:
  :emphasize-lines: 13,26

  # import evolutionary steps
  from posydon.popsyn.binarypopulation import BinaryPopulation
  from posydon.binary_evol.simulationproperties import SimulationProperties
  from posydon.binary_evol.flow_chart import flow_chart
  from posydon.binary_evol.CE.step_CEE import StepCEE
  from posydon.binary_evol.SN.step_SN import StepSN
  from posydon.binary_evol.step_end import step_end
  from posydon.binary_evol.MESA.step_mesa import CO_HeMS_step, MS_MS_step, CO_HMS_RLO_step
  from posydon.binary_evol.DT.step_detached import detached_step
  from posydon.binary_evol.DT.double_CO import DoubleCO

  # import hooks classes
  from posydon.binary_evol.simulationproperties import TimingHooks, StepNamesHooks

  sim_kwargs = dict(
        flow = (flow_chart, {}),
        step_HMS_HMS = (MS_MS_step, {}),
        step_CO_HeMS = (CO_HeMS_step, {}),
        step_CO_HMS_RLO = (CO_HMS_RLO_step, {}),
        step_detached = (detached_step, {}),
        step_CE = (StepCEE, {}),
        step_SN = (StepSN, {}),
        step_dco = (DoubleCO, {}),
        step_end = (step_end, {}),
        # put hooks here
        extra_hooks = [(TimingHooks, {}),(StepNamesHooks, {})]
    )

    sim_prop = SimulationProperties(**sim_kwargs)


The source for ``StepNamesHooks``:

.. code-block:: python

    class StepNamesHooks(EvolveHooks):
        """Add history column 'step_name' to each binary.

        Name of evolutionary step as defined in SimulationProperties.

        >>> pop.to_df(extra_columns=['step_names'])
        """

        def pre_evolve(self, binary):
            """Initialize the step name to match history."""
            binary.step_names = ['initial_cond']
            return binary

        def pre_step(self, binary, step_name):
            """Do not do anything before the step."""
            return binary

        def post_step(self, binary, step_name):
            """Record the step name."""
            binary.step_names.append(step_name)
            len_binary_hist = len(binary.event_history)
            len_step_names = len(binary.step_names)
            diff = len_binary_hist - len_step_names
            if len_binary_hist > len_step_names:
                binary.step_names += [None] * (diff)
            elif len_binary_hist < len_step_names:
                binary.step_names = binary.step_names[-(len_binary_hist - 1):]
            return binary

        def post_evolve(self, binary):
            """Ensure None's are append to step_names to match rows in history."""
            if binary.event == 'END' or binary.event == 'FAILED':
                diff = int(len(binary.event_history) - len(binary.step_names))
                binary.step_names += [None]*diff
            return binary


You can simply combine multiple hooks classes and functions together by passing them
in the list.

.. code-block:: python

    sim_kwargs = dict(
        ...
        extra_hooks = [(TimingHooks, {}),(StepNamesHooks, {}),
                      ('extra_pre_step', print_info_before)],
        ...
    )
