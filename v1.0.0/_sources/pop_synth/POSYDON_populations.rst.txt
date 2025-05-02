#################################
Running populations using POSYDON
#################################

1. Run a basic population
=========================

To run the most basic population with all the pre-set POSYDON defaults,
you only need the following few lines of code. By default the code
runs 100 initial binaries of a constant star formation for a duration of
Hubble time (13.7e9 years).

All the required POSYDON imports:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    from posydon.popsyn.binarypopulation import BinaryPopulation
    from posydon.binary_evol.simulationproperties import SimulationProperties
    from posydon.binary_evol.flow_chart import flow_chart
    from posydon.binary_evol.CE.step_CEE import StepCEE
    from posydon.binary_evol.SN.step_SN import StepSN
    from posydon.binary_evol.step_end import step_end
    from posydon.binary_evol.MESA.step_mesa import CO_HeMS_step, MS_MS_step, CO_HMS_RLO_step
    from posydon.binary_evol.DT.step_detached import detached_step
    from posydon.binary_evol.DT.double_CO import DoubleCO

Define the POSYDON evolutionary steps
~~~~~~~~~~~~~~~~~~~~~~~~

We should define all the evolutionary steps involved. By default, the parameters
for these steps are pre-defined. Then, ``SimulationProperties(**sim_kwargs)``
initializes the population object with the steps defined, in ``sim_prop``. This
initialization step can take a few seconds as it requires loading large data
files.

.. code:: ipython3

    sim_kwargs = dict(
        flow = (flow_chart, {}),
        step_HMS_HMS = (MS_MS_step, {}),
        step_CO_HeMS = (CO_HeMS_step, {}),
        step_CO_HMS_RLO = (CO_HMS_RLO_step, {}),
        step_detached = (detached_step, {}),
        step_CE = (StepCEE, {}),
        step_SN = (StepSN, {}),
        step_dco = (DoubleCO, {}),
        step_end = (step_end, {})
    )

    sim_prop = SimulationProperties(**sim_kwargs)

Initializing the population object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We invoke the population class ``BinaryPopulation``, which is the
population class with the function ``evolve`` that does the work of
evolving the binaries. The population class takes in the simulation
properties to initialize the population.

.. code:: ipython3

    pop = BinaryPopulation(population_properties=sim_prop)

    pop.evolve()

This evolved population can be converted to a pandas dataframe. The
dataframe helps us in studying the evolved population by utilizing the
benefits of the functionalities available for managing dataframes.

.. code:: ipython3

    DF = pop.to_df()

2. Customize the script to run more complicated populations
===========================================================

Now, for most science cases, the default POSYDON parameters would need
to be changed. Below we show an example of how to run a population of
1000 initial binaries with constant star-formation for a duration of
Hubble time.

Define the POSYDON steps and their parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If there is a need to change a parameter in the POSYDON steps that would
be done in the simulation properties. For example, to change the
``interpolation_method`` for the steps that evolve the binaries
according to the MESA grids, we introduce a new string
``interp_str = 'nearest_neighbour'`` and add this string to the
parameter ``interpolation_method`` for the steps ``step_HMS_HMS``,
``step_CO_HeMS``, and ``step_CO_HMS_RLO``.

.. code:: ipython3


    interp_str = 'nearest_neighbour'

    sim_kwargs = dict(
        flow = (flow_chart, {}),
        step_HMS_HMS = (MS_MS_step, dict(interpolation_method=interp_str)),
        step_CO_HeMS = (CO_HeMS_step, dict(interpolation_method=interp_str)),
        step_CO_HMS_RLO = (CO_HMS_RLO_step, dict(interpolation_method=interp_str)),
        step_detached = (detached_step, {}),
        step_CE = (StepCEE, {}),
        step_SN = (StepSN, {}),
        step_dco = (DoubleCO, {}),
        step_end = (step_end, {}),
    )

    sim_prop = SimulationProperties(**sim_kwargs)


Define the properties of the initial generated population
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If some parameters of the initial population are required to be changed,
they are changed as follows in ``kwargs``. Here we change the
``number_of_binaries``.

.. code:: ipython3

    kwargs = {'number_of_binaries' : 1000
             }

Initializing the population object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As before, we invoke the population class, this time with the parameters
of the initial population that we want to include in ``kwargs``.

.. code:: ipython3

    pop = BinaryPopulation(population_properties=sim_prop,
                           **kwargs)

Evolve it! And covert to pandas dataframe
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The boolean parameter ``tqdm`` shows the progress bar.

.. code:: ipython3

    pop.evolve(breakdown_to_df=False, tqdm=True)
    DF = pop.to_df()

To save the population as an HDF5 file use the ``save`` function.

.. code:: ipython3

    pop.save('population_1000.h5', **kwargs )

3. Identify double compact objects
==================================

Now, the evolved population has the complete evolution of all binaries. We
only require the last row from all the binaries as we want to study the
“current” properties of the double compact objects. The example is to
search for binaries not ``disrupted`` in the end, in which ``star_1`` is
a BH and ``star_2`` is an NS.

.. code:: ipython3

    output_cols = ['state','time','event','S1_state','S2_state','S1_mass','S2_mass','orbital_period','eccentricity']
    DF.loc[(DF['S1_state'] == 'BH')&(DF['S2_state'] == 'NS')&(DF['event'] == 'END')&(DF['state'] != 'disrupted')
           ,output_cols]

To look at the evolution of one individual binary, choose proper
quantities to display:

.. code:: ipython3

    index_BHNS = DF.loc[(DF['S1_state'] == 'BH')&(DF['S2_state'] == 'NS')&
                        (DF['event'] == 'END')&
                        (DF['state'] != 'disrupted')].index
    DF.loc[index_BHNS[0], output_cols]

To identify the BHNS binaries that specifically went through a common envelope
phase with a BH and a donor:

.. code:: ipython3

    DF_BHNS = DF.loc[index_BHNS]
    index_BHNS_CE = DF_BHNS.loc[(DF_BHNS["event"]=="oCE2")].index

.. code:: ipython3

    DF_BHNS.loc[index_BHNS_CE[0], output_cols]

4. Identify the X-ray Binaries
==============================

In order to investigate the XRBs from this population, we need to look
into the originally run population again. Again, we only require the
last row from all the binaries as we want to study the “current”
properties of the XRBs. Also, we want to pick out binaries where only
one of the two stars at the end is a compact object (NS or BH), and the
other star is definitively not a compact object (not even a white
dwarf).

.. code:: ipython3


    DF_end = DF[DF['event']=='END']
    star1_is_CO = ( (DF_end["S1_state"] == "NS") | (DF_end["S1_state"] == "BH") ) & (DF_end["S2_state"] != "WD")
    star2_is_CO = ( (DF_end["S2_state"] == "NS") | (DF_end["S2_state"] == "BH") ) & (DF_end["S2_state"] != "WD")
    exactly_one_CO = (star1_is_CO & ~star2_is_CO) | (~star1_is_CO & star2_is_CO)
    only_one_CO = DF_end[exactly_one_CO]


Looking at the final ``state`` of all the binaries, we can figure out
the ones we want to keep. For XRBs, we are interested in states
``detached`` (corresponding to wind accretion) and ``RLO1``/``RLO2``
(corresponding to Roche-lobe overflow).

.. code:: ipython3

    only_one_CO['state'].unique()

Most of our XRB RLO states would be ``RLO2`` (mass-transferred from
``star_2`` to ``star_1``), because the primary (``star_1``) is always
the initially more massive star which would evolve first and form a
compact object. However, there might a rare case where ``star_1``
transferred mass to ``star_2`` making it more massive than ``star_1``,
leading to a compact object for ``star_2``. Therefore, to be sure, we
look for both ``RLO1`` (mass-transferred from ``star_1`` to ``star_2``)
and ``RLO2``.

.. code:: ipython3

    DF_XRB = only_one_CO[(only_one_CO['state']=='detached') |
                         (only_one_CO['state']=='RLO1') |
                         (only_one_CO['state']=='RLO2')]

Now we have our subset of binaries from the total population. To see
them as XRBs, we need some information about their X-ray luminosities.
We can use some simple functions to calculate the corresponding
accretion efficiencies and luminosities.

Importing some packages needed for these calculations.

.. code:: ipython3

    from posydon.utils.common_functions import CO_radius
    import posydon.utils.constants as const
    import numpy as np

Since the accretor can be ``star_1`` or ``star_2``, we need to do a
check to see which star is the compact object. We calculate the
accretion luminosity (``Lacc``) using the relation
``eta * mass_accretion_rate * light_speed^2``, where ``eta`` is the efficiency
at which gravitational potential mass is converted to luminosity.


.. code:: ipython3


    def Lx(S1_state, S1_mass, S1_lg_mdot,S2_state, S2_mass, S2_lg_mdot, lg_mtransfer_rate):

        if np.isnan(lg_mtransfer_rate):
            if S1_state == 'NS':
                eta = 0.1
                acc_lg_mdot = S1_lg_mdot
            elif S1_state == 'BH':
                eta = 0.057
                acc_lg_mdot = S1_lg_mdot
            else:
                print("CO??", S1_state, S2_state)
        else:
            acc_lg_mdot = lg_mtransfer_rate
            if S1_state == 'NS':
                eta = 0.1
            elif S1_state == 'BH':
                eta = 0.057
            else:
                print("CO??", S1_state, S2_state)

        return eta * (10.0**acc_lg_mdot * const.Msun / const.secyer) * const.clight**2


We define a new column ``Lacc`` which stands for the accretion
luminosity and calculate it using the function defined above.

.. code:: ipython3


    DF_XRB['Lacc']=DF_XRB[['S1_state', 'S1_mass', 'S1_lg_mdot',
                     'S2_state', 'S2_mass','S2_lg_mdot',
                     'lg_mtransfer_rate']].apply(lambda x: Lx(x['S1_state'], x['S1_mass'],
                                                              x['S1_lg_mdot'],x['S2_state'],
                                                              x['S2_mass'], x['S2_lg_mdot'],
                                                              x['lg_mtransfer_rate']), axis=1)


Let’s look at the properties of the XRBs, we can use
``pd.set_option('display.max_columns', None)`` to expand the columns.

.. code:: ipython3

    import pandas as pd

    #pd.set_option('display.max_columns', None)

    DF_XRB[['state','time','event','S1_state','S2_state','S1_mass','S2_mass','orbital_period']]
