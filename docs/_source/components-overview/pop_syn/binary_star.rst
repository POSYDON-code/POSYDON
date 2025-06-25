.. _binary-star:


The BinaryStar object
======================

The ``BinaryStar`` object is composed of two ``SingleStar`` objects 
(see :ref:`single-star`) and contains the current and past states of 
the binary. Only parameters in the ``BINARYPROPERTIES`` list are stored in the 
history. The current parameter value of the star object is accessed with, e.g. 
``binary.orbital_period`` and the past history via 
``binary.orbital_period_history``. The two stars are accessed with, e.g. 
(for star 1), ``binary.star_1.mass`` and the past history via 
``binary.star_1.mass_history``.

To use BinaryStar object import it using:

.. code-block:: python

  from posydon.binary_evol.singlestar import SingleStar
  from posydon.binary_evol.binarystar import BinaryStar


Creating a BinaryStar object
----------------------------

BINARYPROPERTIES
~~~~~~~~~~~~~~~~

The binary properties are defined as follows

.. list-table:: BINARYPROPERTIES
  :header-rows: 1
  :widths: 50 150

  * - Properties
    - Descriptions
  * - ``state``
    - The state of the binary, see state options.
  * - ``event``
    - The event of the binary, see event options.
  * - ``time``
    - Age of the binary system in yr.
  * - ``separation``
    - Orbital separation in R_sun.
  * - ``orbital_period``
    - Orbital period in days.
  * - ``eccentricity``
    - Orbital eccentricity.
  * - ``V_sys``
    - Velocity of the center of mass of the binary [Vx, Vy, Vz] in km/s.
  * - ``mass_transfer_case``
    - Mass transfer case, see MT case options.
  * - ``lg_mtransfer_rate``
    - The logarithm of the mass transfer rate in Msun/yr.
  * - ``step_names``
    - Names of the steps in the evolution.
  * - ``step_times``
    - Time spent in the steps in the evolution.
  * - ``rl_relative_overflow_1``
    - The relative overflow of the Roche Lobe of star 1.
  * - ``rl_relative_overflow_2``
    - The relative overflow of the Roche Lobe of star 2.
  * - ``trap_radius``
    - The trapping radius of the binary in R_sun.
  * - ``acc_radius``
    - The accretion radius of the binary in R_sun.
  * - ``t_sync_rad_1``
    - The synchronization time of the radiative zone of star 1 in yr.
  * - ``t_sync_conv_1``
    - The synchronization time of the convective zones of star 1 in yr.
  * - ``t_sync_rad_2``
    - The synchronization time of the radiative zone of star 2 in yr.
  * - ``t_sync_conv_2``
    - The synchronization time of the convective zones of star 2 in yr.
  * - ``nearest_neighbour_distance``
    - The distance to the nearest neighbour for NN interpolation


Additional scalar properties can be added during the evolution.

Since they do not change over time, they are not stored in the history.
These can be requested and will be stored in the output oneline (See the 
:ref:`Synthetic Population<synthetic-population>` and 
:ref:`Population Parameter Guide<pop-params-guide>` for more information).

Additional columns
~~~~~~~~~~~~~~~~~~

The additional columns are defined as follows:

.. list-table:: Additional columns
  :header-rows: 1
  :widths: 50 150

  * - Properties
    - Descriptions
  * - ``interp_class_HMS_HMS``
    - The interpolation class for the HMS_HMS phase.
  * - ``interp_class_CO_HMS_RLO``
    - The interpolation class for the CO_HMS_RLO phase.
  * - ``interp_class_CO_HeMS``
    - The interpolation class for the CO_HeMS phase.
  * - ``interp_class_CO_HeMS_RLO``
    - The interpolation class for the CO_HeMS_RLO phase.
  * - ``mt_history_HMS_HMS``
    -  The mass transfer history for the HMS_HMS phase.
  * - ``mt_history_CO_HMS_RLO``
    - The mass transfer history for the CO_HMS_RLO phase.
  * - ``mt_history_CO_HeMS``
    - The mass transfer history for the CO_HeMS phase.
  * - ``mt_history_CO_HeMS_RLO``
    - The mass transfer history for the CO_HeMS_RLO phase.


State options
~~~~~~~~~~~~~

Binary states are defined according to the following table:

.. list-table:: States
  :header-rows: 1
  :widths: 10 30

  * - State
    - Description
  * - ``initially_single_star``
    - The binary was initially a single star.
  * - ``detached``
    - The stars in the binary are in a detached state.
  * - ``RLO1``
    - The binary is Roche Lobe overflowing, star 1 is overfilling the RL.
  * - ``RLO2``
    - The binary is Roche Lobe overflowing, star 2 is overfilling the RL.
  * - ``contact``
    - The stars in the binary are in contact.
  * - ``disrupted``
    - The binary was disrupted.
  * - ``merged``
    - The stars in the binary merged.
  * - ``initial_RLOF``
    - The binary is in the initial Roche Lobe overflow.
  * - ``maxtime``
    - Max time of the evolution was reached.
  * - ``FAILED``
    - The evolution failed.


Event options
~~~~~~~~~~~~~

Binary events are defined according to the following table:

.. list-table:: Events
  :header-rows: 1
  :widths: 10 30

  * - Event
    - Description
  * - ``ZAMS``
    - Zero Age Main Sequence
  * - ``CC1``
    - Core collapse of star 1.
  * - ``CC2``
    - Core collapse of star 2.
  * - ``oRLO1``
    - The binary is at onset of Roche Lobe overflow, star 1 is overfilling the RL.
  * - ``oRLO2``
    - The binary is at onset of Roche Lobe overflow, star 2 is overfilling the RL.
  * - ``oCE1``
    - The binary is at the onset of Common Envelope initiated by star 1.
  * - ``oCE2``
    - The binary is at the onset of Common Envelope initiated by star 2.
  * - ``oDoubleCE1``
    - | The binary is at the onset of Double Common Envelope initiated by star 1. 
      | Both stars are post main-sequence.
  * - ``oDoubleCE2``
    - | The binary is at the onset of Double Common Envelope initiated by star 2. 
      | Both stars are post main-sequence.
  * - ``CO_contact``
    - The binary reached contact in the compact object phase.
  * - ``redirect_from_ZAMS``
    - | The binary was redirected from ZAMS for a variety of reasons.
      | Only recorded if history_verbose = True
  * - ``redirect_from_CO_HMS_RLO``
    - | The binary was redirected from CO_HMS_RLO for a variety of reasons.
      | Only recorded if history_verbose = True
  * - ``redirect_from_CO_HeMS``
    - | The binary was redirected from CO_HeMS for a variety of reasons.
      | Only recorded if history_verbose = True
  * - ``redirect_from_CO_HeMS_RLO``
    - | The binary was redirected from CO_HeMS_RLO for a variety of reasons.
      | Only recorded if `history_verbose = True`
  * - ``MaxTime_exceeded``
    - The maximum time of the evolution was exceeded.
  * - ``maxtime``
    - The maximum time of the evolution was reached.
  * - ``oMerging1``
    - The binary is at the onset of merging, star 1 is overfilling the RL.
  * - ``oMerging2``
    - The binary is at the onset of merging, star 2 is overfilling the RL.
  * - ``None``
    - No event occurred.
  * - ``ERR``
    - An error occurred in the evolution.
  * - ``END``
    - The binary evolution was stopped.


Mass Transfer case
~~~~~~~~~~~~~~~~~~

The mass transfer cases are stored in `mt_history_GRIDTYPE` and are defined 
according to the following table: TODO: add the table below

.. list-table:: Mass transfer cases
  :header-rows: 1
  :widths: 10 30

  * - Case
    - Description
  * - ``None``
    - The binary is not Roche Lobe overflowing.


TODO: update properties


Basic example
~~~~~~~~~~~~~

The simplest method is to provide the two star objects and `kwargs` of the 
initial binary parameters.

.. code-block:: python

  from posydon.utils.constants import Zsun

  kwargs1 = {'state' : 'H-rich_Core_H_burning',
             'mass' : 20.0,
             'metallicity' : Zsun}

  star_1 = SingleStar(**kwargs1)

  kwargs2 = {'state' : 'H-rich_Core_H_burning',
             'mass' : 10.0,
             'metallicity' : Zsun}

  star_2 = SingleStar(**kwargs2)

  kwargs3 = {'state' : 'detached',
             'event' : None,
             'time' : 0.,
             'orbital_period' : 3.,
             'eccentricity' : 0.}

  binary = BinaryStar(star_1, star_2, **kwargs3)
