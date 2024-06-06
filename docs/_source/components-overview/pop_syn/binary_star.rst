.. _binary-star:


The BinaryStar object
======================

The BinaryStar object is composed of two SingleStar objects and contains the current and past states of the binary. Only parameters in the BINARYPROPERTIES list are stored in the history. The current parameter value of the star object is accessed as, e.g. `binary.orbital_period` while his past history with `binary.orbital_period_history`. The two stars are accesses as, e.g. `binary.star_1.mass` while his past history with `binary.star_1.mass_history`.

To use BinaryStar object import it using:

.. code-block:: python

  from posydon.binary_evol.singlestar import SingleStar
  from posydon.binary_evol.binarystar import BinaryStar


Creating a BinaryStar object
----------------------------

BINARYPROPERTIES
~~~~~~~~~~~~~~~~

The binary properties are defined as follows

.. csv-table:: BINARYPROPERTIES
   :header: "Properties", "Descriptions"
   :widths: 50, 150

   `state`, "The state of the binary, see state options."
   `event`, "The event of the binary, see event options."
   `time`, "Age of the binary system in yr."
   `separation`, "Orbital separation in R_sun."
   `orbital_period`, "Orbital period in days."
   `eccentricity`, "Orbital eccentricity."
   `V_sys`, "Velocity of the centre of mass of the binary [Vx, Vy, Vz] in km/s."
   `mass_transfer_case`, "Mass transfer case, see MT case options."
   `lg_mtransfer_rate`, "the logarithm of the mass transfer rate in Msun/yr."
   `step_names`, "Names of the steps in the evolution."
   `step_times`, 'Time spend in the steps in the evolution.'
   `rl_relative_overflow_1`, "The relative overflow of the Roche Lobe of star 1."
   `rl_relative_overflow_2`, "The relative overflow of the Roche Lobe of star 2."
   `trap_radius`, "The trapping radius of the binary in R_sun."
   `acc_radius`, "The accretion radius of the binary in R_sun."
   `t_sync_rad_1`, "?"
   `t_sync_conv_1`, "?"
   `t_sync_rad_2`, "?"
   `t_sync_conv_2`, "?"
   `nearest_neighbour_distance`, "The distance to the nearest neighbour for NN interpolation."


Additional scalar properties can be added during the evolution.

Since they do not change over time, they are not stored in the history.
These can requested and will be stored in the output oneline (See the :ref:`Synthetic Populatiohn<synthetic-population>` and :ref:`Population Parameter Guide<pop-params-guide>` for more information).

.. csv-table:: Additional columns
  :header:  "Properties", "Descriptions"
  :widths: 50, 150
  'interp_class_HMS_HMS', "Description of interp_class_HMS_HMS."
  'interp_class_CO_HMS_RLO', "Description of interp_class_CO_HMS_RLO."
  'interp_class_CO_HeMS', "Description of interp_class_CO_HeMS."
  'interp_class_CO_HeMS_RLO', "Description of interp_class_CO_HeMS_RLO."
  'mt_history_HMS_HMS', "Description of mt_history_HMS_HMS."
  'mt_history_CO_HMS_RLO', "Description of mt_history_CO_HMS_RLO."
  'mt_history_CO_HeMS', "Description of mt_history_CO_HeMS."
  'mt_history_CO_HeMS_RLO', "Description of mt_history_CO_HeMS_RLO."


State options
~~~~~~~~~~~~~

Binary states are defined according to the following table:

.. csv-table:: States
   :header: "State", "Description"
   :widths: 10, 30

   `detached`, "The stars in the binary are in a detached state."
   `contact`, "The stars in the binary are in contact."
   `RLO1`, "The binary is Roche Lobe overflowing, star 1 is overfilling the RL."
   `RLO2`, "The binary is Roche Lobe overflowing, star 2 is overfilling the RL."
   `CE`, "The binary is in a Common Envelope phase."
   `disrupted`, "The binary was disrupted."
   `MaxTimeChanged`, "Max time of the evolution was reached."

TODO: update properties

Event options
~~~~~~~~~~~~~

Binary events are defined according to the following table:

.. csv-table:: Events
  :header: "State", "Description"
  :widths: 10, 30

  `CC1`, "Core collapse of star 1."
  `CC2`, "Core collapse of star 2."
  `oRLO1`, "The binary is at onset of Roche Lobe overflow, star 1 is overfilling the RL."
  `oRLO2`, "The binary is at onset of Roche Lobe overflow, star 2 is overfilling the RL."
  `oCE`, "The binary is at the onset of Common Envelope."
  `None`, "No event occurred."
  `END`, "The binary evolution was stopped."

TODO: update properties

Mass Transfer case
~~~~~~~~~~~~~~~~~~

The mass transfer cases are stored in `mt_history_GRIDTYPE` and are defined according to the following table: TODO: add the table below

.. csv-table:: Mass transfer cases
  :header: "State", "Description"
  :widths: 10, 30

  `None`, "The binary is not Roche Lobe overflowing."


TODO: update properties


Basic example
~~~~~~~~~~~~~

The simplest method is to provide the two star objects and `kwargs` of the initial binary parameters.

.. code-block:: python

  kwargs1 = {'state' : 'MS',
             'mass' : 20.0,
             'metallicity' : 0.014}

  star_1 = SingleStar(**kwargs1)

  kwargs2 = {'state' : 'MS',
             'mass' : 10.0,
             'metallicity' : 0.014}

  star_2 = SingleStar(**kwargs2)

  kwargs3 = {'state' : 'detached',
             'event' : None,
             'time' : 0.,
             'orbital_period' : 3.,
             'eccentricity' : 0.}


  binary = BinaryStar(star_1, star_2, **kwargs3)
