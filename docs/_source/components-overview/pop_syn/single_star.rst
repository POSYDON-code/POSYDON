.. _single-star:


The SingleStar object
=====================

The Single Star object contains the current and past states of the star.
Only parameters in the ``STARPROPERTIES`` list are stored in the history.
The current parameter value of the star object is accessed as, e.g. ``star.mass`` while its past history with ``star.mass_history``.

To use SingleStar object import it using:

.. code-block:: python

  from posydon.binary_evol.singlestar import SingleStar


Creating a SingleStar object
----------------------------

STARPROPERTIES
~~~~~~~~~~~~~~

The star properties are defined as follows

.. list-table:: STARPROPERTIES
  :header-rows: 1
  :widths: 50 150

  * - Properties
    - Descriptions
  * - ``state``
    - The state of the star, see state options.
  * - ``metallicity``
    - Fractional metal content (Z) of the star.
  * - ``mass``
    - Stellar mass in M_sun.
  * - ``log_R``
    - log10 stellar radius in R_sun.
  * - ``log_L``
    - log10 surface stellar luminosity in L_sun.
  * - ``mdot``
    - Stellar mass-loss rate in Msun/yr.
  * - ``lg_wind_mdot``
    - log10 stellar wind mass-loss rate in Msun/yr.
  * - ``he_core_mass``
    - Helium core mass in M_sun.
  * - ``he_core_radius``
    - Helium core radius in R_sun.
  * - ``c_core_mass``
    - Carbon core mass in M_sun.
  * - ``c_core_radius``
    - Carbon core radius in R_sun.
  * - ``o_core_mass``
    - Oxygen core mass in M_sun.
  * - ``o_core_radius``
    - Oxygen core radius in R_sun.
  * - ``center_h1``
    - Hydrogen central mass fraction abundance.
  * - ``center_he4``
    - Helium central mass fraction abundance.
  * - ``center_c12``
    - Carbon central mass fraction abundance.
  * - ``center_n14``
    - Nitrogen central mass fraction abundance.
  * - ``center_o16``
    - Oxygen central mass fraction abundance.
  * - ``surface_h1``
    - Hydrogen surface mass fraction abundance.
  * - ``surface_he4``
    - Helium surface mass fraction abundance.
  * - ``surface_c12``
    - Carbon surface mass fraction abundance.
  * - ``surface_n14``
    - Nitrogen surface mass fraction abundance.
  * - ``surface_o16``
    - Oxygen surface mass fraction abundance.
  * - ``log_LH``
    - log10 total thermal power from PP and CNO, excluding neutrinos divided by L_sun.
  * - ``log_LHe``
    - log10 total thermal power from triple-alpha, excluding neutrinos divided by L_sun.
  * - ``log_LZ``
    - log10 total burning power excluding LH and LHe and photodisintegrations divided by L_sun.
  * - ``log_Lnuc``
    - log10 total nuclear reaction luminosity (LH + LHe LZ) in L_sun.
  * - ``c12_c12``
    - log10 total luminosity for c12_c12 reaction in L_sun.
  * - ``surf_avg_omega_div_omega_crit``
    - Average surface omega divided by critical omega.
  * - ``total_moment_of_inertia``
    - Total moment of inertia in g*cm^2.
  * - ``log_total_angular_momentum``
    - log10 total angular momentum of the star g*cm^2*s^-1
  * - ``spin``
    - Angular momentum of the star in g*cm^2*s^-1 or dimensionless BH spin.
  * - ``profile``
    - Stellar profile from MESA. [not currently supported for the initial-final interpolator]
  * - ``total_mass_h1``
    - Total Hydrogen mass in M_sun.
  * - ``total_mass_he4``
    - Total Helium mass in M_sun.

Additional scalar properties are added during the evolution depending on which steps the star has undergone. These properties are not stored in the history.

.. list-table:: Additional output
   :header-rows: 1
   :widths: 50 150

  * - Properties
    - Descriptions
  * - ``natal_kick_array``
    - | The natal kick array for the star if it has undergone a SN.
      | contains:

      * velocity
      * theta
      * phi
      * mean anomaly

  * - ``SN_type``
    - The supernova type of the star.
  * - ``f_fb``
    - The fraction of fallback mass.
  * - ``spin_orbit_tilt_first_SN``
    - The spin-orbit tilt after the first SN, if the star has undergone a SN.
  * - ``spin_orbit_tilt_second_SN``
    - The spin-orbit tilt after the second SN, if a second SN has occurred.
  * - ``m_disk_radiated``
    - The mass of the disk radiated in the collapse of the star.
  * - ``m_disk_accreted``
    - The mass of the disk accreted in the collapse of the star.


State options
~~~~~~~~~~~~~

Star states are defined by their burning and surface properties.
These states are combined to describe the stellar state.
We also have additional extra states for objects that are not stars.

.. list-table:: Surface state
  :header-rows: 1
  :widths: 10 30

  * - State
    - Description
  * - ``H-rich``
    - The star has a hydrogen-rich surface.
  * - ``stripped_He``
    - The star has a stripped helium surface.
  * - ``accreted_He``
    - The star has accreted a helium rich layer on its surface.

.. list-table:: Burning state
  :header-rows: 1
  :widths: 10 30

  * - State
    - Description
  * - ``non_burning``
    - The star is not burning.
  * - ``Core_H_burning``
    - The star is burning hydrogen in its core.
  * - ``Shell_H_burning``
    - The star is burning hydrogen in a shell.
  * - ``Core_He_burning``
    - The star is burning helium in its core.
  * - ``Core_He_depleted``
    - The star has a helium depleted core.
  * - ``Shell_He_burning``
    - The star is burning helium in a shell.
  * - ``Core_C_burning``
    - The star is burning carbon in its core.
  * - ``Core_C_depleted``
    - The star has a carbon depleted core.

.. list-table:: Additional States
  :header-rows: 1
  :widths: 10 30

  * - State
    - Description
  * - ``WD``
    - The star is a White Dwarf.
  * - ``NS``
    - The star is a Neutron Star.
  * - ``BH``
    - The star is a Black Hole.
  * - ``massless_remnant``
    - The star exploded or merged. Only its companion is left as a single star.

Basic example
~~~~~~~~~~~~~

The simplest method is to provide `kwargs` of the initial stellar parameters.

.. code-block:: python

  kwargs = {'state' : 'H-rich-Core_H_burning',
            'mass' : 10.0,
            'metallicity' : 0.014}
  SingleStar(**kwargs)

Now, the SingleStar object is ready to be used.
