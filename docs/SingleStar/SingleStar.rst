.. _SingleStar:

######################
The SingleStar object
######################

The Single Star object contains the current and past states of the star.
Only parameters in the STARPROPERTIES list are stored in the history.
The current parameter value of the star object is accessed as, e.g. `star.mass` while
its past history with `star.mass_history`.

To use SingleStar object import it using:

.. code-block:: python

  from posydon.binary_evol.singlestar import SingleStar


Creating a SingleStar object
============================

STARPROPERTIES
--------------

The star properties are defined as follows

.. csv-table:: STARPROPERTIES
   :header: "Properties", "Descriptions"
   :widths: 50, 150

   `state`, "The state of the star, see state options."
   `metallicity`, "Fractional metal content (Z) of the star."
   `mass`, "Stellar mass in M_sun."
   `log_R`, "log10 stellar radius in R_sun."
   `log_L`, "log10 surface stellar luminosity in L_sun."
   `mdot`, "Stellar mass-loss rate in Msun/yr."
   `lg_wind_mdot`, "log10 stellar wind mass-loss rate in Msun/yr."
   `he_core_mass`, "Helium core mass in M_sun."
   `he_core_radius`, "Helium core radius in R_sun."
   `c_core_mass`, "Carbon core mass in M_sun."
   `c_core_radius`, "Carbon core mass in M_sun."
   `o_core_mass`, "Oxygen core mass in M_sun."
   `o_core_radius`,  "Oxygen core radius in R_sun."
   `center_h1`, "Hydrogen central mass fraction abundance."
   `center_he4`, "Helium central mass fraction abundance."
   `center_c12`, "Carbon central mass fraction abundance."
   `center_n14`, "Nitrogen central mass fraction abundance."
   `center_o16`, "Oxygen central mass fraction abundance."
   `surface_h1`, "Hydrogen surface mass fraction abundance."
   `surface_he4`, "Helium surface mass fraction abundance."
   `surface_c12`, "Carbon surface mass fraction abundance."
   `surface_n14`, "Nitrogen surface mass fraction abundance."
   `surface_o16`, "Oxygen surface mass fraction abundance."
   `log_LH`, "log10 total thermal power from PP and CNO, excluding neutrinos devided by L_sun."
   `log_LHe`, "log10 total thermal power from triple-alpha, excluding neutrinos devided by L_sun."
   `log_LZ`, "log10 total burning power excluding LH and LHe and photodisintegrations devided by L_sun."
   `log_Lnuc`, "log10 total nuclear reaction luminosity (LH + LHe LZ) in L_sun."
   `c12_c12`, "log10 total luminosity for c12_c12 reaction in L_sun."
   `surf_avg_omega_div_omega_crit`, "Average surface omega divided by critical omega."
   `total_moment_of_inertia`, "Total moment of inertia in g*cm^2."
   `log_total_angular_momentum`, "log10 total angular momentum of the star g*cm^2*s^-1"
   `spin`, "Angular momentum of the star in g*cm^2*s^-1 or dimensionless BH spin."
   `profile`, "Stellar profile from MESA."

State options
-------------

Star states are defined on their relative position on the HR diagram as follows:

.. csv-table:: State
   :header: "State", "Description"
   :widths: 10, 30

   `PreMS`, "The star is on the pre-main-sequence phase."
   `MS`, "The star is in the Main Sequence phase."
   `PostMS`, "The star is in the post Main Sequence phase."
   `HeMS`, "The star is at the beginning of the He burning phase."
   `PostHeMS`, "The star is evolving through the beginning of the He burning phase."
   `WD`, "The star is a White Dwarf."
   `NS`, "The star is a Neutron Star."
   `BH`, "The star is a Black Hole."

Basic example
-------------
The simplest method is to provide `kwargs` of the initial stellar parameters.

.. code-block:: python

  kwargs = {'state' : 'MS',
            'mass' : 10.0,
            'metallicity' : 0.014}
  SingleStar(**kwargs)

Now, the SingleStar object is ready to be used.
