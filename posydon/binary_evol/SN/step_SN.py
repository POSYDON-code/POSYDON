"""Supernova step.

This step models the end of life of stars by being applied to a binary
object and verifying its state. It performs the collapse prescription
used to initialize the step in the respective star. Depending on the
C and He cores the final state of the star is determined, from the
formation of white dwarfs to electron-capture supernova, Fe core-collapse
supernova, pair pulsation supernova and pair instability supernova.

"""

__authors__ = [
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Jaime Roman Garza <Jaime.Roman@etu.unige.ch>",
    "Emmanouil Zapartas <ezapartas@gmail.com>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Devina Misra <devina.misra@unige.ch>",
    "Zepei Xing <Zepei.Xing@unige.ch>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
    "Tassos Fragos <Anastasios.Fragkos@unige.ch>",
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
]

__credits__ = [
    "Michael Zevin <michael.zevin@ligo.org>",
    "Chase Kimball <charles.kimball@ligo.org",
    "Sam Imperato <samuelimperato2022@u.northwestern.edu>",
]


import os
import warnings
import numpy as np
import scipy as sp
import copy
import pandas as pd

from posydon.utils.data_download import PATH_TO_POSYDON_DATA, data_download
import posydon.utils.constants as const
from posydon.utils.common_functions import is_number
from posydon.utils.common_functions import CO_radius
from posydon.utils.common_functions import (
    orbital_period_from_separation,
    inspiral_timescale_from_separation,
    separation_evol_wind_loss,
    calculate_Patton20_values_at_He_depl,
    THRESHOLD_CENTRAL_ABUNDANCE,
    rotate
)

from posydon.binary_evol.binarystar import BINARYPROPERTIES
from posydon.binary_evol.singlestar import STARPROPERTIES, convert_star_to_massless_remnant
from posydon.binary_evol.SN.profile_collapse import do_core_collapse_BH
from posydon.binary_evol.flow_chart import (STAR_STATES_CO, STAR_STATES_CC,
                                            STAR_STATES_C_DEPLETION)

from posydon.grids.MODELS import MODELS

from pandas import read_csv
from sklearn import neighbors
from scipy.interpolate import interp1d

import json


path_to_Sukhbold_datasets = os.path.join(PATH_TO_POSYDON_DATA,
                                         "Sukhbold+16/")

path_to_Patton_datasets = os.path.join(PATH_TO_POSYDON_DATA,
                                       "Patton+Sukhbold20/")

path_to_Couch_datasets = os.path.join(PATH_TO_POSYDON_DATA,
                                      "Couch+2020/")

MODEL = {
    # core collapse physics
    "mechanism": 'Patton&Sukhbold20-engine',
    "engine": 'N20',
    "PISN": "Marchant+19",
    "ECSN": "Podsiadlowksi+04",
    "conserve_hydrogen_envelope" : False,
    "max_neutrino_mass_loss": 0.5,
    "max_NS_mass": 2.5,
    "use_interp_values": True,
    "use_profiles": True,
    "use_core_masses": True,
    "allow_spin_None" : False,
    "approx_at_he_depletion": False,
    # kick physics
    "kick": True,
    "kick_normalisation": 'one_over_mass',
    "sigma_kick_CCSN_NS": 265.0,
    "sigma_kick_CCSN_BH": 265.0,
    "sigma_kick_ECSN": 20.0,
    # other
    "verbose": False,
}


class StepSN(object):
    """The supernova step in POSYDON.

    Keyword Arguments
    ----------
    mechanism : str
        Mechanism to perform the core-collapse on the star object and
        predict the supernova remnant outcome. Available options are:

        * 'Fryer+12-rapid' : The rapid supernova-engine described in [1]_

        * 'Fryer+12-delayed' : The delayed supernova-engine described in [1]_

        * 'direct' : The pre-supernova mass of the starr is collapsed into the
        remnant baryonic mass.

        * 'direct_he_core' : The pre-supernova He core mass of the starr is
        collapsed into the remnant baryonic mass.

        * 'Sukhbold+16-engine' : Uses the results from [2]_
        to describe the collapse of the star.

        * 'Patton&Sukhbold20-engine': Uses the results from [5]_
        to describe the collapse of the star.

        * 'Couch+20-engine': Uses the results from [6]_
        to describe the collapse of the star.

    engine : str
        Engine used for supernova remnanrt outcome propierties for the
        Sukhbold+16-engineand and Patton&Sukhbold20-engine mechanisms.
        Available options:

        - 'N20'

    PISN : str
        Prescrition to take on the pair-instability supernova.
        Avialable options:

        - 'Marchant+19' : Descripes the pair-instability supernova as [3]_.

    mass_central_BH : double
        Central mass collapsed automatically on black-holes formed by direct
        collapse.

    max_neutrino_mass_loss : double
        Neutrino mass loss during the collapse of the proto neutron-star.

    kick : bool
        If True, the kick velocities are computed corresponding to the
        supernova event, else no kicks are taking into account
        for any supernova outcome.

    kick_normalisation : str
        Renormalise the kick by:
        'one_minus_fallback' : (1-f_fb)
        'one_over_mass' : 1.4/m_BH
        'zero' : 0.
        'one' : 1.
        'NS_one_minus_fallback_BH_one': 1 for BH, (1-f_fb) for NS

    ECSN : str
        Prescription to determine the production of an electron-capture
        supernova.
        Avialable options:

        - 'Tauris+15': Determines the electron capture supernova in terms
        of the CO core mass at pre-supernova, taking the limits from [4]_.

    sigma_kick_CCSN_NS : double
        Standard deviation for a Maxwellian distribution to compute the
        kick velocities from NSs formed by Fe core-collapse
        supernova.

    sigma_kick_CCSN_BH : double
        Standard deviation for a Maxwellian distribution to compute the
        kick velocities from BHs formed by Fe core-collapse
        supernova.

    sigma_kick_ECSN : double
        Standard deviation for a Maxwellian distribution to compute the
        kick velocities from compact-object formed by electron-capture
        supernova.

    max_NS_mass : double
        Maximum neutron-star mass.

    use_interp_values : bool
       The outcome of core collpase was interpolated from a post processed
       MESA grid and stored in the star object in the mesa_step or
       detached_step (default). This option supports only default
       assumptions for all core collase mechanism.

    use_profiles : bool
       Perfrome the core collpase given a MESA profile. To use this option
       a MESA profile must be stored in the star object which is provided
       by nearest neighbor interpolation in the mesa_step or (TODO)
       interpolated in the detached_step.

    use_core_masses : bool
       This option uses the core masses at carbon depletion to determine
       the core collapse outcoume (classical population sythesis
       threatment).

    allow_spin_None : bool
       This option does not determine the spin during core collapse while
       setting other values like in use_core_masses. (used to avoid jumps
       in the spin for interpolator training because of missing profiles)

    approx_at_he_depletion : bool
       This option is relevant only for the mechanism Patton&Sukhbold20-engine.
       In case the core masses at he-depletion are not present in the
       star object, compute them from the history, else (approximation=True)
       approximate it from the core masses at C depletion.

    verbose : bool
        If True, the messages will be prited in the console.

    References
    ----------
    .. [1] Fryer, C. L., Belczynski, K., Wiktorowicz, G., Dominik, M.,
        Kalogera, V., & Holz, D. E. (2012). Compact remnant mass function:
        dependence on the explosion mechanism and metallicity. The
        Astrophysical Journal, 749(1), 91.

    .. [2] Sukhbold, T., Ertl, T., Woosley, S. E., Brown, J. M., & Janka, H. T.
        (2016). Core-collapse supernovae from 9 to 120 solar masses based on
        neutrino-powered explosions. The Astrophysical Journal, 821(1), 38.

    .. [3] Marchant, P., Renzo, M., Farmer, R., Pappas, K. M., Taam, R. E., De
        Mink, S. E., & Kalogera, V. (2019). Pulsational pair-instability
        supernovae in very close binaries. The Astrophysical Journal, 882(1), 36.

    .. [4] Tauris, T. M., Langer, N., & Podsiadlowski, P. (2015).
        Ultra-stripped supernovae: progenitors and fate. Monthly Notices of the
        Royal Astronomical Society, 451(2), 2123-2144.

    .. [5] Patton, R. A. & Sukhbold, T. 2020, MNRAS, 499, 2803. Towards a
        realistic explosion landscape for binary population synthesis

    .. [6] Couch, S. M., Warren, M. L., & O’Connor, E. P. 2020, ApJ, 890, 127.
        Simulating Turbulence-aided Neutrino-driven Core-collapse Supernova
        Explosions in One Dimension

    """

    def __init__(self, **kwargs):
        """Initialize a StepSN instance."""
        # read kwargs to initialize the class
        if kwargs:
            for key in kwargs:
                if key not in MODEL:
                    raise ValueError(key + " is not a valid parameter name!")
            for varname in MODEL:
                default_value = MODEL[varname]
                setattr(self, varname, kwargs.get(varname, default_value))
        else:
            for varname in MODEL:
                default_value = MODEL[varname]
                setattr(self, varname, default_value)

        if self.max_neutrino_mass_loss is None:
            self.max_neutrino_mass_loss = 0

        # Initializing core collapse

        # Available mechanisms for core-collapse supernova
        self.Fryer12_rapid = "Fryer+12-rapid"
        self.Fryer12_delayed = "Fryer+12-delayed"
        self.direct_collapse = "direct"
        self.direct_collapse_hecore = "direct_he_core"
        self.Sukhbold16_engines = "Sukhbold+16-engine"
        self.Patton20_engines = "Patton&Sukhbold20-engine"
        self.Couch20_engines = "Couch+20-engine"

        self.mechanisms = [
            self.Fryer12_rapid,
            self.Fryer12_delayed,
            self.direct_collapse,
            self.direct_collapse_hecore,
            self.Sukhbold16_engines,
            self.Patton20_engines,
            self.Couch20_engines
        ]

        if self.mechanism in self.mechanisms:

            if self.mechanism in [
                self.Fryer12_rapid,
                self.Fryer12_delayed,
                self.direct_collapse,
                self.direct_collapse_hecore,
            ]:
                self.Sukhbold_corecollapse_engine = None

            elif self.mechanism == self.Sukhbold16_engines:
                # set the path to the datasets for each supernova engine
                self.path_to_Sukhbold_datasets = path_to_Sukhbold_datasets
                self.Sukhbold_corecollapse_engine = Sukhbold16_corecollapse(
                    self.engine, self.path_to_Sukhbold_datasets, self.verbose)

            elif self.mechanism == self.Couch20_engines:
                # set the path to the datasets for each supernova engine
                self.path_to_Couch_datasets = path_to_Couch_datasets

                # returns JSON object as
                # a dictionary
                self.Couch_corecollapse_engine = Couch20_corecollapse(
                    turbulence_strength=self.engine,
                    path_engine_dataset=self.path_to_Couch_datasets,
                    verbose=self.verbose)

            elif self.mechanism == self.Patton20_engines:
                self.path_to_Patton_datasets = path_to_Patton_datasets

                def format_data_Patton20(file_name):
                    """Format the Patton&Sukhbold,20 dataset for interpolation.

                    Parameters
                    ----------
                    file_name : str
                        Name of the dataset file.

                    Returns
                    -------
                    CO_core_params : arr
                        Array containing the carbon-oxygen core parameters
                        in a grid of abundance and mass as columns.
                    target_parameter : arr
                        Array with the corresponding value of the target
                        parameter giving the selected dataset.
                    """

                    # Check if interpolation files exist
                    filename = os.path.join(self.path_to_Patton_datasets,
                                            file_name)
                    if not os.path.exists(filename):
                        data_download()

                    # Reading the dataset
                    data = np.loadtxt(filename, skiprows=6, dtype='str')

                    # Extracting the matrix with the values of the target
                    # parameter
                    target_matrix = data[1:].T[1:].T

                    # Formating the target metrix values as a 1D array and
                    # converting the values to float
                    target = target_matrix.astype(float).ravel()

                    # Extracting the values for the CO core parameters
                    M_CO, X_CO = np.meshgrid(data[0][1:], data.T[0][1:])

                    # Stacking the CO core parameters to the corresponding grid
                    # array that defines the injective relation between each
                    # element of CO_core_params to the elements in target
                    CO_core_params = (np.vstack(
                        (X_CO.ravel(), M_CO.ravel())).T).astype(float)

                    return CO_core_params, target

                if self.verbose:
                    print('Loading the train dataset for engine mu4 and M4...')
                CO_core_params_mu4, mu4_target = format_data_Patton20(
                    'Kepler_mu4_table.dat')
                CO_core_params_M4, M4_target = format_data_Patton20(
                    'Kepler_M4_table.dat')

                n_neighbors = 5

                if self.verbose:
                    print('Training the classifier ...')
                self.M4_interpolator = neighbors.KNeighborsRegressor(
                    n_neighbors, weights='distance')
                self.M4_interpolator.fit(CO_core_params_M4, M4_target)

                self.mu4_interpolator = neighbors.KNeighborsRegressor(
                    n_neighbors, weights='distance')
                self.mu4_interpolator.fit(CO_core_params_mu4, mu4_target)
                if self.verbose:
                    print('Done')
        else:
            raise ValueError("Invalid core-collapse mechanism given.")

    def __repr__(self):
        """Get the string representation of the class and any parameters."""
        return "SN step (kick distribution: {})".format(self.kick_distribution)

    def _reset_other_star_properties(self, star):
        """Reset the properties of the star that is not being collapsed."""
        star.lg_mdot = None
        star.lg_system_mdot = None

    def __call__(self, binary):
        """Perform the supernova step on a binary object.

        Parameters
        ----------
        binary : instance of BinaryStar
            The binary to evolve.
        """
        # consistency check
        # self.check()

        # read binary properties of interest
        # do the caclulations
        # update star/binary properties (e.g. period, eccentricity, masses)

        # Check if the binary event is calling correctly the SN_step,
        # this should occour only on the first or second core-collapse
        # CC1 and CC2 respectively.
        print("Manos, calling class SN")
        if binary.event == "CC1":
            # collapse star
            print("Manos, starting CC1")
            self.collapse_star(star=binary.star_1)
            self._reset_other_star_properties(star=binary.star_2)
            print("Manos CC1 outcome", binary, binary.star_1.SN_type, binary.star_1.M4, binary.star_1.mu4, binary.star_2.SN_type, binary.star_2.M4, binary.star_2.mu4)

        elif binary.event == "CC2":
            # collapse star
            print("Manos, starting CC2")
            self.collapse_star(star=binary.star_2)
            self._reset_other_star_properties(star=binary.star_1)
            print("Manos CC2 outcome", binary, binary.star_1.SN_type, binary.star_1.M4, binary.star_1.mu4, binary.star_2.SN_type, binary.star_2.M4, binary.star_2.mu4)
        else:
            raise ValueError("Something went wrong: "
                             "invalid call of supernova step!")

        # do orbital_kick on the binary object
        if self.kick:
            self.orbital_kick(binary=binary)

        # Checks if the binary is not disrupted to compute the
        # inspiral time due to gravitational wave emission
        state1, state2 = binary.star_1.state, binary.star_2.state
        if binary.state == "disrupted" or state1 == "massless_remnant" or state2 == "massless_remnant":
            binary.inspiral_time = np.nan
        elif state1 in STAR_STATES_CO and state2 in STAR_STATES_CO:
            binary.inspiral_time = inspiral_timescale_from_separation(
                binary.star_1.mass,
                binary.star_2.mass,
                binary.separation,
                binary.eccentricity,
            )
        # Cover the case where CC of the companion is immediately followed
        elif state1 in STAR_STATES_CO and state2 in STAR_STATES_C_DEPLETION:
            binary.event = "CC2"
        elif state1 in STAR_STATES_C_DEPLETION and state2 in STAR_STATES_CO:
            binary.event = "CC1"

    def check(self):
        """Check the internal integrity and the values of the parameters."""
        if self.kick_distribution is None:
            raise ValueError("Undefined supernova kick velocity distribution")

    def collapse_star(self, star):
        """Collapse the star object into a compact object.

        This routine supports three options:
        1. use_interp_values : True
           The outcome of core collpase was interpolated from a post processed
           MESA grid and stored in the star object in the mesa_step or
           detached_step (default). This option supports only default
           assumptions for all core collase mechanism.
        2. use_profiles : False
           Perfrome the core collpase given a MESA profile. To use this option
           a MESA profile must be stored in the star object which is provided
           by nearest neighbor interpolation in the mesa_step or (TODO)
           interpolated in the detached_step.
        3. use_core_masses : False
           This option uses the core masses at carbon depletion to determine
           the core collapse outcoume (classical population sythesis
           threatment).
        4. allow_spin_None : False
           This option does not determine the spin during core collapse while
           setting other values like in use_core_masses. (used to avoid jumps
           in the spin for interpolator training because of missing profiles)

        Parameters
        ----------
        star : object
            Star object containing the star properties.


        Returns
        -------
        m_rem : double
            Remnant mass of the compact object in M_sun. This quantity accounts
            for the mass loss thorugh neutrino.

        state : string
            New state of the star object.

        """
        state = star.state
        # after this function is called certain quantities shouldn't be None
        # type objects anymore
        for key in ['m_disk_accreted', 'm_disk_radiated']:
            if getattr(star, key) is None:
                setattr(star, key, np.nan)

        # Verifies if the star is in state state where it can
        # explode
        if state in STAR_STATES_CC:

            # if no profile is avaiable but interpolation quantities are,
            # use those, else continue with or without profile.
            if self.use_interp_values:
                # find MODEL_NAME corresponding to class variable
                print("Manos, entering use_interp_values")
                MODEL_NAME_SEL = None
                for MODEL_NAME, MODEL in MODELS.items():
                    tmp = MODEL_NAME
                    for key, val in MODEL.items():
                        if "use_" in key:
                            continue
                        if getattr(self, key) != val:
                            if self.verbose:
                                print(tmp, 'mismatch:', key, getattr(self, key), val)
                            tmp = None
                            break
                    if tmp is not None:
                        if self.verbose:
                            print('matched to model:', tmp)
                        MODEL_NAME_SEL = tmp

                # check if selected MODEL is supported
                if MODEL_NAME_SEL is None:
                    raise ValueError('Your model assumptions are not'
                                     'supported!')
                elif getattr(star, MODEL_NAME_SEL) is None:
                    # NOTE: this option is needed to do the collapse
                    # for stars evolved with the step_detached or
                    # step_disrupted.
                    # allow to continue with the collapse with profile
                    # or core masses
                    warnings.warn(f'{MODEL_NAME_SEL}: The collapsed star '
                                     'was not interpolated! If use_profiles '
                                     'or use_core_masses is set to True, '
                                     'continue with the collapse.')
                else:
                    # store some properties of the star object
                    # to be used for collapse verification
                    pre_SN_star = copy.deepcopy(star)

                    MODEL_properties = getattr(star, MODEL_NAME_SEL)
                    for key, value in MODEL_properties.items():
                        setattr(star, key, value)

                    if star.state == 'WD':
                        for key in STARPROPERTIES:
                            if key in ["he_core_mass"]:
                                setattr(star, key, star.mass)
                            elif key in ["co_core_mass"]:
                                if star.center_he4 < THRESHOLD_CENTRAL_ABUNDANCE:
                                    setattr(star, key, star.mass)
                                else:
                                    setattr(star, key, 0.)
                            elif key not in ["state", "mass", "spin",
                                             "m_disk_accreted",
                                             "m_disk_radiated", "center_h1",
                                             "center_he4", "center_c12",
                                             "center_n14", "center_o16"]:
                                setattr(star, key, None)

                    else:
                        for key in STARPROPERTIES:
                            if key not in ["state", "mass", "spin",
                                           "m_disk_accreted",
                                           "m_disk_radiated"]:
                                setattr(star, key, None)

                    # check if SN_type matches the predicted CO
                    # and force the SN_type to match the predicted CO.
                    # ie WD is no SN
                    # 1. Check if SN_type and star state match
                    # Non-matching SN_type and star state
                    if not check_SN_CO_match(star):
                        # raise a warning
                        warnings.warn(f'{MODEL_NAME_SEL}: The SN_type '
                                      'does not match the predicted CO! '
                                      'Recalculating the SN_type and CO')
                        # recalculate the SN_type and CO
                        # change some star properties back
                        m_PISN = self.PISN_prescription(pre_SN_star)
                        # no remnant if a PISN happens
                        if pd.isna(m_PISN):
                            star.SN_type = 'PISN'
                            star.state = 'massless_remnant'
                        else:
                            _, _, star.state = self.compute_m_rembar(pre_SN_star, m_PISN)
                            star.SN_type = self.check_SN_type(pre_SN_star.c_core_mass,
                                                            pre_SN_star.he_core_mass,
                                                            pre_SN_star.mass)[-1]

                        # check if the new SN_type matches new SN_type
                        if not check_SN_CO_match(star):
                            # still doesn't match
                            # raise a warning
                            warnings.warn('The SN_type still does not match. '
                                          'Forced the SN type to match the '
                                          'predicted CO.')
                            if star.state == 'WD':
                                star.SN_type = 'WD'
                            elif star.state == 'NS' or star.state == 'BH':
                                star.SN_type = 'CCSN'
                            elif star.state == 'massless_remnant':
                                star.SN_type = 'PISN'
                            else:
                                raise ValueError('Star state not recognized.')

                    del pre_SN_star

                    # No remnant if a PISN happens
                    if star.SN_type == 'PISN':
                        convert_star_to_massless_remnant(star=star)
                        # the mass is set to None
                        # but an orbital kick is still applied.
                        # Since the mass is set to None, this will lead to a disruption
                        # TODO: make it skip the kick caluclation

                    if getattr(star, 'SN_type') != 'PISN':
                        star.log_R = np.log10(CO_radius(star.mass, star.state))
                    return

            # Verifies the selection of core-collapse mechnism to perform
            # the collapse
            if self.mechanism in [
                self.Fryer12_rapid,
                self.Fryer12_delayed,
                self.direct_collapse,
                self.direct_collapse_hecore,
            ]:
                # m_core = star.co_core_mass

                # this flag checks if a profile is available
                profile = star.profile

                # computes the baryonic remnant mass from the
                # PISN and PPISN prescription if the star will
                # experience such event
                m_PISN = self.PISN_prescription(star)

                # the baryonic remnant mass is computed in terms
                # of the core mass.
                m_rembar, star.f_fb, _ = self.compute_m_rembar(star, m_PISN)

                # check if a white dwarf has been born
                if star.SN_type == "WD":
                    star.mass = m_rembar
                    star.state = "WD"
                    star.spin = 0.
                    star.log_R = np.log10(CO_radius(star.mass, star.state))
                    for key in STARPROPERTIES:
                        if key in ["he_core_mass"]:
                            setattr(star, key, star.mass)
                        elif key in ["co_core_mass"]:
                            if star.center_he4 < THRESHOLD_CENTRAL_ABUNDANCE:
                                setattr(star, key, star.mass)
                            else:
                                setattr(star, key, 0.)
                        elif key not in ["state", "mass", "spin",
                                         "m_disk_accreted", "m_disk_radiated",
                                         "center_h1", "center_he4",
                                         "center_c12", "center_n14",
                                         "center_o16"]:
                            setattr(star, key, None)
                    return

                # check if the star was disrupted by the PISN
                if np.isnan(m_rembar):
                    convert_star_to_massless_remnant(star=star)
                    return

                # Computing the gravitational mass of the remnant
                # as in Lattimer & Yahil, 1989
                m_grav = (20.0 / 3.0) * (np.sqrt(1.0 + 0.3 * m_rembar) - 1.0)
                if (m_rembar - m_grav) > self.max_neutrino_mass_loss:
                    m_grav = m_rembar - self.max_neutrino_mass_loss

                # If the profile of the star is available then
                # it will be collapsed to get the information
                # on the compact object spin
                if self.use_profiles and profile is not None:
                    delta_M = m_rembar - m_grav
                    if delta_M > self.max_neutrino_mass_loss:
                        delta_M = self.max_neutrino_mass_loss
                    if m_grav >= self.max_NS_mass:
                        mass_direct_collapse = self.max_NS_mass + delta_M
                        final_BH = do_core_collapse_BH(
                            star=star, mass_collapsing=m_rembar,
                            mass_central_BH=mass_direct_collapse,
                            neutrino_mass_loss=delta_M,
                            max_neutrino_mass_loss=self.max_neutrino_mass_loss,
                            verbose=self.verbose
                        )
                        star.mass = final_BH[0]
                        star.spin = final_BH[1]
                        star.m_disk_accreted = final_BH[2]
                        star.m_disk_radiated = final_BH[3]
                        star.state = "BH"
                    else:
                        star.mass = m_grav
                        star.spin = 0.
                        star.m_disk_accreted = 0.
                        star.m_disk_radiated = 0.
                        star.state = 'NS'

                elif self.use_core_masses:
                    # If the profile is not available the star spin
                    # is used to get the compact object spin
                    star.mass = m_grav
                    if m_grav >= self.max_NS_mass:
                        # see Eq. 14, Fryer, C. L., Belczynski, K., Wiktorowicz,
                        # G., Dominik, M., Kalogera, V., & Holz, D. E. (2012), ApJ, 749(1), 91.

                        # assume the spin value is the AM of the star
                        # convert to CGS units
                        G = const.standard_cgrav
                        c = const.clight
                        Mo = const.Msun
                        star.spin = (10**star.log_total_angular_momentum * c
                                     / (G * (m_grav * Mo) ** 2))
                        if star.spin > 1.0:
                            if self.verbose:
                                print("The spin exceeds 1, capping it to 1...")
                            star.spin = 1.0
                        star.m_disk_accreted = 0.0
                        star.m_disk_radiated = 0.0
                        star.state = "BH"
                    else:
                        star.spin = 0.0
                        star.m_disk_accreted = 0.0
                        star.m_disk_radiated = 0.0
                        star.state = "NS"

                elif self.allow_spin_None:
                    # If the profile is not available and spin can stay
                    # undetermined
                    star.mass = m_grav
                    star.spin = None
                    star.m_disk_accreted = 0.0
                    star.m_disk_radiated = 0.0
                    if m_grav >= self.max_NS_mass:
                        star.state = "BH"
                    else:
                        star.state = "NS"

                else:
                    for key in STARPROPERTIES:
                        setattr(star, key, None)
                    star.state = "ERR"
                    raise ValueError("FAILED core collapse!")

            elif self.mechanism in [self.Sukhbold16_engines,
                                    self.Patton20_engines,
                                    self.Couch20_engines]:
                print("Manos, self.mechanism in Patton20_engines")
                # The final remnant mass and and state
                # is computed by the selected mechanism

                # PISN and PPISN prescription
                m_PISN = self.PISN_prescription(star)

                m_rembar, star.f_fb, state = self.compute_m_rembar(star,
                                                                   m_PISN)
                print("Manos after computer_m_rembar", star.SN_type, star.M4, star.mu4)
                star.state = state

                # check if a white dwarf has been born
                if star.SN_type == "WD":
                    star.mass = m_rembar
                    star.state = "WD"
                    star.spin = 0.
                    star.log_R = np.log10(CO_radius(star.mass, star.state))
                    for key in STARPROPERTIES:
                        if key in ["he_core_mass"]:
                            setattr(star, key, star.mass)
                        elif key in ["co_core_mass"]:
                            if star.center_he4 < THRESHOLD_CENTRAL_ABUNDANCE:
                                setattr(star, key, star.mass)
                            else:
                                setattr(star, key, 0.)
                        elif key not in ["state", "mass", "spin",
                                         "m_disk_accreted", "m_disk_radiated",
                                         "center_h1", "center_he4",
                                         "center_c12", "center_n14",
                                         "center_o16"]:
                            setattr(star, key, None)
                    return

                # check if the star was disrupted by the PISN
                if np.isnan(m_rembar):
                    convert_star_to_massless_remnant(star=star)
                    return

                # Computing the gravitational mass of the remnant
                # as in Lattimer & Yahil, 1989
                m_grav = (20.0 / 3.0) * (np.sqrt(1.0 + 0.3 * m_rembar) - 1.0)
                if (m_rembar - m_grav) > self.max_neutrino_mass_loss:
                    m_grav = m_rembar - self.max_neutrino_mass_loss

                # this flag checks if a profile is available
                profile = star.profile

                if self.use_profiles and profile is not None:
                    print("Manos entering use_profiles", star.SN_type, star.M4, star.mu4)
                    delta_M = m_rembar - m_grav
                    if delta_M > self.max_neutrino_mass_loss:
                        delta_M = self.max_neutrino_mass_loss
                    if m_grav >= self.max_NS_mass and star.state == "BH":
                        mass_direct_collapse = self.max_NS_mass + delta_M
                        final_BH = do_core_collapse_BH(
                            star=star, mass_collapsing=m_rembar,
                            mass_central_BH=mass_direct_collapse,
                            neutrino_mass_loss=delta_M,
                            max_neutrino_mass_loss=self.max_neutrino_mass_loss,
                            verbose=self.verbose
                        )
                        star.mass = final_BH[0]
                        if m_grav != star.mass and self.verbose:
                            print("The star formed a disk during the collapse "
                                  "and lost", round(final_BH[0] - m_rembar, 2),
                                  "M_sun.")
                        star.spin = final_BH[1]
                        star.m_disk_accreted = final_BH[2]
                        star.m_disk_radiated = final_BH[3]
                    elif star.state == "NS":
                        star.mass = m_grav
                        star.m_disk_accreted = 0.0
                        star.m_disk_radiated = 0.0
                        star.spin = 0.0
                    else:
                        for key in STARPROPERTIES:
                            setattr(star, key, None)
                        star.state = "ERR"
                        raise ValueError("Invalid core state", state)

                elif self.use_core_masses:
                    print("Manos entering use_core_masses", star.SN_type, star.M4, star.mu4)
                    star.mass = m_grav
                    if m_grav >= self.max_NS_mass:
                        # see Eq. 14, Fryer, C. L., Belczynski, K., Wiktorowicz,
                        # G., Dominik, M., Kalogera, V., & Holz, D. E. (2012), ApJ, 749(1), 91.

                        # assume the spin value is the AM of the star
                        # convert to CGS units
                        G = const.standard_cgrav
                        c = const.clight
                        Mo = const.Msun
                        star.spin = (10**star.log_total_angular_momentum * c
                                     / (G * (m_grav * Mo) ** 2))
                        if star.spin > 1.0:
                            if self.verbose:
                                print("The spin exceed 1, capping it to 1...")
                            star.spin = 1.0
                        star.m_disk_accreted = 0.0
                        star.m_disk_radiated = 0.0
                        star.state = "BH"
                    else:
                        star.spin = 0.0
                        star.m_disk_accreted = 0.0
                        star.m_disk_radiated = 0.0
                        star.state = "NS"

                elif self.allow_spin_None:
                    # If the profile is not available and spin can stay
                    # undetermined
                    star.mass = m_grav
                    star.spin = None
                    star.m_disk_accreted = 0.0
                    star.m_disk_radiated = 0.0
                    if m_grav >= self.max_NS_mass:
                        star.state = "BH"
                    else:
                        star.state = "NS"

                else:
                    for key in STARPROPERTIES:
                        setattr(star, key, None)
                    star.state = "ERR"
                    raise ValueError("FAILED core collapse!")

        else:
            raise ValueError(f"The star cannot collapse: star state {state}.")

        star.metallicity = star.metallicity_history[-1]

        star.log_R = np.log10(CO_radius(star.mass, star.state))

        for key in STARPROPERTIES:
            if key not in ["state", "mass", "spin", "log_R", "metallicity",
                           "m_disk_accreted", "m_disk_radiated",
                           "co_core_mass"]:
                setattr(star, key, None)
        print("Manos end of collapse_star", star.SN_type, star.M4, star.mu4)

    def PISN_prescription(self, star):
        """Compute baryonic remnant mass for the PPISN and PISN prescription.

        Parameters
        ----------
        star : object
            Star object containing the star properties.

        Returns
        -------
        m_PISN : double
            Maximum stellar mass in M_sun after the PPISN/PISN prescription.

        """
        if self.PISN is None:
            return None

        else:
            # perform the PISN prescription in terms of the
            # He core mass at pre-supernova
            m_He_core = star.he_core_mass
            m_star = star.mass
            if self.PISN == "Marchant+19":
                if m_He_core >= 31.99 and m_He_core <= 61.10:
                    # this is the 8th-order polynomial fit of table 1
                    # value, see COSMIC paper (Breivik et al. 2020)
                    polyfit = (
                        -6.29429263e5
                        + 1.15957797e5 * m_He_core
                        - 9.28332577e3 * m_He_core ** 2.0
                        + 4.21856189e2 * m_He_core ** 3.0
                        - 1.19019565e1 * m_He_core ** 4.0
                        + 2.13499267e-1 * m_He_core ** 5.0
                        - 2.37814255e-3 * m_He_core ** 6.0
                        + 1.50408118e-5 * m_He_core ** 7.0
                        - 4.13587235e-8 * m_He_core ** 8.0
                    )
                    m_PISN = polyfit

                elif m_He_core > 61.10 and m_He_core < 113.29:
                    m_PISN = np.nan

                else:
                    # above the PISN gap we assume direct collapse of the
                    if self.conserve_hydrogen_envelope:
                        m_PISN = m_star
                    else:
                        m_PISN = m_He_core

            elif is_number(self.PISN) and m_He_core > self.PISN:
                m_PISN = self.PISN

            elif is_number(self.PISN) and 0.0 < m_He_core <= self.PISN:
                m_PISN = None

            else:
                raise ValueError(
                    "This choice {} of PISN is not availabe!".format(self.PISN)
                )

        if self.verbose:
            if m_PISN is None:
                print("")
                print("The star did NOT lose any mass because of "
                      "PPIN or PISN.")
            elif not np.isnan(m_PISN):
                print("")
                print(
                    "The star with initial mass {:2.2f}".format(m_He_core),
                    "M_sun went through the PISN routine and lost",
                    "{:2.2f} M_sun.".format(m_He_core - m_PISN),
                    "The new m_rembar mass that will collapse to form a ",
                    "CO object is {:2.2f} M_sun.".format(m_PISN))
            else:
                print("The star was disrupted by the PISN prescription!")

        return m_PISN

    def check_SN_type(self, m_core, m_He_core, m_star):
        """Get the remnant mass, fallback frac., state & SN type of the SN."""
        if self.ECSN == "Tauris+15":
            # Label the supernova type as in Tauris et al. (2015),
            # considering their definition of metal core quivalent
            # to the mass of the CO core the the star object at pre-SN
            min_M_CO_ECSN = 1.37  # Msun from Takahashi et al. (2013)
            max_M_CO_ECSN = 1.43  # Msun from Tauris et al. (2015)

            if m_core < min_M_CO_ECSN:
                # The birth of a white dwarf is assumed
                SN_type = "WD"

                if m_core > 0.:
                    # co_core_mass, note there will be no kick
                    m_rembar = m_core
                elif m_He_core > 0.:
                    m_rembar = m_He_core
                else:
                    # this is catching H-rich_non_burning stars
                    if m_star < 0.5:
                        m_rembar = m_star
                        if ((m_core < 0.)or(m_He_core < 0.)):
                            warnings.warn('Invalid co/He core masses! '
                                          'Setting m_WD=m_star!')
                        else:
                            warnings.warn('co/He core masses are zero! '
                                          'Setting m_WD=m_star!')
                    else:
                        raise ValueError('Invalid co/He core masses!')
                f_fb = 1.0  # no SN the no kick is assumed
                state = "WD"

                return m_rembar, f_fb, state, SN_type

            elif (m_core >= min_M_CO_ECSN) and (m_core <= max_M_CO_ECSN):
                SN_type = "ECSN"
            elif m_core > max_M_CO_ECSN:
                SN_type = "CCSN"
            else:
                raise ValueError(
                    "The SN step was applied for an on object outside the "
                    "domain of electron-capture SN and Fe core-collapse SN."
                )

        elif self.ECSN == 'Podsiadlowksi+04':
            # Limits on He core mass progenitors of ECSN, default on cosmic
            min_M_He_ECSN = 1.4  # Msun from Podsiadlowksi+2004
            max_M_He_ECSN = 2.5  # Msun from Podsiadlowksi+2004

            if m_He_core < min_M_He_ECSN:
                # The birth of a white dwarf is assumed
                SN_type = "WD"

                if m_core > 0.:
                    # co_core_mass, note there will be no kick
                    m_rembar = m_core
                elif m_He_core > 0.:
                    m_rembar = m_He_core
                else:
                    # this is catching H-rich_non_burning stars
                    if m_star < 0.5:
                        m_rembar = m_star
                        if ((m_core < 0.)or(m_He_core < 0.)):
                            warnings.warn('Invalid co/He core masses! '
                                          'Setting m_WD=m_star!')
                        else:
                            warnings.warn('co/He core masses are zero! '
                                          'Setting m_WD=m_star!')
                    else:
                        raise ValueError('Invalid co/He core masses!')
                f_fb = 1.0  # no SN the no kick is assumed
                state = "WD"

                return m_rembar, f_fb, state, SN_type

            elif (m_He_core >= min_M_He_ECSN) and (m_He_core <= max_M_He_ECSN):
                SN_type = "ECSN"
            elif m_He_core > max_M_He_ECSN:
                SN_type = "CCSN"
            else:
                raise ValueError(
                    "The SN step was applied for an on object outside the "
                    "domain of electron-capture SN and Fe core-collapse SN."
                )

        elif self.ECSN is None:
            # Here we consider that any CO core mass less that min_M_CO_ECSN
            # will produce a white dwarf
            min_M_CO_ECSN = 1.37  # Msun from Takahashi et al. (2013)
            if m_core < min_M_CO_ECSN:
                # The birth of a white dwarf is assumed
                SN_type = "WD"

                if m_core > 0.:
                    # co_core_mass, note there will be no kick
                    m_rembar = m_core
                elif m_He_core > 0.:
                    m_rembar = m_He_core
                else:
                    # this is catching H-rich_non_burning stars
                    if m_star < 0.5:
                        m_rembar = m_star
                        if ((m_core < 0.)or(m_He_core < 0.)):
                            warnings.warn('Invalid co/He core masses! '
                                          'Setting m_WD=m_star!')
                        else:
                            warnings.warn('co/He core masses are zero! '
                                          'Setting m_WD=m_star!')
                    else:
                        raise ValueError('Invalid co/He core masses!')
                f_fb = 1.0  # no SN the no kick is assumed
                state = "WD"

                return m_rembar, f_fb, state, SN_type

            else:
                SN_type = "CCSN"

        else:
            raise ValueError("The given ECSN prescription is not available.")

        return None, None, None, SN_type

    def compute_m_rembar(self, star, m_PISN):
        """Compute supernova remnant barionic mass.

        We follow the selected electron-capture and core-collapse mechanisms
        to get the remnant baryonic mass.

        Parameters
        ----------
        star : object
            Star object containing the star properties.

        m_PISN : double
            Maximum stellar mass in M_sun after the PPISN/PISN prescription.

        Returns
        -------
        m_rembar : double
            Barioninc mass of the remnant after the supernova in M_sun. This
            quantity does NOT take into account any neutrino lost, this will be
            taken into account in collapse_star().
        f_fb : double
            Mass fraction falling back onto the compact object created in the
            supernova. The maximum value is 1 and means that all the barionic
            mass is collapsing to form the compact object.
        state : string
            Finall state of the stellar remnant after the supernova.

        """
        if star.state in STAR_STATES_CC:
            m_star = star.mass  # M_sun
            m_core = star.co_core_mass  # M_sun
            m_He_core = star.he_core_mass  # M_sun
        elif star.state_history[-1] in STAR_STATES_CC:
            m_star = star.mass_history[-1]  # M_sun
            m_core = star.co_core_mass_history[-1]  # M_sun
            m_He_core = star.he_core_mass_history[-1]  # M_sun
        else:
            raise ValueError(
                "There are no informations in the evolutionary history"
                "about STAR_STATES_CC."
            )
        if m_core is None or np.isnan(m_core):
            # This should not happen
            raise ValueError("The CO core mass is not correct! CO core = {}".
                             format(m_core))

        # define the collapsing stellar mass: either the H or He core mass
        if self.conserve_hydrogen_envelope:
            m_collapsing = m_star
        else:
            m_collapsing = m_He_core

        m_rembar, f_fb, state, star.SN_type = self.check_SN_type(
            m_core=m_core, m_He_core=m_He_core, m_star=m_star)

        if star.SN_type == "WD":
            return m_rembar, f_fb, state

        # Eq. 15-17 from Fryer, C. L., Belczynski, K., Wiktorowicz,
        # G., Dominik, M., Kalogera, V., & Holz, D. E. (2012), ApJ, 749(1), 91.
        if self.mechanism == self.Fryer12_rapid:
            # Mass of the proto-remnant as Giacobbo N., Mapelli M., 2020, ApJ, 891, 141
            m_proto = 1.1

            if star.SN_type == "ECSN":
                if self.ECSN == 'Podsiadlowksi+04':
                    m_proto = 1.38
                else:
                    m_proto = m_core
                m_fb = 0.0  # as in Giacobbo & Mapelli 2020 for ECSN
                f_fb = 0.0
            elif m_core < 2.5:
                m_fb = 0.2
                f_fb = m_fb / (m_collapsing - m_proto)
            elif m_core >= 2.5 and m_core < 6.0:
                m_fb = 0.286 * m_core - 0.514
                f_fb = m_fb / (m_collapsing - m_proto)
            elif m_core >= 6.0 and m_core < 7.0:
                f_fb = 1.0
                m_fb = f_fb * (m_collapsing - m_proto)
            elif m_core >= 7.0 and m_core < 11.0:
                a = 0.25 - 1.275 / (m_collapsing - m_proto)
                b = -11.0 * a + 1.0
                f_fb = a * m_core + b
                m_fb = f_fb * (m_collapsing - m_proto)
            elif m_core >= 11.0:
                f_fb = 1.0
                m_fb = f_fb * (m_collapsing - m_proto)
            m_rembar = m_proto + m_fb
            state = None

        # Eq. 17-20, from Fryer, C. L., Belczynski, K., Wiktorowicz,
        # G., Dominik, M., Kalogera, V., & Holz, D. E. (2012), ApJ, 749(1), 91.
        elif self.mechanism == self.Fryer12_delayed:
            if m_core < 3.5:
                m_proto = 1.2
            elif m_core >= 3.5 and m_core < 6.0:
                m_proto = 1.3
            elif m_core >= 6 and m_core < 11.0:
                m_proto = 1.4
            elif m_core >= 11.0:
                m_proto = 1.6

            if star.SN_type == "ECSN":
                if self.ECSN == 'Podsiadlowksi+04':
                    m_proto = 1.38
                else:
                    m_proto = m_core
                m_fb = 0.0  # as in Giacobbo & Mapelli 2020 for ECSN
                f_fb = 0.0
            elif m_core < 2.5:
                m_fb = 0.2
                f_fb = m_fb / (m_collapsing - m_proto)
            elif m_core >= 2.5 and m_core < 3.5:
                m_fb = 0.5 * m_core - 1.05
                f_fb = m_fb / (m_collapsing - m_proto)
            elif m_core >= 3.5 and m_core < 11.0:
                a = 0.133 - 0.093 / (m_collapsing - m_proto)
                b = -11.0 * a + 1.0
                f_fb = a * m_core + b
                m_fb = f_fb * (m_collapsing - m_proto)
            elif m_core > 11.0:
                f_fb = 1.0
                m_fb = f_fb * (m_collapsing - m_proto)
            m_rembar = m_proto + m_fb
            state = None

        # direct collapse and f_fb = 1. (no kicks)
        elif self.mechanism == self.direct_collapse:
            m_rembar = m_collapsing
            f_fb = 1.0
            state = None

        # Collapse prescription from the results of
        # Sukhbold, T., Ertl, T., Woosley, S. E., Brown, J. M., & Janka, H. T. (2016). 821(1), 38.
        elif self.mechanism == self.Sukhbold16_engines:

            if star.SN_type == "ECSN":
                if self.ECSN == 'Podsiadlowksi+04':
                    m_proto = 1.38
                else:
                    m_proto = m_core
                m_fb = 0.0
                f_fb = 0.0
                m_rembar = m_proto + m_fb
                state = 'NS'
            else:
                m_rembar, f_fb, state = self.Sukhbold_corecollapse_engine(star,
                                                self.conserve_hydrogen_envelope)

        # Collapse prescription from the results of
        # Couch, S. M., Warren, M. L., & O’Connor, E. P. 2020, ApJ, 890, 127
        elif self.mechanism == self.Couch20_engines:

            if star.SN_type == "ECSN":
                if self.ECSN == 'Podsiadlowksi+04':
                    m_proto = 1.38
                else:
                    m_proto = m_core
                m_fb = 0.0
                f_fb = 0.0
                m_rembar = m_proto + m_fb
                state = 'NS'
            else:
                m_rembar, f_fb, state = self.Couch_corecollapse_engine(star,
                                                self.conserve_hydrogen_envelope)

        elif self.mechanism == self.Patton20_engines:
            if star.SN_type == "ECSN":
                if self.ECSN == 'Podsiadlowksi+04':
                    m_proto = 1.38
                else:
                    m_proto = m_core
                f_fb = 0.0
                m_fb = 0.0
                m_rembar = m_proto + m_fb
                state = 'NS'
            else:
                m_rembar, f_fb, state = self.Patton20_corecollapse(star,
                                                self.engine,
                                                self.conserve_hydrogen_envelope)
        else:
            raise ValueError("Mechanism %s not supported." % self.mechanism)

        # check PISN
        if m_PISN is not None:
            if pd.isna(m_PISN):
                m_rembar = m_PISN
                star.SN_type = "PISN"
            elif m_rembar > m_PISN:
                m_rembar = m_PISN
                star.SN_type = "PPISN"

        return m_rembar, f_fb, state

    def orbital_kick(self, binary):
        """Do the orbital kick.

        This function computes the supernova step of the binary object [1]_,
        [2]_. It checks which binary_state reached the core collapse flag,
        either CC1 or CC2, and runs the step accordingly updating the binary
        object.

        Geometry:
        We work in a right-handed coordinate system. The collapsing helium star,
        here M_he_star, lies on the origin. The companion, here M_companion,
        lies on the negative X axis at rest. The relative velocity of the M_he_star
        with respect to M_companion lies in the X-Y plane, with vY>0.
        The orbital angular momentum vector is in Z direction, which completes the
        right-handed coordinate system.

        psi:
            The angle in the orbital plane between the X axis and the pre-core
            collapse relative velocity. (psi = pi/2 points in Y direction)

        theta :
            The polar angle between the kick velocity and the pre-core collapse
            relative velocity of the M_he_star with respect to M_companion.

        phi :
            The corresponding azimuthal angle such that phi=0 is on the Z axis.

        tilt :
            The angle between pre- and post- supernova orbital angular momentum vectors


        Parameters
        ----------
        binary : object
            Binary object containing the binary properties and the two star
            objects.

        References
        ----------
        .. [1] Kalogera, V. 1996, ApJ, 471, 352

        .. [2] Wong, T.-W., Valsecchi, F., Fragos, T., & Kalogera, V. 2012, ApJ, 747, 111

        """
        # Check that the binary_state is calling correctly the SN_step
        if binary.event != "CC1" and binary.event != "CC2":
            raise ValueError("Something went wrong: invalid call of supernova step!")

        if binary.event == "CC1":
            if binary.star_1.SN_type == "WD":
                # compute the new separaiton prior to reseting the binary prop.
                new_separation = separation_evol_wind_loss(
                    binary.star_1.mass, binary.star_1.mass_history[-1],
                    binary.star_2.mass, binary.separation)
                new_orbital_period = orbital_period_from_separation(
                    new_separation, binary.star_1.mass, binary.star_2.mass
                )
                for key in BINARYPROPERTIES:
                    if key not in ['V_sys', 'nearest_neighbour_distance', 'state']:
                        setattr(binary, key, None)
                    # if key is 'nearest_neighbour_distance':
                    #     setattr(binary, key, ['None', 'None', 'None'])
                binary.separation = new_separation
                if binary.state != "disrupted" and binary.state != "initially_single_star" and binary.state != "merged":
                    binary.state = "detached"

                binary.event = None
                binary.time = binary.time_history[-1]
                binary.eccentricity = binary.eccentricity_history[-1]
                # TODO: in feature we will make the orbital period a callable
                # property
                binary.orbital_period = new_orbital_period
                binary.mass_transfer_case = 'None'
                return

            # load relevant data
            # star1 has already collapsed into a compact object, look in the
            # history of the star to find the properties before supernova
            if binary.star_1.state_history[-1] in STAR_STATES_CC:
                M_he_star = binary.star_1.mass_history[-1]
                # Mcore = binary.star_1.co_core_mass_history[-1]
            else:
                raise ValueError(
                    "There are no informations in the evolutionary history "
                    "about STAR_STATES_CC."
                )
            M_compact_object = binary.star_1.mass
            M_companion = binary.star_2.mass

            # check if a kick is passed, otherwise generate it
            if not binary.star_1.natal_kick_array[0] is None:
                Vkick = binary.star_1.natal_kick_array[0]
            else:
                # Draw a random orbital kick
                # Vkick is the kick velocity with components Vkx, Vky, Vkz in
                # the above coordinate system

                if binary.star_1.SN_type == "ECSN":
                    # Kick for electron-capture SN
                    Vkick = self.generate_kick(
                        star=binary.star_1, sigma=self.sigma_kick_ECSN
                    )
                elif ((binary.star_1.SN_type == "CCSN")
                      or (binary.star_1.SN_type == "PPISN")
                      or (binary.star_1.SN_type == "PISN")):
                    if binary.star_1.state == 'NS':
                        sigma = self.sigma_kick_CCSN_NS
                    elif binary.star_1.state == 'BH':
                        sigma = self.sigma_kick_CCSN_BH
                    else:
                        raise ValueError("CCSN/PPISN/PISN only for NS/BH.")
                    # Kick for core-collapse SN
                    Vkick = self.generate_kick(
                        star=binary.star_1, sigma=sigma
                    )
                elif binary.star_1.SN_type == "WD":
                    # Kick for white dwarfs (allways f_fb = 1 => Vkick = 0)
                    Vkick = 0.0
                else:
                    raise ValueError("The SN type is not ECSN neither CCSN.")

                binary.star_1.natal_kick_array[0] = Vkick

            if not binary.star_1.natal_kick_array[1] is None:
                phi = binary.star_1.natal_kick_array[1]
            else:
                phi = np.random.uniform(0, 2 * np.pi)
                binary.star_1.natal_kick_array[1] = phi

            if not binary.star_1.natal_kick_array[2] is None:
                cos_theta = np.cos(binary.star_1.natal_kick_array[2])
            else:
                cos_theta = np.random.uniform(-1, 1)
                binary.star_1.natal_kick_array[2] = np.arccos(cos_theta)

            # generate random point in the orbit where the kick happens
            if not binary.star_1.natal_kick_array[3] is None:
                mean_anomaly = binary.star_1.natal_kick_array[3]
                # check that ONLY one value is passed and is of type float
                if not isinstance(mean_anomaly, float):
                    raise ValueError(
                        "mean_anomaly must be a single float value.")
            else:
                mean_anomaly = np.random.uniform(0, 2 * np.pi)
                binary.star_1.natal_kick_array[3] = mean_anomaly

        elif binary.event == "CC2":
            if binary.star_2.SN_type == "WD":
                # compute new properties before resting existing binary prop.
                new_separation = separation_evol_wind_loss(
                    binary.star_2.mass, binary.star_2.mass_history[-1],
                    binary.star_1.mass, binary.separation)
                new_orbital_period = orbital_period_from_separation(
                    new_separation, binary.star_1.mass, binary.star_2.mass
                )
                for key in BINARYPROPERTIES:
                    if key not in ['V_sys', 'nearest_neighbour_distance', 'state']:
                        setattr(binary, key, None)
                    # if key is 'nearest_neighbour_distance':
                    #     setattr(binary, key, ['None', 'None', 'None'])
                binary.separation = new_separation
                if binary.state != "disrupted" and binary.state != "initially_single_star" and binary.state != "merged":
                    binary.state = "detached"
                binary.event = None
                binary.time = binary.time_history[-1]
                binary.eccentricity = binary.eccentricity_history[-1]
                # TODO: is the following to be noted?
                # in future we will make the orbital period a callable property
                binary.orbital_period = new_orbital_period
                binary.mass_transfer_case = 'None'
                return

            # load relevant data
            # star1 has already collapsed into a compact object, look in the
            # history of the star to find the properties before supernova
            if binary.star_2.state_history[-1] in STAR_STATES_CC:
                M_he_star = binary.star_2.mass_history[-1]
                # Mcore = binary.star_2.co_core_mass_history[-1]
            else:
                raise ValueError(
                    "There are no informations in the evolutionary history "
                    "about STAR_STATES_CC."
                )

            M_compact_object = binary.star_2.mass
            M_companion = binary.star_1.mass

            # check if a kick is passed, otherwise generate it
            if not binary.star_2.natal_kick_array[0] is None:
                Vkick = binary.star_2.natal_kick_array[0]
            else:
                # Draw a random orbital kick
                # Vkick is the kick velocity with components Vkx, Vky, Vkz in
                # the above coordinate system

                if binary.star_2.SN_type == "ECSN":
                    # Kick for electron-capture SN
                    Vkick = self.generate_kick(star=binary.star_2,
                                               sigma=self.sigma_kick_ECSN)
                elif ((binary.star_2.SN_type == "CCSN")
                      or (binary.star_2.SN_type == "PPISN")
                      or (binary.star_2.SN_type == "PISN")):
                    if binary.star_2.state == 'NS':
                        sigma = self.sigma_kick_CCSN_NS
                    elif binary.star_2.state == 'BH':
                        sigma = self.sigma_kick_CCSN_BH
                    else:
                        raise ValueError("CCSN/PPISN/PISN only for NS/BH.")
                    # Kick for core-collapse SN
                    Vkick = self.generate_kick(star=binary.star_2, sigma=sigma)
                else:
                    raise ValueError("The SN type is not ECSN neither CCSN.")

                binary.star_2.natal_kick_array[0] = Vkick

            if not binary.star_2.natal_kick_array[1] is None:
                phi = binary.star_2.natal_kick_array[1]
            else:
                phi = np.random.uniform(0, 2 * np.pi)
                binary.star_2.natal_kick_array[1] = phi

            if not binary.star_2.natal_kick_array[2] is None:
                cos_theta = np.cos(binary.star_2.natal_kick_array[2])
            else:
                cos_theta = np.random.uniform(-1, 1)
                binary.star_2.natal_kick_array[2] = np.arccos(cos_theta)

            # generate random point in the orbit where the kick happens
            if not binary.star_2.natal_kick_array[3] is None:
                mean_anomaly = binary.star_2.natal_kick_array[3]
                # check that ONLY one value is passed and is of type float
                if not isinstance(mean_anomaly, float):
                    raise ValueError(
                        "mean_anomaly must be a single float value.")
            else:
                mean_anomaly = np.random.uniform(0, 2 * np.pi)
                binary.star_2.natal_kick_array[3] = mean_anomaly

        # update the orbit
        if binary.state == "disrupted" or binary.state == "initially_single_star" or binary.state == "merged":
            #the binary was already disrupted before the SN

            # update the binary object which was disrupted already before the SN
            for key in BINARYPROPERTIES:
                if key not in  ('nearest_neighbour_distance','state'):
                    setattr(binary, key, None)
            #binary.state = "disrupted"
            binary.event = None
            binary.separation = np.nan
            binary.eccentricity = np.nan
            binary.V_sys = np.array([0, 0, 0])
            binary.time = binary.time_history[-1]
            binary.orbital_period = np.nan
            binary.mass_transfer_case = 'None'
            binary.first_SN_already_occurred = True

        else:

            # The binary exist: flag_binary is True if the binary is not disrupted
            flag_binary = True

            # eccentricity before the SN
            epre = binary.eccentricity
            # the orbital semimajor axis is the orbital separation
            Apre = binary.separation
            # Eq 16, Wong, T.-W., Valsecchi, F., Fragos, T., & Kalogera, V. 2012, ApJ, 747, 111
            # for eccentric anomaly
            E_ma = sp.optimize.brentq(
                lambda x: mean_anomaly - x + epre * np.sin(x), 0, 2 * np.pi
            )
            # Eq 15, Wong, T.-W., Valsecchi, F., Fragos, T., & Kalogera, V. 2012, ApJ, 747, 111
            # orbital separation at the time of the exlosion
            rpre = Apre * (1.0 - epre * np.cos(E_ma))

            true_anomaly = 2 * np.arctan(
                np.sqrt((1 + epre) / (1 - epre)) * np.tan(E_ma / 2)
            )
            # load constants in CGS
            G = const.standard_cgrav

            # Convert inputs to CGS
            M_he_star = M_he_star * const.Msun
            M_companion = M_companion * const.Msun
            M_compact_object = M_compact_object * const.Msun
            Apre = Apre * const.Rsun
            Vkick = Vkick * const.km2cm
            rpre = rpre * const.Rsun
            Mtot_pre = M_he_star + M_companion
            Mtot_post = M_compact_object + M_companion

            # get useful quantity
            sin_theta = np.sqrt(1 - (cos_theta ** 2))

            # Eq 1, in Kalogera, V. 1996, ApJ, 471, 352
            # extended to Eq 17 in Wong, T.-W., Valsecchi, F., Fragos, T., & Kalogera, V. 2012, ApJ, 747, 111
            # Vr is velocity of preSN He core relative to M_companion, NOT necessarily
            # in the direction of the Y axis if eccentric
            # Eq from conservation of energy
            Vr = np.sqrt(G * (Mtot_pre) * (2.0 / rpre - 1.0 / Apre))

            # Eq 18, Wong, T.-W., Valsecchi, F., Fragos, T., & Kalogera, V. 2012, ApJ, 747, 111
            # psi is the polar angle of the position vector of the CO with respect
            # to its pre-SN orbital velocity in the companions frame. i.e. angle between Vr and X axis
            # If epre = 0, sin_psi should be 1
            # Eq from setting specific angular momentum r X Vr = sqrt(G*M*A*(1-e**2))
            sin_psi = np.round(
                np.sqrt(G * (Mtot_pre) * (1 - epre ** 2) * Apre)
                / (rpre * Vr), 5)
            cos_psi = np.sqrt(1 - sin_psi ** 2)
            # Allow for -cos_psi (Vr in the -X, +Y quadrant)
            if E_ma > np.pi: cos_psi *= -1

            # Eq 3, in Kalogera, V. 1996, ApJ, 471, 352
            # extended to Eq 13, in Wong, T.-W., Valsecchi, F., Fragos, T., & Kalogera, V. 2012, ApJ, 747, 111
            # get the orbital separation post SN
            # Eq from conservation of energy
            Apost = ((2.0 / rpre)
                    - (((Vkick ** 2) + (Vr ** 2) + (2 * (Vkick * cos_theta) * Vr)) / (G * Mtot_post))
                    ) ** -1


                    # get kicks componets in the coordinate system
            Vkx = Vkick * (sin_theta * np.sin(phi) * sin_psi + cos_theta * cos_psi)
            Vky = Vkick * (-sin_theta * np.sin(phi) * cos_psi + cos_theta * sin_psi)
            Vkz = Vkick * sin_theta * np.cos(phi)


            # Eq 4, in Kalogera, V. 1996, ApJ, 471, 352
            # extended to Eq 14 in Wong, T.-W., Valsecchi, F., Fragos, T., & Kalogera, V. 2012, ApJ, 747, 111
            # get the eccentricity post SN
            # Eq from setting specific angular momentum r X Vr = sqrt(G*M*A*(1-e**2))


            x = ((Vkz ** 2 + (Vky + Vr * sin_psi)** 2)
                * rpre ** 2
                / (G * Mtot_post * Apost))

            # catch negative values, i.e. disrupted binaries
            if 1.-x < 0.:
                epost = np.nan
            else:
                epost = np.sqrt(1 - x)

            # Compute COM velocity, VS, post SN
            # VS_pre in COM frame is 0. So VS_post in COM frame is
            # VS_post - VS_pre in our working frame

            VC0x = M_he_star * Vr * cos_psi / Mtot_pre
            VC0y = M_he_star * Vr * sin_psi / Mtot_pre
            VC0z = 0

            VC1x = M_compact_object * (Vkx + Vr * cos_psi) / Mtot_post
            VC1y = M_compact_object * (Vky + Vr * sin_psi) / Mtot_post
            VC1z = M_compact_object * Vkz / Mtot_post


            VSx = VC1x - VC0x
            VSy = VC1y - VC0y
            VSz = VC1z - VC0z


            # V_sys = np.sqrt(VSx ** 2 + VSy ** 2 + VSz ** 2)

            # Calculate the angle between the pre and post-SN orbital angular momentum vectors
            # Lpre || Z axis
            # Lpost || X axis cross the post SN velocity of the compact object
            # cos(tilt) = Lpre dot Lpost / ||Lpre||||Lpost||
            # For epre=0 (sin_psi=1), reduces to Eq 4, in Kalogera, V. 1996, ApJ, 471, 352

            tilt = np.arccos((Vky + Vr * sin_psi) / np.sqrt( Vkz ** 2 + (Vky + Vr * sin_psi) ** 2 ))

            # Track direction of tilt
            if Vkz < 0: tilt *= -1

            def SNCheck(
                M_he_star,
                M_companion,
                M_compact_object,
                rpre,
                Apost,
                epost,
                Vr,
                Vkick,
                cos_theta,
                verbose,
            ):
                """Check that the binary is not disrupted [1]_, [2]_.

                Parameters
                ----------
                M_he_star : double
                    Helium star mass before the SN in g.
                M_companion : double
                    Companion star mass in g.
                M_compact_object : double
                    Compact object mass left  by the SN in g.
                rpre : double
                    Oribtal separation at the time of the exlosion in cm. If the
                    eccentricity pre SN is 0 this correpond to Apre.
                Apost : double
                    Orbital separtion after the SN in cm.
                epost : double
                    Eccentricity after the SN.
                Vr : double
                    Velocity of pre-SN He core relative to M_companion, directed
                    along the positive y axis in cm/s.
                Vkick : double
                    Kick velocity in cm/s.
                cos_theta : double
                    The cosine of the angle between pre- & post-SN orbital planes.

                Returns
                -------
                flag_binary : bool
                    flag_binary is True if the binary is not disrupted.

                References
                ----------
                .. [1] Willems, B., Henninger, M., Levin, T., et al. 2005, ApJ, 625, 324
                .. [2] Kalogera, V. & Lorimer, D.R. 2000, ApJ, 530, 890

                """
                # flag_binary is True if the binary is not disrupted
                flag_binary = True
                Mtot_pre = M_he_star + M_companion
                Mtot_post = M_compact_object + M_companion

                # Define machine precision (we can probaly lower this number)
                err = const.SNcheck_ERR

                # SNflag1: Eq. 21, Willems, B., Henninger, M., Levin, T., et al. 2005, ApJ, 625, 324 (with typo fixed)
                # from Eq. 10, Flannery, B.P. & van den Heuvel, E.P.J. 1975, A&A, 39, 61
                # Continuity demands post-SN orbit to pass through preSN positions.
                # Updated to work for eccentric orbits,
                # see Eq. 15 in Wong, T.-W., Valsecchi, F., Fragos, T., & Kalogera, V. 2012, ApJ, 747, 111
                SNflag1 = (1 - epost - rpre / Apost <= err) and (
                    rpre / Apost - (1 + epost) <= err
                )

                # SNflag2: Equations 22-23, Willems, B., Henninger, M., Levin, T., et al. 2005, ApJ, 625, 324
                # (see, e.g., Kalogera, V. & Lorimer, D.R. 2000, ApJ, 530, 890)
                tmp1 = 2 - Mtot_pre / Mtot_post * (Vkick / Vr - 1) ** 2
                tmp2 = 2 - Mtot_pre / Mtot_post * (Vkick / Vr + 1) ** 2
                SNflag2 = ((rpre / Apost - tmp1 < err)
                           and (err > tmp2 - rpre / Apost))

                # SNflag3: check that epost does not exeed 1 or is nan
                if epost >= 1.0 or np.isnan(epost):
                    SNflag3 = False
                else:
                    SNflag3 = True

                SNflags = [SNflag1, SNflag2, SNflag3]

                if verbose:
                    print()
                    print("The orbital checks are:", SNflags)
                    print()
                    print("1. Post-SN orbit must pass through pre-SN positions.")
                    print("2. Lower and upper limits on amount of orbital "
                          "contraction or expansion that can take place for a "
                          "given amount of mass loss and a given magnitude of the "
                          "kick velocity.")
                    print("3. Checks that e_post is not larger than 1 or nan.")

                # check if the supernova is valid and doesn't disrupt the system
                if not all(SNflags):
                    flag_binary = False

                return flag_binary

            # check if the binary is disrupted
            flag_binary = SNCheck(M_he_star, M_companion, M_compact_object, rpre,
                                  Apost, epost, Vr, Vkick, cos_theta,
                                  verbose=self.verbose)

            # update the binary object which was bound at least before the SN
            #Check if this is the first SN
            for key in BINARYPROPERTIES:
                if key not in ['nearest_neighbour_distance','event']:
                    setattr(binary, key, None)
            if flag_binary:
                # update the tilt
                if not binary.first_SN_already_occurred:
                    # update the tilt
                    binary.star_1.spin_orbit_tilt_first_SN = tilt
                    binary.star_2.spin_orbit_tilt_first_SN = tilt
                    binary.true_anomaly_first_SN = true_anomaly
                    binary.first_SN_already_occurred = True
                else:
                    if binary.event == 'CC2':
                        # Assume progenitor has aligned with the preSN orbital angular momentum
                        binary.star_2.spin_orbit_tilt_second_SN = tilt
                        binary.star_1.spin_orbit_tilt_second_SN = self.get_combined_tilt(
                            tilt_1 = binary.star_1.spin_orbit_tilt_first_SN,
                            tilt_2 = tilt,
                            true_anomaly_1 = binary.true_anomaly_first_SN,
                            true_anomaly_2 = true_anomaly
                            )
                        binary.true_anomaly_second_SN = true_anomaly
                    elif binary.event == 'CC1':
                        # Assume progenitor has aligned with the preSN orbital angular momentum
                        binary.star_1.spin_orbit_tilt_second_SN = tilt
                        binary.star_2.spin_orbit_tilt_second_SN = self.get_combined_tilt(
                            tilt_1 = binary.star_1.spin_orbit_tilt_first_SN,
                            tilt_2 = tilt,
                            true_anomaly_1 = binary.true_anomaly_first_SN,
                            true_anomaly_2 = true_anomaly
                            )
                        binary.true_anomaly_second_SN = true_anomaly
                    else:
                        raise ValueError("This should never happen!")

                # compute new orbital period before reseting the binary properties
                binary.state = "detached"
                binary.event = None
                binary.separation = Apost / const.Rsun
                binary.eccentricity = epost
                binary.V_sys = np.array([VSx / const.km2cm, VSy / const.km2cm, VSz
                                         / const.km2cm])
                binary.time = binary.time_history[-1]
                # in future we will make the orbital period a callable property
                new_orbital_period = orbital_period_from_separation(
                    binary.separation, binary.star_1.mass, binary.star_2.mass)
                binary.orbital_period = new_orbital_period
                binary.mass_transfer_case = 'None'

            else:
                # update the tilt
                if not binary.first_SN_already_occurred:
                    binary.star_1.spin_orbit_tilt_first_SN = np.nan
                    binary.star_2.spin_orbit_tilt_first_SN = np.nan
                    binary.first_SN_already_occurred = True
                else:
                    binary.star_1.spin_orbit_tilt_second_SN = np.nan
                    binary.star_2.spin_orbit_tilt_second_SN = np.nan


                binary.state = "disrupted"
                binary.event = None
                binary.separation = np.nan
                binary.eccentricity = np.nan
                binary.V_sys = np.array([0, 0, 0])
                binary.time = binary.time_history[-1]
                binary.orbital_period = np.nan
                binary.mass_transfer_case = 'None'

    """
    ##### Generating the CCSN SN kick of a single star #####
    """

    def generate_kick(self, star, sigma):
        """Draw a kick from a Maxwellian distribution.

        We follow Hobbs G., Lorimer D. R., Lyne A. G., Kramer M., 2005, MNRAS, 360, 974
        and choose sigma = 265 km/s

        We rescale the kicks by 1 - f_fb as in Eq. 21 of Fryer, C. L., Belczynski, K., Wiktorowicz,
        G., Dominik, M., Kalogera, V., & Holz, D. E. (2012), ApJ, 749(1), 91.

        Parameters
        ----------
        star : object
            Star object containing the star properties.

        Returns
        -------
        Vkick : double
            Natal orbital kick in km/s.

        """
        if self.kick_normalisation == 'one_minus_fallback':
            # Normalization from Eq. 21, Fryer, C. L., Belczynski, K., Wiktorowicz,
            # G., Dominik, M., Kalogera, V., & Holz, D. E. (2012), ApJ, 749(1), 91.
            norm = (1.0 - star.f_fb)
        elif self.kick_normalisation == 'one_over_mass':
            if star.state == 'BH':
                norm = 1.4/star.mass
            else:
                norm = 1.0
        elif self.kick_normalisation == 'NS_one_minus_fallback_BH_one':
            if star.state == 'BH':
                norm = 1.
            else:
                # Normalization from Eq. 21, Fryer, C. L., Belczynski, K., Wiktorowicz,
                # G., Dominik, M., Kalogera, V., & Holz, D. E. (2012), ApJ, 749(1), 91.
                norm = (1.0 - star.f_fb)
        elif self.kick_normalisation == 'one':
            norm = 1.
        elif self.kick_normalisation == 'zero':
            norm = 0.
        else:
            raise ValueError('kick_normalisation option not supported!')

        if sigma is not None:
            # f_fb = self.compute_m_rembar(star, None)[1]

            Vkick = norm * sp.stats.maxwell.rvs(loc=0., scale=sigma, size=1)[0]
        else:
            Vkick = 0.0

        return Vkick

    def get_combined_tilt(self,tilt_1, tilt_2, true_anomaly_1, true_anomaly_2):
        """Get the combined spin-orbit-tilt after two supernovae, assuming
        the spin as not realigned with the orbital angular momentum after
        SN1

            Parameters
            ----------
            tilt_1: float
                Angle, in radians, through which the orbital plane was tilted
                by SN1
            tilt_2: float
                Angle, in radians, through which the orbital plane was tilted
                by SN2
            true_anomaly_1: float
                Angle, in radians, of the true anomaly at the moment of SN1
            true_anomaly_2: float
                Angle, in radians, of the true anomaly at the moment of SN2


            Returns
            -------
            combined_tilt: float
                Angle, in radians, between the spin and orbital angular momentum
                after SN2
        """
        z_prime = rotate((1,0,0), tilt_1).dot((0,0,1))
        x_prime = rotate(z_prime, true_anomaly_2-true_anomaly_1).dot((1,0,0))

        cos_tilt = np.dot((0,0,1),rotate(x_prime, tilt_2).dot(z_prime))
        combined_tilt = np.arccos(cos_tilt)
        return combined_tilt

    def C_abundance_for_H_stars(self, CO_core_mass):
        """Get the C abundance for a H-star given it's CO core mass."""
        return 0.20/CO_core_mass + 0.15

    def C_abundance_for_He_stars(self, CO_core_mass):
        """Get the C abundance for a He-star given it's CO core mass."""
        return -0.084 * np.log(CO_core_mass) + 0.4

    def get_CO_core_params(self, star, approximation=False):
        """Get the CO core mass and C abundance at the pre-supernova phase.

        If the two parameters are available in the star's profile, perform the
        Patton&Sukhbold,20 core-collapse.

        If the CO core mass is available but not the C abundance then the
        latter is computed from the formulas at Patton&Sukhbold,20.

        Parameters
            ----------
            star : obj
                Star object of a collapsing star containing the MESA profile.
            approximation : bool
                In case the core masses at he-depletion are not present in the
                star object, compute them from the history default behaviour,
                else (approximation=True) approximate it from the core masses
                at C depletion.

            Returns
            -------
            CO_core_mass : float
                Mass of the CO core at He depletion == C core ignition.

            C_core_abundance : float
                C abundance of the CO core  He depletion == C core ignition.
        """
        if approximation:
            CO_core_mass = star.co_core_mass # at C_depletion, which is assumed to be close to He depletion
            if ("H-rich" in star.state) or ("H-rich" in star.state_history[-1]):
                C_core_abundance = self.C_abundance_for_H_stars(CO_core_mass)
            elif ("stripped_He" in star.state) or ("stripped_He" in star.state_history[-1]):
                C_core_abundance = self.C_abundance_for_He_stars(CO_core_mass)
            else:
                raise ValueError("star.state at CC should contain either 'H-rich' or 'stripped_He' ")
        elif ((star.avg_c_in_c_core_at_He_depletion is not None)
              and (star.co_core_mass_at_He_depletion is not None)):
            C_core_abundance = star.avg_c_in_c_core_at_He_depletion
            CO_core_mass = star.co_core_mass_at_He_depletion
        else:
            calculate_Patton20_values_at_He_depl(star)
            C_core_abundance = star.avg_c_in_c_core_at_He_depletion
            CO_core_mass = star.co_core_mass_at_He_depletion

            if (C_core_abundance is None) or (CO_core_mass is None):
                raise ValueError('The history did not contain core masses at'
                                 f' He depletion! {CO_core_mass}'
                                 f' {C_core_abundance}')

        return CO_core_mass, C_core_abundance

    def get_M4_mu4_Patton20(self, CO_core_mass, C_core_abundance):
        """Get the M4 and mu4 using Patton+20."""
        M4 = self.M4_interpolator.predict([[C_core_abundance, CO_core_mass]])
        mu4 = self.mu4_interpolator.predict([[C_core_abundance, CO_core_mass]])

        return M4, mu4

    def Patton20_corecollapse(self, star, engine, conserve_hydrogen_envelope=False):
        """Compute supernova final remnant mass and fallback fraction.

        It uses the results from [1]_. The prediction for the core-collapse
        outcome is performed using the C core mass and its C abundance.
        The criterion by [2]_ is used to determine the final outcome.

        Parameters
        ----------
            star : obj
                Star object of a collapsing star containing the MESA profile.

        Returns
        -------
        m_rem : double
            Remnant mass of the compact object in M_sun.
        f_fb : double
            Fallback mass of the compact object in M_sun.

        References
        ----------
        .. [1] Patton, R. A., & Sukhbold, T. (2020). MNRAS, 499(2), 2803-2816.

        .. [2] Ertl, T., Janka, H. T., Woosley, S. E., Sukhbold, T., & Ugliano,
            M. (2016). ApJ, 818(2), 124.

        """
        Ertl16_k_parameters = {
            'N20': [0.194, 0.0580],
            'S19.8': [0.274, 0.0470],
            'W15': [0.225, 0.0495],
            'W20': [0.284, 0.0393],
            'W18': [0.283, 0.0430],
            'Ertl2020': [0.182, 0.0608],
        }

        if engine not in Ertl16_k_parameters.keys():
            raise ValueError("Engine " + engine + " is not avaiable for the "
                             "Patton&Sukhbold,20 core-collapse prescription, "
                             "please choose one of the following engines to "
                             "compute the collapse: \n" + "\n".join(
                                list(Ertl16_k_parameters.keys())))
        else:

            CO_core_mass, C_core_abundance = self.get_CO_core_params(
                star, self.approx_at_he_depletion)
            M4, mu4 = self.get_M4_mu4_Patton20(CO_core_mass, C_core_abundance)
            star.M4 = M4
            star.mu4 = mu4
            print("Manos M4,m4 storage", star.M4, star.mu4)
            M4 = M4[0]
            mu4 = mu4[0]

            k1 = Ertl16_k_parameters[engine][0]
            k2 = Ertl16_k_parameters[engine][1]

            if CO_core_mass <= 2.5:
                m_rem = 1.25
                f_fb = 0.0
                state = 'NS'

            elif CO_core_mass >= 10.0:
                # Assuming BH formation by direct collapse
                if conserve_hydrogen_envelope:
                    m_rem = star.mass
                else:
                    m_rem = star.he_core_mass
                f_fb = 1.0
                state = 'BH'

            elif ((k1 * (mu4 * M4) + k2) < mu4):
                # The prediction is a failed explosion
                # Assuming BH formation by direct collapse
                if conserve_hydrogen_envelope:
                    m_rem = star.mass
                else:
                    m_rem = star.he_core_mass
                f_fb = 1.0
                state = 'BH'
            else:
                # The prediction is a succesful explosion
                m_rem = M4
                f_fb = 0.0
                state = 'NS'

        return m_rem, f_fb, state


class Sukhbold16_corecollapse(object):
    """Compute supernova final remnant mass, fallback fraction and CO type.

    This consider the nearest neighboor of the He core mass of the star,
    previous to the collapse. Considering a set of data for which the He core
    mass of the compact object projenitos previous the collapse, the final
    remnant mass and final stellar state of the compact object is known.

    Parameters
    ----------
    engine : string
        Engine for the supernova explosion, from the one where used in [1]_.

    path_engine_dataset : string
        Path to the location of the data on initial and final states
        for each engine described in [1]_

    Returns
    -------
    m_rem : double
        Remnant mass of the compact object in M_sun.
    f_fb : double
        Fallback mass of the compact object in M_sun.
    state : string
        Finall state of the stellar remnant after the supernova.

    References
    ----------
    .. [1] Sukhbold, T., Ertl, T., Woosley, S. E., Brown, J. M., & Janka, H. T.
        (2016). Core-collapse supernovae from 9 to 120 solar masses based on
        neutrino-powered explosions. The Astrophysical Journal, 821(1), 38.

    """

    def __init__(self, engine, path_engine_dataset, verbose):
        """Initialize a Sukhbold16_corecollapse instance."""
        self.engines = ['N20', 'S19.8', 'W15', 'W20', 'W18']
        self.engine = engine
        self.path_engine_dataset = path_engine_dataset
        if self.engine in self.engines:
            # path_engine_dataset = path_to_Sukhbold_datasets
            if verbose:
                print(
                    "Class initialisation, load the train dataset for engine "
                    + self.engine
                    + " ..."
                )

            # Check if interpolation files exist
            filename = os.path.join(path_engine_dataset,
                                    "results_" + self.engine + "_table.csv")
            if not os.path.exists(filename):
                data_download()

            Engine_data = read_csv(filename)

            # Selecting only the neutro-stars and black-holes with no fallback
            Engine_data = Engine_data[
                (Engine_data["stellar_state"] == 13)
                | (
                    (Engine_data["stellar_state"] == 14)
                    & (Engine_data["fallback_mass"] == 0)
                )
            ]

            if verbose:
                print("Training the classifier ...")

            # Classifier to assign the type of the remnant after supernova
            # as a function of the He core mass pre-supernova
            # taking the first nearest neighbor
            n_neighbors = 1
            self.stellar_type_classifier = neighbors.KNeighborsClassifier(
                n_neighbors, weights="distance"
            )
            self.stellar_type_classifier.fit(
                np.array(Engine_data["He_c_mass"]).reshape(
                    (len(Engine_data["He_c_mass"]), 1)
                ),
                Engine_data["stellar_state"],
            )
            if verbose:
                print("Done ...\n")
                print("Training the remnant mass interpolator ...")

            # Interpolator to compute the remnant mass
            # as a function of the He core mass pre-supernova
            # and the stellar type of the remnant.
            NS_rem_mass = np.array(
                Engine_data[Engine_data["stellar_state"] == 13]["Rem_mass"]
            )
            NS_He_prog = np.array(
                Engine_data[Engine_data["stellar_state"] == 13]["He_c_mass"]
            )
            self.mass_NS_interpolator = interp1d(NS_He_prog, NS_rem_mass)

            BH_rem_mass = np.array(
                Engine_data[Engine_data["stellar_state"] == 14]["Rem_mass"]
            )
            BH_He_prog = np.array(
                Engine_data[Engine_data["stellar_state"] == 14]["He_c_mass"]
            )
            self.mass_BH_interpolator = interp1d(BH_He_prog, BH_rem_mass)

            if verbose:
                print("Done ...\n")

            # Gets the neutron-star mass in terms of the He core mass
            # if a succesful explotion is predicted
            def extrapolate1d_NS(value, interpolator):
                x = interpolator.x

                if (value >= np.min(x)) and (value <= np.max(x)):
                    result = interpolator(value)
                elif value < np.min(x):
                    result = interpolator(np.min(x))
                elif value > np.max(x):
                    result = interpolator(np.max(x))

                return result

            # Gets the black-hole mass in terms of the He core mass
            # if a unsuccesful explotion is predicted
            def extrapolate1d_BH(value, interpolator):
                x = interpolator.x

                if (value >= np.min(x)) and (value <= np.max(x)):
                    result = interpolator(value)
                elif value < np.min(x):
                    result = interpolator(np.min(x))
                elif value > np.max(x):
                    result = value

                return result

            self.extrapolate_NS = extrapolate1d_NS
            self.extrapolate_BH = extrapolate1d_BH

        else:
            raise ValueError(
                "Engine " + self.engine + " is not avaiable for the"
                "Sukhbold core collapse prescription, please choose"
                "one of the following engines to compute the collapse: ",
                self.engines,
            )

    def __call__(self, star, conserve_hydrogen_envelope=False):
        """Get the mass, fallback franction and state of the remnant."""
        if star.state in STAR_STATES_CC:
            m_star = star.mass  # M_sun
            # m_core = star.co_core_mass  # M_sun
            m_He_core = star.he_core_mass  # M_sun
        elif star.state_history[-1] in STAR_STATES_CC:
            m_star = star.mass_history[-1]  # M_sun
            # m_core = star.co_core_mass_history[-1]  # M_sun
            m_He_core = star.he_core_mass_history[-1]  # M_sun
        else:
            raise ValueError("There are no informations in the evolutionary "
                             "history about STAR_STATES_CC.")
        k_result = int(self.stellar_type_classifier.predict([[m_He_core]])[0])

        if k_result == 13:
            state = "NS"
        elif k_result == 14:
            state = "BH"
        else:
            state = None

        if state == "BH":
            # Assuming BH formation by direct collapse
            if conserve_hydrogen_envelope:
                m_rem = self.extrapolate_BH(m_star, self.mass_BH_interpolator)
            else:
                m_rem = self.extrapolate_BH(m_He_core, self.mass_BH_interpolator)
            f_fb = 1.
        elif state == "NS":
            m_rem = self.extrapolate_NS(m_He_core, self.mass_NS_interpolator)
            f_fb = 0.
        else:
            raise Exception("Need a NS or BH to apply `Sukhbold16_corecollapse`.")

        return float(m_rem), f_fb, state


class Couch20_corecollapse(object):
    """Compute SN final remnant mass, fallback fraction and stellar state.

    This considers the nearest neighboor of the He core mass of the star,
    previous to the collapse. Considering a set of data for which the He core
    mass of the compact object progenitors before the collapse, the final
    remnant mass and final stellar state of the compact object is known.

    Parameters
    ----------
    engine : string
        Engine for the supernova explosion, from the one where used in [1]_.

    path_engine_dataset : string
        Path to the location of the data on initial and final states
        for each engine described in [1]_

    Returns
    -------
    m_rem : double
        Remnant mass of the compact object in M_sun.
    f_fb : double
        Fallback mass of the compact object in M_sun.
    state : string
        Finall state of the stellar remnant after the supernova.

    Notes
    -----
    We need [2]_ data for their cores.

    References
    ----------

    .. [1] Couch, S. M., Warren, M. L., & O’Connor, E. P. 2020, ApJ, 890, 127
        Simulating Turbulence-aided Neutrino-driven Core-collapse Supernova
        Explosions in One Dimension
    .. [2] Sukhbold, T., Ertl, T., Woosley, S. E., Brown, J. M., & Janka, H. T.
        (2016). Core-collapse supernovae from 9 to 120 solar masses based on
        neutrino-powered explosions. The Astrophysical Journal, 821(1), 38.

    """

    def __init__(self, turbulence_strength, path_engine_dataset, verbose):
        """Initialize a Couch20_corecollapse instance."""
        self.turbulence_strength_options = ["1.0", "1.2", "1.23", "1.25",
                                            "1.27", "1.3", "1.4"]
        self.turbulence_strength = turbulence_strength
        self.path_engine_dataset = path_engine_dataset

        if turbulence_strength in self.turbulence_strength_options:
            # path_engine_dataset = path_to_Sukhbold_datasets
            if verbose:
                print(
                    "Class initialisation, load the train dataset for engine "
                    + self.engine
                    + " ..."
                )


            # Check if interpolation files exist
            filename = os.path.join(path_to_Couch_datasets,
                                    'explDatsSTIR2.json')
            if not os.path.exists(filename):
                data_download()

            Couch_data_file = open(filename)
            # Couch_data = json.loads(Couch_data_file)
            Couch_data = json.load(Couch_data_file)
            Couch_data_file.close()
            Couch_data = Couch_data[turbulence_strength]
            # breakpoint()
            # #names = ['MZAMS',['rmax','texp','Eexp']]
            # my_names = ['MZAMS','Eexp']
            # my_formats = ['f8','f8']
            # dtype = dict(names = my_names, formats= my_formats)
            # #dt = np.dtype()
            # #Couch_data_ar = np.array(list(Couch_data.items()), dtype=dtype)
            #
            # #Couch_data_ar = np.array(Couch_data[turbulence_strength])
            # #Couch_data_ar = [np.array(Couch_data[MZAMS])
            #                   for MZAMS in Couch_data]
            # #Couch_data_ar = np.array([(MZAMS,rest["Eexp"]) for (MZAMS,rest)
            #                            in Couch_data.items()], dtype=dtype)
            Couch_MZAMS = []
            Couch_Eexp = []
            Couch_state = []
            for MZAMS, rest in Couch_data.items():
                Couch_MZAMS.append(float(MZAMS))
                Couch_Eexp.append(rest["Eexp"])
                if rest["Eexp"] == 0.0:
                    Couch_state.append(int(14))     # BH
                else:
                    Couch_state.append(int(13))     # NS
            # Couch_MZAMS = np.array(Couch_MZAMS, dtype=dict(
            #     names=my_names[0], formats=my_formats[0]))
            # Couch_Eexp = np.array(Couch_Eexp,dtype=dict(
            #     names=my_names[1], formats= my_formats[1]))
            # Couch_data_ar = np.array()
            # breakpoint()
            # print(Couch_data)

            # we need Sukhbold data for  their cores
            Sukhbold_data = read_csv(
                # path_to_Sukhbold_datasets + "results_N20_table.csv"
                path_to_Couch_datasets + "Sukhbold_Mzams_He_c_core.csv",
                usecols=[0, 1])

            MZAMS = Sukhbold_data["Mzams"]
            He_core_mass = Sukhbold_data["He_c_mass"]
            self.MZAMS_He_core_mass_Sukhbold_interpolator = interp1d(
                MZAMS, He_core_mass)
            # def MZAMS_He_core_mass_Sukhbold_interpolator(MZams):
            #     return Sukhbold_data[Sukhbold_data["Mzams"]==MZams][
            #         "He_c_mass"]

            # # Classifier to assign the He core mass of Sukhbold
            # # as a function of the MZAMS
            # # taking the first nearest neighbor
            # n_neighbors = 1
            # self.stellar_ZAMS_classifier = neighbors.KNeighborsClassifier(
            #     n_neighbors, weights="distance"
            # )
            # self.stellar_ZAMS_classifier.fit(
            #     np.array(Sukhbold_data["Mzams"]).reshape(
            #         (len(Sukhbold_data["Mzams"]), 1)
            #     ),
            #     Sukhbold_data["He_c_mass"],
            # )
            # MZAMS = np.array(
            #     Sukhbold_data["Mzams"]
            # )
            # He_c_mass = np.array(
            #     Sukhbold_data["He_c_mass"]
            # )
            # self.MZAMS_He_core_mass_Sukhbold_interpolator = interp1d(MZAMS,
            #       He_c_mass)

            Couch_He_c_mass = self.MZAMS_He_core_mass_Sukhbold_interpolator(
                Couch_MZAMS)

            if verbose:
                print("Training Couch+20 data...\n")

            # Classifier to assign the type of the remnant after supernova
            # as a function of the He core mass pre-supernova of Sukhbold 2016
            # taking the first nearest neighbor
            n_neighbors = 1
            self.stellar_type_classifier = neighbors.KNeighborsClassifier(
                n_neighbors, weights="distance"
            )
            self.stellar_type_classifier.fit(
                np.array(Couch_He_c_mass).reshape(
                    (len(Couch_He_c_mass), 1)
                ),
                np.array(Couch_state),
            )

            '''
            breakpoint()

            print("Training the remnant mass interpolator ...")

            # Interpolator to compute the remnant mass
            # as a function of the He core mass pre-supernova
            # and the stellar type of the remnant.
            NS_rem_mass = np.array(
                Engine_data[Engine_data["stellar_state"] == 13]["Rem_mass"]
            )
            NS_He_prog = np.array(
                Engine_data[Engine_data["stellar_state"] == 13]["He_c_mass"]
            )
            self.mass_NS_interpolator = interp1d(NS_He_prog, NS_rem_mass)

            BH_rem_mass = np.array(
                Engine_data[Engine_data["stellar_state"] == 14]["Rem_mass"]
            )
            BH_He_prog = np.array(
                Engine_data[Engine_data["stellar_state"] == 14]["He_c_mass"]
            )
            self.mass_BH_interpolator = interp1d(BH_He_prog, BH_rem_mass)

            if verbose:
                print("Done ...\n")

            # Gets the neutron-star mass in terms of the He core mass
            # if a succesful explotion is predicted
            def extrapolate1d_NS(value, interpolator):
                x = interpolator.x

                if (value >= np.min(x)) and (value <= np.max(x)):
                    result = interpolator(value)
                elif value < np.min(x):
                    result = interpolator(np.min(x))
                elif value > np.max(x):
                    result = interpolator(np.max(x))

                return result

            # Gets the black-hole mass in terms of the He core mass
            # if a unsuccesful explotion is predicted
            def extrapolate1d_BH(value, interpolator):
                x = interpolator.x

                if (value >= np.min(x)) and (value <= np.max(x)):
                    result = interpolator(value)
                elif value < np.min(x):
                    result = interpolator(np.min(x))
                elif value > np.max(x):
                    result = value

                return result

            self.extrapolate_NS = extrapolate1d_NS
            self.extrapolate_BH = extrapolate1d_BH
            '''

        else:
            raise ValueError(
                "Turbulence strength " + self.turbulence_strength + " is not "
                "available for the Couch core collapse prescription, please "
                "choose one of the following engines to compute the collapse:",
                self.turbulence_strength_options)

    def __call__(self, star, conserve_hydrogen_envelope=False):
        """Get the mass, fallback fraction and state of the remnant."""
        if star.state in STAR_STATES_CC:
            m_star = star.mass                          # M_sun
            # m_core = star.co_core_mass                  # M_sun
            m_He_core = star.he_core_mass               # M_sun
        elif star.state_history[-1] in STAR_STATES_CC:
            m_star = star.mass_history[-1]              # M_sun
            # m_core = star.co_core_mass_history[-1]      # M_sun
            m_He_core = star.he_core_mass_history[-1]   # M_sun
        else:
            raise ValueError("There are no informations in the evolutionary "
                             "history about STAR_STATES_CC.")
        # single_star_equivalent_ZAMS = \
        #     self.stellar_ZAMS_classifier.predict([[m_He_core]])[0]

        k_result = int(self.stellar_type_classifier.predict([[m_He_core]])[0])

        if k_result == 13:
            state = "NS"
        elif k_result == 14:
            state = "BH"
        else:
            state = None

        if state == "BH":
            # Assuming BH formation by direct collapse
            if conserve_hydrogen_envelope:
                m_rem = m_star
            else:
                m_rem = m_He_core
            # m_rem = self.extrapolate_BH(m_He_core, self.mass_BH_interpolator)
            # TODO: We need to contact Couch et al. to get the remnant masses
            # f_fb = m_rem / m_He_core
            f_fb = 1.
        elif state == "NS":
            # TODO We need to contact Couch et al. to get the remnant masses
            m_rem = 1.4
            # f_fb = m_rem / m_He_core
            f_fb = 0.
        else:
            raise Exception("Need a NS or BH to apply `Sukhbold16_corecollapse`.")

        return float(m_rem), f_fb, state



def check_SN_CO_match(star):
    '''Check if the SN type matches the stellar state of the given star.

    Parameters
    ----------
    star : SingleStar object
        Star object containing the star properties.

    Returns
    -------
    correct_SN_type : bool
        True if the SN type matches the stellar state of the star.
    '''
    # TODO: remove star.state == PISN, because PISN shouldn't be a stellar state
    if star.state == 'PISN':
        star.state = 'massless_remnant'
    correct_SN_type = True
    if star.state == 'WD' and star.SN_type != "WD":
        correct_SN_type = False
    elif (star.state == "NS") and \
            (star.SN_type != 'ECSN' and
            star.SN_type != "CCSN"):
        correct_SN_type = False
    elif (star.state =="BH") and \
            (star.SN_type != "CCSN" and
            star.SN_type != 'PPISN'):
        correct_SN_type = False
    elif (star.state == "massless_remnant" and star.SN_type != 'PISN'):
        correct_SN_type = False
    return correct_SN_type
