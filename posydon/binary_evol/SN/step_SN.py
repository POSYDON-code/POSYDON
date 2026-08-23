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
    "Max Briel <max.briel@unige.ch>",
    "Seth Gossage <seth.gossage@northwestern.edu>",
    "Dimitris Souropanis <dsouropanis@ia.forth.gr>"
]

__credits__ = [
    "Michael Zevin <michael.zevin@ligo.org>",
    "Chase Kimball <charles.kimball@ligo.org",
    "Sam Imperato <samuelimperato2022@u.northwestern.edu>",
]


import copy
import json
import os

import numpy as np
import pandas as pd
import scipy as sp
from pandas import read_csv
from sklearn import neighbors

import posydon.utils.constants as const
from posydon.binary_evol.binarystar import BINARYPROPERTIES
from posydon.binary_evol.flow_chart import (
    STAR_STATES_C_DEPLETED,
    STAR_STATES_CC,
    STAR_STATES_CO,
)
from posydon.binary_evol.singlestar import (
    STARPROPERTIES,
    convert_star_to_massless_remnant,
)
from posydon.binary_evol.SN.profile_collapse import (
    do_core_collapse_BH,
    get_ejecta_element_mass_at_collapse,
)
from posydon.config import PATH_TO_POSYDON_DATA
from posydon.grids.SN_MODELS import DEFAULT_SN_MODEL, get_SN_MODEL_NAME
from posydon.utils.common_functions import (
    CO_radius,
    calculate_Patton20_values_at_He_depl,
    inspiral_timescale_from_separation,
    is_number,
    orbital_period_from_separation,
    rotate,
    separation_evol_wind_loss,
    set_binary_to_failed,
)
from posydon.utils.data_download import data_download
from posydon.utils.interpolators import interp1d
from posydon.utils.limits_thresholds import (
    NEUTRINO_MASS_LOSS_UPPER_LIMIT,
    STATE_NS_STARMASS_UPPER_LIMIT,
    THRESHOLD_CENTRAL_ABUNDANCE,
)
from posydon.utils.posydonerror import ModelError
from posydon.utils.posydonwarning import Pwarn

import kicks


path_to_Sukhbold_datasets = os.path.join(PATH_TO_POSYDON_DATA,
                                         "Sukhbold+16/")

path_to_Patton_datasets = os.path.join(PATH_TO_POSYDON_DATA,
                                       "Patton+Sukhbold20/")

path_to_Couch_datasets = os.path.join(PATH_TO_POSYDON_DATA,
                                      "Couch+2020/")

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

        * 'Maltsev+25-engine': Uses the results from [8]_
        to describe the collapse of the star

    engine : str
        Engine used for supernova remnant outcome propierties for the
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

        - 'Podsiadlowski+04': Determines the electron capture supernova in
        terms of the He core mass at pre-supernova, taking limits from [7]_.
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

    .. [7] Podsiadlowski, P., Langer, N., Poelarends, A. J. T., Rappaport, S.,
        Heger, A., and Pfahl, E. 2004, ApJ, 612, 1044. The Effects of Binary
        Evolution on the Dynamics of Core Collapse and Neutron Star Kicks

    .. [8] K. Maltsev, F.R.N. Schneider, I. Mandel, B. Mueller, A. Heger, F.K. Roepke,
         E. Laplace, 2025,  A&A, 700, A20. Explodability criteria for the neutrino-driven
         supernova mechanism
    """
    DEFAULT_KWARGS = {
        # kick physics
        "kick": True,
        "kick_normalisation": 'one_over_mass',
        "kick_prescription": 'maxwellian',
        "sigma_kick_CCSN_NS": 265.0,
        "mean_kick_CCSN_NS": None,
        "sigma_kick_CCSN_BH": 265.0,
        "mean_kick_CCSN_BH": None,
        "sigma_kick_ECSN": 20.0,
        "mean_kick_ECSN": None,
        # other
        "RNG": None,
        "verbose": False
    }
    # add core collapse physics
    DEFAULT_KWARGS.update(DEFAULT_SN_MODEL)


    def __init__(self, **kwargs):
        """Initialize a StepSN instance."""
        # read kwargs to initialize the class
        if kwargs:
            for key in kwargs:
                if key not in self.DEFAULT_KWARGS:
                    raise ValueError(key + " is not a valid parameter name!")
            for varname in self.DEFAULT_KWARGS:
                setattr(self, varname, kwargs.get(varname, self.DEFAULT_KWARGS[varname]))
            self.RNG = kwargs.get("RNG")
            if self.RNG is None:
                self.RNG = np.random.default_rng()

        else:
            for varname in self.DEFAULT_KWARGS:
                setattr(self, varname, self.DEFAULT_KWARGS[varname])

        # backward compatibility for kick
        if (self.kick_normalisation == 'asym_ej'
            or self.kick_normalisation == 'linear'):
            Pwarn("kick_normalisation 'asym_ej' and 'linear' are "
                "deprecated, use kick_prescription instead. Setting "
                "kick normalization to unity.",
                "DeprecationWarning")
            self.kick_prescription = self.kick_normalisation
            self.kick_normalisation = 'one'

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
        self.Maltsev25_engines = "Maltsev+25-engine"



        self.mechanisms = [
            self.Fryer12_rapid,
            self.Fryer12_delayed,
            self.direct_collapse,
            self.direct_collapse_hecore,
            self.Sukhbold16_engines,
            self.Patton20_engines,
            self.Couch20_engines,
            self.Maltsev25_engines
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

            elif self.mechanism in (self.Patton20_engines, self.Maltsev25_engines):
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
                        data_download(set_name='auxiliary')

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
                    print('Loading the train dataset for engine mu4, M4, Xi, and sc ...')
                CO_core_params_mu4, mu4_target = format_data_Patton20(
                    'Kepler_mu4_table.dat')
                CO_core_params_M4, M4_target = format_data_Patton20(
                    'Kepler_M4_table.dat')
                CO_core_params_Xi, Xi_target = format_data_Patton20(
                    'Kepler_Xi_table.dat')
                CO_core_params_sc, sc_target = format_data_Patton20(
                    'Kepler_sc_table.dat')

                n_neighbors = 5

                if self.verbose:
                    print('Training the classifier ...')
                self.M4_interpolator = neighbors.KNeighborsRegressor(
                    n_neighbors, weights='distance')
                self.M4_interpolator.fit(CO_core_params_M4, M4_target)

                self.mu4_interpolator = neighbors.KNeighborsRegressor(
                    n_neighbors, weights='distance')
                self.mu4_interpolator.fit(CO_core_params_mu4, mu4_target)

                self.Xi_interpolator = neighbors.KNeighborsRegressor(
                    n_neighbors, weights='distance')
                self.Xi_interpolator.fit(CO_core_params_Xi, Xi_target)

                self.sc_interpolator = neighbors.KNeighborsRegressor(
                    n_neighbors, weights='distance')
                self.sc_interpolator.fit(CO_core_params_sc, sc_target)
                if self.verbose:
                    print('Done')
        else:
            raise ValueError("Invalid core-collapse mechanism given.")

    def __repr__(self):
        """Get the string representation of the class and any parameters."""
        return "StepSN:\n" + \
            "\n".join([f"{key} = {getattr(self, key)}" for key in self.__dict__])


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
        if binary.event == "CC1":
            # collapse star
            model_err = self.collapse_star(star=binary.star_1)
            if model_err is not None:
                set_binary_to_failed(binary)
                raise ModelError(model_err)

            self._reset_other_star_properties(star=binary.star_2)
            binary.update_star_states()

        elif binary.event == "CC2":
            # collapse star
            model_err = self.collapse_star(star=binary.star_2)
            if model_err is not None:
                set_binary_to_failed(binary)
                raise ModelError(model_err)

            self._reset_other_star_properties(star=binary.star_1)
            binary.update_star_states()
        else:
            raise ValueError("Something went wrong: "
                             "invalid call of supernova step!")
        # do orbital_kick on the binary object
        if self.kick:
            orbital_kick(
                binary = self.binary,
                sigma_kick_ECSN = self.sigma_kick_ECSN,
                mean_kick_ECSN = self.mean_kick_ECSN,
                sigma_kick_CCSN_NS = self.sigma_kick_CCSN_NS,
                mean_kick_CCSN_NS = self.mean_kick_CCSN_NS,
                sigma_kick_CCSN_BH = self.sigma_kick_CCSN_BH,
                mean_kick_CCSN_BH = self.mean_kick_CCSN_BH,
                RNG = self.RNG,
                verbose = self.verbose,
            )
        else:
            # no kick, but still need to unset the event after CC
            # and update orbital period/separation
            binary.event = None
            new_separation = separation_evol_wind_loss(
                    binary.star_1.mass, binary.star_1.mass_history[-1],
                    binary.star_2.mass, binary.separation)
            new_orbital_period = orbital_period_from_separation(
                    new_separation, binary.star_1.mass, binary.star_2.mass)

            binary.separation = new_separation
            binary.orbital_period = new_orbital_period
            binary.state = "detached"

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
        elif state1 in STAR_STATES_CO and state2 in STAR_STATES_C_DEPLETED:
            binary.event = "CC2"
        elif state1 in STAR_STATES_C_DEPLETED and state2 in STAR_STATES_CO:
            binary.event = "CC1"

        if self.verbose:
            print(f"End of step SN:\n", binary)


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
        error_message : string
            Error message if the collapse fails, else None.

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

            SN_type = ""
            # if no profile is avaiable but interpolation quantities are,
            # use those, else continue with or without profile.
            if self.use_interp_values:
                # find SN_MODEL_NAME corresponding to class variable
                SN_MODEL_NAME_SEL = get_SN_MODEL_NAME(vars(self),
                                                      verbose=self.verbose)

                # check if selected model is supported
                if SN_MODEL_NAME_SEL is None:
                    raise ValueError('Your model assumptions are not'
                                     'supported!')
                elif getattr(star, SN_MODEL_NAME_SEL) is None:
                    # NOTE: this option is needed to do the collapse
                    # for stars evolved with the step_detached or
                    # step_disrupted.
                    # allow to continue with the collapse with profile
                    # or core masses
                    Pwarn(f'{SN_MODEL_NAME_SEL}: The collapsed star was not '
                          'interpolated! If use_profiles or use_core_masses '
                          'is set to True, continue with the collapse.',
                          "InterpolationWarning")
                else:
                    SN_MODEL_properties = getattr(star, SN_MODEL_NAME_SEL)

                    SN_type = self.check_SN_type(m_core=star.co_core_mass,
                                                 m_He_core=star.he_core_mass,
                                                 m_star=star.mass)[3]
                    if self.use_profiles and star.profile is not None:
                        alternative = "Instead use profiles."
                    elif self.use_core_masses:
                        alternative = "Instead use core masses."
                    elif self.allow_spin_None:
                        alternative = "Instead use core mass without spin."
                    else:
                        alternative = ""

                    if SN_MODEL_properties['SN_type'] == "ECSN":
                        # overwrite ECSN in SN MODEL
                        SN_MODEL_properties['SN_type'] = SN_type
                        Pwarn(f"ECSN in SN_MODEL replaced by {SN_type}",
                              "ReplaceValueWarning")

                    if SN_type == "ECSN":
                        # do not use interpolated values for ECSN range instead
                        # behave like use_core_masses=True
                        pass
                    ## star's SN_type mismatches one from the model
                    elif SN_type != SN_MODEL_properties['SN_type']:
                        Pwarn(f"The SN_type does not match the star: {SN_type}"
                              f"!={SN_MODEL_properties['SN_type']}."
                              +alternative, "ApproximationWarning")
                    ## Check if SN_type mismatches the CO_type in the model
                    elif not check_SN_CO_match(SN_MODEL_properties['SN_type'],
                                               SN_MODEL_properties['state']):
                        Pwarn(f"{SN_MODEL_NAME_SEL}: The SN_type does not "
                              "match the predicted CO."+alternative,
                              "ApproximationWarning")
                    ## Check if there is no interpolated remnant mass
                    elif pd.isna(SN_MODEL_properties['mass']):
                        Pwarn(f"There is no interpolated remnant mass."
                              +alternative, "ApproximationWarning")
                    ## Otherwise interpolated values can be used for this SN
                    else:
                        for key, value in SN_MODEL_properties.items():
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

                        # No remnant if a PISN happens
                        if star.SN_type == 'PISN':
                            convert_star_to_massless_remnant(star=star)
                            # the mass is set to None
                            # but an orbital kick is still applied.
                            # Since the mass is set to None, this will lead to
                            # a disruption
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
                if pd.isna(m_rembar):
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
                        # set post-collapse properties/information to store
                        for i in final_BH.keys():
                            setattr(star, i, final_BH[i])

                        # set specific properties manually
                        star.mass = final_BH['M_BH_total']
                        star.spin = final_BH['a_BH_total']
                        star.m_disk_accreted = final_BH['m_disk_accreted']
                        star.m_disk_radiated = final_BH['m_disk_radiated']
                        star.state = "BH"
                    else:
                        star.mass = m_grav
                        star.spin = 0.
                        star.m_disk_accreted = 0.
                        star.m_disk_radiated = 0.
                        star.state = 'NS'
                    star.h1_mass_ej, star.he4_mass_ej = \
                        get_ejecta_element_mass_at_collapse(star,star.mass,verbose=self.verbose)

                elif self.use_core_masses or SN_type == "ECSN":
                    # If the profile is not available the star spin
                    # is used to get the compact object spin
                    star.mass = m_grav
                    if m_grav >= self.max_NS_mass:
                        if SN_type == "ECSN":
                            Pwarn("An ECSN should not form a black hole: "
                                  f"m_grav={m_grav}.",
                                  "InappropriateValueWarning")
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
                    star.h1_mass_ej, star.he4_mass_ej = \
                        np.nan, np.nan

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
                    star.h1_mass_ej, star.he4_mass_ej = \
                        np.nan, np.nan

                else:
                    # This leads to a failed binary
                    for key in STARPROPERTIES:
                        setattr(star, key, None)
                    return "FAILED core collapse!"

            elif self.mechanism in [
                self.Sukhbold16_engines,
                self.Patton20_engines,
                self.Couch20_engines,
                self.Maltsev25_engines
            ]:
                # The final remnant mass and and state
                # is computed by the selected mechanism

                # PISN and PPISN prescription
                m_PISN = self.PISN_prescription(star)

                m_rembar, star.f_fb, state = self.compute_m_rembar(star,
                                                                   m_PISN)
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
                if pd.isna(m_rembar):
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
                        # set post-collapse properties/information to store
                        for i in final_BH.keys():
                            setattr(star, i, final_BH[i])
                        # set specific properties manually
                        star.mass = final_BH['M_BH_total']
                        star.spin = final_BH['a_BH_total']

                        if m_grav != star.mass and self.verbose:
                            print("The star formed a disk during the collapse "
                                  "and lost", round(final_BH['M_BH_total'] - m_rembar, 2),
                                  "M_sun.")

                    elif star.state == "NS":
                        star.mass = m_grav
                        star.m_disk_accreted = 0.0
                        star.m_disk_radiated = 0.0
                        star.spin = 0.0
                    else:
                        # This leads to a failed binary
                        for key in STARPROPERTIES:
                            setattr(star, key, None)
                        return f"FAILED core collapse! (Invalid core state: {state})"

                    star.h1_mass_ej, star.he4_mass_ej = \
                        get_ejecta_element_mass_at_collapse(star,star.mass,verbose=self.verbose)

                elif self.use_core_masses or SN_type == "ECSN":
                    star.mass = m_grav
                    if m_grav >= self.max_NS_mass:
                        if SN_type == "ECSN":
                            Pwarn("An ECSN should not form a black hole: "
                                  f"m_grav={m_grav}.",
                                  "InappropriateValueWarning")
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
                    star.h1_mass_ej, star.he4_mass_ej = \
                        np.nan, np.nan

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
                    star.h1_mass_ej, star.he4_mass_ej = \
                        np.nan, np.nan

                else:
                    # This leads to a failed binary
                    for key in STARPROPERTIES:
                        setattr(star, key, None)
                    return "FAILED core collapse!"

        else:
            # This leads to a failed binary
            return f"The star cannot collapse: star state {state}."

        star.metallicity = star.metallicity_history[-1]

        star.log_R = np.log10(CO_radius(star.mass, star.state))

        for key in STARPROPERTIES:
            if key not in ["state", "mass", "spin", "log_R", "metallicity",
                           "m_disk_accreted", "m_disk_radiated",
                           "co_core_mass"]:
                setattr(star, key, None)

        return

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
            return

        else:
            # perform the PISN prescription in terms of the
            # He core mass at pre-supernova
            m_He_core = star.he_core_mass
            m_CO_core = star.co_core_mass
            m_star = star.mass
            if self.PISN == "Marchant+19":
                if m_He_core >= 31.99 and m_He_core <= 61.10:
                    # this is the 8th-order polynomial fit of table 1
                    # value, see COSMIC paper (Breivik et al. 2020)
                    polyfit = (
                        - 6.29429263e5
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

                elif m_He_core > 61.10 and m_He_core < 124.12:
                    # in Breivik et al. (2020) they qoute the CO core mass
                    # range as 54.48<M_CO-core/Msun<113.29 here, but this
                    # might cause gaps, when switching between core masses,
                    # hence take the He-core masses from table 1 of Marchant
                    # et al. (2019)
                    m_PISN = np.nan

                else:
                    # above the PISN gap we assume direct collapse of the
                    if self.conserve_hydrogen_envelope:
                        m_PISN = m_star
                    else:
                        m_PISN = m_He_core

            elif self.PISN == 'Hendriks+23':
                # Hendriks et al. 2023 PISN prescription
                # 10.1093/mnras/stad2857
                # Shifting PPI and PISN gap
                # works by removing delta_M_PPI from the star
                # and then applying any remnant mass prescription

                delta_M_CO_shift = self.PISN_CO_shift if self.PISN_CO_shift is not None else 0.0
                delta_M_PPI_extra_ML = self.PPI_extra_mass_loss if self.PPI_extra_mass_loss is not None else 0.0

                m_CO_core_PISN_min = 38 + delta_M_CO_shift
                m_CO_core_PISN_max = 114 + delta_M_CO_shift

                if ((m_CO_core >= m_CO_core_PISN_min)
                    and m_CO_core <= m_CO_core_PISN_max):

                    # delta_PPI -> -inf if Z -> 0
                    # limit mass loss to Z = 1e-4 for Z below it.
                    # 1e-4 is the lowest metallicity in the Hendriks et al. 2023
                    if star.metallicity < 1e-4:
                        Z = 1e-4
                    else:
                        Z = star.metallicity
                    # Hendriks et al. 2023 Equation 6
                    # 10.1093/mnras/stad2857
                    delta_M_PPI = (
                        (0.0006 * np.log10(Z * const.Zsun) + 0.0054)
                        * (m_CO_core - delta_M_CO_shift - 34.8)**3
                        - 0.0013 * (m_CO_core - delta_M_CO_shift - 34.8)**2
                        + delta_M_PPI_extra_ML
                    )
                    if self.verbose:
                        print(f"delta_M_PPI: {delta_M_PPI} Msun")
                else:
                    delta_M_PPI = 0.0

                if delta_M_PPI <= 0.0:
                    # no PPI -> use CCSN prescription
                    if self.conserve_hydrogen_envelope:
                        m_PISN = m_star
                    else:
                        m_PISN = m_He_core
                else:
                    # PPI occurs
                    if self.conserve_hydrogen_PPI:
                        m_PISN = m_star - delta_M_PPI
                    else:
                        m_PISN = m_He_core - delta_M_PPI

                    if m_PISN < 0.0:
                        m_PISN = np.nan
                    else:
                        PISN_star = copy.deepcopy(star)
                        PISN_star.mass = m_PISN
                        if PISN_star.he_core_mass > m_PISN:
                            PISN_star.he_core_mass = m_PISN
                        if PISN_star.co_core_mass > m_PISN:
                            PISN_star.co_core_mass = m_PISN
                        m_rembar, _, _ = self.compute_m_rembar(PISN_star, m_PISN)

                        if m_rembar < 10:
                            m_PISN = np.nan
                        else:
                            m_PISN = m_rembar

            elif is_number(self.PISN) and m_He_core > self.PISN:
                m_PISN = self.PISN

            elif is_number(self.PISN) and 0.0 < m_He_core <= self.PISN:
                m_PISN = None

            else:
                raise ValueError("This choice {} of PISN is not available!".format(self.PISN))

        if self.verbose:
            if m_PISN is None:
                print("")
                print("The star did NOT lose any mass because of "
                      "PPIN or PISN.")
            elif not pd.isna(m_PISN):
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
                            Pwarn('Invalid co/He core masses! '
                                          'Setting m_WD=m_star!', "ApproximationWarning")
                        else:
                            Pwarn('co/He core masses are zero! '
                                          'Setting m_WD=m_star!', "ApproximationWarning")
                    else:
                        raise ModelError('Invalid co/He core masses! Cannot complete SN.')
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

        elif self.ECSN == 'Podsiadlowski+04':
            # Limits on He core mass progenitors of ECSN, default on cosmic
            min_M_He_ECSN = 1.4  # Msun from Podsiadlowski+2004
            max_M_He_ECSN = 2.5  # Msun from Podsiadlowski+2004

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
                            Pwarn('Invalid co/He core masses! '
                                          'Setting m_WD=m_star!', "ApproximationWarning")
                        else:
                            Pwarn('co/He core masses are zero! '
                                          'Setting m_WD=m_star!', "ApproximationWarning")
                    else:
                        raise ModelError('Invalid co/He core masses! Cannot complete SN.')
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
                            Pwarn('Invalid co/He core masses! '
                                          'Setting m_WD=m_star!', "ApproximationWarning")
                        else:
                            Pwarn('co/He core masses are zero! '
                                          'Setting m_WD=m_star!', "ApproximationWarning")
                    else:
                        raise ModelError('Invalid co/He core masses! Cannot complete SN.')
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
                "There is no information in the evolutionary history"
                "about STAR_STATES_CC."
            )
        if m_core is None or pd.isna(m_core):
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
                if self.ECSN == 'Podsiadlowski+04':
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
                if self.ECSN == 'Podsiadlowski+04':
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
                if self.ECSN == 'Podsiadlowski+04':
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
                if self.ECSN == 'Podsiadlowski+04':
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
                if self.ECSN == 'Podsiadlowski+04':
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

        elif self.mechanism == self.Maltsev25_engines:
            if star.SN_type == "ECSN":
                if self.ECSN == 'Podsiadlowski+04':
                    m_proto = 1.38
                else:
                    m_proto = m_core
                f_fb = 0.0
                m_fb = 0.0
                m_rembar = m_proto + m_fb
                state = 'NS'
            else:
                m_rembar, f_fb, state = self.Maltsev25_corecollapse(star,
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


    def get_combined_tilt(self, tilt_1, tilt_2, true_anomaly_1, true_anomaly_2):
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
                raise ModelError('The history did not contain core masses at'
                                 f' He depletion! {CO_core_mass}'
                                 f' {C_core_abundance}')

        return CO_core_mass, C_core_abundance

    def get_M4_mu4_Patton20(self, CO_core_mass, C_core_abundance):
        """Get the M4 and mu4 using Patton+20."""

        M4 = self.M4_interpolator.predict([[C_core_abundance, CO_core_mass]])
        mu4 = self.mu4_interpolator.predict([[C_core_abundance, CO_core_mass]])
        Xi = self.Xi_interpolator.predict([[C_core_abundance, CO_core_mass]])
        sc = self.sc_interpolator.predict([[C_core_abundance, CO_core_mass]])

        return M4, mu4, Xi, sc

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
            M4, mu4, Xi, sc = self.get_M4_mu4_Patton20(CO_core_mass, C_core_abundance)
            M4 = M4[0]
            mu4 = mu4[0]
            star.M4 = M4
            star.mu4 = mu4

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

    def Maltsev25_corecollapse(self, star, engine, conserve_hydrogen_envelope=False):
        """Compute supernova final remnant mass and fallback fraction.

        It uses the results from [8]_. The prediction for the core-collapse
        outcome is performed using the C core mass and its C abundance.
        The criterion by [8]_ is used to determine the final outcome.

        Parameters
        ----------
            star : obj
                Star object of a collapsing star containing the MESA profile.
            engine : str
                Engine to use for the core-collapse prescription
                Possible options are: 'M16'
            conserve_hydrogen_envelope : bool
                Whether to assume that the hydrogen envelope is conserved in direct collapse to a BH.

        Returns
        -------
        m_rem : double
            Remnant mass of the compact object in M_sun.
        f_fb : double
            Fallback mass of the compact object in M_sun.
        state : str
            'NS' if the remnant is a neutron star, 'BH' if the remnant is a black hole

        References
        ----------
        .. [8] K. Maltsev, F.R.N. Schneider, I. Mandel, B. Mueller, A. Heger,
            F.K. Roepke, E. Laplace, 2025, A&A, 700, A20. Explodability
            criteria for the neutrino-driven supernova mechanism

        """
        Muller_k_parameters = {
            'M16': [0.005, 0.420] # Section 3.1.1. of [8]_
        }

        if engine not in Muller_k_parameters.keys():
            raise ValueError("Engine " + engine + " is not avaiable for the "
                             "Maltsev+25 core-collapse prescription, "
                             "please choose one of the following engines to "
                             "compute the collapse: \n" + "\n".join(
                                list(Muller_k_parameters.keys())))
        else:

            CO_core_mass, C_core_abundance = self.get_CO_core_params(
                star, self.approx_at_he_depletion)
            M4, mu4, Xi, sc = self.get_M4_mu4_Patton20(CO_core_mass, C_core_abundance)
            M4 = M4[0]
            mu4 = mu4[0]
            Xi = Xi[0]
            sc = sc[0]
            mu4M4 = mu4*M4
            star.M4 = M4
            star.mu4 = mu4
            star.Xi = Xi
            star.sc = sc


            k1 = Muller_k_parameters[engine][0]
            k2 = Muller_k_parameters[engine][1]

            if CO_core_mass <= 2.5:
                m_rem = 1.25
                f_fb = 0.0
                state = 'NS'

            # In the Maltsev prescription, stars with CO core masses above 10 are allowed to explode.
            # However, since this outcome depends on the mass-transfer (MT) history, we handle it
            # in post-processing (for now). For all CO core masses above 10, we assume a failed supernova
            # with fallback = 1 at this stage.
            elif CO_core_mass >= 10.0:
                # Assuming BH formation by direct collapse
                if conserve_hydrogen_envelope:
                    m_rem = star.mass
                else:
                    m_rem = star.he_core_mass
                f_fb = 1.0
                state = 'BH'

            elif (CO_core_mass > 2.5) and (CO_core_mass < 10.0):
                successful_SN = self.explod_crit(Xi, sc, mu4M4, mu4, k1, k2)

                if successful_SN:
                    rem = self.NS_vs_fallbackBH(Xi, CO_core_mass, M4, mu4M4)
                    if rem == 'NS':  # successful SN with NS
                        m_rem = M4
                        f_fb = 0.0
                        state = 'NS'

                    else:  # successful SN but with fallback BH
                        if conserve_hydrogen_envelope:
                            m_rem = star.mass
                        else:
                            m_rem = star.he_core_mass

                        f_fb = 0.99
                        state = 'BH'

                else:
                    if conserve_hydrogen_envelope:
                        m_rem = star.mass
                    else:
                        m_rem = star.he_core_mass

                    f_fb = 1.0
                    state = 'BH'

        return m_rem, f_fb, state

    def NS_vs_fallbackBH(self, comp_val, mco_val, M4_val, mu4M4_val):
        a, b = 1.75, -0.044  # eq. (8) of [8]_
        # conditions for guaranteed NS formation (eq. 7)
        if comp_val <= 0.04 or (comp_val < a*mu4M4_val + b and comp_val <= 0.4) or M4_val/mco_val > 0.6:
            rem = 'NS'
        else:
            # stochastic determination of the remnant type (NS versus fallback-BH)
            rand_number = self.RNG.uniform(0,1)
            if rand_number <= 0.15:  # probability for fallback = 0.15 in Section 3.1.2.
                rem = 'fallback_BH'
            else:
                rem = 'NS'
        return rem

    # implemented from Maltsev+25
    def explod_crit(self, comp_val, sc_val, mu4M4_val, mu4_val, k1, k2):
        ff1, ff2 = [], []
        unclassified = True
        comp_crit1, comp_crit2 = 0.314, 0.544 # compactness
        sc_crit1, sc_crit2 = 0.988, 1.169 # central specific entropy
        mu4M4_crit1, mu4M4_crit2 = 0.247, 0.421 # product of M4 and mu4

        # check whether criterion for failed SN is fulfilled
        if comp_val > comp_crit2 or sc_val > sc_crit2:
            ff2.append(0)
            ff = False
            unclassified = False

        # check whether criterion for successful SN is fulfilled
        if comp_val < comp_crit1 or sc_val < sc_crit1:
            ff1.append(1)
            ff = True
            unclassified = False

        # if there is contradiction or if the progenitor is unclassified based on comp & s_c
        if (len(ff1) > 0 and len(ff2) > 0) or unclassified:

            # final fate classification based on mu4M4
            if mu4M4_val > mu4M4_crit2:
                ff = False
            elif mu4M4_val < mu4M4_crit1:
                ff = True
            # final fate classification based on reversed Ertl criterion
            elif k1 + k2*mu4M4_val - mu4_val > 0:
                ff = False
            else:
                ff = True
        return ff




class Sukhbold16_corecollapse(object):
    """Compute supernova final remnant mass, fallback fraction and CO type.

    This considers the He core mass of the nearest neighbor of the star
    prior to the collapse. Using a set of data for the He core
    mass of the compact object progenitors prior the collapse, the final
    remnant mass and stellar state of the compact object are known.

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
                data_download(set_name='auxiliary')

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
            raise ValueError("There is no information in the evolutionary "
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
            raise ValueError("Need a NS or BH to apply `Sukhbold16_corecollapse`.")

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
            filename = os.path.join(path_to_Couch_datasets, 'explDatsSTIR2.json')
            if not os.path.exists(filename):
                data_download(set_name='auxiliary')

            Couch_data_file = open(filename)
            Couch_data = json.load(Couch_data_file)
            Couch_data_file.close()
            Couch_data = Couch_data[turbulence_strength]

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

            # we need Sukhbold data for their cores
            Sukhbold_data = read_csv(path_to_Couch_datasets + "Sukhbold_Mzams_He_c_core.csv",
                                    usecols=[0, 1])

            MZAMS = Sukhbold_data["Mzams"]
            He_core_mass = Sukhbold_data["He_c_mass"]
            self.MZAMS_He_core_mass_Sukhbold_interpolator = interp1d(MZAMS, He_core_mass)

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
                np.array(Couch_He_c_mass).reshape((len(Couch_He_c_mass), 1)),
                np.array(Couch_state),
            )

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
            raise ValueError("There is no information in the evolutionary "
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
            raise ValueError("Need a NS or BH to apply `Sukhbold16_corecollapse`.")

        return float(m_rem), f_fb, state



def check_SN_CO_match(SN_type, state):
    '''Check if the SN type matches the stellar state of the given star.

    Parameters
    ----------
    SN_type : str
        SN type of the star.
    state : str
        Stellar state of the star.

    Returns
    -------
    correct_SN_type : bool
        True if the SN type matches the stellar state of the star.
    '''
    # TODO: remove star.state == PISN, because PISN shouldn't be a stellar state
    if state == 'PISN':
        state = 'massless_remnant'
    correct_SN_type = True
    if state == 'WD' and SN_type != "WD":
        correct_SN_type = False
    elif (state == "NS") and \
            (SN_type != 'ECSN' and
             SN_type != "CCSN"):
        correct_SN_type = False
    elif (state =="BH") and \
            (SN_type != "CCSN" and
            SN_type != 'PPISN'):
        correct_SN_type = False
    elif (state == "massless_remnant" and SN_type != 'PISN'):
        correct_SN_type = False
    return correct_SN_type
