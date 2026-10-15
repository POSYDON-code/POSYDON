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


import numpy as np
from posydon.grids.SN_MODELS import DEFAULT_SN_MODEL
from posydon.utils.common_functions import (
    set_binary_to_failed,
)
from posydon.utils.posydonerror import ModelError

####################################################################################################
from abc import ABC, abstractmethod

from posydon.binary_evol.SN.explodability.CCSN_explosion_criteria import (
    get_CCSN_prescription,
)
from posydon.binary_evol.SN.explodability.ECSN_prescriptions import (
    get_ECSN_prescription,
)
from posydon.binary_evol.SN.explodability.PISN_prescriptions import (
    get_PISN_prescription,
)
from posydon.binary_evol.SN.kicks.kick_prescriptions import get_kick_prescription
from posydon.binary_evol.SN.remnant_mass.remnant_mass_prescriptions import (
    get_remnant_mass_prescription,
)

####################################################################################################

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
        supernova event, else no kicks are
 taking into account
        for any supernova outcome.

    kick_normalisation : str
        Renormalise the kick by:
        'one_minus_fallback' : (1-f_fb)
        'one_over_mass' : 1.4/m_BH
        'zero' : 0.
        'one' : 1.
        'NS_one_minus_fallback_BH_one': 1 for BH, (1-f_fb) for NS

    no_kick_direct_collapse : bool
        If True, ignore the kick for BHs formed by direct collapse
        (fallback fraction f_fb = 1). Only affects the CCSN/PPISN/PISN
        branch of the orbital kick, where the remnant state is 'BH'.

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
        "kick": True,
        "kick_normalisation": 'one_over_mass',
        "kick_prescription": 'maxwellian',
        "sigma_kick_CCSN_NS": 265.0,
        "mean_kick_CCSN_NS": None,
        "sigma_kick_CCSN_BH": 265.0,
        "mean_kick_CCSN_BH": None,
        "sigma_kick_ECSN": 20.0,
        "mean_kick_ECSN": None,
        "no_kick_direct_collapse": False,
        # other
        "RNG": None,
        "verbose": False
    }

    DEFAULT_KWARGS.update(DEFAULT_SN_MODEL)

    def __init__(self,**kwargs):
        unknown_kwargs = set(kwargs) - set(self.DEFAULT_KWARGS)
        if unknown_kwargs:
            raise ValueError(f"Invalid arguments passed to step_NS: {unknown_kwargs}")
        config = self.DEFAULT_KWARGS.copy()
        config.update(kwargs)
        if config["RNG"] is None:
            config["RNG"] = np.random.default_rng()
        for key, value in config.items():
            setattr(self, key, value)

        self._explosion_checksget_ECSN_prescription = (
            get_ECSN_prescription(**config),     # -> "WD" | "ECSN" | None
            get_CCSN_prescription(**config),     # -> "CCSN"
        )
        self._get_remnant_mass = get_remnant_mass_prescription(**kwargs)
        self._get_kick = get_kick_prescription(**kwargs)

    def __call__(self,binary):
        dying_star,other_star = self._select_dying_star(binary)
        model_err = self.collapse_star(dying_star)
        if model_err is not None:
            set_binary_to_failed(binary)
            raise ModelError(model_err)
        self._reset_other_star_properties(star=other_star)
        binary.update_star_states()
        self.orbital_kick(binary)
        self.finalise_step(binary)

    @staticmethod
    def _select_dying_star(binary):
        if binary.event=="CC1":
            return binary.star_1, binary.star_2
        elif binary.event=="CC2":
            return binary.star_2, binary.star_1
        else:
            raise ValueError(
                "Invalid call of step_SN:"
                "Called during binary event {event}".format(event=binary.event)
                "Step_SN should only be called if binary event is CC1 or CC2"
            )

    def collapse_star(self,star):
        outcome = self.collapse_outcome(star)
        m_rem, f_fb = self._remnant_mass(star, outcome)
        self._update_params(star, outcome, m_rem, f_fb)

    def collapse_outcome(self,star):
        return next(outcome for outcome in (check(star) for check in self._explosion_checks)
                    if outcome is not None)

    def _reset_other_star_properties(self, star):
        """Reset the properties of the star that is not being collapsed."""
        star.lg_mdot = None
        star.lg_system_mdot = None

    def finalise_step(self,binary):
        return
