"""Common envelope evolution step.

Calculates the evolution of a star in a BinaryStar object. It uses
the alpha-perscription and caclulates how much the orbit has to shrink in
order to expel its envelope. If in that separation, one of the star fills its
RL, then the system is considered a merger. Else the common envelope is lost
and we are left with a binary system with the core of the donor star (which
initiated the unstable CEE) and the core of the copmanion. For now it works
with the donor star being either a type of H_giant (with a He core) or
He_giant (with a CO core) and the companion star being either a MS or compact
object star, not contribyting to the CE.

If profiles exist we are able to calculate on the spot the lambda parameter
for the donor star, else we use the default values

Parameters
----------
binary : BinaryStar (An object of class BinaryStar defined in POSYDON)
verbose : Boolean
    In case we want information about the CEE  (the default is False).

"""


__authors__ = [
    "Emmanouil Zapartas <ezapartas@gmail.com>",
    "Devina Misra <devina.misra@unige.ch>",
    "Jaime Roman Garza <Jaime.Roman@etu.unige.ch>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
    "Tassos Fragos <Anastasios.Fragkos@unige.ch>",
]


import numpy as np
from posydon.utils import common_functions as cf
from posydon.utils import constants as const
import warnings

from posydon.binary_evol.binarystar import BINARYPROPERTIES
from posydon.binary_evol.singlestar import STARPROPERTIES
from posydon.utils.common_functions import PATH_TO_POSYDON
from posydon.utils.common_functions import check_state_of_star
from posydon.utils.common_functions import calculate_lambda_from_profile


warnings.simplefilter('always', UserWarning)


MODEL = {"prescription": 'alpha-lambda',
         "common_envelope_efficiency": 1.0,
         "common_envelope_lambda_default": 0.5,
         "common_envelope_option_for_lambda": 'lambda_from_grid_final_values',
         "common_envelope_option_for_HG_star": "optimistic",
         "common_envelope_alpha_thermal": 1.0,
         "core_definition_H_fraction": 0.1,     # with 0.01 no CE BBHs
         "core_definition_He_fraction": 0.1,
         "CEE_tolerance_err": 0.001,
         "verbose": False,
         "common_envelope_option_after_succ_CEE": 'core_not_replaced_noMT',
         # "core_replaced_noMT" for core_definition_H_fraction=0.01
         }


# common_envelope_option_for_lambda:
# 1) 'default_lambda': using for lambda the constant value of
# common_envelope_lambda_default parameter
# 2) 'lambda_from_grid_final_values': using lambda parameter from MESA history
# which was calulated ni the same way as method (5) below
# 3) 'lambda_from_profile_gravitational': calculating the lambda parameter
# from the donor's profile by using the gravitational binding energy from the
# surface to the core (needing "mass", and "radius" as columns in the profile)
# 4) 'lambda_from_profile_gravitational_plus_internal': as above but taking
# into account a factor of common_envelope_alpha_thermal * internal energy too
# in the binding energy (needing also "energy" as column in the profile)
# 5) 'lambda_from_profile_gravitational_plus_internal_minus_recombination':
# as above but not taking into account the recombination energy in the internal
# energy (needing also "y_mass_fraction_He", "x_mass_fraction_H",
# "neutral_fraction_H", "neutral_fraction_He", and "avg_charge_He" as column
# in the profile)
# the mass fraction of an element which is used as threshold to define a core,


STAR_STATE_POST_MS = [
    "H-rich_Core_H_burning",
    "H-rich_Shell_H_burning",
    "H-rich_Core_He_burning",
    "H-rich_Central_He_depleted",
    "H-rich_Central_C_depletion",
    "H-rich_non_burning"
]


STAR_STATE_POST_HeMS = [
    'stripped_He_Core_He_burning',
    'stripped_He_Central_He_depleted',
    'stripped_He_Central_C_depletion',
    'stripped_He_non_burning'
]


class StepCEE(object):
    """Compute supernova final remnant mass, fallback fraction & stellar state.

    This consider the nearest neighboor of the He core mass of the star,
    previous to the collapse. Considering a set of data for which the He core
    mass of the compact object projenitos previous the collapse, the final
    remnant mass and final stellar state of the compact object is known.

    Parameters
    ----------
    verbose : bool
        If True, the messages will be prited in the console.

    Keyword Arguments
    -----------------
    prescription : str
        Prescription to use for computing the prediction of common enevelope
        evolution. Available options are:

        * 'alpha-lambda' : Considers the the alpha-lambda prescription
          described in [1]_ and [2]_ to predict the outcome of the common
          envelope evolution. If the profile of the donor star is available
          then it is used to compute the value of lambda.

    References
    ----------
    .. [1] Webbink, R. F. (1984). Double white dwarfs as progenitors of R
        Coronae Borealis stars and Type I supernovae. The Astrophysical
        Journal, 277, 355-360.

    .. [2] De Kool, M. (1990). Common envelope evolution and double cores of
        planetary nebulae. The Astrophysical Journal, 358, 189-195.
    """

    def __init__(
            self, prescription=MODEL['prescription'],
            common_envelope_efficiency=MODEL['common_envelope_efficiency'],
            common_envelope_lambda_default=MODEL[
                'common_envelope_lambda_default'],
            common_envelope_option_for_lambda=MODEL[
                'common_envelope_option_for_lambda'],
            common_envelope_option_for_HG_star=MODEL[
                'common_envelope_option_for_HG_star'],
            common_envelope_option_after_succ_CEE=MODEL[
                'common_envelope_option_after_succ_CEE'],
            common_envelope_alpha_thermal=MODEL[
                'common_envelope_alpha_thermal'],
            core_definition_H_fraction=MODEL[
                'core_definition_H_fraction'],
            core_definition_He_fraction=MODEL[
                'core_definition_He_fraction'],
            CEE_tolerance_err=MODEL['CEE_tolerance_err'],
            verbose=MODEL['verbose'],
            **kwargs):
        """Initialize a StepCEE instance."""
        # read kwargs to initialize the class
        if kwargs:
            for key in kwargs:
                if key not in MODEL:
                    raise ValueError(key + " is not a valid parameter name!")
            for varname in MODEL:
                default_value = MODEL[varname]
                setattr(self, varname, kwargs.get(varname, default_value))
        else:
            self.prescription = prescription
            self.common_envelope_efficiency = common_envelope_efficiency
            self.common_envelope_lambda_default = \
                common_envelope_lambda_default
            self.common_envelope_option_for_lambda = \
                common_envelope_option_for_lambda
            self.common_envelope_option_for_HG_star = \
                common_envelope_option_for_HG_star
            self.common_envelope_alpha_thermal = common_envelope_alpha_thermal
            self.core_definition_H_fraction = core_definition_H_fraction
            self.core_definition_He_fraction = core_definition_He_fraction
            self.CEE_tolerance_err = CEE_tolerance_err
            self.common_envelope_option_after_succ_CEE = \
                common_envelope_option_after_succ_CEE

        self.verbose = verbose
        self.path_to_posydon = PATH_TO_POSYDON

    def __call__(self, binary):
        """Perform the CEE step for a BinaryStar object."""
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

        # Check for double CE
        if binary.event in ["oDoubleCE1", "oDoubleCE2"]:
            double_CE = True
        else:
            double_CE = False

        if self.verbose:
            print("binary.event : ", binary.event)
            print("common_envelope_efficiency : ",
                  self.common_envelope_efficiency)
            print("common_envelope_option_for_lambda : ",
                  self.common_envelope_option_for_lambda)

        # Check to make sure binary can go through a CE
        if donor_star.state in ['H-rich_Core_H_burning',
                                'stripped_He_Core_He_burning']:
            # system merges
            binary.state = 'merged'
            if binary.event in ["oCE1", "oDoubleCE1"]:
                binary.event = "oMerging1"
            if binary.event in ["oCE2", "oDoubleCE2"]:
                binary.event = "oMerging2"

            return

        if self.common_envelope_option_for_HG_star == "pessimistic":
            # Merging if HG donor
            if donor_star.state in ['H-rich_Shell_H_burning']:
                # system merges
                binary.state = 'merged'
                if binary.event in ["oCE1", "oDoubleCE1"]:
                    binary.event = "oMerging1"
                if binary.event in ["oCE2", "oDoubleCE2"]:
                    binary.event = "oMerging2"

                return

        # Calculate binary's evolution
        if self.prescription == 'alpha-lambda':

            # Determine lambda and associated CE parameters
            lambda1_CE, mc1_i, rc1_i, donor_type = self.calculate_lambda_CE(
                donor_star, verbose=self.verbose)

            if double_CE:
                lambda2_CE, mc2_i, rc2_i, comp_type = self.calculate_lambda_CE(
                    comp_star, verbose=self.verbose)
            else:
                lambda2_CE = None
                mc2_i = comp_star.mass
                rc2_i = 10**comp_star.log_R
                comp_type = "not_giant_companion"

            # Evolve binary
            self.CEE_simple_alpha_prescription(
                binary, donor_star, comp_star, lambda1_CE,
                mc1_i, rc1_i, donor_type, lambda2_CE, mc2_i, rc2_i, comp_type,
                double_CE=double_CE, verbose=self.verbose,
                common_envelope_option_after_succ_CEE=(
                    self.common_envelope_option_after_succ_CEE),
                core_definition_H_fraction=self.core_definition_H_fraction,
                core_definition_He_fraction=self.core_definition_He_fraction)
        else:
            raise ValueError("Invalid common envelope prescription given.")

    def calculate_lambda_CE(self, donor, verbose=False):
        """Calculate lambda_CE from an assumed constant value or from profile.

        If lambda_CE is calculated from donor's profile, we also pass a more
        accurate calculation of the donor core mass for the purposes of CEE.


        Parameters
        ----------
        donor : instance of SingleStar
            The donor star to calculate pre-CE quantities.
        verbose : bool
            In case we want information about the CEE  (the default is False).

        Returns
        -------
        lambda_CE: float
            lambda_CE calculated either as an assumed constant value or
            calculated from profile
        mc1_i: float
            If lambda_CE is calculated from donor's profile, we also pass a
            more accurate calculation of the donor core mass for the purposes
            of CEE.
        rc1_i: float
            The radius of the donor's core

        """
        m1_i = donor.mass
        radius1 = 10. ** donor.log_R

        if donor.state in STAR_STATE_POST_MS:       # "H_Giant":
            donor_type = 'He_core'
            core_element_fraction_definition = self.core_definition_H_fraction
            if core_element_fraction_definition not in [0.01, 0.1, 0.3]:
                raise ValueError("He-core defintion should always be "
                                 "set to a H abundance of 1%, 10%, or 30%")
        elif donor.state in STAR_STATE_POST_HeMS:   # "He_Giant":
            donor_type = 'CO_core'
            core_element_fraction_definition = self.core_definition_He_fraction
            if core_element_fraction_definition != 0.1:
                raise ValueError("CO-core should always be "
                                 "set to a He abudance of 10%")
        else:
            raise ValueError("type = %s of donor of CEE not recognized"
                             % donor.state)

        if self.common_envelope_option_for_lambda == "default_lambda":
            lambda_CE = self.common_envelope_lambda_default
            if donor_type == 'He_core':
                mc1_i = donor.he_core_mass
                rc1_i = donor.he_core_radius
            elif donor_type == "CO_core":
                mc1_i = donor.co_core_mass
                rc1_i = donor.co_core_radius
            else:
                raise ValueError("type = %s of donor of CEE not recognized"
                                 % donor.state)
        elif (self.common_envelope_option_for_lambda
                == "lambda_from_grid_final_values"):
            if donor_type == 'CO_core':
                if core_element_fraction_definition == 0.1:
                    lambda_CE = donor.lambda_CE_pure_He_star_10cent
                    mc1_i = donor.m_core_CE_pure_He_star_10cent
                    rc1_i = donor.r_core_CE_pure_He_star_10cent
                else:
                    raise ValueError("An impossible core_definition was "
                                     "provided in calculate lambda function")
            elif donor_type == 'He_core':
                if core_element_fraction_definition == 0.01:
                    lambda_CE = donor.lambda_CE_1cent
                    mc1_i = donor.m_core_CE_1cent
                    rc1_i = donor.r_core_CE_1cent
                elif core_element_fraction_definition == 0.1:
                    lambda_CE = donor.lambda_CE_10cent
                    mc1_i = donor.m_core_CE_10cent
                    rc1_i = donor.r_core_CE_10cent
                elif core_element_fraction_definition == 0.3:
                    lambda_CE = donor.lambda_CE_30cent
                    mc1_i = donor.m_core_CE_30cent
                    rc1_i = donor.r_core_CE_30cent
                else:
                    raise ValueError("An impossible core_definition was "
                                     "provided in calculate lambda function")
            else:
                raise ValueError("An impossible CE donor-star type was "
                                 "provided in calculate lambda function")

            if verbose:
                print("verbose for lambda_from_grid_final_values")
                print("core_element_fraction_definition :",
                      core_element_fraction_definition)
                print("lambda_CE : ", lambda_CE)
                print("donor state: ", donor.state)
                print("donor_type: ", donor_type)
                print("donor R, rc1_i (for CEE calculation), he_core_radius, "
                      "co_core_radius  = ", radius1, rc1_i,
                      donor.he_core_radius, donor.co_core_radius)
                print("donor M, mc1_i (for CEE calculation), he_core_mass, "
                      "co_core_mass = ", m1_i, mc1_i,
                      donor.he_core_mass, donor.co_core_mass)

        elif donor.profile is None:
            warnings.warn("Profile does not exist -- Proceeding with "
                          "default_lambda alpha-CE prescription")
            # like in the "default_lambda" option
            lambda_CE = self.common_envelope_lambda_default
            if donor_type == 'He_core':
                mc1_i = donor.he_core_mass
                rc1_i = donor.he_core_radius
            elif donor_type == "CO_core":
                mc1_i = donor.co_core_mass
                rc1_i = donor.co_core_radius
            else:
                raise ValueError("type = %s of donor of CEE not recognized"
                                 % donor.state)
        else:
            valid_options = ["lambda_from_profile_gravitational",
                             "lambda_from_profile_gravitational_plus_internal",
                             "lambda_from_profile_gravitational_plus_"
                             + "internal_minus_recombination"]
            if self.common_envelope_option_for_lambda not in valid_options:
                raise ValueError("common_envelope_option_for_lambda not a "
                                 "valid option - Should never occur")

            '''
            lambda_CE, mc1_i, rc1_i = \
                self.calculate_lambda_from_profile_duringstepCEE(
                    donor, core_element_fraction_definition, verbose)
            '''
            lambda_CE, mc1_i, rc1_i = calculate_lambda_from_profile(
                profile=donor.profile, donor_star_state=donor.state,
                m1_i=m1_i, radius1=radius1,
                common_envelope_option_for_lambda=(
                    self.common_envelope_option_for_lambda),
                core_element_fraction_definition=(
                    core_element_fraction_definition),
                ind_core=None,
                common_envelope_alpha_thermal=(
                    self.common_envelope_alpha_thermal),
                tolerance=self.CEE_tolerance_err,
                verbose=verbose)

        return lambda_CE, mc1_i, rc1_i, donor_type

    def CEE_simple_alpha_prescription(
            self, binary, donor, comp_star, lambda1_CE, mc1_i, rc1_i,
            donor_type, lambda2_CE, mc2_i, rc2_i, comp_type, double_CE=False,
            verbose=False, common_envelope_option_after_succ_CEE=MODEL[
                'common_envelope_option_after_succ_CEE'],
            core_definition_H_fraction=MODEL['core_definition_H_fraction'],
            core_definition_He_fraction=MODEL['core_definition_He_fraction']):
        """Apply the alpha-lambda common-envelope prescription.

        It uses energetics to calculate the shrinakge of the orbit
        (and possible merger) of the two components, as well as upadate
        their new masses, during a common-envelope phase.

        Parameters
        ----------
        binary : BinaryStar object
            Instance of the binary to evolve.
        donor : SingleStar object
            The donor star
        comp_star : SingleStar object
            The companion star
        lambda1_CE : float
            Lambda previously calculated for the donor's envelope
            binding energy.
        mc1_i : float
            core mass of the donor
        rc1_i : float
            core radius of the donor
        donor_type : string
            String dictating whether the donor has a H-envelope or He-envelope
        lambda2_CE : float
            Lambda previously calculated for the companion's envelope
            binding energy.
        mc2_i : float
            core mass of the companion
        rc2_i : float
            core radius of the companion
        comp_type : string
            String dictating whether the companion has a
            H-envelope or He-envelope.
        double_CE : bool
            In case we have a double CE situation.
        verbose : bool
            In case we want information about the CEE.
        common_envelope_option_after_succ_CEE: str
            Options are:
            
            1) "core_replaced_noMT"
                he_core_mass/radius (or co_core_mass/radius for CEE of
                stripped_He*) are replaced according to the new core boundary
                used for CEE (based on core_definition_H/He_fraction) but no
                other change in period after succesful ejection at
                alpha-lambda prescription.
            2) "core_not_replaced_noMT"
                he_core_mass/radius (or co_core_mass/radius for CEE of
                stripped_He*) staying as preCEE and no other change in period
                after succesful ejection at alpha-lambda prescription.
            3) "core_not_replaced_stableMT"
                he_core_mass/radius (or co_core_mass/radius for CEE of
                stripped_He*) staying as preCEE and after succesful ejection at
                alpha-lambda prescription, we assume an instantaneous stableMT
                phase (non-conservative, with mass lost from accretor) from the
                donor (or from donor and simultaneously the accretor at
                double_CE), taking away the extra "core" mass as defined by the
                core boundary used for CEE (based on
                core_definition_H/He_fraction).
            4) "core_not_replaced_windloss"
                he_core_mass/radius (or co_core_mass/radius for CEE of
                stripped_He*) staying as preCEE and after succesful ejection
                at alpha-lambda prescription, we assume a instantaneous
                windloss phase from the donor (or from donor and simultaneously
                the accretor at double_CE), taking away the extra "core" mass
                as defined by the core boundary used for CEE (based on
                core_definition_H/He_fraction).
        core_definition_H_fraction: float
            The value of the H abundance to define the envelope-core boundary
            in He_cores.
        core_definition_He_fraction: float
            The value of the He abundance to define the envelope-core boundary
            in CO_cores.

        """
        # Get star properties
        m1_i = donor.mass
        radius1 = 10**donor.log_R
        state1_i = donor.state
        m2_i = comp_star.mass
        radius2 = 10**comp_star.log_R

        if verbose:
            print("Pre CE")
            print("binary event :", binary.event)
            print("double CEE : ", double_CE)
            print("donor state: ", state1_i)
            print("donor mass / core mass = ", m1_i, mc1_i)
            print("companion state: ", comp_star.state)
            print("companion mass / core mass = ", m2_i, mc2_i)

        # Get binary parameters
        period_i = binary.orbital_period
        alpha_CE = self.common_envelope_efficiency

        # calculate evolution of the orbit
        ebind_i = (-const.standard_cgrav / lambda1_CE
                   * (m1_i * const.Msun * (m1_i - mc1_i) * const.Msun)
                   / (radius1 * const.Rsun))
        if double_CE:
            ebind_i += (-const.standard_cgrav / lambda2_CE
                        * (m2_i * const.Msun * (m2_i - mc2_i) * const.Msun)
                        / (radius2 * const.Rsun))

        separation_i = const.Rsun * cf.orbital_separation_from_period(
            period_i, m1_i, m2_i)   # in cgs units
        eorb_i = (-0.5 * const.standard_cgrav * m1_i * const.Msun
                  * m2_i * const.Msun / separation_i)
        # TODO use core masses or total masses?

        eorb_postCEE = eorb_i + ebind_i/alpha_CE

        separation_postCEE = (-0.5 * const.standard_cgrav * mc1_i * const.Msun
                              * mc2_i * const.Msun / eorb_postCEE)

        # Check to make sure final orbital separation is positive
        if not (separation_postCEE > -self.CEE_tolerance_err):
            raise Exception("CEE problem, negative postCEE separation")

        if verbose:
            print("CEE alpha-lambda prescription")
            print("ebind_i", ebind_i)
            print("eorb_i", eorb_i)
            print("eorb_postCEE", eorb_postCEE)
            print("DEorb", eorb_postCEE - eorb_i)
            print("separation_i in Rsun", separation_i/const.Rsun)
            print("separation_postCEE in Rsun", separation_postCEE/const.Rsun)

        # now we check if the roche Lobe of any of the cores that spiralled-in
        # will be filled if reached this final separation
        RL1 = cf.roche_lobe_radius(mc1_i/mc2_i, separation_postCEE/const.Rsun)
        RL2 = cf.roche_lobe_radius(mc2_i/mc1_i, separation_postCEE/const.Rsun)

        if verbose:
            print("donor radius / core radius / RL1:", radius1, rc1_i, RL1)
            print("companion radius / core radius / RL2:", radius2, rc2_i, RL2)

        # Check if binary merges or survives
        if ((rc1_i - RL1) < self.CEE_tolerance_err
                and (rc2_i - RL2) < self.CEE_tolerance_err):
            # system survives CEE, the whole CE is ejected
            # and the new orbital separation for the cores is returned
            if verbose:
                print("system survives CEE, the whole CE is ejected and the "
                      "new orbital separation for the cores is returned")

            # possible change of core masses
            if common_envelope_option_after_succ_CEE in ["core_replaced_noMT"]:
                mc1_f = mc1_i
                rc1_f = rc1_i
                mc2_f = mc2_i
                rc2_f = rc2_i
                if verbose:
                    print("m1 core mass pre CEE (He,CO) / to after CEE : ",
                          donor.he_core_mass, donor.co_core_mass, mc1_f)
                    print("m1 core radius pre CEE (He,CO) / to after CEE : ",
                          donor.he_core_radius, donor.co_core_radius, rc1_f)
                    print("m2 core mass pre CEE (He,CO) / to after CEE : ",
                          comp_star.he_core_mass, comp_star.co_core_mass, mc2_f
                          )
                    print("m2 core radius pre CEE (He,CO) / to after CEE : ",
                          comp_star.he_core_radius,
                          comp_star.co_core_radius, rc2_f)
            elif common_envelope_option_after_succ_CEE in [
                    "core_not_replaced_noMT",
                    "core_not_replaced_stableMT",
                    "core_not_replaced_windloss"]:
                if donor_type == 'He_core':
                    mc1_f = donor.he_core_mass
                    rc1_f = donor.he_core_radius
                elif donor_type == "CO_core":
                    mc1_f = donor.co_core_mass
                    rc1_f = donor.co_core_radius
                if mc1_f > mc1_i:
                    mc1_f = mc1_i
                    warnings.warn(
                        "donor core mass final (after even stable, postCEE MT)"
                        " assumed higher than postCEE core mass. Equalized to "
                        "postCEE mass")
                if not double_CE:
                    mc2_f = mc2_i
                    rc2_f = rc2_i
                else:
                    if comp_type == 'He_core':
                        mc2_f = comp_star.he_core_mass
                        rc2_f = donor.he_core_radius
                    elif comp_type == "CO_core":
                        mc2_f = comp_star.co_core_mass
                        rc2_f = donor.co_core_radius
                    elif comp_type == "not_giant_companion":
                        mc2_f = comp_star.mass
                        rc2_f = 10.**(comp_star.log_R)
                    if mc2_f > mc2_i:
                        mc2_f = mc2_i
                        warnings.warn("accretor's core mass final (after even "
                                      "non-conservative stable, postCEE MT) "
                                      "assumed higher that postCEE core mass. "
                                      "Equialized to postCEE mass")
                if verbose:
                    print("difference between m1 core mass defined by CEE step"
                          " / to the final one as pre CEE : ", mc1_f, mc1_i)
                    print("difference between r1 core mass defined by CEE step"
                          " / to the final one as pre CEE : ", rc1_f, rc1_i)
                    print("difference between m2 core mass defined by CEE step"
                          " / to the final one as pre CEE : ", mc2_f, mc2_i)
                    print("difference between r2 core mass defined by CEE step"
                          " / to the final one as pre CEE : ", rc2_f, rc2_i)
            else:
                raise ValueError(
                    "Not accepted option in common_envelope_option_after_succ_"
                    "CEE = {}, dont know how to proceed".
                    format(common_envelope_option_after_succ_CEE))

            # possible change of orbital period after CEE
            orbital_period_postCEE = cf.orbital_period_from_separation(
                separation_postCEE / const.Rsun, mc1_i, mc2_i)
            if common_envelope_option_after_succ_CEE in [
                    "core_replaced_noMT", "core_not_replaced_noMT"]:
                orbital_period_f = orbital_period_postCEE
            elif common_envelope_option_after_succ_CEE in [
                    "core_not_replaced_stableMT"]:
                # An assumed stable mass transfer case after postCEE with
                # fully non-conservative MT and mass lost from the vicinity
                # of the accretor:
                orbital_period_f = cf.period_change_stabe_MT(
                    orbital_period_postCEE, Mdon_i=mc1_i, Mdon_f=mc1_f,
                    Macc_i=mc2_i, alpha=0.0, beta=1.0)
                if double_CE:
                    # a reverse stable MT assumed to be happening
                    # at the same time
                    orbital_period_f = cf.period_change_stabe_MT(
                        orbital_period_f, Mdon_i=mc2_i, Mdon_f=mc2_f,
                        Macc_i=mc1_i, alpha=0.0, beta=1.0)
                if verbose:
                    print("during the assumed stable MT phase after postCE")
                    print("with 'common_envelope_option_after_succ_CEE' :",
                          common_envelope_option_after_succ_CEE)
                    print("the orbit cahnged from postCEE : ",
                          orbital_period_postCEE)
                    print("to : ", orbital_period_f)
            elif common_envelope_option_after_succ_CEE in [
                    "core_not_replaced_windloss"]:
                # An assumed wind loss after postCEE with mass lost
                # from the vicinity of the donor
                orbital_period_f = cf.period_change_stabe_MT(
                    orbital_period_postCEE, Mdon_i=mc1_i, Mdon_f=mc1_f,
                    Macc_i=mc2_i, alpha=1.0, beta=0.0)
                if double_CE:
                    # a wind mass loss from the 2nd star assumed to be
                    # happening at the same time
                    orbital_period_f = cf.period_change_stabe_MT(
                        orbital_period_f, Mdon_i=mc2_i, Mdon_f=mc2_f,
                        Macc_i=mc1_i, alpha=1.0, beta=0.0)
                if verbose:
                    print("during the assumed windloss phase after postCE")
                    print("with 'common_envelope_option_after_succ_CEE' :",
                          common_envelope_option_after_succ_CEE)
                    print("the orbit cahnged from postCEE : ",
                          orbital_period_postCEE)
                    print("to : ", orbital_period_f)
            separation_f = cf.orbital_separation_from_period(orbital_period_f,
                                                             mc1_f, mc2_f)

            if state1_i == 'stripped_He_Central_He_depleted':
                if donor == binary.star_1:
                    binary.event = 'CC1'
                elif donor == binary.star_2:
                    binary.event = 'CC2'
            else:
                binary.event = None

            # Adjust binary properties
            binary.separation = separation_f
            binary.orbital_period = orbital_period_f
            binary.eccentricity = 0.0
            binary.state = 'detached'
            # binary.event = None

            for key in BINARYPROPERTIES:
                # the binary attributes that are changed in the CE step
                if key not in ["separation", "orbital_period",
                               "eccentricity", "state", "event"]:
                    # the binary attributes that keep the same value from the
                    # previous step
                    if key not in ["time", "V_sys", "mass_transfer_case",
                                   "nearest_neighbour_distance"]:
                        setattr(binary, key, np.nan)  # the rest become np.nan

            stars = [donor]
            star_types = [donor_type]
            core_masses = [mc1_f]
            core_radii = [rc1_f]

            if verbose:
                print("CEE succesfully ejected")
                print("new orbital period = ", binary.orbital_period)
                print("binary event : ", binary.event)
                print("double CEE : ", double_CE)

            if double_CE:
                stars.append(comp_star)
                star_types.append(comp_type)
                core_masses.append(mc2_f)
                core_radii.append(rc2_f)

            # Adjust donor star properties
            for star, star_type, core_mass, core_radius in zip(
                    stars, star_types, core_masses, core_radii):
                star.mass = core_mass
                star.log_R = np.log10(core_radius)
                attributes_changing = [
                                        "mass",
                                        "log_R",
                                        "state"
                                        ]

                if star_type == 'He_core':
                    if common_envelope_option_after_succ_CEE in [
                            "core_replaced_noMT"]:
                        star.surface_h1 = core_definition_H_fraction
                    else:
                        star.surface_h1 = 0.01
                    star.he_core_mass = core_mass
                    star.he_core_radius = core_radius
                    if star.metallicity is None:
                        star.surface_he4 = 1 - star.surface_h1 - 0.0142
                    else:
                        star.surface_he4 = 1.0-star.surface_h1-star.metallicity
                    star.log_LH = -1e99
                    attributes_changing.extend([
                                            "surface_h1",
                                            "surface_he4",
                                            "log_LH",
                                            "log_LHe",
                                            'he_core_mass',
                                            'he_core_radius'
                        ])
                elif star_type == 'CO_core':
                    star.he_core_mass = core_mass
                    star.he_core_radius = core_radius
                    star.co_core_mass = core_mass
                    star.co_core_radius = core_radius
                    star.surface_h1 = 0.0
                    if common_envelope_option_after_succ_CEE in [
                            "core_replaced_noMT"]:
                        star.surface_he4 = core_definition_He_fraction
                    else:
                        star.surface_he4 = 0.1
                    star.log_LH = -1e99
                    star.log_LHe = -1e99
                    attributes_changing.extend([
                                            "surface_h1",
                                            "surface_he4",
                                            "log_LH",
                                            "log_LHe",
                                            'he_core_mass',
                                            'he_core_radius',
                                            'co_core_mass',
                                            'co_core_radius'

                        ])
                elif star_type == "not_giant_companion":
                    continue
                else:
                    raise ValueError("Unrecognized star type:", star_type)

                state_old = star.state
                star.state = check_state_of_star(star)
                if verbose:
                    print("star state of before/ after CE : ",
                          state_old, star.state)

                for key in STARPROPERTIES:
                    if key not in attributes_changing:
                        # the singlestar attributes that are changed
                        # in the CE step

                        # the singlestar attributes that keep the same value
                        # from previous step if not in attributes_changing
                        if key not in [
                                "surface_h1", "surface_he4", "log_LH",
                                "log_LHe", "he_core_mass", "he_core_radius",
                                "co_core_mass", "co_core_radius", "center_h1",
                                "center_he4", "center_c12", "center_o16",
                                "center_n14", "metallicity", "log_LZ",
                                "log_Lnuc", "c12_c12", "center_gamma",
                                "avg_c_in_c_core"]:
                            setattr(star, key, np.nan)  # the rest become NaN's

        else:
            # system merges
            binary.state = 'merged'
            if binary.event in ["oCE1", "oDoubleCE1"]:
                binary.event = "oMerging1"
            if binary.event in ["oCE2", "oDoubleCE2"]:
                binary.event = "oMerging2"

            if verbose:
                print("system merges due to one of the two star's core filling"
                      "its RL")
                print("Rdonor core vs RLdonor core = ", rc1_i, RL1)
                print("Rcompanion vs RLcompanion= ", rc2_i, RL2)

        return
