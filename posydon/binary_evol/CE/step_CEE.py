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
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
]


import numpy as np
import pandas as pd
from posydon.utils import common_functions as cf
from posydon.utils import constants as const
from posydon.binary_evol.binarystar import BINARYPROPERTIES
from posydon.binary_evol.singlestar import STARPROPERTIES
from posydon.config import PATH_TO_POSYDON
from posydon.utils.common_functions import check_state_of_star
from posydon.utils.common_functions import calculate_lambda_from_profile, calculate_Mejected_for_integrated_binding_energy
from posydon.utils.posydonwarning import Pwarn


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
         "mass_loss_during_CEE_merged": False # If False, then no mass loss from this step for a merged star
                                              # If True, then we remove mass according to the alpha-lambda prescription
                                              # assuming a final separation where the inner core RLOF starts.
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
    "H-rich_non_burning",
    "accreted_He_non_burning"
]


STAR_STATE_POST_HeMS = [
    'accreted_He_Core_He_burning',
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
            mass_loss_during_CEE_merged=MODEL['mass_loss_during_CEE_merged'],
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
            self.mass_loss_during_CEE_merged = mass_loss_during_CEE_merged

        self.verbose = verbose
        self.path_to_posydon = PATH_TO_POSYDON

    def __call__(self, binary):
        """Perform the CEE step for a BinaryStar object."""
        # Determine which star is the donor and which is the companion
        if binary.event in ["oCE1", "oDoubleCE1"]:
            donor_star = binary.star_1
            comp_star = binary.star_2
            star_to_merge = "1"
        elif binary.event in ["oCE2", "oDoubleCE2"]:
            donor_star = binary.star_2
            comp_star = binary.star_1
            star_to_merge = "2"
        else:
            raise ValueError("CEE does not apply if `event` is not "
                             "`oCE1`, `oDoubleCE1`, `oCE2`,or `oDoubleCE1`")
        # Check for double CE
        double_CE = binary.event in ["oDoubleCE1", "oDoubleCE2"]

        if self.verbose:
            print("binary.event : ", binary.event)
            print("common_envelope_efficiency : ",
                  self.common_envelope_efficiency)
            print("common_envelope_option_for_lambda : ",
                  self.common_envelope_option_for_lambda)

        # Check to make sure binary can go through a CE
        mergeable_donor = (donor_star.state in [
            'H-rich_Core_H_burning', 'stripped_He_Core_He_burning', 'accreted_He_Core_He_burning'])
        mergeable_HG_donor = (
            self.common_envelope_option_for_HG_star == "pessimistic"
            and donor_star.state in ['H-rich_Shell_H_burning'])
        if mergeable_donor or mergeable_HG_donor:
            # system merges
            binary.state = 'merged'
            binary.event = "oMerging" + star_to_merge
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
                lambda2_CE = np.nan
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
                core_definition_He_fraction=self.core_definition_He_fraction,
                mass_loss_during_CEE_merged=self.mass_loss_during_CEE_merged)
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
            Pwarn("Donor profile does not exist -- proceeding with "
                          "default_lambda alpha-CE prescription", "ApproximationWarning")
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

    def CEE_adjust_post_CE_core_masses(self, donor, mc1_i, rc1_i, donor_type,
            comp_star, mc2_i, rc2_i, comp_type, double_CE, verbose=False):
        """Calculate the post-common-envelope parameters.

        It determines the post-CE parameters based on the core properties.
        Note that these parameters may be updated in a subsequent step 
        depending on assumptions about whether and how the outermost layers
        of the CE core are removed.

        Parameters
        ----------
        donor : SingleStar object
            The donor star
        mc1_i : float
            core mass of the donor (in Msun)
        rc1_i : float
            core radius of the donor (in Rsun)
        donor_type : string
            String dictating whether the donor has a H-envelope or He-envelope
        comp_star : SingleStar object
            The companion star
        mc2_i : float
            core mass of the companion (in Msun)
        rc2_i : float
            core radius of the companion (in Rsun)
        comp_type : string
            String dictating whether the companion has a H-envelope or 
            He-envelope.
        double_CE : bool
            In case we have a double CE situation.
        verbose : bool
            In case we want information about the CEE.

        Returns
        -------
        m1c_f: float
            donor mass after the CE (in Msun)
        r1c_f: float
            donor radius after the CE (in Rsun)
        m2c_f: float
            companion mass after the CE (in Msun)
        r2c_f: float
            companion radius after the CE (in Rsun)
        """
        if donor_type == 'He_core':
            mc1_f = donor.he_core_mass
            rc1_f = donor.he_core_radius
        elif donor_type == "CO_core":
            mc1_f = donor.co_core_mass
            rc1_f = donor.co_core_radius
        else:
            raise ValueError(
                "donor_type in common envelope is not recognized. "
                "donor_type = {}, don't know how to proceed".
                format(donor_type))
        if mc1_f > mc1_i:
            mc1_f = mc1_i
            Pwarn("The final donor core mass (even after stable, postCEE MT)"
                  " is higher than the postCEE core mass. Now equalizing to "
                  "postCEE mass", "ApproximationWarning")
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
                Pwarn("The accretor's final core mass (even after "
                      "non-conservative stable, postCEE MT) is "
                      "higher that postCEE core mass. "
                      "Now equalizing to postCEE mass", "ApproximationWarning")
        if verbose:
            print("difference between m1 core mass defined by CEE step"
                    " / to the final one as pre CEE : ", mc1_f, mc1_i)
            print("difference between r1 core mass defined by CEE step"
                    " / to the final one as pre CEE : ", rc1_f, rc1_i)
            print("difference between m2 core mass defined by CEE step"
                    " / to the final one as pre CEE : ", mc2_f, mc2_i)
            print("difference between r2 core mass defined by CEE step"
                    " / to the final one as pre CEE : ", rc2_f, rc2_i)

        return mc1_f, rc1_f, mc2_f, rc2_f

    def CEE_core_replaced_noMT(self, donor, mc1_i, rc1_i, comp_star, mc2_i, 
                               rc2_i, separation_postCEE, verbose=False):
        """Calculate the post-common-envelope parameters upon exiting a CEE.

        This prescription assumes the he_core_mass/radius (or 
        co_core_mass/radius for CEE of stripped_He*) are replaced according 
        to the new core boundary used for CEE (based on 
        core_definition_H/He_fraction) but no other change in period after 
        succesful ejection at alpha-lambda prescription.

        Parameters
        ----------
        donor : SingleStar object
            The donor star
        mc1_i : float
            core mass of the donor (in Msun)
        rc1_i : float
            core radius of the donor (in Rsun)
        comp_star : SingleStar object
            The companion star
        mc2_i : float
            core mass of the companion (in Msun)
        rc2_i : float
            core radius of the companion (in Rsun)
        separation_postCEE : float
            binary's separation upon exiting the CEE (in cm)
        verbose : bool
            In case we want information about the CEE.

        Returns
        -------
        m1c_f: float
            donor mass after the CE (in Msun)
        r1c_f: float
            donor radius after the CE (in Rsun)
        m2c_f: float
            companion mass after the CE (in Msun)
        r2c_f: float
            companion radius after the CE (in Rsun)
        separation_f : float
            binary's final separation upon exiting the CE (in Rsun)
        orbital_period_f : float
            binary's final orbital period upon exiting the CE (in days)
        merger : bool
            whether the binary merged in the CE
        """
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

        # Check to see if the system has merged
        merger = cf.check_for_RLO(mc1_i, rc1_f, mc2_i, rc2_f,
                        separation_postCEE/const.Rsun, self.CEE_tolerance_err)

        if merger:
            if verbose:
                print("system merges within the CEE")
        else:
            if verbose:
                print("system survives CEE, the whole CE is ejected and the "
                      "new orbital separation for the cores is returned")

        # possible change of orbital period after CEE
        orbital_period_f = cf.orbital_period_from_separation(
            separation_postCEE / const.Rsun, mc1_i, mc2_i)

        # since the binary components' masses do not change, the output
        # separation is just the input separation in different units
        separation_f = separation_postCEE / const.Rsun

        return mc1_f, rc1_f, mc2_f, rc2_f, separation_f, orbital_period_f, merger

    def CEE_core_not_replaced_noMT(self, donor, mc1_i, rc1_i, 
                                   donor_type, comp_star, mc2_i, rc2_i, 
                                   comp_type, double_CE, separation_postCEE,
                                   verbose=False):
        """Calculate the post-common-envelope parameters upon exiting a CEE.

        This prescription assumes the he_core_mass/radius (or 
        co_core_mass/radius for CEE of stripped_He*) staying as preCEE and no 
        other change in period after succesful ejection at alpha-lambda 
        prescription.

        Parameters
        ----------
        donor : SingleStar object
            The donor star
        mc1_i : float
            core mass of the donor (in Msun)
        rc1_i : float
            core radius of the donor (in Rsun)
        comp_star : SingleStar object
            The companion star
        mc2_i : float
            core mass of the companion (in Msun)
        rc2_i : float
            core radius of the companion (in Rsun)
        separation_postCEE : float
            binary's separation upon exiting the CEE (in cm)
        verbose : bool
            In case we want information about the CEE.

        Returns
        -------
        m1c_f: float
            donor mass after the CE (in Msun)
        r1c_f: float
            donor radius after the CE (in Rsun)
        m2c_f: float
            companion mass after the CE (in Msun)
        r2c_f: float
            companion radius after the CE (in Rsun)
        separation_f : float
            binary's final separation upon exiting the CE (in Rsun)
        orbital_period_f : float
            binary's final orbital period upon exiting the CE (in days)
        merger : bool
            whether the binary merged in the CE
        """

        # First find the post-CE parameters for each star
        mc1_f, rc1_f, mc2_f, rc2_f = self.CEE_adjust_post_CE_core_masses(
            donor, mc1_i, rc1_i, donor_type, comp_star, mc2_i, rc2_i, 
            comp_type, double_CE, verbose=verbose)

        # Adopt the final orbital period to be the same as the initial period
        orbital_period_f = cf.orbital_period_from_separation(
            separation_postCEE / const.Rsun, mc1_i, mc2_i)

        # Calculate the orbital separation accordingly
        separation_f = cf.orbital_separation_from_period(orbital_period_f,
            mc1_f, mc2_f)

        # Check to see if the system has merged
        merger = cf.check_for_RLO(mc1_f, rc1_f, mc2_f, rc2_f,
            separation_f, self.CEE_tolerance_err)

        if merger:
            if verbose:
                print("system merges within the CEE")
        else:
            if verbose:
                print("system survives CEE, the whole CE is ejected and the "
                      "new orbital separation for the cores is returned")

        return mc1_f, rc1_f, mc2_f, rc2_f, separation_f, orbital_period_f, merger

    def CEE_core_not_replaced_stableMT(self, donor, mc1_i, rc1_i, 
                                       donor_type, comp_star, mc2_i, rc2_i, 
                                       comp_type, double_CE, 
                                       separation_postCEE, verbose=False):
        """Calculate the post-common-envelope parameters upon exiting a CEE.

        This prescription assumes the he_core_mass/radius (or 
        co_core_mass/radius for CEE of stripped_He*) stays as preCEE and after 
        a succesful ejection using the alpha-lambda prescription, we assume an 
        instantaneous stableMT phase (non-conservative, with mass lost from 
        accretor) from the donor (or from donor and simultaneously the accretor
        at double_CE), taking away the extra "core" mass as defined by the core
        boundary used for CEE (based on core_definition_H/He_fraction).

        Parameters
        ----------
        donor : SingleStar object
            The donor star
        mc1_i : float
            core mass of the donor (in Msun)
        rc1_i : float
            core radius of the donor (in Rsun)
        comp_star : SingleStar object
            The companion star
        mc2_i : float
            core mass of the companion (in Msun)
        rc2_i : float
            core radius of the companion (in Rsun)
        separation_postCEE : float
            binary's separation upon exiting the CEE (in cm)
        verbose : bool
            In case we want information about the CEE.

        Returns
        -------
        m1c_f: float
            donor mass after the CE (in Msun)
        r1c_f: float
            donor radius after the CE (in Rsun)
        m2c_f: float
            companion mass after the CE (in Msun)
        r2c_f: float
            companion radius after the CE (in Rsun)
        separation_f : float
            binary's final separation upon exiting the CE (in Rsun)
        orbital_period_f : float
            binary's final orbital period upon exiting the CE (in days)
        merger : bool
            whether the binary merged in the CE
        """
        # First find the post-CE parameters for each star
        mc1_f, rc1_f, mc2_f, rc2_f = self.CEE_adjust_post_CE_core_masses(
            donor, mc1_i, rc1_i, donor_type, comp_star, mc2_i, rc2_i,
            comp_type, double_CE, verbose=verbose)

        # calculate the orbital period after CEE (before adjustment)
        orbital_period_postCEE = cf.orbital_period_from_separation(
            separation_postCEE / const.Rsun, mc1_i, mc2_i)

        # First, check if merger happens before stable MT phase. Note this uses
        # the 1% core values to calculate the He-core radii, but it uses the 
        # post-CEE masses
        merger = cf.check_for_RLO(mc1_i, rc1_f, mc2_i, rc2_f,
            separation_postCEE / const.Rsun, self.CEE_tolerance_err)

        if merger:
            if verbose:
                print("system merges within the CEE, prior to stable MT")
            return mc1_i, rc1_f, mc2_i, rc2_f, \
                separation_postCEE / const.Rsun, orbital_period_postCEE, merger
        else:
            if verbose:
                print("system survives initial CEE detachment")
                print("donor core mass / core radius:", mc1_f, rc1_i)
                print("companion core mass / core radius:", mc2_f, rc2_i)


        # An assumed stable mass transfer case after postCEE with
        # fully non-conservative MT and mass lost from the vicinity
        # of the accretor:
        orbital_period_f = cf.period_change_stable_MT(
            orbital_period_postCEE, Mdon_i=mc1_i, Mdon_f=mc1_f,
            Macc_i=mc2_i, alpha=0.0, beta=1.0)

        # TODO: This seems wrong!
        if double_CE:
            # a reverse stable MT assumed to be happening
            # at the same time
            orbital_period_f = cf.period_change_stable_MT(
                orbital_period_f, Mdon_i=mc2_i, Mdon_f=mc2_f,
                Macc_i=mc1_i, alpha=0.0, beta=1.0)

        if verbose:
            print("during the assumed stable MT phase after postCE")
            print("with 'common_envelope_option_after_succ_CEE' :",
                "core_not_replaced_stableMT")
            print("the orbit changed from postCEE : ", orbital_period_postCEE)
            print("to : ", orbital_period_f)

        # Calculate the post-CEE separation
        separation_f = cf.orbital_separation_from_period(orbital_period_f,
            mc1_f, mc2_f)

        # Check once more to see if the system has merged during stable MT
        merger = cf.check_for_RLO(mc1_f, rc1_f, mc2_f, rc2_f, separation_f, 
            self.CEE_tolerance_err)

        if merger:
            if verbose:
                print("system merges within the CEE")
        else:
            if verbose:
                print("system survives CEE, the whole CE is ejected and the "
                      "new orbital separation for the cores is returned")

        return mc1_f, rc1_f, mc2_f, rc2_f, separation_f, orbital_period_f, merger

    def CEE_core_not_replaced_windloss(self, donor, mc1_i, rc1_i,
                                       donor_type, comp_star, mc2_i, rc2_i,
                                       comp_type, double_CE,
                                       separation_postCEE, verbose=False):
        """Calculate the post-common-envelope parameters upon exiting a CEE.

        This prescription assumes the he_core_mass/radius (or 
        co_core_mass/radius for CEE of stripped_He*) staying as preCEE and 
        after succesful ejection at alpha-lambda prescription, we assume a 
        instantaneous windloss phase from the donor (or from donor and 
        simultaneously the accretor at double_CE), taking away the extra 
        "core" mass as defined by the core boundary used for CEE (based on
        core_definition_H/He_fraction).

        Parameters
        ----------
        donor : SingleStar object
            The donor star
        mc1_i : float
            core mass of the donor (in Msun)
        rc1_i : float
            core radius of the donor (in Rsun)
        comp_star : SingleStar object
            The companion star
        mc2_i : float
            core mass of the companion (in Msun)
        rc2_i : float
            core radius of the companion (in Rsun)
        separation_postCEE : float
            binary's separation upon exiting the CEE (in cm)
        verbose : bool
            In case we want information about the CEE.

        Returns
        -------
        m1c_f: float
            donor mass after the CE (in Msun)
        r1c_f: float
            donor radius after the CE (in Rsun)
        m2c_f: float
            companion mass after the CE (in Msun)
        r2c_f: float
            companion radius after the CE (in Rsun)
        separation_f : float
            binary's final separation upon exiting the CE (in Rsun)
        orbital_period_f : float
            binary's final orbital period upon exiting the CE (in days)
        merger : bool
            whether the binary merged in the CE
        """
        # First find the post-CE parameters for each star
        mc1_f, rc1_f, mc2_f, rc2_f = self.CEE_adjust_post_CE_core_masses(
            donor, mc1_i, rc1_i, donor_type, comp_star, mc2_i, rc2_i,
            comp_type, double_CE, verbose=verbose)

        # calculate the orbital period after CEE (before adjustment)
        orbital_period_postCEE = cf.orbital_period_from_separation(
            separation_postCEE / const.Rsun, mc1_i, mc2_i)

        # First, check if merger happens before stable MT phase. Note this uses
        # the 1% core values to calculate the He-core radii, but it uses the 
        # post-CEE masses
        merger = cf.check_for_RLO(mc1_i, rc1_f, mc2_i, rc2_f,
            separation_postCEE / const.Rsun, self.CEE_tolerance_err)

        if merger:
            if verbose:
                print("system merges within the CEE, prior to stable MT")
            return mc1_i, rc1_f, mc2_i, rc2_f, \
                separation_postCEE / const.Rsun, orbital_period_postCEE, merger
        else:
            if verbose:
                print("system survives initial CEE detachment")
                print("donor core mass / core radius:", mc1_f, rc1_i)
                print("companion core mass / core radius:", mc2_f, rc2_i)


        # An assumed wind loss after postCEE with mass lost
        # from the vicinity of the donor
        orbital_period_f = cf.period_change_stable_MT(
            orbital_period_postCEE, Mdon_i=mc1_i, Mdon_f=mc1_f,
            Macc_i=mc2_i, alpha=1.0, beta=0.0)

        if double_CE:
            # a wind mass loss from the 2nd star assumed to be
            # happening at the same time
            orbital_period_f = cf.period_change_stable_MT(
                orbital_period_f, Mdon_i=mc2_i, Mdon_f=mc2_f,
                Macc_i=mc1_i, alpha=1.0, beta=0.0)
        if verbose:
            print("during the assumed windloss phase after postCE")
            print("with 'common_envelope_option_after_succ_CEE' :",
                "core_not_replaced_windloss")
            print("the orbit changed from postCEE : ", orbital_period_postCEE)
            print("to : ", orbital_period_f)

        # Calculate the post-CEE separation
        separation_f = cf.orbital_separation_from_period(orbital_period_f,
            mc1_f, mc2_f)

        # Check once more to see if the system has merged during stable MT
        merger = cf.check_for_RLO(mc1_f, rc1_f, mc2_f, rc2_f, separation_f,
            self.CEE_tolerance_err)

        if merger:
            if verbose:
                print("system merges within the CEE")
        else:
            if verbose:
                print("system survives CEE, the whole CE is ejected and the "
                      "new orbital separation for the cores is returned")

        return mc1_f, rc1_f, mc2_f, rc2_f, separation_f, orbital_period_f, merger

    def CEE_simple_alpha_prescription(
            self, binary, donor, comp_star, lambda1_CE, mc1_i, rc1_i,
            donor_type, lambda2_CE, mc2_i, rc2_i, comp_type, double_CE=False,
            verbose=False, common_envelope_option_after_succ_CEE=MODEL[
                'common_envelope_option_after_succ_CEE'],
            core_definition_H_fraction=MODEL['core_definition_H_fraction'],
            core_definition_He_fraction=MODEL['core_definition_He_fraction'],
            mass_loss_during_CEE_merged=MODEL['mass_loss_during_CEE_merged']):
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
        mass_loss_during_CEE_merged: Boolean
             If False, then no mass loss from this step for a merged CEE
             If True, then we remove mass according to the alpha-lambda prescription
              assuming a final separation where the inner core(s) RLOF starts, and the same lambda(s)

        """
        # Get star properties
        m1_i = donor.mass
        radius1 = 10**donor.log_R
        state1_i = donor.state
        m2_i = comp_star.mass
        radius2 = 10**comp_star.log_R
        # he_core_radius1 = donor.he_core_radius

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
        if pd.isna(lambda1_CE):
            ebind_i = 0.0
        else:
            ebind_i = (-const.standard_cgrav / lambda1_CE
                       * (m1_i * const.Msun * (m1_i - mc1_i) * const.Msun)
                       / (radius1 * const.Rsun))
        if (double_CE and pd.notna(lambda2_CE)):
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
        if separation_postCEE < -self.CEE_tolerance_err:
            raise ValueError("CEE problem, negative postCEE separation")

        if verbose:
            print("CEE alpha-lambda prescription")
            print("ebind_i", ebind_i)
            print("eorb_i", eorb_i)
            print("eorb_postCEE", eorb_postCEE)
            print("DEorb", eorb_postCEE - eorb_i)
            print("separation_i in Rsun", separation_i/const.Rsun)
            print("separation_postCEE in Rsun", separation_postCEE/const.Rsun)


        # Calculate the post-CE binary properties
        if common_envelope_option_after_succ_CEE == "CEE_core_replaced_noMT":
            mc1_f, rc1_f, mc2_f, rc2_f, separation_f, orbital_period_f, merger \
                    = self.CEE_core_replaced_noMT(donor, mc1_i, rc1_i, comp_star,
                                             mc2_i, rc2_i, separation_postCEE,
                                             verbose)
        elif common_envelope_option_after_succ_CEE == "CEE_core_not_replaced_noMT":
            mc1_f, rc1_f, mc2_f, rc2_f, separation_f, orbital_period_f, merger \
                    = self.CEE_core_not_replaced_noMT(donor, mc1_i, rc1_i,
                                                 donor_type, comp_star, mc2_i,
                                                 rc2_i, comp_type, double_CE,
                                                 separation_postCEE, verbose)
        elif common_envelope_option_after_succ_CEE == "CEE_core_not_replaced_stableMT":
            mc1_f, rc1_f, mc2_f, rc2_f, separation_f, orbital_period_f, merger \
                    = self.CEE_core_not_replaced_stableMT(donor, mc1_i, rc1_i,
                                                     donor_type, comp_star,
                                                     mc2_i, rc2_i, comp_type,
                                                     double_CE,
                                                     separation_postCEE,
                                                     verbose)
        elif common_envelope_option_after_succ_CEE == "CEE_core_not_replaced_windloss":
            mc1_f, rc1_f, mc2_f, rc2_f, separation_f, orbital_period_f, merger \
                    = self.CEE_core_not_replaced_windloss(donor, mc1_i, rc1_i,
                                                     donor_type, comp_star,
                                                     mc2_i, rc2_i, comp_type,
                                                     double_CE,
                                                     separation_postCEE,
                                                     verbose)
        else:
            raise ValueError(
                "Not accepted option in common_envelope_option_after_succ_"
                "CEE = {}, dont know how to proceed".
                format(common_envelope_option_after_succ_CEE))

        # Adjust stellar and binary parameters depending on whether the system
        # mergers in the CE or not
        if merger:
            # Calculate the amount of mass lost during the merger
            if mass_loss_during_CEE_merged:
                Pwarn("NOTE:This is legacy code from a previous "
                      "implementation. It is not currently supported "
                      "and may lead to unexpected behavior", 
                      "UnsupportedModelWarning")
                Mejected_donor, Mejected_comp = \
                    self.CEE_adjust_mass_loss_during_CEE_merged(donor,
                        comp_star, double_CE, verbose=verbose)
            else:
                Mejected_donor = 0.0
                Mejected_comp = 0.0

            # Adjust the binary and single star properties due to merger
            self.CEE_adjust_binary_upon_merger(binary, donor, comp_star, m1_i,
                m2_i, donor_type, Mejected_donor, Mejected_comp, verbose)
        else:
            # Adjust the binary and single star properties due to ejection
            self.CEE_adjust_binary_upon_ejection(binary, donor, mc1_f, rc1_f,
                donor_type, comp_star, mc2_f, rc2_f, comp_type, double_CE,
                separation_f, orbital_period_f,
                common_envelope_option_after_succ_CEE,
                core_definition_H_fraction,
                core_definition_He_fraction,
                verbose)

        return



    def CEE_adjust_binary_upon_ejection(self, binary, donor, mc1_f, rc1_f,
                donor_type, comp_star, mc2_f, rc2_f, comp_type, double_CE,
                separation_f, orbital_period_f, 
                common_envelope_option_after_succ_CEE, 
                core_definition_H_fraction, 
                core_definition_He_fraction,
                verbose=False):
        """Update the binary and component stars upon exiting a CEE.

        The binary's parameters (orbital period, separation, state, etc.) are
        updated along with the donor's (and in the case of a double CE the
        companion's as well) parameters are also updated. Note that certain 
        parameters are set to np.nan if they are irrelevant for the CE.

        Parameters
        ----------
        binary: BinaryStar object
            The binary system
        donor : SingleStar object
            The donor star
        mc1_f : float
            final core mass of the donor (in Msun)
        rc1_f : float
            final core radius of the donor (in Rsun)
        donor_type : string
            descriptor for the stellar type of the donor
        comp_star : SingleStar object
            The companion star
        mc2_f : float
            core mass of the companion (in Msun)
        rc2_f : float
            core radius of the companion (in Rsun)
        comp_type : string
            descriptor for the stellar type of the companion
        double_CE : bool
            whether the CEE is a double CE or not
        separation_f : float
            binary's separation upon exiting the CEE (in Rsun)
        orbital_period_f : float
            binary's orbital period upon exiting the CEE (in days)
        common_envelope_option_after_succ_CEE : string
            which type of post-common envelope evolution is used to remove
            the final layers around the helium core
        core_definition_H_fraction : float
            the fractional abundance of H defining the He core (0.3, 0.1, or 
            0.01)
        core_definition_He_fraction : float
            the fractional abundance of He defining the CO core (typically 0.1)
        verbose : bool
            In case we want information about the CEE.
        """
        # Adjust binary properties
        binary.separation = separation_f
        binary.orbital_period = orbital_period_f
        binary.eccentricity = 0.0
        binary.state = 'detached'
        binary.event = None

        if verbose:
            print("CEE succesfully ejected")
            print("new orbital period = ", binary.orbital_period)
            print("binary event : ", binary.event)
            print("double CEE : ", double_CE)

        # TODO: Check this!
        # Adjust the binary event to CC1/CC2 if star will explode
        # if state1_i == 'stripped_He_Central_He_depleted':
        #     if donor == binary.star_1:
        #         binary.event = 'CC1'
        #     elif donor == binary.star_2:
        #         binary.event = 'CC2'


        # Set binary values that are unchanged in CE step to np.nan
        for key in BINARYPROPERTIES:

            # the binary attributes that are changed in the CE step
            if key in ["separation", "orbital_period",
                       "eccentricity", "state", "event"]:
                continue

            # the binary attributes that keep the same value from the
            # previous step
            if key in ["time", "V_sys", "mass_transfer_case",
                       "nearest_neighbour_distance"]:
                continue

            # the rest become np.nan
            setattr(binary, key, np.nan)

        # Create list of stars that need updated parameters. If normal CE,
        # then only donor has updated properties
        stars = [donor]
        star_types = [donor_type]
        core_masses = [mc1_f]
        core_radii = [rc1_f]

        # If a double CE, then the companion's properties are updated, too
        if double_CE:
            stars.append(comp_star)
            star_types.append(comp_type)
            core_masses.append(mc2_f)
            core_radii.append(rc2_f)

        # Adjust stellar properties
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
                #TODO: update 'total_mass_h1' and 'total_mass_he4'
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
                #TODO: update 'total_mass_h1' and 'total_mass_he4'
            elif star_type == "not_giant_companion":
                continue
            else:
                raise ValueError("Unrecognized star type:", star_type)

            # Update state of star
            # TODO shouldn't stellar states be updated after the parameters are updated?
            state_old = star.state
            star.state = check_state_of_star(star)
            if verbose:
                print("star state of before/ after CE : ",
                        state_old, star.state)

            # Set star values that are unchanged in CE step to np.nan
            for key in STARPROPERTIES:

                # the singlestar attributes that are changed
                # in the CE step
                if key in attributes_changing:
                    continue

                # the singlestar attributes that keep the same value
                # from previous step if not in attributes_changing
                if key in [
                        "surface_h1", "surface_he4", "log_LH",
                        "log_LHe", "he_core_mass", "he_core_radius",
                        "co_core_mass", "co_core_radius", "center_h1",
                        "center_he4", "center_c12", "center_o16",
                        "center_n14", "metallicity", "log_LZ",
                        "log_Lnuc", "c12_c12", "center_gamma",
                        "avg_c_in_c_core"]:
                    continue

                # the rest become NaN's
                setattr(star, key, np.nan)
        return

    def CEE_adjust_mass_loss_during_CEE_merged(self, donor, comp_star, 
                                               double_CE, verbose=False):
        """ Calculate the amount of mass lost during a stellar merger in a CEE

        The binary's parameters (orbital period, separation, state, etc.) are
        updated along with the donor's (and in the case of a double CE the
        companion's as well) parameters are also updated. Note that certain 
        parameters are set to np.nan if they are irrelevant for the CE.

        Parameters
        ----------
        binary: BinaryStar object
            The binary system
        donor : SingleStar object
            The donor star
        mc1_f : float
            final core mass of the donor (in Msun)
        rc1_f : float
            final core radius of the donor (in Rsun)
        donor_type : string
            descriptor for the stellar type of the donor
        comp_star : SingleStar object
            The companion star
        mc2_f : float
            core mass of the companion (in Msun)
        rc2_f : float
            core radius of the companion (in Rsun)
        comp_type : string
            descriptor for the stellar type of the companion
        double_CE : bool
            whether the CEE is a double CE or not
        separation_f : float
            binary's separation upon exiting the CEE (in Rsun)
        orbital_period_f : float
            binary's orbital period upon exiting the CEE (in days)
        common_envelope_option_after_succ_CEE : string
            which type of post-common envelope evolution is used to remove
            the final layers around the helium core
        core_definition_H_fraction : float
            the fractional abundance of H defining the He core (0.3, 0.1, or 
            0.01)
        core_definition_He_fraction : float
            the fractional abundance of He defining the CO core (typically 0.1)
        verbose : bool
            In case we want information about the CEE.
        """
        # we calculate the ejected mass from part of the common envelope, using
        # a_f = separation_postCEE so that one of the cores (or MS star) is 
        # filling its inner Roche lobe, and assuming that 
        # lambda(Menvelope) ~ lamda(Mejected) 
        # although Mejected < Menvelope (e.g. see Fig1 of Dewi+Tauris2000)

        # First, make sure we have information to calculate ejected mass
        if not double_CE and donor.profile is None:
            Mejected_donor = 0.0
            Mejected_comp = 0.0
            Pwarn("mass_loss_during_CEE_merged == True, but no profile found "
                            "for the donor star. Proceeding with no partial mass ejection.", "ApproximationWarning")
            if verbose:
                print("Mejecta_donor = 0 in Msun compared to its initial "
                      "envelope =",  m1_i  - mc1_i)
                print("Mejecta_comp = 0 in Msun compared to its initial "
                      "envelope =",  m2_i  - mc2_i)
            return Mejected_donor, Mejected_comp
        elif double_CE and (donor.profile is None or comp_star.profile is None):
            Mejected_donor = 0.0
            Mejected_comp = 0.0
            Pwarn("mass_loss_during_CEE_merged == True, but not profile found "
                            "for the donor or companion star in double_CE. Proceeding with no partial mass ejection.", 
                            "ApproximationWarning")
            if verbose:
                print("Mejecta_donor = 0 in Msun compared to its initial "
                      "envelope =",  m1_i  - mc1_i)
                print("Mejecta_comp = 0 in Msun compared to its initial "
                      "envelope =",  m2_i  - mc2_i)
            return Mejected_donor, Mejected_comp
        else:
            pass

        separation_for_inner_RLO1 = rc1_i / cf.roche_lobe_radius(mc1_i, mc2_i, a_orb=1)
        separation_for_inner_RLO2 = rc2_i / cf.roche_lobe_radius(mc2_i, mc1_i, a_orb=1)

        separation_before_merger = max( separation_for_inner_RLO1, separation_for_inner_RLO2 ) * const.Rsun

        if verbose:
            print("separation_before_merger (for the calculation of Mejected "
                  "for the merger): ", separation_before_merger/ const.Rsun, 
                  "in Rsun")
            print("which is the max of RLO1 or RLO2 of the inner cores: ", 
                  separation_for_inner_RLO1 ,  separation_for_inner_RLO2, 
                  "in Rsun")

        E_orb_used_up_to_inner_RLOF = alpha_CE * ( - const.standard_cgrav * m1_i * const.Msun * m2_i * const.Msun/(2. * separation_i) + \
                    const.standard_cgrav * mc1_i  * const.Msun * mc2_i  * const.Msun/(2. * separation_before_merger) )

        # We assume "lambda_from_profile_gravitational_plus_internal_minus_recombination"
        if not double_CE:
            Mejected_donor = \
                calculate_Mejected_for_integrated_binding_energy(donor.profile,
                    E_orb_used_up_to_inner_RLOF, mc1_i, rc1_i, m1_i, radius1)

        else: # in double_CE

            # Assuming that the ratio of orbital energy used for the partial 
            # ejection of each (common) envelope is the same as the ratio of 
            # the initial envelope masses:

            WF = (m1_i  - mc1_i)/ (m2_i - mc2_i  +  m1_i  - mc1_i) # weight factor for 1:
                                                # = Mdonor,envelope / Mcomp,envelope
            Eorb_for_partial_ej_1 = E_orb_used_up_to_inner_RLOF * WF
            Mejected_donor = calculate_Mejected_for_integrated_binding_energy(donor.profile, Eorb_for_partial_ej_1, mc1_i, rc1_i, m1_i, radius1)

            Eorb_for_partial_ej_2 = E_orb_used_up_to_inner_RLOF * (1. - WF)
            Mejected_comp = calculate_Mejected_for_integrated_binding_energy(comp_star.profile, Eorb_for_partial_ej_2, mc2_i, rc2_i, m2_i, radius2)

        if verbose:
            print("Mejecta_donor = ", Mejected_donor, "in Msun compared to its initial envelope =",  m1_i  - mc1_i)
            print("Mejecta_comp = ", Mejected_comp, "in Msun compared to its initial envelope =",  m2_i  - mc2_i)
        '''
        if not ( (Mejected_donor <= m1_i  - mc1_i) and (Mejected_comp <= m2_i  - mc2_i) ):
            raise Exception("Mejected_donor in double CEE with mass_loss_during_CEE_merged is found more than Menvelope")
        if not (Mejected_donor >= 0.):
            raise Exception("The root of the equation in double CEE with mass_loss_during_CEE_merged is found negative")
        '''

        if (Mejected_donor > m1_i  - mc1_i) or (Mejected_comp > m2_i  - mc2_i):
            Mejected_donor = (m1_i  - mc1_i) -0.01 # at least this value of envelope is left.
            Mejected_comp = (m2_i  - mc2_i) -0.01
            Pwarn("M_ejected of at least one star in double CEE is found to be more than the initial envelope. "
                            "Reducing both to their initial_envelope - 0.01 Msun", "ApproximationWarning")

        if verbose:
            print("The mass loss during merging_CEE is: ", mass_loss_during_CEE_merged)
            print("so m1_i , m1_f = ", m1_i, donor.mass)
            print("so m2_i , m2_f = ", m2_i, comp_star.mass)


        return Mejected_donor, Mejected_comp

    def CEE_adjust_binary_upon_merger(self, binary, donor, comp_star, m1_i,
                                      m2_i, donor_type, Mejected_donor,
                                      Mejected_comp, verbose=False):
        """Update the binary and component stars upon merging within a CEE.

        The binary's state and event are updated along with the donor and 
        companion star masses and radii corresponding to a merger event.

        Parameters
        ----------
        binary: BinaryStar object
            The binary system
        donor : SingleStar object
            The donor star
        comp_star : SingleStar object
            The companion star
        m1_i : float
            mass of the donor upon entering a CE (in Msun)
        m2_i : float
            mass of the companion upon entering a CE (in Msun)
        donor_type : string
            descriptor for the stellar type of the donor
        Mejected_donor : float
            how much mass is ejected from the donor upon merger (in Msun)
        Mejected_comp : float
            how much mass is ejected from the companion upon merger (in Msun)
        verbose : bool
            In case we want information about the CEE.
        """
        # system merges
        binary.state = 'merged'
        if binary.event in ["oCE1", "oDoubleCE1"]:
            binary.event = "oMerging1"
        if binary.event in ["oCE2", "oDoubleCE2"]:
            binary.event = "oMerging2"

        if verbose:
            print("system merges due to one of the two star's core filling"
                    "its RL")
            # print("Rdonor core vs RLdonor core = ", core_radius1, RL1)
            # print("Rcompanion vs RLcompanion= ", rc2_i, RL2)

        donor.mass = m1_i - Mejected_donor
        donor.log_R = np.nan
        comp_star.mass = m2_i - Mejected_comp
        comp_star.log_R = np.nan
        if donor_type == 'CO_core':
            donor.he_core_mass = m1_i - Mejected_donor
            donor.he_core_radius = np.nan
            comp_star.he_core_mass = m2_i - Mejected_comp
            comp_star.he_core_radius = np.nan

        return
