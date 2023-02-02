"""Merging and isolated evolution step."""


__authors__ = [
    "Emmanouil Zapartas <ezapartas@gmail.com>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>"
]


import os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize
from scipy.optimize import root

from posydon.utils.data_download import PATH_TO_POSYDON_DATA
from posydon.binary_evol.binarystar import BINARYPROPERTIES
from posydon.binary_evol.singlestar import STARPROPERTIES
from posydon.interpolation import GRIDInterpolator
from posydon.interpolation.data_scaling import DataScaler
from posydon.utils.common_functions import (
    bondi_hoyle,
    orbital_period_from_separation,
    orbital_separation_from_period,
    roche_lobe_radius,
    check_state_of_star,
    PchipInterpolator2
)
from posydon.binary_evol.flow_chart import (STAR_STATES_CC)
import posydon.utils.constants as const
from posydon.binary_evol.step_detached import detached_step
from posydon.binary_evol.step_isolated import isolated_step

import warnings

from posydon.binary_evol.flow_chart import (STAR_STATES_ALL,
    STAR_STATES_CO,
    STAR_STATES_H_RICH,
    STAR_STATES_HE_RICH,
    STAR_STATES_NOT_CO
    )



LIST_ACCEPTABLE_STATES_FOR_HMS = ["H-rich_Core_H_burning"]
LIST_ACCEPTABLE_STATES_FOR_HeMS = ["stripped_He_Core_He_burning"]

LIST_ACCEPTABLE_STATES_FOR_POSTMS = STAR_STATES_H_RICH.copy()
[LIST_ACCEPTABLE_STATES_FOR_POSTMS.remove(x) for x in LIST_ACCEPTABLE_STATES_FOR_HMS]

LIST_ACCEPTABLE_STATES_FOR_POSTHeMS = STAR_STATES_HE_RICH.copy()
[LIST_ACCEPTABLE_STATES_FOR_POSTHeMS.remove(x) for x in LIST_ACCEPTABLE_STATES_FOR_HeMS]



class merged_step(isolated_step):
        """
        Prepare a merging star to do an an isolated_step)
        """

    def __init__(self,
        grid_name_Hrich=None,
        grid_name_stripedHe=None,
        path=PATH_TO_POSYDON_DATA,
        merger_critical_rot = 0.4,
        rel_mass_lost_HMS_HMS = 0.1,
        *args, **kwargs):

        super().__init__(
        grid_name_Hrich=grid_name_Hrich,
        grid_name_stripedHe=grid_name_stripedHe,
        *args,
        **kwargs)

    def merged_star_properties(star_base,comp):
        """
        Make assumptions about the core/total mass of the star of a merged product.

        Similar to the table of merging in BSE

        star_base: Single Star
            is our base star that engulfs its companions. The merged star will have this star as a base
        comp: Single Star
            is the star that is engulfed
        """
        #by default the stellar attributes that keep the same value from the
        merged_star = star_base

        s1 = star_base.state
        s2 = comp.state


        def mass_weighted_avg(star1=star_base,star2=comp, abundance_name="center_h1", mass_weight_name="mass"):
            A1 = getattr(star1, abundance_name)
            A2 = getattr(star2, abundance_name)
            M1 = getattr(star1, mass_weight_name)
            M2 = getattr(star2, mass_weight_name)
            return (A1*M1 + A2*M2 ) / (M1+M2)

        if ( s1 is in LIST_ACCEPTABLE_STATES_FOR_HMS
            and s2 is in LIST_ACCEPTABLE_STATES_FOR_HMS):
                #these stellar attributes change value
                merged_star.mass = (star_base.mass + comp.mass) * (1.-rel_mass_lost_HMS_HMS)

                merged_star.center_h1 = mass_weighted_avg()
                merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4")
                merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12")
                merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14")
                merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16")

                #TODO: should I check if the abundaces above end up in ~1 (?)

                for key in STARPROPERTIES:
                    # these stellar attributes become np.nan
                    for substring in ["log_", "lg_", "surface","surf_", "conv_", "lambda", "profile"]:
                        if (substring in key) :
                            setattr(merged_star, key, np.nan)
                        if key in [ "c12_c12", "center_gamma",
                                   "avg_c_in_c_core", "total_moment_of_inertia", "spin", "envelope_binding_energy"]:
                            setattr(merged_star, key, np.nan)


        elif (s1 is in LIST_ACCEPTABLE_STATES_FOR_POSTMS:
            and s2 is in LIST_ACCEPTABLE_STATES_FOR_HMS):

                merged_star.mass = star_base.mass + comp.mass

                for key in STARPROPERTIES:
                    # these stellar attributes become np.nan
                    for substring in ["log_", "lg_", "surface","surf_", "conv_", "lambda", "profile"]:
                        if (substring in key) :
                            setattr(merged_star, key, np.nan)
                        if key in [ "c12_c12", "center_gamma",
                                   "avg_c_in_c_core", "total_moment_of_inertia", "spin"]:
                            setattr(merged_star, key, np.nan)

                merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!

        # as above but opposite stars
        elif (s1 is in LIST_ACCEPTABLE_STATES_FOR_HMS:
            and s2 is in LIST_ACCEPTABLE_STATES_FOR_POSTMS):

                merged_star = comp
                merged_star.mass = star_base.mass + comp.mass

                for key in STARPROPERTIES:
                    # these stellar attributes become np.nan
                    for substring in ["log_", "lg_", "surface","surf_", "conv_", "lambda", "profile"]:
                        if (substring in key) :
                            setattr(merged_star, key, np.nan)
                        if key in [ "c12_c12", "center_gamma",
                                   "avg_c_in_c_core", "total_moment_of_inertia", "spin"]:
                            setattr(merged_star, key, np.nan)

                merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!


        elif (s1 is in LIST_ACCEPTABLE_STATES_FOR_POSTMS:
            and s2 is in LIST_ACCEPTABLE_STATES_FOR_POSTMS):

                # add total and core masses
                for key in ["mass", "he_core_mass", "c_core_mass", "o_core_mass", "co_core_mass"]:
                    current = getattr(merged_star, key) + getattr(comp, key)
                    setattr(merged_star, key,current)

                # weighted central abundances if merging cores. Else only from star_base
                if star_base.co_core_mass == 0 and comp.co_core_mass == 0: # two stars with Helium cores
                    merged_star.center_h1 = mass_weighted_avg(mass_weight_name="he_core_mass")
                    merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight_name="he_core_mass")
                    merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight_name="he_core_mass")
                    merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight_name="he_core_mass")
                    merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight_name="he_core_mass")
                elif (star_base.co_core_mass > 0 and comp.co_core_mass == 0): # star_base with CO core and the comp has a He core
                    continue # the central abundances are kept as the ones of star_base
                elif (comp.co_core_mass > 0 and star_base.co_core_mass == 0): # comp with CO core and the star_base has a He core
                    merged_star.center_h1 = comp.center_h1
                    merged_star.center_he4 = comp.center_he4
                    merged_star.center_c12 = comp.center_c12
                    merged_star.center_n14 = comp.center_n14
                    merged_star.center_o16 = comp.center_o16
                elif (star_base.co_core_mass > 0 and comp.co_core_mass > 0):
                    merged_star.center_h1 = mass_weighted_avg(mass_weight_name="he_core_mass")
                    merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight_name="he_core_mass")
                    merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight_name="he_core_mass")
                    merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight_name="he_core_mass")
                    merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight_name="he_core_mass")
                else:
                    warnings.warn("weird compbination of CO core masses during merging")

                for key in STARPROPERTIES:
                    # these stellar attributes become np.nan
                    for substring in ["log_", "lg_", "surface","surf_", "conv_", "lambda", "profile"]:
                        if (substring in key) :
                            setattr(merged_star, key, np.nan)
                        if key in [ "c12_c12", "center_gamma",
                                   "avg_c_in_c_core", "total_moment_of_inertia", "spin"]:
                            setattr(merged_star, key, np.nan)
                merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!

        elif (s1 is in LIST_ACCEPTABLE_STATES_FOR_POSTMS:
            and s2 is in LIST_ACCEPTABLE_STATES_FOR_HeMS):

                # add total and core masses
                for key in ["mass", "he_core_mass", "c_core_mass", "o_core_mass", "co_core_mass"]:
                    current = getattr(merged_star, key) + getattr(comp, key)
                    setattr(merged_star, key,current)

                for key in STARPROPERTIES:
                    # these stellar attributes become np.nan
                    for substring in ["log_", "lg_", "surface","surf_", "conv_", "lambda", "profile"]:
                        if (substring in key) :
                            setattr(merged_star, key, np.nan)
                        if key in [ "c12_c12", "center_gamma",
                                   "avg_c_in_c_core", "total_moment_of_inertia", "spin"]:
                            setattr(merged_star, key, np.nan)
                merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!

        # as above but opposite stars
        elif (s1 is in  LIST_ACCEPTABLE_STATES_FOR_HeMS :
            and s2 is in LIST_ACCEPTABLE_STATES_FOR_POSTMS):

                merged_star = comp

                # add total and core masses
                for key in ["mass", "he_core_mass", "c_core_mass", "o_core_mass", "co_core_mass"]:
                    current = getattr(merged_star, key) + getattr(star_base, key)
                    setattr(merged_star, key,current)

                for key in STARPROPERTIES:
                    # these stellar attributes become np.nan
                    for substring in ["log_", "lg_", "surface","surf_", "conv_", "lambda", "profile"]:
                        if (substring in key) :
                            setattr(merged_star, key, np.nan)
                        if key in [ "c12_c12", "center_gamma",
                                   "avg_c_in_c_core", "total_moment_of_inertia", "spin"]:
                            setattr(merged_star, key, np.nan)
                merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!

        elif (s1 is in LIST_ACCEPTABLE_STATES_FOR_POSTMS:
            and s2 is in LIST_ACCEPTABLE_STATES_FOR_POSTHeMS):

                # add total and core masses
                for key in ["mass", "he_core_mass", "c_core_mass", "o_core_mass", "co_core_mass"]:
                    current = getattr(merged_star, key) + getattr(comp, key)
                    setattr(merged_star, key,current)

                # weighted central abundances if merging cores. Else only from star_base
                if star_base.co_core_mass == 0 and comp.co_core_mass == 0: # two stars with Helium cores
                    merged_star.center_h1 = mass_weighted_avg(mass_weight_name="he_core_mass")
                    merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight_name="he_core_mass")
                    merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight_name="he_core_mass")
                    merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight_name="he_core_mass")
                    merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight_name="he_core_mass")
                elif (star_base.co_core_mass > 0 and comp.co_core_mass == 0): # star_base with CO core and the comp has a He core
                    continue # the central abundances are kept as the ones of star_base
                elif (comp.co_core_mass > 0 and star_base.co_core_mass == 0): # comp with CO core and the star_base has a He core
                    merged_star.center_h1 = comp.center_h1
                    merged_star.center_he4 = comp.center_he4
                    merged_star.center_c12 = comp.center_c12
                    merged_star.center_n14 = comp.center_n14
                    merged_star.center_o16 = comp.center_o16
                elif (star_base.co_core_mass > 0 and comp.co_core_mass > 0):
                    merged_star.center_h1 = mass_weighted_avg(mass_weight_name="co_core_mass")
                    merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight_name="co_core_mass")
                    merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight_name="co_core_mass")
                    merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight_name="co_core_mass")
                    merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight_name="co_core_mass")
                else:
                    warnings.warn("weird compbination of CO core masses during merging")

                for key in STARPROPERTIES:
                    # these stellar attributes become np.nan
                    for substring in ["log_", "lg_", "surface","surf_", "conv_", "lambda", "profile"]:
                        if (substring in key) :
                            setattr(merged_star, key, np.nan)
                        if key in [ "c12_c12", "center_gamma",
                                   "avg_c_in_c_core", "total_moment_of_inertia", "spin"]:
                            setattr(merged_star, key, np.nan)
                merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!

        # as above but the opposite stars
        elif (s1 is in LIST_ACCEPTABLE_STATES_FOR_POSTHeMS:
            and s2 is in LIST_ACCEPTABLE_STATES_FOR_POSTMS):

                merged_star = comp
                # add total and core masses
                for key in ["mass", "he_core_mass", "c_core_mass", "o_core_mass", "co_core_mass"]:
                    current = getattr(merged_star, key) + getattr(star_base, key)
                    setattr(merged_star, key,current)

                # weighted central abundances if merging cores. Else only from star_base
                if star_base.co_core_mass == 0 and comp.co_core_mass == 0: # two stars with Helium cores
                    merged_star.center_h1 = mass_weighted_avg(mass_weight_name="he_core_mass")
                    merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight_name="he_core_mass")
                    merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight_name="he_core_mass")
                    merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight_name="he_core_mass")
                    merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight_name="he_core_mass")
                elif (star_base.co_core_mass > 0 and comp.co_core_mass == 0): # star_base with CO core and the comp has a He core
                    continue # the central abundances are kept as the ones of star_base
                elif (comp.co_core_mass > 0 and star_base.co_core_mass == 0): # comp with CO core and the star_base has a He core
                    merged_star.center_h1 = comp.center_h1
                    merged_star.center_he4 = comp.center_he4
                    merged_star.center_c12 = comp.center_c12
                    merged_star.center_n14 = comp.center_n14
                    merged_star.center_o16 = comp.center_o16
                elif (star_base.co_core_mass > 0 and comp.co_core_mass > 0):
                    merged_star.center_h1 = mass_weighted_avg(mass_weight_name="co_core_mass")
                    merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight_name="co_core_mass")
                    merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight_name="co_core_mass")
                    merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight_name="co_core_mass")
                    merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight_name="co_core_mass")
                else:
                    warnings.warn("weird compbination of CO core masses during merging")

                for key in STARPROPERTIES:
                    # these stellar attributes become np.nan
                    for substring in ["log_", "lg_", "surface","surf_", "conv_", "lambda", "profile"]:
                        if (substring in key) :
                            setattr(merged_star, key, np.nan)
                        if key in [ "c12_c12", "center_gamma",
                                   "avg_c_in_c_core", "total_moment_of_inertia", "spin"]:
                            setattr(merged_star, key, np.nan)
                merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!


            elif (s1 is in STAR_STATES_HE_RICH:
                and s2 is in STAR_STATES_HE_RICH):

                    # add total and core masses
                    for key in ["mass", "he_core_mass", "c_core_mass", "o_core_mass", "co_core_mass"]:
                        current = getattr(merged_star, key) + getattr(comp, key)
                        setattr(merged_star, key,current)

                    # weighted central abundances if merging cores. Else only from star_base
                    if star_base.co_core_mass == 0 and comp.co_core_mass == 0: # two stars with Helium cores
                        merged_star.center_h1 = mass_weighted_avg(mass_weight_name="he_core_mass")
                        merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight_name="he_core_mass")
                        merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight_name="he_core_mass")
                        merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight_name="he_core_mass")
                        merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight_name="he_core_mass")
                    elif (star_base.co_core_mass > 0 and comp.co_core_mass == 0): # star_base with CO core and the comp has a He core
                        continue # the central abundances are kept as the ones of star_base
                    elif (comp.co_core_mass > 0 and star_base.co_core_mass == 0): # comp with CO core and the star_base has a He core
                        merged_star.center_h1 = comp.center_h1
                        merged_star.center_he4 = comp.center_he4
                        merged_star.center_c12 = comp.center_c12
                        merged_star.center_n14 = comp.center_n14
                        merged_star.center_o16 = comp.center_o16
                    elif (star_base.co_core_mass > 0 and comp.co_core_mass > 0):
                        merged_star.center_h1 = mass_weighted_avg(mass_weight_name="co_core_mass")
                        merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight_name="co_core_mass")
                        merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight_name="co_core_mass")
                        merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight_name="co_core_mass")
                        merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight_name="co_core_mass")
                    else:
                        warnings.warn("weird compbination of CO core masses during merging")

                    for key in STARPROPERTIES:
                        # these stellar attributes become np.nan
                        for substring in ["log_", "lg_", "surface","surf_", "conv_", "lambda", "profile"]:
                            if (substring in key) :
                                setattr(merged_star, key, np.nan)
                            if key in [ "c12_c12", "center_gamma",
                                       "avg_c_in_c_core", "total_moment_of_inertia", "spin"]:
                                setattr(merged_star, key, np.nan)
                    merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!

            elif (s1 is in STAR_STATES_HE_RICH:
                and s2 is in ["WD"]):

                    # add total and core masses
                    for key in ["mass", "he_core_mass", "c_core_mass", "o_core_mass", "co_core_mass"]:
                        current = getattr(merged_star, key) + getattr(comp, "mass")
                        setattr(merged_star, key,current)

                    # weighted central abundances if merging cores. Else only from star_base
                    if (comp.co_core_mass > 0 and star_base.co_core_mass == 0): # comp with CO core and the star_base has a He core
                        merged_star.center_h1 = comp.center_h1
                        merged_star.center_he4 = comp.center_he4
                        merged_star.center_c12 = comp.center_c12
                        merged_star.center_n14 = comp.center_n14
                        merged_star.center_o16 = comp.center_o16
                    elif (star_base.co_core_mass > 0 and comp.co_core_mass > 0):
                        merged_star.center_h1 = mass_weighted_avg(mass_weight_name="co_core_mass")
                        merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight_name="co_core_mass")
                        merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight_name="co_core_mass")
                        merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight_name="co_core_mass")
                        merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight_name="co_core_mass")
                    else:
                        warnings.warn("weird compbination of CO core masses during merging")

                    for key in STARPROPERTIES:
                        # these stellar attributes become np.nan
                        for substring in ["log_", "lg_", "surface","surf_", "conv_", "lambda", "profile"]:
                            if (substring in key) :
                                setattr(merged_star, key, np.nan)
                            if key in [ "c12_c12", "center_gamma",
                                       "avg_c_in_c_core", "total_moment_of_inertia", "spin"]:
                                setattr(merged_star, key, np.nan)
                    merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!

            elif (s1 is in STAR_STATES_HE_RICH:
                and s2 is in ["NS", "BH"]):
                    merged_star = comp

            else:
                print("Combination of merging star states not expected: ", s1, s2)

        # ad hoc spin of merged star to be used in the detached step
        merged_star.surf_avg_omega_div_omega_crit = merger_critical_rot

        return merged_star, None


    def __call__(self,binary):

        if binary.state == "merged":
            if binary.event == 'oMerging1':
                binary.star_1,binary.star_2 = merged_star_properties(binary.star_1,binary.star_2)
            elif binary.event == 'oMerging2':
                binary.star_2,binary.star_1 = merged_star_properties(binary.star_2,binary.star_1)
            else:
                raise ValueError("binary.state='merged' but binary.event != 'oMerging1/2'")
        else:
            raise ValueError("step_merging initated but binary.state != 'merged'")

        binary.event == None

        super().__call__(binary)
