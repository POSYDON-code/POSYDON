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
from posydon.binary_evol.DT.step_detached import detached_step
from posydon.binary_evol.DT.step_isolated import IsolatedStep

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


def convert_star_to_massless_remnant(star):
    star.state="massless_remnant"
    for key in STARPROPERTIES:
        if key is not "state":
            setattr(star, key, np.nan)
    return star


class MergedStep(IsolatedStep):
    """
    Prepare a merging star to do an an IsolatedStep
    """

    def __init__(
        self,
        grid_name_Hrich=None,
        grid_name_strippedHe=None,
        path=PATH_TO_POSYDON_DATA,
        merger_critical_rot = 0.4,
        rel_mass_lost_HMS_HMS = 0.1,
        list_for_matching_HMS = [["mass", "center_h1", "log_R", "he_core_mass"],
                              [20.0, 1.0, 2.0, 10.0],
                                 # [[m_min_H, m_max_H], [0, None]],
                              ["log_min_max" , "min_max", "min_max", "min_max"] ],
        list_for_matching_postMS = [["mass", "center_he4", "he_core_mass"],
                                 [20.0, 1.0, 10.0],
                                 #[[m_min_H, m_max_H], [0, None]],
                                 ["log_min_max" , "min_max",  "min_max"] ],
        list_for_matching_HeStar = [["he_core_mass", "center_he4"],
                                 [10.0, 1.0],
                                 #[[m_min_He, m_max_He], [0, None]],
                                 ["min_max" , "min_max"]  ],
        *args,
        **kwargs
    ):

        self.merger_critical_rot = merger_critical_rot
        self.rel_mass_lost_HMS_HMS = rel_mass_lost_HMS_HMS

        super().__init__(
        grid_name_Hrich=grid_name_Hrich,
        grid_name_strippedHe=grid_name_strippedHe,
        list_for_matching_HMS = list_for_matching_HMS,
        list_for_matching_postMS = list_for_matching_postMS,
        list_for_matching_HeStar = list_for_matching_HeStar,
        *args,
        **kwargs)


    def __call__(self,binary):
        merged_star_properties = self.merged_star_properties
        print("call step_merged")
        print("BEFORE", binary.star_1.state,binary.star_2.state,binary.state, binary.event)
        print(binary.star_1.mass,binary.star_2.mass,binary.star_1.center_he4,binary.star_2.center_he4, binary.star_1.surface_he4,binary.star_2.surface_he4)
        if binary.state == "merged":
            if binary.event == 'oMerging1':
                binary.star_1,binary.star_2 = merged_star_properties(binary.star_1,binary.star_2)
            elif binary.event == 'oMerging2':
                binary.star_2,binary.star_1 = merged_star_properties(binary.star_2,binary.star_1)
            else:
                raise ValueError("binary.state='merged' but binary.event != 'oMerging1/2'")
        else:
            raise ValueError("step_merging initiated but binary.state != 'merged'")

        binary.event = None
        print("AFTER",binary.star_1.state,binary.star_2.state, binary.state, binary.event)
        print(binary.star_1.mass,binary.star_2.mass,binary.star_1.center_he4,binary.star_2.center_he4, binary.star_1.surface_he4,binary.star_2.surface_he4)

        super().__call__(binary)

    def merged_star_properties(self,star_base,comp):
        """
        Make assumptions about the core/total mass, and abundances of the star of a merged product.

        Similar to the table of merging in BSE

        star_base: Single Star
            is our base star that engulfs its companion. The merged star will have this star as a base
        comp: Single Star
            is the star that is engulfed
        """
        #by default the stellar attributes that keep the same value from the
        #merged_star = copy.copy(star_base)
        merged_star = star_base

        s1 = star_base.state
        s2 = comp.state

        def mass_weighted_avg(star1=star_base,star2=comp, abundance_name="center_h1", mass_weight1="mass", mass_weight2=None):
            A1 = getattr(star1, abundance_name)
            A2 = getattr(star2, abundance_name)
            #print(2.1, A1, A2)

            #print(2.15, mass_weight1, mass_weight2)
            if mass_weight1 == "H-rich_envelope_mass":
                M1 = getattr(star1, "mass") - getattr(star1, "he_core_mass")
                #print(2.16, getattr(star1, "mass"), getattr(star1, "he_core_mass"), M1)
            elif mass_weight1 == "He-rich_envelope_mass":
                M1 = getattr(star1, "he_core_mass") - getattr(star1, "co_core_mass")
            else:
                M1 = getattr(star1, mass_weight1)

            if mass_weight2 is None:
                mass_weight2 = mass_weight1
            #print(2.17, mass_weight2)
            if mass_weight2 == "H-rich_envelope_mass":
                M2 = getattr(star2, "mass") - getattr(star2, "he_core_mass")
            elif mass_weight2 == "He-rich_envelope_mass":
                M2 = getattr(star2, "he_core_mass") - getattr(star2, "co_core_mass")
            else:
                M2 = getattr(star2, mass_weight2)
            #print(2.2, M1, M2)
            print(abundance_name, A1,A2,M1, M2)
            print("avg", (A1*M1 + A2*M2 ) / (M1+M2))
            return (A1*M1 + A2*M2 ) / (M1+M2)

        # MS + MS
        if ( s1 in LIST_ACCEPTABLE_STATES_FOR_HMS
            and s2 in LIST_ACCEPTABLE_STATES_FOR_HMS):
                #these stellar attributes change value
                merged_star.mass = (star_base.mass + comp.mass) * (1.-self.rel_mass_lost_HMS_HMS)

                #TODO for key in ["center_h1", "center_he4", "center_c12", "center_n14","center_o16"]:
                merged_star.center_h1 = mass_weighted_avg()
                merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4")
                merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12")
                merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14")
                merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16")
                #TODO: should I check if the abundaces above end up in ~1 (?)

                # weigheted mixing on the surface abundances
                merged_star.surface_h1 = mass_weighted_avg(abundance_name = "surface_h1")
                merged_star.surface_he4 = mass_weighted_avg(abundance_name = "surface_he4")
                merged_star.surface_c12 = mass_weighted_avg(abundance_name = "surface_c12")
                merged_star.surface_n14 = mass_weighted_avg(abundance_name = "surface_n14")
                merged_star.surface_o16 = mass_weighted_avg(abundance_name = "surface_o16")

                for key in STARPROPERTIES:
                    # these stellar attributes become np.nan
                    for substring in ["log_", "lg_", "surf_", "conv_", "lambda", "profile"]:
                        if (substring in key) :
                            setattr(merged_star, key, np.nan)
                        if key in [ "c12_c12", "center_gamma",
                                   "avg_c_in_c_core", "total_moment_of_inertia", "spin", "envelope_binding_energy"]:
                            setattr(merged_star, key, np.nan)
                massless_remnant = convert_star_to_massless_remnant(comp)

        #postMS + MS
        elif (s1 in LIST_ACCEPTABLE_STATES_FOR_POSTMS
            and s2 in LIST_ACCEPTABLE_STATES_FOR_HMS):

                merged_star.mass = star_base.mass + comp.mass #TODO: in step_CEE we need to eject part of the (common) envelope
                #print(1.1)
                # weigheted mixing on the surface abundances of the whole comp with the envelope of star_base
                merged_star.surface_h1 = mass_weighted_avg(abundance_name = "surface_h1", mass_weight1="H-rich_envelope_mass", mass_weight2="mass")
                merged_star.surface_he4 = mass_weighted_avg(abundance_name = "surface_he4", mass_weight1="H-rich_envelope_mass", mass_weight2="mass")
                merged_star.surface_c12 = mass_weighted_avg(abundance_name = "surface_c12", mass_weight1="H-rich_envelope_mass", mass_weight2="mass")
                merged_star.surface_n14 = mass_weighted_avg(abundance_name = "surface_n14", mass_weight1="H-rich_envelope_mass", mass_weight2="mass")
                merged_star.surface_o16 = mass_weighted_avg(abundance_name = "surface_o16", mass_weight1="H-rich_envelope_mass", mass_weight2="mass")

                for key in STARPROPERTIES:
                    # these stellar attributes become np.nan
                    for substring in ["log_", "lg_", "surf_", "conv_", "lambda", "profile"]:
                        if (substring in key) :
                            setattr(merged_star, key, np.nan)
                        if key in [ "c12_c12", "center_gamma",
                                   "avg_c_in_c_core", "total_moment_of_inertia", "spin"]:
                            setattr(merged_star, key, np.nan)
                #print(1.2)
                merged_star.log_LHe = star_base.log_LHe
                merged_star.log_LZ = star_base.log_LZ

                merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!
                massless_remnant = convert_star_to_massless_remnant(comp)
                #print(1.3)
        # as above but opposite stars
        elif (s1 in LIST_ACCEPTABLE_STATES_FOR_HMS
            and s2 in LIST_ACCEPTABLE_STATES_FOR_POSTMS):

                merged_star = comp
                merged_star.mass = star_base.mass + comp.mass

                merged_star.surface_h1 = mass_weighted_avg(abundance_name = "surface_h1", mass_weight2="H-rich_envelope_mass", mass_weight1="mass")
                merged_star.surface_he4 = mass_weighted_avg(abundance_name = "surface_he4", mass_weight2="H-rich_envelope_mass", mass_weight1="mass")
                merged_star.surface_c12 = mass_weighted_avg(abundance_name = "surface_c12", mass_weight2="H-rich_envelope_mass", mass_weight1="mass")
                merged_star.surface_n14 = mass_weighted_avg(abundance_name = "surface_n14", mass_weight2="H-rich_envelope_mass", mass_weight1="mass")
                merged_star.surface_o16 = mass_weighted_avg(abundance_name = "surface_o16", mass_weight2="H-rich_envelope_mass", mass_weight1="mass")

                for key in STARPROPERTIES:
                    # these stellar attributes become np.nan
                    for substring in ["log_", "lg_", "surf_", "conv_", "lambda", "profile"]:
                        if (substring in key) :
                            setattr(merged_star, key, np.nan)
                        if key in [ "c12_c12", "center_gamma",
                                   "avg_c_in_c_core", "total_moment_of_inertia", "spin"]:
                            setattr(merged_star, key, np.nan)

                merged_star.log_LHe = comp.log_LHe
                merged_star.log_LZ = comp.log_LZ

                merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!
                massless_remnant = convert_star_to_massless_remnant(star_base)

        #postMS + postMS
        elif (s1 in LIST_ACCEPTABLE_STATES_FOR_POSTMS
            and s2 in LIST_ACCEPTABLE_STATES_FOR_POSTMS):

                # add total and core masses
                for key in ["mass", "he_core_mass", "c_core_mass", "o_core_mass", "co_core_mass"]:
                    current = getattr(merged_star, key) + getattr(comp, key)
                    setattr(merged_star, key,current)

                # weighted central abundances if merging cores. Else only from star_base
                if star_base.co_core_mass == 0 and comp.co_core_mass == 0: # two stars with Helium cores
                    merged_star.center_h1 = mass_weighted_avg(mass_weight1="he_core_mass")
                    merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight1="he_core_mass")
                    merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight1="he_core_mass")
                    merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight1="he_core_mass")
                    merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight1="he_core_mass")
                elif (star_base.co_core_mass > 0 and comp.co_core_mass == 0): # star_base with CO core and the comp has a He core
                    pass # the central abundances are kept as the ones of star_base
                elif (comp.co_core_mass > 0 and star_base.co_core_mass == 0): # comp with CO core and the star_base has a He core
                    merged_star.center_h1 = comp.center_h1
                    merged_star.center_he4 = comp.center_he4
                    merged_star.center_c12 = comp.center_c12
                    merged_star.center_n14 = comp.center_n14
                    merged_star.center_o16 = comp.center_o16
                elif (star_base.co_core_mass > 0 and comp.co_core_mass > 0):
                    merged_star.center_h1 = mass_weighted_avg(mass_weight1="co_core_mass") #TODO : maybe he_core_mass makes more sense?
                    merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight1="co_core_mass")
                    merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight1="co_core_mass")
                    merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight1="co_core_mass")
                    merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight1="co_core_mass")
                else:
                    warnings.warn("weird compbination of CO core masses during merging")

                # weigheted mixing on the surface abundances based on the envelopes of the two stars
                merged_star.surface_h1 = mass_weighted_avg(abundance_name = "surface_h1", mass_weight1="H-rich_envelope_mass", mass_weight2="H-rich_envelope_mass")
                merged_star.surface_he4 = mass_weighted_avg(abundance_name = "surface_he4", mass_weight1="H-rich_envelope_mass", mass_weight2="H-rich_envelope_mass")
                merged_star.surface_c12 = mass_weighted_avg(abundance_name = "surface_c12", mass_weight1="H-rich_envelope_mass", mass_weight2="H-rich_envelope_mass")
                merged_star.surface_n14 = mass_weighted_avg(abundance_name = "surface_n14", mass_weight1="H-rich_envelope_mass", mass_weight2="H-rich_envelope_mass")
                merged_star.surface_o16 = mass_weighted_avg(abundance_name = "surface_o16", mass_weight1="H-rich_envelope_mass", mass_weight2="H-rich_envelope_mass")


                for key in STARPROPERTIES:
                    # these stellar attributes become np.nan
                    for substring in ["log_", "lg_", "surf_", "conv_", "lambda", "profile"]:
                        if (substring in key) :
                            setattr(merged_star, key, np.nan)
                        if key in [ "c12_c12", "center_gamma",
                                   "avg_c_in_c_core", "total_moment_of_inertia", "spin"]:
                            setattr(merged_star, key, np.nan)
                merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!
                massless_remnant = convert_star_to_massless_remnant(comp)

        #postMS + HeMSStar
        elif (s1 in LIST_ACCEPTABLE_STATES_FOR_POSTMS
            and s2 in LIST_ACCEPTABLE_STATES_FOR_HeMS):

                # add total and core masses
                for key in ["mass", "he_core_mass", "c_core_mass", "o_core_mass", "co_core_mass"]:
                    current = getattr(merged_star, key) + getattr(comp, key)
                    setattr(merged_star, key,current)

                # weighted central abundances if merging cores. Else only from star_base
                if star_base.co_core_mass == 0 and comp.co_core_mass == 0: # two stars with only Helium cores, not CO cores
                    merged_star.center_h1 = mass_weighted_avg(mass_weight1="he_core_mass")
                    merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight1="he_core_mass")
                    merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight1="he_core_mass")
                    merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight1="he_core_mass")
                    merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight1="he_core_mass")
                elif (star_base.co_core_mass > 0 and comp.co_core_mass == 0): # star_base with CO core and the comp has just a He core (is a HeMS star)
                    pass # the central abundances are kept as the ones of star_base
                else:
                    warnings.warn("weird compbination of CO core masses during merging")

                for key in STARPROPERTIES:
                    # these stellar attributes become np.nan
                    for substring in ["log_", "lg_", "surf_", "conv_", "lambda", "profile"]:
                        if (substring in key) :
                            setattr(merged_star, key, np.nan)
                        if key in [ "c12_c12", "center_gamma",
                                   "avg_c_in_c_core", "total_moment_of_inertia", "spin"]:
                            setattr(merged_star, key, np.nan)
                merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!
                massless_remnant = convert_star_to_massless_remnant(comp)

        # as above but opposite stars
        elif (s1 in  LIST_ACCEPTABLE_STATES_FOR_HeMS
            and s2 in LIST_ACCEPTABLE_STATES_FOR_POSTMS):

                merged_star = comp

                # add total and core masses
                for key in ["mass", "he_core_mass", "c_core_mass", "o_core_mass", "co_core_mass"]:
                    current = getattr(merged_star, key) + getattr(star_base, key)
                    setattr(merged_star, key,current)

                # weighted central abundances if merging cores. Else only from star_base
                if star_base.co_core_mass == 0 and comp.co_core_mass == 0: # two stars with only Helium cores, not CO cores
                    merged_star.center_h1 = mass_weighted_avg(mass_weight1="he_core_mass")
                    merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight1="he_core_mass")
                    merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight1="he_core_mass")
                    merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight1="he_core_mass")
                    merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight1="he_core_mass")
                elif (star_base.co_core_mass == 0 and comp.co_core_mass > 0): # star_base is the HeMS Star and comp has a CO core
                    pass # the central abundances are kept as the ones of star_base
                else:
                    warnings.warn("weird compbination of CO core masses during merging")

                for key in STARPROPERTIES:
                    # these stellar attributes become np.nan
                    for substring in ["log_", "lg_", "surf_", "conv_", "lambda", "profile"]:
                        if (substring in key) :
                            setattr(merged_star, key, np.nan)
                        if key in [ "c12_c12", "center_gamma",
                                   "avg_c_in_c_core", "total_moment_of_inertia", "spin"]:
                            setattr(merged_star, key, np.nan)
                merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!
                massless_remnant = convert_star_to_massless_remnant(star_base)

        #postMS + HeStar that is not in HeMS
        elif (s1 in LIST_ACCEPTABLE_STATES_FOR_POSTMS
            and s2 in LIST_ACCEPTABLE_STATES_FOR_POSTHeMS):

                # add total and core masses
                for key in ["mass", "he_core_mass", "c_core_mass", "o_core_mass", "co_core_mass"]:
                    current = getattr(merged_star, key) + getattr(comp, key)
                    setattr(merged_star, key,current)

                # weighted central abundances if merging cores. Else only from star_base
                if star_base.co_core_mass == 0 and comp.co_core_mass == 0: # two stars with Helium cores
                    merged_star.center_h1 = mass_weighted_avg(mass_weight1="he_core_mass")
                    merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight1="he_core_mass")
                    merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight1="he_core_mass")
                    merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight1="he_core_mass")
                    merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight1="he_core_mass")
                elif (star_base.co_core_mass > 0 and comp.co_core_mass == 0): # star_base with CO core and the comp has a He core
                    pass # the central abundances are kept as the ones of star_base
                elif (comp.co_core_mass > 0 and star_base.co_core_mass == 0): # comp with CO core and the star_base has a He core
                    merged_star.center_h1 = comp.center_h1
                    merged_star.center_he4 = comp.center_he4
                    merged_star.center_c12 = comp.center_c12
                    merged_star.center_n14 = comp.center_n14
                    merged_star.center_o16 = comp.center_o16
                elif (star_base.co_core_mass > 0 and comp.co_core_mass > 0):
                    merged_star.center_h1 = mass_weighted_avg(mass_weight1="co_core_mass")
                    merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight1="co_core_mass")
                    merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight1="co_core_mass")
                    merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight1="co_core_mass")
                    merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight1="co_core_mass")
                else:
                    warnings.warn("weird compbination of CO core masses during merging")

                for key in STARPROPERTIES:
                    # these stellar attributes become np.nan
                    for substring in ["log_", "lg_", "surf_", "conv_", "lambda", "profile"]:
                        if (substring in key) :
                            setattr(merged_star, key, np.nan)
                        if key in [ "c12_c12", "center_gamma",
                                   "avg_c_in_c_core", "total_moment_of_inertia", "spin"]:
                            setattr(merged_star, key, np.nan)
                merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!
                massless_remnant = convert_star_to_massless_remnant(comp)

        # as above but the opposite stars
        elif (s1 in LIST_ACCEPTABLE_STATES_FOR_POSTHeMS
            and s2 in LIST_ACCEPTABLE_STATES_FOR_POSTMS):

                merged_star = comp
                # add total and core masses
                for key in ["mass", "he_core_mass", "c_core_mass", "o_core_mass", "co_core_mass"]:
                    current = getattr(merged_star, key) + getattr(star_base, key)
                    setattr(merged_star, key,current)

                # weighted central abundances if merging cores. Else only from star_base
                if star_base.co_core_mass == 0 and comp.co_core_mass == 0: # two stars with Helium cores
                    merged_star.center_h1 = mass_weighted_avg(mass_weight1="he_core_mass")
                    merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight1="he_core_mass")
                    merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight1="he_core_mass")
                    merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight1="he_core_mass")
                    merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight1="he_core_mass")
                elif (star_base.co_core_mass > 0 and comp.co_core_mass == 0): # star_base with CO core and the comp has a He core
                    pass # the central abundances are kept as the ones of star_base
                elif (comp.co_core_mass > 0 and star_base.co_core_mass == 0): # comp with CO core and the star_base has a He core
                    merged_star.center_h1 = comp.center_h1
                    merged_star.center_he4 = comp.center_he4
                    merged_star.center_c12 = comp.center_c12
                    merged_star.center_n14 = comp.center_n14
                    merged_star.center_o16 = comp.center_o16
                elif (star_base.co_core_mass > 0 and comp.co_core_mass > 0):
                    merged_star.center_h1 = mass_weighted_avg(mass_weight1="co_core_mass")
                    merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight1="co_core_mass")
                    merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight1="co_core_mass")
                    merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight1="co_core_mass")
                    merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight1="co_core_mass")
                else:
                    warnings.warn("weird compbination of CO core masses during merging")

                for key in STARPROPERTIES:
                    # these stellar attributes become np.nan
                    for substring in ["log_", "lg_", "surf_", "conv_", "lambda", "profile"]:
                        if (substring in key) :
                            setattr(merged_star, key, np.nan)
                        if key in [ "c12_c12", "center_gamma",
                                   "avg_c_in_c_core", "total_moment_of_inertia", "spin"]:
                            setattr(merged_star, key, np.nan)
                merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!
                massless_remnant = convert_star_to_massless_remnant(star_base)

        # HeStar + HeStar
        elif (s1 in STAR_STATES_HE_RICH
            and s2 in STAR_STATES_HE_RICH):

                # add total and core masses
                for key in ["mass", "he_core_mass", "c_core_mass", "o_core_mass", "co_core_mass"]:
                    current = getattr(merged_star, key) + getattr(comp, key)
                    setattr(merged_star, key,current)

                # weighted central abundances if merging cores. Else only from star_base
                if star_base.co_core_mass == 0 and comp.co_core_mass == 0: # two stars with Helium cores
                    merged_star.center_h1 = mass_weighted_avg(mass_weight1="he_core_mass")
                    merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight1="he_core_mass")
                    merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight1="he_core_mass")
                    merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight1="he_core_mass")
                    merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight1="he_core_mass")
                elif (star_base.co_core_mass > 0 and comp.co_core_mass == 0): # star_base with CO core and the comp has a He core
                    pass # the central abundances are kept as the ones of star_base
                elif (comp.co_core_mass > 0 and star_base.co_core_mass == 0): # comp with CO core and the star_base has a He core
                    merged_star.center_h1 = comp.center_h1
                    merged_star.center_he4 = comp.center_he4
                    merged_star.center_c12 = comp.center_c12
                    merged_star.center_n14 = comp.center_n14
                    merged_star.center_o16 = comp.center_o16
                elif (star_base.co_core_mass > 0 and comp.co_core_mass > 0):
                    merged_star.center_h1 = mass_weighted_avg(mass_weight1="co_core_mass")
                    merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight1="co_core_mass")
                    merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight1="co_core_mass")
                    merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight1="co_core_mass")
                    merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight1="co_core_mass")
                else:
                    warnings.warn("weird compbination of CO core masses during merging")

                # weigheted mixing on the surface abundances based on the He-rich envelopes of the two stars
                merged_star.surface_h1 = mass_weighted_avg(abundance_name = "surface_h1", mass_weight1="He-rich_envelope_mass", mass_weight2="He-rich_envelope_mass")
                merged_star.surface_he4 = mass_weighted_avg(abundance_name = "surface_he4", mass_weight1="He-rich_envelope_mass", mass_weight2="He-rich_envelope_mass")
                merged_star.surface_c12 = mass_weighted_avg(abundance_name = "surface_c12", mass_weight1="He-rich_envelope_mass", mass_weight2="He-rich_envelope_mass")
                merged_star.surface_n14 = mass_weighted_avg(abundance_name = "surface_n14", mass_weight1="He-rich_envelope_mass", mass_weight2="He-rich_envelope_mass")
                merged_star.surface_o16 = mass_weighted_avg(abundance_name = "surface_o16", mass_weight1="He-rich_envelope_mass", mass_weight2="He-rich_envelope_mass")

                for key in STARPROPERTIES:
                    # these stellar attributes become np.nan
                    for substring in ["log_", "lg_", "surf_", "conv_", "lambda", "profile"]:
                        if (substring in key) :
                            setattr(merged_star, key, np.nan)
                        if key in [ "c12_c12", "center_gamma",
                                   "avg_c_in_c_core", "total_moment_of_inertia", "spin"]:
                            setattr(merged_star, key, np.nan)
                merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!
                massless_remnant = convert_star_to_massless_remnant(comp)

        # Star + WD
        elif (s1 in STAR_STATES_NOT_CO
            and s2 in ["WD"]):

                #WD is considered a stripped CO core

                # add total and core masses
                for key in ["mass", "he_core_mass", "c_core_mass", "o_core_mass", "co_core_mass"]:
                    current = getattr(merged_star, key) + getattr(comp, "mass")
                    setattr(merged_star, key,current)

                # weighted central abundances if merging cores. Else only from star_base
                if (comp.co_core_mass > 0 and star_base.co_core_mass == 0): # comp with CO core and the star_base has not
                    merged_star.center_h1 = comp.center_h1
                    merged_star.center_he4 = comp.center_he4
                    merged_star.center_c12 = comp.center_c12
                    merged_star.center_n14 = comp.center_n14
                    merged_star.center_o16 = comp.center_o16
                elif (star_base.co_core_mass > 0 and comp.co_core_mass > 0):
                    merged_star.center_h1 = mass_weighted_avg(mass_weight1="co_core_mass")
                    merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight1="co_core_mass")
                    merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight1="co_core_mass")
                    merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight1="co_core_mass")
                    merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight1="co_core_mass")
                else:
                    warnings.warn("weird compbination of CO core masses during merging")

                for key in STARPROPERTIES:
                    # these stellar attributes become np.nan
                    for substring in ["log_", "lg_", "surf_", "conv_", "lambda", "profile"]:
                        if (substring in key) :
                            setattr(merged_star, key, np.nan)
                        if key in [ "c12_c12", "center_gamma",
                                   "avg_c_in_c_core", "total_moment_of_inertia", "spin"]:
                            setattr(merged_star, key, np.nan)
                merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!
                massless_remnant = convert_star_to_massless_remnant(comp)

        # Star + NS/BH
        elif (s1 in STAR_STATES_NOT_CO
            and s2 in ["NS", "BH"]):
                merged_star = comp
                # TODO: potentially flag a Thorne-Zytkov object
                massless_remnant = convert_star_to_massless_remnant(star_base)
        else:
            print("Combination of merging star states not expected: ", s1, s2)

        # ad hoc spin of merged star to be used in the detached step
        merged_star.surf_avg_omega_div_omega_crit = self.merger_critical_rot

        #print("3", merged_star.state,massless_remnant.state)
        return merged_star, massless_remnant
