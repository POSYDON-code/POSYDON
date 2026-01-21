"""Merging and isolated evolution step."""


__authors__ = [
    "Emmanouil Zapartas <ezapartas@gmail.com>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Max Briel <max.briel@gmail.com>",
]


import numpy as np

from posydon.binary_evol.DT.step_isolated import IsolatedStep
from posydon.binary_evol.flow_chart import (
    STAR_STATES_H_RICH,
    STAR_STATES_HE_RICH,
    STAR_STATES_NOT_CO,
)
from posydon.binary_evol.singlestar import (
    STARPROPERTIES,
    convert_star_to_massless_remnant,
)
from posydon.config import PATH_TO_POSYDON_DATA
from posydon.utils.common_functions import check_state_of_star
from posydon.utils.constants import Zsun, zams_table
from posydon.utils.posydonerror import ModelError
from posydon.utils.posydonwarning import Pwarn

LIST_ACCEPTABLE_STATES_FOR_HMS = ["H-rich_Core_H_burning"]
LIST_ACCEPTABLE_STATES_FOR_HeMS = ["stripped_He_Core_He_burning"]

LIST_ACCEPTABLE_STATES_FOR_POSTMS = STAR_STATES_H_RICH.copy()
[LIST_ACCEPTABLE_STATES_FOR_POSTMS.remove(x) for x in LIST_ACCEPTABLE_STATES_FOR_HMS]

LIST_ACCEPTABLE_STATES_FOR_POSTHeMS = STAR_STATES_HE_RICH.copy()
[LIST_ACCEPTABLE_STATES_FOR_POSTHeMS.remove(x) for x in LIST_ACCEPTABLE_STATES_FOR_HeMS]


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
        rel_mass_lost_HMS_HMS = 0.1, # [0-1) or "Glebbeek+2013" (dependent on mass ratio q; Glebbeek E., Gaburov E., Portegies Zwart S., Pols O. R., 2013, MNRAS, 434, 3497)
        HMS_HMS_merging_rejuvenation = True, # if True then new abundances based on total_mass_h1 or he4 
                                             # ("Schneider+2016"; i.e. Schneider, F. R. N., Podsiadlowski, P., Langer, N., Castro, N., & Fossati, L. 2016, MNRAS,457, 2355
                                             # https://ui.adsabs.harvard.edu/abs/2016MNRAS.457.2355S/abstract)
        list_for_matching_HMS = [
                #["mass", "center_h1", "he_core_mass"],
                #[20.0, 1.0, 10.0],
                ["mass", "total_mass_h1", "he_core_mass"],
                [20.0, 20.0, 10.0],
                ["log_min_max", "min_max", "min_max"],
                #[m_min_H, m_max_H], [0, None]
                [None, None], [0, None]
            ],
        list_for_matching_postMS = [
                ["mass", "center_he4", "he_core_mass"],
                [20.0, 1.0, 10.0],
                ["log_min_max", "min_max", "min_max"],
                #[m_min_H, m_max_H], [0, None]
                [None, None], [0, None]
            ],
        list_for_matching_HeStar = [
                ["he_core_mass", "center_he4"],
                [10.0, 1.0],
                ["min_max" , "min_max"],
                #[[m_min_He, m_max_He], [0, None]],
                [None, None], [0, None]
            ],
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

        if self.verbose:
            print("Before Merger", binary.star_1.state, binary.star_2.state,
                  binary.state, binary.event)
            print("M1 , M2, he_core_mass1, he_core_mass2, co_core_mass1, co_core_mass2: ",
                  binary.star_1.mass, binary.star_2.mass,
                  binary.star_1.he_core_mass, binary.star_2.he_core_mass, binary.star_1.co_core_mass, binary.star_2.co_core_mass)
            print("star_1.center_h1, star_2.center_h1, star_1.center_he4, star_1.center_c12, star_2.center_c12, "
                  "star_2.center_he4, star_1.surface_he4, star_2.surface_he4: ",
                  binary.star_1.center_h1, binary.star_2.center_h1, binary.star_1.center_he4,
                  binary.star_2.center_he4, binary.star_1.center_c12, binary.star_2.center_c12, binary.star_1.surface_he4,
                  binary.star_2.surface_he4)

        if binary.state == "merged":
            if binary.event == 'oMerging1':
                binary.star_1, binary.star_2 = merged_star_properties(binary.star_1, binary.star_2)
            elif binary.event == 'oMerging2':
                binary.star_2, binary.star_1 = merged_star_properties(binary.star_2, binary.star_1)
            else:
                raise ValueError("binary.state='merged' but binary.event != 'oMerging1/2'")

        ## assume that binaries in RLO with two He-rich stars always merge
        elif binary.star_1.state in STAR_STATES_HE_RICH and binary.star_2.state in STAR_STATES_HE_RICH:
            binary.state = "merged"
            if binary.event == 'oRLO1':
                binary.star_1, binary.star_2 = merged_star_properties(binary.star_1, binary.star_2)
            elif binary.event == 'oRLO2':
                binary.star_2, binary.star_1 = merged_star_properties(binary.star_2, binary.star_1)
            else:
                raise ValueError("step_merged initiated for He stars but RLO not initiated")
        else:
            raise ValueError("step_merged initiated but binary is not in valid merging state!")

        binary.event = None

        if self.verbose:
            print("After Merger", binary.star_1.state, binary.star_2.state,
                  binary.state, binary.event)
            if binary.event == 'oMerging1':
                print("M_merged , he_core_mass merged, co_core_mass merged: ", binary.star_1.mass,
                  binary.star_1.he_core_mass, binary.star_1.co_core_mass)
                print("star_1.center_h1, star_1.center_he4, star_1.center_c12, star_1.surface_he4: ",
                  binary.star_1.center_h1, binary.star_1.center_he4, binary.star_1.center_c12, binary.star_1.surface_he4)
            elif binary.event=='oMerging2':
                print("M_merged , he_core_mass merged, co_core_mass merged: ", binary.star_2.mass,
                  binary.star_2.he_core_mass, binary.star_2.co_core_mass)
                print("star_2.center_h1, star_2.center_he4, star_2.center_c12, star_2.surface_he4: ",
                  binary.star_2.center_h1, binary.star_2.center_he4, binary.star_2.center_c12, binary.star_2.surface_he4)

        super().__call__(binary)

    def merged_star_properties(self, star_base, comp):
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
            if mass_weight1 == "H-rich_envelope_mass":
                M1 = getattr(star1, "mass") - getattr(star1, "he_core_mass")
            elif mass_weight1 == "He-rich_envelope_mass":
                M1 = getattr(star1, "he_core_mass") - getattr(star1, "co_core_mass")
            else:
                M1 = getattr(star1, mass_weight1)

            if mass_weight2 is None:
                mass_weight2 = mass_weight1
            if mass_weight2 == "H-rich_envelope_mass":
                M2 = getattr(star2, "mass") - getattr(star2, "he_core_mass")
            elif mass_weight2 == "He-rich_envelope_mass":
                M2 = getattr(star2, "he_core_mass") - getattr(star2, "co_core_mass")
            else:
                M2 = getattr(star2, mass_weight2)

            try:
            	mass_weighted_avg_value=(A1*M1+A2*M2)/(M1+M2)
            except TypeError:
            	mass_weighted_avg_value= np.nan

            if self.verbose:
                print(abundance_name, mass_weight1, mass_weight2)
                print("A_base, M_base_abundance, A_comp, M_comp_abundanc", A1, M1, A2, M2)
                print("mass_weighted_avg= ", mass_weighted_avg_value)

            return mass_weighted_avg_value

        # MS + MS
        if ( s1 in LIST_ACCEPTABLE_STATES_FOR_HMS
            and s2 in LIST_ACCEPTABLE_STATES_FOR_HMS):

                q =  comp.mass/star_base.mass
                # first we calculate the mass loss percentage:
                if self.rel_mass_lost_HMS_HMS == "Glebbeek+2013":
                    # Eq. 4  of Glebbeek+2013
                    rel_mass_lost_HMS_HMS = (0.3*q)/((1.+q)**2.)
                else:
                    rel_mass_lost_HMS_HMS = self.rel_mass_lost_HMS_HMS


                if HMS_HMS_merging_rejuvenation :
                    X_average1 = star_base.total_mass_h1 / star_base.mass
                    X_average2 = comp.total_mass_h1 / comp.mass

                    Zinitial_div_Zsun = star_base.metallicity
                    if Zinitial_div_Zsun in zams_table.keys():
                        Y_initial = zams_table[Zinitial_div_Zsun]
                    else:
                        raise KeyError(f"{Zinitial_div_Zsun} is a not defined metallicity")
                    Z_initial = Zinitial_div_Zsun*Zsun
                    X_initial = 1.0 - Y_initial - Z_initial

                    merged_star.mass = (star_base.mass + comp.mass) * (1.-rel_mass_lost_HMS_HMS)

                    X_average_merged = (star_base.total_mass_h1 + comp.total_mass_h1 - (1.-rel_mass_lost_HMS_HMS) * X_initial) / merged_star.mass
                    Y_average_merged = (star_base.total_mass_he4 + comp.total_mass_he4 - (1.-rel_mass_lost_HMS_HMS) * Y_initial) / merged_star.mass

                    merged_star.total_mass_h1 = X_average_merged * merged_star.mass
                    merged_star.center_h1 = X_average_merged # this is only consistent with the merged_star.total_mass_h1 above if we assume full mixing during merger
                    merged_star.surface_h1 = X_average_merged # this is only consistent with the merged_star.total_mass_h1 above if we assume full mixing during merger
                    #in any case, the track_match after the step_merged will decide which parameter will be used for matching (center_h1, center_he4 or total_mass_h1)

                    merged_star.total_mass_he4 = Y_average_merged * merged_star.mass
                    merged_star.center_he4 = Y_average_merged # this is only consistent with the merged_star.total_mass_he4 above if we assume full mixing during merger
                    merged_star.surface_he4 = Y_average_merged # this is only consistent with the merged_star.total_mass_he4 above if we assume full mixing during merger
                    #in any case, the track_match after the step_merged will decide which parameter will be used for matching (center_h1, center_he4 or total_mass_h1)

                    for key in STARPROPERTIES:
                    # these stellar attributes become np.nan as
                        if key in [ "center_c12", "center_n14", "center_o16"]:
                            setattr(merged_star, key, np.nan)

                else:
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

                    # The change of masses occurs after the calculation of weighted averages
                    merged_star.mass = (star_base.mass + comp.mass) * (1.-rel_mass_lost_HMS_HMS)

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

                # weigheted mixing on the surface abundances of the whole comp with the envelope of star_base
                merged_star.surface_h1 = mass_weighted_avg(abundance_name = "surface_h1", mass_weight1="H-rich_envelope_mass", mass_weight2="mass")
                merged_star.surface_he4 = mass_weighted_avg(abundance_name = "surface_he4", mass_weight1="H-rich_envelope_mass", mass_weight2="mass")
                merged_star.surface_c12 = mass_weighted_avg(abundance_name = "surface_c12", mass_weight1="H-rich_envelope_mass", mass_weight2="mass")
                merged_star.surface_n14 = mass_weighted_avg(abundance_name = "surface_n14", mass_weight1="H-rich_envelope_mass", mass_weight2="mass")
                merged_star.surface_o16 = mass_weighted_avg(abundance_name = "surface_o16", mass_weight1="H-rich_envelope_mass", mass_weight2="mass")

                # The change of masses occurs after the calculation of weighted averages
                merged_star.mass = star_base.mass + comp.mass

                for key in STARPROPERTIES:
                    # these stellar attributes become np.nan
                    for substring in ["log_", "lg_", "surf_", "conv_", "lambda", "profile"]:
                        if (substring in key) :
                            setattr(merged_star, key, np.nan)
                        if key in [ "c12_c12", "total_moment_of_inertia", "spin"]:
                            setattr(merged_star, key, np.nan)

                merged_star.log_LHe = star_base.log_LHe
                merged_star.log_LZ = star_base.log_LZ

                merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!
                massless_remnant = convert_star_to_massless_remnant(comp)

        # as above but opposite stars
        elif (s1 in LIST_ACCEPTABLE_STATES_FOR_HMS
            and s2 in LIST_ACCEPTABLE_STATES_FOR_POSTMS):

                merged_star = comp

                merged_star.surface_h1 = mass_weighted_avg(abundance_name = "surface_h1", mass_weight2="H-rich_envelope_mass", mass_weight1="mass")
                merged_star.surface_he4 = mass_weighted_avg(abundance_name = "surface_he4", mass_weight2="H-rich_envelope_mass", mass_weight1="mass")
                merged_star.surface_c12 = mass_weighted_avg(abundance_name = "surface_c12", mass_weight2="H-rich_envelope_mass", mass_weight1="mass")
                merged_star.surface_n14 = mass_weighted_avg(abundance_name = "surface_n14", mass_weight2="H-rich_envelope_mass", mass_weight1="mass")
                merged_star.surface_o16 = mass_weighted_avg(abundance_name = "surface_o16", mass_weight2="H-rich_envelope_mass", mass_weight1="mass")

                # The change of masses occurs after the calculation of weighted averages
                merged_star.mass = star_base.mass + comp.mass

                for key in STARPROPERTIES:
                    # these stellar attributes become np.nan
                    for substring in ["log_", "lg_", "surf_", "conv_", "lambda", "profile"]:
                        if (substring in key) :
                            setattr(merged_star, key, np.nan)
                        if key in [ "c12_c12", "total_moment_of_inertia", "spin"]:
                            setattr(merged_star, key, np.nan)

                merged_star.log_LHe = comp.log_LHe
                merged_star.log_LZ = comp.log_LZ

                merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!
                massless_remnant = convert_star_to_massless_remnant(star_base)

        #postMS + postMS
        elif (s1 in LIST_ACCEPTABLE_STATES_FOR_POSTMS
            and s2 in LIST_ACCEPTABLE_STATES_FOR_POSTMS):

                # weighted central abundances if merging cores. Else only from star_base
                if star_base.co_core_mass == 0 and comp.co_core_mass == 0: # two stars with Helium cores
                    merged_star.center_h1 = mass_weighted_avg(mass_weight1="he_core_mass")
                    merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight1="he_core_mass")
                    merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight1="he_core_mass")
                    merged_star.avg_c_in_c_core = mass_weighted_avg(abundance_name = "avg_c_in_c_core", mass_weight1="he_core_mass")
                    merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight1="he_core_mass")
                    merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight1="he_core_mass")
                    setattr(merged_star, "center_gamma", np.nan)
                elif (star_base.co_core_mass > 0 and comp.co_core_mass == 0): # star_base with CO core and the comp has a He core
                    pass # the central abundances are kept as the ones of star_base
                elif (comp.co_core_mass > 0 and star_base.co_core_mass == 0): # comp with CO core and the star_base has a He core
                    merged_star.center_h1 = comp.center_h1
                    merged_star.center_he4 = comp.center_he4
                    merged_star.center_c12 = comp.center_c12
                    merged_star.avg_c_in_c_core = comp.avg_c_in_c_core
                    merged_star.center_n14 = comp.center_n14
                    merged_star.center_o16 = comp.center_o16
                    merged_star.center_gamma =  comp.center_gamma
                elif (star_base.co_core_mass > 0 and comp.co_core_mass > 0):
                    merged_star.center_h1 = mass_weighted_avg(mass_weight1="co_core_mass") #TODO : maybe he_core_mass makes more sense?
                    merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight1="co_core_mass")
                    merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight1="co_core_mass")
                    merged_star.avg_c_in_c_core = mass_weighted_avg(abundance_name = "avg_c_in_c_core", mass_weight1="co_core_mass")
                    merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight1="co_core_mass")
                    merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight1="co_core_mass")
                    setattr(merged_star, "center_gamma", np.nan)
                else:
                    Pwarn("weird compbination of CO core masses during merging", "EvolutionWarning")

                # weigheted mixing on the surface abundances based on the envelopes of the two stars
                merged_star.surface_h1 = mass_weighted_avg(abundance_name = "surface_h1", mass_weight1="H-rich_envelope_mass", mass_weight2="H-rich_envelope_mass")
                merged_star.surface_he4 = mass_weighted_avg(abundance_name = "surface_he4", mass_weight1="H-rich_envelope_mass", mass_weight2="H-rich_envelope_mass")
                merged_star.surface_c12 = mass_weighted_avg(abundance_name = "surface_c12", mass_weight1="H-rich_envelope_mass", mass_weight2="H-rich_envelope_mass")
                merged_star.surface_n14 = mass_weighted_avg(abundance_name = "surface_n14", mass_weight1="H-rich_envelope_mass", mass_weight2="H-rich_envelope_mass")
                merged_star.surface_o16 = mass_weighted_avg(abundance_name = "surface_o16", mass_weight1="H-rich_envelope_mass", mass_weight2="H-rich_envelope_mass")

                # add total and core masses after calculations of weighter average
                for key in ["mass", "he_core_mass", "c_core_mass", "o_core_mass", "co_core_mass"]:
                    current = getattr(star_base, key) + getattr(comp, key)
                    setattr(merged_star, key,current)

                for key in STARPROPERTIES:
                    # these stellar attributes become np.nan
                    for substring in ["log_", "lg_", "surf_", "conv_", "lambda", "profile"]:
                        if (substring in key) :
                            setattr(merged_star, key, np.nan)
                        if key in [ "c12_c12", "total_moment_of_inertia", "spin"]:
                            setattr(merged_star, key, np.nan)
                merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!
                massless_remnant = convert_star_to_massless_remnant(comp)

        #postMS + HeMSStar
        elif (s1 in LIST_ACCEPTABLE_STATES_FOR_POSTMS
            and s2 in LIST_ACCEPTABLE_STATES_FOR_HeMS):

                # weighted central abundances if merging cores. Else only from star_base
                if star_base.co_core_mass == 0 and comp.co_core_mass == 0: # two stars with only Helium cores, not CO cores
                    merged_star.center_h1 = mass_weighted_avg(mass_weight1="he_core_mass")
                    merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight1="he_core_mass")
                    merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight1="he_core_mass")
                    merged_star.avg_c_in_c_core = mass_weighted_avg(abundance_name = "avg_c_in_c_core", mass_weight1="he_core_mass")
                    merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight1="he_core_mass")
                    merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight1="he_core_mass")
                    setattr(merged_star, "center_gamma", np.nan)
                elif (star_base.co_core_mass > 0 and comp.co_core_mass == 0): # star_base with CO core and the comp has just a He core (is a HeMS star)
                    pass # the central abundances are kept as the ones of star_base
                else:
                    Pwarn("weird compbination of CO core masses during merging", "EvolutionWarning")
                # surface abundances from star_base

                # add total and core masses after abundance mass weighted calculations
                for key in ["mass", "he_core_mass", "c_core_mass", "o_core_mass", "co_core_mass"]:
                    current = getattr(star_base, key) + getattr(comp, key)
                    setattr(merged_star, key,current)

                for key in STARPROPERTIES:
                    # these stellar attributes become np.nan
                    for substring in ["log_", "lg_", "surf_", "conv_", "lambda", "profile"]:
                        if (substring in key) :
                            setattr(merged_star, key, np.nan)
                        if key in [ "c12_c12", "total_moment_of_inertia", "spin"]:
                            setattr(merged_star, key, np.nan)
                merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!
                massless_remnant = convert_star_to_massless_remnant(comp)

        # as above but opposite stars
        elif (s1 in  LIST_ACCEPTABLE_STATES_FOR_HeMS
            and s2 in LIST_ACCEPTABLE_STATES_FOR_POSTMS):

                merged_star = comp

                # weighted central abundances if merging cores. Else only from comp
                if star_base.co_core_mass == 0 and comp.co_core_mass == 0: # two stars with only Helium cores, not CO cores
                    merged_star.center_h1 = mass_weighted_avg(mass_weight1="he_core_mass")
                    merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight1="he_core_mass")
                    merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight1="he_core_mass")
                    merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight1="he_core_mass")
                    merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight1="he_core_mass")
                elif (star_base.co_core_mass == 0 and comp.co_core_mass > 0): # star_base is the HeMS Star and comp has a CO core
                    pass # the central abundances are kept as the ones of comp
                else:
                    Pwarn("weird compbination of CO core masses during merging", "EvolutionWarning")
                # surface abundances from comp

                # add total and core masses after weighted averages above
                for key in ["mass", "he_core_mass", "c_core_mass", "o_core_mass", "co_core_mass"]:
                    current = getattr(star_base, key) + getattr(comp, key)
                    setattr(merged_star, key,current)

                for key in STARPROPERTIES:
                    # these stellar attributes become np.nan
                    for substring in ["log_", "lg_", "surf_", "conv_", "lambda", "profile"]:
                        if (substring in key) :
                            setattr(merged_star, key, np.nan)
                        if key in [ "c12_c12", "total_moment_of_inertia", "spin"]:
                            setattr(merged_star, key, np.nan)

                merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!
                massless_remnant = convert_star_to_massless_remnant(star_base)

        #postMS + HeStar that is not in HeMS
        elif (s1 in LIST_ACCEPTABLE_STATES_FOR_POSTMS
            and s2 in LIST_ACCEPTABLE_STATES_FOR_POSTHeMS):

                # weighted central abundances if merging cores. Else only from star_base
                if star_base.co_core_mass == 0 and comp.co_core_mass == 0: # two stars with Helium cores
                    merged_star.center_h1 = mass_weighted_avg(mass_weight1="he_core_mass")
                    merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight1="he_core_mass")
                    merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight1="he_core_mass")
                    merged_star.avg_c_in_c_core = mass_weighted_avg(abundance_name = "avg_c_in_c_core", mass_weight1="he_core_mass")
                    merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight1="he_core_mass")
                    merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight1="he_core_mass")
                    setattr(merged_star, "center_gamma", np.nan)
                elif (star_base.co_core_mass > 0 and comp.co_core_mass == 0): # star_base with CO core and the comp has a He core
                    pass # the central abundances are kept as the ones of star_base
                elif (comp.co_core_mass > 0 and star_base.co_core_mass == 0): # comp with CO core and the star_base has a He core
                    merged_star.center_h1 = comp.center_h1
                    merged_star.center_he4 = comp.center_he4
                    merged_star.center_c12 = comp.center_c12
                    merged_star.avg_c_in_c_core = comp.avg_c_in_c_core
                    merged_star.center_n14 = comp.center_n14
                    merged_star.center_o16 = comp.center_o16
                    merged_star.center_gamma =  comp.center_gamma
                elif (star_base.co_core_mass > 0 and comp.co_core_mass > 0):
                    merged_star.center_h1 = mass_weighted_avg(mass_weight1="co_core_mass")
                    merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight1="co_core_mass")
                    merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight1="co_core_mass")
                    merged_star.avg_c_in_c_core = mass_weighted_avg(abundance_name = "avg_c_in_c_core", mass_weight1="co_core_mass")
                    merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight1="co_core_mass")
                    merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight1="co_core_mass")
                    setattr(merged_star, "center_gamma", np.nan)
                else:
                    Pwarn("weird compbination of CO core masses during merging", "EvolutionWarning")

                # add total and core masses
                for key in ["mass", "he_core_mass", "c_core_mass", "o_core_mass", "co_core_mass"]:
                    current = getattr(star_base, key) + getattr(comp, key)
                    setattr(merged_star, key,current)

                for key in STARPROPERTIES:
                    # these stellar attributes become np.nan
                    for substring in ["log_", "lg_", "surf_", "conv_", "lambda", "profile"]:
                        if (substring in key) :
                            setattr(merged_star, key, np.nan)
                        if key in [ "c12_c12", "total_moment_of_inertia", "spin"]:
                            setattr(merged_star, key, np.nan)
                merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!
                massless_remnant = convert_star_to_massless_remnant(comp)

        # as above but the opposite stars
        elif (s1 in LIST_ACCEPTABLE_STATES_FOR_POSTHeMS
            and s2 in LIST_ACCEPTABLE_STATES_FOR_POSTMS):

                merged_star = comp

                # weighted central abundances if merging cores. Else only from star_base
                if star_base.co_core_mass == 0 and comp.co_core_mass == 0: # two stars with Helium cores
                    merged_star.center_h1 = mass_weighted_avg(mass_weight1="he_core_mass")
                    merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight1="he_core_mass")
                    merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight1="he_core_mass")
                    merged_star.avg_c_in_c_core = mass_weighted_avg(abundance_name = "avg_c_in_c_core", mass_weight1="he_core_mass")
                    merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight1="he_core_mass")
                    merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight1="he_core_mass")
                    setattr(merged_star, "center_gamma", np.nan)
                elif (star_base.co_core_mass > 0 and comp.co_core_mass == 0): # star_base with CO core and the comp has a He core
                    merged_star.center_h1 = star_base.center_h1
                    merged_star.center_he4 = star_base.center_he4
                    merged_star.center_c12 = star_base.center_c12
                    merged_star.center_n14 = star_base.center_n14
                    merged_star.center_o16 = star_base.center_o16
                elif (comp.co_core_mass > 0 and star_base.co_core_mass == 0): # comp with CO core and the star_base has a He core
                    pass # the central abundances are kept as the ones of comp
                elif (star_base.co_core_mass > 0 and comp.co_core_mass > 0):
                    merged_star.center_h1 = mass_weighted_avg(mass_weight1="co_core_mass")
                    merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight1="co_core_mass")
                    merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight1="co_core_mass")
                    merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight1="co_core_mass")
                    merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight1="co_core_mass")
                else:
                    Pwarn("weird compbination of CO core masses during merging", "EvolutionWarning")

                # add total and core masses
                for key in ["mass", "he_core_mass", "c_core_mass", "o_core_mass", "co_core_mass"]:
                    current = getattr(star_base, key) + getattr(comp, key)
                    setattr(merged_star, key,current)

                for key in STARPROPERTIES:
                    # these stellar attributes become np.nan
                    for substring in ["log_", "lg_", "surf_", "conv_", "lambda", "profile"]:
                        if (substring in key) :
                            setattr(merged_star, key, np.nan)
                        if key in [ "c12_c12", "total_moment_of_inertia", "spin"]:
                            setattr(merged_star, key, np.nan)
                merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!
                massless_remnant = convert_star_to_massless_remnant(star_base)

        # HeStar + HeStar
        elif (s1 in STAR_STATES_HE_RICH
            and s2 in STAR_STATES_HE_RICH):

                # weighted central abundances if merging cores. Else only from star_base
                if star_base.co_core_mass == 0 and comp.co_core_mass == 0: # two stars with Helium cores
                    merged_star.center_h1 = mass_weighted_avg(mass_weight1="he_core_mass")
                    merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight1="he_core_mass")
                    merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight1="he_core_mass")
                    merged_star.avg_c_in_c_core = mass_weighted_avg(abundance_name = "avg_c_in_c_core", mass_weight1="he_core_mass")
                    merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight1="he_core_mass")
                    merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight1="he_core_mass")
                    setattr(merged_star, "center_gamma", np.nan)
                elif (star_base.co_core_mass > 0 and comp.co_core_mass == 0): # star_base with CO core and the comp has a He core
                    pass # the central abundances are kept as the ones of star_base
                elif (comp.co_core_mass > 0 and star_base.co_core_mass == 0): # comp with CO core and the star_base has a He core
                    merged_star.center_h1 = comp.center_h1
                    merged_star.center_he4 = comp.center_he4
                    merged_star.center_c12 = comp.center_c12
                    merged_star.avg_c_in_c_core = comp.avg_c_in_c_core
                    merged_star.center_n14 = comp.center_n14
                    merged_star.center_o16 = comp.center_o16
                    merged_star.center_gamma =  comp.center_gamma
                elif (star_base.co_core_mass > 0 and comp.co_core_mass > 0):
                    merged_star.center_h1 = mass_weighted_avg(mass_weight1="co_core_mass")
                    merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight1="co_core_mass")
                    merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight1="co_core_mass")
                    merged_star.avg_c_in_c_core = mass_weighted_avg(abundance_name = "avg_c_in_c_core", mass_weight1="co_core_mass")
                    merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight1="co_core_mass")
                    merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight1="co_core_mass")
                    setattr(merged_star, "center_gamma", np.nan)
                else:
                    Pwarn("weird compbination of CO core masses during merging", "EvolutionWarning")

                # add total and core masses
                for key in ["mass", "he_core_mass", "c_core_mass", "o_core_mass", "co_core_mass"]:
                    current = getattr(star_base, key) + getattr(comp, key)
                    setattr(merged_star, key,current)

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
                        if key in [ "c12_c12", "total_moment_of_inertia", "spin"]:
                            setattr(merged_star, key, np.nan)
                merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!
                massless_remnant = convert_star_to_massless_remnant(comp)

        # Star + WD
        elif (s1 in STAR_STATES_NOT_CO
            and s2 in ["WD"]):

                # WD is considered a stripped CO core
                # WD would always be the comp, it cannot be the engulfing star (so no need to do the opposite stars case below)

                # weighted central abundances if merging cores. Else only from star_base
                if (comp.co_core_mass is not None and star_base.co_core_mass == 0): # comp with CO core and the star_base has not
                    merged_star.center_h1 = comp.center_h1
                    merged_star.center_he4 = comp.center_he4
                    merged_star.center_c12 = comp.center_c12
                    merged_star.avg_c_in_c_core = comp.avg_c_in_c_core
                    merged_star.center_n14 = comp.center_n14
                    merged_star.center_o16 = comp.center_o16
                    merged_star.center_gamma =  comp.center_gamma
                elif (comp.co_core_mass is not None and star_base.co_core_mass > 0):
                    merged_star.center_h1 = mass_weighted_avg(mass_weight1="co_core_mass")
                    merged_star.center_he4 = mass_weighted_avg(abundance_name = "center_he4", mass_weight1="co_core_mass")
                    merged_star.center_c12 = mass_weighted_avg(abundance_name = "center_c12", mass_weight1="co_core_mass")
                    merged_star.avg_c_in_c_core = mass_weighted_avg(abundance_name = "avg_c_in_c_core", mass_weight1="co_core_mass")
                    merged_star.center_n14 = mass_weighted_avg(abundance_name = "center_n14", mass_weight1="co_core_mass")
                    merged_star.center_o16 = mass_weighted_avg(abundance_name = "center_o16", mass_weight1="co_core_mass")
                    setattr(merged_star, "center_gamma", np.nan)
                else:
                    Pwarn("weird compbination of CO core masses during merging", "EvolutionWarning")

                # add total and core masses
                for key in ["mass", "he_core_mass", "c_core_mass", "o_core_mass", "co_core_mass"]:
                    current = getattr(star_base, key) + getattr(comp, "mass")
                    setattr(merged_star, key,current)

                for key in STARPROPERTIES:
                    # these stellar attributes become np.nan
                    for substring in ["log_", "lg_", "surf_", "conv_", "lambda", "profile"]:
                        if (substring in key) :
                            setattr(merged_star, key, np.nan)
                        if key in [ "c12_c12", "total_moment_of_inertia", "spin"]:
                            setattr(merged_star, key, np.nan)
                merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!
                massless_remnant = convert_star_to_massless_remnant(comp)

        # Star + NS/BH
        elif (s1 in STAR_STATES_NOT_CO
            and s2 in ["NS", "BH"]):
                # TODO: potentially flag a Thorne-Zytkov object
                massless_remnant = convert_star_to_massless_remnant(star_base)

                ## in this case, want CO companion object to stay the same, and base star to be assigned massless remnant
                return massless_remnant, comp
        else:
            raise ModelError(f"Combination of merging star states not expected: {s1} {s2}")
        # ad hoc spin of merged star to be used in the detached step
        merged_star.surf_avg_omega_div_omega_crit = self.merger_critical_rot

        return merged_star, massless_remnant
