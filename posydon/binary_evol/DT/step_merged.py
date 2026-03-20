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
        rel_mass_lost_HMS_HMS = 0.1,
        list_for_matching_HMS = [
                ["mass", "center_h1", "he_core_mass"],
                [20.0, 1.0, 10.0],
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
            print("Before Merger:\n"
                  f"_____________\n")
            print(f"M1 = {binary.star_1.mass}\n"
                  f"M2 = {binary.star_2.mass}\n"
                  f"he_core_mass1 = {binary.star_1.he_core_mass}\n"
                  f"he_core_mass2 = {binary.star_2.he_core_mass}\n"
                  f"co_core_mass1 = {binary.star_1.co_core_mass}\n"
                  f"co_core_mass2 = {binary.star_2.co_core_mass}\n")
            print(f"star_1.center_h1 = {binary.star_1.center_h1}\n"
                  f"star_2.center_h1 = {binary.star_2.center_h1}\n"
                  f"star_1.center_he4 = {binary.star_1.center_he4}\n"
                  f"star_2.center_he4 = {binary.star_2.center_he4}\n"
                  f"star_1.center_c12 = {binary.star_1.center_c12}\n"
                  f"star_2.center_c12 = {binary.star_2.center_c12}\n"
                  f"star_1.surface_he4 = {binary.star_1.surface_he4}\n"
                  f"star_2.surface_he4 = {binary.star_2.surface_he4}\n"
                  f"_____________\n")

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
            print("After Merger:\n"
                  f"_____________\n"
                  f"star_1.state = {binary.star_1.state}\n"
                  f"star_2.state = {binary.star_2.state}\n"
                  f"binary.state = {binary.state}\n"
                  f"binary.event = {binary.event}\n")
            if binary.event == 'oMerging1':
                print(f"M_merged = {binary.star_1.mass}\n"
                      f"he_core_mass merged = {binary.star_1.he_core_mass}\n"
                      f"co_core_mass merged = {binary.star_1.co_core_mass}\n")
                print(f"star_1.center_h1 = {binary.star_1.center_h1}\n"
                      f"star_1.center_he4 = {binary.star_1.center_he4}\n"
                      f"star_1.center_c12 = {binary.star_1.center_c12}\n"
                      f"star_1.surface_he4 = {binary.star_1.surface_he4}\n")
            elif binary.event=='oMerging2':
                print(f"M_merged = {binary.star_2.mass}\n"
                      f"he_core_mass merged = {binary.star_2.he_core_mass}\n"
                      f"co_core_mass merged = {binary.star_2.co_core_mass}\n")
                print(f"star_2.center_h1 = {binary.star_2.center_h1}\n"
                      f"star_2.center_he4 = {binary.star_2.center_he4}\n"
                      f"star_2.center_c12 = {binary.star_2.center_c12}\n"
                      f"star_2.surface_he4 = {binary.star_2.surface_he4}\n")

        super().__call__(binary)

    def mass_weighted_avg(self, star1, star2, abundance_name, mass_weight1="mass", mass_weight2='mass'):
        """Compute the mass-weighted average of an abundance between two stars.

        Parameters
        ----------
        star1 : SingleStar
            Primary star
        star2 : SingleStar
            Companion star
        abundance_name : str
            Name of the SingleStar attribute to average.
        mass_weight1 : str
            Mass attribute to use as weight for star1.  Special values
            ``"H-rich_envelope_mass"`` and ``"He-rich_envelope_mass"`` are
            computed on the fly; any other value is looked up directly on the
            star object.
        mass_weight2 : str or None
            Mass attribute for star2.
        """
        # give NaN if the abundance is not defined for one of the stars
        A1 = getattr(star1, abundance_name, np.nan)
        A2 = getattr(star2, abundance_name, np.nan)
        if A1 is None:
            A1 = np.nan
        if A2 is None:
            A2 = np.nan

        if mass_weight1 == "mass":
            M1 = getattr(star1, "mass", np.nan)
        elif mass_weight1 == "H-rich_envelope_mass":
            M1 = getattr(star1, "mass", np.nan) - getattr(star1, "he_core_mass", np.nan)
        elif mass_weight1 == "He-rich_envelope_mass":
            M1 = getattr(star1, "he_core_mass", np.nan) - getattr(star1, "co_core_mass", np.nan)
        else:
            M1 = getattr(star1, mass_weight1, np.nan)

        if M1 is None:
            M1 = np.nan

        if mass_weight2 == "mass":
            M2 = getattr(star2, "mass", np.nan)
        elif mass_weight2 == "H-rich_envelope_mass":
            M2 = getattr(star2, "mass", np.nan) - getattr(star2, "he_core_mass", np.nan)
        elif mass_weight2 == "He-rich_envelope_mass":
            M2 = getattr(star2, "he_core_mass", np.nan) - getattr(star2, "co_core_mass", np.nan)
        else:
            M2 = getattr(star2, mass_weight2, np.nan)

        if M2 is None:
            M2 = np.nan

        den = M1 + M2
        if not (np.isfinite(A1) and np.isfinite(A2) and np.isfinite(M1) and np.isfinite(M2)):
            mass_weighted_avg_value = np.nan
        elif den == 0:
            mass_weighted_avg_value = np.nan
        else:
            mass_weighted_avg_value = (A1 * M1 + A2 * M2) / den

        if self.verbose:
            print(f"Mass weighted average for {abundance_name}:\n"
                  f"weight 1: {mass_weight1}\n"
                  f"weight 2: {mass_weight2}\n"
                  f"A_base = {A1}\n"
                  f"M_base_abundance = {M1}\n"
                  f"A_comp = {A2}\n"
                  f"M_comp_abundance = {M2}\n"
                  f"mass_weighted_avg = {mass_weighted_avg_value}\n")

        return mass_weighted_avg_value

    def _apply_nan_attributes(self, star, extra_nan_keys=None):
        """Set to np.nan the attributes of the star that are not expected to be meaningful after a merger event.

        Parameters set to np.nan are:
        - all attributes containing the substrings
            {"log_", "lg_", "surf_", "conv_", "lambda", "profile"}
        - all attributes in the set:
            {"c12_c12", "total_moment_of_inertia", "spin",}

        Parameters
        ----------
        star: SingleStar
            The star for which the attributes will be set to np.nan
        extra_nan_keys: set of str
            Set of keys to be set to np.nan in addition to the default ones.
        """
        if extra_nan_keys is None:
            extra_nan_keys = set()

        _NAN_SUBSTRINGS = ("log_", "lg_", "surf_", "conv_", "lambda", "profile")
        _NAN_KEYS = set(["c12_c12", "total_moment_of_inertia", "spin"])

        for key in STARPROPERTIES:
            if any(sub in key for sub in _NAN_SUBSTRINGS) or (key in extra_nan_keys) or (key in _NAN_KEYS):
                setattr(star, key, np.nan)

    def star_type(self, state):
        if state in LIST_ACCEPTABLE_STATES_FOR_HMS:
            return "HMS"
        elif state in LIST_ACCEPTABLE_STATES_FOR_POSTMS:
            return "postMS"
        elif state in LIST_ACCEPTABLE_STATES_FOR_HeMS:
            return "HeMS"
        elif state in LIST_ACCEPTABLE_STATES_FOR_POSTHeMS:
            return "postHeMS"
        elif state in ["WD", "NS", "BH"]:
            return state
        
    def mix_cores(self, merged_star, star_base, comp):

        parameters_to_mix = ["center_h1", "center_he4", "center_c12",  
                             "center_n14", "center_o16", "avg_c_in_c_core"]

        # both have CO cores
        # can trigger: postMS-postMS, postMS-postHeMS
        if (star_base.co_core_mass > 0 and comp.co_core_mass > 0):

            #TODO : maybe he_core_mass makes more sense?
            mass_weight1='co_core_mass'
            mass_weight2='co_core_mass'
            setattr(merged_star, "center_gamma", np.nan)
        
        # both have Helium cores
        # can trigger: postMS-postMS, postMS-HeMS, HeMS-postMS, postMS-postHeMS
        elif (star_base.co_core_mass == 0 and comp.co_core_mass == 0) & \
             (star_base.he_core_mass > 0 and comp.he_core_mass > 0):
            
            mass_weight1 = 'he_core_mass'
            mass_weight2 = 'he_core_mass'
            setattr(merged_star, "center_gamma", np.nan)

        # star_base He core and companion CO core
        # can trigger: postMS-postMS, HeMS-postMS, postMS-postHeMS
        elif (star_base.co_core_mass == 0 and comp.co_core_mass > 0):

            for abundance_name in parameters_to_mix:
                setattr(merged_star, abundance_name, getattr(comp, abundance_name))

            setattr(merged_star, "center_gamma", comp.center_gamma)
            return
        
        # star_base CO core and comp He core
        # can trigger: postMS-postMS, postMS-HeMS, postMS-postHeMS
        elif (star_base.co_core_mass > 0 and comp.co_core_mass == 0) & \
             (comp.he_core_mass > 0):
            # the central abundances are kept as the ones of star_base
            return
        # both have H cores
        elif not (star_base.he_core_mass and comp.he_core_mass):
            mass_weight1 = "mass"
            mass_weight2 = "mass" 

        else:
            Pwarn("Unexpected combination of CO core masses during merging", "EvolutionWarning")
            return

        for abundance_name in parameters_to_mix:
            setattr(merged_star, abundance_name, self.mass_weighted_avg(star_base,
                                                                        comp,
                                                                        abundance_name = abundance_name,
                                                                        mass_weight1 = mass_weight1,
                                                                        mass_weight2 = mass_weight2))

    def mix_envelope(self, merged_star, star_base, comp):
        
        mass_weight = [None, None]
        for i, star in enumerate([star_base, comp]):
            if self.star_type(star) == "HeMS":
                mass_weight[i] = "He-rich_envelope_mass"
            elif self.star_type(star) == "HMS":
                mass_weight[i] = "mass"
            # post HMS
            else:
                mass_weight[i] = "H-rich_envelope_mass"

        parameters_to_mix = ["surface_h1", "surface_he4", "surface_c12", 
                             "surface_n14", "surface_o16"]

        # Weighted mixing on the surface abundances based on the envelopes of the two stars
        for abundance_name in parameters_to_mix:
            setattr(merged_star, abundance_name, self.mass_weighted_avg(star_base,
                                                                        comp,
                                                                        abundance_name=abundance_name,
                                                                        mass_weight1=mass_weight[0],
                                                                        mass_weight2=mass_weight[1]))
        
    def set_post_merger(self, merged_star, star_base, comp):

        if merged_star.merger_type == "HMS-HMS":
            # The change of masses occurs after the calculation of weighted averages
            merged_star.mass = (star_base.mass + comp.mass) * (1.-self.rel_mass_lost_HMS_HMS)

            # Set parameters that are not expected to be meaningful after a merger to np.nan
            extra_nan_keys = set(["center_gamma", "avg_c_in_c_core", "envelope_binding_energy"])

        elif merged_star.merger_type == "postMS-HMS" or \
             merged_star.merger_type == "HMS-postMS":
            # The change of masses occurs after the calculation of weighted averages
            # Note that the he core mass is unchanged, which means that
            # adding the mass changes the envelope mass only.
            merged_star.mass = star_base.mass + comp.mass

            # Set parameters that are not expected to be meaningful after a merger to np.nan
            extra_nan_keys = None

            # Set burning luminosities and star state.
            merged_star.log_LHe = star_base.log_LHe
            merged_star.log_LZ = star_base.log_LZ
        
        else:
            # add total and core masses
            for key in ["mass", "he_core_mass", "c_core_mass", "o_core_mass", "co_core_mass"]:
                current = getattr(star_base, key) + getattr(comp, "mass")
                setattr(merged_star, key,current)

            # Set parameters that are not expected to be meaningful after a merger to np.nan
            extra_nan_keys = None

        self._apply_nan_attributes(merged_star, extra_nan_keys=extra_nan_keys)
        merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!
        massless_remnant = convert_star_to_massless_remnant(comp)

        return massless_remnant
        
    def merged_star_properties(self, star_base, comp):

        """Set the properties of the merged star after a merger event.

        Generally, the Roche lobe overflowing star becomes star_base, except
        when the companion is further evolved than star_base, in which
        case the comp becomes the base for the merged star.

        Abundances are mass-weighted averages of the two stars, with the
        weights depending on the type of merger and the abundance considered.

        star_base: Single Star
            The star that engulfs its companion.
            (generally the base for the merged_star)
        comp: Single Star
            The star that is engulfed by star_base.
        """

        # By default the stellar attributes are kept from star_base.
        merged_star = star_base
        state1 = star_base.state
        state2 = comp.state
        merged_star.merger_type = "-".join([self.star_type(state1),
                                            self.star_type(state2)])

        # MS + MS
        if (merged_star.merger_type == "HMS-HMS"):
            # Mix central and surface abundances of the two stars with
            # mass-weighted averages based on the total masses of the two stars.
            self.mix_envelope(merged_star, star_base, comp)
            self.mix_cores(merged_star, star_base, comp)
            massless_remnant = self.set_post_merger(merged_star, star_base, comp)

        # postMS + MS
        elif (merged_star.merger_type == "postMS-HMS"):
            # weighted mixing on the surface abundances of the whole
            # companion with the envelope of star_base
            self.mix_envelope(merged_star, star_base, comp)
            massless_remnant = self.set_post_merger(merged_star, star_base, comp)

        # MS + postMS (the opposite of above)
        elif (merged_star.merger_type == "HMS-postMS"):
            merged_star = comp
            self.mix_envelope(merged_star, star_base, comp)
            massless_remnant = self.set_post_merger(merged_star, comp, star_base)

        # postMS + postMS
        elif (merged_star.merger_type == "postMS-postMS"):
            # Weighted central abundances if merging cores.
            self.mix_cores(merged_star, star_base, comp)
            self.mix_envelope(merged_star, star_base, comp)
            massless_remnant = self.set_post_merger(merged_star, star_base, comp)

        # postMS + HeMS
        elif (merged_star.merger_type == "postMS-HeMS"):
            # Keep surface abundances from star_base.
            self.mix_cores(merged_star, star_base, comp)
            massless_remnant = self.set_post_merger(merged_star, star_base, comp)

        # HeMS + postMS (the opposite of above)
        elif (merged_star.merger_type == "HeMS-postMS"):
            # surface abundances from comp
            merged_star = comp
            self.mix_cores(merged_star, star_base, comp)
            massless_remnant = self.set_post_merger(merged_star, comp, star_base)

        # postMS + postHeMS
        elif (merged_star.merger_type == "postMS-postHeMS"):
            self.mix_cores(merged_star, star_base, comp)
            massless_remnant = self.set_post_merger(merged_star, star_base, comp)

        # postHeMS + postMS (the opposite of above)
        elif (merged_star.merger_type == "postHeMS-postMS"):
            merged_star = comp
            self.mix_cores(merged_star, star_base, comp)
            massless_remnant = self.set_post_merger(merged_star, comp, star_base)

        # He star + He star
        elif (merged_star.merger_type == "HeMS-HeMS"):
            # weigheted mixing on the surface abundances based on the He-rich envelopes of the two stars
            self.mix_envelope(merged_star, star_base, comp)
            self.mix_cores(merged_star, star_base, comp)
            massless_remnant = self.set_post_merger(merged_star, star_base, comp)

        # Star + WD
        elif ("-WD" in merged_star.merger_type):
            self.mix_cores(merged_star, star_base, comp)
            massless_remnant = self.set_post_merger(merged_star, star_base, comp)

        # Star + NS/BH
        elif ("-NS" in merged_star.merger_type or "-BH" in merged_star.merger_type):
            # TODO: potentially flag a Thorne-Zytkov object
            massless_remnant = convert_star_to_massless_remnant(star_base)
            ## in this case, want CO companion object to stay the same, and base star to be assigned massless remnant
            return massless_remnant, comp
        
        else:
            raise ModelError(f"Combination of merging star states not expected:\n"
                             f"Star 1: {state1}\n"
                             f"Star 2: {state2}\n")
        
        # ad hoc spin of merged star to be used in the detached step
        merged_star.surf_avg_omega_div_omega_crit = self.merger_critical_rot
        return merged_star, massless_remnant
    