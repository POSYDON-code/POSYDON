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

    def __init__(self,
                 merger_critical_rot = 0.4,
                 rel_mass_lost_HMS_HMS = 0.1,
                 *args,
                 **kwargs):

        self.merger_critical_rot = merger_critical_rot
        self.rel_mass_lost_HMS_HMS = rel_mass_lost_HMS_HMS

        super().__init__(*args, **kwargs)

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

        #binary.event = None

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

        s1 = star_base.state
        s2 = comp.state

        # MS + MS
        if ( s1 in LIST_ACCEPTABLE_STATES_FOR_HMS
            and s2 in LIST_ACCEPTABLE_STATES_FOR_HMS):
                # Mix central and surface abundances of the two stars with
                # mass-weighted averages based on the total masses of the two stars.
                parameters_to_mix = ["center_h1", "center_he4",  "center_c12",  "center_n14",  "center_o16",
                                    "surface_h1", "surface_he4", "surface_c12", "surface_n14", "surface_o16"]

                for abundance_name in parameters_to_mix:
                    #TODO: should I check if the abundaces above end up in ~1 (?)
                    setattr(merged_star, abundance_name, self.mass_weighted_avg(star_base,
                                                                                comp,
                                                                                abundance_name=abundance_name,
                                                                                mass_weight1="mass",
                                                                                mass_weight2='mass'))

                # The change of masses occurs after the calculation of weighted averages
                merged_star.mass = (star_base.mass + comp.mass) * (1.-self.rel_mass_lost_HMS_HMS)

                # Set parameters that are not expected to be meaningful after a merger to np.nan
                self._apply_nan_attributes(merged_star,
                                           extra_nan_keys=set(["center_gamma",
                                                               "avg_c_in_c_core",
                                                               "envelope_binding_energy"]))

                # Set companion to a massless remnant.
                massless_remnant = convert_star_to_massless_remnant(comp)

        # postMS + MS
        elif (s1 in LIST_ACCEPTABLE_STATES_FOR_POSTMS
            and s2 in LIST_ACCEPTABLE_STATES_FOR_HMS):

                # weighted mixing on the surface abundances of the whole
                # companion with the envelope of star_base
                parameters_to_mix = ["surface_h1", "surface_he4", "surface_c12", "surface_n14", "surface_o16"]
                for abundance_name in parameters_to_mix:
                    setattr(merged_star, abundance_name, self.mass_weighted_avg(star_base,
                                                                                comp,
                                                                                abundance_name=abundance_name,
                                                                                mass_weight1="H-rich_envelope_mass",
                                                                                mass_weight2="mass"))

                # The change of masses occurs after the calculation of weighted averages
                # Note that the he core mass is unchanged, which means that
                # adding the mass changes the envelope mass only.
                merged_star.mass = star_base.mass + comp.mass

                # Set parameters that are not expected to be meaningful after a merger to np.nan
                self._apply_nan_attributes(merged_star)

                # Set burning luminosities and star state.
                merged_star.log_LHe = star_base.log_LHe
                merged_star.log_LZ = star_base.log_LZ
                merged_star.state = check_state_of_star(merged_star, star_CO=False)  # TODO for sure this needs testing!

                massless_remnant = convert_star_to_massless_remnant(comp)

        # MS + postMS (the opposite of above)
        elif (s1 in LIST_ACCEPTABLE_STATES_FOR_HMS
            and s2 in LIST_ACCEPTABLE_STATES_FOR_POSTMS):

                merged_star = comp

                parameters_to_mix = ["surface_h1", "surface_he4", "surface_c12", "surface_n14", "surface_o16"]
                for abundance_name in parameters_to_mix:
                    setattr(merged_star, abundance_name, self.mass_weighted_avg(star_base,
                                                                                comp,
                                                                                abundance_name=abundance_name,
                                                                                mass_weight1="mass",
                                                                                mass_weight2="H-rich_envelope_mass",))

                # The change of masses occurs after the calculation of weighted averages
                merged_star.mass = star_base.mass + comp.mass

                # Set parameters that are not expected to be meaningful after a merger to np.nan
                self._apply_nan_attributes(merged_star)

                # Set burning luminosities and star state.
                merged_star.log_LHe = comp.log_LHe
                merged_star.log_LZ = comp.log_LZ
                merged_star.state = check_state_of_star(merged_star, star_CO=False)

                massless_remnant = convert_star_to_massless_remnant(star_base)

        # postMS + postMS
        elif (s1 in LIST_ACCEPTABLE_STATES_FOR_POSTMS
            and s2 in LIST_ACCEPTABLE_STATES_FOR_POSTMS):
                # Weighted central abundances if merging cores.

                parameters_to_mix = ["center_h1", "center_he4",  "center_c12",  "center_n14",  "center_o16",
                                     'avg_c_in_c_core']

                # two stars with Helium cores
                if star_base.co_core_mass == 0 and comp.co_core_mass == 0:
                    mass_weight1 = 'he_core_mass'
                    mass_weight2 = 'he_core_mass'

                    for abundance_name in parameters_to_mix:
                        setattr(merged_star, abundance_name, self.mass_weighted_avg(star_base,
                                                                                    comp,
                                                                                    abundance_name = abundance_name,
                                                                                    mass_weight1 = mass_weight1,
                                                                                    mass_weight2 = mass_weight2))

                    setattr(merged_star, "center_gamma", np.nan)

                # star_base with CO core and the comp has a He core
                elif (star_base.co_core_mass > 0 and comp.co_core_mass == 0):
                    pass # the central abundances are kept as the ones of star_base

                # Companion with CO core and the star_base has a He core
                elif (comp.co_core_mass > 0 and star_base.co_core_mass == 0):

                    for abundance_name in parameters_to_mix:
                        setattr(merged_star, abundance_name, getattr(comp, abundance_name))

                    setattr(merged_star, "center_gamma", comp.center_gamma)

                # both stars have CO cores
                elif (star_base.co_core_mass > 0 and comp.co_core_mass > 0):

                    #TODO : maybe he_core_mass makes more sense?
                    mass_weight1='co_core_mass'
                    mass_weight2='co_core_mass'

                    for abundance_name in parameters_to_mix:
                        setattr(merged_star, abundance_name, self.mass_weighted_avg(star_base,
                                                                                    comp,
                                                                                    abundance_name=abundance_name,
                                                                                    mass_weight1=mass_weight1,
                                                                                    mass_weight2=mass_weight2))

                    setattr(merged_star, "center_gamma", np.nan)

                else:
                    Pwarn("Unexpected combination of CO core masses during merging", "EvolutionWarning")

                additional_parameter_to_mix = ['surface_h1', 'surface_he4', 'surface_c12', 'surface_n14', 'surface_o16']

                # Weighted mixing on the surface abundances based on the envelopes of the two stars
                for abundance_name in additional_parameter_to_mix:
                    setattr(merged_star, abundance_name, self.mass_weighted_avg(star_base,
                                                                                comp,
                                                                                abundance_name=abundance_name,
                                                                                mass_weight1="H-rich_envelope_mass",
                                                                                mass_weight2="H-rich_envelope_mass"))

                # Add total and core masses after calculations of weighted average
                for key in ["mass", "he_core_mass", "c_core_mass", "o_core_mass", "co_core_mass"]:
                    current = getattr(star_base, key) + getattr(comp, key)
                    setattr(merged_star, key,current)

                # Set parameters that are not expected to be meaningful after a merger to np.nan
                self._apply_nan_attributes(merged_star)

                # TODO for sure this needs testing!
                merged_star.state = check_state_of_star(merged_star, star_CO=False)
                massless_remnant = convert_star_to_massless_remnant(comp)

        # postMS + HeMS
        elif (s1 in LIST_ACCEPTABLE_STATES_FOR_POSTMS
            and s2 in LIST_ACCEPTABLE_STATES_FOR_HeMS):
                # Keep surface abundances from star_base.

                parameters_to_mix = ["center_h1", "center_he4",  "center_c12",  "center_n14",  "center_o16", "avg_c_in_c_core"]

                # Weighted central abundances if merging cores.
                # two stars with only Helium cores, not CO cores
                if star_base.co_core_mass == 0 and comp.co_core_mass == 0:

                    for abundance_name in parameters_to_mix:
                        setattr(merged_star, abundance_name, self.mass_weighted_avg(star_base,
                                                                                    comp,
                                                                                    abundance_name=abundance_name,
                                                                                    mass_weight1="he_core_mass",
                                                                                    mass_weight2="he_core_mass"))

                    setattr(merged_star, "center_gamma", np.nan)

                # star_base with CO core and the comp has just a He core (is a HeMS star)
                # Keep central abundances as those from star_base.
                elif (star_base.co_core_mass > 0 and comp.co_core_mass == 0):
                    pass

                else:
                    Pwarn("Unexpected combination of CO core masses during merging", "EvolutionWarning")

                # add total and core masses after abundance mass weighted calculations
                for key in ["mass", "he_core_mass", "c_core_mass", "o_core_mass", "co_core_mass"]:
                    current = getattr(star_base, key) + getattr(comp, key)
                    setattr(merged_star, key, current)

                # Set parameters that are not expected to be meaningful after a merger to np.nan
                self._apply_nan_attributes(merged_star)

                # TODO for sure this needs testing!
                merged_star.state = check_state_of_star(merged_star, star_CO=False)
                massless_remnant = convert_star_to_massless_remnant(comp)

        # HeMS + postMS (the opposite of above)
        elif (s1 in  LIST_ACCEPTABLE_STATES_FOR_HeMS
            and s2 in LIST_ACCEPTABLE_STATES_FOR_POSTMS):
                # surface abundances from comp
                merged_star = comp
                parameters_to_mix = ["center_h1", "center_he4",  "center_c12",  "center_n14",  "center_o16", "avg_c_in_c_core"]

                # Weighted central abundances if merging cores. Else only from comp
                if star_base.co_core_mass == 0 and comp.co_core_mass == 0: # two stars with only Helium cores, not CO cores

                    for abundance_name in parameters_to_mix:
                        setattr(merged_star, abundance_name, self.mass_weighted_avg(star_base,
                                                                                    comp,
                                                                                    abundance_name=abundance_name,
                                                                                    mass_weight1="he_core_mass",
                                                                                    mass_weight2="he_core_mass"))

                    setattr(merged_star, "center_gamma", np.nan)

                # star_base is the HeMS Star and comp has a CO core
                # the central abundances are kept from comp
                elif (star_base.co_core_mass == 0 and comp.co_core_mass > 0):
                    pass
                else:
                    Pwarn("Unexpected combination of CO core masses during merging", "EvolutionWarning")

                # add total and core masses after weighted averages above
                for key in ["mass", "he_core_mass", "c_core_mass", "o_core_mass", "co_core_mass"]:
                    current = getattr(star_base, key) + getattr(comp, key)
                    setattr(merged_star, key, current)

                # Set parameters that are not expected to be meaningful after a merger to np.nan
                self._apply_nan_attributes(merged_star)

                # TODO for sure this needs testing!
                merged_star.state = check_state_of_star(merged_star, star_CO=False)
                massless_remnant = convert_star_to_massless_remnant(star_base)

        # postMS + postHeMS
        elif (s1 in LIST_ACCEPTABLE_STATES_FOR_POSTMS
            and s2 in LIST_ACCEPTABLE_STATES_FOR_POSTHeMS):

                parameters_to_mix = ["center_h1", "center_he4",  "center_c12",  "center_n14",  "center_o16", "avg_c_in_c_core"]

                # Weighted central abundances if merging cores
                # Two stars with Helium cores
                if star_base.co_core_mass == 0 and comp.co_core_mass == 0:

                    for abundance_name in parameters_to_mix:
                        setattr(merged_star, abundance_name,
                                self.mass_weighted_avg(star_base,
                                                       comp,
                                                       abundance_name=abundance_name,
                                                       mass_weight1="he_core_mass",
                                                       mass_weight2="he_core_mass"))

                    setattr(merged_star, "center_gamma", np.nan)

                # Star_base with CO core and the comp has a He core
                # the central abundances are kept as the ones of star_base
                elif (star_base.co_core_mass > 0 and comp.co_core_mass == 0):
                    pass

                # comp with CO core and the star_base has a He core
                elif (comp.co_core_mass > 0 and star_base.co_core_mass == 0):

                    for abundance_name in parameters_to_mix:
                        setattr(merged_star, abundance_name, getattr(comp, abundance_name))

                    setattr(merged_star, "center_gamma", comp.center_gamma)

                # both stars have CO cores
                elif (star_base.co_core_mass > 0 and comp.co_core_mass > 0):

                    for abundance_name in parameters_to_mix:
                        setattr(merged_star, abundance_name,
                                self.mass_weighted_avg(star_base,
                                                       comp,
                                                       abundance_name=abundance_name,
                                                       mass_weight1="co_core_mass",
                                                       mass_weight2="co_core_mass"))

                    setattr(merged_star, "center_gamma", np.nan)

                else:
                    Pwarn("Unexpected combination of CO core masses during merging", "EvolutionWarning")

                # add total and core masses
                for key in ["mass", "he_core_mass", "c_core_mass", "o_core_mass", "co_core_mass"]:
                    current = getattr(star_base, key) + getattr(comp, key)
                    setattr(merged_star, key, current)

                # Set parameters that are not expected to be meaningful after a merger to np.nan
                self._apply_nan_attributes(merged_star)

                # TODO for sure this needs testing!
                merged_star.state = check_state_of_star(merged_star, star_CO=False)
                massless_remnant = convert_star_to_massless_remnant(comp)

        # postHeMS + postMS (the opposite of above)
        elif (s1 in LIST_ACCEPTABLE_STATES_FOR_POSTHeMS
            and s2 in LIST_ACCEPTABLE_STATES_FOR_POSTMS):

                merged_star = comp

                parameters_to_mix = ["center_h1", "center_he4",  "center_c12",  "center_n14",  "center_o16", "avg_c_in_c_core"]

                # weighted central abundances if merging cores
                # two stars with Helium cores
                if star_base.co_core_mass == 0 and comp.co_core_mass == 0:

                    for abundance_name in parameters_to_mix:
                        setattr(merged_star, abundance_name, self.mass_weighted_avg(star_base,
                                                                                    comp,
                                                                                    abundance_name=abundance_name,
                                                                                    mass_weight1="he_core_mass",
                                                                                    mass_weight2="he_core_mass"))

                    setattr(merged_star, "center_gamma", np.nan)

                # star_base with CO core and the comp has a He core
                elif (star_base.co_core_mass > 0 and comp.co_core_mass == 0):

                    for abundance_name in parameters_to_mix:
                        setattr(merged_star, abundance_name, getattr(star_base, abundance_name))

                    setattr(merged_star, "center_gamma", star_base.center_gamma)

                # comp with CO core and the star_base has a He core
                # the central abundances are kept as the ones of comp
                elif (comp.co_core_mass > 0 and star_base.co_core_mass == 0):
                    pass

                # both stars have CO cores
                elif (star_base.co_core_mass > 0 and comp.co_core_mass > 0):

                    for abundance_name in parameters_to_mix:
                        setattr(merged_star, abundance_name, self.mass_weighted_avg(star_base, comp, abundance_name=abundance_name, mass_weight1="co_core_mass", mass_weight2="co_core_mass"))

                    setattr(merged_star, "center_gamma", np.nan)

                else:
                    Pwarn("Unexpected combination of CO core masses during merging", "EvolutionWarning")

                # add total and core masses
                for key in ["mass", "he_core_mass", "c_core_mass", "o_core_mass", "co_core_mass"]:
                    current = getattr(star_base, key) + getattr(comp, key)
                    setattr(merged_star, key,current)

                # Set parameters that are not expected to be meaningful after a merger to np.nan
                self._apply_nan_attributes(merged_star)

                # TODO for sure this needs testing!
                merged_star.state = check_state_of_star(merged_star, star_CO=False)
                massless_remnant = convert_star_to_massless_remnant(star_base)

        # He star + He star
        elif (s1 in STAR_STATES_HE_RICH
            and s2 in STAR_STATES_HE_RICH):

                parameters_to_mix = ["center_h1", "center_he4",  "center_c12",  "center_n14",  "center_o16", "avg_c_in_c_core"]
                # weighted central abundances if merging cores. Else only from star_base
                # two stars with Helium cores
                if star_base.co_core_mass == 0 and comp.co_core_mass == 0:

                    for abundance_name in parameters_to_mix:
                        setattr(merged_star, abundance_name, self.mass_weighted_avg(star_base,
                                                                                     comp,
                                                                                     abundance_name=abundance_name,
                                                                                     mass_weight1="he_core_mass",
                                                                                     mass_weight2="he_core_mass"))

                    setattr(merged_star, "center_gamma", np.nan)

                # star_base with CO core and the comp has a He core
                # the central abundances are kept as the ones of star_base
                elif (star_base.co_core_mass > 0 and comp.co_core_mass == 0):
                    pass

                # comp with CO core and the star_base has a He core
                elif (comp.co_core_mass > 0 and star_base.co_core_mass == 0):

                    for abundance_name in parameters_to_mix:
                        setattr(merged_star, abundance_name, getattr(comp, abundance_name))

                    setattr(merged_star, "center_gamma", comp.center_gamma)

                # both stars have CO cores
                elif (star_base.co_core_mass > 0 and comp.co_core_mass > 0):

                    for abundance_name in parameters_to_mix:
                        setattr(merged_star, abundance_name, self.mass_weighted_avg(star_base,
                                                                                    comp,
                                                                                    abundance_name=abundance_name,
                                                                                    mass_weight1="co_core_mass",
                                                                                    mass_weight2="co_core_mass"))

                    setattr(merged_star, "center_gamma", np.nan)

                else:
                    Pwarn("Unexpected combination of CO core masses during merging", "EvolutionWarning")

                # add total and core masses
                for key in ["mass", "he_core_mass", "c_core_mass", "o_core_mass", "co_core_mass"]:
                    current = getattr(star_base, key) + getattr(comp, key)
                    setattr(merged_star, key,current)

                # weigheted mixing on the surface abundances based on the He-rich envelopes of the two stars
                additional_parameter_to_mix = ['surface_h1', 'surface_he4', 'surface_c12', 'surface_n14', 'surface_o16']
                for abundance_name in additional_parameter_to_mix:
                    setattr(merged_star, abundance_name, self.mass_weighted_avg(star_base,
                                                                                comp,
                                                                                abundance_name=abundance_name,
                                                                                mass_weight1="He-rich_envelope_mass",
                                                                                mass_weight2="He-rich_envelope_mass"))

                # Set parameters that are not expected to be meaningful after a merger to np.nan
                self._apply_nan_attributes(merged_star)

                # TODO for sure this needs testing!
                merged_star.state = check_state_of_star(merged_star, star_CO=False)
                massless_remnant = convert_star_to_massless_remnant(comp)

        # Star + WD
        elif (s1 in STAR_STATES_NOT_CO
            and s2 in ["WD"]):

                parameters_to_mix = ["center_h1", "center_he4",  "center_c12",  "center_n14",  "center_o16", "avg_c_in_c_core"]
                # WD is considered a stripped CO core
                # WD would always be the comp, it cannot be the engulfing star (so no need to do the opposite stars case below)

                # weighted central abundances if merging cores. Else only from star_base
                # comp with CO core and the star_base has not
                if (comp.co_core_mass is not None and star_base.co_core_mass == 0):

                    for abundance_name in parameters_to_mix:
                        setattr(merged_star, abundance_name, getattr(comp, abundance_name))

                    setattr(merged_star, "center_gamma", comp.center_gamma)


                # Star_base with CO core and the comp is a WD (so has a CO core)
                elif (comp.co_core_mass is not None and star_base.co_core_mass > 0):

                    for abundance_name in parameters_to_mix:
                        setattr(merged_star, abundance_name, self.mass_weighted_avg(star_base,
                                                                                     comp,
                                                                                     abundance_name=abundance_name,
                                                                                     mass_weight1="co_core_mass",
                                                                                     mass_weight2="co_core_mass"))

                    setattr(merged_star, "center_gamma", np.nan)

                else:
                    Pwarn("Unexpected combination of CO core masses during merging", "EvolutionWarning")

                # add total and core masses
                for key in ["mass", "he_core_mass", "c_core_mass", "o_core_mass", "co_core_mass"]:
                    current = getattr(star_base, key) + getattr(comp, "mass")
                    setattr(merged_star, key,current)

                # Set parameters that are not expected to be meaningful after a merger to np.nan
                self._apply_nan_attributes(merged_star)

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
