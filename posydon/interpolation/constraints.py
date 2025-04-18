"""Module to enforce constraints on interpolated quantities.

Constraints are split up into 3 types:

 1. Astrophysical laws that must hold
 2. Inequalities between values interpolated upon which must hold
 3. Sums of values interpolated upon which must satisfy some inequality

For each type of constraint, different parameter dictionaries are required.

For type 1 constraints a dictionary with the following fields is required:

    type - the type of constraint
    name - the name of the constraint
    dependent - the field for which the constraint must be calculated
    independent - the field names which are used to compute the dependent
    function - a function that computes the dependent from the independent

For type 2 constraints a dictionary with the following fields is required:

    type - the type of constraint
    name - the name of the constraint
    fields - a multi dimensional array containing the field names in the
             desired descending order (if usage with default enforcer is
             wanted)
    log - a boolean denoting whether or not the largest (first) field in
          each constraint list is in log scale (by default false)
    function - a function that enforces the desired constraint (by default
               inequality is corrected)

For type 3 constraints a dictionary with the following fields is required:

    type - the type of constraint
    name - the name of the constraint
    sum - an array of fields whose sum needs to satisfy some inequality
    constraint - either a real number or fieldname that the sum needs to be
                 constrained by
"""


__authors__ = [
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Philipp Moura Srivastava <philipp.msrivastava@gmail.com>",
    "Emmanouil Zapartas <ezapartas@gmail.com>",
]


import math
import numpy as np
from posydon.utils.common_functions import (stefan_boltzmann_law,
                                            orbital_separation_from_period)


# toggle this flag to enable/disable constraints (used for debugging)
INTERPOLATION_CONSTRAINTS_ON = True

# which quantities are excluded from the interpolator
# QUANTITIES_EXCLUDED_FROM_INTERP = ["S1_log_Teff", "S2_log_Teff"]
QUANTITIES_EXCLUDED_FROM_INTERP = []

EPS = 1.0e-9   # to avoid numerical precision errors


def boltzman_constraint_func(log_L, log_R):
    """Boltzmann law constraint."""
    return np.log10(stefan_boltzmann_law(10.0 ** log_L, 10.0 ** log_R))


def kepler_3rd_law_func(period_days, mass1, mass2):
    """Kepler's 3rd law constraint."""
    return orbital_separation_from_period(period_days, mass1, mass2)


def Lnuc_func(log_LH, log_LHe, log_LZ):
    """Total nuclear luminosity should be the sum of H, He and the rest."""
    return np.log10(10.0**log_LH + 10.0**log_LHe + 10.0**log_LZ)


def lg_system_mdot_func(lg_system_mdot, lg_mtransfer_rate):
    """System mdot must be less than the total mass-transfer rate."""
    return min(lg_system_mdot, lg_mtransfer_rate)


def xfer_fraction_func(lg_system_mdot_2, lg_mtransfer_rate):
    """Infer the xfer_fraction from the system and mass-transfer rate."""
    if lg_system_mdot_2 < -98.0 or lg_mtransfer_rate < -98.0:
        return 1.0
    if np.isfinite(lg_system_mdot_2 + lg_mtransfer_rate):
        assert(lg_system_mdot_2 <= lg_mtransfer_rate)
        return 1.0 - 10.0 ** (lg_system_mdot_2 - lg_mtransfer_rate)
    return np.nan


def correct_inequality(fv_dict, key_sets, star, log=False):
    """Apply an "unequality constraint" (type 2) if needed.

    Parameters
    ----------
    fv_dict : dict
        The final values dictionary.
    key_sets : list of lists
        List of collections of keys connected with >= inequalities (descending
        order). E.g., if one constraint is `a >= b >= c`, then key_sets is
        [[a, b, c], ...].
    star : int
        For which star (1 or 2), the constraint will be applied. E.g.,
        `star_{}_mass` will be either `star_1_mass` or `star_2_mass`.
    log : bool
        If True, then the largest parameter (first key), is given in log10
        scale, while all the rest are in linear scale.

    Returns
    -------
    dict
        The sanitized version of `fv_dict`.

    """
    for key_set in key_sets:

        vals = [fv_dict[key] for key in key_set]

        if log:         # taking care of logged variables
            vals[0] = 10**vals[0]

        if vals != sorted(vals, reverse=True):
            for indx, (val, key) in enumerate(zip(vals, key_set)):
                if(indx != 0 and vals[indx - 1] < val):
                    fv_dict[key] = vals[indx - 1]
                    vals[indx] = fv_dict[key]

    return fv_dict


def correct_sum(fv_dict, sum_keys, constraint):
    """Normalize quantities should their sum does not have the expected value.

    Parameters
    ----------
    fv_dict : dict
        The final values dictionary.
    sum_keys : array-like
        Which quantities should add up to a specific value.
    constraint : float or str
        If float, then the sum of the quantities should be eqaul to that.
        If str, the fieldname of the quanity that the sum is constrained by.

    Returns
    -------
    dict
        The sanitized version of the final values dictionary.

    """
    keys_sum = sum([fv_dict[key] for key in sum_keys])

    target = fv_dict[constraint] if isinstance(constraint, str) else constraint

    if abs(keys_sum - target) > EPS:
        for key in sum_keys:
            fv_dict[key] = constraint * (fv_dict[key] / keys_sum)

    return fv_dict


def get_thickness_n_radius(fv_dict, constraint, keys):
    """Corrects thickness and radius constraint violation of one is present.

    Parameters
    ----------
    fv_dict : dict
        The final values dictionary
    constraint : dict
        The dictionary specifying the constraint to be applied
    keys : list
        A list of formatted key names

    Returns
    -------
    dict
        The sanitized version of the final values dictionary.

    """
    r = fv_dict[keys[2]]
    t = fv_dict[keys[1]]
    R = 10**fv_dict[keys[0]]
    if not np.isfinite(r+t+R):      # if any of those is nan
        return fv_dict

    if(not (t > 0 and r >= 0.5 * t and r + 0.5 * t <= R)):

        centroid = (R / 3, R / 2)   # centroid of convex hull

        if(centroid[1] == r):
            fv_dict[keys[2]] = 0.5 * R
            fv_dict[keys[1]] = R
        else:
            slope = ((r - centroid[1])/(t - centroid[0]))

            coeff_mat_one = np.array([[-1, slope],
                                      [1, -0.5]])

            dep_vec_one = np.array([-centroid[1] + slope * centroid[0], 0])

            coeff_mat_two = np.array([[-1, slope],
                                      [1, 0.5]])

            dep_vec_two = np.array([-centroid[1] + slope * centroid[0], R])

            inters = [np.linalg.solve(coeff_mat_one, dep_vec_one),
                      np.linalg.solve(coeff_mat_two, dep_vec_two)]

            valid_inters = [inter for inter in inters
                            if inter[1] > 0
                            and inter[0] >= 0.5 * inter[1] - EPS
                            and inter[0] + 0.5 * inter[1] <= R + EPS]

            assert len(valid_inters) > 0

            min_indx = np.argmin([math.sqrt((vi[0] - r)**2 + (vi[1] - t)**2)
                                  for vi in valid_inters])

            fv_dict[keys[1]] = valid_inters[min_indx][1]
            fv_dict[keys[2]] = valid_inters[min_indx][0]

    return fv_dict


CONSTRAINTS = [
    {
        "type": 1,
        "name": "Boltzmann",
        "dependent": "S{}_log_Teff",
        "independent": ["S{}_log_L", "S{}_log_R"],
        "function": boltzman_constraint_func
    },
    {
        "type": 1,
        "name": "Kepler's 3rd law",
        "dependent": "binary_separation",
        "independent": ["period_days", "star_1_mass", "star_2_mass"],
        "function": kepler_3rd_law_func
    },
    {
        "type": 1,
        "name": "Total nuclear luminosity",
        "dependent": "S{}_log_Lnuc",
        "independent": ["S{}_log_LH", "S{}_log_LHe", "S{}_log_LZ"],
        "function": Lnuc_func
    },
    {   # IMPORTANT: this should precede the `xfer_fraction` constraint
        "type": 1,
        "name": "lg_system_mdot <= lg_mtransfer_rate",
        "dependent": "lg_system_mdot_{}",
        "independent": ["lg_system_mdot_{}", "lg_mtransfer_rate"],
        "function": lg_system_mdot_func
    },
    {
        "type": 1,
        "name": "xfer_fraction",
        "dependent": "xfer_fraction",
        "independent": ["lg_system_mdot_2", "lg_mtransfer_rate"],
        "function": xfer_fraction_func
    },
    {
        "type": 2,
        "name": "Element Mass Inequalities Constraint (2.2)",
        "fields": [["star_{}_mass", "S{}_he_core_mass", "S{}_co_core_mass"]],
        "log": False
    },
    {
        "type": 2,
        "name": "Element Radius Inequalities Constraint (2.3)",
        "fields": [["S{}_log_R", "S{}_he_core_radius", "S{}_co_core_radius"]],
        "log": True
    },
    {
        "type": 2,
        "name": "2.4",
        "fields": [
            ["star_{}_mass", "S{}_m_core_CE_1cent"],
            ["star_{}_mass", "S{}_m_core_CE_10cent"],
            ["star_{}_mass", "S{}_m_core_CE_30cent"],
            ["star_{}_mass", "S{}_m_core_CE_pure_He_star_10cent"]
        ],
        "log": False
    },
    {
        "type": 2,
        "name": "2.5",
        "fields": [
            ["S{}_log_R", "S{}_r_core_CE_1cent"],
            ["S{}_log_R", "S{}_r_core_CE_10cent"],
            ["S{}_log_R", "S{}_r_core_CE_30cent"],
            ["S{}_log_R", "S{}_r_core_CE_pure_He_star_10cent"]
        ],
        "log": True
    },
    {
        "type": 2,
        "name": "2.6",
        "fields": [["star_{}_mass", "S{}_mass_conv_reg_fortides"]],
        "log": False
    },
    {
        "type": 2,
        "name": "2.7",
        "fields": [["S{}_log_R",
                    "S{}_thickness_conv_reg_fortides",
                    "S{}_radius_conv_reg_fortides"]],
        "log": False,
        "function": get_thickness_n_radius
    },
    {
        "type": 2,
        "name": "2.9",
        "fields": [
            ["star_{}_mass", "S{}_direct_mass"],
            ["star_{}_mass", "S{}_Fryer+12-rapid_mass"],
            ["star_{}_mass", "S{}_Fryer+12-delayed_mass"],
            ["star_{}_mass", "S{}_Sukhbold+16-engineN20_mass"],
            ["star_{}_mass", "S{}_Patton&Sukhbold20-engineN20_mass"]
        ],
        "log": False
    },
    {
        "type": 2,
        "name": "2.11",
        "fields": [["lg_mtransfer_rate", "lg_system_mdot_{}"]],
        "log": False
    },
    {
        "type": 2,
        "name": "Mass in Hydrogen/Helium (2.12)",
        "fields": [
            ["star_{}_mass", "S{}_total_mass_h1"],
            ["star_{}_mass", "S{}_total_mass_he4"]
        ],
        "log": False
    },
    {
        "type": 3,
        "name": "Center Sum Equal One",
        "sum_fields": ["S{}_center_h1", "S{}_center_he4", "S{}_center_c12",
                       "S{}_center_n14", "S{}_center_o16", "S{}_center_other"],
        "constraint": 1.0
    },
    {
        "type": 3,
        "name": "Surface Sum Equal One",
        "sum_fields": ["S{}_surface_h1", "S{}_surface_he4",
                       "S{}_surface_c12", "S{}_surface_n14",
                       "S{}_surface_o16", "S{}_surface_other"],
        "constraint": 1.0
    }
]


def check_order_of_constraints():
    """Report possible issues with the order of application of constraints.

    TODO: this needs to be updated for the new types of constraints.
    """
    changed_so_far, used_so_far = set(), set()

    for constraint in CONSTRAINTS:
        # find all colnames for the independent variables (including per star)
        indep_cols = set([col.format(star) for star in [1, 2]
                         for col in constraint.independent])
        # same with dependent variables
        dep_cols = set([constraint.dependent.format(star) for star in [1, 2]])

        # track quantities relevant to this constraint, already encountered
        already_changed = changed_so_far.intersection(indep_cols)
        about_to_change = used_so_far.intersection(dep_cols)

        # if any, notify the user
        if len(already_changed) > 0:
            print("In constraint {}, independent quanitities ({}) may have "
                  "been modified by a previous constraint.".format(
                      constraint.name, already_changed))
        if len(about_to_change) > 0:
            print("Constraint {} may modify previously used dependent "
                  "quanitities ({}).".format(constraint.name, about_to_change))

        # you've seen these now...
        changed_so_far.update(dep_cols)
        used_so_far.upate(indep_cols)


def find_constraints_to_apply(out_keys):
    """Find out which constraints can be applied.

    Find out which constraints can be applied based on the out_keys
    specified by the user. This function assumes that the keys for star 1
    have a corresponding key in star 2 when applicable.

    Parameters
    ----------
    out_keys : list of strings
        The keys specified in the IFInterpolator

    Returns
    -------
    list of dicts
        A list of constraints that can and will be applied

    """
    constraints_to_apply = []

    for constraint in CONSTRAINTS:

        should_apply = True

        if constraint["type"] == 1:

            for field in constraint["independent"]:
                if field.format(1) not in out_keys:
                    should_apply = False

        elif constraint["type"] == 2:

            new_fields = []

            for fields in constraint["fields"]:
                add_subconstraint = True

                for field in fields:
                    if field.format(1) not in out_keys:
                        add_subconstraint = False

                if add_subconstraint:
                    new_fields.append(fields)

            if len(new_fields) == 0:
                should_apply = False
            else:
                constraint["fields"] = new_fields

        elif constraint["type"] == 3:

            for field in constraint["sum_fields"]:
                if field.format(1) not in out_keys:
                    should_apply = False

        if should_apply:
            constraints_to_apply.append(constraint)

    return constraints_to_apply


def sanitize_interpolated_quantities(fvalues, constraints, verbose=False):
    """Apply constraints in final values returned from the IF interpolator.

    Parameters
    ----------
    fvalues : dict
        Dictionary containing the final values returned from the initial-final
        interpolator evaluation method.
    constraints : list of dicts
        A list of constraints that are to be applied. Can be computed using the
        find_constraints_to_apply function
    verbose: bool
        If True, report the actions.

    Returns
    -------
    dict
        The dictionary of sanitized final values.

    """
    if not INTERPOLATION_CONSTRAINTS_ON:
        return fvalues

    sanitized = fvalues.copy()

    stars_present = [1, 2] if ("S1_log_R" in sanitized.keys()
                               and "S2_log_R" in sanitized.keys()) else [1]

    for constraint in constraints:
        if verbose:
            print("Enforcing constraint {}".format(constraint["name"]))

        if constraint["type"] == 1:
            y_colname = constraint["dependent"]
            per_star = "{}" in y_colname
            both_stars = per_star and y_colname.format(2) in sanitized.keys()
            stars = [1, 2] if (per_star and both_stars) else [1]
            for star in stars:
                x_cols = [x_colname.format(star)
                          for x_colname in constraint["independent"]]
                if verbose:
                    print("    {} -> {}".format(
                        x_cols, y_colname.format(star)))
                x_vals = [sanitized[x_col] for x_col in x_cols]
                if np.any(~np.isfinite(x_vals)):
                    y_value = np.nan
                else:
                    y_value = constraint["function"](*x_vals)
                sanitized[y_colname.format(star)] = y_value
        elif constraint["type"] == 2:
            stars = stars_present
            for star in stars:
                # creating fieldnames
                key_sets = [[key.format(star) for key in key_set]
                            for key_set in constraint["fields"]]

                function = constraint.get("function")
                if function is None:
                    sanitized = correct_inequality(sanitized, key_sets, star,
                                                   constraint["log"])
                else:
                    for key_set in key_sets:
                        sanitized = function(sanitized, constraint, key_set)
        elif constraint["type"] == 3:
            stars = stars_present
            for star in stars:
                sum_keys = [key.format(star)
                            for key in constraint["sum_fields"]]

                sanitized = correct_sum(sanitized, sum_keys,
                                        constraint["constraint"])

    return sanitized
