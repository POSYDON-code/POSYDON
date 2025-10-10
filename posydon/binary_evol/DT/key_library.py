"""Dictionary containing data column name (key) translations between POSYDON h5
file PSyGrid data names (items) and MESA data names (keys).

"""

__authors__ = [
    "Seth Gossage <seth.gossage@northwestern.edu>"
]

# translates MESA stellar data column names into 
# what we use for property names in interpolators.
# Most are the same, a few are different.
DEFAULT_TRANSLATION = {
    "time": "time",
    "orbital_period": "porb",
    "eccentricity": "ecc",
    "separation": "sep",
    "state": None,
    "event": None,
    "rl_relative_overflow_1": "rl_relative_overflow_1",
    "rl_relative_overflow_2": "rl_relative_overflow_2",
    "lg_mtransfer_rate": "lg_mtransfer_rate",
    "V_sys": None,
    "mass": "mass",
    "log_R": "log_R",
    "R": "R",
    "lg_mdot": "mdot",
    "log_L": "log_L",
    "lg_wind_mdot": "mdot",
    "lg_system_mdot": "lg_mdot",
    "he_core_mass": "he_core_mass",
    "he_core_radius": "he_core_radius",
    "c_core_mass": "c_core_mass",
    "c_core_radius": "c_core_radius",
    "o_core_mass": "o_core_mass",
    "o_core_radius": "o_core_radius",
    "center_h1": "center_h1",
    "center_he4": "center_he4",
    "center_c12": "center_c12",
    "center_o16": "center_o16",
    "center_n14": "center_n14",
    "surface_h1": "surface_h1",
    "surface_he4": "surface_he4",
    "surface_c12": "surface_c12",
    "surface_n14": "surface_n14",
    "surface_o16": "surface_o16",
    "center_gamma": "center_gamma",
    "log_LH": "log_LH",
    "log_LHe": "log_LHe",
    "log_LZ": "log_LZ",
    "log_Lnuc": "log_Lnuc",
    "c12_c12": "c12_c12",
    "avg_c_in_c_core": "avg_c_in_c_core",
    "surf_avg_omega_div_omega_crit": "surf_avg_omega_div_omega_crit",
    "surf_avg_omega": "omega",
    "total_moment_of_inertia": "inertia",
    "log_total_angular_momentum": "log_total_angular_momentum",
    "profile": None,
    "metallicity": None,
    "spin": "spin_parameter",
    "conv_env_top_mass": "conv_env_top_mass",
    "conv_env_bot_mass": "conv_env_bot_mass",
    "conv_env_top_radius": "conv_env_top_radius",
    "conv_env_bot_radius": "conv_env_bot_radius",
    "conv_env_turnover_time_g": "conv_env_turnover_time_g",
    "conv_env_turnover_time_l_b": "conv_env_turnover_time_l_b",
    "conv_env_turnover_time_l_t": "conv_env_turnover_time_l_t",
    "envelope_binding_energy": "envelope_binding_energy",
    "mass_conv_reg_fortides": "mass_conv_reg_fortides",
    "thickness_conv_reg_fortides": "thickness_conv_reg_fortides",
    "radius_conv_reg_fortides": "radius_conv_reg_fortides",
    "lambda_CE_1cent": "lambda_CE_1cent",
    "lambda_CE_10cent": "lambda_CE_10cent",
    "lambda_CE_30cent": "lambda_CE_30cent",
    "co_core_mass": "co_core_mass",
    "co_core_radius": "co_core_radius",
    "lambda_CE_pure_He_star_10cent": "lambda_CE_pure_He_star_10cent",
    "trap_radius": "trap_radius",
    "acc_radius": "acc_radius",
    "t_sync_rad_1": "t_sync_rad_1",
    "t_sync_conv_1": "t_sync_conv_1",
    "t_sync_rad_2": "t_sync_rad_2",
    "t_sync_conv_2": "t_sync_conv_2",
    "mass_transfer_case": None,
    "nearest_neighbour_distance": None,
    "total_mass_h1": "total_mass_h1",
    "total_mass_he4": "total_mass_he4",
}

# keys from POSYDON h5 file translated to MESA data column names
DEFAULT_TRANSLATED_KEYS = (
    'age',
    'mass',
    'mdot',
    'inertia',
    'conv_mx1_top_r',
    'conv_mx1_bot_r',
    'surface_h1',
    'center_h1',
    'mass_conv_reg_fortides',
    'thickness_conv_reg_fortides',
    'radius_conv_reg_fortides',
    'log_Teff',
    'surface_he3',
    'surface_he4',
    'center_he4',
    'avg_c_in_c_core',
    'log_LH',
    'log_LHe',
    'log_LZ',
    'log_Lnuc',
    'c12_c12',
    'center_c12',
    'he_core_mass',
    'log_L',
    'log_R',
    'c_core_mass',
    'o_core_mass',
    'co_core_mass',
    'c_core_radius',
    'o_core_radius',
    'co_core_radius',
    'spin_parameter',
    'log_total_angular_momentum',
    'center_n14',
    'center_o16',
    'surface_n14',
    'surface_o16',
    'conv_env_top_mass',
    'conv_env_bot_mass',
    'conv_env_top_radius',
    'conv_env_bot_radius',
    'conv_env_turnover_time_g',
    'conv_env_turnover_time_l_b',
    'conv_env_turnover_time_l_t',
    'envelope_binding_energy',
    'lambda_CE_1cent',
    'lambda_CE_10cent',
    'lambda_CE_30cent',
    'lambda_CE_pure_He_star_10cent',
    'center_gamma',
    'total_mass_h1',
    'total_mass_he4'
)

KEYS_POSITIVE = (
    'mass_conv_reg_fortides',
    'thickness_conv_reg_fortides',
    'radius_conv_reg_fortides'
)

DEFAULT_PROFILE_KEYS = (
    'radius',
    'mass',
    'logRho',
    'energy',
    'x_mass_fraction_H',
    'y_mass_fraction_He',
    'z_mass_fraction_metals',
    'neutral_fraction_H',
    'neutral_fraction_He',
    'avg_charge_He'
)

"""
# states to ID an HMS star
# TODO: build these from the other lists to (hopefully) 
#       ensure consistency
STAR_STATES_FOR_HMS_MATCHING = ["H-rich_Core_H_burning",
                                "accreted_He_Core_H_burning"]

# states to ID a postMS star
STAR_STATES_FOR_postMS_MATCHING = [
    "H-rich_Shell_H_burning",
    "H-rich_Core_He_burning",
    "H-rich_Central_He_depleted",
    "H-rich_Core_C_burning",
    "H-rich_Central_C_depletion",
    "H-rich_non_burning",
    "accreted_He_non_burning"]

# states to ID an He star
STAR_STATES_FOR_strippedHe_MATCHING = [
    'accreted_He_Core_He_burning',
    'stripped_He_Core_He_burning',
    'stripped_He_Shell_He_burning',     # includes stars burning C in core, does this exist?
    'stripped_He_Central_He_depleted',  # includes stars burning C in core
    'stripped_He_Central_C_depletion',
    'stripped_He_non_burning'
    ]
"""
