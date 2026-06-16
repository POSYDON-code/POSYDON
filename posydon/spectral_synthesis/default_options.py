"""Default spectal synthesis and spectral grid parameters parameters."""


__authors__ = [
    "Eirini Kasdagli <kasdaglie@ufl.edu>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
]


default_grid_kwargs = {
    # Default grid options:

    # The main grid option:
    'main_grid':'sg-C3K-demo.h5',#'sg-C3K-coarse.h5',"sg-CAP18-coarse.h5",
    # Secondary grid option to capture the failures of the main grid.
    'secondary_grid': 'sg-BOSZ-MARCS-v1.h5',
    'stripped_grid': "sg-Gotberg23.h5",
    'bstar_grid' :'sg-BSTAR2006-ufl-vtb2.h5',# "sg-BSTAR2006-ufl-vtb2.h5",#,'sg-C3K-coarse.h5'
    'ostar_grid': "sg-OSTAR2002-ufl.h5",
    'WR_grid' : "sg-PoWR-WNL-H20-new.h5",
    'WNL_grid' :"sg-PoWR-MW-WNL.h5",
    'WNE_grid' : "sg-PoWR-MW-WNE.h5",
    'WC_grid' : "sg-PoWR-MW-WC.h5",

    # Setting the wavelengths range:
    'lam_min': 110.,
    'lam_max': 100000,

    # Number of points in the wavelength range:
    'lam_res': 300000,
    # Performance variables:
    'cache_limit': 512,
    'filters': ['U', 'B', 'V']
}
default_kwargs = {
    'max_number_of_binaries': None,
    #int
    #
    'save_population_data': True,
    'save_data': True,
    #False,True

    #The desired path for the output file.
    'output_file_path':None,
    #The default option is ./

    # Accepts a list of states from state_list, or None to include all.
    # e.g. ['detached', 'merged']
    'include_states': None,       # None = include all binary states
    ## Accepts a list of spectral types from libraries, or None to include all.
    # e.g. ['stripped_grid', 'ostar_grid']
    'include_spectral_types': None,  # None = include all spectral types

    # Ostar Temp cutoff:
    'ostar_temp_cut_off': 28000,
    #Bstar Temp cutoff:
    'bstar_temp_cut_off': 15000,
 }
