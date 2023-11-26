"""Default spectal synthesis and spectral grid parameters parameters."""


__authors__ = [
    "Eirini Kasdagli <kasdaglie@ufl.edu>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
]


default_grid_kwargs = {
    # Default grid options:

    # The main grid option:
    'main_grid':'sg-C3K-coarse.h5', #"sg-CAP18-coarse.h5",

    # Secondary grid option to capture the failures of the main grid.
    'secondary_grid': "sg-BSTAR2006-ufl-vtb2.h5",
    'stripped_grid': "sg-Gotberg18.h5",
    'bstar_grid' : "sg-BSTAR2006-ufl-vtb2.h5",#,'sg-C3K-coarse.h5'
    'ostar_grid': "sg-OSTAR2002-ufl.h5",


    # Setting the wavelengths range:
    'lam_min': 105.,
    'lam_max': 50000.,

    # Number of points in the wavelength range:
    'lam_res': 300000,
    # Performance variables:
    'cache_limit': 256,
    'filters': ['U', 'B', 'V']
}
default_kwargs = {
    'max_number_of_binaries': None,
    #int
    'save_data': False,
    #False,True

    #The desired path for the output file.
    'output_file_path':None,
    #The default option is ./

    # Ostar Temp cutoff:
    'ostar_temp_cut_off': 28000,
    #Bstar Temp cutoff:
    'bstar_temp_cut_off': 15000,
 }
