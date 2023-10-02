"""Default spectal synthesis and spectral grid parameters parameters."""


__authors__ = [
    "Eirini Kasdagli <kasdaglie@ufl.edu>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
]


default_grid_kwargs = {
    """Write docstring here."""
    # Default grid options:

    # The main grid option:
    'main_grid': "sg-CAP18-coarse.h5",

    # Secondary grid option to capture the failures of the main grid.
    'secondary_grid': "sg-BSTAR2006-medium.h5",
    'stripped_grid': "sg-Gotberg18.h5",
    'ostar_grid': "sg-OSTAR2002-medium.h5",

    # Setting the wavelengths range:
    'lam_min': 3000.,
    'lam_max': 7000.,

    # Number of points in the wavelength range:
    'lam_res': 2000,
    # Performance variables:
    'cache_limit': 256,
    # Specific Temp cutoff:
    'ostar_temp_cut_off': 27000,
    'filters': ['U', 'B', 'V']
}
default_kwargs = {

 }
