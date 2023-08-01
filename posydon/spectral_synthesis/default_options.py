
default_kwargs = {
    #Default grid options:
    #The main grid option: 
    'main_grid_file' : "sg-CAP18-coarse.h5",
    #Secondary grid option to capture the failures of the main grid.
    'secondary_grid_file' :"sg-BSTAR2006-medium.h5",
    'stripped_grid_file' : "sg-Gotberg18.h5",
    'ostar_grid_file' : "sg-OSTAR2002-medium.h5",

    #Setting the wavelengths range:
    'lam_min':3000.,
    'lam_max':7000.,
    #Number of points in the wavelength range:
    'lam_res':2000,
    


    #Performance variables: 
    'cache_limit' : 256


}