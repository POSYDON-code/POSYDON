from posydon.spectral_synthesis.spectral_tools import population_data
from posydon.spectral_synthesis.spectral_grids import spectral_grids
import copy
import numpy as np
import astropy.constants as con
import astropy.units as unt

grid_keys = [
    'main_grid',
    'secondary_grid',
    'ostar_grid',
    'stripped_grid',
]
Zo = 0.0142
class population_cmd():
    
    def __init__(self,file,**kwargs):
        ######Initializing the parameters ##################
        for key, value in kwargs.items():
            setattr(self, key, value)

        #population_file = kwargs.get('population_file')
        population_file = file
        time = kwargs.get('time')
        self.failed_stars = 0 #int, the stars that failed durin the spectra process 
        self.missing_stars = 0 #The stars that are missing due to POSYDON, some detached cases for now. 
        #Creating lists/numpy arrays for investigation for the failled and succesful stars. 
        self.stars_fails = []
        self.stars_run   = []
        self.stars_grid_fail = []
        self.population = population_data(population_file,time)
        #Creating readable arrays for the stars objects. 
        self.total_binaries = len(self.population)

        #Initializing the spectral_grids object and parameters used.
        #To do put an option for changing the wavelength 
        self.grids = spectral_grids()
        self.scaling_factor = kwargs.get('scaling_factor')
        self.filters = self.grids.filters
        self.photgrids = self.grids.photgrids
    #Making a V-(B-V) diagram! Working on a similar logic that has been used for the spectra!

    
    
    def calc_colours(self,star):
   
        scale = self.scaling_factor
        F_obs = {}

        if "stripped" in star.state:
            Z = star.metallicity*Zo
            M = star.mass/con.M_sun
            M_min =  self.grids.spectral_grids['stripped_grid'].axis_x_min['M_init']
            M_max = self.grids.spectral_grids['stripped_grid'].axis_x_max['M_init']
            if M < M_min or M > M_max:
                self.failed_stars +=1
                return None
            x = {'M_init':M,'Z': Z}
            F_obs = self.grids.photogrid_flux('stripped_grid',1*unt.m/unt.m,**x)
            return F_obs
    
        Fe_H = star.Fe_H
        Z_Zo = star.metallicity
        Teff = copy.copy(star.get_Teff(self.grids.T_max,self.grids.T_min))
        logg = copy.copy(star.get_logg(self.grids.logg_max,self.grids.logg_min))
        x = {'Teff':Teff ,'log(g)': logg,'[Fe/H]': Fe_H,'Z/Zo':Z_Zo}

        if Teff is None or logg is None: 
            return None 

        #For warmer stars we go to the ostar regime:
        if Teff > 30000:
            try:

                F_obs = self.grids.photogrid_flux('ostar_grid',star.R**2*scale**-2,**x)
                return F_obs
            except:
                return None
        try: 
            F_obs= self.grids.photogrid_flux('main_grid',star.R**2*scale**-2,**x)
            return F_obs
        except LookupError:
            if Teff > 20000:
                logg = max(logg, 4.0)
            elif Teff > 12000:
                logg = max(logg, 3.0)
            elif Teff > 8000:
                logg = max(logg, 2.0)
            elif Teff > 6000:
                logg = max(logg, 1.0)
            
            F_obs = self.grids.photogrid_flux('main_grid',star.R**2*scale**-2,**x)
            return F_obs
  
            
    def colour_mag(self,star):
        mags = {}
        F_obs = self.calc_colours(star)

        if F_obs is None:
            return None 

        for filter in F_obs:
            mags[filter] = -2.5*np.log10(F_obs[filter])
        return mags
  
        


    def population_mag(self,num_binaries = None ):
        if num_binaries is None:
            num_binaries = len(self.population)
        elif num_binaries > len(self.population):
            raise Exception('The number of binaries exceeds the number of binaries in the population!')
        V = []
        B_V = []
        L = []
        population = self.population\
        
        for i in range(num_binaries):
            binary = population[i]
            for j in (0,1):
                newstar = binary[j]
                magnitude = self.colour_mag(newstar)
                if magnitude is not None:
                    V.append(magnitude['V'].value)
                    B_V.append(magnitude['B'].value-magnitude['V'].value)
                    L.append(newstar.L/con.L_sun)
        return B_V,V,L

    def sub_population_mag(self,number_of_binaries,state):
        V = []
        B_V = []
        L = []
        stars = []
        for binary_number in range(number_of_binaries):
            for i in (1,2):
                newstar = star(binary_number,i)
                #if state in newstar.state:
                if newstar.binary_state in state:
                    magnitude = self.colour_mag(newstar)
                    if magnitude is not None:
                        stars.append(newstar)
                        V.append(magnitude['V'].value)
                        B_V.append(magnitude['B'].value-magnitude['V'].value)
                        L.append(newstar.L/con.L_sun)
        return B_V,L,V,stars











    """
    def create_spectrum_single(self,star,**kwargs):
        scale = self.scaling_factor
        if "stripped" in star.state:
            M = star.mass/con.M_sun
            M_min = self.grids.specgrid_stripped.axis_x_min['M_init']
            M_max = self.grids.specgrid_stripped.axis_x_max['M_init']
            if M < M_min or M > M_max:
                self.failed_stars +=1
                return None
            return self.grids.stripped_grid_flux(star)
        
        Fe_H = star.Fe_H
        Z_Zo = star.metallicity
        Teff = copy.copy(star.get_Teff(self.grids.T_max,self.grids.T_min))
        logg = copy.copy(star.get_logg(self.grids.logg_max,self.grids.logg_min))
        x = {'Teff':Teff ,'log(g)': logg,'[Fe/H]': Fe_H,'Z/Zo':Z_Zo}

        if Teff == None or logg == None:
            self.failed_stars +=1 
            return None
        
        #For temperature higher than then lo_limits of Teff in main grids we calculate the spectra using the Ostar girds. 
        if Teff >= 30000:
            #Setting the acceptable boundaries for the ostar grid in the logg.
            logg_min = self.specgrid_ostar.axis_x_min['log(g)']
            logg_max = self.specgrid_ostar.axis_x_max['log(g)']
            if logg > logg_min and logg<logg_max: 
                    try:
                        Flux = self.ostar_grid_flux(**x)
                        return Flux*star.R**2*scale**-2
                    except LookupError:
                        self.failed_stars +=1
                        return None
            else:
                self.failed_stars +=1
                return None
            
        try:
            #normal_start = datetime.datetime.now()
            Flux = self.main_grid_flux(**x) 
            #normal_end = datetime.datetime.now()
            #print( 'The spectral normal time is: {time}'.format(time = normal_end - normal_start ))            
            return Flux*star.R**2*scale**-2
        except LookupError:
            try:
                # lo limits for the secondary grid: 
                #logg_min = grids.specgrid_secondary.axis_x_min['log(g)']
                #Teff_min = grids.specgrid_secondary.axis_x_min['Teff']
                #Calculating all the exeption that the first grid gave with the secondary. 
                #if Teff>15000.0 and logg > logg_min:
                F = self.secondary_grid_flux(**x)
                return F*star.R**2*scale**-2
            except LookupError:
                # if and else statements that fix the grid voids.  
                if Teff > 20000:
                    logg = max(logg, 4.0)
                elif Teff > 12000:
                    logg = max(logg, 3.0)
                elif Teff > 8000:
                    logg = max(logg, 2.0)
                elif Teff > 6000:
                    logg = max(logg, 1.0)
                try:
                    F = self.main_grid_flux(**x)
                    return F*star.R**2*scale**-2 
                except LookupError:
                    self.failed_stars +=1
                    return None



    def add_spectra(self,spectrum_1,spectrum_2):
        if spectrum_1 is not None and spectrum_2 is not None:
            return spectrum_1 + spectrum_2
        elif spectrum_1 is None:
            return spectrum_2
        elif spectrum_2 is None:
            return spectrum_1
        else:
            raise Exception("Spectrum 1 and Spectrum 2 have not allowed values.")
    
    def create_spectrum_binary(self,binary):
    # check if binary has two or any non-degenerate stars
        star1 = binary[0]
        star2 = binary[1]
        if star1.state not in ['NS','BH','massless_remnant']:
            spectrum_1 = self.create_spectrum_single(star1)
        else:
            spectrum_1 = None

        if star2.state not in ['NS','BH','massless_remnant']:
            spectrum_2 = self.create_spectrum_single(star2)
        else: 
            spectrum_2 = None
        return self.add_spectra(spectrum_1, spectrum_2)


    def create_spectrum_population(self,num_binaries = None ):
        if num_binaries is None:
            num_binaries = len(self.population)
        elif num_binaries > len(self.population):
            raise Exception('The number of binaries exceeds the number of binaries in the population!')
        
        create_spectrum_population = np.zeros(len(self.grids.lam_c))*unt.m/unt.m*0
        failed_binaries = 0
        population = self.population
        for i in range(num_binaries):
            binary = population[i]
            spec_start = datetime.datetime.now()
            spectrum = self.create_spectrum_binary(binary)
            spec_end = datetime.datetime.now()
           # print('The {i} binary time was: {time}'.format(i=i,time = spec_end - spec_start ))
            if spectrum is not None:
                create_spectrum_population += spectrum
        return create_spectrum_population , self.grids.lam_c
       

    def create_spectrum_subpopulation(self,state,num_binaries = None):
        if num_binaries is None:
            num_binaries = len(self.population)
        elif num_binaries > len(self.population):
            raise Exception('The number of binaries exceeds the number of binaries in the population!')
        failed_binaries = 0
        create_spectrum_population = np.zeros(len(self.grids.lam_c))*unt.m/unt.m*0
        population = self.population

        for i in range(num_binaries):
            binary = population[i]
            star1 = binary[0]
            if star1.binary_state == state:
                spectrum = self.create_spectrum_binary(binary)
                if spectrum is not None:
                    self.create_spectrum_population += spectrum
        return create_spectrum_population , self.grids.lam_c
"""