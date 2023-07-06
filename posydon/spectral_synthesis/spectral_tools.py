#This is tools and useful funtion that are using in the spectral synthesis and the creation of CMD diagrams combined. 
import numpy as np
import astropy.constants as con
import astropy.units as unt

def array_data(N,last_binary_line,first_binary_line):
    total_binaries = N
    mass_arr = np.zeros((total_binaries,2))
    radius_arr = np.zeros((total_binaries,2))
    L_arr = np.zeros((total_binaries,2))
    state_arr = np.zeros((total_binaries,2),dtype=object)

    for i in range(N):
        try:
            mass_arr[i][:] = np.array([ last_binary_line.S1_mass[i], last_binary_line.S2_mass[i]])
            radius_arr[i][:] = np.array([ 10**last_binary_line.S1_log_R[i], 10**last_binary_line.S2_log_R[i]])
            L_arr[i][:] = np.array([ 10**last_binary_line.S1_log_L[i], 10**last_binary_line.S2_log_L[i]])
            state_arr[i][:] = np.array([ str(last_binary_line.S1_state[i]), str(last_binary_line.S2_state[i])])
        except:
            mass_arr[i][:] = np.array(None,None)
            radius_arr[i][:] = np.array(None,None)
            L_arr[i][:] =  np.array(None,None)
            state_arr[i][:] = np.array(["Failed","Failed"])
        finally:
            if "stripped" in state_arr[i][0]:
                mass_arr[i][0] = first_binary_line.S1_mass[i]
            elif "stripped" in state_arr[i][1]:
                mass_arr[i][1] = first_binary_line.S2_mass[i]
    return mass_arr,radius_arr,L_arr,state_arr

class star():
    
    def __init__(self,binary_number,number,mass,state,R,L,metallicity):
        self.binary_number = binary_number
        self.number = number
        self.mass = mass*con.M_sun
        self.state = state
        self.R = R*con.R_sun
        self.L = L*con.L_sun
        self.logg = None
        self.Teff = None
        self.metallicity = metallicity
        self.Fe_H = np.log10(metallicity)
        
        
    def get_logg(self,max_logg,min_logg): 
        self.logg = np.log10(con.G*self.mass/self.R**2/(unt.cm/unt.s**2))
        if self.logg < max_logg and self.logg > min_logg:
            return self.logg
        else:
            return None
    
    def get_Teff(self,max_Teff,min_Teff):
        self.Teff = (self.L/(4*np.pi*self.R**2*con.sigma_sb))**0.25/unt.K
        if self.Teff < max_Teff and self.Teff > min_Teff:
            return self.Teff
        else:
            return None
        
    def set_metallicity(self,new_Fe_H,max_Fe_H,min_Fe_H):
        if new_Fe_H < max_Fe_H and new_Fe_H > min_Fe_H:
            self.Fe_H = new_Fe_H
        else: 
            raise TypeError
    #This is for the stars that we can't calculate their spectra so we change their logg and thus we need a new radius for them and we set and return. 
    def new_radius(self,new_logg):
        R = np.sqrt(con.G*self.mass/10**(new_logg)/unt.cm *unt.s**2).decompose()
        return R