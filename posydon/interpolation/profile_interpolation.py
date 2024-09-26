"""Module for performing initial-final profile interpolation
"""

__authors__ = [
    "Elizabeth Teng <elizabethteng@u.northwestern.edu>"
]

import pickle

# POSYDON
from posydon.grids.psygrid import PSyGrid
from posydon.interpolation.IF_interpolation import IFInterpolator
from posydon.utils.posydonwarning import Pwarn

# Math and ML
import os
import random
import numpy as np
import pandas as pd
try:
    import tensorflow as tf
except ImportError:
    raise ImportError('tensorflow is not installed. Please run `pip install .[ml]` in the POSYDON base directory')
tf.get_logger().setLevel('ERROR')
from tensorflow.keras import layers, losses, models, optimizers, backend, utils
from sklearn.decomposition import PCA
from scipy.interpolate import interp1d

class CompileData:

    def __init__(self, train_path, test_path, hms_s2=False,
                 profile_names=['radius','logRho', 'x_mass_fraction_H',
                                'y_mass_fraction_He','z_mass_fraction_metals',
                                'omega','energy']):
        """Extracts profile data from '.h5' grid files and saves to file.
        Args:
            train_path (str) : path/name of '.h5' file for training data.
            test_path (str) : path/name of '.h5' file for testing data.
            hms_s2 (Boolean) : option to get profiles of star 2 in HMS-HMS grid
            profile_names (array-like) : list of profile quantities to extract.
        """
        self.names = profile_names
        
        # extract testing data
        print("extracting testing data")
        test = PSyGrid(test_path)  # load PSyGrid object for testing grid
        self.test_scalars = pd.DataFrame()
        self.test_profiles = []
        testing_failed = []
        
        for i in range(len(test)):
            try:
                scalars,profiles = self.scrape(test,i,hms_s2)
                self.test_scalars = pd.concat([self.test_scalars,
                                              pd.DataFrame(scalars)],ignore_index=True)
                self.test_profiles.append(profiles)
            except:
                testing_failed.append(i)
                pass
                
        if len(testing_failed)>0:
            Pwarn(f"{len(testing_failed)} binaries failed",
                  "IncompletenessWarning")

        # extract training data
        print("extracting training data")
        train = PSyGrid(train_path)  # load PSyGrid object for training grid
        self.scalars = pd.DataFrame()
        self.profiles = []
        training_failed = []
        
        for i in range(len(train)):
            try:
                scalars,profiles = self.scrape(train,i,hms_s2)
                self.scalars = pd.concat([self.scalars,
                                         pd.DataFrame(scalars)],ignore_index=True)
                self.profiles.append(profiles)
            except:
                training_failed.append(i)
                pass       
            
        if 'omega' in self.names:
            self.names.append('norm_omega')
            
        if len(training_failed)>0:
            Pwarn(f"{len(training_failed)} training binaries failed",
                  "IncompletenessWarning")

    def scrape(self,grid,ind,hms_s2):
        """Extracts profile data from one MESA run.
        Args:
            grid (obj) : PSyGrid object.
            ind (int) : index of run to be scraped.
            hms_s2 (Boolean) : option to get profiles of star 2 in HMS-HMS grid
        Returns:
            scalars (array-like) : dictionary containing initial 
                                   m1, m2, p, mass transfer class, 
                                   final star_state, final total_mass.
            profiles (array-like) : all N specified profiles, shape (N,200).
        """
        # open individual run as a DataFrame
        
        if hms_s2==False: # default star 1 profile information
            df = pd.DataFrame(grid[ind]['final_profile1'])
            mass_key = "star_1_mass"
            state_key = "S1_state"
        else: # requesting star 2 information for HMS-HMS grid
            df = pd.DataFrame(grid[ind]['final_profile2'])
            mass_key = "star_2_mass"
            state_key = "S2_state"
            
        df = df.sort_values(by="mass")
        df = df.reset_index(drop=True)

        # grab input values, final star 1 state, final star 1 mass
        total_mass = df["mass"].iloc[-1]
        
        scalars = {"m1":grid.initial_values["star_1_mass"][ind],
                   "m2":grid.initial_values["star_2_mass"][ind],
                   "p":grid.initial_values["period_days"][ind],
                   "MT_class":grid.final_values["interpolation_class"][ind],
                   "star_state":grid.final_values[state_key][ind],
                   "total_mass":total_mass}

        # grab output vectors, interpolate to normalize
        profiles=np.zeros([1+len(self.names),200])
        for i,prof in enumerate(self.names):
            if prof in df.columns:
            
                f = interp1d(df['mass']/total_mass,df[prof],
                             fill_value="extrapolate")
                profile_new = f(np.linspace(0,1,200))
                profiles[i] = profile_new
            else:
                Pwarn(f"{prof} profile not saved in grid, will not be "
                      "included in file", "InappropriateValueWarning")
        
        if 'omega' in self.names:
            profiles[-1]= profiles[self.names.index('omega')]/  \
                          (grid.final_values['S1_surf_avg_omega'][ind]/ \
                           grid.final_values['S1_surf_avg_omega_div_omega_crit'][ind])
        
        if not np.isnan(grid.final_values['S1_surf_avg_omega_div_omega_crit'][ind]):
            return scalars, profiles
        else:
            Pwarn("nan in final values, binary will not be"
                  "included in file","InappropriateValueWarning")

    def save(self, filename):
        """Save extracted profile data.
        Args:
            filename (str) : path/name of '.pkl' file where the data will be saved.
        """
        SAVE_ATTRS = self.__dict__
        DONT_SAVE = []
        myattrs = {key: SAVE_ATTRS[key]
                   for key in SAVE_ATTRS if key not in DONT_SAVE}
        with open(filename, 'wb') as f:
            pickle.dump(myattrs, f)

            
class ProfileInterpolator:
    
    def __init__(self,seed_value):
        """Interfaces with other classes, trains models and predicts profiles.
        Args:
            seed_value (int) : random seed to ensure consistent results
        """
        if seed_value is not None:
            os.environ['PYTHONHASHSEED']=str(seed_value)
            random.seed(seed_value)
            np.random.seed(seed_value)
            tf.random.set_seed(seed_value)
            utils.set_random_seed(seed_value)

    def load_profiles(self,filename,valid_split=0.2):
        """Load and process extracted profile data.
        Args:
            filename (str) : path/name of '.pkl' file to be loaded.
            valid_split (float) : percentage of training data used for validation set
        """
        with open(filename, 'rb') as f:
            myattrs = pd.read_pickle(f)
            for key in myattrs:
                setattr(self, key, myattrs[key])
                
        # processing
        self.test_profiles = np.array(self.test_profiles)
        test_linear_initial = np.transpose([self.test_scalars["m1"],
                                             self.test_scalars["m2"],
                                             self.test_scalars["p"]]) 
        self.test_initial = np.log10(np.array(test_linear_initial))
        
        linear_initial = np.transpose([self.scalars["m1"],
                                       self.scalars["m2"],
                                       self.scalars["p"]])
        
        # random split for training and validation data (default 80/20)
        binaries = np.arange(len(self.profiles))
        np.random.shuffle(binaries)    
        split = int(len(self.profiles)*valid_split) # index at which to split data
    
        self.valid_profiles = np.array(self.profiles)[binaries[:split]]      
        self.valid_initial = np.log10(linear_initial)[binaries[:split]]  
        self.valid_scalars = self.scalars.iloc[binaries[:split]]
        
        self.profiles = np.array(self.profiles)[binaries[split:]]
        self.initial = np.log10(linear_initial)[binaries[split:]]
        self.scalars = self.scalars.iloc[binaries[split:]]

    def train(self,IF_interpolator,density_epochs=1000,density_patience=200,
             comp_bounds_epochs=500,comp_bounds_patience=50,loss_history=False,hms_s2=False,
             depth=12,width=256,depthn=12,widthn=256,lr=0.0001):
        """Trains models for density, H mass fraction, and He mass fraction profile models. 
        Args:
            IF_interpolator (str) : path to '.pkl' file for IF interpolator.
            density_epochs (int) : number of epochs used to train density profile model
            density_patience (int) : patience parameter for NN callback in density profile model
            comp_bounds_epochs (int) : number of epochs used to train composition profiles model
            comp_bounds_patience (int) : patience parameter for NN callback in composition profiles model
            loss_history (Boolean) : option to return training and validation loss histories
            hms_s2 (Boolean) : option to do profiles of star 2 in HMS-HMS grid
            depth (int) : depth of neural network for principal component weights
            width (int) : width of neural network for principal component weights
            depthn (int) : depth of neural network for normalizing value
            widthn (int) : width of neural network for normalizing value
            lr (float) : learning rate for neural network training
        Returns:
            self.comp.loss_history (array-like) : training and validation loss history for composition profiles
            self.dens.loss_history (array-like) : training and validation loss history for density profiles
            
        """        
        # instantiate and train composition (H and He mass fraction) profiles model
        self.comp = Composition(self.initial, 
                                self.profiles[:,self.names.index("x_mass_fraction_H")], 
                                self.profiles[:,self.names.index("y_mass_fraction_He")], 
                                self.scalars["star_state"], 
                                self.valid_initial, 
                                self.valid_profiles[:,self.names.index("x_mass_fraction_H")], 
                                self.valid_profiles[:,self.names.index("y_mass_fraction_He")], 
                                self.valid_scalars["star_state"], 
                                IF_interpolator,
                                comp_bounds_epochs,comp_bounds_patience,hms_s2)

        # instantiate and train density profile model
        self.dens = Density(self.initial,
                            self.profiles[:,self.names.index("logRho")],
                            self.scalars["MT_class"],
                            self.valid_initial,
                            self.valid_profiles[:,self.names.index("logRho")],
                            self.valid_scalars["MT_class"],
                            IF_interpolator,hms_s2=hms_s2,
                            depth=depth, width=width,
                            depthn=depthn, widthn=widthn)
        self.dens.train(prof_epochs=density_epochs,prof_patience=density_patience,lr=lr)
        
        if loss_history==True:
            return self.comp.loss_history, self.dens.loss_history

    def predict(self,inputs):
        """Predict density, H mass fraction, and He mass fraction profiles from inputs.
        Args:
            inputs (array-like) : log-space initial conditions of N binaries to predict, shape (N,3).
        Returns:
            mass_coords (array-like) : linear-scale mass enclosed profile coordinates.
            density_profiles (array-like) : log-scale density profile coordinates.
            h_profiles (array-like) : H mass fraction profile coordinates
            he_profiles (array-like) : He mass fraction profile coordinates
        """
        mass_coords, density_profiles = self.dens.predict(inputs)
        mass_coords, h_profiles, he_profiles = self.comp.predict(inputs)
        
        return mass_coords, density_profiles, h_profiles, he_profiles
        
                    
    def save(self, filename):
        """Save complete profiles interpolation model.
        Args:
            filename (str) : path/name of '.pkl' file where the model will be saved.
        """
        SAVE_ATTRS = self.__dict__
        DONT_SAVE = []
        myattrs = {key: SAVE_ATTRS[key]
                   for key in SAVE_ATTRS if key not in DONT_SAVE}

        with open(filename, 'wb') as f:
            pickle.dump(myattrs, f)

    def load(self, filename):
        """Load interpolation model, which can be used for predictions.
        Args:
            filename (str) : path/name of '.pkl' file to be loaded.
        """
        with open(filename, 'rb') as f:
            myattrs = pickle.load(f)
            for key in myattrs:
                setattr(self, key, myattrs[key])
    
    def mono_decrease(self,profiles):
        """Enforce monotonicity in profiles such as density that must monotonically decrease
        Args:
            profiles (array-like) : N profiles to be post-processed for monotonicity
        Returns:
            profiles_mono (array-like) : post-processed profiles
        """
        
        def mono_renorm(arr):
            arr_copy = arr.copy()
            # starting from surface, force dropping points up
            for i in range(1,len(arr_copy)):
                if arr_copy[-i]>arr_copy[-i-1]:
                    arr_copy[-i-1]=arr_copy[-i]
            # cut off points above center value
            return np.where(arr_copy>arr[0], arr[0], arr_copy)

        profiles_mono=profiles.copy()

        for i in range(len(profiles)):
            if len(np.where(profiles[i][1:]-profiles[i][:-1]>0)[0]>0):
                profiles_mono[i] = mono_renorm(profiles[i])

        return profiles_mono

    
class Density:

    def __init__(self,initial,profiles,mt,valid_initial,
                 valid_profiles,valid_mt,IF_interpolator,n_comp=8, hms_s2=False,
                 depth=12,width=256,depthn=12,widthn=256):
        """Creates and trains density profile model.
        Args:
            initial (array-like) : log-space initial conditions for training data.
            profiles (array-like) : final density profiles for training data. 
            valid_initial (array-like) : log-space initial conditions for validation data.
            valid_profiles (array-like) : final density profiles for validation data.
            IF_interpolator (string) : path to .pkl file for IF interpolator for central density, final mass values
            n_comp (int) : number of PCA components. 
            hms_s2 (Boolean) : option to do profiles of star 2 in HMS-HMS grid
        """
        self.n_comp = n_comp
        self.hms_s2 = hms_s2
        
        # process training data 
        self.initial = initial  # initial conditions in log space
        self.mt = mt
        self.surf_val = profiles[:,-1] 
        self.center_val = profiles[:,0]
        
        # process validation data
        self.valid_initial = valid_initial
        self.valid_mt = valid_mt
        self.valid_surf_val = valid_profiles[:,-1] 
        self.valid_center_val = valid_profiles[:,0]

        # process profiles for modeling
        profiles_norm = (profiles-self.surf_val[:,np.newaxis])\
                            /(self.center_val-self.surf_val)[:,np.newaxis]  # minmax normalized profiles
        valid_profiles_norm = (valid_profiles-self.valid_surf_val[:,np.newaxis])\
                            /(self.valid_center_val-self.valid_surf_val)[:,np.newaxis]
           
        self.pca = PCA(n_components=self.n_comp).fit(profiles_norm) # perform PCA
        self.pca_weights = self.pca.transform(profiles_norm)
        self.valid_pca_weights = self.pca.transform(valid_profiles_norm)
        
        # instantiate model
        self.prof_models = {}
        for mt in self.mt.unique():
            model = models.Sequential()
            model.add(layers.Dense(width,input_dim=3,activation='relu'))
            for i in range(depth-1):
                model.add(layers.Dense(width,input_dim=width,activation='relu'))           
            model.add(layers.Dense(n_comp,input_dim=width,activation=None))
            self.prof_models[mt] = model
                
        self.model_norm = models.Sequential()
        self.model_norm.add(layers.Dense(widthn,input_dim=3,activation="relu"))
        for i in range(depthn-1):
            self.model_norm.add(layers.Dense(widthn,input_dim=widthn,activation='relu'))
        self.model_norm.add(layers.Dense(1,input_dim=widthn,activation=None))
        
        self.model_IF = IFInterpolator()  # instantiate POSYDON initial-final interpolator object
        self.model_IF.load(filename=IF_interpolator)

    def train(self,loss=losses.MeanSquaredError(),prof_epochs=1000,prof_patience=200):
        """Trains NN models. 
        Args: 
            loss (object) : loss function for training.
            prof_epochs (int) : number of epochs used to train neural network
            prof_patience (int) : patience parameter for callback in neural network
        """
        print("training on PCA weights...")
        
        self.loss_history = {}
        for mt in self.mt.unique():
            if mt=="not_converged" or mt=="initial_MT":
                pass
            else:
                print(mt)
                inds = np.where(self.mt==mt)[0]
                valid_inds = np.where(self.valid_mt==mt)[0]
                self.prof_models[mt].compile(optimizers.Adam(clipnorm=1,learning_rate=lr),loss=loss)
                callback = tf.keras.callbacks.EarlyStopping(monitor="loss",patience=prof_patience)
                history = self.prof_models[mt].fit(self.initial[inds],self.pca_weights[inds],
                                              epochs=prof_epochs,callbacks=[callback],verbose=0,
                                              validation_data=(self.valid_initial[valid_inds],
                                                               self.valid_pca_weights[valid_inds]))
                self.loss_history[mt] = np.array([history.history['loss'],history.history['val_loss']])
        
        self.model_norm.compile(optimizers.Adam(clipnorm=1,learning_rate=lr),loss=loss)
        callback = tf.keras.callbacks.EarlyStopping(monitor="loss",patience=40)
        self.model_norm.fit(self.initial,self.surf_val,epochs=300,callbacks=[callback],
                            validation_data = (self.valid_initial,self.valid_surf_val))
                                     
        print("done training")
    
    def predict(self,inputs):
        """Predicts profile for n sets of given inputs, in array of shape (n,3).
        Args: 
            inputs (array-like) : log-space initial conditions of N binaries to predict, shape (N,3).
        Returns:
            mass_coords (array-like) : linear-scale mass enclosed profile coordinates.
            density_profiles (array_like) : log-scale density profile coordinates.
        """
        # predict PCA weights
        regress_prof = lambda x: self.model_prof(x)
        pca_weights_pred = regress_prof(inputs).numpy()
                
        # predict surface density
        regress_rho = lambda x: self.model_rho(x)
        min_rho = regress_rho(inputs).numpy()[:,0]
        
        # IF interpolate final mass, center density 
        if self.hms_s2==False:
            m_ind = self.model_IF.interpolators[0].out_keys.index("star_1_mass")
            center_ind = self.model_IF.interpolators[0].out_keys.index('S1_log_center_Rho')
        else:
            m_ind = self.model_IF.interpolators[0].out_keys.index("star_2_mass")
            center_ind = self.model_IF.interpolators[0].out_keys.index('S2_log_center_Rho')
        max_rho = self.model_IF.interpolators[0].test_interpolator(10**inputs)[:,center_ind]                    
        pred_mass = self.model_IF.interpolators[0].test_interpolator(10**inputs)[:,m_ind]
            
        # reconstruct profile
        norm_prof = self.pca.inverse_transform(pca_weights_pred*self.scaling)
        density_profiles = norm_prof*(max_rho[:,np.newaxis]-min_rho[:,np.newaxis]) \
                           + min_rho[:,np.newaxis]
        
        # construct mass enclosed profile coordinates
        mass_coords = np.linspace(0,1,200)*pred_mass[:,np.newaxis] 
        
        return mass_coords,density_profiles

    
class Composition:
    
    def __init__(self,initial,h_profiles,he_profiles,star_state,
                 valid_initial,valid_h_profiles,valid_he_profiles,valid_star_state,
                 IF_interpolator,training_epochs=500, training_patience=50,hms_s2=False):
        """Creates and trains H mass fraction and He mass fraction profiles model.
        Args:
            initial (array-like) : log-space initial conditions for training data.
            h_profiles (array-like) : final H mass fraction profiles for training data. 
            he_profiles (array-like) : final He mass fraction profiles for training data. 
            star_state (array-like) : final star 1 state for training data.
            valid_initial (array-like) : log-space initial conditions for validation data.
            valid_h_profiles (array-like) : final H mass fraction profiles for validation data.
            valid_he_profiles (array-like) : final He mass fraction profiles for validation data.
            valid_star_state (array-like) : final star 1 state for testing data.
            IF_interpolator (string) : path to .pkl file for IF interpolator
            training_epochs (int) : number of epochs used to train neural networks
            training_patience (int) : patience parameter for callback in neural networks
            hms_s2 (Boolean) : option to do profiles of star 2 in HMS-HMS grid
        """
        self.hms_s2 = hms_s2
        
        self.initial = initial
        self.h_profiles = h_profiles
        self.he_profiles = he_profiles
        self.star_state = star_state
        self.valid_initial = valid_initial
        self.valid_h_profiles = valid_h_profiles
        self.valid_he_profiles = valid_he_profiles
        self.valid_star_state = valid_star_state

        # load IF interpolator
        self.model_IF = IFInterpolator()  # instantiate POSYDON initial-final interpolator object
        self.model_IF.load(filename=IF_interpolator)
        self.interp = self.model_IF.interpolators[0]
        
        if self.hms_s2==False: # default to star 1
            self.c_h_ind = self.interp.out_keys.index("S1_center_h1")
            self.s_h_ind = self.interp.out_keys.index("S1_surface_h1")
            self.c_he_ind = self.interp.out_keys.index("S1_center_he4")
            self.s_he_ind = self.interp.out_keys.index("S1_surface_he4")
        else: # if interpolating star 2 profiles in HMS-HMS grid
            self.c_h_ind = self.interp.out_keys.index("S2_center_h1")
            self.s_h_ind = self.interp.out_keys.index("S2_surface_h1")
            self.c_he_ind = self.interp.out_keys.index("S2_center_he4")
            self.s_he_ind = self.interp.out_keys.index("S2_surface_he4")
            
        # final star states
        self.star_states = ['None','WD','NS','BH',
                             'stripped_He_non_burning',
                             'stripped_He_Core_He_burning',
                             'stripped_He_Core_C_burning',
                             'stripped_He_Central_He_depleted',
                             'stripped_He_Central_C_depletion',
                             'H-rich_non_burning',
                             'H-rich_Central_He_depleted',
                             'H-rich_Central_C_depletion',
                             'H-rich_Shell_H_burning',
                             'H-rich_Core_H_burning',
                             'H-rich_Core_He_burning',
                             'H-rich_Core_C_burning']
        
        # sort training data into classes by final star state
        self.sort_ind = {}
        self.valid_sort_ind = {}
        for name in self.star_states:
            self.sort_ind[name] = np.where(self.star_state==name)[0]
            self.valid_sort_ind[name] = np.where(self.valid_star_state==name)[0]
        
        # create and train models for profile boundaries
        self.bounds_models, self.loss_history = self.learn_bounds(training_epochs,training_patience)
        
    def learn_bounds(self,training_epochs,training_patience):
        """Creates and trains NNs to predict boundary points for each star 1 state.
        Args:
            training_epochs (int) : number of epochs used to train neural networks
            training_patience (int) : patience parameter for callback in neural networks
        Returns:
            b_models (array-like) : dictionary containing boundary models.
            loss_history (array-like) : training and validation loss histories
        """
        b_models = {}
        loss_history = {}
        self.empty = [] # keep track of which star states have no data
        for state in ['stripped_He_Core_He_burning',
                      'stripped_He_Core_C_burning',
                      'stripped_He_Central_He_depleted',
                      'stripped_He_Central_C_depletion',
                      'H-rich_Central_He_depleted',
                      'H-rich_Central_C_depletion',
                      'H-rich_Shell_H_burning',
                      'H-rich_Core_H_burning',
                      'H-rich_Core_He_burning',
                      'H-rich_Core_C_burning']:
            
            # identify input/output training data for class
            indices = self.sort_ind[state]
            inputs = self.initial[indices]
            h_prof = self.h_profiles[indices]
            he_prof = self.he_profiles[indices]

            # identify input/output testing data for class
            valid_indices = self.valid_sort_ind[state]
            valid_inputs = self.valid_initial[valid_indices]
            valid_h_prof = self.valid_h_profiles[valid_indices]
            valid_he_prof = self.valid_he_profiles[valid_indices]
            
            # calculate boundary points for training and testing data 
            if state in ['H-rich_Shell_H_burning',
                         'H-rich_Core_H_burning',
                         'H-rich_Core_He_burning',
                         'H-rich_Core_C_burning']:  # profiles with 2 boundary points
                outs = 2
                bounds = []   # collect information on parameters for given H and He profile shapes
                nonflat = []  # ensures that training data only has the correct "non-flat" shape
                              # to avoid issues with calculating the boundary points
                valid_bounds = []
                valid_nonflat = []
                # calculate first and last points in each H profile with large increases
                # (He profiles have the same boundary points)
                for i in range(len(indices)):
                    try:
                        diff = np.where(h_prof[i][1:]-h_prof[i][:-1]>0.002)[0]
                        bounds.append(np.array([diff[0]+1,diff[-1]+1])/200)
                        nonflat.append(i)
                    except:
                        bounds.append([np.nan,np.nan])
                for i in range(len(valid_indices)):
                    try:
                        diff = np.where(valid_h_prof[i][1:]-valid_h_prof[i][:-1]>0.002)[0]
                        valid_bounds.append(np.array([diff[0]+1,diff[-1]+1])/200)
                        valid_nonflat.append(i)
                    except:
                        valid_bounds.append([np.nan,np.nan])
                        
            elif state in ['H-rich_Central_He_depleted',
                           'H-rich_Central_C_depletion']: # H profiles with 1 boundary point
                outs = 3
                # calculate the point in each H profile with the largest increase
                hbound = (np.argmax(h_prof[:,1:]-h_prof[:,:-1],axis=1)+1)/200
                valid_hbound = (np.argmax(valid_h_prof[:,1:]-valid_h_prof[:,:-1],axis=1)+1)/200
                max_He = np.max(he_prof,axis=1) # maximum He mass fraction
                valid_max_He = np.max(valid_he_prof,axis=1)
                dep_He = (np.argmax(he_prof[:,1:]-he_prof[:,:-1],axis=1)+1)/200 # location of Helium depletion zone
                valid_dep_He = (np.argmax(valid_he_prof[:,1:]-valid_he_prof[:,:-1],axis=1)+1)/200
                bounds = np.transpose([hbound,dep_He,max_He])
                valid_bounds = np.transpose([valid_hbound,valid_dep_He,valid_max_He])
            
            elif state in ['stripped_He_Core_He_burning',
                         'stripped_He_Core_C_burning',
                         'stripped_He_Central_He_depleted',
                         'stripped_He_Central_C_depletion']:
                outs = 2
                bounds = []
                nonflat = []  # ensures that training data only has the correct "non-flat" shape
                              # to avoid issues with calculating the boundary points
                valid_bounds = []
                valid_nonflat = []
                # calculate first and points in each profile with large increases
                for i in range(len(indices)):
                    try:
                        diff = np.where(he_prof[i][1:]-he_prof[i][:-1]>0.002)[0]
                        bounds.append(np.array([diff[0]+1,diff[-1]+1])/200)
                        nonflat.append(i)
                    except:
                        bounds.append([np.nan,np.nan])
                for i in range(len(valid_indices)):
                    try:
                        diff = np.where(valid_he_prof[i][1:]-valid_he_prof[i][:-1]>0.002)[0]
                        valid_bounds.append(np.array([diff[0]+1,diff[-1]+1])/200)
                        valid_nonflat.append(i)
                    except:
                        valid_bounds.append([np.nan,np.nan])
                
            # instantiate and train model on bounds
            model = models.Sequential([
                    layers.Dense(10,input_dim=3,activation=None),
                    layers.Dense(10,input_dim=10,activation='relu'),
                    layers.Dense(10,input_dim=10,activation='relu'),
                    layers.Dense(10,input_dim=10,activation='tanh'),
                    layers.Dense(10,input_dim=10,activation='tanh'),
                    layers.Dense(10,input_dim=10,activation='tanh'),
                    layers.Dense(outs,input_dim=10,activation="sigmoid")])

            model.compile(optimizers.Adam(clipnorm=1),loss=losses.MeanSquaredError())
            callback = tf.keras.callbacks.EarlyStopping(monitor='loss', patience=training_patience)
            
            if len(indices)==0:
                Pwarn(f"no training data available for {state}",
                      "InappropriateValueWarning")
                loss_history[state]=np.nan
                self.empty.append(state)
            
            elif state in ['H-rich_Central_He_depleted',
                         'H-rich_Central_C_depletion']:
                history = model.fit(inputs,np.array(bounds),
                                    epochs=training_epochs,verbose=0,callbacks=[callback],
                                    validation_data=(valid_inputs,
                                                     np.array(valid_bounds)))
                loss_history[state] = np.array([history.history['loss'],history.history['val_loss']])

            else:
                history = model.fit(inputs[nonflat],np.array(bounds)[nonflat],
                                    epochs=training_epochs,verbose=0,callbacks=[callback],
                                    validation_data=(valid_inputs[valid_nonflat],
                                                     np.array(valid_bounds)[valid_nonflat]))
                loss_history[state] = np.array([history.history['loss'],history.history['val_loss']])
                        
            b_models[state] = model
            print(f"finished {state}")
        return b_models, loss_history
            
    def predict_single(self,initial,center_H,surface_H,center_He,surface_He,star_state):
        """Predict a profile for a single binary
        Args:
            initial (array-like) : initial position of binary
            center_H (float) : center H mass fraction 
            surface_H (float) : surface H mass fraction
            center_He (float) : center He mass fraction 
            surface_He (float) : surface He mass fraction
            star_state (str) : final star state
        Returns:
            H (array-like) : predicted H mass fraction profile 
            He (array-like) : predicted He mass fraction profile
        """
        if "stripped_He" in star_state:
            H = np.zeros(200) # stripped Helium stars have no Hydrogen
        
        if star_state == "stripped_He_non_burning":
            He = np.ones(200)*surface_He # non-burning stars have a flat He profile
        
        if star_state in ['stripped_He_Core_He_burning',
                 'stripped_He_Core_C_burning',
                 'stripped_He_Central_He_depleted',
                 'stripped_He_Central_C_depletion']:
            # predicting profile shape parameters 
            b = self.bounds_models[star_state](tf.convert_to_tensor([initial])).numpy()[0]                        
            He = np.ones(200) * center_He
            He[int(b[1]*200):] = surface_He
            f_He = interp1d(b,[center_He,surface_He],fill_value="extrapolate")
            # use f to construct the profile points in the shell burning region
            He[int(b[0]*200):int(b[1]*200)] = f_He(np.linspace(0,1,200)[int(b[0]*200):int(b[1]*200)])
            
        if star_state in ["H-rich_non_burning","None"]: # these states have flat H and He profiles
            H = np.ones(200)*surface_H
            He = np.ones(200)*surface_He
            
        if star_state in ["WD","NS","BH"]: # profiles for compact objects are arbitrary -- return nans
            H = np.ones(200)*np.nan
            He = np.ones(200)*np.nan
        
        # construct step-shaped profile - H and He profiles have symmetrical shapes
        if star_state in ["H-rich_Central_He_depleted",
                  "H-rich_Central_C_depletion"]:
            # predicting profile shape parameters 
            b = self.bounds_models[star_state](tf.convert_to_tensor([initial])).numpy()[0]            
            H = np.ones(200)*center_H
            H[int(b[0]*200):] = surface_H
            He = np.ones(200) * center_He
            He[int(b[1]*200):] = b[2]
            He[int(b[0]*200):] = surface_He
        
        # construct shell-shaped profile - H and He profiles have symmetrical shapes
        if star_state in ["H-rich_Shell_H_burning",
                  "H-rich_Core_H_burning",
                  "H-rich_Core_He_burning",
                  "H-rich_Core_C_burning"]:
            # predicting profile shape parameters 
            b = self.bounds_models[star_state](tf.convert_to_tensor([initial])).numpy()[0]                        
            new = np.ones([2,200]) * np.array([center_H,center_He])[:,np.newaxis]
            new[:,int(b[1]*200):] = np.array([surface_H,surface_He])[:,np.newaxis]
            f_H = interp1d(b,[center_H,surface_H],fill_value="extrapolate")
            f_He = interp1d(b,[center_He,surface_He],fill_value="extrapolate")
            # use f to construct the profile points in the shell burning region
            new[0][int(b[0]*200):int(b[1]*200)] = f_H(np.linspace(0,1,200)[int(b[0]*200):int(b[1]*200)])
            new[1][int(b[0]*200):int(b[1]*200)] = f_He(np.linspace(0,1,200)[int(b[0]*200):int(b[1]*200)])
            H = new[0]
            He = new[1]
        
        # if there is no training data in the state, return nans
        if star_state in self.empty:
            H = np.ones(200)*np.nan
            He = np.ones(200)*np.nan
            
        return H, He
        
    def predict(self,inputs):
        """Predict H mass fraction profiles from inputs.
        Args:
            inputs (array-like) : log-space initial conditions of N binaries to predict, shape (N,3).
        Returns:
            mass_coords (array-like) : linear-scale mass enclosed profile coordinates.
            h_profiles (array_like) : H mass fraction profile coordinates.
            he_profiles (array_like) : He mass fraction profile coordinates.
        """
        # IF interpolate H mass fraction values at center, surface; final star state
        center_h_vals = self.interp.test_interpolator(10**inputs)[:,self.c_h_ind]
        surface_h_vals = self.interp.test_interpolator(10**inputs)[:,self.s_h_ind]
        center_he_vals = self.interp.test_interpolator(10**inputs)[:,self.c_he_ind]
        surface_he_vals = self.interp.test_interpolator(10**inputs)[:,self.s_he_ind]
        
        # IF interpolate final masses, final star states 
        if self.hms_s2==False:
            m_ind = self.model_IF.interpolators[0].out_keys.index("star_1_mass")
            star_state_vals = self.interp.test_classifiers(10**inputs)['S1_state']   
        else:
            m_ind = self.model_IF.interpolators[0].out_keys.index("star_2_mass")
            star_state_vals = self.interp.test_classifiers(10**inputs)['S2_state']   
        pred_mass = self.interp.test_interpolator(10**inputs)[:,m_ind]
        
        # generate predicted profiles
        pred_profiles = []
        for i in range(len(inputs)):
            pred_H,pred_He = self.predict_single(inputs[i],center_h_vals[i],surface_h_vals[i],
                                                 center_he_vals[i],surface_he_vals[i],star_state_vals[i])
            pred_profiles.append([pred_H,pred_He])

        # generate mass enclosed profile coordinates
        mass_coords = np.linspace(0,1,200)*pred_mass[:,np.newaxis] 
        h_profiles = np.array(pred_profiles)[:,0]
        he_profiles = np.array(pred_profiles)[:,1]
            
        return mass_coords, h_profiles, he_profiles
