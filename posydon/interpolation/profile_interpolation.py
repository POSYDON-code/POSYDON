"""Module for performing final profile interpolation
"""

__authors__ = [
    "Elizabeth Teng <elizabethteng@u.northwestern.edu>"
]

import pickle
import warnings

# POSYDON
from posydon.grids.psygrid import PSyGrid
from posydon.interpolation.IF_interpolation import IFInterpolator

# Math and ML
import numpy as np
import pandas as pd
import tensorflow as tf
tf.get_logger().setLevel('ERROR')
from tensorflow.keras import layers, losses, models, optimizers
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
                self.test_scalars = self.test_scalars.append(scalars,ignore_index=True)
                self.test_profiles.append(profiles)
            except:
                testing_failed.append(i)
                pass
        warnings.warn(f"{len(testing_failed)} binaries failed")

        # extract training data
        print("extracting training data")
        train = PSyGrid(train_path)  # load PSyGrid object for training grid
        self.scalars = pd.DataFrame()
        self.profiles = []
        training_failed = []
        
        for i in range(len(train)):
            try:
                scalars,profiles = self.scrape(train,i,hms_s2)
                self.scalars = self.scalars.append(scalars,ignore_index=True)
                self.profiles.append(profiles)
            except:
                training_failed.append(i)
                pass
        warnings.warn(f"{len(training_failed)} training binaries failed")
        

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
        total_mass = grid.final_values[mass_key][ind]
        scalars = {"m1":grid.initial_values["star_1_mass"][ind],
                   "m2":grid.initial_values["star_2_mass"][ind],
                   "p":grid.initial_values["period_days"][ind],
                   "MT_class":grid.final_values["interpolation_class"][ind],
                   "star_state":grid.final_values[state_key][ind],
                   "total_mass":total_mass}

        # grab output vectors, interpolate to normalize
        profiles=np.zeros([len(self.names),200])
        for i,prof in enumerate(self.names):
            if prof in df.columns:
            
                f = interp1d(df['mass']/total_mass,df[prof],
                             fill_value="extrapolate")
                profile_new = f(np.linspace(0,1,200))
                profiles[i] = profile_new
            else:
                warnings.warn(f"{prof} profile not saved in grid, will not be included in file")

        return scalars, profiles

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
    
    def __init__(self):
        """Interfaces with other classes, trains models and predicts profiles.
        """

    def load_profiles(self,filename):
        """Load extracted profile data.
        Args:
            filename (str) : path/name of '.pkl' file to be loaded.
        """
        with open(filename, 'rb') as f:
            myattrs = pd.read_pickle(f)
            for key in myattrs:
                setattr(self, key, myattrs[key])
        self.profiles = np.array(self.profiles)
        self.test_profiles = np.array(self.test_profiles)
    
    def train(self,IF_interpolator,density_epochs=3000,density_patience=200,
             comp_bounds_epochs=500,comp_bounds_patience=50,loss_history=False):
        """Trains models for density, H mass fraction, and He mass fraction profile models. 
        Args:
            IF_interpolator (str) : path to '.pkl' file for IF interpolator.
            density_epochs (int) : number of epochs used to train density profile model
            density_patience (int) : patience parameter for NN callback in density profile model
            comp_bounds_epochs (int) : number of epochs used to train composition profiles model
            comp_bounds_patience (int) : patience parameter for NN callback in composition profiles model
            loss_history (Boolean) : option to return training and validation loss histories
        Returns:
            self.comp.loss_history (array-like) : training and validation loss history for composition profiles
            self.dens.loss_history (array-like) : training and validation loss history for density profiles
            
        """
        # processing
        linear_initial = np.transpose([
            self.scalars["m1"],self.scalars["m2"],self.scalars["p"]])
        initial = np.log10(np.array(linear_initial))
        final_m1 = self.scalars["final_m1"].astype(np.float64)

        test_linear_initial = np.transpose([self.test_scalars["m1"],
                                             self.test_scalars["m2"],
                                             self.test_scalars["p"]]) 
        test_initial = np.log10(np.array(test_linear_initial))
        test_total_mass = self.test_scalars["total_mass"].astype(np.float64)
        
        
        valid_initial = test_intial  #TODO: update these
        valid_scalars = self.test_scalars
        valid_profiles = self.test_profiles
        
        
        # instantiate and train H mass fraction profile model
        h_ind = self.names.index("x_mass_fraction_H")
        he_ind = self.names.index("y_mass_fraction_He")
        self.comp = Composition(initial, self.profiles[:,h_ind], self.profiles[:,he_ind], self.scalars["s1_state"], 
                     valid_initial, valid_profiles[:,h_ind], valid_profiles[:,he_ind], 
                     valid_scalars["star_state"], IF_interpolator,
                     comp_bounds_epochs,comp_bounds_patience)

        # instantiate and train density profile model
        dens_ind = self.names.index("logRho")
        self.dens = Density(initial,
                       self.profiles[:,dens_ind],
                       valid_initial,
                       valid_profiles[:,dens_ind],
                       IF_interpolator)
        self.dens.train(prof_epochs=density_epochs,prof_patience=density_patience)
        
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
        """Load interpolation model, which can only be used for predictions.
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
            # force rising points down
            for i in range(1,len(arr_copy)):
                if arr_copy[i]>arr_copy[i-1]:
                    arr_copy[i]=arr_copy[i-1]
            # cut off points below surface value
            return np.where(arr_copy<arr[-1], arr[-1], arr_copy)
        
        profiles_copy=profiles.copy()
        
        for i in range(len(profiles)):
            if len(np.where(profiles[i][1:]-profiles[i][:-1]>0)[0]>0):
                profiles_copy[i] = mono_renorm(profiles[i])
                
        return profiles_copy
    
    
class Density:

    def __init__(self,initial,profiles,valid_initial,
                 valid_profiles,IF_interpolator,n_comp=8):
        """Creates and trains density profile model.
        Args:
            initial (array-like) : log-space initial conditions for training data.
            profiles (array-like) : final density profiles for training data. 
            valid_initial (array-like) : log-space initial conditions for validation data.
            valid_profiles (array-like) : final density profiles for validation data.
            IF_interpolator (string) : path to .pkl file for IF interpolator for central density, final mass values
            n_comp (int) : number of PCA components. 
        """
        self.n_comp = n_comp
        
        # process training data 
        self.initial = initial  # initial conditions in log space
        self.rho_min = np.min(profiles,axis=1) 
        rho_max = np.max(profiles,axis=1)
        profiles_norm = (profiles-self.rho_min[:,np.newaxis])\
                            /(rho_max-self.rho_min)[:,np.newaxis]  # minmax normalized profiles
        self.pca = PCA(n_components=self.n_comp).fit(profiles_norm) # perform principal component analysis
        pca_weights_unscaled = self.pca.transform(profiles_norm)
        self.scaling = np.std(pca_weights_unscaled,axis=0)
        self.pca_weights = pca_weights_unscaled/self.scaling  # scaled PCA weights
        
        # process testing data
        self.valid_initial = valid_initial
        self.valid_rho_min = np.min(valid_profiles,axis=1)
        valid_rho_max = np.max(valid_profiles,axis=1)
        valid_profiles_norm = (valid_profiles-self.valid_rho_min[:,np.newaxis])\
                            /(valid_rho_max-self.valid_rho_min)[:,np.newaxis]
        valid_pca_weights_unscaled = self.pca.transform(valid_profiles_norm)
        self.valid_pca_weights = valid_pca_weights_unscaled/self.scaling
        
        # instantiate models
        self.model_prof = models.Sequential([
            layers.Dense(15,input_dim=3,activation=None),
            layers.Dense(15,input_dim=15,activation="relu"),
            layers.Dense(15,input_dim=15,activation="relu"),
            layers.Dense(15,input_dim=15,activation="tanh"),
            layers.Dense(15,input_dim=15,activation="tanh"),
            layers.Dense(10,input_dim=10,activation="tanh"),
            layers.Dense(10,input_dim=10,activation="tanh"),
            layers.Dense(n_comp,input_dim=10,activation=None)])
                
        self.model_rho = models.Sequential([
            layers.Dense(10,input_dim=3,activation='relu'),
            layers.Dense(10,input_dim=10,activation='relu'),
            layers.Dense(10,input_dim=10,activation='tanh'),
            layers.Dense(10,input_dim=10,activation='tanh'),
            layers.Dense(10,input_dim=10,activation='tanh'),
            layers.Dense(10,input_dim=10,activation='tanh'),
            layers.Dense(1,activation=None)])
        
        self.model_IF = IFInterpolator()  # instantiate POSYDON initial-final interpolator object
        self.model_IF.load(filename=IF_interpolator)
                
    def train(self,loss=losses.MeanSquaredError(),prof_epochs=3000,prof_patience=200):
        """Trains NN models. 
        Args: 
            loss (object) : loss function for training.
            prof_epochs (int) : number of epochs used to train neural network
            prof_patience (int) : patience parameter for callback in neural network
        """
        print("training on PCA weights...")
        
        self.model_prof.compile(optimizers.Adam(clipnorm=1),loss=loss)
        callback = tf.keras.callbacks.EarlyStopping(monitor="loss",patience=prof_patience)
        history = self.model_prof.fit(self.initial,self.pca_weights,
                                      epochs=prof_epochs,callbacks=[callback],verbose=0,
                                      validation_data=(self.valid_initial,
                                                       self.valid_pca_weights))
        
        self.model_rho.compile(optimizers.Adam(clipnorm=1),loss=loss)
        callback = tf.keras.callbacks.EarlyStopping(monitor="loss",patience=40)
        history = self.model_rho.fit(self.initial,self.rho_min,
                                     epochs=500, callbacks=[callback], verbose=0,
                                     validation_data=(self.valid_initial,
                                                      self.valid_rho_min))
        
        self.loss_history = np.array([history.history['loss'],history.history['val_loss']])
                                     
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
        
        # IF interpolate center density
        center_ind = self.model_IF.interpolators[0].out_keys.index('S1_log_center_Rho')
        max_rho = self.model_IF.interpolators[0].test_interpolator(10**inputs)[:,center_ind]                    
             
        # reconstruct profile
        norm_prof = self.pca.inverse_transform(pca_weights_pred*self.scaling)
        density_profiles = norm_prof*(max_rho[:,np.newaxis]-min_rho[:,np.newaxis]) \
                           + min_rho[:,np.newaxis]
        
        # IF interpolate final mass, construct mass enclosed profile coordinates
        m1_ind = self.model_IF.interpolators[0].out_keys.index("star_1_mass")
        pred_mass = self.model_IF.interpolators[0].test_interpolator(10**inputs)[:,m1_ind]
        mass_coords = np.linspace(0,1,200)*pred_mass[:,np.newaxis] 
        
        return mass_coords,density_profiles

    
class Composition:
    
    def __init__(self,initial,h_profiles,he_profiles,s1,
                 valid_initial,valid_h_profiles,valid_he_profiles,valid_s1,
                 IF_interpolator,training_epochs=500, training_patience=50):
        """Creates and trains H mass fraction and He mass fraction profiles model.
        Args:
            initial (array-like) : log-space initial conditions for training data.
            h_profiles (array-like) : final H mass fraction profiles for training data. 
            he_profiles (array-like) : final He mass fraction profiles for training data. 
            s1 (array-like) : final star 1 state for training data.
            valid_initial (array-like) : log-space initial conditions for validation data.
            valid_h_profiles (array-like) : final H mass fraction profiles for validation data.
            valid_he_profiles (array-like) : final He mass fraction profiles for validation data.
            valid_s1 (array-like) : final star 1 state for testing data.
            IF_interpolator (string) : path to .pkl file for IF interpolator
            training_epochs (int) : number of epochs used to train neural networks
            training_patience (int) : patience parameter for callback in neural networks

        """
        self.initial = initial
        self.h_profiles = h_profiles
        self.he_profiles = he_profiles
        self.s1 = s1
        self.valid_initial = valid_initial
        self.valid_h_profiles = valid_h_profiles
        self.valid_he_profiles = valid_he_profiles
        self.valid_s1 = valid_s1

        # load IF interpolator
        self.model_IF = IFInterpolator()  # instantiate POSYDON initial-final interpolator object
        self.model_IF.load(filename=IF_interpolator)
        self.c_h_ind = self.model_IF.interpolators[0].out_keys.index("S1_center_h1")
        self.s_h_ind = self.model_IF.interpolators[0].out_keys.index("S1_surface_h1")
        self.c_he_ind = self.model_IF.interpolators[0].out_keys.index("S1_center_he4")
        self.s_he_ind = self.model_IF.interpolators[0].out_keys.index("S1_surface_he4")
        self.interp = self.model_IF.interpolators[0]
        
        # star 1 states
        self.class_names = ['None','WD','NS','BH',
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
        
        # sort training data into classes by s1 state
        self.sort_ind = {}
        self.valid_sort_ind = {}
        for name in class_names:
            self.sort_ind[name] = np.where(self.s1==name)[0]
            self.valid_sort_ind[name] = np.where(self.valid_s1==name)[0]
        
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
                bounds = []
                nonflat = []  # ensures that training data only has the correct "non-flat" shape
                              # to avoid issues with calculating the boundary points
                valid_bounds = []
                valid_nonflat = []
                # calculate first and points in each profile with large increases
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
                # calculate the point in each profile with the largest increase
                hbound = (np.argmax(h_prof[:,1:]-h_prof[:,:-1],axis=1)+1)/200
                valid_hbound = (np.argmax(valid_h_prof[:,1:]-valid_h_prof[:,:-1],axis=1)+1)/200
                max_He = np.max(he_prof,axis=1)
                valid_max_He = np.max(valid_he_prof,axis=1)
                dep_He = (np.argmax(he_prof[:,1:]-he_prof[:,:-1],axis=1)+1)/200
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
                warnings.warn(f"no training data available for s1 state {state}. model will return random results")
                
            
            elif state in ['H-rich_Central_He_depleted',
                         'H-rich_Central_C_depletion']:
                history = model.fit(inputs,np.array(bounds),
                                    epochs=training_epochs,verbose=0,callbacks=[callback],
                                    validation_data=(valid_inputs,
                                                     np.array(valid_bounds)))
            
            else:
                history = model.fit(inputs[nonflat],np.array(bounds)[nonflat],
                                    epochs=training_epochs,verbose=0,callbacks=[callback],
                                    validation_data=(valid_inputs[valid_nonflat],
                                                     np.array(valid_bounds)[valid_nonflat]))
                        
            b_models[state] = model
            loss_history[state] = np.array([history.history['loss'],history.history['val_loss']])
            print(f"finished {state}")
        return b_models, loss_history
            
    def predict_single(self,initial,center_H,surface_H,center_He,surface_He,s1):
        """Predict a profile for a single binary
        Args:
            initial (array-like) : initial position of binary
            center_H (float) : center H mass fraction 
            surface_H (float) : surface H mass fraction
            center_He (float) : center He mass fraction 
            surface_He (float) : surface He mass fraction
            s1 (str) : final star 1 state
        Returns:
            H (array-like) : predicted H mass fraction profile 
            He (array-like) : predicted He mass fraction profile
        """
        if "stripped_He" in s1:
            H = np.zeros(200)
        
        if s1 == "stripped_He_non_burning":
            He = np.ones(200)*surface_He
        
        if s1 in ['stripped_He_Core_He_burning',
                 'stripped_He_Core_C_burning',
                 'stripped_He_Central_He_depleted',
                 'stripped_He_Central_C_depletion']:
            b = self.bounds_models[s1](tf.convert_to_tensor([initial])).numpy()[0]                        
            He = np.ones(200) * center_He
            He[int(b[1]*200):] = surface_He
            f_He = interp1d(b,[center_He,surface_He],fill_value="extrapolate")
            # use f to construct the profile points in the shell burning region
            He[int(b[0]*200):int(b[1]*200)] = f_He(np.linspace(0,1,200)[int(b[0]*200):int(b[1]*200)])
            
        if s1 in ["H-rich_non_burning","None"]:
            H = np.ones(200)*surface_H
            He = np.ones(200)*surface_He
            
        if s1 in ["WD","NS","BH"]:
            H = np.ones(200)*np.nan
            He = np.ones(200)*np.nan
        
        # construct step-shaped profile
        if s1 in ["H-rich_Central_He_depleted",
                  "H-rich_Central_C_depletion"]:
            b = self.bounds_models[s1](tf.convert_to_tensor([initial])).numpy()[0]            
            H = np.ones(200)*center_H
            H[int(b[0]*200):] = surface_H
            He = np.ones(200) * center_He
            He[int(b[1]*200):] = b[2]
            He[int(b[0]*200):] = surface_He
        
        # construct shell-shaped profile
        if s1 in ["H-rich_Shell_H_burning",
                  "H-rich_Core_H_burning",
                  "H-rich_Core_He_burning",
                  "H-rich_Core_C_burning"]:
            b = self.bounds_models[s1](tf.convert_to_tensor([initial])).numpy()[0]                        
            new = np.ones([2,200]) * np.array([center_H,center_He])[:,np.newaxis]
            new[:,int(b[1]*200):] = np.array([surface_H,surface_He])[:,np.newaxis]
            f_H = interp1d(b,[center_H,surface_H],fill_value="extrapolate")
            f_He = interp1d(b,[center_He,surface_He],fill_value="extrapolate")
            # use f to construct the profile points in the shell burning region
            new[0][int(b[0]*200):int(b[1]*200)] = f_H(np.linspace(0,1,200)[int(b[0]*200):int(b[1]*200)])
            new[1][int(b[0]*200):int(b[1]*200)] = f_He(np.linspace(0,1,200)[int(b[0]*200):int(b[1]*200)])
            H = new[0]
            He = new[1]
            
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
        # IF interpolate H mass fraction values at center, surface; star 1 state
        center_h_vals = self.interp.test_interpolator(10**inputs)[:,self.c_h_ind]
        surface_h_vals = self.interp.test_interpolator(10**inputs)[:,self.s_h_ind]
        center_he_vals = self.interp.test_interpolator(10**inputs)[:,self.c_he_ind]
        surface_he_vals = self.interp.test_interpolator(10**inputs)[:,self.s_he_ind]
        s1_vals = self.interp.test_classifiers(10**inputs)['S1_state']        
        
        # generate predicted profiles
        pred_profiles = []
        for i in range(len(inputs)):
            pred_H,pred_He = self.predict_single(inputs[i],center_h_vals[i],surface_h_vals[i],
                                                 center_he_vals[i],surface_he_vals[i],s1_vals[i])
            pred_profiles.append([pred_H,pred_He])
           
        # IF interpolate final masses, generate mass enclosed profile coordinates
        m1_ind = self.model_IF.interpolators[0].out_keys.index("star_1_mass")
        pred_mass = self.model_IF.interpolators[0].test_interpolator(10**inputs)[:,m1_ind]
        mass_coords = np.linspace(0,1,200)*pred_mass[:,np.newaxis] 
        h_profiles = np.array(pred_profiles)[:,0]
        he_profiles = np.array(pred_profiles)[:,1]
            
        return mass_coords, h_profiles, he_profiles