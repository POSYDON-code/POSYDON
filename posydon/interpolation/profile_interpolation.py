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

    def __init__(self, train_path, valid_path,
                 profile_names=['radius','logRho', 'x_mass_fraction_H',
                                'y_mass_fraction_He','z_mass_fraction_metals',
                                'omega','energy']):
        """Extracts profile data from '.h5' grid files and saves to file.
        Args:
            train_path (str) : path/name of '.h5' file for training data.
            valid_path (str) : path/name of '.h5' file for testing data.
            profile_names (array-like) : list of profile quantities to extract.
        """
        self.names = profile_names
        
        # extract testing data
        print("extracting testing data")
        valid = PSyGrid(valid_path)  # load PSyGrid object for testing grid
        self.valid_scalars = pd.DataFrame()
        self.valid_profiles = []
        testing_failed = []
        
        try:
            print(valid.final_values["S1_state"][0])
        except:
            print("grid does not have S1_state information")


        for i in range(len(valid)):
            try:
                scalars,profiles = self.scrape(valid,i)
                self.valid_scalars = self.valid_scalars.append(scalars,ignore_index=True)
                self.valid_profiles.append(profiles)
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
                scalars,profiles = self.scrape(train,i)
                self.scalars = self.scalars.append(scalars,ignore_index=True)
                self.profiles.append(profiles)
            except:
                training_failed.append(i)
                pass
        warnings.warn(f"{len(training_failed)} training binaries failed")
        

    def scrape(self,grid,ind):
        """Extracts profile data from one MESA run.
        Args:
            grid (obj) : PSyGrid object.
            ind (int) : index of run to be scraped.
        Returns:
            scalars (array-like) : dictionary containing initial 
                                   m1, m2, p and final s1_state, final_m1.
            profiles (array-like) : all N specified profiles, shape (N,200).
        """
        # open individual run as a DataFrame
        df = pd.DataFrame(grid[ind]['final_profile1'])
        df = df.sort_values(by="mass")
        df = df.reset_index(drop=True)

        # grab input values, final star 1 state, final star 1 mass
        final_m1 = grid.final_values["star_1_mass"][ind]
        scalars = {"m1":grid.initial_values["star_1_mass"][ind],
                   "m2":grid.initial_values["star_2_mass"][ind],
                   "p":grid.initial_values["period_days"][ind],
                   "s1_state":grid.final_values["S1_state"][ind],
                   "final_m1":final_m1}

        # grab output vectors, interpolate to normalize
        profiles=np.zeros([len(self.names),200])
        for i,prof in enumerate(self.names):
            if prof in df.columns:
            
                f = interp1d(df['mass']/final_m1,df[prof],
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
    
    def __init__(self,load_interpolator=None):
        """Interfaces with other classes, trains models and predicts profiles.
        Args:
            load_interpolator (str) : optional path/name of interpolator file to be loaded
        """

        if load_interpolator is not None:
            self.load(load_interpolator)
        
    def load_profiles(self,filename):
        """Load extracted profile data.
        Args:
            filename (str) : path/name of '.pkl' file to be loaded.
        Returns:
            self.profiles (array-like) : M profiles for each of N training binaries, shape (N,M,200).
            self.scalars (array-like) : DataFrame containing initial and final conditions for N training binaries.
            self.valid_profiles (array-like) : M profiles for each of N testing binaries, shape (N,M,200).
            self.valid_scalars (array-like) : DataFrame containing initial and final conditions for N testing binaries.
        """
        with open(filename, 'rb') as f:
            myattrs = pd.read_pickle(f)
            for key in myattrs:
                setattr(self, key, myattrs[key])
        self.profiles = np.array(self.profiles)
        self.valid_profiles = np.array(self.valid_profiles)
    
    def train(self,IF_interpolator):
        """Trains models for density and H mass fraction profile models. 
        Args:
            IF_interpolator (str) : path to '.pkl' file for IF interpolator.
        """
        # processing
        linear_initial = np.transpose([
            self.scalars["m1"],self.scalars["m2"],self.scalars["p"]])
        initial = np.log10(np.array(linear_initial))
        final_m1 = self.scalars["final_m1"].astype(np.float64)

        valid_linear_initial = np.transpose([self.valid_scalars["m1"],
                                             self.valid_scalars["m2"],
                                             self.valid_scalars["p"]]) 
        valid_initial = np.log10(np.array(valid_linear_initial))
        valid_final_m1 = self.valid_scalars["final_m1"].astype(np.float64)
        
        # instantiate and train H mass fraction profile model
        h_ind = self.names.index("x_mass_fraction_H")
        self.h = X_H(initial, self.profiles[:,h_ind], self.scalars["s1_state"], 
                     valid_initial, self.valid_profiles[:,h_ind], 
                     self.valid_scalars["s1_state"], IF_interpolator)

        # instantiate and train density profile model
        dens_ind = self.names.index("logRho")
        self.dens = Density(initial,
                       self.profiles[:,dens_ind],
                       valid_initial,
                       self.valid_profiles[:,dens_ind],
                       IF_interpolator)
        self.dens.train()        

    def predict(self,inputs):
        """Predict density and H mass fraction profiles from inputs.
        Args:
            inputs (array-like) : log-space initial conditions of N binaries to predict, shape (N,3).
        Returns:
            mass_coords (array-like) : linear-scale mass enclosed profile coordinates.
            density_profiles (array_like) : log-scale density profile coordinates.
            h_profiles (array_like) : H mass fraction profile coordinates
        """
        mass_coords, density_profiles = self.dens.predict(inputs)
        mass_coords, h_profiles = self.h.predict(inputs)
        
        return mass_coords, density_profiles, h_profiles
        
                    
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
            np.where(arr_copy<arr[-1], arr[-1], arr_copy)
        
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
            valid_initial (array-like) : log-space initial conditions for testing data.
            valid_profiles (array-like) : final density profiles for testing data.
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
        weights_unscaled = self.pca.transform(profiles_norm)
        self.scaling = np.std(weights_unscaled,axis=0)
        self.weights = weights_unscaled/self.scaling  # scaled PCA weights
        
        # process testing data
        self.valid_initial = valid_initial
        self.valid_rho_min = np.min(valid_profiles,axis=1)
        valid_rho_max = np.max(valid_profiles,axis=1)
        valid_profiles_norm = (valid_profiles-self.valid_rho_min[:,np.newaxis])\
                            /(valid_rho_max-self.valid_rho_min)[:,np.newaxis]
        valid_weights_unscaled = self.pca.transform(valid_profiles_norm)
        self.valid_weights = valid_weights_unscaled/self.scaling
        
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
                
    def train(self,loss=losses.MeanSquaredError()):
        """Trains NN models. 
        Args: 
            loss (object) : loss function for training.
        """
        print("training on PCA weights...")
        
        self.model_prof.compile(optimizers.Adam(clipnorm=1),loss=loss)
        callback = tf.keras.callbacks.EarlyStopping(monitor="loss",patience=200)
        history = self.model_prof.fit(self.initial,self.weights,
                                      epochs=3000,callbacks=[callback],verbose=0,
                                      validation_data=(self.valid_initial,
                                                       self.valid_weights))
        
        self.model_rho.compile(optimizers.Adam(clipnorm=1),loss=loss)
        callback = tf.keras.callbacks.EarlyStopping(monitor="loss",patience=40)
        history = self.model_rho.fit(self.initial,self.rho_min,
                                     epochs=500, callbacks=[callback], verbose=0,
                                     validation_data=(self.valid_initial,
                                                      self.valid_rho_min))
                                     
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
        weights_pred = regress_prof(inputs).numpy()
                
        # predict surface density
        regress_rho = lambda x: self.model_rho(x)
        min_rho = regress_rho(inputs).numpy()[:,0]
        
        # IF interpolate center density
        center_ind = self.model_IF.interpolators[0].out_keys.index('S1_log_center_Rho')
        max_rho = self.model_IF.interpolators[0].test_interpolator(10**inputs)[:,center_ind]                    
             
        # reconstruct profile
        norm_prof = self.pca.inverse_transform(weights_pred*self.scaling)
        density_profiles = norm_prof*(max_rho[:,np.newaxis]-min_rho[:,np.newaxis]) \
                           + min_rho[:,np.newaxis]
        
        # IF interpolate final mass, construct mass enclosed profile coordinates
        m1_ind = self.model_IF.interpolators[0].out_keys.index("star_1_mass")
        pred_mass = self.model_IF.interpolators[0].test_interpolator(10**inputs)[:,m1_ind]
        mass_coords = np.linspace(0,1,200)*pred_mass[:,np.newaxis] 
        
        return mass_coords,density_profiles

    
class X_H:
    
    def __init__(self,initial,profiles,s1,valid_initial, 
                 valid_profiles,valid_s1,IF_interpolator):
        """Creates and trains H mass fraction profile model.
        Args:
            initial (array-like) : log-space initial conditions for training data.
            profiles (array-like) : final H mass fraction profiles for training data. 
            s1 (array-like) : final star 1 state for training data.
            valid_initial (array-like) : log-space initial conditions for testing data.
            valid_profiles (array-like) : final H mass fraction profiles for testing data.
            valid_s1 (array-like) : final star 1 state for testing data.
            IF_interpolator (string) : path to .pkl file for IF interpolator for central density, final mass values.
        """
        self.initial = initial
        self.profiles = profiles
        self.s1 = s1
        self.valid_initial = valid_initial
        self.valid_profiles = valid_profiles
        self.valid_s1 = valid_s1

        # load IF interpolator
        self.model_IF = IFInterpolator()  # instantiate POSYDON initial-final interpolator object
        self.model_IF.load(filename=IF_interpolator)
        self.c_ind = self.model_IF.interpolators[0].out_keys.index("S1_center_h1")
        self.s_ind = self.model_IF.interpolators[0].out_keys.index("S1_surface_h1")
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
        self.sort_ind = {class_name:[] for class_name in self.class_names}
        for i in range(len(self.s1)):
            state = self.s1[i]
            self.sort_ind[state].append(i)
        
        # sort testing data into classes by s1 state
        self.valid_sort_ind = {class_name:[] for class_name in self.class_names}
        for i in range(len(self.valid_s1)):
            state = self.valid_s1[i]
            self.valid_sort_ind[state].append(i)
        
        # create and train models for profile boundaries
        self.bounds_models = self.learn_bounds()
        
    def learn_bounds(self):
        """Creates and trains NNs to predict boundary points for each star 1 state.
        Returns:
            b_models (array-like) : dictionary containing boundary models.
        """
        b_models = {}
        for state in ['WD','NS','BH',
                      'H-rich_Central_He_depleted',
                      'H-rich_Central_C_depletion',
                      'H-rich_Shell_H_burning',
                      'H-rich_Core_H_burning',
                      'H-rich_Core_He_burning',
                      'H-rich_Core_C_burning']:
            
            # identify input/output training data for class
            indices = self.sort_ind[state]
            inputs = self.initial[indices]
            prof = self.profiles[indices]

            # identify input/output testing data for class
            valid_indices = self.valid_sort_ind[state]
            valid_inputs = self.valid_initial[valid_indices]
            valid_prof = self.valid_profiles[valid_indices]
            
            # calculate boundary points for training and testing data 
            if "burning" in state:  # these classes' profile shapes have 2 boundary points
                outs = 2
                bounds = []
                nonflat = []  # TODO: explain this
                valid_bounds = []
                valid_nonflat = []
                # calculate first and points in each profile with large increases
                for i in range(len(indices)):
                    try:
                        diff = np.where(prof[i][1:]-prof[i][:-1]>0.002)[0]
                        bounds.append(np.array([diff[0]+1,diff[-1]+1])/200)
                        nonflat.append(i)
                    except:
                        bounds.append([np.nan,np.nan])
                for i in range(len(valid_indices)):
                    try:
                        diff = np.where(valid_prof[i][1:]-valid_prof[i][:-1]>0.002)[0]
                        valid_bounds.append(np.array([diff[0]+1,diff[-1]+1])/200)
                        valid_nonflat.append(i)
                    except:
                        valid_bounds.append([np.nan,np.nan])
                        
            else: # the rest of the profile shapes have 1 boundary point
                outs = 1
                # calculate the point in each profile with the largest increase
                bounds = (np.argmax(prof[:,1:]-prof[:,:-1],axis=1)+1)/200
                valid_bounds = (np.argmax(valid_prof[:,1:]-valid_prof[:,:-1],axis=1)+1)/200
            
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
            callback = tf.keras.callbacks.EarlyStopping(monitor='loss', patience=50)
            
            if len(indices)==0:
                warnings.warn(f"no training data available for s1 state {state}. model will return random results")
                
            elif "burning" in state:
                history = model.fit(inputs[nonflat],np.array(bounds)[nonflat],
                                    epochs=500,verbose=0,callbacks=[callback],
                                    validation_data=(valid_inputs[valid_nonflat],
                                                     np.array(valid_bounds)[valid_nonflat]))
            else:
                history = model.fit(inputs,np.array(bounds),
                                    epochs=500,verbose=0,callbacks=[callback],
                                    validation_data=(valid_inputs,
                                                     np.array(valid_bounds)))
            
            b_models[state] = model
            print(f"finished {state}")
        return b_models
            
    def predict_single(self,initial,center,surface,s1):
        """Predict a profile for a single binary
        Args:
            initial (array-like) : initial position of binary
            center (float) : center H mass fraction 
            surface (float) : surface H mass fraction
            s1 (str) : final star 1 state
        """
        if "stripped_He" in s1:
            return np.zeros(200)
        
        if s1 in ["H-rich_non_burning","None"] or surface-center<0.002:
            return np.ones(200)*surface
        
        # construct step-shaped profile
        if s1 in ["H-rich_Central_He_depleted",
                  "H-rich_Central_C_depletion",
                  "WD","NS","BH"]:
            b = self.bounds_models[s1](tf.convert_to_tensor([initial])).numpy()[0]            
            new = np.ones(200)*center
            new[int(b*200):] = surface
            return new
        
        # construct shell-shaped profile
        if s1 in ["H-rich_Shell_H_burning",
                  "H-rich_Core_H_burning",
                  "H-rich_Core_He_burning",
                  "H-rich_Core_C_burning"]:
            b = self.bounds_models[s1](tf.convert_to_tensor([initial])).numpy()[0]                        
            new = np.ones(200)*center
            new[int(b[1]*200):] = surface
            # "f" defines the shape of the H abundance profile in the shell burning region
            f = interp1d(b,[center,surface],fill_value="extrapolate")
            # use f to construct the profile points in the shell burning region
            new[int(b[0]*200):int(b[1]*200)] = f(np.linspace(0,1,200)[int(b[0]*200):int(b[1]*200)])
            return new
        
    def predict(self,inputs):
        """Predict H mass fraction profiles from inputs.
        Args:
            inputs (array-like) : log-space initial conditions of N binaries to predict, shape (N,3).
        Returns:
            mass_coords (array-like) : linear-scale mass enclosed profile coordinates.
            h_profiles (array_like) : H mass fraction profile coordinates.
        """
        # IF interpolate H mass fraction values at center, surface; star 1 state
        center_vals = self.interp.test_interpolator(10**inputs)[:,self.c_ind]
        surface_vals = self.interp.test_interpolator(10**inputs)[:,self.s_ind]
        s1_vals = self.interp.test_classifiers(10**inputs)['S1_state']        
        
        # generate predicted profiles
        pred_profiles = []
        for i in range(len(inputs)):
            prediction = self.predict_single(inputs[i],center_vals[i],surface_vals[i],s1_vals[i])
            pred_profiles.append(prediction)
        h_profiles = np.array(pred_profiles)  
           
        # IF interpolate final masses, generate mass enclosed profile coordinates
        m1_ind = self.model_IF.interpolators[0].out_keys.index("star_1_mass")
        pred_mass = self.model_IF.interpolators[0].test_interpolator(10**inputs)[:,m1_ind]
        mass_coords = np.linspace(0,1,200)*pred_mass[:,np.newaxis] 
            
        return mass_coords, h_profiles