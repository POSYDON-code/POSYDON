
.. _ProfileInterpolator:

###########################
Profile Interpolation
###########################


We demonstrate the profile interpolator which delivers predictions 
on internal stellar structure at the end of MESA evolution. 
To use the profile interpolator we firstimport the ``ProfileInterpolator`` 
object from the POSYDON library.

.. code-block:: python

  from posydon.interpolation.profile_interpolation import ProfileInterpolator, CompileData


Extracting Profile Data from .h5 Files for Training Interpolators
===============================

To extract profile data from existing grids, we pass arguments 
``train_path`` and ``valid_path`` pointing to .h5 files 
containing instances of the ``PSyGrid`` class from which the
training and testing grids can be loaded. These can be constructed 
by running ones own simulations or by loading a precomputed simulation 
stored in the data directory of the POSYDON repository.
The ``profile_names`` argument denotes which profiles will be extracted 
from the grid and saved to ``filename`` using the function ``save``. 

.. code-block:: python

  compiler = CompileData(train_data = '/path/to/training/grid.h5',
                         valid_data = '/path/to/testing/grid.h5',
                         profile_names=['radius','logRho',
                         'x_mass_fraction_H','y_mass_fraction_He',
                         'z_mass_fraction_metals','omega','energy'])
  compiler.save(filename = '/path/for/profiles/datafile.pkl')
  
  
Loading a Pretrained Interpolator
===============================

To load a pretrained interpolator we need to pass the optional
``load_interpolator`` argument to the ``ProfileInterpolator`` 
instance which specifies the path to a .pkl file where the 
pretrained interpolator can be loaded from. 

.. code-block:: python

  model = ProfileInterpolator(load_interpolator = "/path/to/profile/interpolator.pkl")


Training the Interpolator
=========================

The interpolator is trained on profiles stored in a .pkl file generated with the 
aforementioned ``CompileData`` class. This file is passed using the argument 
``filename``. The profile interpolation is anchored to POSYDON's IF interpolation, 
and so we pass the name of the interpolator file to the ``train`` function using the ``IF_interpolator`` argument. 

.. code-block:: python

  model = ProfileInterpolator(load_interpolator=None)
  model.load_profiles(filename = "/path/to/profiles/datafile.pkl")
  model.train(IF_interpolator = "/path/to/IFinterpolator.pkl")
  

Using the Interpolator
======================

Once the interpolator has been trained or loaded from a .pkl file it can be used
to predict profiles for sets of initial conditions passed through argument ``inputs``.
These initial conditions must be in log space and in shape (N,3) for N binaries. 
The order of the coordinates is star 1 mass, star 2 mass, period. The prediction 
function returns three arrays containing the profiles' coordinates along mass 
enclosed, log density, and Hydrogen mass fraction. The Density and H mass fraction 
profiles have the same mass coordinates. 

.. code-block:: python

  mass_coords, density_profiles, h_profiles = model.predict(inputs)

Finally a trained interpolator can be easily saved by specifying a path to a .pkl file
where the interpolator will be saved to.

.. code-block:: python

 model.save("path/for/profile/interpolator.pkl")
   
Evaluating on Testing Data
==========================

To evaluate the interpolator on the testing grid that was used as validation data in 
training, we can pull the testing data out of the ``ProfileInterpolator`` class as follows:

.. code-block:: python

  valid_initial = np.log10(np.transpose([model.valid_scalars["m1"],
                                         model.valid_scalars['m2'],
                                         model.valid_scalars['p']]))
                                         
  valid_mass_coords = np.transpose(np.array(model.valid_scalars["final_mass"])*np.linspace(0,1,200)[:,np.newaxis])
  valid_density_profiles = model.valid_profiles[:,model.names.index("logRho")]
  valid_H_profiles = model.valid_profiles[:,model.names.index("x_mass_fraction_H")]
                                    