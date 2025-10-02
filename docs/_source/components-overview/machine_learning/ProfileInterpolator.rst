.. _profile-interpolation:

###########################
Profile Interpolation
###########################

We demonstrate the profile interpolator which delivers predictions on internal
stellar structure at the end of MESA evolution.

To use the profile interpolator we first import the ``ProfileInterpolator`` and
``CompileData`` objects from the POSYDON library.

.. code-block:: python

  from posydon.interpolation.profile_interpolation import ProfileInterpolator, CompileData


Extracting Profile Data from :samp:`h5` Files for Training Interpolators
========================================================================

To extract profile data from existing grids, we pass arguments ``train_path``
and ``test_path`` pointing to :samp:`h5` files containing instances of the
``PSyGrid`` class from which the training and testing grids can be loaded. These
can be constructed by running one's own simulations or by loading a precomputed
simulation stored in the data directory of the POSYDON repository. The
``profile_names`` argument denotes which profiles will be extracted from the
grid and saved to ``filename`` using the function ``save``. The Boolean
argument``hms_s2`` should be set to True only to extract star 2 profiles from
the HMS-HMS grid.

.. code-block:: python

  compiler = CompileData(train_path = '/path/to/training/grid.h5',
                         test_path = '/path/to/testing/grid.h5',
                         hms_s2 = False,
                         profile_names=['radius','logRho',
                         'x_mass_fraction_H','y_mass_fraction_He',
                         'z_mass_fraction_metals','omega','energy'])
  compiler.save(filename = '/path/for/profiles/datafile.pkl')


Loading a Pretrained Interpolator
=================================

To load a pretrained interpolator we instantiate a ``ProfileInterpolator``
object instance and then use the ``load`` function which specifies the path
to a :samp:`pkl` file from which the pretrained interpolator can be loaded.

.. code-block:: python

  model = ProfileInterpolator()
  model.load(filename = "/path/to/profile/interpolator.pkl")


Training the Interpolator
=========================

The interpolator is trained on profiles stored in a :samp:`pkl` file (passed
with argument ``filename``) generated with the ``CompileData`` class. The
function ``load_profiles`` also randomly splits the training data into a
training set and validation set based on the argument ``valid_split``, which
specifies the percentage of the training data that should be used for
validation. The profile interpolation is anchored to POSYDON's IF
interpolation, and so we pass the name of the interpolator file to the
``train`` function using the ``IF_interpolator`` argument. The four following
arguments specify the numbers of epochs and level of patience for callback that
will be used for training the neural networks in the density and composition
models. The ``loss_history`` argument provides an option to return the loss
histories for the composition and density model training. Once again, the
Boolean argument ``hms_s2`` should be set to True to build an interpolator for
star 2 profiles from the HMS-HMS grid.

.. code-block:: python

  model = ProfileInterpolator()
  model.load_profiles(filename = "/path/to/profiles/datafile.pkl",
                      valid_split=0.2)
  model.train(IF_interpolator = "/path/to/IFinterpolator.pkl",
              density_epochs=3000,
              density_patience=200,
              comp_bounds_epochs=500,
              comp_bounds_patience=50,
              loss_history=False,
              hms_s2=False)


Using the Interpolator
======================

Once the interpolator has been trained or loaded from a :samp:`pkl` file it can
be used to predict profiles for sets of initial conditions passed through
argument ``inputs``. These initial conditions must be in log space and in shape
(N,3) for N binaries. The order of the coordinates is star 1 mass, star 2 mass,
period. The prediction function returns four arrays containing the profiles'
coordinates along mass enclosed, log density, Hydrogen mass fraction, and
Helium mass fraction. All profiles share the same coordinates.

.. code-block:: python

  mass_coords, density_profiles, h_profiles, he_profiles = model.predict(inputs)

Finally a trained interpolator can be easily saved by specifying a path to a
:samp:`pkl` file where the interpolator will be saved to.

.. code-block:: python

 model.save(filename = "path/for/profile/interpolator.pkl")

Evaluating on Testing Data
==========================

To evaluate the interpolator on the testing grid, we can pull the testing data
out of the ``ProfileInterpolator`` class as follows:

.. code-block:: python

  test_initial = model.test_initial
  test_mass_coords = np.transpose(np.array(model.test_scalars["total_mass"])*np.linspace(0,1,200)[:,np.newaxis])
  test_density_profiles = model.test_profiles[:,model.names.index("logRho")]
  test_H_profiles = model.test_profiles[:,model.names.index("x_mass_fraction_H")]
  test_He_profiles = model.test_profiles[:,model.names.index("y_mass_fraction_He")]

