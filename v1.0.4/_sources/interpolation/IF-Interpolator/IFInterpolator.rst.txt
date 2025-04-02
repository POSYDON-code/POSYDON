
.. _IFInterpolator:

###########################
Initial-Final Interpolation
###########################


We showcase the initial-final interpolator which plays a critical role in the
evolving binary populations. To use the initial-final interpolator we first
import the ``IFInterpolator`` object from the POSYDON library.

.. code-block:: python

  # importing interpolator
  from posydon.interpolation.IF_interpolation import IFInterpolator

Loading a Pretrained Interpolator
===============================

To load a pretrained interpolator we need to
pass in the ``filename`` argument into the ``IFInterpolator``
constructor which specifies the path to a .pkl file where
the pretrained interpolator can be loaded from. POSYDON provides
various pretrained models whose corresponding .pkl files
can be found in the data directory of the POSYDON repository.

.. code-block:: python

  model = IFInterpolator(filename = "path/to/file.pkl")


Training the Interpolator
=========================

The interpolator can be trained using an instance of the ``PSyGrid`` 
class which can be constructed by running ones own simulations or
by loading a simulation from an h5 file stored in the data directory
of the POSYDON repository. For more details on the ``PSyGrid`` class
please visit the ``PSyGrid`` documentation. The ``IFInterpolator``
class relies on the ``BaseIFInterpolator`` class to perform the interpolation 
so parameters to construct instances of the ``BaseIFInterpolator`` classes are
required to construct the ``IFInterpolator``. These parameters are:

1. ``interp_method``: the interpolation method to be used can be either
``linear``, ``1NN``, or a list specifying which interpolation method to
be used for each type of track. If ``interp_classes`` is specified and this
parameter is not a list then the interpolator will use the specified method
for all classes.

2. ``interp_classes``: a list of classes that the simulation tracks can
fall into. Usually specified as the mass transfer type. This only needs
be specified if ``interp_method`` is a list.

3. ``class_method``: the classification method to be used, either ``kNN`` or
``1NN``.

4. ``in_keys``: the keys to be used as the input to the interpolator, by default
these are ``star_1_mass``, ``star_2_mass``, and ``period_days``.

5. ``out_keys``: the keys for which the interpolator is supposed to provide values,
by default all keys are used.

6. ``in_scaling``: The scalings for the input keys, by default these scalings are
optimized through Monte Carlo Cross Validation.

7. ``out_scaling``: The scalings for the output keys, by default these scalings
are optimized through Monte Carlo Cross Validation.

8. ``c_keys``: A list of strings specifying which classifiers are to be trained

9. ``c_key``: A string specifying by which class the interpolator should interpolate
binaries. Only to be specified in the MCInterpolator case.

For most applications specifying only the first four parameters is recommended.

.. code-block:: python

    from posydon.grids.psygrid import PSyGrid
    from posydon.interpolation.IF_interpolation import IFInterpolator

    grid = PSyGrid("path/to/h5/file.h5") # loading grid from h5 file

    interp = IFInterpolator(grid = grid, interpolators = [
        { 
            "interp_method": ["linear", "linear", "linear"], 
            "interp_classes": ["no_MT", "stable_MT", "unstable_MT"],
            "out_keys": first,
            "class_method": "kNN",
            "c_keys": ["interpolation_class"],
            "c_key": "interpolation_class"
        }, 
        { 
            "interp_method": ["linear", "linear", "linear", "linear"], 
            "interp_classes": ["None", "BH", "WD", "NS"],
            "out_keys": second,
            "class_method": "kNN",
            "c_keys": ['S1_direct_state'],
            "c_key": 'S1_direct_state'
        },
        { 
            "interp_method": ["linear", "linear", "linear", "linear"], 
            "interp_classes": ["None", "BH", "WD", "NS"],
            "out_keys": third,
            "class_method": "kNN",
            "c_keys": ['S1_Fryer+12-rapid_state'],
            "c_key": 'S1_Fryer+12-rapid_state'
        },
        { 
            "interp_method": ["linear", "linear", "linear", "linear"], 
            "interp_classes": ["None", "BH", "WD", "NS"],
            "out_keys": fourth,
            "class_method": "kNN",
            "c_keys": ['S1_Fryer+12-delayed_state'],
            "c_key": 'S1_Fryer+12-delayed_state'
        },
        { 
            "interp_method": ["linear", "linear", "linear", "linear"], 
            "interp_classes": ["None", "BH", "WD", "NS"],
            "out_keys": fifth,
            "class_method": "kNN",
            "c_keys": ['S1_Sukhbold+16-engineN20_state'],
            "c_key": 'S1_Sukhbold+16-engineN20_state'
        },
        { 
            "interp_method": ["linear", "linear", "linear", "linear"], 
            "interp_classes": ["None", "BH", "WD", "NS"],
            "out_keys": sixth,
            "class_method": "kNN",
            "c_keys": ['S1_Patton&Sukhbold20-engineN20_state'],
            "c_key": 'S1_Patton&Sukhbold20-engineN20_state'
        }
    ]) # constructing IFInterpolator
    
    interp.train() # training interpolator


Using the Interpolator
======================

Once the interpolator has been trained or loaded from a .pkl file it can be used
to accomplish various tasks which most commonly are to classify a track into its class
given an input vector and or to approximate a final vector given an input vector.

.. code-block:: python

    from posydon.binary_evol.binarystar import BinaryStar
    from posydon.binary_evol.singlestar import SingleStar


    binary = BinaryStar(**binary_params,
                        star_1=SingleStar(**star1_params),
                        star_2=SingleStar(**star2_params)) # creating binary, refer to BinaryStar documentation

    interpolation, classification = interp.evaluate(binary) # evaluating returns a tuple of dictionaries


Finally a trained interpolator can be easily saved by specifying a path to a .pkl file
where the interpolator will be saved to.

.. code-block:: python

   model.save("path/to/file.pkl") # saving interpolator


