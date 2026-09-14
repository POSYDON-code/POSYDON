.. _initial-final-interp:

############################################
Initial-Final Classification & Interpolation
############################################

We showcase the initial-final interpolator which plays a critical role in the
evolving binary populations. To use the initial-final interpolator we first
import the ``IFInterpolator`` object from the POSYDON library.

.. code-block:: python

  # importing interpolator
  from posydon.interpolation.IF_interpolation import IFInterpolator

.. note::

  Familiarity with the overall system, which can be gained by referencing
  section 3 of `2411.02376 <https://arxiv.org/abs/2411.02376>`_, is
  recommended to fully understand the interpolator and the parameters
  described below.

Loading a Pretrained Interpolator
==================================

To load a pretrained interpolator, first construct an ``IFInterpolator`` in
loading mode by passing ``load = True``, then call its ``load`` method with
the path to the :samp:`pkl` file. POSYDON provides various pretrained models
whose corresponding :samp:`pkl` files can be found in the data directory of
the POSYDON repository.

.. code-block:: python

  model = IFInterpolator(load = True)
  model = model.load(filename = "path/to/file.pkl")

``load`` accepts the following optional arguments:

1. ``filename``: the path to the :samp:`pkl` file the interpolator should be
   loaded from.

2. ``sn_model``: the supernova model whose classifiers should be kept after
   loading (default ``"SN_MODEL_v2_01"``). To keep the interpolator small,
   only the classifiers associated with this supernova model are retained,
   along with the general-purpose ``interpolation_class`` and ``mt_history``
   classifiers.

3. ``nearest_neighbor_mode``: a boolean (default ``False``) that, when
   ``True``, forces the interpolator to always fall back to nearest-neighbor
   interpolation instead of the default Delaunay-triangulation-based
   interpolation.

``load`` returns the fully-populated interpolator object (it does not modify
``model`` in place), so make sure to assign its return value as shown above.

Training the Interpolator
=========================

Training the interpolator requires two ``PSyGrid`` instances passed together
as a list: a regularly-sampled training grid and a separate validation grid
used to select hyperparameters. Each can be constructed by running one's own
simulations or by loading a simulation from an :samp:`h5` file stored in the
data directory of the POSYDON repository. For more details on the
``PSyGrid`` class please visit the ``PSyGrid`` documentation.

The ``IFInterpolator`` constructor accepts the following parameters:

1. ``grids``: a list of two ``PSyGrid`` objects, ``[training_grid,
   validation_grid]``. The first is used to fit the classifiers and
   interpolants; the second is used to evaluate different hyperparameter and
   normalization choices during training.

2. ``in_keys``: the keys to be used as the input to the interpolator, e.g.
   ``star_1_mass``, ``star_2_mass``, and ``period_days``.

3. ``out_keys``: a dictionary whose keys are the names of the discrete
   classification schemes (classifiers) to be trained, and whose values are
   lists of the continuous quantities to be interpolated using the
   predictions of the corresponding classifier. This dictionary **must**
   contain ``"interpolation_class"`` as one of its keys, since it is used to
   determine the mass-transfer type of a binary. Any classifier keys tied to
   a specific supernova engine should follow the naming convention
   ``S1_<sn_model>_CO_interpolation_class``, where ``<sn_model>`` matches the
   ``sn_model`` string later passed to ``evaluate`` (e.g. the default
   ``sn_model = "SN_MODEL_v2_01"`` corresponds to the key
   ``S1_SN_MODEL_v2_01_CO_interpolation_class``).

4. ``max_k``: the maximum number of neighbors ``k`` considered when
   optimizing the kNN classifier for each classification scheme against the
   validation grid.

5. ``load``: a boolean (default ``False``) that puts the constructor into
   loading mode, skipping the training-grid setup described above. Only used
   in conjunction with the ``load`` method described in the previous
   section.

6. ``nearest_neighbor_mode``: a boolean (default ``False``) that forces the
   interpolator to always use nearest-neighbor interpolation rather than
   Delaunay-triangulation-based interpolation.

Unlike earlier versions of the interpolator, individual input/output
normalization schemes, interpolation methods (e.g. ``linear`` vs ``1NN``),
and classification methods no longer need to be specified by hand: calling
``train()`` automatically searches over the available normalization schemes
and values of ``k`` and selects the combination that performs best against
the validation grid for each classification scheme.

Provided below is an example instantiation of the ``IFInterpolator`` class.

.. code-block:: python

    from posydon.grids.psygrid import PSyGrid
    from posydon.interpolation.IF_interpolation import IFInterpolator

    training_grid = PSyGrid("path/to/training/h5/file.h5")
    validation_grid = PSyGrid("path/to/validation/h5/file.h5")

    out_keys = {
        "interpolation_class": first,
        "S1_direct_state": second,
        "S1_Fryer+12-rapid_state": third,
        "S1_Fryer+12-delayed_state": fourth,
        "S1_Sukhbold+16-engineN20_state": fifth,
        "S1_Patton&Sukhbold20-engineN20_state": sixth,
    } # each value (first, ..., sixth) is a list of continuous output keys

    interp = IFInterpolator(
        grids = [training_grid, validation_grid],
        in_keys = ["star_1_mass", "star_2_mass", "period_days"],
        out_keys = out_keys,
        max_k = 20
    ) # constructing IFInterpolator

    interp.train() # training interpolator

After training, ``interp.stats(_print = True)`` can be called to report the
percentage of validation points, per classification scheme, whose initial
conditions fell outside the convex hull of the training grid (and therefore
fell back to nearest-neighbor interpolation).

Using the Interpolator
======================

Once the interpolator has been trained or loaded from a :samp:`pkl` file it
can be used to accomplish various tasks which most commonly are to classify a
track into its class given an input vector and/or to approximate a final
vector given an input vector.

.. code-block:: python

    from posydon.binary_evol.binarystar import BinaryStar
    from posydon.binary_evol.singlestar import SingleStar


    binary = BinaryStar(**binary_params,
                        star_1=SingleStar(**star1_params),
                        star_2=SingleStar(**star2_params)) # creating binary, refer to BinaryStar documentation

    interpolated_values, classes, meta_data = interp.evaluate(binary) # evaluating returns a tuple of three items

``evaluate`` accepts either a single ``BinaryStar`` instance or a numpy array
of shape ``(N, 3)`` containing ``[star_1_mass, star_2_mass, period_days]``
for batch evaluation, along with an optional ``sn_model`` string (default
``"SN_MODEL_v2_01"``) selecting which supernova engine's compact-object
classification should be used. It returns a tuple of three items:

1. ``interpolated_values``: a numpy array of the interpolated continuous
   quantities for each evaluated binary.

2. ``classes``: a numpy array where each row contains two predicted class
   labels — the mass-transfer type (``interpolation_class``) and the
   compact-object type predicted for the chosen ``sn_model``.

3. ``meta_data``: a list of dictionaries, one per evaluated binary,
   containing diagnostic information about the interpolation (e.g. the
   barycentric ``weights``, the neighboring initial conditions ``ics``, the
   queried initial condition ``ic``, and the raw ``interpolated`` values).

Finally a trained interpolator can be easily saved by specifying a path to
a :samp:`pkl` file where the interpolator will be saved to.

.. code-block:: python

   model.save("path/to/file.pkl") # saving interpolator