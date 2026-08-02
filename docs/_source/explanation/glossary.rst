.. _glossary:

#####################
Glossary of Concepts
#####################

POSYDON's documentation uses a number of specialized terms. This glossary
defines the most common ones so that the rest of the pages share a common
vocabulary. The terms are grouped by topic.

Stellar and binary evolution
-----------------------------

``SingleStar``
    The class that stores the properties and evolutionary history of one star.

``BinaryStar``
    The class that stores the properties and evolutionary history of a single
    binary system, together with its two ``SingleStar`` components.

``HMS``
    Hydrogen main-sequence star. Shorthand used in grid run types (e.g.
    ``HMS-HMS`` for a binary of two hydrogen-rich main-sequence stars).

``DT``
    "Detached" phase: the phase of a binary in which neither star fills its
    Roche lobe and there is no mass transfer. The evolution of a detached
    binary (wind-driven mass loss, orbital shrinkage, the two stars evolving)
    is handled by the ``DT`` (detached) routines.

``CE``
    "Common envelope": a phase in which one star dives deep into the envelope
    of its binary companion, dramatically shrinking the orbit and ejecting the
    envelope. Common-envelope evolution is physically uncertain and in POSYDON
    it is modeled with a prescription (see the code module
    ``binary_evol.CE``).

``SN``
    "Supernova": the explosive end of a star. The ``binary_evol.SN`` module
    determines the remnant (neutron star or black hole) and applies the
    relevant natal kick to the binary.

``Natal kick``
    The recoil velocity given to a compact remnant (and hence to the binary) at
    the moment of a supernova.

``Roche lobe``
    The region of space around a star within which matter is gravitationally
    bound to it. When a star fills its Roche lobe it starts transferring mass to
    its companion.

``termination flag``
    The label stored by MESA (and kept in a ``PSyGrid``'s final values) that
    says how a simulation ended, e.g. that the star reached a supernova, that
    the binary came into contact, or that the run was stopped for another
    reason. Termination flags are used to select which runs to "rerun" (see
    :ref:`run-hpc`).

Grids, data and processing
-----------------------------

``PSyGrid``
    The object that stores the processed data of a MESA grid (initial values,
    final values, downsampled histories, final profiles) and the metadata about
    how it was run. See :ref:`psygrid`.

``grid``
    A collection of MESA simulations spanning a set of initial conditions. See
    :ref:`grid-types`.

``run type``
    The type of binaries a grid models, e.g. ``HMS-HMS``. See :ref:`grid-types`.

``downsampled history``
    The evolution of a simulated star or binary stored at a reduced time
    resolution to save disk space.

``metallicity set``
    The collection of metallicities covered by a dataset. DR2 provides data for
    eight metallicities; see :ref:`data-releases`.

``interpolation``
    Estimating the value of a *continuous* property at an arbitrary point of the
    grid from its neighbours, rather than re-running the simulation.

``classification``
    Deciding a *discrete* outcome (e.g. "black hole" vs "neutron star", or
    whether an object fills its Roche lobe).

POSYDON uses both interpolation and classification as part of its
machine-learning layer; together they are the fast surrogates that make
population synthesis feasible. See :ref:`machine-learning`.

``initial-final interpolator``
    An interpolator that maps the initial conditions of an object or binary to
    its final conditions. Used to step binaries quickly during population
    synthesis.

``profile interpolator``
    An interpolator that reconstructs the interior profile (structure) of a
    star. Experimental in POSYDON; see the ``ProfileInterpolator`` reference.

``population synthesis``
    The process of evolving a large number of stars / binaries (a population)
    quickly using the trained surrogates, then aggregating them into a
    synthetic population. See :ref:`workflow`.

The pieces of a population model
----------------------------------

``flow chart``
    The state machine that defines which step a binary moves to next, given its
    current "state" and "event". Configuring the flow chart is how you customize
    the evolution of your binaries. See :ref:`flow-chart`.

``step``
    One unit of evolution in the flow chart. A step takes a ``BinaryStar`` and
    advances it through a particular phase (e.g. a detached phase step, a
    common-envelope step, a supernova step).

``user_modules``
    A directory where you can drop your own code and import it within the
    POSYDON namespace without modifying the core. See :ref:`user_modules`.

``post-processing pipeline``
    The sequence of steps that extracts the important data from raw MESA output
    and stores it in a ``PSyGrid``. See :ref:`processing-pipeline`.

.. seealso::
   :ref:`workflow` places all of these terms in the context of the entire
   pipeline, from MESA grids to an analyzed synthetic population.
