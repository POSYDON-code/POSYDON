.. _binary-population:

The Binary Population Object
============================

The BinaryPopulation class is the main object for handling the evolution of
binary populations in POSYDON. `posydon-popsyn` and `PopulationRunner` build upon
it to provide a more user-friendly interface for simulating populations.

However, its features are useful if you want to adapt part of the population synthesis,
for example, to test a specific population of binaries with custom parameters.


The initialisation of the BinaryPopulation object
-------------------------------------------------

To create a BinaryPopulation object, you need to provide a set of parameters that define the population's characteristics. 
The easiest way is to load in the parameters from the default population `ini`` file.

```python

import os
import shutil
from posydon.config import PATH_TO_POSYDON

path_to_params = os.path.join(PATH_TO_POSYDON, "posydon/popsyn/population_params_default.ini")
shutil.copyfile(path_to_params, './population_params.ini')

```

Here's an example of how to initialise a BinaryPopulation object:

```python

from posydon.popsyn.binarypopulation import BinaryPopulation
from posydon.popsyn.io import binarypop_kwargs_from_ini

# Load the parameters from the ini file
ini_kwargs = binarypop_kwargs_from_ini('./population_params.ini')

# Create the BinaryPopulation object
population = BinaryPopulation(**ini_kwargs)

print(population.number_of_binaries)
```

This code will create a BinaryPopulation object with the parameters defined in the `population_params.ini` file and print the number of binaries in the population.
You can adapt the ``ini_kwargs`` dictionary to change the parameters of the population, such as the number of binaries and the metallicity.


You can also skip the read in from the ini file and create the BinaryPopulation object directly. 
This performs the same task as the previous example, but automatically.

```python

from posydon.popsyn.binarypopulation import BinaryPopulation

population = BinaryPopulation.from_ini('./population_params.ini')

print(population.number_of_binaries)
```

On initialisation, the BinaryPopulation class will initialise the following attributes:

.. list-table:: Core Attributes
   :widths: 30 70
   :header-rows: 1

   * - Attribute
     - Description
   * - ``kwargs``
     - Dictionary containing all initialization parameters, starting with defaults and updated with user inputs
   * - ``number_of_binaries``
     - Total number of binaries in the population
   * - ``population_properties``
     - SimulationProperties instance containing the evolution steps and configuration
   * - ``metallicity``
     - Metallicity value for the population [Z/Z_sun]
   * - ``metallicities``
     - List of metallicity values, defaults to [metallicity] for single metallicity populations
   * - ``history_verbose``
     - Boolean flag controlling verbosity of evolution history output
   * - ``entropy``
     - Entropy value used for random number generation seeding
   * - ``RNG``
     - NumPy random number generator instance for reproducible sampling
   * - ``manager``
     - PopulationManager instance that handles binary creation, evolution, and data management
   * - ``to_df``
     - Method reference to manager's to_df method for converting binaries to DataFrame
   * - ``to_oneline_df``
     - Method reference to manager's to_oneline_df method for converting binaries to summary DataFrame
   * - ``find_failed``
     - Method reference to manager's find_failed method for identifying failed binary evolutions

During the initialisation, the BinaryPopulation creates a PopulationManager and SImulationPropeteries isntance

Several variables are only set when running with MPI or job arrays, which are used for parallel processing of the population synthesis.
These will not be used in a standard run of the BinaryPopulation class, but are set when using the `posydon-popsyn` or `PopulationRunner` classes.

.. list-table:: Conditional Attributes (MPI/Job Array)
   :widths: 30 70
   :header-rows: 1

   * - Attribute
     - Description
   * - ``comm``
     - MPI communicator object (only when running with MPI)
   * - ``JOB_ID``
     - Job array ID for parallel processing (only when using job arrays)
   * - ``rank``
     - Process rank for parallel processing (set when using MPI or job arrays)
   * - ``size``
     - Total number of processes (set when using MPI or job arrays)



Evolving the population
------------------------

With the population parameters defined, you can evolve the population using the `evolve` method.
This method will 

1. Sample the initial system parameters (single or binary systems).
2. Evolve each binary system through its evolutionary steps.

Here's an example of how to evolve the population:

```python

population.evolve()

```

Additional kwargs can be passed to the `evolve` method to control the evolution process, such as:

.. list-table:: Additional evolve kwargs
   :widths: 30 50 20
   :header-rows: 1

   * - Parameter
     - Description
     - Default
   * - ``indices``
     - Custom binary indices to use instead of range(number_of_binaries). If running with MPI, indices are split between processes
     - None
   * - ``breakdown_to_df``
     - Convert binaries to dataframe and remove from memory after evolution to save RAM
     - True
   * - ``tqdm``
     - Whether to show a progress bar during evolution
     - False
   * - ``from_hdf``
     - Whether to load the population from an HDF5 file instead of evolving it
     - False
   * - ``optimize_ram``
     - Enable RAM optimization by processing binaries in batches. Uses dump_rate from the ini file
     - True
   * - ``ram_per_cpu``
     - Amount of RAM per CPU for batch size calculation (in GB)
     - None
   * - ``temp_directory``
     - Directory path for storing temporary batch files during evolution
     - "batches"


Depending on the parameters, ``evolve`` will create a temporary directory to 
store the batches of binaries during the evolution process.
Within this folder, a batch will write ``dump_rate`` binaries to a temporary file:
``{dump_rate}_evolution.batch``. At the end of the evolution, these files will be merged into a single HDF5 file:
``evolution.combined``.

If you're running with SLURM or MPI, all different processes will write to the 
same folder, with different batch indicators for each process: ``{dump_rate}_evolution.batch.{rank}``.

.. note::

    The merging of different processes is different from the merging of batches within a single process!
    We recommend running this with the `PopulationRunner` class, which will handle 
    the batch writing and merging the output of multiple processes automatically into a single HDF5 file.

When evolving a population, you can read the starting conditions from an HDF5 file
or sample the initial conditions from the given population parameters.
See :ref:`pop-params-guide` for more details about the population parameters file.


Accessing the evolved population
---------------------------------

Depending on the initialisation parameters, the evolved population can be accessed in different ways.

1. If not written to file with ``breakdown_to_df=True``, the population is stored in memory as a list of BinaryStar objects.
   You can access the individual binaries using the `manager` attribute:

   ```python
   first_binary = population.manager.binaries[0]
   print(first_binary)
   ```

   Additionally you can show turn the binary into a history DataFrame or create a oneline summary DataFrame:

   ```python
   history = population.to_df()
   print(history)

   oneline = population.to_oneline_df()
   print(oneline)
   ```

2. If ``breakdown_to_df=True``, the population is removed from memory and written to the population file.
   You can access the population with the normal ``Population`` class.
   Make sure the file name has ``.h5`` extension, as this is required for the Population class to read the file correctly.

   ```python
    from posydon.popsyn.synthetic_population import Population
    population = Population('./population.h5')
    ```


BinaryGenerator class
---------------------

The BinaryGenerator class is a helper class for generating binary systems based on the population parameters.
It can be used to create a population of binaries with specific characteristics, such as mass ratios,
metallicity, and initial conditions.
Please see the :class:`BinaryGenerator <posydon.popsyn.binarypopulation.BinaryGenerator>` documentation for more details on how to use this class.

