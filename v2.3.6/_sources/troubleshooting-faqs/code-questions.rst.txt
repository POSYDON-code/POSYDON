.. _code-usage:

Code Usage Questions
--------------------

Welcome to the code usage FAQ for POSYDON. This section is designed to answer questions related to how to use various features and functionalities of POSYDON. If your question isn't addressed here, consider reaching out to our :ref:`support channels <contact_info>` or checking our :ref:`general FAQ <installation-issues>`.

Frequently Asked Questions
~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **How do I start a basic population synthesis simulation with POSYDON?**
    Start with the :ref:`binary population synthesis guide <binary-pop-syn>` to learn how to run POSYDON.

2. **I'm getting an error when using the** ``where`` **parameter in the** ``select`` **function**
    If you're trying to select based on values in an unsupported column, you'll receive an error message like the one below.
    Only specific columns are supported for selection. These are limited to string columns, the indices, and the column names.
    Please perform the selection in-memory instead.

    .. code-block:: python

        population.history.select('S1_mass > 10')

    The above code will produce the following error:

    .. code-block::

        ValueError: The passed where expression S1_mass > 10
        contains an invalid variable reference. All of the variable references
        must be a reference to an axis (e.g. 'index' or 'columns'), or a data_column.
        The currently defined references are: index, columns, state, event, step_names, S1_state, S2_state


3. **How should I tune my memory usage for a population synthesis run?**
    A population run requires a bare minimum of 8GB of memory per CPU. The DR2 POSYDON grids need to be loaded in completely. At the moment, the binary grids are loaded in at startup, but single star models are loaded when needed.
    However, 8GB allows you to run a small population, but it restricts the number of binaries you can keep in memory and requires a :code:`dump_rate < 1000` to keep the memory usage low, which slows down the simulation.

    As such, 9GB or 10GB per CPU is a better starting point. This allows you to keep more binaries in memory and increases the speed of the population synthesis run and work well with a :code:`dump_rate` of 5.000 binaries.

    In general, the :code:`dump_rate` should be at least 500 for populations of 100,000 binaries or more.
    Setting a very low :code:`dump_rate` for larger populations will cause many I/O operations to be called (we are dumping out and saving our binary star data more frequently) and can potentially introduce issues during the reading, writing,
    and merging of output files. Keeping the ``dump_rate`` at a relatively high number, but not so much that it causes excessive memory usage will allow your population to run faster by limiting I/O operations. (However, of course setting a
    ``dump_rate`` that is greater than the number of binaries in your population will effectivley counteract the purpose of this parameter, as it means you would only dump your binaries into a file after you have run all of them.)


4. **What should the walltime and job array size be for my population synthesis run?**
    The :code:`walltime` and job array size are dependent on the number of binaries you want to simulate and the memory usage of the simulation.
    The job array size should be set such that the number of binaries per job is at least 1000, since there's a minimum overhead per job due to loading the grids.

    The :code:`walltime` depends on the number of binaries per job, where each binary takes about 1-2 seconds to run.
    For example, 100,000 binaries split over 100 jobs (per metallicity) means that every job runs 1,000 binaries. This will take around 33 minutes per job. So a :code:`walltime` of :code:`00:45:00` is reasonable.

    The balance between :code:`walltime` and the size of the job array is important.
    If the :code:`walltime` is too long, it might be worth increasing the job array size to decrease the time per job and allow the population synthesis to finish faster.
    But if the :code:`walltime` is too short, the job array size should be decreased, since each job has an initial overhead that is not dependent on the number of binaries in the job.

    .. note::
        The processing time increases if you make the :code:`dump_rate` too low due to many I/O operations.

5. **I am unable to open HDF5 files created by POSYDON. What should I do?**
    If you're on a Mac, there might be an issue with the HDF5 installation.
    Make sure you have the :code:`hdf5` and :code:`pytables` packages installed through conda in your environment with :code:`conda install hdf5 pytables` before running POSYDON!
    Although they are dependencies of POSYDON, sometimes they are not installed correctly on Mac.

6. **Can I run POSYDON on an HPC facility?**
    Absolutely! Refer to `our HPC guide <../tutorials-examples/population-synthesis/pop_syn.ipynb>`_ for detailed instructions on running POSYDON in an HPC environment.

7. **Help, I'm stuck! Where can I get support?**
    Please check if your question has been answered already on the the POSYDON `Github Discussions <https://github.com/POSYDON-code/POSYDON/discussions>`_.
    If it has not yet been answered, please make a new discussion with your question and a detailed description of your problem.

    Are you unable to post on Github? Please check `our email group <https://groups.google.com/g/posydon-users>`_ or email us at posydon-users [at] googlegroups.com

8. **How can I stay updated with the latest features and updates?**
    You can regularly visit our `official website <https://posydon.org>`_ for news and updates.

9. **I've come across a FAILED binary. What does this mean?**
     A :code:`FAILED` binary has encountered an error during the simulation because POSYDON was unable to evolve it. This can be due to a variety of reasons:

        -  The evolutionary state of the binary is not represented in the currently supported stellar evolution grids. For example, we do not have a grid for Roche lobe overflow between two helium stars.
        -  The binary has masses outside the grid range. For example, the HMS-HMS grid does not contain binaries with a secondary mass below 0.5.
        -  The binary could not be matched to a single star or a binary due to a too large matching error, preventing further evolution.

10. **What approximations does POSYDON make?**
     This is a complex question and the best answer is provided in the POSYDON papers: `Fragos et al. (2023) <https://ui.adsabs.harvard.edu/abs/2023ApJS..264...45F/abstract>`_ and `Andrews et al. (submitted) <https://ui.adsabs.harvard.edu/abs/2024arXiv241102376A/abstract>`_.


Additional Resources
~~~~~~~~~~~~~~~~~~~~
1. **Examples and Tutorials**: Learn by doing!
2. **API Reference**: Dive deep into the functionality provided by POSYDON with our :ref:`API Reference <modules>`.
3. **Github Discussions**: Engage with the community, ask questions, and share your experiences on our `GitHub Discussions <https://github.com/POSYDON-code/POSYDON/discussions>`_ page.

Still Have Questions?
~~~~~~~~~~~~~~~~~~~~~

If your question remains unanswered, we're here to help! Reach out to our community through the :ref:`support channels <contact_info>` or consider checking our :ref:`general installation FAQ <installation-issues>` for non-usage related questions.

Your feedback helps us improve the code and documentation. If you think a common question should be added here, don't hesitate to suggest it!
