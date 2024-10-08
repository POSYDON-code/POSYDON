.. _code-usage:

Code Usage Questions
--------------------

Welcome to the code usage FAQ for POSYDON. This section is designed to answer questions related to how to use various features and functionalities of POSYDON. If your question isn't addressed here, consider reaching out to our :ref:`support channels <contact_info>` or checking our :ref:`general FAQ <installation-issues>`.

Frequently Asked Questions
~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **How do I start a basic population synthesis simulation with POSYDON?**
    - Answer: Start with the :ref:`binary population synthesis guide <binary-pop-syn>` to learn how to run POSYDON.

2. **I'm getting an error when using the** ``where`` **parameter in the** ``select`` **function**

    If you're trying to select based on values in an unsupported column, you'll receive an error message like below.
    Only specific columns are supported for selection. These are limited to string columns, the indices, and the column names.
    Please perform the selection in-memory instead.

    .. code-block:: python

        population.history.select('S1_mass > 10')

    The above code will produce the following error:
    .. code-block:: python
        
        ValueError: The passed where expression: S1_mass > 10
            contains an invalid variable reference
            all of the variable references must be a reference to
            an axis (e.g. 'index' or 'columns'), or a data_column
            The currently defined references are: index,columns,state,event,step_names,S1_state,S2_state


3. **How should I tune my memory usage for a population synthesis run?**

    A population run requires at a bare minimum of 4GB of memory per CPU.
    However, this restricts the number of binaries you can keep in memory and requires a :code:`dump_rate < 1000` to keep the memory usage low, which slows down the simulation.
    
    As such, 5GB per CPU is a better starting point. This allows you to keep more binaries in memory and increases the speed of the population synthesis run.
    The figure below can be used to fine-tune your memory usage.
    It shows the memory usage of two large population synthesis runs with a 7GB and a 8GB limit with a :code:`dump_rate` of 8000 and 10.000 binaries respectively.

    .. image:: ./large_pop_runs_memory.png
        :alt: Memory usage of two large population synthesis runs.
        :width: 700px
    
    The memory usage of the 8000 :code:`dump_rate` run is stable at around 6GB, while the 10.000 :code:`dump_rate` run is stable at around 6.8GB.


4. **What should the walltime and job array size be for my population synthesis run?**

    The :code:`walltime`` and job array size are dependent on the number of binaries you want to simulate and the memory usage of the simulation.
    The job array size should be set such that the number of binaries per job is at least 1000, since there's a minimum overhead per job due to loading the grids.
    
    The :code:`walltime` depends on the number of binaries per job, where each binary takes about 1-2 seconds to run.
    For example, with 100.000 binaries split over 100 jobs (per metallicity), means that every job runs 1.000 binaries. This will take around 33 minutes per job. So a :code:`walltime` of `00:45:00` is reasonable.

    The balance between :code:`walltime` and the size of the job array is important.
    If the :code:`walltime` is too long, it might be worth increasing the job array size to decrease the time per job and allowing the population synthesis to finish faster. 
    But if the :code:`walltime` is too short, the job array size should be decreased, since each job has an initial overhead that is not dependent on the number of binaries in the job.

    .. note::
        The processing time increases if you make the `dump_rate` too low due to many I/O operations.

2. **Which parameters can I customize for my simulations?**
    - Answer: POSYDON allows customization of ... (TODO list or briefly describe customizable parameters).

3. **I'm getting an error message when trying to run a simulation. What does it mean?**
    - Answer: Error messages typically indicate ... (TODO provide guidance on interpreting common error messages or point to a troubleshooting guide).

4. **How can I visualize the results from my simulation?**
    - Answer: You can utilize the experimental visualization libraries by ... (TODO provide a brief overview or link to visualization documentation).

5. **How do I use the new machine learning features in POSYDON?**
    - Answer: To leverage the ML features, ensure you've installed the necessary dependencies. Then, you can ... (TODO brief overview or link to ML features documentation).

6. **Are there any examples or tutorials available?**
    - Answer: Yes, you can check our :ref:`roadmap <roadmap>` for tutorials related to different POSYDON components, including population synthesis, creating core datasets, and running your own MESA grids with POSYDON.

7. **Can I run POSYDON on an HPC facility?**
    - Answer: Absolutely! Ensure you've installed ``mpi4py`` as outlined in :ref:`our installation guide <installation-guide>`. Refer to `our HPC guide <../tutorials-examples/population-synthesis/pop_syn.ipynb>`_ for detailed instructions on running POSYDON in an HPC environment.

8. **How can I stay updated with the latest features and updates?**
    - Answer: You can regularly visit our `official website <https://posydon.org>`_ for news and updates. Also, consider subscribing to our mailing list.

Additional Resources
~~~~~~~~~~~~~~~~~~~~

1. **User Guide**: For detailed instructions on all features of POSYDON, visit our comprehensive :ref:`roadmap <roadmap>`.
 
2. **API Reference**: Dive deep into the functionality provided by POSYDON with our :ref:`API Reference <modules>`.

3. **Examples and Tutorials**: Learn by doing! Visit :ref:`our roadmap page <roadmap>` for hands-on learning.

Still Have Questions?
~~~~~~~~~~~~~~~~~~~~~

If your query remains unanswered, we're here to help! Reach out to our community through the :ref:`support channels <contact_info>` or consider checking our :ref:`general installation FAQ <installation-issues>` for non-usage related questions.

Your feedback helps us improve. If you think a common question should be added here, don't hesitate to suggest it!
