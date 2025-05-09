.. _pipeline_additions:

##################
Pipeline additions
##################

All the additions can run after any of the steps and will run in parallel to
the next step. They are run the same way as the :ref:`steps <pipeline_steps>`.

.. _pipeline_plots:

Creating plots
--------------

Plots can be created after each step with the data available from the previous
step. Hence, each corresponding :samp:`csv` file is called 
:samp:`step_?_plots.csv`, where the question mark will be the number of the step 
the plots take the final :samp:`PSyGrid` object from. All the :samp:`csv` files 
for plotting have the same structure:

.. code-block::

    path_to_grid,grid_type,quantities_to_plot,path_to_plot,plot_extension

Beside the grid, it states the type and takes a list of quantities to plot. All
final quantities supported for a :ref:`2D plot <plot_2d>`, as third dimension
can be specified. Additionally, you can put a :samp:`LOG10_` in front of each
of them to switch on plotting in log-scale. Beside that there are predefined
plots. Finally, the path to the directory, where the plots should get stored,
and the extension of the image files (those need to be valid extension for
`mathplotlib <https://matplotlib.org/>`_) are given. There is one additional
extension :samp:`multipage-pdf`, which will create a PDF, where several plots
are stored as pages in a single PDF.

.. table:: Basic predefined plots

    ========================  ====================  ========================  ======  ======  ======
    quantities_to_plot        'term_flag'           'zvar'                    'zmin'  'zmax'  'zlog'
    ========================  ====================  ========================  ======  ======  ======
    'combined_TF12'           'combined_TF12'       None                      None    None    False
    'termination_flag_1'      'termination_flag_1'  'lg_mtransfer_rate'       -8      -1      False
    'termination_flag_2'      'termination_flag_2'  None                      None    None    False
    'termination_flag_3'      'termination_flag_3'  None                      None    None    False
    'termination_flag_4'      'termination_flag_4'  None                      None    None    False
    'rl_relative_overflow_1'  'debug'               'rl_relative_overflow_1'  -0.5    0.5     False
    'rl_relative_overflow_2'  'debug'               'rl_relative_overflow_2'  -0.5    0.5     False
    'lg_mtransfer_rate'       'debug'               'lg_mtransfer_rate'       -8      -1      False
    ========================  ====================  ========================  ======  ======  ======

After :ref:`Step3: calculating extra values from detailed data 
<pipeline_step3>`, the supernova model quantities get available, too.

.. table:: Supernova predefined plots

    ============================  ====================  ============================  ======  ======  ======
    quantities_to_plot            'term_flag'           'zvar'                        'zmin'  'zmax'  'zlog'
    ============================  ====================  ============================  ======  ======  ======
    'S1_MODEL??_CO_type'          'S1_MODEL01_CO_type'  'S1_MODEL??_CO_type'          None    None    False
    'S1_MODEL??_SN_type'          'S1_MODEL01_SN_type'  'S1_MODEL??_SN_type'          None    None    False
    'S1_MODEL??_mass'             'termination_flag_1'  'S1_MODEL??_mass'             0.      2.      True
    'S1_MODEL??_spin'             'termination_flag_1'  'S1_MODEL??_spin'             0.      1.      False
    'S1_MODEL??_m_disk_radiated'  'termination_flag_1'  'S1_MODEL??_m_disk_radiated'  0.      3.      False
    ============================  ====================  ============================  ======  ======  ======

After :ref:`Step4: training of the interpolators <pipeline_step4>`, the
interpolators can be used.

.. table:: Predefined plots with variable `QUANTITY`

    ============================  ====================  ============  ======  ======  ======
    quantities_to_plot            'term_flag'           'zvar'        'zmin'  'zmax'  'zlog'
    ============================  ====================  ============  ======  ======  ======
    '`QUANTITY`'                  'termination_flag_1'  '`QUANTITY`'  None    None    False
    'LOG10\_ `QUANTITY`'          'termination_flag_1'  '`QUANTITY`'  None    None    True
    'INTERP\_ERROR\_ `QUANTITY`'  None                  '`QUANTITY`'  0       0.1     False
    ============================  ====================  ============  ======  ======  ======

.. _pipeline_checks:

Doing checks
------------

After each step one can perform some checks with the data from that step.

.. table:: Currently supported checks

    ==============  ===========
    Check           Description
    ==============  ===========
    'failure_rate'  calculates the failure rate of the grid
    'CO_type'       gets counts of compact object types
    'SN_type'       gets counts of supernova types
    ==============  ===========

