.. _pipeline:

########################
Post-processing pipeline
########################

The post-processing is done in steps with the possibility to get things run in
parallel after each step. The pipeline_steps describe each step in detail, while
the pipeline_additions describe plots and checks that can be done after a
step.

.. toctree::

   pipeline/pipeline_steps

.. toctree::

   pipeline/pipeline_additions

We recommend running the pipeline with an :ref:`ini file <pipeline_ini>` to keep
a better overview of large grids and being able to post-process them with a few
commands.
