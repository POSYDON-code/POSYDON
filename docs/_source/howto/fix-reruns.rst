.. _fix-reruns:

Fixing and Rerunning FAILED binaries
====================================

You will inevitably encounter ``FAILED`` binaries in a population run. This page
explains what a ``FAILED`` binary is, why it happens, and how to rerun the
affected subset rather than repeating the whole run.

What does FAILED mean?
----------------------
A binary is marked ``FAILED`` when POSYDON was unable to evolve it. There are
three common reasons:

- The evolutionary state of the binary is not represented in the currently
  supported stellar grids. For example, there is no grid for Roche-lobe
  overflow between two helium stars.
- The binary has masses outside the range of the grid. For example, the
  ``HMS-HMS`` grid does not contain binaries with a secondary mass below 0.5
  solar masses.
- The binary could not be matched to a single star or a binary with a small
  enough matching error, preventing further evolution.

A FAILED binary is not necessarily a bug in your run; it usually means the
binary lies outside the coverage or assumptions of the supplied grids.

How to rerun a subset
---------------------

Rather than restarting the whole population, you can rerun just the binaries
you are interested in. Two companion concepts in the code base support this:

* ``PSyGrid.rerun`` exports selected runs of a grid so they can be re-run with a
  new MESA setup. You can either write your own logic that picks the indices you
  want, or select by termination flag. See also the full :ref:`psygrid` page.

* For a population run, the *step rerun* workflow lets you redo a particular
  step rather than the full pipeline. See the ``just_step_1`` and
  ``step_rerun`` notebooks under ``tutorials-examples/generating-datasets/``.

When the failure was just a transient numerical issue or a bad set of MESA
parameters, changing the setup slightly and rerunning the flagged subset will
often turn a FAILED binary into a successful one.

Related
---------------

* :ref:`psygrid` for operating on ``PSyGrid`` data and the ``rerun`` function.
* :ref:`processing-pipeline` for how reruns fit into the post-processing
  workflow.
