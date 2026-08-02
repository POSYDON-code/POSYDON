.. _add-step-user-module:

Adding a Custom Step or User Module
====================================

POSYDON is designed so that you can add your own evolutionary steps (a
custom supernova or common-envelope prescription) and your own flow charts
without touching the core code. This page summarizes how to implement one of
these custom pieces and where to find the full worked example.

Two kinds of customization
---------------------------

There is a useful distinction to keep in mind:

* a **custom step** is a piece of physics that is inserted into the evolution
  of a binary. It is a class that gets called on the ``BinaryStar`` when the
  flow chart routes that binary to it.
* a custom **flow chart** (or "flow") is the routing logic that decides, given
  a binary's state and event, which step runs next.

If your new physics needs its own step, you will generally write *both a
custom step and* a custom flow (or modify the existing flow) so the binaries
are routed to it.

Implementing a custom step
------------------------------------

Place your module in the :ref:`user_modules <user_modules>` directory so you
can import it inside the POSYDON namespace without changing the core code.
Your step must define a class with:

* an ``__init__`` method that takes ``verbose`` as an argument; and
* a ``__call__`` method that takes a ``BinaryStar`` object as an argument.

The physics lives in the ``__call__`` method, which should alter the binary's
properties according to your prescription. Be sure to update the binary's
``state`` and ``event`` attributes so the flow knows how to handle the
outgoing binaries.

Adjusting the flow
----------------------

Next, update the ``flow`` so that POSYDON directs the relevant binaries to your
new step. Because binary ``flow`` is a dictionary describing the state machine,
you can create a custom flow for your model. See the
:ref:`flow chart <flow-chart>` documentation to understand this object and how
to build your own.

Full tutorial
------------------

If you would like to see a complete example, follow the notebook *Creating
custom steps and flow chart for the evolution of a binary*:
``tutorials-examples/population-synthesis/custom_step_and_flow.ipynb``. The
companion ``custom_step_and_flow.py`` shows the implementation as a script.

Related
--------

* :ref:`user_modules` -- the directory where custom code is placed.
* :ref:`flow-chart` -- the routing object and how to build one.
