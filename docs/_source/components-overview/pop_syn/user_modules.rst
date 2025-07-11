.. _user_modules:

User Modules
------------

To allow users to add code easily within the POSYDON namespace without modifying the code core we added a directory called :samp:`user_modules`. There users can put their own code and simply import it as any other function within the posydon namespace. This feature is designed for users who would like to implement their own common envelope or supernova prescription. Implementing your own module requires essentially two steps. 

First, the class defining your new step needs to be constructed. This class ought to have an :samp:`__init__` method that takes in :samp:`verbose` as an argument and a :samp:`__call__` method that takes in a :samp:`BinaryStar` object as an argument. The physics describing the new step will be included in the :samp:`call` method within which the binary's parameters will be altered according to your prescription. Be sure to adjust the binary :samp:`state` and :samp:`event` so the flow knows how to handle your output binaries.

Second, the binary ``flow`` needs to be adjusted so POSYDON's population synthesis code knows to direct the relevant binaries to the newly constructed step. This can be done by creating a custom user module, similar to how it is done for a custom step. Read our :ref:`flow chart <flow-chart>` documentation to understand what this object is and how you can create your own.

For a complete example of how to implement a new user module, we include a tutorial describing the process: :ref:`/tutorials-examples/population-synthesis/custom_step_and_flow.ipynb`.