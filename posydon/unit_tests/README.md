# POSYDON: unit tests

Here we collect our unit tests for the POSYDON code. We are using the [python `unittest` module](https://docs.python.org/3/library/unittest.html).

The tests in each `TESTFILE` can be run via

    python -m unittest TESTFILE

If you like to run all unit tests you can use

    python -m unittest $PATH_TO_POSYDON/posydon/unit_tests/*/*.py

## File structure

As already visible from the commands how to run all unit tests we have in the `unit_tests` directory individual directories to combine several tests into one suite. It is a natural choice to reuse the topics used in the code itself. For example all the unit tests for files inside `$PATH_TO_POSYDON/posydon/utils` are put in `$PATH_TO_POSYDON/posydon/unit_tests/utils`. In those suits, there should be a file for each code file, which is simply extended by the prefix `test_`, e.g. `test_posydonerror.py` contains the unit tests of `posydonerror.py`. If the code file is in a subdirectory, may add this to the prefix. *I'd discourage to create deep directory structures inside the `unit_tests`.*

## How to add a new unit test

There is a template file `$PATH_TO_POSYDON/posydon/unit_tests/test_template.py`, which can be copied. Please rename/replace the copy to fulfill the [File structure](#file-structure). Please update the doc-string, the author list, the imported modules, and the test classes.

### Conventions

To keep things in a common style, here are some conventions to follow. You can never make too many tests as long as they are unique!

#### Doc string

Each test file should have a doc-string. This should at least contain the information which part of the code will be tested.

#### Author list

Like in code files we add an author list to this files to let others know whom to ask when doing changes in this file.

#### Imports

We always need to import the `unittest` module. Beside that we always need to import the module we'd like to test. *I recommend to import it with the local name `totest`.*

There might be additional modules/functions needed to do the testing. *I recommend to not import stuff, which is already imported in the module, you like to test, instead use the reference in the tested module, e.g. if the module to test imports X, access it via `totest.X` in the unit test, to ensure using the same code, which is not subject to this test.*

#### Test classes

According to [python `unittest` module](https://docs.python.org/3/library/unittest.html), there are conventions to follow, when creating test classes and their functions. It should be noted, that there is no guarantee for the order the tests are run in. Thus, each test should be considered a stand alone.

#### Check the elements of a module

The first thing to check is the completeness and types of the module elements. *I suggest to use the build-in function `dir`.* All variables should be tested for their type via `assertIsInstance`. All classes and functions can be tested for their type by making us of `isclass` and `isroutine` form the `inspect` module, which should return `True`. *I suggest to put all those checks into one test class with functions for each element.*

#### Check values

Beside the existence and the type of a variable, should verify the integrity of it's default value. *I suggest to put all those value checks into one test class with a function for each variable.*  Depending on the kind of variable different tests can be done, e.g. a fixed value for a constant, a range of a changeable variable, or needed items of a list ...

#### Check functions

Function in the module need checks according to their functionality. This should coincide with the doc-string of each function. *I suggest to have one class with all the function tests and a test function for each function in the module, which gets tested.* If the functions change variables, please create a `setUp` and/or `tearDown` function.

Functions may include prints statements. In such cases it is useful to redirect the data stream going there into a variable to be able to validate the output (there is a similar context manager to redirect `stderr`). Here an example code:

    from io import StringIO
    from contextlib import redirect_stdout
    ...
    with redirect_stdout(StringIO()) as print_out:
        totest.function_with_print_statements()
    self.assertEqual("Text of first print statement.\nText of second print statement.\n", print_out.getvalue())

#### Check classes

Each class inside a module should be get its components checked like a module itself. *I suggest to have a test class for each class in the tested module and the test of each class function should get an own test function.* Again, `setUp` and/or `tearDown` functions should be used to ensure that all tests run under the same conditions.

### Calling the unit test

Each test file should include the following lines at the end to allow the to run it as standalone:

    if __name__ == "__main__":
        unittest.main()

## How to update a unit test

First, check what is already there. Second, edit or add new stuff following the [conventions](#conventions).