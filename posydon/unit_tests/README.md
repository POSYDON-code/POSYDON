# POSYDON: unit tests

Here we collect our unit tests for the POSYDON code. We are using the [pytest package](https://docs.pytest.org/en/stable/index.html).

Before using it you need to install pytest or consider to upgrade it

    pip install -U pytest

The tests in each `TESTFILE` can be run via

    pytest TESTFILE

If you like to run all unit tests you can use

    pytest $PATH_TO_POSYDON/posydon/unit_tests/

## File structure

As already visible from the commands how to run all unit tests we have a the `unit_tests` directory. In there are individual directories to combine several tests into one suite. It is a natural choice to reuse the topics used in the code itself. For example all the unit tests for files inside `$PATH_TO_POSYDON/posydon/utils` are put in `$PATH_TO_POSYDON/posydon/unit_tests/utils`. In those suites, there should be a file for each code file, which is simply extended by the prefix `test_`, e.g. `test_posydonerror.py` contains the unit tests of `posydonerror.py`. If the code file is in a subdirectory, may add this to the prefix. *I'd discourage to create deep directory structures inside the `unit_tests`.*

## How to add a new unit test

There is a template file `$PATH_TO_POSYDON/posydon/unit_tests/test_template.py`, which can be copied. Please rename/replace the copy to fulfill the [File structure](#file-structure). Please update the doc-string, the author list, the imported modules, and the test classes.

### Conventions

To keep things in a common style, here are some conventions to follow. You can never make too many tests as long as they are unique!

#### Doc string

Each test file should have a doc-string. This should at least contain the information which part of the code will be tested.

#### Author list

Like in code files we add an author list to those files to let others know whom to ask when doing changes in such a file.

#### Imports

We always need to import the module we'd like to test. *I recommend to import it with the local name `totest`.*

There might be additional modules/functions needed to do the testing. *I recommend to not import stuff, which is already imported in the module, you like to test, instead use the reference in the tested module, e.g. if the module to test imports X, access it via `totest.X` in the unit test, to ensure using the same code, which is not subject to this test. For convenience you might define some aliases for imported stuff, e.g. `np = totest.np`.*

#### Test functions

For a small single test you can simply define a single test function. The name should start with the prefix `test`. *I recommend to use the prefix `test_`.* The actual test is usually done by an assert statement. Usually, there will be several tests collected in a [class](#test-classes), thus single functions will be rare.

#### Test classes

To group tests it is useful to create classes which contain several test functions. All test classes should have the prefix `Test`. *I recommend to continue with a capital letter after the prefix.* In each class there should be test functions, again having a prefix `test`. *I recommend to use the prefix `test_`.* It should be noted, that the test functions inside a class require a parameter containing the class object. *I recommend to use the standard `self`.* It should be noted, that there is no guarantee for the order the tests are run in. Thus, each test should be considered a stand alone.

##### Check the elements of a module

The first thing to check is the completeness and types of the module elements. *I suggest to use the build-in function `dir`.* All variables should be tested for their type via `isinstance`. All classes and functions can be tested for their type by making us of `isclass` and `isroutine` form the `inspect` module, which should return `True`. *I suggest to put all those checks into one test class with functions for each element.*

##### Check values

Beside the existence and the type of a variable, we should verify the integrity of it's default value. *I suggest to put all those value checks into one test class with a function for each variable.*  Depending on the kind of variable different tests can be done, e.g. a fixed value for a constant, a range of a changeable variable, or needed items of a list ...

##### Check functions

Functions in the module need checks according to their functionality. This should coincide with the doc-string of each function. *I suggest to have one class with all the function tests and a test function for each function in the module, which gets tested.* If the functions need variables, you may [use fixtures](#using-fixtures).

Functions may include prints statements. To capture outputs to `stdout` and `stderr` pytest has global defined fixtures `capsys` (and `capsysbinary` for bytes data instead of usual text; addtionally, there are `capfd` and `capfdbinary` to capture the file descriptors `1` and `2`).

To access the captured content you simply call the `readouterr` function of the fixture. It will return a `(out, err)` namedtuple. Here an example how to get the stdout content, captured so far:

    capsys.readouterr().out

It should be noted, that the call of `readouterr` will clear the buffer, hence if you like to get both `out` and `err` at the same time, you need to store the namedtuple and access the components from the stored version, e.g.

    captured_output = capsys.readouterr()
    captured_output.out
    captured_output.err

In case you like to access what is captured to far but keep it in the buffer, you'd need to reprint the part you read, e.g.

    captured_output = capsys.readouterr()
    captured_output.out
    print(captured_output.out)

Using those fixtures as arguments to a function will capture all outputs of that function. To exclude some parts from getting its output captured you need to put that into the `disabled` context of the used fixture. It should be noted, that only the capturing is disabled, but the fixture object is still available, e.g.

    with capsys.disabled():
        print(capsys.readouterr().out)

will print the collected prints to stdout, at this moment all together (and clears the buffer).

##### Check classes

Each class inside a module should be get its components checked like a module itself. *I suggest to have a test class for each class in the tested module and the test of each class function should get an own test function.* Again, [fixtures](#using-fixtures) can be used to ensure that all tests run under the same conditions. A commonly useful fixture for a class test is an object of this class initialized with the defaults.

#### Using fixtures

You can define [fixtures](https://docs.pytest.org/en/stable/how-to/fixtures.html) at any level. You need to decorate them with the [`pytest.fixture` function](https://docs.pytest.org/en/stable/reference/reference.html#pytest-fixture), which needs to be imported, via

    @pytest.fixture

*I suggest to import the function via `from pytest import fixture` and than call it via `@fixture`.*

Fixtures replace the `setUp` and `tearDown`. To use a fixture to prepare something before a test, you can simply write it as a function and the variable will contain the returned value.

For cleaning things up after a test, instead of having a final return, you separate setUp and tearDown with a `yield` statement, which ensure that all before is executed when the fixture is requested and the stuff after when it get deleted (usually at the end of the test function). For chains of fixtures it should be noted, that the clean up happens in the reverse order to the creation, because the innermost fixture will get deleted first.

There are some useful [predefined fixtures](https://docs.pytest.org/en/stable/reference/fixtures.html). One was already introduced earlier, `capsys` can be used to interact with captured output. Another very useful one is `tmp_path`, which is a path object to a directory, which is created for each test function and will be subject to removal afterwards. Hence, whenever you need to read/write files for the test, it should be done in there. In most cases, `tmp_path` will not be used in your test function, but you use it in a self-defined fixture, which e.g. creates a file for the function. A third predefined fixture of interest is `monkeypatch`. This can be used to replace the objects inside the tested code just for the test. More [predefined fixtures can be found online](https://docs.pytest.org/en/stable/reference/fixtures.html).

#### Catching raised errors

Pytest has the [context manager `pytest.raises`](https://docs.pytest.org/en/stable/reference/reference.html#pytest-raises) to catch raised errors. You use it like other context managers via a `with` statement. Beside the expected exception, you can specify a `match`, which will be checked against the error message as a regular expression, e.g.:

    with pytest.raises(TypeError, match="Error message"):
        raise TypeError("Error message")

The context object can be used to check for more details of the raised error (this is useful to overcome some limitations of escape characters in regular expressions). Here an example how to check the error message from the context object:

    with pytest.raises(TypeError) as error_info:
        raise TypeError("Error message")
    assert error_info.value.args[0] == "Error message"

It should be noted that the error type will match subclasses successfully.

#### Catching warnings

Usually, pytest will catch all warnings and print them at the end of all tests. If your test will cause a warning which you don't like to have displayed, you can filter the warnings caught by pytest. To filter all warnings in a function or class you can decorate it with a filter, e.g. `@pytest.mark.filterwarnings("ignore:WARNINGTEXT")`. There are more things you can do on [warnings in pytest](https://docs.pytest.org/en/stable/how-to/capture-warnings.html), but you should use that only were needed. But you should be careful with the pytest warning catching, because it overwrites some parts of the python warnings, which even interferes badly with our POSYDON warnings (especially the filter changes). By using the `pytest.warns` context you can capture and check for warnings the same way as for [errors](#catching-raised-errors).

### Do not test functions called inside or need for a function

In more complex cases, the function the test acts on will call other functions. Here it can be useful to replace inner functions with `monkeypatch`.

Some functions may require input generated from other functions, here it can be useful to call the other function within a try and return in case it fails to let the test fail on the function, which makes problems not a function, which needs its output. *I recommened therefore to first write a unit test for the functions which do not other functions.*

### Check that it can fail

Whenever you write a test, it should succeed on the existing code. On the other hand, you should also try out introducing a temporary change that you expect to lead to a failure to ensure the test is doing its job.

## How to update a unit test

First, check what is already there. Second, edit or add new stuff following the [conventions](#conventions) and test that it works by let it fail once.

## Check code coverage

There is a plug-in for `pytest` to measure the coverage of executed code, called [pytest-cov](https://pypi.org/project/pytest-cov). It can be installed via

    pip install pytest-cov

and run by including the `cov` option in the `pytest` call

    pytest TESTFILE --cov=MODULE

Strictly speaking it runs the `pytest` inside of [coverage](https://coverage.readthedocs.io). Thus you can run it although via

    coverage run -m pytest TESTFILE

*I suggest to use the `--branch` option, which is `--cov-branch` in pytest.* The latter version to run the test has the advantage that you can make use of all the options provided by `coverage`, while running it through the plug-in will give you only access to the options implemented in there. On the other hand, running it with the plug-in will exclude the test code from coverage and give you the coverage report together with the report of the tests, while running it via `coverage` the report in the file `.coverage` needs to readout via

    coverage report

*I suggest to use the `-m` option to show uncovered lines, which is `--cov-report term-missing` in pytest.*

We should aim for 100% coverage. If there is code which should be excluded from the coverage, please mark it with `# pragma: no cover`. The main usage is for code, which under normal conditions should never run, e.g. in POSYDONwarnings is a part to act, when the python `sys` module fails.

You can also enforce 100% coverage in the Github action performing the tests by adding `--cov-fail-under=100` to the pytest command, for example:

    python -m pytest posydon/unit_tests/ --cov=posydon.utils --cov-branch --cov-report term-missing --cov-fail-under=100 

This particular line ensures that the tests in `posydon/unit_tests/` run 100% of the code in `posydon.utils`, and that the check will otherwise fail.

