# POSYDON: unit tests

Here we collect our unit tests for the POSYDON code. We are using the [python `unittest` module](https://docs.python.org/3/library/unittest.html).

The tests in each `TESTFILE` can be run via

	python -m unittest TESTFILE

If you like to run all unit tests you can use

    python -m unittest $PATH_TO_POSYDON/posydon/unit_tests/*/*.py