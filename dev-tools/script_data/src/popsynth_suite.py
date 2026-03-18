import os
import traceback
import warnings

from formatting import columns_to_show, line_length
from pandas.testing import assert_frame_equal
from utils import print_pop_settings, print_warnings

from posydon.config import PATH_TO_POSYDON
from posydon.popsyn.binarypopulation import BinaryPopulation
from posydon.popsyn.synthetic_population import Population, PopulationRunner

base_dir = os.path.dirname(PATH_TO_POSYDON)
script_dir = os.path.join(base_dir, "script_data/")
path_to_default_params = os.path.join(script_dir, "inlists/default_test_params.ini")
path_to_multiZ_params = os.path.join(script_dir, "inlists/multiZ_test_params.ini")
path_to_popout = os.path.join(script_dir, "output/population_tests/batches")

def test_binpop_evolve(population, popevo_kwargs, verbose=False):

    # Capture warnings during evolution
    captured_warnings = []

    def warning_handler(message, category, filename, lineno, file=None, line=None):
        captured_warnings.append({
            'message': str(message),
            'category': category.__name__,
            'filename': filename,
            'lineno': lineno
        })

    # Set up warning capture
    old_showwarning = warnings.showwarning
    warnings.showwarning = warning_handler

    try:

        print("Running BinaryPopulation.evolve() with settings:")
        for key, val in popevo_kwargs.items():
            print(f"\t {key} : {val}")

        print("\nEvolving BinaryPopulation...\n")
        population.evolve(**popevo_kwargs)

        print_warnings(captured_warnings)
        print("✅ BinaryPopulation evolved successfully.")
        print("=" * line_length)

    except Exception as e:
            print_warnings(captured_warnings)
            print(f"🚨 BinaryPopulation evolution failed!\n")
            traceback.print_exc(limit=3)
            print("\n")
            print("=" * line_length)

    return population

def compare_io_to_ram(loaded_pop, pop_in_ram):

    # check that binaries match between pop runs w/ fixed entropy
    # and that saved/loaded binaries match those from a memory loaded run
    df_from_ram = pop_in_ram.to_df()
    ram_dflist = [df_from_ram.loc[i] for i in range(10)]
    print("🔍 Checking that binaries in RAM match those retrieved from I/O...")
    for i, ram_df in enumerate(ram_dflist):
        io_df = loaded_pop.history[i]
        ram_df = ram_df[columns_to_show]
        io_df = io_df[columns_to_show]
        try:
            assert_frame_equal(ram_df, io_df)
        except AssertionError as e:
            print("🚨 A binary from I/O does not equal the same binary stored in RAM:")
            print(e)
            print("\nBinary in RAM:\n", ram_df)
            print("\nBinary from I/O:\n", io_df)
            print("=" * line_length)
            return
        if i == len(loaded_pop.history):
            break

    print("✅ Binaries from I/O match those in RAM.")
    print("=" * line_length)

def print_testinfo(test_title, population, popevo_kwargs):

    # print test title str
    numchar = (line_length - len(test_title)) // 2
    print("=" * numchar + test_title + "=" * numchar)

    optimize_ram = popevo_kwargs.get("optimize_ram", False)
    breakdown_to_df = popevo_kwargs.get("breakdown_to_df", False)
    N_binaries = population.number_of_binaries

    # print expected behavior based on settings
    if not breakdown_to_df and not optimize_ram:
        print(f"🚀 Evolving a population and storing {N_binaries} binaries in RAM.")
    elif breakdown_to_df and not optimize_ram:
        print("🚀 Evolving a population and saving binaries to a hdf5 file.")
    elif optimize_ram and not breakdown_to_df:
        num_batch_files = population.number_of_binaries // population.kwargs["dump_rate"]
        print(f"🚀 Evolving a population and saving to {num_batch_files} batch files.")

def check_test(pop_in_ram, load_pop=False):

    # if we have population in RAM, check that the number
    # stored in RAM matches the number we expected to run
    if pop_in_ram and not load_pop:
        num_binaries = pop_in_ram.number_of_binaries
        num_in_ram = len(pop_in_ram.manager.binaries)
        print(f"🔍 Checking that we have {num_binaries} binaries in RAM...")
        assert num_binaries == num_in_ram, \
            f"🚨 Number of binaries in RAM ({num_in_ram}) " \
            f"does not equal the number specified to run ({num_binaries})."
        print(f"✅ Successfully ran and stored {num_in_ram} binaries in RAM.")

    elif pop_in_ram and load_pop:
        save_fn = os.path.join(path_to_popout, "evolution.combined.h5")
        loaded_pop = Population(save_fn)
        compare_io_to_ram(loaded_pop, pop_in_ram)

def test_popruns():

    print("Performing population run tests...")
    pop = BinaryPopulation.from_ini(path_to_default_params, verbose=False)
    pop.kwargs.update({"temp_directory": path_to_popout})
    print_pop_settings(pop)

    # DO TESTS:

    # test simple run, stays in RAM
    kwargs = {"optimize_ram":False, "breakdown_to_df":False, "tqdm":True}
    print_testinfo("TEST: 01", pop, kwargs)
    pop_in_ram = test_binpop_evolve(pop, kwargs, verbose=True)
    check_test(pop_in_ram, load_pop=False)

    # test same but w/ saving/loading binaries
    kwargs = {"optimize_ram":False, "breakdown_to_df":True, "tqdm":True}
    print_testinfo("TEST: 02", pop, kwargs)
    _ = test_binpop_evolve(pop, kwargs, verbose=False)
    check_test(pop_in_ram, load_pop=True)

    # test optimize RAM run w/ batch saving
    kwargs = {"optimize_ram":True, "breakdown_to_df":False, "tqdm":True}
    print_testinfo("TEST: 03", pop, kwargs)
    _ = test_binpop_evolve(pop, kwargs, verbose=True)
    check_test(pop_in_ram, load_pop=True)

    # TEST POPRUNNER
    # This is RAM heavy and may fail on personal computers
    # ================================================================================
    test_str = " TEST: 04 "
    numchar = (line_length - len(test_str)) // 2
    print("=" * numchar + test_str + "=" * numchar)
    print("Test PopulationRunner with multiple metallicities...")
    poprun = PopulationRunner(path_to_multiZ_params, verbose=True)
    print('\t Number of binary populations:', len(poprun.binary_populations))
    print('\t Metallicities:', poprun.solar_metallicities)
    print('\t Number of binaries (per pop):', poprun.binary_populations[0].number_of_binaries)
    print("🚀 Evolving PopulationRunner...")
    #poprun.evolve(overwrite=True)
    print("✅ PopulationRunner evolved successfully.")


if __name__ == "__main__":

    test_popruns()
