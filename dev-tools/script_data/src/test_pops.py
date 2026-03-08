import os
import traceback
import warnings

from formatting import columns_to_show, line_length
from pandas.testing import assert_frame_equal
from utils import print_pop_settings, print_warnings

from posydon.config import PATH_TO_POSYDON
from posydon.popsyn.binarypopulation import BinaryPopulation
from posydon.popsyn.synthetic_population import Population, PopulationRunner

path_to_default_params = os.path.join(PATH_TO_POSYDON, "dev-tools/script_data/inlists/test_population_params.ini")
path_to_multiZ_params = os.path.join(PATH_TO_POSYDON, "dev-tools/script_data/inlists/test_multiZ_population_params.ini")

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
            print(f"🚨 Binary Population Evolution Failed!\n")
            traceback.print_exc(limit=3)
            print("\n")
            print("=" * line_length)

    return population

def test_popruns():

    print("Performing population run tests...")

    pop = BinaryPopulation.from_ini(path_to_default_params, verbose=False)
    print_pop_settings(pop)

    # test simple run, stays in RAM
    print("🚀 Evolving a population and storing binaries in RAM.")
    kwargs = {"optimize_ram":False, "breakdown_to_df":False, "tqdm":True}
    pop_in_ram = test_binpop_evolve(pop, kwargs, verbose=True)

    # test same but w/ saving/loading binaries
    print("🚀 Evolving a population and saving binaries to a hdf5 file.")
    kwargs = {"optimize_ram":False, "breakdown_to_df":True, "tqdm":True}
    _ = test_binpop_evolve(pop, kwargs, verbose=False)
    loaded_pop = Population("batches/evolution.combined.h5")

    # check that binaries match between pop runs w/ fixed entropy
    # and that saved/loaded binaries match those from a memory loaded run
    df_from_ram = pop_in_ram.to_df()
    ram_dflist = [df_from_ram.loc[i] for i in range(10)]
    print("⚔️ Checking that binaries in RAM match those retrieved from I/O...")
    for i, ram_df in enumerate(ram_dflist):
        io_df = loaded_pop.history[i]
        ram_df = ram_df[columns_to_show]
        io_df = io_df[columns_to_show]
        try:
            assert_frame_equal(ram_df, io_df)
        except AssertionError as e:
            print("🚨 A binary from I/O does not equal the same binary stored in RAM:")
            print(e)
            print("Binary in RAM:\n", ram_df)
            print("Binary from I/O:\n", io_df)
            return
        if i == len(loaded_pop.history):
            break

    print("✅ Binaries from I/O match those in RAM.")
    print("=" * line_length)

    # TEST POPRUNNER
    # This is RAM heavy and may fail on personal computers
    # ================================================================================
    print("Test PopulationRunner with multiple metallicities...")
    poprun = PopulationRunner(path_to_multiZ_params, verbose=True)
    print('Number of binary populations:', len(poprun.binary_populations))
    print('Metallicities:', poprun.solar_metallicities)
    print('Number of binaries (per pop):', poprun.binary_populations[0].number_of_binaries)
    print("🚀 Evolving PopulationRunner...")
    poprun.evolve(overwrite=True)
    print("✅ PopulationRunner evolved successfully.")


if __name__ == "__main__":

    test_popruns()
