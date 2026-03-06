import os
import shutil
import traceback
import tracemalloc
import warnings

import pandas as pd
from binaries_suite import write_binary_to_screen
from pandas.testing import assert_frame_equal
from tabulate import tabulate

from posydon.config import PATH_TO_POSYDON
from posydon.popsyn.binarypopulation import BinaryPopulation
from posydon.popsyn.synthetic_population import Population, PopulationRunner

line_length = 140
path_to_default_params = os.path.join(PATH_TO_POSYDON, "dev-tools/script_data/test_population_params.ini")

def print_pop_settings(population):

    print("\nPopulation settings:")

    ignore_kwargs = ["extra_columns", "only_select_columns", "scalar_names",
                     "include_S1", "S1_kwargs", "include_S2", "S2_kwargs",
                     "population_properties", "warnings_verbose", "history_verbose",
                      "error_checking_verbose", "use_MPI", "read_samples_from_file",
                      "RANK", "size", "optimize_ram", "ram_per_cpu",
                      "dump_rate", "tqdm", "temp_directory", "breakdown_to_df"]

    for key, val in population.kwargs.items():
        if key in ignore_kwargs:
            continue
        else:
            print(f"\t {key} : {val}")

    print("\n")

def print_warnings(captured_warnings):
    # Show warnings if any were captured
    if captured_warnings:
        print(f"⚠️  {len(captured_warnings)} warning(s) raised during evolution:")
        for i, warning in enumerate(captured_warnings[:3], 1):  # Show max 3 warnings
            print(f"   {i}. {warning['category']}: {warning['message']}")
        if len(captured_warnings) > 3:
            print(f"   ... and {len(captured_warnings) - 3} more warning(s)")
        elif len(captured_warnings) <= 3:
            for i in range(4-len(captured_warnings)):
                print("")
    else:
        print(f"No warning(s) raised during evolution\n\n")


def test_binpop_evolve(popevo_kwargs, verbose=False):

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

    #pop = BinaryPopulation.from_ini(path_to_default_params, verbose=False)
    #pop.evolve(breakdown_to_df=False, optimize_ram=False, tqdm=False)

    try:
        pop = BinaryPopulation.from_ini(path_to_default_params, verbose=False)
        if verbose:
            print_pop_settings(pop)

        print("Running BinaryPopulation.evolve() with settings:")
        for key, val in popevo_kwargs.items():
            print(f"\t {key} : {val}")

        print("\nEvolving BinaryPopulation...\n")
        pop.evolve(**popevo_kwargs)

        print_warnings(captured_warnings)
        print("Done!")
        print("=" * line_length)

    except Exception as e:
            #print("=" * line_length)
            print_warnings(captured_warnings)
            print(f"🚨 Binary Population Evolution Failed!\n")
            traceback.print_exc(limit=3)
            print("\n")
            print("=" * line_length)

    return pop

def test_popruns():

    kwargs = {"optimize_ram":False, "breakdown_to_df":False, "tqdm":False}
    pop1 = test_binpop_evolve(kwargs, verbose=True)

    # test saving
    kwargs = {"optimize_ram":False, "breakdown_to_df":True, "tqdm":False}
    pop2 = test_binpop_evolve(kwargs, verbose=False)

    loaded_pop = Population("batches/evolution.combined.h5")
    #print(loaded_pop.history)

    df2 = pop1.to_df()
    dflist = [df2.loc[i] for i in range(10)]

    for i, df in enumerate(dflist):
        df1 = df#.drop(columns=["step_times"])
        df2 = loaded_pop.history[i]#.drop(columns=["step_times"])
        cols = ['time', 'step_names', 'state', 'event', 'S1_state', 'S2_state', 'S1_mass', 'S2_mass', 'orbital_period']
        df1 = df1[cols]
        df2 = df2[cols]
        try:
            assert_frame_equal(df1, df2)
            print(i)
        except AssertionError as e:
            print(e)
            print(df1)
            print(df2)
            break
        if i == len(loaded_pop.history):
            break

    print("Done!")

    #print(pop.manager.binaries)
    #write_binary_to_screen(pop.manager.binaries[0])



    # ================================================================================

    #poprun = PopulationRunner(my_ini_filename, verbose=False)
    #print('Number of binary populations:',len(poprun.binary_populations))
    #print('Metallicities:', poprun.solar_metallicities)
    #print('Number of binaries:', poprun.binary_populations[0].number_of_binaries)
    #print("Evolving PopulationRunner...")
    #poprun.evolve(overwrite=True)
    #print("Done!")


if __name__ == "__main__":

    #tracemalloc.start()

    test_popruns()

    #current, peak = tracemalloc.get_traced_memory()
    #print(f"Current: {current/1e6:.1f} MB, Peak: {peak/1e6:.1f} MB")
