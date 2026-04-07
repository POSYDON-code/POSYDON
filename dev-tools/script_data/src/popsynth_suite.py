import argparse
import os
import shutil
import subprocess
import traceback
import warnings

from formatting import LINE_LENGTH, columns_to_show
from pandas.testing import assert_frame_equal
from utils import print_pop_settings, print_warnings

from posydon.config import PATH_TO_POSYDON
from posydon.popsyn.binarypopulation import BinaryPopulation
from posydon.popsyn.synthetic_population import Population, PopulationRunner

script_dir = os.path.dirname(os.path.abspath(__file__))

def test_binpop_evolve(population, popevo_kwargs, verbose=False):

    # this function runs a pop.evolve test given a set of kwargs

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
        print("=" * LINE_LENGTH)

    except Exception as e:
            print_warnings(captured_warnings)
            print(f"🚨 BinaryPopulation evolution failed!\n")
            traceback.print_exc(limit=3)
            print("\n")
            print("=" * LINE_LENGTH)

    return population

def compare_io_to_ram(loaded_pop, pop_in_ram):

    # this function compares binaries saved/loaded to any stored in RAM

    # check that binaries match between pop runs w/ fixed entropy
    # and that saved/loaded binaries match those from a memory loaded run
    N = pop_in_ram.number_of_binaries
    df_from_ram = pop_in_ram.to_df()
    ram_dflist = [df_from_ram.loc[i] for i in range(N)]
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
            print("=" * LINE_LENGTH)
            return
        if i == len(loaded_pop.history):
            break

    print("✅ Binaries from I/O match those in RAM.")
    print("=" * LINE_LENGTH)

def print_testinfo(test_title, population, popevo_kwargs):

    # prints some info about the tests

    # print test title str
    numchar = (LINE_LENGTH - len(test_title)) // 2
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

def check_test(pop_in_ram, out_path, load_pop=False):

    # checks that a test went OK

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
        save_fn = os.path.join(out_path, "batches", "evolution.combined.h5")
        loaded_pop = Population(save_fn)
        compare_io_to_ram(loaded_pop, pop_in_ram)

def test_popruns(ini_path, multiz_path, out_path, verbose):

    # primary function

    print("Performing population run tests...")
    print(f"Reading inlist: {ini_path}")
    pop = BinaryPopulation.from_ini(ini_path, verbose=verbose)
    pop.kwargs.update({"temp_directory": os.path.join(out_path, "batches")})
    print_pop_settings(pop)

    # DO TESTS:

    # test simple run, stays in RAM
    kwargs = {"optimize_ram":False, "breakdown_to_df":False, "tqdm":True}
    print_testinfo("TEST: 01", pop, kwargs)
    pop_in_ram = test_binpop_evolve(pop, kwargs, verbose=verbose)
    check_test(pop_in_ram, out_path, load_pop=False)

    # test same but w/ saving/loading binaries
    kwargs = {"optimize_ram":False, "breakdown_to_df":True, "tqdm":True}
    print_testinfo("TEST: 02", pop, kwargs)
    _ = test_binpop_evolve(pop, kwargs, verbose=verbose)
    check_test(pop_in_ram, out_path, load_pop=True)

    # test optimize RAM run w/ batch saving
    kwargs = {"optimize_ram":True, "breakdown_to_df":False, "tqdm":True}
    print_testinfo("TEST: 03", pop, kwargs)
    _ = test_binpop_evolve(pop, kwargs, verbose=verbose)
    check_test(pop_in_ram, out_path, load_pop=True)

    # TEST POPRUNNER
    # This is can be RAM heavy (may fail esp. on personal computers)
    # Using flush on print here since we are running subprocceses and want them to 
    # show in order with shell stdout.
    # ================================================================================
    os.chdir(out_path)
    test_str = " TEST: 04 "
    numchar = (LINE_LENGTH - len(test_str)) // 2
    print("=" * numchar + test_str + "=" * numchar, flush=True)
    print("Test PopulationRunner with multiple metallicities...", flush=True)
    print(f"Reading inlist: {multiz_path}", flush=True)
    poprun = PopulationRunner(multiz_path, verbose=True)
    print('\t Number of binary populations:', len(poprun.binary_populations), flush=True)
    print('\t Metallicities:', poprun.solar_metallicities, flush=True)
    print('\t Number of binaries (per pop):', poprun.binary_populations[0].number_of_binaries, flush=True)
    print("🚀 Evolving PopulationRunner...", flush=True)
    poprun.evolve(overwrite=True)
    print("✅ PopulationRunner evolved successfully.", flush=True)
    print("=" * LINE_LENGTH, flush=True)

    # TEST PIPELINE
    # This is can also be RAM heavy
    # ================================================================================
    test_str = " TEST: 05 "
    numchar = (LINE_LENGTH - len(test_str)) // 2
    print("=" * numchar + test_str + "=" * numchar, flush=True)
    print("Test posydon-popsyn pipeline for multiple metallicities...", flush=True)
    shutil.copy(os.path.join(script_dir, "setup_poprun.sh"), out_path)
    subprocess.run(["bash", "setup_poprun.sh", multiz_path], check=True)
    # mimic SLURM job array env vars, as if jobs submitted with --job_array=1
    # this is needed to test merge_metallicity.py, which looks for jobs per task ID
    # to merge.
    os.environ["SLURM_ARRAY_JOB_ID"] = "0"
    os.environ["SLURM_ARRAY_TASK_MIN"] = "0"
    os.environ["SLURM_ARRAY_TASK_ID"] = "0"
    os.environ["SLURM_ARRAY_TASK_COUNT"] = "1"

    for metallicity in poprun.solar_metallicities:
        subprocess.run(["echo", f"🚀 Running pipeline for metallicity {metallicity}..."])
        subprocess.run(["python", "run_metallicity.py", str(metallicity)], check=True)
        subprocess.run(["python", "merge_metallicity.py", str(metallicity)], check=True)

    print("✅ Successfully evolved multiple populations with posydon-popsyn.", flush=True)



if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='Evolve test binary populations for POSYDON branch validation.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--verbose', '-v', action='store_true', default=False,
                        help='Enable verbose output')
    parser.add_argument('--output', '-o', type=str, required=True,
                        help='Path to save population synthesis output')
    parser.add_argument('--ini', type=str, default=None,
                        help='Path to params ini file (auto-detected if not given)')
    parser.add_argument('--multiz', type=str, default=None,
                        help='Path to params ini file with multiple metallicities (auto-detected if not given)')
    args = parser.parse_args()

    test_popruns(ini_path=args.ini,
                 multiz_path=args.multiz,
                 out_path=args.output,
                 verbose=args.verbose)
