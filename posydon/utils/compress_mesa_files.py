"""Functions for bin/compress-mesa to handle the compression of files created
with MESA

"""

__authors__ = [
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Philipp Moura Srivastava <philipp.msrivastava@gmail.com>",
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
]

import os
import sys
import shutil
import random
import argparse
from tqdm import tqdm
from posydon.utils.posydonwarning import Pwarn

def _parse_commandline():
    """Parse the arguments given on the command line

        Returns
        -------
        Namespace
            All the passed arguments from the command line or their defaults.

    """
    parser = argparse.ArgumentParser(description="Compressing MESA files")
    parser.add_argument("-td", "--test_dir",
                        type=str,
                        help="the path to where the testing directory should "
                             "be set up instead of compressing files",
                        default=None)
    parser.add_argument("-dsr", "--dsr",
                        type=float,
                        help="downsampling rate when creating test directory",
                        default=0.01)
    parser.add_argument("-v", "--verbose",
                        help="enable outputs",
                        default=False,
                        action='store_true')
    parser.add_argument("-d", "--debug",
                        help="enable debugging outputs",
                        default=False,
                        action='store_true')
    parser.add_argument("mesa_dir",
                        type=str,
                        help="the path to the directory containing "
                             "MESA-generated data",
                        default=None)
    args = parser.parse_args()
    return args


def textsize(filesize, floatfmt=".3g", base=1024, threshold=1000):
    """Get a human-readable file size in string.

    Parameters
    ----------
    filesize : int or float
        The size of the file, or directory, ... (in bytes)
    floatfmt : str
        The format for the float before the size descriptor (e.g., 2.34K).
    base : int
        The base of the units. Typically 1024 but 1000 can be used as well.
    threshold : int or float
        The threshold for using the next unit. For example, if base is 1024
        and threshold is 1000, then 1000 bytes will be returned as 0.98K.

    Returns
    -------
    str
        The filesize in string format using b, K, M, ...

    """
    if base not in [1000, 1024]:
        raise ValueError(f"base={base} should be 1000 or 1024")
    if threshold <= 0:
        raise ValueError(f"threshold={threshold} should be larger than 0")
    if base < threshold:
        raise ValueError(f"threshold={threshold} should be smaller or equal "\
                         f"to base={base}")
    if filesize < 0:
        return "-" + textsize(-filesize, floatfmt=floatfmt, base=base,\
                              threshold=threshold)

    units = ["b", "K", "M", "G", "T", "P", "E", "Z", "Y"]
    unit_values = [base**i for i in range(len(units))]

    for unit, unit_value in zip(units, unit_values):
        if filesize < unit_value * threshold:
            quantity = filesize / unit_value
            return f"{quantity:{floatfmt}}{unit}"

    return f"{filesize:.3g} bytes"


def set_up_test(args):
    """Set up a testing directory in the requested directory. It copies data
    from the mesa_dir into the testing directory.

    Parameters (keys in `args`)
    ---------------------------
    mesa_dir : string
        The directory where the MESA tracks are stored.
    test_dir : string
        The directory where the test directory is to be set up.
    dsr : float
        Downsampling rate when creating test directory. The test directory will
        contain a random sample of the runs from mesa_dir, downsampled by the
        factor set here.

    """
    if args.test_dir is None:
        raise NameError("--test_dir needs to be specified for set_up_test")
    elif not os.path.isdir(args.test_dir):
        raise NotADirectoryError(f"Directory {args.test_dir} does not exist.")
    if args.mesa_dir is None:
        raise NameError("mesa_dir needs to be specified for set_up_test")
    elif not os.path.isdir(args.mesa_dir):
        raise NotADirectoryError(f"Directory {args.mesa_dir} does not exist.")

    for folder in os.listdir(args.mesa_dir):
        if os.path.isdir(os.path.join(args.mesa_dir, folder)):
            is_mesa_run = False
            sub_dir = os.listdir(os.path.join(args.mesa_dir, folder))
            track_dirs = []

            for _f in sub_dir:
                if "_grid_index_" in _f:
                    is_mesa_run = True
                track_dirs.append(os.path.join(args.mesa_dir, folder, _f))

            if is_mesa_run:     # checking if directory is a mesa run
                os.mkdir(os.path.join(args.test_dir, folder))
                # choosing which tracks to copy over
                inds = random.sample(list(range(len(track_dirs))),
                                     int(len(track_dirs) * args.dsr))
                for ind in inds:
                    if os.path.isdir(track_dirs[ind]):
                        shutil.copytree(track_dirs[ind], os.path.join(
                            args.test_dir,
                            folder,
                            os.path.split(track_dirs[ind])[1]))
                    else:
                        shutil.copy(track_dirs[ind], os.path.join(
                            args.test_dir,
                            folder,
                            os.path.split(track_dirs[ind])[1]))

    print(f"Created Test Directory at {args.test_dir}.")


def get_size(start_path="."):
    """Gets the size of a directory and selects MESA files for compression and
    removal.

    Parameters
    ----------
    start_path : string
        The directory root to start the file system walk.

    Returns
    -------
    total_size : int
        The size in bytes.
    remove_files : list
        List of files to remove.
    compress_files : list
        List of files to compress.
    n_runs : int
        Number of MESA run directories.
    n_remove_files : int
        Number of files to remove.
    n_compress_files : int
        Number of files to compress.

    """
    total_size = 0
    remove_files = []
    compress_files = []
    n_runs = 0
    n_remove_files = 0
    n_compress_files = 0
    for dirpath, _, filenames in os.walk(start_path):
        if "_grid_index_" in dirpath:   # checking if directory is mesa run
            new_remove_files = []
            new_compress_files = []
            if "_grid_index_" in os.path.basename(dirpath):
                n_runs += 1
        else:
            new_remove_files = None
            new_compress_files = None
        for filename in filenames:
            filepath = os.path.join(dirpath, filename)
            # skip if it is symbolic link
            if not os.path.islink(filepath):
                total_size += os.path.getsize(filepath)
            # check for files in mesa run, whether to remove or compress it
            if new_remove_files is not None:
                name, ext = os.path.splitext(filename)
                if name == "core":
                    # remove core dump files
                    new_remove_files.append(filename)
                elif ext in [".data", ".mod", ".txt"]:
                    # compress .data, .mod, .txt files
                    new_compress_files.append(filename)
        if ((new_remove_files is not None) and
            (len(new_remove_files)>0)):
            remove_files.append((dirpath, new_remove_files))
            n_remove_files += len(new_remove_files)
        if ((new_compress_files is not None) and
            (len(new_compress_files)>0)):
            compress_files.append((dirpath, new_compress_files))
            n_compress_files += len(new_compress_files)
    return (total_size, remove_files, compress_files, n_runs,
            n_remove_files, n_compress_files)


def compress_dir(args):
    """Compresses a directory containing tracks evolved with MESA.

    Parameters (keys in `args`)
    ---------------------------
    verbose : bool
        Enable/Disable additional output.
    mesa_dir : string
        The directory where the MESA tracks are stored.

    """
    if args.mesa_dir is None:
        raise NameError("mesa_dir needs to be specified for set_up_test")
    elif not os.path.isdir(args.mesa_dir):
        raise NotADirectoryError(f"Directory {args.mesa_dir} does not exist.")

    og_size, to_remove, to_compress, n_runs, n_remove_files, n_compress_files\
        = get_size(args.mesa_dir)

    if args.verbose:
        print("remove", n_remove_files, "core dump files in", len(to_remove),\
              "directories of", n_runs, "MESA runs")
    for folder, files in tqdm(to_remove):
        for remove_file in files:
            if os.path.isfile(os.path.join(folder, remove_file)):
                if args.debug:
                    print("remove:", os.path.join(folder, remove_file))
                try:
                    os.remove(os.path.join(folder, remove_file))
                except: #pragma: no cover
                    print("Could not remove:", remove_file, "in", folder)
            else: #pragma: no cover
                raise FileNotFoundError(f"{os.path.join(folder, remove_file)}"
                                        " is not a file.")

    if args.verbose:
        print("compress", n_compress_files, "files in", len(to_compress),\
              "directories of", n_runs, "MESA runs")
    for folder, files in tqdm(to_compress):
        for compress_file in files:
            if os.path.isfile(os.path.join(folder, compress_file)):
                if args.debug:
                    print("compress:", os.path.join(folder, compress_file))
                os.system(f"gzip -1 {os.path.join(folder, compress_file)}")
            else: #pragma: no cover
                raise FileNotFoundError(f"{os.path.join(folder, remove_file)}"
                                        " is not a file.")

    new_size, to_remove, to_compress, n_runs, n_remove_files, n_compress_files\
        = get_size(args.mesa_dir)
    if args.verbose:
        print("")
        print("Compressed MESA tracks")
        print(f"Original size {textsize(og_size)} | "\
              f"Compressed size {textsize(new_size)}")
    if len(to_remove)>0: #pragma: no cover
    	Pwarn("Still files to remove: {}".format(to_remove),
    	      "IncompletenessWarning")
    if len(to_compress)>0: #pragma: no cover
    	Pwarn("Still files to compress: {}".format(to_compress),
    	      "IncompletenessWarning")


def _compress_MESA():
    """Run the compression of MESA files
    
    """
    args = _parse_commandline()
    if args.test_dir is not None:
        set_up_test(args)
    else:
        compress_dir(args)
