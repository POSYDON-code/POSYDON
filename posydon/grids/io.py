"""Provides the functions that handle input/output of Job data.

Example
-------
    print("Example using new script (Simone's runs)")
    path = "/home/konstantinos/pCloudDrive/"
           "Shared/POSYDON/BinaryGrids/example_posydon_grid"

    grid_reader = GridReader(path, verbose=True)

    print("READING THE PROFILE OF STAR1 IN THE SECOND BINARY OF THE GRID.")
    profile = np.genfromtxt(grid_reader.runs[1].final_profile1_path,
                            skip_header=5, names=True)
    print("Columns:")
    print(profile.dtype.names)
    print("Array:")
    print(profile)

Notes
-----
The input path can be a list of paths or a string. The paths can be given with
wildcard characters to incorporate multiple folders. E.g.,

path = ["/home/mygrid_part1", "/home/mygrid_part2"]

or

path = ["home/mygrid_part*"]

Runs with the same initial parameters (folder name without the `_grid_index_*`
part) are substituted. See `psygrid` module docstring for more information.


"""


__authors__ = [
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
    "Devina Misra <devina.misra@unige.ch>",
    "Scott Coughlin <scottcoughlin2014@u.northwestern.edu>",
    "Emmanouil Zapartas <ezapartas@gmail.com>",
    "Tassos Fragos <Anastasios.Fragkos@unige.ch>",
]


import os
import glob
import gzip

from posydon.utils.gridutils import read_MESA_data_file
from posydon.utils.posydonwarning import Pwarn


POSYDON_FORMAT_OPTIONS = {
    # subfolders in the grid parent folder that are unnecessary
    "ignored folders": ["make", "star1", "star2", "binary", "data",
                        "column_lists", "new_data", "template",
                        ".ipynb_checkpoints"],
    # which files contain useful metadata concerning the grid
    "grid metadata": ["grid_test.csv", "grid_test.csv.gz",
                      "grid.csv", "grid.csv.gz"],
    # which files contain useful metadata concerning individual grids
    "run metadata": ["inlist_grid_points", "summary.txt", "out.txt",
                     "inlist_grid_point.gz", "summary.txt.gz", "out.txt.gz"]
}


# Note that the `.gz` versions of these filenames will also be searched
BINARY_OUTPUT_FILE = "out.txt"
SINGLE_OUTPUT_FILE = "out_star1_formation_step0.txt"

# Possible extensions of EEP files (including zipped versions)
EEP_FILE_EXTENSIONS = [".data.eep", ".data.clean.eep",
                       ".data.eep.gz", ".data.clean.eep.gz"]


class RunReader:
    """Class for reading MESA output files of an individual run."""

    def __init__(self, path, fmt="posydon", binary=True,
                 verbose=False, verbose_maxlines=10):
        """Initialize object and immediately read the data.

        Parameters
        ----------
        path : str
            Path of the file or folder containing the MESA run output.
        fmt : str
            Format specifier. Linked to data-types, filename conventions, etc.
        binary : bool
            Whether the run belongs to a binary grid.
        verbose : bool
            If True, it prints reports of performed actions.
        verbose_maxlines : int
            In case of `verbose=True`, it sets the maximum number of lines of
            metadata files to be printed. If None, all lines will be printed.
            If 0, nothing will be printed.

        """
        self.fmt = fmt
        self.path = path
        self.verbose = verbose
        self.verbose_maxlines = verbose_maxlines
        self.binary = binary

        # used only if the run is saved in separate files
        self.history1_path = None
        self.history2_path = None
        self.binary_history_path = None
        self.final_profile1_path = None
        self.final_profile2_path = None
        self.initial_profile1_path = None
        self.initial_profile2_path = None
        self.final_star1_path = None
        self.final_star2_path = None
        self.out_txt_path = None
        self.metadata_files = []

        if fmt in ["posydon", "posydon_single"]:
            self.read_posydon_format()
        else:
            raise ValueError("Format {} not supported.".format(fmt))

    def read_posydon_format(self):
        """Read run structure, expecting the POSYDON collaboration format."""
        if self.verbose:
            print("Reading run in {}".format(self.path))

        def joined_exists(folder, filename, allow_gzip=True):
            """Return the joined path of `folder` and `filename`, if exists."""
            fullpath = os.path.join(folder, filename)
            fullpath_gz = fullpath + ".gz"
            if os.path.exists(fullpath):
                return fullpath
            if allow_gzip and os.path.exists(fullpath_gz):
                return fullpath_gz
            return None

        files = os.listdir(self.path)
        for file in files:
            fullpath = os.path.join(self.path, file)
            if file in ["binary_history.data", "binary_history.data.gz"]:
                self.binary_history_path = fullpath
            elif file in ["final_star1.mod", "final_star1.mod.gz"]:
                self.final_star1_path = fullpath
            elif file in ["final_star2.mod", "final_star2.mod.gz"]:
                self.final_star2_path = fullpath
            elif ((file in [BINARY_OUTPUT_FILE, BINARY_OUTPUT_FILE+".gz"]) and
                  (self.binary)):
                self.out_txt_path = fullpath
            elif ((file in [SINGLE_OUTPUT_FILE, SINGLE_OUTPUT_FILE+".gz"]) and
                  not self.binary):
                self.out_txt_path = fullpath
            elif file in POSYDON_FORMAT_OPTIONS["run metadata"]:
                self.metadata_files.append(fullpath)

        if joined_exists(self.path, 'LOGS1'):
            if os.path.isdir(os.path.join(self.path, 'LOGS1')):
                self.history1_path = joined_exists(
                    self.path, 'LOGS1/history.data', allow_gzip=True)
                self.final_profile1_path = joined_exists(
                    self.path, 'LOGS1/final_profile.data', allow_gzip=True)
                self.initial_profile1_path = joined_exists(
                    self.path, 'LOGS1/initial_profile.data', allow_gzip=True)

        if joined_exists(self.path, 'LOGS2'):
            if os.path.isdir(os.path.join(self.path, 'LOGS2')):
                self.history2_path = joined_exists(
                    self.path, 'LOGS2/history.data', allow_gzip=True)
                self.final_profile2_path = joined_exists(
                    self.path, 'LOGS2/final_profile.data', allow_gzip=True)
                self.initial_profile2_path = joined_exists(
                    self.path, 'LOGS2/initial_profile.data', allow_gzip=True)

        if ((self.history1_path is None) and (self.history2_path is None) and
            (self.binary_history_path is None) and
            (self.final_profile1_path is None) and
            (self.final_profile2_path is None) and
            (self.initial_profile1_path is None) and
            (self.initial_profile2_path is None) and
            (self.final_star1_path is None) and (self.final_star2_path is None)
            and (self.out_txt_path is None) and (len(self.metadata_files)==0)):
            Pwarn("No relevant files found in {}".format(self.path),
                  "MissingFilesWarning")

        if self.verbose:
            self.report()

    def report(self):
        """Report what data or metadata were found."""
        def isfound(var):
            return "Yes" if var is not None else "No"

        print("-" * 80)
        print("DATA FOUND")
        print("-" * 80)
        print("MESA screen output       :", isfound(self.out_txt_path))
        print("History of Star 1        :", isfound(self.history1_path))
        print("History of Star 2        :", isfound(self.history2_path))
        print("Binary history           :", isfound(self.binary_history_path))
        print("Final profile of Star 1  :", isfound(self.final_profile1_path))
        print("Final profile of Star 2  :", isfound(self.final_profile2_path))
        print("Initial profile of Star 1:",
              isfound(self.initial_profile1_path))
        print("Initial profile of Star 2:",
              isfound(self.initial_profile2_path))
        print("Final model of Star 1    :", isfound(self.final_star1_path))
        print("Final model of Star 2    :", isfound(self.final_star2_path))
        print("-" * 80)
        print("METADATA FILES: {}".format(len(self.metadata_files)))
        print("-" * 80)
        for metafile in self.metadata_files:
            print_meta_contents(metafile, max_lines=self.verbose_maxlines)


class GridReader:
    """Class for quick reading of MESA grid files."""

    def __init__(self, path, fmt="posydon", binary=True,
                 verbose=False, verbose_maxlines=10):
        """Initialize object and immediately load the grid file structure.

        Parameters
        ----------
        path : str
            Path of the file or folder containing the MESA runs.
        fmt : str
            Format specifier. Linked to data-types, filename conventions, etc.
        binary : bool
            Whether the grid(s) are binary.
        verbose : bool
            If True, it prints reports of performed actions.
        verbose_maxlines : int
            In case of `verbose=True`, it sets the maximum number of lines of
            metadata files to be printed. If None, all lines will be printed.
            If 0, nothing will be printed.

        """
        self.fmt = fmt
        self.path = path
        self.verbose = verbose
        self.verbose_maxlines = verbose_maxlines
        self.metadata_files = []
        self.runs = []
        self.binary = binary

        if fmt == "posydon":
            self.input_folders = self._get_input_folders()
            self.read_posydon_format()
        else:
            raise ValueError("Format {} not supported.".format(fmt))

    def _get_input_folders(self):
        """Get a list of all the folders with MESA runs."""
        if isinstance(self.path, (str, bytes)):
            folder_paths = [self.path]
        else:
            try:
                folder_paths = list(self.path)  # raise TypeErr if not iterable
            except TypeError as err:
                raise ValueError("MESA path must be a string or list") from err

        # if wildcard's in paths, expand them
        folders = dict()  # key=initial params, value=fullpath
        out_name = BINARY_OUTPUT_FILE if self.binary else SINGLE_OUTPUT_FILE
        for folder_path in folder_paths:
            if self.verbose:
                print("Searching for MESA runs in `{}`".format(folder_path))

            search_string = os.path.join(folder_path, "**/" + out_name)
            search_string_gz = search_string + ".gz"
            out_files_txt = glob.glob(search_string, recursive=True)
            out_files_gz = glob.glob(search_string_gz, recursive=True)
            out_files = list(out_files_txt) + list(out_files_gz)
            for out_file in out_files:
                fullpath = os.path.dirname(os.path.abspath(out_file))
                params_part = initial_values_from_dirname(fullpath)
                if params_part in folders:
                    Pwarn("Run in {} substitutes run in {}".format(fullpath,
                              folders[params_part]), "ReplaceValueWarning")
                folders[params_part] = fullpath

            # discover metadata files
            self.metadata_files = []
            for meta_path in POSYDON_FORMAT_OPTIONS["grid metadata"]:
                search_string = os.path.join(folder_path, meta_path)
                search_string_gz = search_string + ".gz"
                meta_files_txt = glob.glob(search_string)
                meta_files_gz = glob.glob(search_string_gz)
                meta_files = list(meta_files_txt) + list(meta_files_gz)
                self.metadata_files.extend(meta_files)

        if len(folders) == 0:
            raise ValueError("No folders found in {}".format(self.path))
        if self.verbose:
            print("Found {} MESA runs.".format(len(folders)))

        return folders.values()

    def read_posydon_format(self):
        """Read grid structure, expecting the POSYDON collaboration format."""
        if self.verbose:
            print("Reading grid in {}".format(self.path))

        for fullpath in self.input_folders:
            if fullpath not in POSYDON_FORMAT_OPTIONS["ignored folders"]:
                new_run = RunReader(fullpath, fmt=self.fmt,
                                    binary=self.binary, verbose=False)
                if new_run.out_txt_path is None:
                    Pwarn("Folder "+fullpath+" has no stdout file.",
                          "MissingFilesWarning")
                self.runs.append(new_run)

    def infer_history_columns(self, BH_cols, H1_cols, H2_cols):
        """Infer the columns found in the grid's history data files.

        Parameters
        ----------
        BH_cols : array-like
            Which columns to consider from `binary_history`.
        H1_cols : array-like
            Which columns to consider from `history1`.
        H2_cols : array-like
            Which columns to consider from `history2`.

        Returns
        -------
        list
            The names of the columns in history data. Columns in `history1` and
            `history2` are prefixed with `star1_` and `star2_` respectively.

        """
        BH_from, H1_from, H2_from = None, None, None

        for run in self.runs:
            if BH_from is None:
                BH_from = run.binary_history_path
            if H1_from is None:
                H1_from = run.history1_path
            if H2_from is None:
                H2_from = run.history2_path
            if not (BH_from is None or H1_from is None or H2_from is None):
                break

        BH, H1, H2 = [read_MESA_data_file(path, columns)
                      for path, columns in zip([BH_from, H1_from, H2_from],
                                               [BH_cols, H1_cols, H2_cols])]

        BH_names = [] if BH is None else list(BH.dtype.names)
        H1_names = [] if H1 is None else list(H1.dtype.names)
        H2_names = [] if H2 is None else list(H2.dtype.names)

        H1_names = ["S1_" + h1_name for h1_name in H1_names]
        H2_names = ["S2_" + h2_name for h2_name in H2_names]

        return BH_names + H1_names + H2_names


def print_meta_contents(path, max_lines=None, max_chars_per_line=80):
    """Print parts of metadata files for inspection.

    Parameters
    ----------
    path : str
        Path of file containing the metadata.
    max_lines : int or None
        Maximum number of lines to print. If None (default), print all of them.
    max_chars_per_line: int or "warp"
        If integer, it is the maximum number of character to be printed
        (e.g., useful when printing MESA output). If "warp", no truncation.

    """
    if max_lines == 0:
        return

    print("CONTENTS OF METADATA FILE: {}".format(path))
    print("-" * 80)

    if path.endswith(".gz"):
        f = gzip.open(path, "rt")
    else:
        f = open(path, "r")
    for i, line in enumerate(f):
        if max_lines is not None and i >= max_lines:
            break
        line = line.rstrip("\n")
        if max_chars_per_line != "warp":
            line = line[:max_chars_per_line]
        print(line)
    f.close()
    print("-" * 80)


def read_initial_values(mesa_dir):
    """Read grid point values given the MESA run directory."""
    path = os.path.join(mesa_dir, "inlist_grid_points")
    if not os.path.exists(path):
        return None

    initial_values = {}
    if path.endswith(".gz"):
        f = gzip.open(path, "rt")
    else:
        f = open(path, "r")
    for line in f:
        if "=" not in line:
            continue
        fields = line.strip().split("=")
        if len(fields) != 2:
            Pwarn("Multiple `=`, skipping line in {}.".format(path),
                  "InappropriateValueWarning")
            continue
        varname = fields[0].strip()
        valueparts = fields[1].split("d")
        value = float(valueparts[0].strip())
        if len(valueparts)>1:
            value = value*10**float(valueparts[1].strip())
        if varname == "m1":
            varname = "star_1_mass"
        elif varname == "m2":
            varname = "star_2_mass"
        elif varname == "initial_period_in_days":
            varname = "period_days"
        initial_values[varname] = value
    f.close()
    return initial_values


def initial_values_from_dirname(mesa_dir):
    """Use the name of the directory for inferring the main initial values."""
    dirname = str(os.path.basename(os.path.normpath(mesa_dir)))
    if "initial_mass" in dirname:                           # single-star grid
        if "v1/" in dirname: # version 1 dirnames don't contain initial_z
            variable_names = ["initial_mass"]
        else:
            variable_names = ["initial_mass", "initial_z"]
    else:                                                   # binary-star grid
        if "v1/" in dirname: # version 1 dirnames don't contain initial_z
            variable_names = ["m1", "m2", "initial_period_in_days"]
        else:
            variable_names = ["m1", "m2", "initial_period_in_days", "initial_z"]
        for variable_name in variable_names:
            assert variable_name in dirname

    values = [dirname.split(variable_name+"_")[1].split("_")[0]
              for variable_name in variable_names]
    return tuple(values)
