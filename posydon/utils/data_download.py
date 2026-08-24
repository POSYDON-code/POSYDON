"""Functions for bin/get-posydon-data to handle the download from Zenodo

"""

__authors__ = [
    "Jeff J Andrews <jeffrey.andrews@northwestern.edu>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
]

import argparse
import hashlib
import os
import tarfile
import textwrap
import urllib.request

import progressbar
from tqdm import tqdm

from posydon.config import PATH_TO_POSYDON_DATA
from posydon.utils.common_functions import convert_metallicity_to_string
from posydon.utils.datasets import COMPLETE_SETS, ZENODO_COLLECTION
from posydon.utils.posydonwarning import Pwarn

# grid directories below PATH_TO_POSYDON_DATA, in which the *_Zsun.h5 files
# of the grid data sets are expected (see simulationproperties.py)
_GRID_DIRS = ['single_HMS', 'single_HeMS', 'HMS-HMS', 'HMS-HMS_RLO',
              'CO-HMS_RLO', 'CO-HeMS', 'CO-HeMS_RLO']
# binary grid directories, which contain pre-trained interpolators, and the
# interpolation methods, whose interpolators are distributed with the grids
_INTERP_GRID_DIRS = ['HMS-HMS', 'HMS-HMS_RLO', 'CO-HMS_RLO', 'CO-HeMS',
                     'CO-HeMS_RLO']
_INTERP_METHODS = ['1NN_1NN', 'linear3c_kNN']


def _parse_commandline():
    """Parse the arguments given on the command-line

        Returns
        -------
        Namespace
            All the passed arguments from the commoand line or their defaults.

    """
    defined_sets = list(COMPLETE_SETS.keys()) + list(ZENODO_COLLECTION.keys())
    parser = argparse.ArgumentParser(description="Downloading POSYDON data "
                                                 "from Zenodo")
    parser.add_argument('dataset',
                        help="Name of the dataset to download (default: DR2)",
                        nargs='?',
                        default='DR2')
    parser.add_argument('-l', '--listedsets',
                        help="list the datasets: 'complete' shows the full "
                             "dataset able to run POSYDON, 'individual' lists "
                             "the datasets on zenodo, which might need others "
                             "to run population synthesis (default: complete)",
                        nargs='?',
                        const='complete',
                        choices=['complete', 'individual'])
    parser.add_argument('-n', '--nomd5check',
                        help="do not confirm md5 checksum (default: False)",
                        default=False,
                        action='store_true')
    parser.add_argument('-f', '--force',
                        help="download the data even if they seem to be "
                             "already installed (default: False)",
                        default=False,
                        action='store_true')
    parser.add_argument('-v', '--verbose',
                        help="run in Verbose Mode (default: False)",
                        default=False,
                        action='store_true')
    args = parser.parse_args()
    if args.dataset not in defined_sets:
        raise parser.error("unknown dataset, use -l to show defined sets")
    return args

class ProgressBar():
    def __init__(self):
        self.pbar = None
        self.widgets = [progressbar.Bar(marker="#",left="[",right="]"),
                        progressbar.Percentage(), " | ",
                        progressbar.FileTransferSpeed(), " | ",
                        progressbar.DataSize(), " / ",
                        progressbar.DataSize(variable="max_value"), " | ",
                        progressbar.ETA()]

    def __call__(self, block_num, block_size, total_size):
        if not self.pbar:
            self.pbar=progressbar.ProgressBar(widgets=self.widgets,
                                              max_value=total_size)
            self.pbar.start()

        downloaded = block_num * block_size
        if downloaded < total_size:
            self.pbar.update(downloaded)
        else:
            self.pbar.finish()

def list_datasets(individual_sets=False, verbose=False):
    """Print a list of available datasets

        Parameters
        ----------
        individual_sets : boolean (default: False)
            Show the individual sets or only the complete sets.
        verbose : boolean (default: False)
            Enables verbose output.

    """
    if individual_sets:
        print("Defined individual sets are:")
        for dataset in ZENODO_COLLECTION:
            prefix = f"  - '{dataset}': "
            indent = " "*len(prefix)
            wrapper = textwrap.TextWrapper(initial_indent=prefix, width=80,
                                           subsequent_indent=indent)
            print(wrapper.fill(ZENODO_COLLECTION[dataset]['title']))
            if verbose:
                wrapper = textwrap.TextWrapper(initial_indent=indent, width=80,
                                               subsequent_indent=indent)
                print(wrapper.fill(ZENODO_COLLECTION[dataset]['description']))
                print(wrapper.fill("more information at "
                                   +ZENODO_COLLECTION[dataset]['url']))
    else:
        print("Defined complete sets are:")
        for set_name,complete_set in COMPLETE_SETS.items():
            print(f"  - '{set_name}' consisting of:")
            for dataset in complete_set:
                prefix = f"    - '{dataset}': "
                indent = " "*len(prefix)
                wrapper = textwrap.TextWrapper(initial_indent=prefix, width=80,
                                               subsequent_indent=indent)
                print(wrapper.fill(ZENODO_COLLECTION[dataset]['title']))
                if verbose:
                    wrapper = textwrap.TextWrapper(initial_indent=indent,
                                                   width=80,
                                                   subsequent_indent=indent)
                    print(wrapper.fill(
                        ZENODO_COLLECTION[dataset]['description']))
                    print(wrapper.fill("more information at "
                                       +ZENODO_COLLECTION[dataset]['url']))

def _expected_paths(dataset):
    """Get the paths, relative to PATH_TO_POSYDON_DATA, whose presence
    indicates that a data set is installed.

        Parameters
        ----------
        dataset : string
            Name of the data set in ZENODO_COLLECTION.

        Returns
        -------
        list of strings or None
            Relative paths created by extracting the data set. None, if
            they cannot be determined.

    """
    if dataset.startswith('DR2_grids_') and dataset.endswith('Zsun'):
        suffix = dataset[len('DR2_grids_'):-len('Zsun')]
        try:
            z_str = convert_metallicity_to_string(float(suffix))
        except ValueError:
            return None
        return [os.path.join(grid_dir, z_str + "_Zsun.h5")
                for grid_dir in _GRID_DIRS] \
               + [os.path.join(grid_dir, "interpolators", interp_method,
                               z_str + "_Zsun.pkl")
                  for grid_dir in _INTERP_GRID_DIRS
                  for interp_method in _INTERP_METHODS]
    elif dataset == 'auxiliary':
        return ["SFR/IllustrisTNG.npz", "SFR/Zavala+21.txt",
                "selection_effects/pdet_grid.hdf5", "Sukhbold+16",
                "Patton+Sukhbold20", "Couch+2020"]
    elif dataset == 'DR1_for_v2.0.0-pre1':
        z_str = convert_metallicity_to_string(1.)
        return [os.path.join(grid_dir, z_str + "_Zsun.h5")
                for grid_dir in _GRID_DIRS] \
               + ["Sukhbold+16", "Patton+Sukhbold20", "Couch+2020"]
    elif dataset == 'DR1-super_Eddington':
        # NOTE: it replaces the solar metallicity CO-* grids of
        #       'DR1_for_v2.0.0-pre1', hence these two data sets cannot be
        #       distinguished by their extracted files: use 'force' to
        #       switch between them.
        z_str = convert_metallicity_to_string(1.)
        return [os.path.join(grid_dir, z_str + "_Zsun.h5")
                for grid_dir in ['CO-HMS_RLO', 'CO-HeMS', 'CO-HeMS_RLO']]
    return None

def _dataset_installed(dataset):
    """Check whether a data set seems to be already installed.

        Parameters
        ----------
        dataset : string
            Name of the data set in ZENODO_COLLECTION.

        Returns
        -------
        boolean
            True, if all expected paths of the data set exist.

    """
    expected_paths = _expected_paths(dataset)
    if not expected_paths:
        return False
    return all(os.path.exists(os.path.join(PATH_TO_POSYDON_DATA, path))
               for path in expected_paths)

def _md5_of_file(filepath):
    """Calculate the MD5 checksum of a file without loading it into memory.

        Parameters
        ----------
        filepath : string
            Path to the file.

        Returns
        -------
        string
            The hexadecimal MD5 checksum of the file.

    """
    md5 = hashlib.md5()
    with open(filepath, "rb") as file_to_check:
        for chunk in iter(lambda: file_to_check.read(65536), b""):
            md5.update(chunk)
    return md5.hexdigest()

def download_one_dataset(dataset='DR2_1Zsun', MD5_check=True, verbose=False,
                         force=False):
    """Download a data set from Zenodo if it is not installed yet.

        Parameters
        ----------
        dataset : string (default: 'DR2_1Zsun')
            Name of the data set to be in COMPLETE_SETS or ZENODO_COLLECTION.
        MD5_check : boolean (default: True)
            Use the MD5 check to make sure data is not corrupted.
        verbose : boolean (default: False)
            Enables verbose output.
        force : boolean (default: False)
            Download the data even if they seem to be already installed.

    """
    if not isinstance(dataset, str):
        raise TypeError("'dataset' should be a string.")
    if dataset not in ZENODO_COLLECTION:
        raise KeyError(f"The dataset '{dataset}' is not defined.")

    # skip data sets, which seem to be installed already
    if not force and _dataset_installed(dataset):
        print(f"POSYDON data '{dataset}' is already present, skipping.")
        return

    # First, generate filename and make sure the path exists
    data_url = ZENODO_COLLECTION[dataset]['data']
    if data_url is None:
        raise ValueError(f"The dataset '{dataset}' has no publication yet.")
    original_md5 = ZENODO_COLLECTION[dataset]['md5']
    if original_md5 is None:
        MD5_check = False
        Pwarn("MD5 undefined, skip MD5 check.", "ReplaceValueWarning")
    filename = os.path.basename(data_url)
    directory = os.path.dirname(PATH_TO_POSYDON_DATA)
    filepath = os.path.join(directory, filename)
    partpath = filepath + ".part"
    if not os.path.isdir(os.path.dirname(filepath)):
        raise NotADirectoryError("PATH_TO_POSYDON_DATA does not refer to a "
                                 "valid directory.")

    # handle leftovers of previous interrupted downloads: an incomplete
    # download gets removed and restarted, whereas a complete archive is
    # verified below and extracted instead of being downloaded again
    use_existing_archive = os.path.exists(filepath)
    if os.path.exists(partpath):
        print(f"Removing incomplete download '{filename}.part'...")
        os.remove(partpath)

    # download the data (unless a complete archive exists already) and
    # verify its integrity; a corrupted leftover archive gets replaced by a
    # fresh download instead of aborting
    verified = not MD5_check
    while True:
        if use_existing_archive:
            if verbose:
                print(f"Verifying existing archive '{filename}'...")
        else:
            print(f"Downloading POSYDON data '{dataset}' from Zenodo to "
                  +directory)
            urllib.request.urlretrieve(data_url, partpath, ProgressBar())
            os.replace(partpath, filepath)

        # Compare original MD5 with freshly calculated
        if MD5_check:
            try:
                verified = (_md5_of_file(filepath) == original_md5)
                if verified and verbose:
                    print("MD5 verified.")
            except OSError:
                verified = True
                print('Failed to read the tar.gz file for MD5 verification, '
                      'cannot guarantee file integrity (this error seems to '
                      'happen only on macOS).')
        if verified:
            break
        os.remove(filepath)
        if use_existing_archive:
            use_existing_archive = False
            print("The existing archive did not pass the MD5 verification, "
                  "downloading it again.")
        else:
            raise ValueError("MD5 verification failed!.")

    # extract each file
    print(f"Extracting POSYDON data '{dataset}' from tar file...")
    with tarfile.open(filepath) as tar:
        for member in tqdm(iterable=tar.getmembers(),
                           total=len(tar.getmembers())):
            tar.extract(member=member, path=directory)

    # remove tar files after extracted
    if os.path.exists(filepath):
        if verbose:
            print('Removed downloaded tar file.')
        os.remove(filepath)

def data_download(set_name='DR2', MD5_check=True, verbose=False, force=False):
    """Download data files from Zenodo if they are not installed yet.

        Parameters
        ----------
        set_name : string (default: 'DR2')
            Name of the data set to be in COMPLETE_SETS or ZENODO_COLLECTION.
        MD5_check : boolean (default: True)
            Use the MD5 check to make sure data is not corrupted.
        verbose : boolean (default: False)
            Enables verbose output.
        force : boolean (default: False)
            Download the data even if they seem to be already installed.

    """
    if not isinstance(set_name, str):
        raise TypeError("'set_name' should be a string.")
    # Check whether the set is in the complete sets or just a single dataset.
    if set_name in COMPLETE_SETS:
        for dataset in COMPLETE_SETS[set_name]:
            download_one_dataset(dataset=dataset, MD5_check=MD5_check,
                                 verbose=verbose, force=force)
    elif set_name in ZENODO_COLLECTION:
        if verbose:
            print("You are downloading a single data set, which might not "
                  "contain all the data needed.")
        download_one_dataset(dataset=set_name, MD5_check=MD5_check,
                             verbose=verbose, force=force)
    else:
        raise KeyError(f"The dataset '{set_name}' is not defined.")

def _get_posydon_data():
    """Run the data download or list the datasets

    """
    args = _parse_commandline()
    if args.listedsets == 'complete':
        list_datasets(individual_sets=False, verbose=args.verbose)
    elif args.listedsets == 'individual':
        list_datasets(individual_sets=True, verbose=args.verbose)
    else:
        data_download(set_name=args.dataset, MD5_check=not args.nomd5check,
                      verbose=args.verbose, force=args.force)
