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
import shutil
import tarfile
import textwrap
import time
import urllib.error
import urllib.request
from http.client import IncompleteRead

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
# number of download attempts and waiting time in seconds before the first
# retry, which gets doubled for each further retry
_DOWNLOAD_ATTEMPTS = 3
_RETRY_WAIT = 10
# number of fresh downloads to retry after a failed MD5 verification
_MD5_RETRIES = 1
# errors caused by a connection dropping during the download, which are worth
# a retry; an unreachable server raises a URLError and is not retried
_CONNECTION_ERRORS = (urllib.error.ContentTooShortError, IncompleteRead,
                      ConnectionResetError, ConnectionAbortedError,
                      TimeoutError)
# temporary HTTP errors, which are worth a retry: too many requests (429) and
# server errors (5xx); other HTTP errors (e.g. 404) are not retried
_RETRY_HTTP_CODES = [429] + list(range(500, 600))


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
                             "already installed or an archive is left over "
                             "(default: False)",
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
    """Get the paths, relative to PATH_TO_POSYDON_DATA, created by extracting
    a data set. Returns None, if they cannot be determined."""
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
        return ["SFR/IllustrisTNG.h5", "SFR/Zavala+21.txt",
                "selection_effects/pdet_grid.hdf5", "Sukhbold+16",
                "Patton+Sukhbold20", "Couch+2020"]
    return None

def _dataset_installed(dataset):
    """Check whether all expected paths of a data set exist."""
    expected_paths = _expected_paths(dataset)
    if not expected_paths:
        return False
    return all(os.path.exists(os.path.join(PATH_TO_POSYDON_DATA, path))
               for path in expected_paths)

def _md5_of_file(filepath):
    """Calculate the MD5 checksum of a file in chunks, to keep memory low."""
    md5 = hashlib.md5()
    with open(filepath, "rb") as file_to_check:
        for chunk in iter(lambda: file_to_check.read(65536), b""):
            md5.update(chunk)
    return md5.hexdigest()

def _archive_readable(filepath):
    """Check whether all members of a tar archive can be read."""
    try:
        with tarfile.open(filepath) as tar:
            tar.getmembers()
    except (tarfile.TarError, EOFError, OSError):
        return False
    return True

def _archive_verified(filepath, md5=None, verbose=False):
    """Check the integrity of a downloaded archive. Without an expected
    MD5 checksum, only check that the archive can be read completely."""
    if md5 is None:
        return _archive_readable(filepath)
    try:
        verified = (_md5_of_file(filepath) == md5)
    except OSError:
        # an unreadable archive cannot pass the MD5 check
        return False
    if verified and verbose:
        print("MD5 verified.")
    return verified

def _remote_size(data_url):
    """Get the size of a file on the server in bytes via a HEAD request.
    Returns None, if the server does not tell."""
    request = urllib.request.Request(data_url, method="HEAD")
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            size = response.headers.get("Content-Length")
        return None if size is None else int(size)
    except (OSError, ValueError):
        # the download itself reports problems with the server
        return None

def _check_disk_space(directory, required):
    """Raise an OSError, if the directory has less than the required number
    of bytes free."""
    free = shutil.disk_usage(directory).free
    if free < required:
        raise OSError(f"Not enough disk space in {directory}: "
                      f"{required/1e9:.1f} GB needed, but only "
                      f"{free/1e9:.1f} GB free.")

def _download_with_retries(data_url, partpath):
    """Download a file and retry with an increasing waiting time, if the
    connection drops or the server has a temporary error. The partial
    download gets removed after each failed attempt."""
    attempt = 1
    while True:
        try:
            urllib.request.urlretrieve(data_url, partpath, ProgressBar())
            return
        except (urllib.error.HTTPError,) + _CONNECTION_ERRORS as e:
            # Zenodo does not support resuming a download, hence restart it
            if os.path.exists(partpath):
                os.remove(partpath)
            permanent_http_error = (isinstance(e, urllib.error.HTTPError)
                                    and e.code not in _RETRY_HTTP_CODES)
            if permanent_http_error or attempt == _DOWNLOAD_ATTEMPTS:
                raise
            wait = _RETRY_WAIT * 2**(attempt-1)
            print(f"\nDownload interrupted ({e}), retrying in {wait} seconds "
                  f"(attempt {attempt+1} of {_DOWNLOAD_ATTEMPTS})...")
            time.sleep(wait)
            attempt += 1

def _archive_path(dataset):
    """Get the path to store the archive of a data set at, next to
    PATH_TO_POSYDON_DATA."""
    data_url = ZENODO_COLLECTION[dataset]['data']
    if data_url is None:
        raise ValueError(f"The dataset '{dataset}' has no publication yet.")
    directory = os.path.dirname(PATH_TO_POSYDON_DATA)
    if not os.path.isdir(directory):
        raise NotADirectoryError("PATH_TO_POSYDON_DATA does not refer to a "
                                 "valid directory.")
    return os.path.join(directory, os.path.basename(data_url))

def _expected_md5(dataset, MD5_check=True):
    """Get the MD5 checksum to verify the archive of a data set with. Returns
    None, if the MD5 check is skipped."""
    if not MD5_check:
        return None
    md5 = ZENODO_COLLECTION[dataset]['md5']
    if md5 is None:
        Pwarn("MD5 undefined, skip MD5 check.", "ReplaceValueWarning")
    return md5

def _clean_up_leftovers(filepath, md5=None, verbose=False, force=False):
    """Handle the leftovers of a previous, interrupted run.

    An incomplete download gets removed. A complete archive is kept to be
    extracted instead of being downloaded again, unless it is corrupted or
    a fresh download is forced.
    """
    partpath = filepath + ".part"
    if os.path.exists(partpath):
        print("Removing incomplete download "
              f"'{os.path.basename(partpath)}'...")
        os.remove(partpath)
    if os.path.exists(filepath) and force:
        print("Removing existing archive "
              f"'{os.path.basename(filepath)}' to download it again...")
        os.remove(filepath)
    if os.path.exists(filepath):
        if verbose:
            print("Verifying existing archive "
                  f"'{os.path.basename(filepath)}'...")
        if not _archive_verified(filepath, md5, verbose):
            os.remove(filepath)
            print("The existing archive did not pass the verification, "
                  "downloading it again.")

def _download_archive(dataset, filepath, md5=None, verbose=False):
    """Download and verify the archive of a data set.

    A corrupted download gets retried _MD5_RETRIES times. Before, it is
    checked that there is enough disk space for the archive and its
    extracted content, which is at least as large as the archive.
    """
    partpath = filepath + ".part"
    data_url = ZENODO_COLLECTION[dataset]['data']
    size = _remote_size(data_url)
    if size is not None:
        _check_disk_space(os.path.dirname(filepath), 2*size)
    for attempt in range(1+_MD5_RETRIES):
        if attempt > 0:
            print("The download did not pass the verification, "
                  "downloading it again.")
        print(f"Downloading POSYDON data '{dataset}' from Zenodo to "
              +os.path.dirname(filepath))
        _download_with_retries(data_url, partpath)
        os.replace(partpath, filepath)
        if _archive_verified(filepath, md5, verbose):
            return
        os.remove(filepath)
    raise ValueError(("MD5" if md5 else "Archive")+" verification failed!")

def _extract_archive(dataset, filepath, verbose=False):
    """Extract the archive of a data set and remove it afterwards."""
    # the extracted content is at least as large as the archive
    _check_disk_space(os.path.dirname(filepath), os.path.getsize(filepath))
    print(f"Extracting POSYDON data '{dataset}' from tar file...")
    with tarfile.open(filepath) as tar:
        for member in tqdm(tar.getmembers()):
            tar.extract(member=member, path=os.path.dirname(filepath))
    os.remove(filepath)
    if verbose:
        print('Removed downloaded tar file.')

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
            Download the data even if they seem to be already installed or
            an archive is left over.

    """
    if not isinstance(dataset, str):
        raise TypeError("'dataset' should be a string.")
    if dataset not in ZENODO_COLLECTION:
        raise KeyError(f"The dataset '{dataset}' is not defined.")
    filepath = _archive_path(dataset)

    # 1. skip installed datasets; a leftover archive indicates an
    #    interrupted extraction, hence the data set is incomplete
    if (not force
        and not os.path.exists(filepath)
        and _dataset_installed(dataset)
            ):
        print(f"POSYDON data '{dataset}' is already present, skipping.")
        return

    # 2. handle the leftovers of a previous, interrupted run
    md5 = _expected_md5(dataset, MD5_check)
    _clean_up_leftovers(filepath, md5, verbose, force)

    # 3. download the archive, unless a verified one is left over
    if not os.path.exists(filepath):
        _download_archive(dataset, filepath, md5, verbose)

    # 4. extract the archive + install the dataset
    _extract_archive(dataset, filepath, verbose)

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
            Download the data even if they seem to be already installed or
            an archive is left over.

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
