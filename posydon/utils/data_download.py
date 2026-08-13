"""Functions for bin/get-posydon-data to handle the download from Zenodo

"""

__authors__ = [
    "Jeff J Andrews <jeffrey.andrews@northwestern.edu>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
]

import argparse
import hashlib
import math
import os
import tarfile
import textwrap
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

import progressbar
from tqdm import tqdm

from posydon.config import PATH_TO_POSYDON_DATA
from posydon.utils.datasets import COMPLETE_SETS, ZENODO_COLLECTION
from posydon.utils.posydonwarning import Pwarn


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
                    wrapper = textwrap.TextWrapper(initial_indent=indent, width=80,
                                                   subsequent_indent=indent)
                    print(wrapper.fill(ZENODO_COLLECTION[dataset]['description']))
                    print(wrapper.fill("more information at "
                                       +ZENODO_COLLECTION[dataset]['url']))

def _download_chunk_range(url, start, end, chunk_id, filepath, progress_callback):
    """Downloads a single block range from Zenodo and updates a progress tracking file."""
    headers = {"Range": f"bytes={start}-{end}"}
    req = urllib.request.Request(url, headers=headers)

    # Write directly to the assigned segment offset using standard file seeks
    with urllib.request.urlopen(req) as response:
        with open(filepath, "r+b") as f:
            f.seek(start)
            block_size = 1024 * 1024
            while True:
                buffer = response.read(block_size)
                if not buffer:
                    break
                f.write(buffer)
                # Fallback hook for your custom progress tracker context
                progress_callback(len(buffer))

def download_one_dataset(dataset='DR2_1Zsun', MD5_check=True, verbose=False, num_threads=8):
    """Download a data set from Zenodo if they do not exist.

        Parameters
        ----------
        dataset : string (default: 'DR2_1Zsun')
            Name of the data set to be in COMPLETE_SETS or ZENODO_COLLECTION.
        MD5_check : boolean (default: True)
            Use the MD5 check to make sure data is not corrupted.
        verbose : boolean (default: False)
            Enables verbose output.

    """
    if not isinstance(dataset, str):
        raise TypeError("'dataset' should be a string.")
    if dataset not in ZENODO_COLLECTION:
        raise KeyError(f"The dataset '{dataset}' is not defined.")

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

    if not os.path.isdir(os.path.dirname(filepath)):
        raise NotADirectoryError("PATH_TO_POSYDON_DATA does not refer to a valid directory.")
    if os.path.exists(filepath):
        raise FileExistsError(f"POSYDON data already exists at {filepath}.")

    print(f"Downloading POSYDON data '{dataset}' via {num_threads} threads...")

    # Fetch total dataset size via a HEAD request to Zenodo
    req = urllib.request.Request(data_url, method='HEAD')
    with urllib.request.urlopen(req) as resp:
        total_size = int(resp.headers.get('Content-Length', 0))

    if total_size == 0:
        raise ValueError("Could not retrieve file size from Zenodo.")

    # Create an empty file with the full size payload to safely write out of order
    with open(filepath, "wb") as f:
        f.truncate(total_size)

    # Setup clean block progress tracking via tqdm
    chunk_size = math.ceil(total_size / num_threads)

    with tqdm(total=total_size, unit='B', unit_scale=True, desc="Downloading") as pbar:
        def update_progress(bytes_written):
            pbar.update(bytes_written)

        # Fire parallel thread workers
        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            futures = []
            for i in range(num_threads):
                start = i * chunk_size
                end = min(total_size - 1, start + chunk_size - 1)
                futures.append(
                    executor.submit(_download_chunk_range, data_url, start, end, i, filepath, update_progress)
                )

            # Wait for all thread pools to successfully resolve
            for future in as_completed(futures):
                future.result()

    # Compare original MD5 with freshly calculated
    if MD5_check:
        print("Verifying MD5 Checksum...")
        try:
            # Use chunking, 4096 bytes at a time to not overload RAM with large files.
            md5_hash = hashlib.md5()
            with open(filepath, "rb") as f:
                for byte_block in iter(lambda: f.read(4096), b""):
                    md5_hash.update(byte_block)
            md5_returned = md5_hash.hexdigest()

            if original_md5 == md5_returned:
                if verbose:
                    print("MD5 verified.")
            else:
                os.remove(filepath)
                raise ValueError("MD5 verification failed!.")
        except Exception as e:
            print(f'Failed MD5 verification check: {e}')

    # Extract file archives
    print(f"Extracting POSYDON data '{dataset}' from tar file...")
    with tarfile.open(filepath) as tar:
        for member in tqdm(iterable=tar.getmembers(), total=len(tar.getmembers())):
            tar.extract(member=member, path=directory)

    if os.path.exists(filepath):
        if verbose:
            print('Removed downloaded tar file.')
        os.remove(filepath)

def data_download(set_name='DR2', MD5_check=True, verbose=False):
    """Download data files from Zenodo if they do not exist.

        Parameters
        ----------
        set_name : string (default: 'DR2')
            Name of the data set to be in COMPLETE_SETS or ZENODO_COLLECTION.
        MD5_check : boolean (default: True)
            Use the MD5 check to make sure data is not corrupted.
        verbose : boolean (default: False)
            Enables verbose output.

    """
    if not isinstance(set_name, str):
        raise TypeError("'set_name' should be a string.")
    # Check whether the set is in the complete sets or just a single dataset.
    if set_name in COMPLETE_SETS:
        for dataset in COMPLETE_SETS[set_name]:
            download_one_dataset(dataset=dataset, MD5_check=MD5_check,
                                 verbose=verbose)
    elif set_name in ZENODO_COLLECTION:
        if verbose:
            print("You are downloading a single data set, which might not "
                  "contain all the data needed.")
        download_one_dataset(dataset=set_name, MD5_check=MD5_check,
                             verbose=verbose)
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
                      verbose=args.verbose)
