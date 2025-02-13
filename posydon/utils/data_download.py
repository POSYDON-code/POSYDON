"""Functions handling the download of data from Zenodo.

"""

__authors__ = [
    "Jeff J Andrews <jeffrey.andrews@northwestern.edu>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Matthias Kruckow <Matthias.Kruckow@unige.ch>",
]

import os
import urllib.request
import hashlib
import progressbar
import tarfile
from tqdm import tqdm
from posydon.config import PATH_TO_POSYDON_DATA
from posydon.utils.datasets import COMPLETE_SETS, ZENODO_COLLECTION

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

def download_one_dataset(dataset='v1_for_v2.0.0-pre1', MD5_check=True,
                         verbose=False):
    """Download a data set from Zenodo if they do not exist.

        Parameters
        ----------
        dataset : string (default: 'v1_for_v2.0.0-pre1')
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

    # First, generate filename and make sure the path does not exist
    data_url = ZENODO_COLLECTION[dataset]['data']
    original_md5 = ZENODO_COLLECTION[dataset]['md5']
    filename = os.path.basename(data_url)
    directory = os.path.dirname(PATH_TO_POSYDON_DATA)
    filepath = os.path.join(directory, filename)
    if not os.path.isdir(os.path.dirname(filepath)):
        raise NotADirectoryError("PATH_TO_POSYDON_DATA does not refer to a "
                                 "valid directory.")
    if os.path.exists(filepath):
        raise FileExistsError(f"POSYDON data already exists at {filename}.")

    # Download the data
    print(f"Downloading POSYDON data from Zenodo to directory={directory}")
    urllib.request.urlretrieve(data_url, filepath, ProgressBar())

    # Compare original MD5 with freshly calculated
    if MD5_check:
        try:
            with open(filepath, "rb") as file_to_check:
                # read contents of the file
                data = file_to_check.read()

            # pipe contents of the file through
            md5_returned = hashlib.md5(data).hexdigest()

            if original_md5 == md5_returned:
                if verbose:
                    print("MD5 verified.")
            else:
                # Delete file - we cannot rely upon that data
                os.remove(filepath)

                # Raise value error
                raise ValueError("MD5 verification failed!.")
        except:
            print('Failed to read the tar.gz file for MD5 verification, '
                  'cannot guarantee file integrity (this error seems to '
                  'happen only on macOS).')

    # extract each file
    print('Extracting POSYDON data from tar file...')
    with tarfile.open(filepath) as tar:
        for member in tqdm(iterable=tar.getmembers(),
                           total=len(tar.getmembers())):
            tar.extract(member=member, path=directory)

    # remove tar files after extracted
    if os.path.exists(filepath):
        if verbose:
            print('Removed downloaded tar file.')
        os.remove(filepath)
        return

    return

def data_download(set_name='v1', MD5_check=True, verbose=False):
    """Download data files from Zenodo if they do not exist.

        Parameters
        ----------
        set_name : string (default: 'v1')
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
