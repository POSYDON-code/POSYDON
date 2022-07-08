"""Functions handling the download of data from Zenodo."""

__authors__ = [
    "Jeff J Andrews <jeffrey.andrews@northwestern.edu>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
]

import os
import urllib.request
import hashlib
import progressbar
import tarfile
from tqdm import tqdm

# get path to data, if not provided use the working directory
PATH_TO_POSYDON_DATA = os.environ.get("PATH_TO_POSYDON_DATA",'./')
file = os.path.join(PATH_TO_POSYDON_DATA, "POSYDON_data.tar.gz")
PATH_TO_POSYDON_DATA = os.path.join(PATH_TO_POSYDON_DATA, 'POSYDON_data/')

data_url = "https://zenodo.org/record/6655751/files/POSYDON_data.tar.gz"
original_md5 = "8873544d9a568ebb85bccffbf1bdcd99"

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
                                              maxval=total_size)
            self.pbar.start()

        downloaded = block_num * block_size
        if downloaded < total_size:
            self.pbar.update(downloaded)
        else:
            self.pbar.finish()

def data_download(file=file, MD5_check=True, verbose=False):
    """Download data files from Zenodo if they do not exist.

       Parameters
       ----------
       file : string
            - Filename of the data file to be downloaded
        MD5_check : boolean
            - Use the MD5 checksum to make sure data is not corrupted
        verbose : boolean
            - verbose output

    """
    # First, make sure the path does not exist
    if os.path.exists(file):
        if verbose:
            print('POSYDON data alraedy exists at', file)
        return

    # Split the file into its directory and filename
    directory, filename = os.path.split(file)

    # Download the data
    print('Downloading POSYDON data from Zenodo '
          f'to PATH_TO_POSYDON_DATA={PATH_TO_POSYDON_DATA}')
    urllib.request.urlretrieve(data_url, file, ProgressBar())

    # Open,close, read file and calculate MD5 on its contents
    with open(file, 'rb') as file_to_check:

        # read contents of the file
        data = file_to_check.read()

    # Compare original MD5 with freshly calculated
    if MD5_check:

        # pipe contents of the file through
        md5_returned = hashlib.md5(data).hexdigest()

        if original_md5 == md5_returned:
            if verbose:
                print("MD5 verified.")
        else:
            # Delete file - we cannot rely upon that data
            os.remove(file)

            # Raise value error
            raise ValueError("MD5 verification failed!.")

    # extract each file
    print('Extracting POSYDON data from tar file...')
    with tarfile.open(file) as tar:
        for member in tqdm(iterable=tar.getmembers(), total=len(tar.getmembers())):
            tar.extract(member=member, path=directory)

    if os.path.exists(file):
        if verbose:
            print('Removed downloaded tar file.')
        os.remove(file)
        return

    return
