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
from posydon.config import PATH_TO_POSYDON_DATA

file = os.path.join(PATH_TO_POSYDON_DATA, "../POSYDON_data.tar.gz")
data_url = "https://zenodo.org/record/14205146/files/POSYDON_data.tar.gz"
original_md5 = "cf645a45b9b92c2ad01e759eb1950beb"

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

def data_download(file=file, MD5_check=True, verbose=False):
    """Download data files from Zenodo if they do not exist.

       Parameters
       ----------
       file : string
            - Filename of the data file to be downloaded
        MD5_check : boolean
            - Use the MD5 check to make sure data is not corrupted
        verbose : boolean
            - verbose output

    """
    # First, make sure the path does not exist
    if os.path.exists(file):
        raise FileExistsError(f'POSYDON data already exists at {file}')

    # Second, check to make sure PATH_TO_POSYDON_DATA is defined
    if "PATH_TO_POSYDON_DATA" not in os.environ:
        raise NameError('You must define the PATH_TO_POSYDON_DATA environment '
                        'variable before downloading POSYDON datasets')

    # Split the file into its directory and filename
    directory, filename = os.path.split(file)

    # Download the data
    print('Downloading POSYDON data from Zenodo '
          f'to PATH_TO_POSYDON_DATA={PATH_TO_POSYDON_DATA}')
    urllib.request.urlretrieve(data_url, file, ProgressBar())

    # Compare original MD5 with freshly calculated
    if MD5_check:
        try:
            with open(file, "rb") as file_to_check:
                # read contents of the file
                data = file_to_check.read()

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
        except:
            print('Failed to read the tar.gz file for MD5 verification, '
                  'cannot guarantee file integrity (this error seems to '
                  'happen only on macOS).')

    # extract each file
    print('Extracting POSYDON data from tar file...')
    with tarfile.open(file) as tar:
        for member in tqdm(iterable=tar.getmembers(), total=len(tar.getmembers())):
            tar.extract(member=member, path=directory)

    # remove tar files after extracted
    if os.path.exists(file):
        if verbose:
            print('Removed downloaded tar file.')
        os.remove(file)
        return

    return
