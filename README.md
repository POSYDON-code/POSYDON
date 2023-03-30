# POSYDON

POSYDON is a next-generation single and binary-star population synthesis code incorporating full stellar structure and evolution modeling with the use of the MESA code (https://docs.mesastar.org).

POSYDON is being developed by a collaborative team of astrophysicists and computer scientists led by Principal Investigators Tassos Fragos (Université de Genève) and Vicky Kalogera (Northwestern University). The code is modular in many aspects and the user can specify initial population properties and adopt choices that determine how stellar evolution proceeds. Populations are simulated with the use of MESA evolutionary tracks for single, non-interacting, and interacting binaries organized in grids. Machine-learning methods are incorporated and applied on the grids for classification and various interpolation calculations, and the development of irregular grids guided by active learning, for computational efficiency.

## Python distribution and virtual environments
We recommend to install the Python distribution Anaconda and to install POSYDON in a virtual environment. Specifically, we recommend using the installation we have provided which can be accessed through conda-forge. On Linux, the new conda environment can be created (we have named our environment `posydon-example`, but you can choose any name), the conda-forge channel added, and the required library installation can all be completed in one line:
```
conda create --name posydon-example -c posydon -c conda-forge posydon
```
On Mac (or if you have problems with the above command on Linux), these steps likely need to be separately run:
```
conda create -n posydon-conda python=3.7
conda activate posydon-conda
conda config --add channels conda-forge
conda config --add channels posydon
conda config --set channel_priority false
conda install posydon
```

Now, you can source the environment with
```
conda activate posydon
```

### Cloning POSYDON
Clone the repository in a local directory, e.g. '/home/POSYDON/', with
```
git clone https://github.com/POSYDON-code/POSYDON.git
git checkout origin/main
```
The directory will contain the following structure:
```
posydon/
README.md
setup.py
...
```

## Installing POSYDON


### Exporting the global path
Export the path to the cloned POSYDON code (you can add this line to your .bashrc or .bash_profile.), e.g.  
```
export PATH_TO_POSYDON=/home/POSYDON/
```

### Installing POSYDON

From the cloned POSYDON directory execute the commands to install POSYDON and the `mpi4py` dependency
```
pip install -e .
conda install mpi4py
```

### Downloading POSYDON data

#### [OPTION 1] Zenodo (latest official version)
Export the path to where you want to clone the data, e.g. `/home/`, and download the data with the following commands
```
export PATH_TO_POSYDON_DATA=/home/
get-posydon-data
```

A new folder named `POSYDON_data` will be created in the indicated location. Make sure you use a path where you have writing permissions (in some Linux distributions you may use `/home/[username]/`).

#### [OPTION 2] git LFS submodule (latest development version)
Inside the cloned POSYDON repository run the following commands to
install and initialise git LFS, and download the data.
```
export PATH_TO_POSYDON_DATA=$PATH_TO_POSYDON/data/
conda install git-lfs
git lfs install
git submodule init
git submodule update data/POSYDON_data
```

## Installing POSYDON documentations modules

These modules are needed in order to compile the documentation
```
pip install -e .[doc]
```
To compile the documentation and open the html page do the following
```
cd docs/
make html
open _build/html/index.html
```
