"""Setup the posydon package."""

from __future__ import print_function
import glob
import sys
import os.path
sys.path.insert(0, os.path.dirname(__file__))

import versioneer

cmdclass = {}


# VERSIONING

__version__ = versioneer.get_version()
cmdclass.update(versioneer.get_cmdclass())


# TOGGLE WRAPPING C/C++ OR FORTRAN

WRAP_C_CPP_OR_FORTRAN = False

if WRAP_C_CPP_OR_FORTRAN:
    from distutils.command.sdist import sdist

    try:
        from numpy.distutils.core import setup, Extension
    except ImportError:
        raise ImportError("Building fortran extensions requires numpy.")

    cmdclass["sdist"] = sdist
else:
    from setuptools import setup, find_packages


# DOCUMENTATION

# import sphinx commands
try:
    from sphinx.setup_command import BuildDoc
except ImportError:
    pass
else:
    cmdclass["build_sphinx"] = BuildDoc

# read description
with open("README.md", "rb") as f:
    longdesc = "f.read().decode().strip()"


# DEPENDENCIES
setup_requires = [
    'setuptools >= 76.0.0',
]
if 'test' in sys.argv:
    setup_requires += [
        'pytest-runner',
    ]


# These pretty common requirement are commented out. Various syntax types
# are all used in the example below for specifying specific version of the
# packages that are compatbile with your software.
# TODO NOTE: before the v2.0.0 code release, we should froze the versions
# the correct way to do this is to make sure that they are available on
# conda and pip for all platforms we support (see prerequisites doc page).
install_requires = [
    'numpy >= 1.24.2, < 2.0.0',
    'scipy >= 1.10.1, <= 1.14.1',
    'iminuit >= 2.21.3, <= 2.30.1',
    'configparser >= 5.3.0, <= 7.1.0',
    'astropy >= 5.2.2, <= 6.1.6',
    'pandas >= 2.0.0, <= 2.2.3',
    'scikit-learn == 1.2.2',
    'matplotlib >=  3.9.0, <= 3.9.2',
    'matplotlib-label-lines >= 0.5.2, <= 0.7.0',
    'h5py >= 3.8.0, <= 3.12.1',
    'psutil >= 5.9.4, <= 6.1.0',
    'tqdm >= 4.65.0, <= 4.67.0',
    'tables >= 3.8.0, <= 3.10.1',
    'progressbar2 >= 4.2.0, <= 4.5.0', # for downloading data
    'hurry.filesize >= 0.9, <= 0.9',
    'python-dotenv >= 1.0.0, <= 1.0.1',
]

tests_require = [
    "pytest >= 7.3.1",
    "pytest-cov >= 4.0.0",
]

# For documentation
extras_require = {
    # to build documentation
    "doc": [
        "ipython",
        "sphinx >= 8.2.2",
        "numpydoc",
        "sphinx_rtd_theme",
        "sphinxcontrib_programoutput",
        "PSphinxTheme",
        "nbsphinx",
        "pandoc",
    ],
    # for experimental visualization features, e.g. VDH diagrams
    "vis": ["PyQt5 >= 5.15.9, <= 5.15.11"],
    # for profile macjhine learning features, e.g. profile interpolation
    "ml": ["tensorflow >= 2.13.0"],
    # for running population synthesis on HPC facilities
    "hpc": ["mpi4py >= 3.0.3"],
}

# RUN SETUP

packagenames = find_packages()

# Executables go in a folder called bin
scripts = glob.glob(os.path.join("bin", "*"))

PACKAGENAME = "posydon"
DISTNAME = "posydon"
AUTHOR = "POSYDON Collaboration"
AUTHOR_EMAIL = "posydon.team@gmail.com"
LICENSE = "GPLv3+"
DESCRIPTION = "POSYDON the Next Generation of Population Synthesis"
GITHUBURL = "https://github.com/POSYDON-code/POSYDON"

# Additional included files via include_package_data are defined in MANIFEST.in

setup(
    name=DISTNAME,
    provides=[PACKAGENAME],
    version=__version__,
    description=DESCRIPTION,
    long_description=longdesc,
    long_description_content_type="text/markdown",
    ext_modules=[wrapper] if WRAP_C_CPP_OR_FORTRAN else [],
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    license=LICENSE,
    packages=packagenames,
    include_package_data=True,
    cmdclass=cmdclass,
    url=GITHUBURL,
    scripts=scripts,
    setup_requires=setup_requires,
    install_requires=install_requires,
    tests_require=tests_require,
    extras_require=extras_require,
    python_requires=">=3.11, <3.12",
    use_2to3=False,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Intended Audience :: End Users/Desktop",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.11",
        "Operating System :: POSIX",
        "Operating System :: Unix",
        "Operating System :: MacOS",
        "Natural Language :: English",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3+)",
    ],
)
