"""Setup the posydon package."""

from __future__ import print_function
import glob
import sys
import versioneer
import os.path

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
if 'test' in sys.argv:
    setup_requires = [
        'setuptools',
        'pytest-runner',
    ]
else:
    setup_requires = []


# These pretty common requirement are commented out. Various syntax types
# are all used in the example below for specifying specific version of the
# packages that are compatbile with your software.
install_requires = [
    'numpy == 1.24.2',
    'scipy == 1.10.1',
    'iminuit == 2.21.3',
    'configparser == 5.3.0',
    'astropy == 5.2.2',
    'pandas == 2.0.0',
    'scikit-learn == 1.2.2',
    'matplotlib ==  3.7.1',
    'matplotlib-label-lines == 0.5.2',
    'PyQt5 == 5.15.9',
    'h5py == 3.8.0',
    'psutil == 5.9.4',
    'tqdm == 4.65.0',
    'tables == 3.8.0',
    'progressbar2 == 4.2.0',
    'hurry.filesize == 0.9',
]

tests_require = [
    "pytest == 7.3.1",
    "pytest-cov >= 4.0.0",
]

# For documentation
extras_require = {
    "doc": [
        "ipython",
        "sphinx == 6.1.3",
        "numpydoc",
        "sphinx_rtd_theme",
        "sphinxcontrib_programoutput",
        "PSphinxTheme",
        "tensorflow == 2.12.0" # for profile interpolation
    ],
    "hpc": ["mpi4py == 3.0.3"],
}

# RUN SETUP

packagenames = find_packages()

# Executables go in a folder called bin
scripts = glob.glob(os.path.join("bin", "*"))

PACKAGENAME = "posydon"
DISTNAME = "posydon"
AUTHOR = "POSYDON Collaboration"
AUTHOR_EMAIL = "scottcoughlin2014@u.northwestern.edu"
LICENSE = "GPLv3+"
DESCRIPTION = "POSYDON the Next Generation of Population Synthesis"
GITHUBURL = "https://github.com/POSYDON-code/POSYDON"

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
    python_requires=">3.5, <4",
    use_2to3=False,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Intended Audience :: Science/Research",
        "Intended Audience :: End Users/Desktop",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics",
        "Operating System :: POSIX",
        "Operating System :: Unix",
        "Operating System :: MacOS",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3+)",
    ],
)
