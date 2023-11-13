"""Setup the posydon package."""

from __future__ import print_function
import glob
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
    'numpy == 1.19.1',
    'scipy == 1.5.2',
    'iminuit == 1.4.9',
    'configparser == 5.0.0',
    'astropy == 4.0.1',
    'pandas == 1.3.0',
    'scikit-learn == 0.21.3',
    'matplotlib ==  3.5.0',
    'more-itertools == 9.1.0',
    'matplotlib-label-lines == 0.3.8',
    'PyQt5 == 5.15.3',
    'h5py == 3.7.0',
    'psutil == 5.6.7',
    'tqdm == 4.48.2',
    'tables == 3.6.1',
    'python_utils == 3.5.2',
    'progressbar2 == 4.0.0',
]

tests_require = [
    "pytest == 6.2.2",
    "pytest-cov >= 2.4.0",
]

# For documenation
extras_require = {
    "doc": [
        "ipython",
        "sphinx <= 4.2.0",
        "numpydoc",
        "sphinx_rtd_theme",
        "sphinxcontrib_programoutput",
        "PSphinxTheme",
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
