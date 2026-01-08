"""Minimal setup.py for POSYDON package.

Most configuration has been moved to pyproject.toml.
This file only handles versioneer and optional sphinx documentation builds.
"""

import os.path
import sys

sys.path.insert(0, os.path.dirname(__file__))

from setuptools import setup

import versioneer

# Get version from versioneer
__version__ = versioneer.get_version()
cmdclass = versioneer.get_cmdclass()

# Optional: Add sphinx documentation build command
try:
    from sphinx.setup_command import BuildDoc
    cmdclass["build_sphinx"] = BuildDoc
except ImportError:
    pass

# Minimal setup call - metadata is in pyproject.toml
setup(
    version=__version__,
    cmdclass=cmdclass,
)
