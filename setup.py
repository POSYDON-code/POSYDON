"""Minimal setup.py for POSYDON package.

All configuration is in pyproject.toml.
This file only handles optional sphinx documentation builds.
"""

from setuptools import setup

# Optional: Add sphinx documentation build command
cmdclass = {}
try:
    from sphinx.setup_command import BuildDoc
    cmdclass["build_sphinx"] = BuildDoc
except ImportError:
    pass

# Minimal setup call - all metadata including version is in pyproject.toml
# Version is automatically determined by setuptools-scm from git tags
setup(cmdclass=cmdclass if cmdclass else None)
