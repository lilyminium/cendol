"""
cendol
Tools for parametrizing a residue-based polymer force field
"""

# Add imports here
from .cendol import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
