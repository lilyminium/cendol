"""
cendol
Tools for parametrizing a residue-based polymer force field
"""

# Add imports here
from .combine import join_monomers
from .fragment import fragment_capped_monomers

# Handle versioneer
from ._version import get_versions

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
