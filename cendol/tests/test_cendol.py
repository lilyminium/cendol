"""
Unit and regression test for the cendol package.
"""

# Import package, test suite, and other packages as needed
import cendol
import pytest
import sys


def test_cendol_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "cendol" in sys.modules
