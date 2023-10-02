#! /usr/bin/env python3

##############################################################################
#  Adapted from: DendroPy Phylogenetic Computing Library (2010 Jeet Sukumaran
#  and Mark T. Holder)
##############################################################################

"""
geotaxsel testing suite.
"""

import os
import re


def get_test_file_names():
    """Get list of test file names."""
    path = os.path.dirname(__file__)
    default_tests = _get_test_file_names_from_dir(path, "geotaxsel.test.")
    return default_tests


def _get_test_file_names_from_dir(dirpath, pref):
    files = os.listdir(dirpath)
    t = []
    pat = re.compile(r"^test.*\.py$")
    for f in files:
        if pat.match(f):
            rp = pref + f[:-3]  # [:-3] to strip ".py"
            t.append(rp)
    return t
