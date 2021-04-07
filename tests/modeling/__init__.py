
from __future__ import print_function

import os
import sys

import matplotlib

from pycsamt.utils._csamtpylog import csamtpylog

from difflib import unified_diff



if os.name == "posix" and 'DISPLAY' not in os.environ:
    
    print("MATPLOTLIB: No Display found, using non-interactive svg backend", file=sys.stderr)
    matplotlib.use('svg')
    import matplotlib.pyplot as plt
    
    pyCSAMT_TEST_HAS_DISPLAY = False
else:
    #matplotlib.use('svg')
    import matplotlib.pyplot as plt
    pyCSAMT_TEST_HAS_DISPLAY  = True
    plt.ion()
    
csamtpylog.get_csamtpy_logger(__name__).info("Testing using matplotlib backend {}".\
                                          format(matplotlib.rcParams['backend']))
def diff_files(after, before, ignores=None):
    """
    compare two files using diff
    :param ignores:
    :param before:
    :param after:
    :return: the number count of different lines
    """

    with open(before) as f2p:
        before_lines = f2p.readlines()
    with open(after) as f1p:
        after_lines = f1p.readlines()

    before_lines = [line.strip() for line in before_lines]
    after_lines = [line.strip() for line in after_lines]

    if ignores:
        for ignored_term in ignores:
            before_lines = [line for line in before_lines if ignored_term not in line]
            after_lines = [line for line in before_lines if ignored_term not in line]

    msg = "Comparing {} and {}:\n".format(before, after)

    lines = [line for line in unified_diff(
        before_lines,
        after_lines,
        fromfile="baseline ({})".format(before),
        tofile="test ({})".format(after),
        n=0)]


    if lines:
        msg += "  Found differences:\n\t" + "\n\t".join(lines)
        is_identical = False
    else:
        msg += " NO differences found."
        is_identical = True

    return is_identical, msg  
 
def reset_matplotlib():
    
    interactive = matplotlib.rcParams['interactive']
    backend = matplotlib.rcParams['backend']
    matplotlib.rcdefaults()  # reset the rcparams to default
    matplotlib.rcParams['backend'] = backend
    matplotlib.rcParams['interactive'] = interactive
    logger = csamtpylog().get_csamtpy_logger(__name__)
    
    logger.info("Testing using matplotlib backend {}".format(matplotlib.rcParams['backend']))
    