import os
import sys
from contextlib import contextmanager


@contextmanager
def nullify_output(suppress_stdout=True, suppress_stderr=True):
    """
    suppress stdout and stderr messages using context manager.
    https://www.codeforests.com/2020/11/05/python-suppress-stdout-and-stderr/
    """
    stdout = sys.stdout
    stderr = sys.stderr
    devnull = open(os.devnull, "w")
    try:
        if suppress_stdout:
            sys.stdout = devnull
        if suppress_stderr:
            sys.stderr = devnull
        yield
    finally:
        if suppress_stdout:
            sys.stdout = stdout
        if suppress_stderr:
            sys.stderr = stderr
