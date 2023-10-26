import os
from contextlib import contextmanager

@contextmanager
def pushd(path):
    cwd = os.getcwd()
    if not os.path.isdir(path):
        os.makedirs(path)
    os.chdir(path)
    yield
    os.chdir(cwd)

