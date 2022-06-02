"""
    conftest.py for pmd.

    Read more about conftest.py under:
    - https://docs.pytest.org/en/stable/fixture.html
    - https://docs.pytest.org/en/stable/writing_plugins.html
"""

import os
import shutil

import pytest


@pytest.fixture
def data_path(tmp_path, request):
    '''
    Fixture responsible for searching a folder with the same name of test
    module and, if available, moving all contents to a temporary directory so
    tests can use them freely.
    '''
    filename = request.module.__file__
    test_dir, _ = os.path.splitext(filename)
    d = tmp_path / "data"

    if os.path.isdir(test_dir):
        shutil.copytree(test_dir, d)

    return d
