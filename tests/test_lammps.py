import os
import shutil
import pytest

from pmd.core import GAFF2, Lammps
from pmd.core.Procedure import Equilibration, Minimization


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


def test_lammps_write(data_path, tmp_path):
    expected_output = data_path / 'lmp.in'
    d = tmp_path / "result"
    actual_output = d / "lmp.in"

    lmp = Lammps(read_data_from='data.lmps', force_field=GAFF2())
    lmp.add_procedure(Minimization())
    lmp.add_procedure(Equilibration(Teq=300, Peq=1, Tmax=800, Pmax=49346.163))
    lmp.write_lammps(d)

    assert actual_output.read_text() == expected_output.read_text()
