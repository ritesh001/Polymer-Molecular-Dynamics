from __future__ import unicode_literals

import os
from distutils import dir_util

from pytest import fixture

from pmd.core import GAFF2, Lammps
from pmd.core.Procedure import Equilibration, Minimization


@fixture
def datadir(tmpdir, request):
    '''
    Fixture responsible for searching a folder with the same name of test
    module and, if available, moving all contents to a temporary directory so
    tests can use them freely.
    '''
    filename = request.module.__file__
    test_dir, _ = os.path.splitext(filename)

    if os.path.isdir(test_dir):
        dir_util.copy_tree(test_dir, str(tmpdir))

    return tmpdir


def test_lammps_write(datadir, tmp_path):
    expected_output = datadir.join('lmp.in')
    lmp = Lammps(read_data_from='data.lmps', force_field=GAFF2())
    lmp.add_procedure(Minimization())
    lmp.add_procedure(Equilibration(Teq=300, Peq=1, Tmax=800, Pmax=49346.163))

    # https://stackoverflow.com/questions/36070031/creating-a-temporary-directory-in-pytest
    d = tmp_path / "tmp"
    lmp.write_lammps(d)
    p = d / "lmp.in"

    assert p.read() == expected_output.read()
