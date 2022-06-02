import os
import shutil

import pytest

from pmd.core import (GAFF2, Lammps, Minimization, Equilibration,
                      TgMeasurement, MSDMeasurement, TensileDeformation,
                      ShearDeformation)


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


@pytest.mark.parametrize(
    "procedures, lmp_input_fname",
    [
        (Minimization(), 'lmp_Minimization.in'),
        ([Minimization()], 'lmp_Minimization.in'),
        ([
            Minimization(),
            Equilibration(Teq=300, Peq=1, Tmax=800, Pmax=49346.163)
        ], 'lmp_Equilibration.in'),
        ([
            Minimization(),
            Equilibration(Teq=300, Peq=1, Tmax=800, Pmax=49346.163),
            TgMeasurement(
                Tinit=600, Tfinal=100, Tinterval=20, step_duration=1000000)
        ], 'lmp_TgMeasurement.in'),
        ([
            Minimization(),
            Equilibration(Teq=300, Peq=1, Tmax=800, Pmax=49346.163),
            MSDMeasurement(T=300,
                           group='type 1',
                           create_block_every=10000000,
                           duration=200000000,
                           dump_image=True,
                           reset_timestep_before_run=True)
        ], 'lmp_MSDMeasurement.in'),
        ([
            Minimization(),
            Equilibration(Teq=300, Peq=1, Tmax=800, Pmax=49346.163),
            TensileDeformation(duration=10**7,
                               erate=10**-6,
                               T=300,
                               P=1,
                               reset_timestep_before_run=True)
        ], 'lmp_TensileDeformation.in'),
        ([
            Minimization(),
            Equilibration(Teq=300, Peq=1, Tmax=800, Pmax=49346.163),
            ShearDeformation(duration=10**7,
                             erate=10**-6,
                             T=300,
                             reset_timestep_before_run=True)
        ], 'lmp_ShearDeformation.in'),
    ],
)
def test_lammps_write(data_path, tmp_path, procedures, lmp_input_fname):
    d = tmp_path / "result"
    expected_output = data_path / lmp_input_fname
    actual_output = d / lmp_input_fname

    lmp = Lammps(read_data_from='data.lmps',
                 force_field=GAFF2(),
                 procedures=procedures,
                 lmp_input_fname=lmp_input_fname)
    lmp.write_lammps(d)

    assert actual_output.read_text() == expected_output.read_text()
