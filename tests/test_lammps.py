import pytest

from pmd.core import (GAFF2, Equilibration, Lammps, Minimization,
                      MSDMeasurement, ShearDeformation, TensileDeformation,
                      TgMeasurement)


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


def test_lammps_add_procedure(data_path, tmp_path):
    d = tmp_path / "result"
    expected_output = data_path / 'lmp_Minimization.in'
    actual_output = d / 'lmp_Minimization.in'

    lmp = Lammps(read_data_from='data.lmps',
                 force_field=GAFF2(),
                 lmp_input_fname='lmp_Minimization.in')
    lmp.add_procedure(Minimization())
    lmp.write_lammps(d)

    assert actual_output.read_text() == expected_output.read_text()
