import pytest

from pmd.core import (EMC, NPT, NVT, PSP, Equilibration, Lammps, Minimization,
                      MSDMeasurement, ShearDeformation, TensileDeformation,
                      TgMeasurement)


@pytest.mark.parametrize(
    "procedures, builder, lmp_input_fname",
    [
        (Minimization(), PSP('gaff2-gasteiger'), 'lmp_Minimization.in'),
        ([Minimization()], PSP('gaff2-gasteiger'), 'lmp_Minimization.in'),
        ([
            Minimization(),
            Equilibration(Teq=300, Peq=1, Tmax=800, Pmax=49346.163)
        ], PSP('opls-lbcc'), 'lmp_Equilibration.in'),
        ([
            Minimization(),
            Equilibration(Teq=300, Peq=1, Tmax=800, Pmax=49346.163),
            TgMeasurement(
                Tinit=600, Tfinal=100, Tinterval=20, step_duration=1000000)
        ], EMC('opls-aa'), 'lmp_TgMeasurement.in'),
        ([
            Minimization(),
            Equilibration(Teq=300, Peq=1, Tmax=800, Pmax=49346.163),
            MSDMeasurement(T=300,
                           group='type 1',
                           create_block_every=10000000,
                           duration=200000000,
                           dump_image=True,
                           reset_timestep_before_run=True)
        ], EMC('pcff'), 'lmp_MSDMeasurement.in'),
        ([
            Minimization(),
            Equilibration(Teq=300, Peq=1, Tmax=800, Pmax=49346.163),
            TensileDeformation(duration=10**7,
                               erate=10**-6,
                               T=300,
                               P=1,
                               reset_timestep_before_run=True)
        ], EMC('trappe'), 'lmp_TensileDeformation.in'),
        ([
            Minimization(),
            Equilibration(Teq=300, Peq=1, Tmax=800, Pmax=49346.163),
            ShearDeformation(duration=10**7,
                             erate=10**-6,
                             T=300,
                             reset_timestep_before_run=True)
        ], PSP('gaff2-gasteiger'), 'lmp_ShearDeformation.in'),
        ([
            Minimization(),
            Equilibration(Teq=300, Peq=1, Tmax=800, Pmax=49346.163),
            NPT(Tinit=300, Tfinal=300, Pinit=1, Pfinal=1, duration=10**7),
            NVT(Tinit=300, Tfinal=300, duration=10**8)
        ], PSP('gaff2-gasteiger'), 'lmp_NPT_NVT.in'),
    ],
)
def test_lammps_write(data_path, tmp_path, procedures, builder,
                      lmp_input_fname):
    d = tmp_path / "result"
    expected_output = data_path / lmp_input_fname
    actual_output = d / lmp_input_fname

    # procedures provided through the constructor
    lmp = Lammps(read_data_from='data.lmps',
                 get_functional_form_from=builder,
                 procedures=procedures,
                 lmp_input_fname=lmp_input_fname)
    lmp.write_lammps(d)

    assert actual_output.read_text() == expected_output.read_text()

    # procedures provided through the add_procedure function
    lmp = Lammps(read_data_from='data.lmps',
                 get_functional_form_from=builder,
                 lmp_input_fname=lmp_input_fname)
    lmp.add_procedure(procedures)
    lmp.write_lammps(d)

    assert actual_output.read_text() == expected_output.read_text()
