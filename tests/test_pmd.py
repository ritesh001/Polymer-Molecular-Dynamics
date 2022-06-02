import pytest

from pmd.core import (GAFF2, Equilibration, Lammps, Minimization, Pmd, System,
                      Torque)
from pmd.entry import load

# TODO: Add System into tests


@pytest.fixture
def test_data(data_path):
    return {
        'config_file': data_path / 'config.yaml',
        'lmp_file': data_path / 'lmp.in',
        'job_file': data_path / 'job.pbs',
    }


def test_pmd_create(tmp_path, test_data):
    d = tmp_path / "result"
    actual_output = d / 'config.yaml'

    system = System(smiles='*CC*',
                    density=0.5,
                    force_field=GAFF2(),
                    natoms_total=2500,
                    natoms_per_chain=150)

    # Equilibrate the system
    lmp = Lammps(read_data_from=system)
    lmp.add_procedure(Minimization())
    lmp.add_procedure(Equilibration(Teq=300, Peq=1, Tmax=800, Pmax=49346.163))

    # Setup for the Torque scheduling system's job file
    job = Torque(run_lammps=lmp,
                 jobname='PE_equilibration',
                 project='GT-rramprasad3-CODA20',
                 nodes=1,
                 ppn=24,
                 walltime='24:00:00')

    # Create all the files in the PE_equilibration folder
    run = Pmd(lammps=lmp, job=job)
    run.create(output_dir=d, save_config=True)

    assert actual_output.read_text() == test_data['config_file'].read_text()


def test_pmd_load(tmp_path, test_data):
    d = tmp_path / "result"
    actual_lmp_output = d / 'lmp.in'
    actual_job_output = d / 'job.pbs'

    Pmd.load_config(test_data['config_file'], d)

    assert actual_lmp_output.read_text() == test_data['lmp_file'].read_text()
    assert actual_job_output.read_text() == test_data['job_file'].read_text()


def test_pmd_load_cli(tmp_path, test_data):
    d = tmp_path / 'result'
    actual_lmp_output = d / 'lmp.in'
    actual_job_output = d / 'job.pbs'
    load.main([str(test_data['config_file']), '-o', str(d)])

    assert actual_lmp_output.read_text() == test_data['lmp_file'].read_text()
    assert actual_job_output.read_text() == test_data['job_file'].read_text()