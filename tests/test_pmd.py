from pmd.core import (GAFF2, Equilibration, Lammps, Minimization, Pmd, System,
                      Torque)


def test_pmd_create(data_path, tmp_path):
    d = tmp_path / "result"
    expected_output = data_path / 'config.yaml'
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

    assert actual_output.read_text() == expected_output.read_text()
