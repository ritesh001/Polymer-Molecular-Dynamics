import pmd

if __name__ == '__main__':
    # Build a Polyethlyene system
    system = pmd.System(smiles='*CC*',
                        density=0.5,
                        force_field=pmd.GAFF2(),
                        natoms_total=2500,
                        natoms_per_chain=150)

    # Equilibrate the system
    lmp = pmd.Lammps(read_data_from=system)
    lmp.add_procedure(pmd.Minimization())
    lmp.add_procedure(
        pmd.Equilibration(Teq=300, Peq=1, Tmax=800, Pmax=49346.163))

    # Setup for the Torque scheduling system's job file
    job = pmd.Torque(run_lammps=lmp,
                     jobname='PE_equilibration',
                     project='GT-rramprasad3-CODA20',
                     nodes=1,
                     ppn=24,
                     walltime='24:00:00')

    # Create all the files in the PE_equilibration folder
    run = pmd.Pmd(system=system, lammps=lmp, job=job)
    run.create(output_dir='PE_equilibration', save_config=True)
