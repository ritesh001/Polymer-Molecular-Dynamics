import pandas as pd
import sys
import pmd

if __name__ == '__main__':
    # Get system ID from command line argument
    try:
        system_id = sys.argv[1]
    except BaseException:
        print('System ID is not provided, please pass in ID as an argument')
        exit()

    # Get system-dependent parameters from the csv file
    try:
        df = pd.read_csv('Sample_list_of_SMILES.csv')
        smiles = df.loc[df['ID'] == int(system_id), 'SMILES'].item()
    except BaseException:
        print('Having trouble getting info from the csv file')
        exit()

    system = pmd.System(smiles=smiles,
                        density=0.8,
                        force_field='gaff2',
                        natoms_total=5000,
                        natoms_per_chain=150)
    system.write_data(output_dir=system_id)

    lmp = pmd.Lammps(system)
    lmp.add_procedure(pmd.Minimization())
    lmp.add_procedure(
        pmd.Equilibration(Tfinal=600, Pfinal=1, Tmax=1000, Pmax=49346.163))
    lmp.add_procedure(
        pmd.TgMeasurement(Tinit=600,
                          Tfinal=100,
                          Tinterval=20,
                          step_duration=1000000))
    lmp.write_input(output_dir=system_id, lmp_input_fname='lmp.in')

    job = pmd.Job(jobname=system_id,
                  project='GT-rramprasad3-CODA20',
                  nodes=2,
                  ppn=24,
                  walltime='72:00:00')
    job.write_pbs(output_dir=system_id, pbs_fname='job.pbs')
