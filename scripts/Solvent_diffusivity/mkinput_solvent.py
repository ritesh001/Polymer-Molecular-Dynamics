import sys

import pandas as pd

import pmd


def get_system_info(csv_file, id):
    print('ID:', id)
    df = pd.read_csv(csv_file)
    row = df['ID'] == int(id)
    smiles = df.loc[row, 'SMILES'].item()
    smiles_solvent = df.loc[row, 'SMILES_solvent'].item()
    ratio = df.loc[row, 'Ratio'].item()
    return smiles, smiles_solvent, ratio


if __name__ == '__main__':
    # Get system ID from command line argument
    try:
        system_id = sys.argv[1]
    except IndexError:
        print('System ID is not provided, please pass in ID as an argument')
        exit()

    # Get system-dependent parameters from the csv file
    try:
        system_info_csv = 'solvent_diffusivity.csv'
        smiles, solvent, ratio = get_system_info(system_info_csv, system_id)
    except FileNotFoundError:
        raise FileNotFoundError(
            'Having trouble getting info from the csv file')

    system = pmd.SolventSystem(smiles=smiles,
                               solvent_smiles=solvent,
                               ru_nsolvent_ratio=ratio,
                               force_field=pmd.GAFF2(charge_method='am1bcc'),
                               density=0.8,
                               natoms_total=5000,
                               natoms_per_chain=150)

    lmp = pmd.Lammps(read_data_from=system)
    lmp.add_procedure(pmd.Minimization())
    lmp.add_procedure(
        pmd.Equilibration(Teq=300, Peq=1, Tmax=600, Pmax=49346.163))
    lmp.add_procedure(
        pmd.NPT(Tinit=300,
                Tfinal=300,
                Pinit=1,
                Pfinal=1,
                duration=10000000,
                reset_timestep_before_run=True))
    lmp.add_procedure(
        pmd.MSDMeasurement(T=300,
                           group=system.solvent_group,
                           create_block_every=10000000,
                           duration=200000000,
                           dump_image=True,
                           reset_timestep_before_run=True))

    job = pmd.Torque(run_lammps=lmp,
                     jobname=system_id,
                     project='GT-rramprasad3-CODA20',
                     nodes=3,
                     ppn=24,
                     walltime='72:00:00')

    run = pmd.Pmd(system, lmp, job)
    run.create(system_id, save_config=True)
