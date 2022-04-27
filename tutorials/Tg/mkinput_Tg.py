import pandas as pd
import os
import sys


def get_smiles(csv_file, index):
    print('ID:', index)
    df = pd.read_csv(csv_file)
    return df.loc[df['ID'] == int(index), 'SMILES'].item()


def mk_csv(smiles, natoms_per_chain, natoms_total, output_fname):
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    natoms = mol.GetNumAtoms(onlyExplicit=0) - 2
    length = round(natoms_per_chain / natoms)
    nchains = round(natoms_total / (natoms * length + 2))

    print('SMILES:', smiles)
    print('Natom_per_RU:', natoms)
    print('length:', length)
    print('Nchains:', nchains)

    with open(output_fname, 'w') as f:
        f.write('ID,smiles,Len,Num,NumConf,Loop,LeftCap,RightCap\n')
        f.write('Poly,{},{},{},1,False,[*][H],[*][H]'.format(
            smiles, length, nchains))


def run_psp(output_fname, density):
    import psp.AmorphousBuilder as ab

    input_df = pd.read_csv(output_fname, low_memory=False)
    amor = ab.Builder(input_df,
                      ID_col='ID',
                      SMILES_col='smiles',
                      Length='Len',
                      NumConf='NumConf',
                      LeftCap='LeftCap',
                      RightCap='RightCap',
                      Loop='Loop',
                      density=density,
                      box_type='c',
                      BondInfo=False)
    amor.Build()
    amor.get_gaff2(output_fname='amor_gaff2.lmps',
                   atom_typing='antechamber',
                   swap_dict={
                       'ns': 'n',
                       'nt': 'n',
                       'nv': 'nh'
                   })


def run_pmd(jobname):
    import pmd

    lmp = pmd.Lammps(data_fname='amorphous_models/amor_gaff2.lmps',
                     force_field='gaff2')
    lmp.add_procedure(pmd.Minimization(min_style='cg'))
    lmp.add_procedure(
        pmd.Equilibration(Tfinal=600, Pfinal=1, Tmax=800, Pmax=49346.163))
    lmp.add_procedure(
        pmd.TgMeasurement(Tinit=600, Tfinal=100, Tinterval=25, step=1000000))
    lmp.write_input(output_dir=".")

    job = pmd.Job(jobname=jobname,
                  project='GT-rramprasad3-CODA20',
                  nodes=2,
                  ppn=24,
                  walltime='48:00:00')
    job.write_pbs(output_dir=".")


def build_dir(output_dir):
    try:
        os.mkdir(output_dir)
    except OSError:
        pass


if __name__ == '__main__':
    # Get system ID from command line argument
    try:
        system_id = sys.argv[1]
    except BaseException:
        print('System ID is not provided, please pass in ID as an argument')
        exit()

    # Get system-dependent parameters from the csv file
    try:
        smiles = get_smiles('Sample_list_of_SMILES.csv', system_id)
    except BaseException:
        print('Having trouble getting info from the csv file')
        exit()

    # Define system-independent parameters
    natoms_per_chain = 150
    natoms_total = 5000
    density = 0.8
    output_fname = 'input.csv'

    # Build a directory with the name of id and cd into it
    build_dir(system_id)
    previous_dir = os.getcwd()
    os.chdir(system_id)

    # Make input files
    mk_csv(smiles, natoms_per_chain, natoms_total, output_fname)
    run_psp(output_fname, density)
    run_pmd(system_id)

    # Change directory back to the original
    os.chdir(previous_dir)
