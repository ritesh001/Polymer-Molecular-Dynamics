import os
import sys

import pandas as pd


def get_smiles(csv_file, index):
    print('ID:', index)
    df = pd.read_csv(csv_file)
    return df.loc[df['ID'] == int(index), 'smiles_string'].item()


def mk_csv(smiles, natoms_per_chain, nchains, output_fname):
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    natoms = mol.GetNumAtoms(onlyExplicit=0) - 2
    length = round(natoms_per_chain / natoms)

    print('SMILES:', smiles)
    print('Natom_per_RU:', natoms)
    print('length:', length)

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
    # amor.get_gaff2()
    amor.get_gaff2(output_fname='amor_gaff2.lmps',
                   atom_typing='antechamber',
                   swap_dict={
                       'ns': 'n',
                       'nt': 'n',
                       'nv': 'nh'
                   })


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
        smiles = get_smiles('random_500.csv', system_id)
    except BaseException:
        print('Having trouble getting info from the csv file')
        exit()

    # Define system-independent parameters
    natoms_per_chain = 150
    nchains = 27
    density = 0.85
    output_fname = 'input.csv'

    # Build a directory with the name of id and cd into it
    build_dir(system_id)
    previous_dir = os.getcwd()
    os.chdir(system_id)

    # Make input files
    mk_csv(smiles, natoms_per_chain, nchains, output_fname)
    run_psp(output_fname, density)

    # Change directory back to the original
    os.chdir(previous_dir)
