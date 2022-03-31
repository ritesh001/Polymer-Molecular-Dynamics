import pandas as pd
import os

def get_smiles(csv_file, index):
    print('ID:', index)
    df = pd.read_csv(csv_file)
    return df.loc[df['ID'] == int(index), 'smiles_string'].item()


def mk_csv(smiles, natoms_per_chain, nchains, output_fname):
    from rdkit import Chem
    
    mol = Chem.MolFromSmiles(smiles)
    natoms = mol.GetNumAtoms(onlyExplicit=0)-2
    length = round(natoms_per_chain/natoms)

    print('SMILES:', smiles)
    print('Natom_per_RU:', natoms)
    print('length:', length)

    with open(output_fname, 'w') as f:
        f.write('ID,smiles,Len,Num,NumConf,Loop,LeftCap,RightCap\n')
        f.write('Poly,{},{},{},1,False,[*][H],[*][H]'.format(smiles, length, nchains))


def run_psp(output_fname, density):
    import psp.AmorphousBuilder as ab

    input_df = pd.read_csv(output_fname, low_memory=False)
    amor = ab.Builder(
        input_df,
        ID_col='ID',
        SMILES_col='smiles',
        Length='Len',
        NumConf='NumConf',
        LeftCap = 'LeftCap',
        RightCap = 'RightCap',
        Loop='Loop',
        density=density,
        box_type='c',
        BondInfo=False
    )
    amor.Build()
    # amor.get_gaff2()
    amor.get_gaff2(output_fname='amor_gaff2.lmps', atom_typing='antechamber', swap_dict={'ns': 'n', 'nt': 'n', 'nv': 'nh'})


if __name__ == '__main__':
    # get SMILES string from a csv file
    smiles = get_smiles('../random_500.csv', os.path.basename(os.getcwd()))

    natoms_per_chain = 150
    nchains = 27
    density = 0.85
    output_fname = 'input.csv'
    mk_csv(smiles, natoms_per_chain, nchains, output_fname)
    run_psp(output_fname, density)
