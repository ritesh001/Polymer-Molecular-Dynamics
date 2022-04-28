import pandas as pd
import os, sys


def get_system_info(csv_file, id):
    print('ID:', id)
    df = pd.read_csv(csv_file)
    row = df['ID'] == int(id)
    smiles = df.loc[row, 'SMILES'].item()
    smiles_solvent = df.loc[row, 'SMILES_solvent'].item()
    ratio = df.loc[row, 'Ratio'].item()
    return smiles, smiles_solvent, ratio


def mk_psp_input(smiles, solvent, ratio, natoms_per_chain, natoms_target,
                 psp_input_fname):
    from rdkit import Chem

    # Get the number of atoms of a repeating unit and determine the polymer chain length
    mol = Chem.MolFromSmiles(smiles)
    natoms = mol.GetNumAtoms(onlyExplicit=0) - 2
    length = round(natoms_per_chain / natoms)

    # Get the number of atoms of a solvent molecule
    mol_solvent = Chem.MolFromSmiles(solvent)
    natoms_solvent = mol_solvent.GetNumAtoms(onlyExplicit=0)

    # Calculate number of polymer chains and solvents based on target total number of atoms
    ntotal_singlechain = ratio * length * natoms_solvent + (length * natoms +
                                                            2)
    nchains = round(natoms_target / ntotal_singlechain)
    nsolvents = round(ratio * length * nchains)

    # Create PSP input file
    with open(psp_input_fname, 'w') as f:
        f.write('ID,smiles,Len,Num,NumConf,Loop,LeftCap,RightCap\n')
        f.write('Sol,{},{},{},1,False\n'.format(solvent, 1, nsolvents))
        f.write('Poly,{},{},{},1,False,[*][H],[*][H]'.format(
            smiles, length, nchains))

    print('--------Polymer Stats--------')
    print('Polymer SMILES:', smiles)
    print('Polymer length:', length)
    print('Polymer number:', nchains)
    print('')
    print('--------Solvent Stats--------')
    print('Solvent SMILES:', solvent)
    print('Solvent number:', nsolvents)
    print('')
    print('--------System Stats--------')
    print('Target Nsolvents/Nrepeatunits', ratio)
    print('Final Nsolvents/Nrepeatunits', nsolvents / (length * nchains))
    print('Total number of atoms:',
          nsolvents * natoms_solvent + (length * natoms + 2) * nchains)
    print('')


def run_pmd(jobname):
    import pmd

    lmp = pmd.Lammps(data_fname='amorphous_models/amor_gaff2.lmps',
                     force_field='gaff2')
    lmp.add_procedure(pmd.Minimization())
    lmp.add_procedure(
        pmd.Equilibration(Tfinal=300, Pfinal=1, Tmax=600, Pmax=49346.163))
    lmp.add_procedure(
        pmd.NPT(Tinit=300,
                Tfinal=300,
                Pinit=1,
                Pfinal=1,
                duration=10000000,
                reset_timestep=True))
    lmp.add_procedure(
        pmd.NVT(Tinit=300,
                Tfinal=300,
                duration=200000000,
                reset_timestep=False))
    lmp.write_input(output_dir=".", lmp_input_fname='lmp.in')

    job = pmd.Job(jobname=jobname,
                  project='GT-rramprasad3-CODA20',
                  nodes=2,
                  ppn=24,
                  walltime='48:00:00')
    job.write_pbs(output_dir=".", pbs_fname='job.pbs')


def run_psp(psp_input_fname, initial_density):
    import psp.AmorphousBuilder as ab

    print('--------Starting PSP--------')
    input_df = pd.read_csv(psp_input_fname, low_memory=False)
    amor = ab.Builder(input_df,
                      ID_col="ID",
                      SMILES_col="smiles",
                      Length='Len',
                      NumConf='NumConf',
                      LeftCap="LeftCap",
                      RightCap="RightCap",
                      Loop='Loop',
                      density=initial_density,
                      box_type='c',
                      BondInfo=False)
    amor.Build()
    amor.get_gaff2(output_fname='amor_gaff2.lmps',
                   atom_typing='antechamber',
                   am1bcc_charges=False,
                   swap_dict={'ns': 'n'})


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
        system_info_csv = 'solvent_diffusivity.csv'
        smiles, solvent, ratio = get_system_info(system_info_csv, system_id)
    except BaseException:
        print('Having trouble getting info from the csv file')
        exit()

    # Define system-independent parameters
    natoms_per_chain = 150
    natoms_target = 5000
    initial_density = 0.85
    psp_input_fname = 'input.csv'

    # Build a directory with the name of id and cd into it
    build_dir(system_id)
    previous_dir = os.getcwd()
    os.chdir(system_id)

    # Make input files
    mk_psp_input(smiles, solvent, ratio, natoms_per_chain, natoms_target,
                 psp_input_fname)
    run_pmd(system_id)
    run_psp(psp_input_fname, initial_density)

    # Change directory back to the original
    os.chdir(previous_dir)
