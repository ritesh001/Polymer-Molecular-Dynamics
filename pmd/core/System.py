import re
from typing import Union
from rdkit import Chem

from pmd.core.ForceField import GAFF2, OPLS, ForceField
from pmd.util import Util


class System:
    '''Template object to contain System initialization settings

    Attributes:
        smiles (str): SMILES string of the polymer (use * as connecting point)

        density (float): Density of the system

        force_field (ForceField): Force field (One of `GAFF2` and `OPLS`)

        natoms_total (int): Total number of atoms of the system

        natoms_per_chain (int): Number of atoms per polymer chain, one of this
                                attribute, `mw_per_chain`, and `ru_per_chain`
                                has to be provided but not both (providing both
                                will result in an error); default: `None`

        mw_per_chain (int): Molecular weight of the polymer, one of this 
                            attribute, `natoms_per_chain`, and `ru_per_chain`
                            has to be provided but not both (providing both 
                            will result in an error); default: `None`

        ru_per_chain (int): Number of repeating unit per polymer chain, one of
                            this attribute, `natoms_per_chain`, and 
                            `mw_per_chain` has to be provided but not both
                            (providing both will result in an error); default:
                            `None`

        data_fname (str): File name of the output data file, which will be read in by 
                          LAMMPS [read_data](https://docs.lammps.org/read_data.html) 
                          command; default: `data.lmps`
    '''

    def __init__(self,
                 smiles: str,
                 force_field: ForceField,
                 density: float,
                 natoms_total: int,
                 *,
                 natoms_per_chain: int = None,
                 mw_per_chain: int = None,
                 ru_per_chain: int = None,
                 data_fname: str = 'data.lmps'):

        chain_length_options = [natoms_per_chain, mw_per_chain, ru_per_chain]
        num_given_options = sum(option is not None
                                for option in chain_length_options)
        if num_given_options == 0:
            raise ValueError('One of natoms_per_chain, mw_per_chain, and '
                             'ru_per_chain has to be provided')
        elif num_given_options > 1:
            raise ValueError('Only one of natoms_per_chain, mw_per_chain, and '
                             'ru_per_chain can be provided')

        # convert * to [*]
        stars_no_bracket = re.findall(r'(?<!\[)\*(?!\])', smiles)
        if len(stars_no_bracket) == 2:
            smiles = smiles.replace('*', '[*]')

        self._smiles = smiles
        self._density = density
        self._force_field = force_field
        self._natoms_total = natoms_total
        self._mw_per_chain = mw_per_chain
        self._natoms_per_chain = natoms_per_chain
        self._ru_per_chain = ru_per_chain
        self._data_fname = data_fname

        self._calculate_system_spec()

    def __repr__(self) -> str:
        return type(self).__name__

    @property
    def smiles(self) -> str:
        return self._smiles

    @property
    def data_fname(self) -> str:
        return self._data_fname

    @property
    def force_field(self) -> Union[OPLS, GAFF2]:
        return self._force_field

    @smiles.setter
    def smiles(self, smiles: str):
        self._smiles = smiles
        self._calculate_system_spec()

    def _calculate_system_spec(self):
        mol = Chem.MolFromSmiles(self._smiles)
        natoms_per_RU = mol.GetNumAtoms(onlyExplicit=0) - 2
        if self._natoms_per_chain:
            self._length = round(self._natoms_per_chain / natoms_per_RU)
        elif self._mw_per_chain:
            mw_per_RU = Chem.Descriptors.ExactMolWt(mol)
            self._length = round(self._mw_per_chain / mw_per_RU)
        else:
            self._length = self._ru_per_chain
        self._nchains = round(self._natoms_total /
                              (natoms_per_RU * self._length + 2))

        print('-----------------------System Stats-----------------------')
        print('SMILES:', self._smiles)
        print('Natom_per_RU:', natoms_per_RU)
        print('length:', self._length)
        print('Nchains:', self._nchains)

    def write_data(self, output_dir: str = '.', cleanup: bool = True) -> None:
        '''Method to make LAMMPS data file (which contains coordinates and force 
        field parameters)

        Parameters:
        output_dir (str): Directory for the generated LAMMPS data file
                          ; default: `.`

        cleanup (bool): Whether to clean up files other than the LAMMPS data 
                        file PSP generated
        
        Returns:
            None
        '''

        psp_input_data = {
            'ID': ['Poly'],
            'smiles': [self._smiles],
            'Len': [self._length],
            'Num': [self._nchains],
            'NumConf': [1],
            'Loop': [False],
            'LeftCap': ['[*][H]'],
            'RightCap': ['[*][H]']
        }

        _run_psp(psp_input_data, self._density, self._force_field,
                 self._data_fname, output_dir, cleanup)


class SolventSystem(System):
    '''Template object to contain System with solvents initialization settings

    Attributes:
        smiles (str): SMILES string of the polymer (use * as connecting point)

        solvent_smiles (str): SMILES string of the solvent
        
        ru_nsolvent_ratio (float): The ratio of total number of repeating units
                                   in the system and total number of solvent
                                   molecules

        density (float): Density of the system

        force_field (Force Field): Force field (One of `GAFF2` and `OPLS`)

        natoms_total (int): Total number of atoms of the system

        natoms_per_chain (int): Number of atoms per polymer chain, one of this
                                attribute, `mw_per_chain`, and `ru_per_chain`
                                has to be provided but not both (providing both
                                will result in an error); default: `None`

        mw_per_chain (int): Molecular weight of the polymer, one of this 
                            attribute, `natoms_per_chain`, and `ru_per_chain`
                            has to be provided but not both (providing both 
                            will result in an error); default: `None`

        ru_per_chain (int): Number of repeating unit per polymer chain, one of
                            this attribute, `natoms_per_chain`, and 
                            `mw_per_chain` has to be provided but not both
                            (providing both will result in an error); default:
                            `None`

        data_fname (str): File name of the output data file, which will be read in by 
                          LAMMPS [read_data](https://docs.lammps.org/read_data.html) 
                          command; default: `data.lmps`
    '''

    def __init__(self,
                 smiles: str,
                 solvent_smiles: str,
                 ru_nsolvent_ratio: float,
                 force_field: ForceField,
                 density: float,
                 natoms_total: int,
                 *,
                 natoms_per_chain: int = None,
                 mw_per_chain: int = None,
                 ru_per_chain: int = None,
                 data_fname: str = 'data.lmps'):

        self._solvent_smiles = solvent_smiles
        self._ru_nsolvent_ratio = ru_nsolvent_ratio

        super().__init__(smiles,
                         force_field,
                         density,
                         natoms_total,
                         natoms_per_chain=natoms_per_chain,
                         mw_per_chain=mw_per_chain,
                         ru_per_chain=ru_per_chain,
                         data_fname=data_fname)

    @property
    def solvent_group(self):
        return f'molecule <= {self._nsolvents}'

    def _calculate_system_spec(self):
        # Get the number of atoms of a repeating unit and determine the polymer
        # chain length
        mol = Chem.MolFromSmiles(self._smiles)
        natoms_per_ru = mol.GetNumAtoms(onlyExplicit=0) - 2
        if self._natoms_per_chain:
            self._length = round(self._natoms_per_chain / natoms_per_ru)
        elif self._mw_per_chain:
            mw_per_ru = Chem.Descriptors.ExactMolWt(mol)
            self._length = round(self._mw_per_chain / mw_per_ru)
        else:
            self._length = self._ru_per_chain

        # Get the number of atoms of a solvent molecule
        mol_solvent = Chem.MolFromSmiles(self._solvent_smiles)
        natoms_solvent = mol_solvent.GetNumAtoms(onlyExplicit=0)

        # Calculate number of polymer chains and solvents based on target total
        # number of atoms
        natoms_total_onechain = (self._ru_nsolvent_ratio * self._length *
                                 natoms_solvent) + (
                                     self._length * natoms_per_ru + 2)
        self._nchains = round(self._natoms_total / natoms_total_onechain)
        self._nsolvents = round(self._ru_nsolvent_ratio * self._length *
                                self._nchains)

        print('--------Polymer Stats--------')
        print('Polymer SMILES:', self._smiles)
        print('Polymer length:', self._length)
        print('Polymer number:', self._nchains)
        print('')
        print('--------Solvent Stats--------')
        print('Solvent SMILES:', self._solvent_smiles)
        print('Solvent number:', self._nsolvents)
        print('')
        print('--------System Stats--------')
        print('Target Nsolvents/Nrepeatunits', self._ru_nsolvent_ratio)
        print('Final Nsolvents/Nrepeatunits',
              self._nsolvents / (self._length * self._nchains))
        print(
            'Total number of atoms:', self._nsolvents * natoms_solvent +
            (self._length * natoms_per_ru + 2) * self._nchains)
        print('')

    def write_data(self, output_dir: str = '.', cleanup: bool = True) -> None:
        '''Method to make LAMMPS data file (which contains coordinates and force 
        field parameters)

        Parameters:
        output_dir (str): Directory for the generated LAMMPS data file
                          ; default: `.`

        cleanup (bool): Whether to clean up files other than the LAMMPS data 
                        file PSP generated
        
        Returns:
            None
        '''

        try:
            import numpy as np
        except ImportError:
            raise ImportError(
                'SolventSystem\'s write_data function requires numpy to '
                'function properly, please install numpy')

        psp_input_data = {
            'ID': ['Sol', 'Poly'],
            'smiles': [self._solvent_smiles, self._smiles],
            'Len': [1, self._length],
            'Num': [self._nsolvents, self._nchains],
            'NumConf': [1, 1],
            'Loop': [False, False],
            'LeftCap': [np.nan, '[*][H]'],
            'RightCap': [np.nan, '[*][H]']
        }

        _run_psp(psp_input_data, self._density, self._force_field,
                 self._data_fname, output_dir, cleanup)


def _run_psp(input_data: dict, density: float, force_field: ForceField,
             data_fname: str, output_dir: str, cleanup: bool) -> None:
    try:
        import psp.AmorphousBuilder as ab
        import pandas as pd
    except ImportError:
        raise ImportError('System\'s write_data function requires PSP to '
                          'function properly, please install PSP')

    print('--------Creating the system, this may take a while--------')
    try:
        with Util.HiddenPrints():
            amor = ab.Builder(pd.DataFrame(data=input_data),
                              density=density,
                              OutDir=output_dir)
            amor.Build()

            if isinstance(force_field, OPLS):
                amor.get_opls(
                    output_fname=data_fname,
                    lbcc_charges=force_field.charge_method == 'cm1a-lbcc')
            elif isinstance(force_field, GAFF2):
                amor.get_gaff2(
                    output_fname=data_fname,
                    atom_typing='antechamber',
                    am1bcc_charges=force_field.charge_method == 'am1bcc',
                    swap_dict={
                        'ns': 'n',
                        'nt': 'n',
                        'nv': 'nh'
                    })
        print('--------------- System successfully created---------------')
    finally:
        if cleanup:
            import os, shutil
            force_field_dname = ['ligpargen'] if isinstance(
                force_field, OPLS) else ['pysimm']
            dnames = ['molecules', 'packmol'] + force_field_dname
            for dir in dnames:
                try:
                    shutil.rmtree(os.path.join(output_dir, dir))
                except FileNotFoundError:
                    pass

            fnames = ['amor_model.data', 'amor_model.vasp']
            for file in fnames:
                try:
                    os.remove(os.path.join(output_dir, file))
                except FileNotFoundError:
                    pass

            fnames = ['output_MB.csv', 'molecules.csv']
            for file in fnames:
                try:
                    os.remove(file)
                except FileNotFoundError:
                    pass
