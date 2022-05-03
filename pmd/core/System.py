from ast import Str
from typing import TypeVar

System = TypeVar("System", bound="System")


class System:
    '''Template object to contain System initialization settings

    Attributes:
        smiles (str): SMILES string of the polymer (use * as connecting point)

        density (float): Density of the system

        force_field (str): Force field (One of `gaff2` and `opls`)

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
                 force_field: str,
                 density: float,
                 natoms_total: int,
                 natoms_per_chain: int = None,
                 mw_per_chain: int = None,
                 ru_per_chain: int = None,
                 data_fname: str = 'data.lmps'):

        if (force_field != 'opls' and force_field != 'gaff2'):
            raise Exception('Force field options are opls and gaff2')

        chain_length_options = [natoms_per_chain, mw_per_chain, ru_per_chain]
        num_given_options = sum(option is not None
                                for option in chain_length_options)
        if num_given_options == 0:
            raise Exception('One of natoms_per_chain, mw_per_chain, and '
                            'ru_per_chain has to be provided')
        elif num_given_options > 1:
            raise Exception('Only one of natoms_per_chain, mw_per_chain, and '
                            'ru_per_chain can be provided')

        self._smiles = smiles
        self._density = density
        self._force_field = force_field
        self._natoms_total = natoms_total
        self._mw_per_chain = mw_per_chain
        self._natoms_per_chain = natoms_per_chain
        self._ru_per_chain = ru_per_chain
        self._data_fname = data_fname

    def get_data_fname(self) -> str:
        return self._data_fname

    def get_force_field(self) -> str:
        return self._force_field

    def update_smiles(self, smiles: str) -> System:
        '''Method to update the SMILES of this system

        Parameters:
            smiles (str): Updated SMILES string
        
        Returns:
            system (System): System instance itself (builder design pattern)
        '''
        self._smiles = smiles
        return self

    def update_natoms_total(self, natoms_total: int) -> System:
        '''Method to update total number of atoms of the system

        Parameters:
            natoms_total (int): Updated total number of atoms
        
        Returns:
            system (System): System instance itself (builder design pattern)
        '''
        self._natoms_total = natoms_total
        return self

    def update_force_field(self, force_field: Str) -> System:
        '''Method to update force field of the system

        Parameters:
            force_field (str): Updated force field
        
        Returns:
            system (System): System instance itself (builder design pattern)
        '''
        self._force_field = force_field
        return self

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
            from rdkit import Chem
        except:
            raise Exception('System\'s write_data function requires RDKit to '
                            'function properly, please install RDKit')
        try:
            import psp.AmorphousBuilder as ab
            import pandas as pd
        except:
            raise Exception('System\'s write_data function requires PSP to '
                            'function properly, please install PSP')
        import os

        mol = Chem.MolFromSmiles(self._smiles)
        natoms_per_RU = mol.GetNumAtoms(onlyExplicit=0) - 2
        if self._natoms_per_chain:
            length = round(self._natoms_per_chain / natoms_per_RU)
        elif self._mw_per_chain:
            mw_per_RU = Chem.Descriptors.ExactMolWt(mol)
            length = round(self._mw_per_chain / mw_per_RU)
        else:
            length = self._ru_per_chain
        nchains = round(self._natoms_total / (natoms_per_RU * length + 2))

        print('--------System Stats--------')
        print('SMILES:', self._smiles)
        print('Natom_per_RU:', natoms_per_RU)
        print('length:', length)
        print('Nchains:', nchains)

        psp_csv_fname = 'psp_input.csv'
        psp_csv_fpath = os.path.join(output_dir, psp_csv_fname)
        with open(psp_csv_fpath, 'w') as f:
            f.write('ID,smiles,Len,Num,NumConf,Loop,LeftCap,RightCap\n')
            f.write('Poly,{},{},{},1,False,[*][H],[*][H]'.format(
                self._smiles, length, nchains))

        input_df = pd.read_csv(psp_csv_fpath, low_memory=False)
        amor = ab.Builder(input_df, density=self._density, OutDir=output_dir)
        amor.Build()

        if self.force_field == 'opls':
            amor.get_opls(output_fname=self._data_fname)
        elif self.force_field == 'gaff2':
            amor.get_gaff2(output_fname=self._data_fname,
                           atom_typing='antechamber',
                           swap_dict={
                               'ns': 'n',
                               'nt': 'n',
                               'nv': 'nh'
                           })

        if cleanup:
            import shutil
            dnames = ['molecules', 'packmol', 'pysimm']
            for dir in dnames:
                try:
                    shutil.rmtree(os.path.join(output_dir, dir))
                except BaseException:
                    print('problem removing {} folder during cleanup'.format(
                        dir))

            fnames = ['amor_model.data', 'amor_model.vasp']
            for file in fnames:
                try:
                    os.remove(os.path.join(output_dir, file))
                except BaseException:
                    print(
                        'problem removing {} file during cleanup'.format(file))
            try:
                os.remove("output_MB.csv")
            except BaseException:
                print('problem removing output_MB.csv during cleanup')
