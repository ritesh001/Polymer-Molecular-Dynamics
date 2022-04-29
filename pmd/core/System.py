class System:
    '''Template object to contain System initialization settings

    Attributes:
        smiles (str): SMILES string of the polymer (use * as connecting point)

        density (float): Density of the system

        force_field (str): Force field (One of `gaff2` and `opls`)

        natoms_total (int): Total number of atoms of the system

        natoms_per_chain (int): Number of atoms of the polymer, one of this 
                                attribute and `mw_per_chain` has to be
                                provided but not both (providing both will
                                result in an error); default: `None`

        mw_per_chain (int): Molecular weight of the polymer, one of this 
                            attribute and `natoms_per_chain` has to be
                            provided but not both (providing both will
                            result in an error); default: `None`

        data_fname(str): File name of the data file, which will be read in by LAMMPS 
                         [read_data](https://docs.lammps.org/read_data.html) command,
                         only enter it at instantialization if you will not use the
                         `write_data` function
    '''

    def __init__(self,
                 smiles: str,
                 density: float,
                 force_field: str,
                 natoms_total: int,
                 natoms_per_chain: int = None,
                 mw_per_chain: int = None,
                 data_fname: str = ''):

        if (force_field != 'opls' and force_field != 'gaff2'):
            raise Exception('Force field options are opls and gaff2')

        if natoms_per_chain and mw_per_chain:
            raise Exception(
                'Only one of natoms_per_chain and mw_per_chain can be provided'
            )
        elif not natoms_per_chain and not mw_per_chain:
            raise Exception(
                'One of natoms_per_chain or mw_per_chain has to be provided')

        self._smiles = smiles
        self._mw_per_chain = mw_per_chain
        self._natoms_per_chain = natoms_per_chain
        self._natoms_total = natoms_total
        self._density = density

        # These can be accessed by Lammps objects
        self.data_fname = data_fname
        self.force_field = force_field

    def write_data(self,
                   output_dir: str,
                   data_fname: str = 'amor_data.lmps',
                   cleanup: bool = True) -> None:
        '''Method to make LAMMPS data file (which contains coordinates and force 
        field parameters)

        Parameters:
        output_dir (str): Directory for the generated LAMMPS data file

        data_fname (str): File name of the data file, which will be read in by LAMMPS 
                          [read_data](https://docs.lammps.org/read_data.html) command

        cleanup (bool): Whether to clean up files other than the LAMMPS data file PSP
                        generated
        
        Returns:
            None
        '''

        self.data_fname = data_fname
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
            amor.get_opls(output_fname=data_fname)
        elif self.force_field == 'gaff2':
            amor.get_gaff2(output_fname=data_fname,
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
