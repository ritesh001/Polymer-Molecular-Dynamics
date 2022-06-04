import glob
import os
import shutil
import pyemc

from pmd.util import build_dir, Pmdlogging, HiddenPrints


class Builder:

    def __init__(self, force_field: str) -> None:
        self._force_field = force_field

    def write_data(self) -> None:
        raise NotImplementedError


class EMC(Builder):

    def __init__(self, force_field: str) -> None:
        super().__init__(force_field)

    @build_dir
    def write_data(self,
                   output_dir,
                   output_prefix='system',
                   tmp_ff='opls-aa',
                   terminator='*[H]',
                   cleanup=True):

        previous_dir = os.getcwd()
        os.chdir(output_dir)

        # Write .esh file required to run EMC
        tmp_eshfile = '{}.esh'.format(output_prefix)
        RU_mw = Descriptors.ExactMolWt(MolFromSmiles(self.smiles))
        chainlength = int(self.mw / RU_mw)
        with open(tmp_eshfile, 'w') as f:
            f.write('#!/usr/bin/env emc_setup.pl\n')
            f.write('ITEM OPTIONS\n')
            f.write('replace true\n')
            f.write('field {}\n'.format(self._force_field))
            f.write('density {}\n'.format(self.density))
            f.write('ntotal {}\n'.format(self.ntotal))
            # f.write('emc_execute true\n')
            f.write('ITEM END\n')
            f.write('\n')
            f.write('ITEM GROUPS\n')
            f.write('RU {},1,RU:2\n'.format(self.smiles))
            f.write('terminator {},1,RU:1,1,RU:2\n'.format(terminator))
            f.write('ITEM END\n')
            f.write('\n')
            f.write('ITEM CLUSTERS\n')
            f.write('poly alternate 1\n')
            f.write('ITEM END\n')
            f.write('\n')
            f.write('ITEM POLYMERS\n')
            f.write('poly\n')
            f.write('1 RU,{},terminator,2\n'.format(chainlength))
            f.write('ITEM END\n')

        pyemc.setup(tmp_eshfile)
        pyemc.build('build.emc')

        # Clean up all EMC generated files except for the data file
        if cleanup:
            fnames = ['build.emc']
            fnames += glob.glob('*.esh')
            fnames += glob.glob('*.gz')
            fnames += glob.glob('*.in')
            fnames += glob.glob('*.vmd')
            fnames += glob.glob('*.params')
            for fname in fnames:
                try:
                    os.remove(fname)
                except Exception:
                    print(f'problem removing {fname} during cleanup')

        os.chdir(previous_dir)


class PSP(Builder):

    def __init__(self, force_field: str) -> None:
        super().__init__(force_field)

    def write_data(self, input_data: dict, density: float, data_fname: str,
                   output_dir: str, cleanup: bool) -> None:

        try:
            import psp.AmorphousBuilder as ab
        except ImportError:
            raise ImportError('PSP builder requires PSP to be installed to'
                              'function properly, please install PSP')

        Pmdlogging.info('Creating the system, this may take a while...')
        try:
            with HiddenPrints():
                amor = ab.Builder(pd.DataFrame(data=input_data),
                                  density=density,
                                  OutDir=output_dir)
                amor.Build()

                if isinstance(self._force_field, OPLS):
                    amor.get_opls(
                        output_fname=data_fname,
                        lbcc_charges=force_field.charge_method == 'cm1a-lbcc')
                elif isinstance(self._force_field, GAFF2):
                    amor.get_gaff2(
                        output_fname=data_fname,
                        atom_typing='antechamber',
                        am1bcc_charges=force_field.charge_method == 'am1bcc',
                        swap_dict={
                            'ns': 'n',
                            'nt': 'n',
                            'nv': 'nh'
                        })
            Pmdlogging.info(
                f'System file - {data_fname} successfully created in {output_dir}'
            )
        finally:
            if cleanup:
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
