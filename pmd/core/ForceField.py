from io import TextIOWrapper
from typing import TypeVar

ForceField = TypeVar("ForceField", bound="ForceField")

OPLS_CHARGE_METHOD_OPTIONS = ['cm1a-lbcc', 'cm1a']
GAFF2_CHARGE_METHOD_OPTIONS = ['gasteiger', 'am1bcc']


class ForceField():

    def __init__(self, charge_method: str) -> None:
        self._charge_method = charge_method

    def __repr__(self) -> str:
        return type(self).__name__

    @property
    def charge_method(self) -> str:
        return self._charge_method

    def write_settings(self, f: TextIOWrapper) -> None:
        raise NotImplementedError


class OPLS(ForceField):
    '''Template OPLS object to contain OPLS force field settings

    Attributes:
        charge_method (str): Charge method; has to be one of `"cm1a-lbcc"` or
                             `"cm1a"`; default: `"cm1a-lbcc"`
    '''

    def __init__(self, charge_method: str = 'cm1a-lbcc') -> None:
        if charge_method not in OPLS_CHARGE_METHOD_OPTIONS:
            raise ValueError(f'Invalid OPLS charge method, valid options are '
                             f'{", ".join(OPLS_CHARGE_METHOD_OPTIONS)}')

        super().__init__(charge_method)

    def write_settings(self, f: TextIOWrapper):
        f.write(f'{"pair_style":<15} lj/cut/coul/long 9.0\n')
        f.write(f'{"pair_modify":<15} mix geometric tail yes\n')
        f.write(f'{"kspace_style":<15} pppm 1e-4\n')
        f.write(f'{"bond_style":<15} harmonic\n')
        f.write(f'{"angle_style":<15} harmonic\n')
        f.write(f'{"dihedral_style":<15} opls\n')
        f.write(f'{"improper_style":<15} cvff\n')
        f.write(f'{"special_bonds":<15} lj/coul 0.0 0.0 0.5\n')


class GAFF2(ForceField):
    '''Template GAFF2 object to contain GAFF2 force field settings

    Attributes:
        charge_method (str): Charge method; has to be one of `"gasteiger"` or
                             `"am1bcc"`; default: `"gasteiger"`
    '''

    def __init__(self, charge_method: str = 'gasteiger') -> None:
        if charge_method not in GAFF2_CHARGE_METHOD_OPTIONS:
            raise ValueError(f'Invalid OPLS charge method, valid options are '
                             f'{", ".join(GAFF2_CHARGE_METHOD_OPTIONS)}')

        super().__init__(charge_method)

    def write_settings(self, f: TextIOWrapper):
        f.write(f'{"pair_style":<15} lj/cut/coul/long 12.0 12.0\n')
        f.write(f'{"pair_modify":<15} mix arithmetic\n')
        f.write(f'{"kspace_style":<15} pppm 1e-4\n')
        f.write(f'{"bond_style":<15} harmonic\n')
        f.write(f'{"angle_style":<15} harmonic\n')
        f.write(f'{"dihedral_style":<15} fourier\n')
        f.write(f'{"improper_style":<15} cvff\n')
        f.write(f'{"special_bonds":<15} amber\n')
