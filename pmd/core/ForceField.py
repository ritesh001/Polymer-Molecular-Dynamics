from abc import ABC, abstractmethod
from io import TextIOWrapper


class ForceField(ABC):

    @abstractmethod
    def write_settings(self, f: TextIOWrapper):
        pass

    @property
    @abstractmethod
    def charge_method(self):
        pass


class OPLS(ForceField):

    def __init__(self, charge_method: str = 'lbcc') -> None:
        self._charge_method = charge_method

    @property
    def charge_method(self):
        return self._charge_method

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

    def __init__(self, charge_method: str = 'am1bcc') -> None:
        self._charge_method = charge_method

    @property
    def charge_method(self):
        return self._charge_method

    def write_settings(self, f: TextIOWrapper):
        f.write(f'{"pair_style":<15} lj/cut/coul/long 12.0 12.0\n')
        f.write(f'{"pair_modify":<15} mix arithmetic\n')
        f.write(f'{"kspace_style":<15} pppm 1e-4\n')
        f.write(f'{"bond_style":<15} harmonic\n')
        f.write(f'{"angle_style":<15} harmonic\n')
        f.write(f'{"dihedral_style":<15} fourier\n')
        f.write(f'{"improper_style":<15} cvff\n')
        f.write(f'{"special_bonds":<15} amber\n')
