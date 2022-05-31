import sys

from .core import (GAFF2, NPT, NVT, OPLS, Equilibration, Lammps, Minimization,
                   MSDMeasurement, Pmd, ShearDeformation, Slurm, SolventSystem,
                   System, TensileDeformation, TgMeasurement, Torque)
from .postprocessing import (calculate_Tg, read_lammpstrj,
                             read_lammpstrj_by_type)

if sys.version_info[:2] >= (3, 8):
    from importlib.metadata import PackageNotFoundError, version
else:
    from importlib_metadata import PackageNotFoundError, version

try:
    dist_name = __name__
    __version__ = version(dist_name)
except PackageNotFoundError:
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError

__all__ = [
    Lammps, Pmd, Torque, Slurm, System, SolventSystem, OPLS, GAFF2,
    Minimization, Equilibration, TgMeasurement, MSDMeasurement,
    TensileDeformation, ShearDeformation, NVT, NPT, calculate_Tg,
    read_lammpstrj, read_lammpstrj_by_type
]
