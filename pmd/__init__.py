from pmd.core.Lammps import Lammps
from pmd.core.Pmd import Pmd
from pmd.core.Job import Torque, Slurm
from pmd.core.System import System, SolventSystem
from pmd.core.ForceField import OPLS, GAFF2
from pmd.core.Procedure import Minimization, Equilibration, TgMeasurement, MSDMeasurement, Deformation, NVT, NPT
from pmd.postprocessing.Analysis import calculate_Tg
from pmd.postprocessing.TrajectoryReader import read_lammpstrj, read_lammpstrj_by_type

import sys

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
