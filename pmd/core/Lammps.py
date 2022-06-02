import os
from typing import List, Optional, TypeVar, Union

# These have to be written explicitly for typing
from pmd.core.ForceField import ForceField
from pmd.core.Procedure import Procedure
from pmd.core.System import System
from pmd.util import Pmdlogging, build_dir, validate_options

Lammps = TypeVar("Lammps", bound="Lammps")

DATA_SOURCE_OPTIONS = ('read_data_from', 'read_restart_from')


class Lammps:
    '''Template object to contain LAMMPS initialization settings

    Attributes:
        read_data_from (System | str): System object that the data file will
                                 be read from. This can also be your data
                                 file name string if you do not generate your
                                 system via PMD; one of this attribute and
                                 `read_restart_from` has to be provided but not
                                 both (providing both will result in an error)
                                 ; default: `None`

        read_restart_from (Lammps | str): Lammps object that the last restart
                                    file created will be read from. This can
                                    also be your previous lammps file name
                                    string if you do not have previous Lammps
                                    object; one of this attribute and
                                    `read_data_from` has to be provided but
                                    not both (providing both will result in an
                                    error); default: `None`

        force_foeld (ForceField): Only needed if `read_data_from` or
                                  `read_restart_from` is provided as a file
                                  name string. This is needed for specifying
                                  potential form for LAMMPS input file
                                  ; default: `None`

        atom_style (str): LAMMPS
                          [atom_style](https://docs.lammps.org/atom_style.html)
                          to use during simulation; default: `"full"`

        units (str): LAMMPS [units](https://docs.lammps.org/units.html) to use
                     during simulation; default: `"real"`

        timestep (float): LAMMPS
                          [timestep](https://docs.lammps.org/timestep.html) to
                          use during simulation; default: `1.0` (in unit of fs
                          if `units` is `"real"`)

        neighbor_skin (float): LAMMPS
                            [neighbor](https://docs.lammps.org/neighbor.html)
                            skin size to use during the simulation; default:
                            `2.0 Angstrom`

        neighbor_every (int): LAMMPS
                            [neighbor](https://docs.lammps.org/neighbor.html)
                            list checking frequency to use during the
                            simulation; default: `1`

        thermo (int): LAMMPS [thermo](https://docs.lammps.org/thermo.html)
                      to use during simulation; default: `1000 timestep`

        lmp_input_fname (str): Name of the LAMMPS input file; default:
                               `"lmp.in"`
    '''

    def __init__(self,
                 read_data_from: Optional[Union[System, str]] = None,
                 read_restart_from: Optional[Union[Lammps, str]] = None,
                 force_field: Optional[ForceField] = None,
                 procedures: Optional[Union[Procedure,
                                            List[Procedure]]] = None,
                 atom_style: str = 'full',
                 units: str = 'real',
                 timestep: int = 1,
                 neighbor_skin: float = 2.0,
                 neighbor_every: int = 1,
                 thermo: int = 1000,
                 lmp_input_fname: str = 'lmp.in'):

        self._read_data_from = read_data_from
        self._read_restart_from = read_restart_from
        self._force_field = force_field
        self._atom_style = atom_style
        self._units = units
        self._timestep = timestep
        self._neighbor_skin = neighbor_skin
        self._neighbor_every = neighbor_every
        self._thermo = thermo
        self._lmp_input_fname = lmp_input_fname
        if not procedures:
            procedures = []
        elif isinstance(procedures, Procedure):
            procedures = [procedures]
        self._procedures = procedures

        # reassign data source and force field if objects are provided
        if isinstance(read_data_from, System):
            self._read_data_from = read_data_from.data_fname
            self._force_field = read_data_from.force_field
        elif isinstance(read_restart_from, Lammps):
            # TODO: implement this
            # self._read_restart_from =
            # read_restart_from.last_restart_fname
            self._force_field = read_restart_from.force_field

        # Make sure only 1 data source option is given
        validate_options(self, DATA_SOURCE_OPTIONS)

    def __repr__(self) -> str:
        return type(self).__name__

    @property
    def lmp_input_fname(self) -> str:
        return self._lmp_input_fname

    @property
    def force_field(self) -> Optional[ForceField]:
        return self._force_field

    def add_procedure(self, procedure: Union[Procedure,
                                             List[Procedure]]) -> Lammps:
        '''Method to add simulation procedure
        Parameters:
            procedure (Procedure): One of `Minimization`, `Equilibration`,
            `NPT`, `NVT`, `MSDMeasurement`, `TgMeasurement`, and `Deformation`

        Returns:
            Lammps (Lammps): Lammps instance itself (builder design pattern)
        '''
        if not isinstance(procedure, list):
            procedure = [procedure]

        for p in procedure:
            self._procedures.append(p)
        return self

    @build_dir
    def write_lammps(self, output_dir: str = '.') -> None:
        '''Method to make LAMMPS input files
        Parameters:
            output_dir (str): Directory for the generated LAMMPS input file
                              ; default: `"."`

        Returns:
            None
        '''

        Pmdlogging.info('Creating Lammps input file...')

        # Write LAMMPS input file
        with open(os.path.join(output_dir, self._lmp_input_fname), 'w') as f:

            f.write('# LAMMPS input file generated by PMD package\n')
            f.write('\n')

            f.write('### Initialization\n')
            f.write(f'{"atom_style":<15} full\n')
            f.write(f'{"units":<15} {self._units}\n')
            f.write('\n')

            if self._force_field:
                self._force_field.write_settings(f)
            f.write('\n')

            if self._read_data_from:
                f.write(f'{"read_data":<15} {self._read_data_from}\n')
            # TODO: add last_restart_fname
            # elif self._read_restart_from:
            #   f.write(f'{"read_restart":<15} '
            #           f'{self._read_restart_from.last_restart_fname}\n')
            f.write('\n')

            f.write(f'{"neighbor":<15} {self._neighbor_skin} bin\n')
            f.write(f'{"neigh_modify":<15} delay 0 every '
                    f'{self._neighbor_every} check yes\n')
            f.write('\n')

            f.write(
                f'{"thermo_style":<15} custom step temp density vol press ke '
                f'pe ebond evdwl ecoul elong\n')
            f.write(f'{"thermo":<15} {self._thermo}\n')
            f.write(f'{"timestep":<15} {self._timestep}\n')

            for procedure in self._procedures:
                f.write('\n')
                f.write('\n')
                procedure.write_before_run(f)
                procedure.write_lammps(f)
                procedure.write_after_run(f)
                Pmdlogging.info(f'Adding procedure - {procedure} '
                                f'to {self._lmp_input_fname}')

        Pmdlogging.info(f'Lammps input file - {self._lmp_input_fname} '
                        f'successfully created in {output_dir}')
