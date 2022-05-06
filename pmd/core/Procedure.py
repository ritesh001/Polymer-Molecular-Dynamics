from abc import ABC, abstractmethod
from io import TextIOWrapper


class Procedure(ABC):

    @abstractmethod
    def write_lammps(self):
        pass


class Minimization(Procedure):
    '''Perform an energy minimization of the system, by iteratively adjusting
    atom coordinates. Iterations are terminated when one of the stopping 
    criteria is satisfied. At that point the configuration will hopefully be in
    local potential energy minimum.

    Attributes:
        min_style (str): Minimization algorithm, see 
                         [here](https://docs.lammps.org/min_style.html) for all
                         options; default: `cg`   
                         
        etol (float): Stopping tolerance for energy (unitless); default: 
                      `10**(-6)`
        
        ftol (float): Stopping tolerance for force (force units); default: 
                      `10**(-8)`
        
        maxiter (int): Max iterations of minimizer; default: `10**5`
        
        maxeval (int): Max number of force/energy evaluations; default: `10**7`
    '''

    def __init__(self,
                 min_style: str = 'cg',
                 etol: float = 10**(-6),
                 ftol: float = 10**(-8),
                 maxiter: int = 10**5,
                 maxeval: int = 10**7):
        self._min_style = min_style
        self._etol = etol
        self._ftol = ftol
        self._maxiter = maxiter
        self._maxeval = maxeval

    def write_lammps(self, f: TextIOWrapper):
        f.write('### Minimization\n')
        f.write(f'{"min_style":<15} {self._min_style}\n')
        f.write(f'{"minimize":<15} {self._etol} {self._ftol} '
                f'{self._maxiter} {self._maxeval}\n')
        f.write('\n')
        f.write('\n')


class Equilibration(Procedure):
    '''Perform a 21-step amorphous polymer equilibration process. Ref: Abbott, 
    Hart, and Colina, Theoretical Chemistry Accounts, 132(3), 1-19, 2013.
       
    Attributes:
        Teq (float): Target equilibration temperature; default: `300`
        
        Peq (float): Target equilibration pressure; default: `1`
        
        Tmax (float): Maximum temperature during the equilibration; default: 
                      `600`
        
        Pmax (float): Maximum pressure during the equilibration; default: 
                      `50000`
        
        Tdamp (str): Damping parameter for thermostats; default: `$(100.0*dt)`
        
        Pdamp (str): Damping parameter for barostats; default: `$(100.0*dt)`
        
        dump_fname (str): Name of the dump file; default: `equil.lammpstrj`
        
        reset_timestep_before_run (bool): Whether to reset timestep after the 
                                          procedure; default: `True`
    '''

    def __init__(self,
                 Teq: float = 300,
                 Peq: float = 1,
                 Tmax: float = 600,
                 Pmax: float = 50000,
                 Tdamp: str = '$(100.0*dt)',
                 Pdamp: str = '$(100.0*dt)',
                 dump_fname: str = 'equil.lammpstrj',
                 reset_timestep_before_run: bool = True):
        self._Tdamp = Tdamp
        self._Pdamp = Pdamp
        self._dump_fname = dump_fname
        self._reset_timestep_before_run = reset_timestep_before_run
        self._eq_totaltime = 0
        self._eq_step = [
            ['nvt', 50000, Tmax],
            ['nvt', 50000, Teq],
            ['npt', 50000, Teq, 0.02 * Pmax],
            ['nvt', 50000, Tmax],
            ['nvt', 100000, Teq],
            ['npt', 50000, Teq, 0.6 * Pmax],
            ['nvt', 50000, Tmax],
            ['nvt', 100000, Teq],
            ['npt', 50000, Teq, Pmax],
            ['nvt', 50000, Tmax],
            ['nvt', 100000, Teq],
            ['npt', 5000, Teq, 0.5 * Pmax],
            ['nvt', 5000, Tmax],
            ['nvt', 10000, Teq],
            ['npt', 5000, Teq, 0.1 * Pmax],
            ['nvt', 5000, Tmax],
            ['nvt', 10000, Teq],
            ['npt', 5000, Teq, 0.01 * Pmax],
            ['nvt', 5000, Tmax],
            ['nvt', 10000, Teq],
            ['npt', 800000, Teq, Peq],
        ]

        for i in self._eq_step:
            self._eq_totaltime += i[1]

    def write_lammps(self, f: TextIOWrapper):
        f.write('### Equilibration\n')
        if self._reset_timestep_before_run:
            f.write(f'{"reset_timestep":<15} 0\n')
        f.write('\n')
        f.write(f'{"dump":<15} dump_eq all custom 10000 {self._dump_fname} '
                f'id mol type q xs ys zs ix iy iz\n')
        f.write(f'{"restart":<15} {self._eq_totaltime} equilibrated.restart\n')
        f.write('\n')

        for n, i in enumerate(self._eq_step):
            if i[0] == 'nvt':
                f.write(f'{"fix":<15} step{n + 1} all nvt temp '
                        f'{i[2]} {i[2]} {self._Tdamp}\n')
            elif i[0] == 'npt':
                f.write(f'{"fix":<15} step{n + 1} all npt temp {i[2]} {i[2]} '
                        f'{self._Tdamp} iso {i[3]} {i[3]} {self._Pdamp}\n')
            f.write(f'{"run":<15} {i[1]}\n')
            f.write(f'{"unfix":<15} step{n + 1}\n')
            f.write('\n')
        f.write(f'{"undump":<15} dump_eq\n')
        f.write('\n')
        f.write('\n')


class NPT(Procedure):
    '''Perform the simulation under NPT ensemble (via Nose-Hoover thermostat
    and barostat).
       
    Attributes:
        duration (int): Duration of this NPT procedure (timestep unit)
    
        Tinit (float): Initial temperature
        
        Tfinal (float): Final temperature
        
        Pinit (float): Initial pressure
        
        Pfinal (float): Final pressure
        
        Tdamp (str): Damping parameter for thermostats; default: `$(100.0*dt)`
        
        Pdamp (str): Damping parameter for barostats; default: `$(100.0*dt)`
        
        dump_fname (str): Name of the dump file; default: `npt.lammpstrj`
        
        reset_timestep_before_run (bool): Whether to reset timestep after the 
                                          procedure; default: `False`
    '''

    def __init__(self,
                 duration: int,
                 Tinit: float,
                 Tfinal: float,
                 Pinit: float,
                 Pfinal: float,
                 Tdamp: str = '$(100.0*dt)',
                 Pdamp: str = '$(100.0*dt)',
                 dump_fname: str = 'npt.lammpstrj',
                 reset_timestep_before_run: bool = False):
        self._duration = duration
        self._Tinit = Tinit
        self._Tfinal = Tfinal
        self._Pinit = Pinit
        self._Pfinal = Pfinal
        self._Tdamp = Tdamp
        self._Pdamp = Pdamp
        self._dump_fname = dump_fname
        self._reset_timestep_before_run = reset_timestep_before_run

    def write_lammps(self, f: TextIOWrapper):
        f.write('### NPT simulation\n')
        if self._reset_timestep_before_run:
            f.write(f'{"reset_timestep":<15} 0\n')
        f.write('\n')
        f.write(f'{"dump":<15} dump_npt all custom 10000 {self._dump_fname} '
                f'id mol type q xs ys zs ix iy iz\n')
        f.write(f'{"restart":<15} {self._duration} npt.restart\n')
        f.write('\n')
        f.write(
            f'{"fix":<15} fNPT all npt temp {self._Tinit} {self._Tfinal} '
            f'{self._Tdamp} iso {self._Pinit} {self._Pfinal} {self._Pdamp}\n')
        f.write(f'{"run":<15} {self._duration}\n')
        f.write(f'{"unfix":<15} fNPT\n')
        f.write('\n')
        f.write(f'{"undump":<15} dump_npt\n')
        f.write('\n')
        f.write('\n')


class NVT(Procedure):
    '''Perform the simulation under NVT ensemble (via Nose-Hoover thermostat).
       
    Attributes:
        duration (int): Duration of this NVT procedure (timestep unit)
    
        Tinit (float): Initial temperature
        
        Tfinal (float): Final temperature
        
        Tdamp (str): Damping parameter for thermostats; default: `$(100.0*dt)`
        
        dump_fname (str): Name of the dump file; default: `nvt.lammpstrj`
        
        reset_timestep_before_run (bool): Whether to reset timestep after the 
                                          procedure; default: `False`
    '''

    def __init__(self,
                 duration: int,
                 Tinit: float,
                 Tfinal: float,
                 Tdamp: str = '$(100.0*dt)',
                 dump_fname: str = 'nvt.lammpstrj',
                 reset_timestep_before_run: bool = False):
        self._duration = duration
        self._Tinit = Tinit
        self._Tfinal = Tfinal
        self._Tdamp = Tdamp
        self._dump_fname = dump_fname
        self._reset_timestep_before_run = reset_timestep_before_run

    def write_lammps(self, f: TextIOWrapper):
        f.write('### NVT simulation\n')
        if self._reset_timestep_before_run:
            f.write(f'{"reset_timestep":<15} 0\n')
        f.write('\n')
        f.write(f'{"dump":<15} dump_nvt all custom 10000 {self._dump_fname} '
                f'id mol type q xs ys zs ix iy iz\n')
        f.write(f'{"dump":<15} dump_image all image {self._duration} '
                f'image.*.jpg type type\n')
        f.write(f'{"restart":<15} {self._duration} nvt.restart\n')
        f.write('\n')
        f.write(f'{"fix":<15} fNVT all nvt temp {self._Tinit} '
                f'{self._Tfinal} {self._Tdamp}\n')
        f.write(f'{"run":<15} {self._duration}\n')
        f.write(f'{"unfix":<15} fNVT\n')
        f.write('\n')
        f.write(f'{"undump":<15} dump_nvt\n')
        f.write(f'{"undump":<15} dump_image\n')
        f.write('\n')
        f.write('\n')


class TgMeasurement(Procedure):
    '''Perform glass transition temperature measurement of the system, 
    by iteratively cooling the system and equilibrate.
       
    Attributes:    
        Tinit (float): Initial temperature of the cooling process; default: 
                       `500`
        
        Tfinal (float): Final temperature of the cooling process; default: 
                        `100`
        
        Tinterval (float): Temperature interval of the cooling process
                           ; default: `20`
        
        step_duration (int): Duration of each temperature step 
                             (timestep unit); default: `1000000`
        
        pressure (float): Pressure during the cooling process; default: `1`
        
        Tdamp (str): Damping parameter for thermostats; default: `$(100.0*dt)`
        
        Pdamp (str): Damping parameter for barostats; default: `$(100.0*dt)`
        
        dump_fname (str): Name of the dump file; default: 
                          `Tg_measurement.lammpstrj`
        
        result_fname (str): Name of the result file; default: 
                            `temp_vs_density.txt`
    '''

    def __init__(self,
                 Tinit: float = 500,
                 Tfinal: float = 100,
                 Tinterval: float = 20,
                 step_duration: int = 1000000,
                 pressure: float = 1,
                 Tdamp: str = '$(100.0*dt)',
                 Pdamp: str = '$(100.0*dt)',
                 dump_fname: str = 'Tg_measurement.lammpstrj',
                 result_fname: str = 'temp_vs_density.txt'):
        self._Tinit = Tinit
        self._Tfinal = Tfinal
        self._Tinterval = Tinterval
        self._step_duration = step_duration
        self._pressure = pressure
        self._Tdamp = Tdamp
        self._Pdamp = Pdamp
        self._dump_fname = dump_fname
        self._result_fname = result_fname

    def write_lammps(self, f: TextIOWrapper):
        f.write('### Production - Tg measurement\n')
        f.write(f'{"dump":<15} dump_Tg all custom 10000 {self._dump_fname} '
                f'id mol type q xs ys zs ix iy iz\n')
        f.write(f'{"restart":<15} {self._step_duration} production.restart\n')
        f.write(f'{"variable":<15} Rho equal density\n')
        f.write(f'{"variable":<15} Temp equal temp\n')
        f.write(
            f'{"fix":<15} fDENS all ave/time '
            f'{int(self._step_duration / 100 / 4)} {100} '
            f'{self._step_duration} v_Temp v_Rho file {self._result_fname}\n')
        f.write('\n')

        f.write(f'{"label":<15} loop\n')
        f.write(f'{"variable":<15} a loop '
                f'{int((self._Tinit - self._Tfinal) / self._Tinterval + 1)}\n')
        f.write(f'{"variable":<15} b equal '
                f'{self._Tinit}-{self._Tinterval}*($a-1)\n')
        f.write(f'{"fix":<15} fNPT all npt temp $b $b {self._Tdamp} iso '
                f'{self._pressure} {self._pressure} {self._Pdamp}\n')
        f.write(f'{"run":<15} {self._step_duration}\n')
        f.write(f'{"unfix":<15} fNPT\n')
        f.write(f'{"next":<15} a\n')
        f.write(f'{"jump":<15} SELF loop\n')
        f.write(f'{"variable":<15} a delete\n')
        f.write('\n')
        f.write(f'{"undump":<15} dump_Tg\n')
        f.write('\n')
        f.write('\n')
