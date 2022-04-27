from abc import ABC, abstractmethod
from io import TextIOWrapper


class Procedure(ABC):

    @abstractmethod
    def write_input(self):
        pass


class Minimization(Procedure):
    '''Perform an energy minimization of the system, by iteratively adjusting atom coordinates. 
    Iterations are terminated when one of the stopping criteria is satisfied. At that point the 
    configuration will hopefully be in local potential energy minimum.
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

    def write_input(self, f: TextIOWrapper):
        f.write('### Minimization\n')
        f.write('{:<15} {}\n'.format('min_style', self._min_style))
        f.write('{:<15} {} {} {} {}\n'.format('minimize', self._etol,
                                              self._ftol, self._maxiter,
                                              self._maxeval))
        f.write('{:<15} 0\n'.format('reset_timestep'))
        f.write('\n')
        f.write('\n')


class Equilibration(Procedure):
    '''Perform a 21-step amorphous polymer equilibration process. 
       Ref: Abbott, Hart, and Colina, Theoretical Chemistry Accounts, 132(3), 1-19, 2013.
    '''

    def __init__(self,
                 Tfinal=300,
                 Pfinal=1,
                 Tmax=1000,
                 Pmax=50000,
                 Tdamp='$(100.0*dt)',
                 Pdamp='$(100.0*dt)',
                 reset_timestep=True):
        self._Tfinal = Tfinal
        self._Pfinal = Pfinal
        self._Tmax = Tmax
        self._Pmax = Pmax
        self._Tdamp = Tdamp
        self._Pdamp = Pdamp
        self._reset_timestep = reset_timestep
        self._eq_totaltime = 0
        self._eq_step = [
            ['nvt', 50000, Tmax],
            ['nvt', 50000, Tfinal],
            ['npt', 50000, Tfinal, 0.02 * Pmax],
            ['nvt', 50000, Tmax],
            ['nvt', 100000, Tfinal],
            ['npt', 50000, Tfinal, 0.6 * Pmax],
            ['nvt', 50000, Tmax],
            ['nvt', 100000, Tfinal],
            ['npt', 50000, Tfinal, Pmax],
            ['nvt', 50000, Tmax],
            ['nvt', 100000, Tfinal],
            ['npt', 5000, Tfinal, 0.5 * Pmax],
            ['nvt', 5000, Tmax],
            ['nvt', 10000, Tfinal],
            ['npt', 5000, Tfinal, 0.1 * Pmax],
            ['nvt', 5000, Tmax],
            ['nvt', 10000, Tfinal],
            ['npt', 5000, Tfinal, 0.01 * Pmax],
            ['nvt', 5000, Tmax],
            ['nvt', 10000, Tfinal],
            ['npt', 800000, Tfinal, Pfinal],
        ]

        for i in self._eq_step:
            self._eq_totaltime += i[1]

    def write_input(self, f: TextIOWrapper):
        f.write('### Equilibration\n')
        f.write(
            '{:<15} dump_eq all custom 10000 equil.lammpstrj id mol type q xs ys zs ix iy iz\n'
            .format('dump'))
        f.write('{:<15} {} equilibrated.restart\n'.format(
            'restart', self._eq_totaltime))
        f.write('\n')

        for n, i in enumerate(self._eq_step):
            if i[0] == 'nvt':
                f.write('{:<15} step{} all nvt temp {} {} {}\n'.format(
                    'fix', n + 1, i[2], i[2], self._Tdamp))
            elif i[0] == 'npt':
                f.write('{:<15} step{} all npt temp {} {} {} iso {} {} {}\n'.
                        format('fix', n + 1, i[2], i[2], self._Tdamp, i[3],
                               i[3], self._Pdamp))
            f.write('{:<15} {}\n'.format('run', i[1]))
            f.write('{:<15} step{}\n'.format('unfix', n + 1))
            f.write('\n')
        f.write('{:<15} dump_eq\n'.format('undump'))
        if self._reset_timestep:
            f.write('{:<15} 0\n'.format('reset_timestep'))
        f.write('\n')
        f.write('\n')


class TgMeasurement(Procedure):
    '''Perform glass transition temperature measurement of the system, 
    by iteratively cooling the system and equilibrate.
    '''

    def __init__(self,
                 Tinit=500,
                 Tfinal=100,
                 Tinterval=20,
                 step=1000000,
                 pressure=1,
                 Tdamp='$(100.0*dt)',
                 Pdamp='$(100.0*dt)'):
        self._Tinit = Tinit
        self._Tfinal = Tfinal
        self._Tinterval = Tinterval
        self._step = step
        self._pressure = pressure
        self._Tdamp = Tdamp
        self._Pdamp = Pdamp

    def write_input(self, f: TextIOWrapper):
        f.write('### Production - Tg measurement\n')
        f.write(
            '{:<15} dump_Tg all custom 10000 production.lammpstrj id mol type q xs ys zs ix iy iz\n'
            .format('dump'))
        f.write('{:<15} {} production.restart\n'.format('restart', self._step))
        f.write('{:<15} Rho equal density\n'.format('variable'))
        f.write('{:<15} Temp equal temp\n'.format('variable'))
        f.write(
            '{:<15} fDENS all ave/time {} {} {} v_Temp v_Rho file temp_vs_density.txt\n'
            .format('fix', int(self._step / 100 / 4), 100, self._step))
        f.write('\n')

        f.write('{:<15} loop\n'.format('label'))
        f.write('{:<15} a loop {}\n'.format(
            'variable',
            int((self._Tinit - self._Tfinal) / self._Tinterval + 1)))
        f.write('{:<15} b equal {}-{}*($a-1)\n'.format('variable', self._Tinit,
                                                       self._Tinterval))
        f.write('{:<15} fNPT all npt temp $b $b {} iso {} {} {}\n'.format(
            'fix', self._Tdamp, self._pressure, self._pressure, self._Pdamp))
        f.write('{:<15} {}\n'.format('run', self._step))
        f.write('{:<15} fNPT\n'.format('unfix'))
        f.write('{:<15} a\n'.format('next'))
        f.write('{:<15} SELF loop\n'.format('jump'))
        f.write('{:<15} a delete\n'.format('variable'))
        f.write('\n')
        f.write('{:<15} dump_Tg\n'.format('undump'))
        f.write('\n')
        f.write('\n')
