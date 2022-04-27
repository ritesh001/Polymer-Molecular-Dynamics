from abc import ABC, abstractmethod
from io import TextIOWrapper


class Procedure(ABC):

    @abstractmethod
    def write_input(self):
        pass


class Minimization(Procedure):

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
