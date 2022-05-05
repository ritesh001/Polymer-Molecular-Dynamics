from pmd.core.Lammps import Lammps
from pmd.util import Util
from abc import ABC, abstractmethod


class Job(ABC):

    @abstractmethod
    def write_job(self, output_dir: str):
        pass


class Torque(Job):
    '''Template Torque object to contain Torque job initialization settings

    Attributes:
        jobname (str): Job name
        project (str): Project name
        nodes (int): Number of nodes
        ppn (int): Number of processors (CPU)
        walltime (str): Job walltime
        gpus (int): Number of processors (GPU)
        job_fname (str): Name of the PBS input file; default: `job.pbs`
    '''

    def __init__(self,
                 run_lammps: Lammps,
                 jobname: str,
                 project: str,
                 nodes: int,
                 ppn: int,
                 walltime: str,
                 gpus: int = None,
                 job_fname: str = 'job.pbs'):
        self._run_lammps = run_lammps
        self._jobname = jobname
        self._project = project
        self._nodes = nodes
        self._ppn = ppn
        self._walltime = walltime
        self._gpus = gpus
        self._job_fname = job_fname

    def write_job(self, output_dir: str = '.') -> None:
        '''Method to make the Torque job scheduler input file

        Parameters:
            output_dir (str): Directory for the Torque input file; default: `.`

        Returns:
            None
        '''

        Util.build_dir(output_dir)
        with open(output_dir + '/' + self._job_fname, 'w') as f:
            f.write('#PBS -A {}\n'.format(self._project))
            f.write('#PBS -q inferno\n')
            f.write('#PBS -N {}\n'.format(self._jobname))
            if self._gpus:
                f.write(
                    '#PBS -l nodes={}:ppn={}:gpus={}:RTX6000:default\n'.format(
                        self._nodes, self._ppn, self._gpus))
            else:
                f.write('#PBS -l nodes={}:ppn={}\n'.format(
                    self._nodes, self._ppn))
            f.write('#PBS -l walltime={}\n'.format(self._walltime))
            f.write('#PBS -j oe\n')
            f.write('#PBS -o out.$PBS_JOBID\n')
            f.write('\n')
            f.write('cd $PBS_O_WORKDIR\n')
            if self._gpus:
                f.write(
                    'module load gcc/8.3.0 mvapich2/2.3.2 lammps-gpu/29Oct20\n'
                )
                f.write('mpirun -np {} lmp -sf gpu -pk gpu {} -in {}\n'.format(
                    self._nodes * self._ppn, self._gpus,
                    self._run_lammps.lmp_input_fname))
            else:
                f.write(
                    'module load intel/19.0.5 mvapich2/2.3.4 lammps/09Jan20\n')
                f.write('mpirun -np {} lmp -in {}\n'.format(
                    self._nodes * self._ppn, self._run_lammps.lmp_input_fname))


class Slurm(Job):
    '''Template Slurm object to contain Slurm job initialization settings

    Attributes:
        jobname (str): Job name
        nodes (int): Number of nodes
        ntasks_per_node (int): Number of processors (CPU)
        time (str): Job time
        gpus (int): Number of processors (GPU)
        job_fname (str): Name of the Slurm input file; default: `job.sh`
    '''

    def __init__(self,
                 run_lammps: Lammps,
                 jobname: str,
                 nodes: int,
                 ntasks_per_node: int,
                 time: str,
                 gpus: int = 0,
                 job_fname: str = 'job.sh'):
        self._run_lammps = run_lammps
        self._jobname = jobname
        self._nodes = nodes
        self._ntasks_per_node = ntasks_per_node
        self._time = time
        self._gpus = gpus
        self._job_fname = job_fname

    def write_job(self, output_dir: str = '.') -> None:
        '''Method to make the Slurm job scheduler input file

        Parameters:
            output_dir (str): Directory for the Slurm input file; default: `.`

        Returns:
            None
        '''
        Util.build_dir(output_dir)
        with open(output_dir + '/' + self._job_fname, 'w') as f:
            f.write('#!/bin/bash')
            f.write('#SBATCH --job-name={}\n'.format(self._jobname))
            f.write('#SBATCH -o out.o%j \n')
            f.write('#SBATCH -e err.e%j \n')
            f.write('#SBATCH --nodes={}\n'.format(self._nodes))
            if self._gpus:
                f.write('#SBATCH --gpus={}\n'.format(self._gpus))
            else:
                f.write('#SBATCH --ntasks-per-node={}\n'.format(
                    self._ntasks_per_node))
            f.write('#SBATCH --time={}\n'.format(self._time))
            if self._gpus:
                # TODO
                print('Have not implemented GPU Slurm yet')
            else:
                f.write('module load intel/18.0.2 impi/18.0.2 lammps/9Jan20\n')
                f.write('ibrun lmp -in {}\n'.format(
                    self._run_lammps.lmp_input_fname))
