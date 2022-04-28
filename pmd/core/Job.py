from pmd.util import Util


class Job:
    '''Template object to contain job initialization settings

    Attributes:
        jobname (str): Job name
        project (str): Project name
        nodes (int): Number of nodes
        ppn (int): Number of processors (CPU)
        walltime (str): Job walltime
        gpus (int): Number of processors (GPU)
    '''

    def __init__(self,
                 jobname: str,
                 project: str,
                 nodes: int,
                 ppn: int,
                 walltime: str,
                 gpus: int = 0):
        self.jobname = jobname
        self.project = project
        self.nodes = nodes
        self.ppn = ppn
        self.walltime = walltime
        self.gpus = gpus

    def write_pbs(self, output_dir: str, pbs_fname: str) -> None:
        '''Method to make the PBS job scheduler input file

        Parameters:
            output_dir (str): Directory for the PBS input file
            pbs_fname (str): Name of the PBS input file 

        Returns:
            None
        '''

        Util.build_dir(output_dir)
        with open(output_dir + '/' + pbs_fname, 'w') as f:
            f.write('#PBS -A {}\n'.format(self.project))
            f.write('#PBS -q inferno\n')
            f.write('#PBS -N {}\n'.format(self.jobname))
            if self.gpus == 0:
                f.write('#PBS -l nodes={}:ppn={}\n'.format(
                    self.nodes, self.ppn))
            else:
                f.write('#PBS -l nodes={}:ppn={}:gpus={}:RTX6000\n'.format(
                    self.nodes, self.ppn, self.gpus))
            f.write('#PBS -l walltime={}\n'.format(self.walltime))
            f.write('#PBS -j oe\n')
            f.write('#PBS -o out.$PBS_JOBID\n')
            f.write('\n')
            f.write('cd $PBS_O_WORKDIR\n')
            if self.gpus == 0:
                f.write(
                    'module load intel/19.0.5 mvapich2/2.3.4 lammps/09Jan20\n')
            else:
                f.write(
                    'module load gcc/8.3.0 mvapich2/2.3.2 lammps-gpu/29Oct20\n'
                )
            f.write('mpirun -np {} lmp -in lmp.in\n'.format(
                int(self.nodes * self.ppn)))

    def write_slurm(self, output_dir: str, slurm_fname: str) -> None:
        '''Method to make the Slurm job scheduler input file

        Parameters:
            output_dir (str): Directory for the Slurm input file
            slurm_fname (str): Name of the Slurm input file

        Returns:
            None
        '''

        Util.build_dir(output_dir)
        # Todo
