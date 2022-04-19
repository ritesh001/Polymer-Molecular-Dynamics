from pmd.util import Util


class Job:
    '''pmd.core.Job.Job
    Template object to contain job initialization settings
    Attributes:
        jobname: str
            Job name
        project: str
            Project name
        nodes: int
            Number of nodes
        ppn: int
            Number of processors (CPU)
        gpus: int
            Number of processors (GPU)
        walltime: str
            Job walltime
        LAMMPS_EXEC: str
            Directory of the LAMMPS executable file
    '''
    def __init__(self, jobname, project, nodes, ppn, walltime, gpus=0):
        self.jobname = jobname
        self.project = project
        self.nodes = nodes
        self.ppn = ppn
        self.walltime = walltime
        self.gpus = gpus

    def write_pbs(self, output_dir):
        Util.build_dir(output_dir)
        pbs_fname = 'job.pbs'
        with open(output_dir + '/' + pbs_fname, 'w') as f:
            f.write('#PBS -A {}\n'.format(self.project))
            f.write('#PBS -q inferno\n')
            f.write('#PBS -N {}\n'.format(self.jobname))
            if self.gpus == 0:
                f.write('#PBS -l nodes={}:ppn={}\n'.format(self.nodes, self.ppn))
            else:
                f.write('#PBS -l nodes={}:ppn={}:gpus={}:RTX6000\n'.format(self.nodes, self.ppn, self.gpus))
            f.write('#PBS -l walltime={}\n'.format(self.walltime))
            f.write('#PBS -j oe\n')
            f.write('#PBS -o out.$PBS_JOBID\n')
            f.write('\n')
            f.write('cd $PBS_O_WORKDIR\n')
            if self.gpus == 0:
                f.write('module load intel/19.0.5 mvapich2/2.3.4 lammps/09Jan20\n')
            else:
                f.write('module load gcc/8.3.0 mvapich2/2.3.2 lammps-gpu/29Oct20\n')
            f.write('mpirun -np {} lmp -in lmp.in\n'.format(
                int(self.nodes * self.ppn)))
