# Tutorial to generate input files for equilibration

Before you start, make sure PSP and PMD are installed.

As an example, say we want to equilibrate a simple polyethylene system

## Build the system
Run: `python mkinput.py`

The `mkinput.py` defines the `System`, `Lammps`, and `Job` objects as well as creates all the required data and input files for LAMMPS equilibration simulation.

This should generate a folder called `PE_equilibration` with all the necessary files for running MD simulation on your cluster.

## Run the simulation
For Ramprasad group users, simply upload the folder on PACE and `cd` to the directory and do
```bash
$ qsub job.pbs
```

Alternatively, you can run it locally on Gaanam clusters by
```bash
$ mpirun -np 4 /home/appl/utility/lammps/src/lmp_mpi -in lmp.in
```
