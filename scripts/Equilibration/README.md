# Tutorial to generate input files for equilibration

Before you start, make sure PSP and PMD are installed.

As an example, say we want to equilibrate a simple polyethylene system

Run: `python mkinput.py`

The `mkinput.py` defines the `System`, `Lammps`, and `Job` objects as well as creates all the required data and input files for LAMMPS equilibration simulation.

This should generate a folder called `simple equilibration` with all the necessary files for running MD simulation on PACE. Simply upload the folder on PACE and `cd` to the directory and do `qsub job.pbs`.
