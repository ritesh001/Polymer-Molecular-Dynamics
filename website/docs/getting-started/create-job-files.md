---
sidebar_position: 4
title: Creating job files
---

To create a job scheduler input file, first we need to define a Torque or Slurm object depending on your cluster's scheduling system.

## Define the job for running LAMMPS

As an example, let us instantiate a Torque object:

```python
import pmd

job = pmd.Torque(run_lammps=lmp, jobname='Your-job-name', project='Your-project-id',
                 nodes=2, ppn=24, walltime='24:00:00')
```

You will first need to define a Lammps object (e.g. `lmp` shown above) and provide it as the `run_lammps` argument input. This is because the job would need to know the Lammps' input file name to work.

## Define the job for running Python

You can also create jobs that run Python script

```python
import pmd

# highlight-next-line
job = pmd.Torque(run_python='Your-python-script.py',
                 jobname='Your-job-name',
                 project='Your-project-id',
                 nodes=2, ppn=24, walltime='24:00:00')
```

And of course, you will need to make sure the Python script is in the same directory as the generated job file.
