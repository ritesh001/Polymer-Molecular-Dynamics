---
sidebar_position: 1
title: Tg calculation
---

## Starting the simulation

The following script sets up the **TgMeasurement procedure** for polyethylene Tg calculation in a MD simulation:

```python
import pmd

system = pmd.System(smiles='*CC*', density=0.8,
                    natoms_total=5000, natoms_per_chain=150,
                    builder=pmd.EMC('opls-aa'))

lmp = pmd.Lammps(read_data_from=system)
lmp.add_procedure(pmd.Minimization())
lmp.add_procedure(pmd.Equilibration(Teq=600, Peq=1, Tmax=1000, Pmax=49346.163))
# highlight-next-line
lmp.add_procedure(pmd.TgMeasurement(Tinit=600, Tfinal=100, Tinterval=20, step_duration=10**6))

job = pmd.Torque(run_lammps=lmp, jobname='*CC*', project='GT-rramprasad3-CODA20',
                 nodes=2, ppn=24, walltime='72:00:00')

run = pmd.Pmd(system=system, lammps=lmp, job=job)
run.create()
```

## Result Analysis

The TgMeasurement procedure calcualte the average density at each temperature on the fly and dump the result to the `temp_vs_density.txt` file. Once the simulation is done, use the **Analysis module** to get the T<sub>g</sub> value:

```python
import pmd

Tg = pmd.calculate_Tg('temp_vs_density.txt')
# Do whatever you want with the Tg here (e.g. put the value into a csv file)
```

You can also do the same thing through command line

```bash
$ pmd-analyze -p Tg 'temp_vs_density.txt'  # this calls the calculate_Tg function
```

The `calculate_Tg` function not only returns the T<sub>g</sub> value of the system, but also make a result plot with the fitting lines like below:

![Tg result](/img/guides/temp_vs_density.png)
