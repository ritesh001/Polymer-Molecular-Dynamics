---
sidebar_position: 1
title: Installation
---

The main package used here is the **Polymer Molecular Dynamics (PMD) package**, which can build LAMMPS input files, PBS/Slurm job files, and analyze results for various MD calculations.

However, to run MD simulations, we also need the initial configuration and force field parameters of the system. The System objects in PMD utilizes the **Polymer Structure Predictor (PSP) package** to prepare those.

Belows are the installation instruction on these two packages:

## PMD

To install PMD, simply do

```bash
$ pip install pmd
```

## PSP

:::tip Tip
If you're using PMD only for creating LAMMPS input files or PBS/Slurm job files, PSP is not required.
:::

### For general users

1. Follow the installation instruction on [PSP GitHub Repo](https://github.com/Ramprasad-Group/PSP).

2. Unfortunately, there are quite a few packages needed to be installed for PSP to work, make sure that you install all the packages required for PSP.

### For Ramprasad group users

1. On Gaanam or Tyrion clusters, you can use the python version at
   `/home/modules/anaconda3/bin/python3`, which has PSP already installed.

2. To fully utilize PSP, you need the following PATHs to your `~/.bashrc`

```bash title="~/.bashrc"
# PACKMOL PATH
export PACKMOL_EXEC='/home/appls/utility/packmol/packmol'

# Pysimm PATHs
export PYTHONPATH=$PYTHONPATH:'/home/appls/utility/pysimm'
export PATH=$PATH:'/home/appls/utility/pysimm/bin'

# AmberTools PATH
export ANTECHAMBER_EXEC='/home/modules/anaconda3/bin/antechamber'

# OBabel PATH
export PATH=$PATH:'/home/modules/anaconda3/bin'

# BOSS PATH
export BOSSdir='/home/appls/utility/boss'
```
